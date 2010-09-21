/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/preview.h>
#include <mitsuba/render/triaccel.h>
#include <mitsuba/render/imageproc_wu.h>

MTS_NAMESPACE_BEGIN

void PreviewWorker::serialize(Stream *stream, InstanceManager *manager) const {
	Log(EError, "Serialization is not supported!");
}

ref<WorkUnit> PreviewWorker::createWorkUnit() const {
	return new RectangularWorkUnit();
}

ref<WorkResult> PreviewWorker::createWorkResult() const {
	/* Only pixel data, nothing else */
	return new ImageBlock(Vector2i(m_blockSize, m_blockSize), 
		0, false, false, false, false);
}

void PreviewWorker::prepare() {
	m_scene = static_cast<Scene *>(getResource("scene"));
	m_kdtree = m_scene->getKDTree();
	m_shapes = &m_kdtree->getShapes();
}

void PreviewWorker::process(const WorkUnit *workUnit, WorkResult *workResult, 
	const bool &stop) {
	if (m_coherent)
		processCoherent(workUnit, workResult, stop);
	else
		processIncoherent(workUnit, workResult, stop);
}

void PreviewWorker::processIncoherent(const WorkUnit *workUnit, WorkResult *workResult, 
	const bool &stop) {
	const RectangularWorkUnit *rect = static_cast<const RectangularWorkUnit *>(workUnit);
	ImageBlock *block = static_cast<ImageBlock *>(workResult);

	block->setOffset(rect->getOffset());
	block->setSize(rect->getSize());

	const int sx = rect->getOffset().x, sy = block->getOffset().y;
	const int ex = sx + rect->getSize().x, ey = sy + rect->getSize().y;

	/* Some local variables */
	int pos = 0;
	Intersection its;
	Spectrum value, bsdfVal;
	Vector toIts;
	Ray primary, secondary;
	int numRays = 0;

	for (int y=sy; y<ey; ++y) {
		for (int x=sx; x<ex; ++x) {
			/* Generate a camera ray without normalization */
			primary = Ray(m_cameraO, m_cameraTL + m_cameraDx*x + m_cameraDy * y);

			++numRays;
			if (!m_kdtree->rayIntersect(primary, its)) {
				block->setPixel(pos++, m_scene->LeBackground(primary));
				continue;
			}

			if (its.shape->isLuminaire())
				value = its.Le(-primary.d);
			else
				value = Spectrum(0.0f);

			toIts = its.p - m_vpl.its.p;
			secondary = Ray(m_vpl.its.p, toIts, 0.001, 0.999);
			++numRays;
			if (m_kdtree->rayIntersect(secondary)) {
				block->setPixel(pos++, value);
				continue;
			}

			Float length = toIts.length();
			toIts /= length;

			BSDFQueryRecord rr(its, -its.toLocal(toIts));
			rr.wi = normalize(rr.wi);
			bsdfVal = its.shape->getBSDF()->fCos(rr);
			length = std::max(length, m_minDist);

			if (m_vpl.type == ESurfaceVPL) {
				BSDFQueryRecord bRec(m_vpl.its, m_vpl.its.toLocal(toIts));
				bRec.quantity = EImportance;
				value += m_vpl.P * bsdfVal * m_vpl.its.shape->getBSDF()->fCos(bRec) / (length*length);
			} else {
				EmissionRecord eRec(m_vpl.luminaire, 
					ShapeSamplingRecord(m_vpl.its.p, m_vpl.its.shFrame.n), toIts);
				eRec.type = EmissionRecord::EPreview;
				value += m_vpl.P * bsdfVal * m_vpl.luminaire->f(eRec) 
					* ((m_vpl.luminaire->getType() & Luminaire::EOnSurface ?
					dot(m_vpl.its.shFrame.n, toIts) : (Float) 1)
					/ (length*length));
			}
			block->setPixel(pos++, value);
		}
	}
	block->setExtra(numRays);
}

void PreviewWorker::processCoherent(const WorkUnit *workUnit, WorkResult *workResult, 
	const bool &stop) {
#if defined(MTS_HAS_COHERENT_RT)
	const RectangularWorkUnit *rect = static_cast<const RectangularWorkUnit *>(workUnit);
	ImageBlock *block = static_cast<ImageBlock *>(workResult);

	block->setOffset(rect->getOffset());
	block->setSize(rect->getSize());

	/* Some constants */
	const int sx = rect->getOffset().x, sy = block->getOffset().y;
	const int ex = sx + rect->getSize().x, ey = sy + rect->getSize().y;
	const int width = rect->getSize().x;
	const SSEVector MM_ALIGN16 xOffset(0.0f, 1.0f, 0.0f, 1.0f);
	const SSEVector MM_ALIGN16 yOffset(0.0f, 0.0f, 1.0f, 1.0f);
	const int pixelOffset[] = {0, 1, width, width+1};
	const __m128 clamping = _mm_set1_ps(1/(m_minDist*m_minDist));

	const __m128 camTL[3] = {
		 _mm_set1_ps(m_cameraTL.x),
		 _mm_set1_ps(m_cameraTL.y),
		 _mm_set1_ps(m_cameraTL.z)
	}; 
	const __m128 camDx[3] = {
		 _mm_set1_ps(m_cameraDx.x),
		 _mm_set1_ps(m_cameraDx.y),
		 _mm_set1_ps(m_cameraDx.z)
	}; 
	const __m128 camDy[3] = {
		 _mm_set1_ps(m_cameraDy.x),
		 _mm_set1_ps(m_cameraDy.y),
		 _mm_set1_ps(m_cameraDy.z)
	}; 
	const __m128 lumPos[3] = {
		_mm_set1_ps(m_vpl.its.p.x),
		_mm_set1_ps(m_vpl.its.p.y),
		_mm_set1_ps(m_vpl.its.p.z)
	};
	const __m128 lumDir[3] = {
		_mm_set1_ps(m_vpl.its.shFrame.n.x),
		_mm_set1_ps(m_vpl.its.shFrame.n.y),
		_mm_set1_ps(m_vpl.its.shFrame.n.z)
	};

	/* Some local variables */
	int pos = 0;
	int numRays = 0;
	RayPacket4 MM_ALIGN16 primRay4, secRay4;
	Intersection4 MM_ALIGN16 its4, secIts4;
	RayInterval4 MM_ALIGN16 itv4, secItv4;
	SSEVector MM_ALIGN16 nSecD[3], cosThetaLight, invLengthSquared;
	Spectrum emitted[4], direct[4];
	Intersection its;
	LuminaireSamplingRecord lRec;
	Vector wo, wi;
	its.hasUVPartials = false;

	bool diffuseVPL = false, vplOnSurface = false;
	Spectrum vplWeight;

	if (m_vpl.type == ESurfaceVPL && (m_diffuseSources || m_vpl.its.shape->getBSDF()->getType() == BSDF::EDiffuseReflection)) {
		diffuseVPL = true;
		vplOnSurface = true;
		vplWeight = m_vpl.its.shape->getBSDF()->getDiffuseReflectance(m_vpl.its) * m_vpl.P / M_PI;
	} else if (m_vpl.type == ELuminaireVPL) {
		vplOnSurface = m_vpl.luminaire->getType() & Luminaire::EOnSurface;
		diffuseVPL = m_vpl.luminaire->getType() & Luminaire::EDiffuseDirection;
		EmissionRecord eRec(m_vpl.luminaire, 
			ShapeSamplingRecord(m_vpl.its.p, m_vpl.its.shFrame.n), m_vpl.its.shFrame.n);
		vplWeight = m_vpl.P * m_vpl.luminaire->f(eRec);
	}

	primRay4.o[0].ps = _mm_set1_ps(m_cameraO.x);
	primRay4.o[1].ps = _mm_set1_ps(m_cameraO.y);
	primRay4.o[2].ps = _mm_set1_ps(m_cameraO.z);
	secItv4.mint.ps = _mm_set1_ps(0.001);

	/* Work on 2x2 sub-blocks */
	for (int y=sy; y<ey; y += 2, pos += width) {
		for (int x=sx; x<ex; x += 2, pos += 2) {
			/* Generate camera rays without normalization */
			const __m128
				xPixel = _mm_add_ps(xOffset.ps, _mm_set1_ps((float) x)),
				yPixel = _mm_add_ps(yOffset.ps, _mm_set1_ps((float) y));

			primRay4.d[0].ps = _mm_add_ps(camTL[0], _mm_add_ps(
				_mm_mul_ps(xPixel, camDx[0]), _mm_mul_ps(yPixel, camDy[0])));
			primRay4.d[1].ps = _mm_add_ps(camTL[1], _mm_add_ps(
				_mm_mul_ps(xPixel, camDx[1]), _mm_mul_ps(yPixel, camDy[1])));
			primRay4.d[2].ps = _mm_add_ps(camTL[2], _mm_add_ps(
				_mm_mul_ps(xPixel, camDx[2]), _mm_mul_ps(yPixel, camDy[2])));

			primRay4.dRcp[0].ps = _mm_div_ps(SSEConstants::one.ps, primRay4.d[0].ps);
			primRay4.dRcp[1].ps = _mm_div_ps(SSEConstants::one.ps, primRay4.d[1].ps);
			primRay4.dRcp[2].ps = _mm_div_ps(SSEConstants::one.ps, primRay4.d[2].ps);

			/* Ray coherence test */
			const int primSignsX = _mm_movemask_ps(primRay4.d[0].ps);
			const int primSignsY = _mm_movemask_ps(primRay4.d[1].ps);
			const int primSignsZ = _mm_movemask_ps(primRay4.d[2].ps);

			const bool primCoherent =
				   (primSignsX == 0 || primSignsX == 0xF)
				&& (primSignsY == 0 || primSignsY == 0xF)
				&& (primSignsZ == 0 || primSignsZ == 0xF);

			/* Trace the primary rays */
			its4.t = SSEConstants::p_inf;
			if (EXPECT_TAKEN(primCoherent)) {
				primRay4.signs[0][0] = primSignsX ? 1 : 0;
				primRay4.signs[1][0] = primSignsY ? 1 : 0;
				primRay4.signs[2][0] = primSignsZ ? 1 : 0;
				m_kdtree->rayIntersectPacket(primRay4, itv4, its4);
			} else {
				m_kdtree->rayIntersectPacketIncoherent(primRay4, itv4, its4);
			}
			numRays += 4;

			/* Generate secondary rays */
			secRay4.o[0].ps = _mm_add_ps(primRay4.o[0].ps, _mm_mul_ps(its4.t.ps, primRay4.d[0].ps));
			secRay4.o[1].ps = _mm_add_ps(primRay4.o[1].ps, _mm_mul_ps(its4.t.ps, primRay4.d[1].ps));
			secRay4.o[2].ps = _mm_add_ps(primRay4.o[2].ps, _mm_mul_ps(its4.t.ps, primRay4.d[2].ps));
			secRay4.d[0].ps = _mm_sub_ps(lumPos[0], secRay4.o[0].ps);
			secRay4.d[1].ps = _mm_sub_ps(lumPos[1], secRay4.o[1].ps);
			secRay4.d[2].ps = _mm_sub_ps(lumPos[2], secRay4.o[2].ps);

			/* Normalization */
			const __m128 
				lengthSquared = _mm_add_ps(_mm_add_ps(
					_mm_mul_ps(secRay4.d[0].ps, secRay4.d[0].ps),
					_mm_mul_ps(secRay4.d[1].ps, secRay4.d[1].ps)),
					_mm_mul_ps(secRay4.d[2].ps, secRay4.d[2].ps)),
				invLength = _mm_rsqrt_ps(lengthSquared);
	
			invLengthSquared.ps = _mm_min_ps(_mm_rcp_ps(lengthSquared), clamping);

			nSecD[0].ps = _mm_mul_ps(secRay4.d[0].ps, invLength);
			nSecD[1].ps = _mm_mul_ps(secRay4.d[1].ps, invLength);
			nSecD[2].ps = _mm_mul_ps(secRay4.d[2].ps, invLength);

			secRay4.dRcp[0].ps = _mm_div_ps(SSEConstants::one.ps, secRay4.d[0].ps);
			secRay4.dRcp[1].ps = _mm_div_ps(SSEConstants::one.ps, secRay4.d[1].ps);
			secRay4.dRcp[2].ps = _mm_div_ps(SSEConstants::one.ps, secRay4.d[2].ps);
			
			cosThetaLight.ps = _mm_sub_ps(_mm_setzero_ps(),
				_mm_add_ps(_mm_add_ps(
					_mm_mul_ps(nSecD[0].ps, lumDir[0]),
					_mm_mul_ps(nSecD[1].ps, lumDir[1])),
					_mm_mul_ps(nSecD[2].ps, lumDir[2])));
			secItv4.maxt.ps = _mm_set1_ps(0.999);

			/* Shading (scalar) --- this is way too much work and should be 
			   rewritten to be smarter in special cases */
			for (int idx=0; idx<4; ++idx) {
				if (EXPECT_NOT_TAKEN(its4.t.f[idx] == std::numeric_limits<float>::infinity())) {
					/* Don't trace a secondary ray */
					secItv4.maxt.f[idx] = 0;
					emitted[idx] = m_scene->LeBackground(Ray(
						Point(primRay4.o[0].f[idx], primRay4.o[1].f[idx], primRay4.o[2].f[idx]),
						Vector(primRay4.d[0].f[idx], primRay4.d[1].f[idx], primRay4.d[2].f[idx])
					));
					memset(&direct[idx], 0, sizeof(Spectrum));
					continue;
				}
				const unsigned int primIndex = its4.primIndex.i[idx];
				const Shape *shape = (*m_shapes)[its4.shapeIndex.i[idx]];
				const BSDF *bsdf = shape->getBSDF();


				if (EXPECT_TAKEN(primIndex != KNoTriangleFlag)) {
					const TriMesh *mesh = static_cast<const TriMesh *>(shape);
					const Vertex *vb = mesh->getVertexBuffer();
					const Triangle &t = mesh->getTriangles()[primIndex];
					const Vertex &v0 = vb[t.idx[0]], &v1 = vb[t.idx[1]], &v2 = vb[t.idx[2]];
					const Float beta  = its4.u.f[idx],
								gamma = its4.v.f[idx],
								alpha = 1.0f - beta - gamma;
					its.shFrame.n = normalize(v0.n * alpha + v1.n * beta + v2.n * gamma);
					its.uv = v0.uv * alpha + v1.uv * beta + v2.uv * gamma;

					if (EXPECT_NOT_TAKEN(bsdf->getType() != BSDF::EDiffuseReflection || !diffuseVPL)) {
						its.dpdu = v0.dpdu * alpha + v1.dpdu * beta + v2.dpdu * gamma;
						its.dpdv = v0.dpdv * alpha + v1.dpdv * beta + v2.dpdv * gamma;
					}
				} else {
					Ray ray(
						Point(primRay4.o[0].f[idx], primRay4.o[1].f[idx], primRay4.o[2].f[idx]),
						Vector(primRay4.d[0].f[idx], primRay4.d[1].f[idx], primRay4.d[2].f[idx])
					);
					shape->rayIntersect(ray, its);
				}

				wo.x = nSecD[0].f[idx]; wo.y = nSecD[1].f[idx]; wo.z = nSecD[2].f[idx];

				if (EXPECT_TAKEN(!shape->isLuminaire())) {
					memset(&emitted[idx], 0, sizeof(Spectrum));
				} else {
					lRec.d = -Vector(primRay4.d[0].f[idx], primRay4.d[1].f[idx], primRay4.d[2].f[idx]);
					lRec.sRec.n = its.shFrame.n;
					emitted[idx] = shape->getLuminaire()->Le(lRec);
				}

				if (EXPECT_TAKEN(bsdf->getType() == BSDF::EDiffuseReflection && diffuseVPL)) {
					/* Fast path */
					direct[idx] = (bsdf->getDiffuseReflectance(its) * vplWeight)
						* (std::max((Float) 0.0f, dot(wo, its.shFrame.n))
						* (vplOnSurface ? (std::max(cosThetaLight.f[idx], (Float) 0.0f) * INV_PI) : INV_PI)
						* invLengthSquared.f[idx]);
				} else {
					wi.x = -primRay4.d[0].f[idx];
					wi.y = -primRay4.d[1].f[idx];
					wi.z = -primRay4.d[2].f[idx];
					its.p.x = secRay4.o[0].f[idx];
					its.p.y = secRay4.o[1].f[idx];
					its.p.z = secRay4.o[2].f[idx];
					its.shFrame.s = normalize(its.dpdu - its.shFrame.n
						* dot(its.shFrame.n, its.dpdu));
					its.shFrame.t = cross(its.shFrame.n, its.shFrame.s);
					const Float ctLight = cosThetaLight.f[idx];
					wi = normalize(wi);

					its.wi = its.toLocal(wi);
					wo = its.toLocal(wo);

					if (!diffuseVPL) {
						if (m_vpl.type == ESurfaceVPL) {
							BSDFQueryRecord bRec(m_vpl.its, m_vpl.its.toLocal(wi));
							bRec.quantity = EImportance;
							vplWeight = m_vpl.its.shape->getBSDF()->fCos(bRec) * m_vpl.P;
						} else {
							EmissionRecord eRec(m_vpl.luminaire, 
								ShapeSamplingRecord(m_vpl.its.p, m_vpl.its.shFrame.n), wi);
							eRec.type = EmissionRecord::EPreview;
							vplWeight = m_vpl.luminaire->f(eRec) * m_vpl.P;
						}
					}

					if (EXPECT_TAKEN(ctLight >= 0)) {
						direct[idx] = (bsdf->fCos(BSDFQueryRecord(its, wo)) * vplWeight
							* ((vplOnSurface ? std::max(ctLight, (Float) 0.0f) : 1.0f) * invLengthSquared.f[idx]));
					} else {
						memset(&direct[idx], 0, sizeof(Spectrum));
					}
				}
				++numRays;
			}

			/* Shoot the secondary rays */
			const int secSignsX = _mm_movemask_ps(secRay4.d[0].ps);
			const int secSignsY = _mm_movemask_ps(secRay4.d[1].ps);
			const int secSignsZ = _mm_movemask_ps(secRay4.d[2].ps);

			const bool secCoherent =
				   (secSignsX == 0 || secSignsX == 0xF)
				&& (secSignsY == 0 || secSignsY == 0xF)
				&& (secSignsZ == 0 || secSignsZ == 0xF);

			/* Shoot the secondary rays */
			secIts4.t = SSEConstants::p_inf;
			if (EXPECT_TAKEN(secCoherent)) {
				secRay4.signs[0][0] = secSignsX ? 1 : 0;
				secRay4.signs[1][0] = secSignsY ? 1 : 0;
				secRay4.signs[2][0] = secSignsZ ? 1 : 0;
				m_kdtree->rayIntersectPacket(secRay4, secItv4, secIts4);
			} else {
				m_kdtree->rayIntersectPacketIncoherent(secRay4, secItv4, secIts4);
			}

			for (int idx=0; idx<4; ++idx) {
				if (EXPECT_TAKEN(secIts4.t.f[idx] == std::numeric_limits<float>::infinity()))
					block->setPixel(pos+pixelOffset[idx], direct[idx]+emitted[idx]);
				else
					block->setPixel(pos+pixelOffset[idx], emitted[idx]);
			}
		}
	}
	block->setExtra(numRays);
#else
	Log(EError, "Coherent raytracing support was not compiled into this binary!");
#endif
}

ref<WorkProcessor> PreviewWorker::clone() const {
	return new PreviewWorker(m_blockSize, m_cameraO, m_cameraTL, 
		m_cameraDx, m_cameraDy, m_vpl, m_minDist, m_coherent,
		m_diffuseSources, m_diffuseReceivers);
}

MTS_IMPLEMENT_CLASS(PreviewWorker, false, WorkProcessor)
MTS_NAMESPACE_END
