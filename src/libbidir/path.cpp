/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/bidir/path.h>

MTS_NAMESPACE_BEGIN

void Path::initialize(const Scene *scene, Float time,
		ETransportMode mode, MemoryPool &pool) {
	release(pool);
	m_vertices.push_back(pool.allocVertex());
	m_vertices[0]->makeEndpoint(scene, time, mode);
}

int Path::randomWalk(const Scene *scene, Sampler *sampler,
		int nSteps, int rrStart, ETransportMode mode,
		MemoryPool &pool) {
	/* Determine the relevant edge and vertex to start the random walk */
	PathVertex *curVertex  = m_vertices[m_vertices.size()-1],
	           *predVertex = m_vertices.size() < 2 ? NULL :
	                         m_vertices[m_vertices.size()-2];
	PathEdge *predEdge     = m_edges.empty() ? NULL :
	                         m_edges[m_edges.size()-1];
	Spectrum throughput(1.0f);

	for (int i=0; i<nSteps || nSteps == -1; ++i) {
		PathVertex *succVertex = pool.allocVertex();
		PathEdge *succEdge = pool.allocEdge();

		if (!curVertex->sampleNext(scene, sampler, predVertex, predEdge, succEdge,
				succVertex, mode, rrStart != -1 && i >= rrStart, &throughput)) {
			pool.release(succVertex);
			pool.release(succEdge);
			return i;
		}

		append(succEdge, succVertex);

		predVertex = curVertex;
		curVertex = succVertex;
		predEdge = succEdge;
	}

	return nSteps;
}

int Path::randomWalkFromPixel(const Scene *scene, Sampler *sampler,
		int nSteps, const Point2i &pixelPosition, int rrStart, MemoryPool &pool) {

	PathVertex *v1 = pool.allocVertex(), *v2 = pool.allocVertex();
	PathEdge *e0 = pool.allocEdge(), *e1 = pool.allocEdge();

	/* Use a special sampling routine for the first two sensor vertices so that
	   the resulting subpath passes through the specified pixel position */
	int t = vertex(0)->sampleSensor(scene,
		sampler, pixelPosition, e0, v1, e1, v2);

	if (t < 1) {
		pool.release(e0);
		pool.release(v1);
		return 0;
	}

	append(e0, v1);

	if (t < 2) {
		pool.release(e1);
		pool.release(v2);
		return 1;
	}

	append(e1, v2);

	PathVertex *predVertex = v1, *curVertex = v2;
	PathEdge *predEdge = e1;
	Spectrum throughput(1.0f);

	for (; t<nSteps || nSteps == -1; ++t) {
		PathVertex *succVertex = pool.allocVertex();
		PathEdge *succEdge = pool.allocEdge();

		if (!curVertex->sampleNext(scene, sampler, predVertex, predEdge, succEdge,
				succVertex, ERadiance, rrStart != -1 && t >= rrStart, &throughput)) {
			pool.release(succVertex);
			pool.release(succEdge);
			return t;
		}

		append(succEdge, succVertex);

		predVertex = curVertex;
		curVertex = succVertex;
		predEdge = succEdge;
	}

	return nSteps;
}


std::pair<int, int> Path::alternatingRandomWalkFromPixel(const Scene *scene, Sampler *sampler,
		Path &emitterPath, int nEmitterSteps, Path &sensorPath, int nSensorSteps,
		const Point2i &pixelPosition, int rrStart, MemoryPool &pool) {
	/* Determine the relevant edges and vertices to start the random walk */
	PathVertex *curVertexS  = emitterPath.vertex(0),
	           *curVertexT  = sensorPath.vertex(0),
	           *predVertexS = NULL, *predVertexT = NULL;
	PathEdge   *predEdgeS  = NULL, *predEdgeT = NULL;

	PathVertex *v1 = pool.allocVertex(), *v2 = pool.allocVertex();
	PathEdge *e0 = pool.allocEdge(), *e1 = pool.allocEdge();

	/* Use a special sampling routine for the first two sensor vertices so that
	   the resulting subpath passes through the specified pixel position */
	int t = curVertexT->sampleSensor(scene,
		sampler, pixelPosition, e0, v1, e1, v2);

	if (t >= 1) {
		sensorPath.append(e0, v1);
	} else {
		pool.release(e0);
		pool.release(v1);
	}

	if (t == 2) {
		sensorPath.append(e1, v2);
		predVertexT = v1;
		curVertexT = v2;
		predEdgeT = e1;
	} else {
		pool.release(e1);
		pool.release(v2);
		curVertexT = NULL;
	}

	Spectrum throughputS(1.0f), throughputT(1.0f);

	int s = 0;
	do {
		if (curVertexT && (t < nSensorSteps || nSensorSteps == -1)) {
			PathVertex *succVertexT = pool.allocVertex();
			PathEdge *succEdgeT = pool.allocEdge();

			if (curVertexT->sampleNext(scene, sampler, predVertexT,
					predEdgeT, succEdgeT, succVertexT, ERadiance,
					rrStart != -1 && t >= rrStart, &throughputT)) {
				sensorPath.append(succEdgeT, succVertexT);
				predVertexT = curVertexT;
				curVertexT = succVertexT;
				predEdgeT = succEdgeT;
				t++;
			} else {
				pool.release(succVertexT);
				pool.release(succEdgeT);
				curVertexT = NULL;
			}
		} else {
			curVertexT = NULL;
		}

		if (curVertexS && (s < nEmitterSteps || nEmitterSteps == -1)) {
			PathVertex *succVertexS = pool.allocVertex();
			PathEdge *succEdgeS = pool.allocEdge();

			if (curVertexS->sampleNext(scene, sampler, predVertexS,
					predEdgeS, succEdgeS, succVertexS, EImportance,
					rrStart != -1 && s >= rrStart, &throughputS)) {
				emitterPath.append(succEdgeS, succVertexS);
				predVertexS = curVertexS;
				curVertexS = succVertexS;
				predEdgeS = succEdgeS;
				s++;
			} else {
				pool.release(succVertexS);
				pool.release(succEdgeS);
				curVertexS = NULL;
			}
		} else {
			curVertexS = NULL;
		}
	} while (curVertexS || curVertexT);

	return std::make_pair(s, t);
}

void Path::reverse() {
	std::reverse(m_vertices.begin(), m_vertices.end());
	std::reverse(m_edges.begin(), m_edges.end());
}

void Path::release(MemoryPool &pool) {
	for (size_t i=0; i<m_vertices.size(); ++i)
		pool.release(m_vertices[i]);
	for (size_t i=0; i<m_edges.size(); ++i)
		pool.release(m_edges[i]);
	m_vertices.clear();
	m_edges.clear();
}

void Path::clone(Path &target, MemoryPool &pool) const {
	target.release(pool);

	for (size_t i=0; i<m_vertices.size(); ++i)
		target.append(m_vertices[i]->clone(pool));
	for (size_t i=0; i<m_edges.size(); ++i)
		target.append(m_edges[i]->clone(pool));
}

void Path::append(const Path &path) {
	for (size_t i=0; i<path.vertexCount(); ++i)
		m_vertices.push_back(path.vertex(i));
	for (size_t i=0; i<path.edgeCount(); ++i)
		m_edges.push_back(path.edge(i));
}

void Path::append(const Path &path, size_t start, size_t end, bool reverse) {
	size_t count = end - start;

	for (size_t i=start; i<end; ++i) {
		m_vertices.push_back(path.vertex(i));
		if (i+1 < end)
			m_edges.push_back(path.edge(i));
	}

	if (reverse) {
		std::reverse(m_vertices.end() - count, m_vertices.end());
		std::reverse(m_edges.end() - (count - 1), m_edges.end());
	}
}

void Path::release(size_t start, size_t end, MemoryPool &pool) {
	for (size_t i=start; i<end; ++i) {
		pool.release(m_vertices[i]);
		if (i+1 < end)
			pool.release(m_edges[i]);
	}
}

bool Path::operator==(const Path &path) const {
	if (m_vertices.size() != path.vertexCount() ||
		m_edges.size() != path.edgeCount())
		return false;
	for (size_t i=0; i<m_vertices.size(); ++i)
		if (*path.vertex(i) != *m_vertices[i])
			return false;
	for (size_t i=0; i<m_edges.size(); ++i)
		if (*path.edge(i) != *m_edges[i])
			return false;
	return true;
}

Float Path::miWeight(const Scene *scene, const Path &emitterSubpath,
		const PathEdge *connectionEdge, const Path &sensorSubpath,
		int s, int t, bool sampleDirect, bool lightImage) {
	int k = s+t+1, n = k+1;

	const PathVertex
			*vsPred = emitterSubpath.vertexOrNull(s-1),
			*vtPred = sensorSubpath.vertexOrNull(t-1),
			*vs = emitterSubpath.vertex(s),
			*vt = sensorSubpath.vertex(t);

	/* pdfImp[i] and pdfRad[i] store the area/volume density of vertex
	   'i' when sampled from the adjacent vertex in the emitter
	   and sensor direction, respectively. */

	Float ratioEmitterDirect = 0.0f, ratioSensorDirect = 0.0f;
	Float *pdfImp      = (Float *) alloca(n * sizeof(Float)),
		  *pdfRad      = (Float *) alloca(n * sizeof(Float));
	bool  *connectable = (bool *)  alloca(n * sizeof(bool)),
		  *isNull      = (bool *)  alloca(n * sizeof(bool));

	/* Keep track of which vertices are connectable / null interactions */
	int pos = 0;
	for (int i=0; i<=s; ++i) {
		const PathVertex *v = emitterSubpath.vertex(i);
		connectable[pos] = v->isConnectable();
		isNull[pos] = v->isNullInteraction() && !connectable[pos];
		pos++;
	}

	for (int i=t; i>=0; --i) {
		const PathVertex *v = sensorSubpath.vertex(i);
		connectable[pos] = v->isConnectable();
		isNull[pos] = v->isNullInteraction() && !connectable[pos];
		pos++;
	}

	if (k <= 3)
		sampleDirect = false;

	EMeasure vsMeasure = EArea, vtMeasure = EArea;
	if (sampleDirect) {
		/* When direct sampling is enabled, we may be able to create certain
		   connections that otherwise would have failed (e.g. to an
		   orthographic camera or a directional light source) */
		const AbstractEmitter *emitter = (s > 0 ? emitterSubpath.vertex(1) : vt)->getAbstractEmitter();
		const AbstractEmitter *sensor = (t > 0 ? sensorSubpath.vertex(1) : vs)->getAbstractEmitter();

		EMeasure emitterDirectMeasure = emitter->getDirectMeasure();
		EMeasure sensorDirectMeasure  = sensor->getDirectMeasure();

		connectable[0]   = emitterDirectMeasure != EDiscrete && emitterDirectMeasure != EInvalidMeasure;
		connectable[1]   = emitterDirectMeasure != EInvalidMeasure;
		connectable[k-1] = sensorDirectMeasure != EInvalidMeasure;
		connectable[k]   = sensorDirectMeasure != EDiscrete && sensorDirectMeasure != EInvalidMeasure;

		/* The following is needed to handle orthographic cameras &
		   directional light sources together with direct sampling */
		if (t == 1)
			vtMeasure = sensor->needsDirectionSample() ? EArea : EDiscrete;
		else if (s == 1)
			vsMeasure = emitter->needsDirectionSample() ? EArea : EDiscrete;
	}

	/* Collect importance transfer area/volume densities from vertices */
	pos = 0;
	pdfImp[pos++] = 1.0;

	for (int i=0; i<s; ++i)
		pdfImp[pos++] = emitterSubpath.vertex(i)->pdf[EImportance]
			* emitterSubpath.edge(i)->pdf[EImportance];

	pdfImp[pos++] = vs->evalPdf(scene, vsPred, vt, EImportance, vsMeasure)
		* connectionEdge->pdf[EImportance];

	if (t > 0) {
		pdfImp[pos++] = vt->evalPdf(scene, vs, vtPred, EImportance, vtMeasure)
			* sensorSubpath.edge(t-1)->pdf[EImportance];

		for (int i=t-1; i>0; --i)
			pdfImp[pos++] = sensorSubpath.vertex(i)->pdf[EImportance]
				* sensorSubpath.edge(i-1)->pdf[EImportance];
	}

	/* Collect radiance transfer area/volume densities from vertices */
	pos = 0;
	if (s > 0) {
		for (int i=0; i<s-1; ++i)
			pdfRad[pos++] = emitterSubpath.vertex(i+1)->pdf[ERadiance]
				* emitterSubpath.edge(i)->pdf[ERadiance];

		pdfRad[pos++] = vs->evalPdf(scene, vt, vsPred, ERadiance, vsMeasure)
			* emitterSubpath.edge(s-1)->pdf[ERadiance];
	}

	pdfRad[pos++] = vt->evalPdf(scene, vtPred, vs, ERadiance, vtMeasure)
		* connectionEdge->pdf[ERadiance];

	for (int i=t; i>0; --i)
		pdfRad[pos++] = sensorSubpath.vertex(i-1)->pdf[ERadiance]
			* sensorSubpath.edge(i-1)->pdf[ERadiance];

	pdfRad[pos++] = 1.0;


	/* When the path contains specular surface interactions, it is possible
	   to compute the correct MI weights even without going through all the
	   trouble of computing the proper generalized geometric terms (described
	   in the SIGGRAPH 2012 specular manifolds paper). The reason is that these
	   all cancel out. But to make sure that that's actually true, we need to
	   convert some of the area densities in the 'pdfRad' and 'pdfImp' arrays
	   into the projected solid angle measure */
	for (int i=1; i <= k-3; ++i) {
		if (i == s || !(connectable[i] && !connectable[i+1]))
			continue;

		const PathVertex *cur = i <= s ? emitterSubpath.vertex(i) : sensorSubpath.vertex(k-i);
		const PathVertex *succ = i+1 <= s ? emitterSubpath.vertex(i+1) : sensorSubpath.vertex(k-i-1);
		const PathEdge *edge = i < s ? emitterSubpath.edge(i) : sensorSubpath.edge(k-i-1);

		pdfImp[i+1] *= edge->length * edge->length / std::abs(
			(succ->isOnSurface() ? dot(edge->d, succ->getGeometricNormal()) : 1) *
			(cur->isOnSurface()  ? dot(edge->d, cur->getGeometricNormal())  : 1));
	}

	for (int i=k-1; i >= 3; --i) {
		if (i-1 == s || !(connectable[i] && !connectable[i-1]))
			continue;

		const PathVertex *cur = i <= s ? emitterSubpath.vertex(i) : sensorSubpath.vertex(k-i);
		const PathVertex *succ = i-1 <= s ? emitterSubpath.vertex(i-1) : sensorSubpath.vertex(k-i+1);
		const PathEdge *edge = i <= s ? emitterSubpath.edge(i-1) : sensorSubpath.edge(k-i);

		pdfRad[i-1] *= edge->length * edge->length / std::abs(
			(succ->isOnSurface() ? dot(edge->d, succ->getGeometricNormal()) : 1) *
			(cur->isOnSurface()  ? dot(edge->d, cur->getGeometricNormal())  : 1));
	}

	int emitterRefIndirection = 2, sensorRefIndirection = k-2;

	/* One more array sweep before the actual useful work starts -- phew! :)
	   "Collapse" edges/vertices that were caused by BSDF::ENull interactions.
	   The BDPT implementation is smart enough to connect straight through those,
	   so they shouldn't be treated as Dirac delta events in what follows */
	for (int i=1; i <= k-3; ++i) {
		if (!connectable[i] || !isNull[i+1])
			continue;

		int start = i+1, end = start;
		while (isNull[end+1])
			++end;

		if (!connectable[end+1]) {
			/// The chain contains a non-ENull interaction
			isNull[start] = false;
			continue;
		}

		const PathVertex *before = i     <= s ? emitterSubpath.vertex(i) : sensorSubpath.vertex(k-i);
		const PathVertex *after  = end+1 <= s ? emitterSubpath.vertex(end+1) : sensorSubpath.vertex(k-end-1);

		Vector d = before->getPosition() - after->getPosition();
		Float lengthSquared = d.lengthSquared();
		d /= std::sqrt(lengthSquared);

		Float geoTerm = std::abs(
			(before->isOnSurface() ? dot(before->getGeometricNormal(), d) : 1) *
			(after->isOnSurface()  ? dot(after->getGeometricNormal(),  d) : 1)) / lengthSquared;

		pdfRad[start-1] *= pdfRad[end] * geoTerm;
		pdfRad[end] = 1;
		pdfImp[start] *= pdfImp[end+1] * geoTerm;
		pdfImp[end+1] = 1;

		/* When an ENull chain starts right after the emitter / before the sensor,
		   we must keep track of the reference vertex for direct sampling strategies. */
		if (start == 2)
			emitterRefIndirection = end + 1;
		else if (end == k-2)
			sensorRefIndirection = start - 1;

		i = end;
	}

	double initial = 1.0f;

	/* When direct sampling strategies are enabled, we must
	   account for them here as well */
	if (sampleDirect) {
		/* Direct connection probability of the emitter */
		const PathVertex *sample = s>0 ? emitterSubpath.vertex(1) : vt;
		const PathVertex *ref = emitterRefIndirection <= s
			? emitterSubpath.vertex(emitterRefIndirection) : sensorSubpath.vertex(k-emitterRefIndirection);
		EMeasure measure = sample->getAbstractEmitter()->getDirectMeasure();

		if (connectable[1] && connectable[emitterRefIndirection])
			ratioEmitterDirect = ref->evalPdfDirect(scene, sample, EImportance,
				measure == ESolidAngle ? EArea : measure) / pdfImp[1];

		/* Direct connection probability of the sensor */
		sample = t>0 ? sensorSubpath.vertex(1) : vs;
		ref = sensorRefIndirection <= s ? emitterSubpath.vertex(sensorRefIndirection)
			: sensorSubpath.vertex(k-sensorRefIndirection);
		measure = sample->getAbstractEmitter()->getDirectMeasure();

		if (connectable[k-1] && connectable[sensorRefIndirection])
			ratioSensorDirect = ref->evalPdfDirect(scene, sample, ERadiance,
				measure == ESolidAngle ? EArea : measure) / pdfRad[k-1];

		if (s == 1)
			initial /= ratioEmitterDirect;
		else if (t == 1)
			initial /= ratioSensorDirect;
	}

	double weight = 1, pdf = initial;

	/* With all of the above information, the MI weight can now be computed.
	   Since the goal is to evaluate the power heuristic, the absolute area
	   product density of each strategy is interestingly not required. Instead,
	   an incremental scheme can be used that only finds the densities relative
	   to the (s,t) strategy, which can be done using a linear sweep. For
	   details, refer to the Veach thesis, p.306. */
	for (int i=s+1; i<k; ++i) {
		double next = pdf * (double) pdfImp[i] / (double) pdfRad[i],
		       value = next;

		if (sampleDirect) {
			if (i == 1)
				value *= ratioEmitterDirect;
			else if (i == sensorRefIndirection)
				value *= ratioSensorDirect;
		}


		int tPrime = k-i-1;
		if (connectable[i] && (connectable[i+1] || isNull[i+1]) && (lightImage || tPrime > 1))
			weight += value*value;

		pdf = next;
	}

	/* As above, but now compute pdf[i] with i<s (this is done by
	   evaluating the inverse of the previous expressions). */
	pdf = initial;
	for (int i=s-1; i>=0; --i) {
		double next = pdf * (double) pdfRad[i+1] / (double) pdfImp[i+1],
		       value = next;

		if (sampleDirect) {
			if (i == 1)
				value *= ratioEmitterDirect;
			else if (i == sensorRefIndirection)
				value *= ratioSensorDirect;
		}

		int tPrime = k-i-1;
		if (connectable[i] && (connectable[i+1] || isNull[i+1]) && (lightImage || tPrime > 1))
			weight += value*value;

		pdf = next;
	}

	return (Float) (1.0 / weight);
}

void Path::collapseTo(PathEdge &target) const {
	BDAssert(m_edges.size() > 0);

	target.pdf[ERadiance] = 1.0f;
	target.pdf[EImportance] = 1.0f;
	target.weight[ERadiance] = Spectrum(1.0f);
	target.weight[EImportance] = Spectrum(1.0f);
	target.d = m_edges[0]->d;
	target.medium = m_edges[0]->medium;
	target.length = 0;

	for (size_t i=0; i<m_edges.size(); ++i) {
		const PathEdge *edge = m_edges[i];
		target.weight[ERadiance] *= edge->weight[ERadiance];
		target.weight[EImportance] *= edge->weight[EImportance];
		target.pdf[ERadiance] *= edge->pdf[ERadiance];
		target.pdf[EImportance] *= edge->pdf[EImportance];
		target.length += edge->length;
		if (target.medium != edge->medium)
			target.medium = NULL;
	}

	for (size_t i=0; i<m_vertices.size(); ++i) {
		const PathVertex *vertex = m_vertices[i];
		BDAssert(vertex->isSurfaceInteraction() &&
			vertex->componentType == BSDF::ENull);
		target.weight[ERadiance] *= vertex->weight[ERadiance];
		target.weight[EImportance] *= vertex->weight[EImportance];
		target.pdf[ERadiance] *= vertex->pdf[ERadiance];
		target.pdf[EImportance] *= vertex->pdf[EImportance];
	}
}

std::string Path::toString() const {
	std::ostringstream oss;
	oss << "Path[" << endl;

	if (m_vertices.size() == m_edges.size() + 1) {
		for (size_t i=0; i<m_vertices.size(); ++i) {
			oss << "  Vertex " << i << " => " << indent(m_vertices[i]->toString());

			if (i<m_edges.size()) {
				oss << "," << endl;
				oss << "  Edge " << i << " => " << indent(m_edges[i]->toString());
			}

			if (i+1 < m_vertices.size())
				oss << ",";
			oss << endl;
		}
	} else if (m_edges.size() == m_vertices.size() + 1) {
		for (size_t i=0; i<m_edges.size(); ++i) {
			oss << "  Edge " << i << " => " << indent(m_edges[i]->toString());

			if (i<m_vertices.size()) {
				oss << "," << endl;
				oss << "  Vertex " << i << " => " << indent(m_vertices[i]->toString());
			}

			if (i+1 < m_edges.size())
				oss << ",";
			oss << endl;
		}
	} else {
		SLog(EError, "Unknown path configuration!");
	}

	oss << "]";
	return oss.str();
}

std::string Path::summarize() const {
	std::ostringstream oss;

	BDAssert(m_vertices.size() == m_edges.size() + 1);

	for (size_t i=0; i<m_vertices.size(); ++i) {
		const PathVertex *vertex = m_vertices[i];
		switch (vertex->type) {
			case PathVertex::EEmitterSupernode: oss << "E"; break;
			case PathVertex::ESensorSupernode: oss << "C"; break;
			case PathVertex::EEmitterSample: oss << "e"; break;
			case PathVertex::ESensorSample: oss << "c"; break;
			case PathVertex::ESurfaceInteraction: oss << "S"; break;
			case PathVertex::EMediumInteraction: oss << "M"; break;
			default:
				SLog(EError, "Unknown vertex typ!");
		}
		if (!vertex->isConnectable())
			oss << "d";

		if (i<m_edges.size()) {
			const PathEdge *edge = m_edges[i];
			if (edge->medium == NULL)
				oss << " - ";
			else
				oss << " = ";
		}
	}
	return oss.str();
}

MTS_NAMESPACE_END
