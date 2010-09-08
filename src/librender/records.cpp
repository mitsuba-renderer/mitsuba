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

#include <mitsuba/render/records.h>
	
MTS_NAMESPACE_BEGIN

void Intersection::computePartials(const RayDifferential &ray) {
	/* Compute the texture coordinates partials wrt. 
		changes in the screen-space position. Based on PBRT */
	if (hasUVPartials)
		return;
	hasUVPartials = true;

	if (!ray.hasDifferentials) {
		dudx = dvdx = dudy = dvdy = 0.0f;
		return;
	}

	/* Offset of the plane passing through the surface */
	const Float d = -dot(geoFrame.n, p);

	const Float txRecip = dot(geoFrame.n, ray.rx.d),
				tyRecip = dot(geoFrame.n, ray.ry.d);

	if (EXPECT_NOT_TAKEN(txRecip == 0 || tyRecip == 0)) {
		dudx = dvdx = dudy = dvdy = 0.0f;
		return;
	}

	/* Ray distances traveled */
	const Float tx = -(dot(geoFrame.n, ray.rx.o) + d) / 
		txRecip;
	const Float ty = -(dot(geoFrame.n, ray.ry.o) + d) / 
		tyRecip;

	/* Auxilary intersection point of the adjacent rays */
	Point px = ray.rx(tx), py = ray.ry(ty);

	/* Calculate the U and V partials by solving two out
		of a set of 3 equations in an overconstrained system */
	Float A[2][2], Bx[2], By[2], x[2];
	int axes[2];

	if (std::abs(geoFrame.n.x) > std::abs(geoFrame.n.y)
		&& std::abs(geoFrame.n.x) > std::abs(geoFrame.n.z)) {
		axes[0] = 1; axes[1] = 2;
	} else if (std::abs(geoFrame.n.y) > std::abs(geoFrame.n.z)) {
		axes[0] = 0; axes[1] = 2;
	} else {
		axes[0] = 0; axes[1] = 1;
	}

	A[0][0] = dpdu[axes[0]];
	A[0][1] = dpdv[axes[0]];
	A[1][0] = dpdu[axes[1]];
	A[1][1] = dpdv[axes[1]];
	Bx[0] = px[axes[0]] - p[axes[0]];
	Bx[1] = px[axes[1]] - p[axes[1]];
	By[0] = py[axes[0]] - p[axes[0]];
	By[1] = py[axes[1]] - p[axes[1]];

	if (solveLinearSystem2x2(A, Bx, x)) {
		dudx = x[0]; dvdx = x[1];
	} else {
		dudx = 1; dvdx = 0;
	}

	if (solveLinearSystem2x2(A, By, x)) {
		dudy = x[0]; dvdy = x[1];
	} else {
		dudy = 0; dudy = 1;
	}
}

std::string Intersection::toString() const {
	if (!isValid())
		return "Intersection[invalid]";
	std::ostringstream oss;
	oss << "Intersection[" << std::endl
		<< "  p = " << p.toString() << "," << std::endl
		<< "  wi = " << wi.toString() << "," << std::endl
		<< "  t = " << t << "," << std::endl
		<< "  geoFrame = " << indent(geoFrame.toString()) << "," << std::endl
		<< "  shFrame = " << indent(shFrame.toString()) << "," << std::endl
		<< "  uv = " << uv.toString() << "," << std::endl
		<< "  dpdu = " << dpdu.toString() << "," << std::endl
		<< "  dpdv = " << dpdv.toString() << "," << std::endl
		<< "  shape = " << indent(((Object *)shape)->toString()) << std::endl
		<< "]";
	return oss.str();
}


/// Return a string representation
std::string RadianceQueryRecord::toString() const {
	std::ostringstream oss;
	oss << "RadianceQueryRecord[" << std::endl
		<< "  type = { ";
	if (type & EEmittedRadiance) oss << "emitted ";
	if (type & ESubsurfaceRadiance) oss << "subsurface ";
	if (type & EDirectRadiance) oss << "direct ";
	if (type & EIndirectRadiance) oss << "indirect ";
	if (type & ECausticRadiance) oss << "caustic ";
	if (type & EInscatteredDirectRadiance) oss << "inscatteredDirect ";
	if (type & EInscatteredIndirectRadiance) oss << "inscatteredIndirect ";
	if (type & EDistance) oss << "distance ";
	if (type & EOpacity) oss << "opacity ";
	if (type & EIntersection) oss << "intersection ";
	oss << "}," << std::endl
		<< "  depth = " << depth << "," << std::endl
		<< "  its = " << indent(its.toString()) << std::endl
		<< "  alpha = " << alpha << "," << std::endl
		<< "  extra = " << extra << "," << std::endl
		<< "]" << std::endl;
	return oss.str();
}

std::string BSDFQueryRecord::toString() const {
	std::ostringstream oss;
	oss << "BSDFQueryRecord[" << std::endl
		<< "  wi = " << wi.toString() << "," << std::endl
		<< "  wo = " << wo.toString() << "," << std::endl
		<< "  sample = " << sample.toString() << "," << std::endl
		<< "  quantity = " << quantity << "," << std::endl
		<< "  typeMask = " << typeMask << "," << std::endl
		<< "  sampledType = " << sampledType << "," << std::endl
		<< "  component = " << component << "," << std::endl
		<< "  sampledComponent = " << sampledComponent << "," << std::endl
		<< "]";
	return oss.str();
}
	
std::string ShapeSamplingRecord::toString() const {
	std::ostringstream oss;
	oss << "ShapeSamplingRecord[" << std::endl
		<< "  p = " << p.toString() << "," << std::endl
		<< "  n = " << n.toString() << std::endl
		<< "]";
	return oss.str();
}

std::string MediumSamplingRecord::toString() const {
	std::ostringstream oss;
	oss << "MediumSamplingRecord[" << std::endl
		<< "  t = " << t << "," << std::endl
		<< "  p = " << p.toString() << "," << std::endl
		<< "  sigmaA = " << sigmaA.toString() << "," << std::endl
		<< "  sigmaS = " << sigmaS.toString() << "," << std::endl
		<< "  pdf = " << pdf << "," << std::endl
		<< "  medium = " << indent(((Object *) medium)->toString()) << std::endl
		<< "]";
	return oss.str();
}

std::string LuminaireSamplingRecord::toString() const {
	std::ostringstream oss;
	oss << "LuminaireSamplingRecord[" << std::endl
		<< "  sRec = " << indent(sRec.toString()) << "," << std::endl
		<< "  d = " << d.toString() << "," << std::endl
		<< "  pdf = " << pdf << "," << std::endl
		<< "  Le = " << Le.toString() << "," << std::endl
		<< "  luminaire = " << ((luminaire == NULL ) ? "null" : indent(((Object *) luminaire)->toString()).c_str()) << std::endl
		<< "]";
	return oss.str();
}

std::string EmissionRecord::toString() const {
	std::ostringstream oss;
	oss << "EmissionRecord[" << std::endl
		<< "  sRec = " << indent(sRec.toString()) << "," << std::endl
		<< "  d = " << d.toString() << "," << std::endl
		<< "  P = " << P.toString() << "," << std::endl
		<< "  pdfArea = " << pdfArea << "," << std::endl
		<< "  pdfDir = " << pdfDir << "," << std::endl
		<< "  luminaire = " << ((luminaire == NULL ) ? "null" : indent(((Object *) luminaire)->toString()).c_str()) << std::endl
		<< "]";
	return oss.str();
}
	
void operator<<(const ETransportQuantity &quantity, std::ostream &os) {
	switch (quantity) {
		case EImportance: os << "importance"; break;
		case ERadiance:   os << "radiance"; break;
		default: os << "invalid"; break;
	};
}

MTS_NAMESPACE_END
