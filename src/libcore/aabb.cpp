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

#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

Point AABB::getCorner(uint8_t corner) const {
	return Point(corner & 1 ? max.x : min.x,
			corner & 2 ? max.y : min.y,
			corner & 4 ? max.z : min.z);
}


bool AABB::overlaps(const BSphere &sphere) const {
	Float distance = 0;
	for (int i=0; i<3; ++i) {
		if (sphere.center[i] < min[i]) {
			Float d = sphere.center[i]-min[i];
			distance += d*d;
		} else if (sphere.center[i] > max[i]) {
			Float d = sphere.center[i]-max[i];
			distance += d*d;
		}
	}
	return distance < sphere.radius*sphere.radius;
}

BSphere AABB::getBSphere() const {
	Point3 center = getCenter();
	return BSphere(center, (center - max).length());
}

MTS_NAMESPACE_END
