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

#include <mitsuba/core/aabb.h>

MTS_NAMESPACE_BEGIN

Point AABB::getCorner(uint8_t corner) const {
	return Point(corner & 1 ? max.x : min.x,
			corner & 2 ? max.y : min.y,
			corner & 4 ? max.z : min.z);
}

Point AABB::getMidPoint() const {
	return Point(
		min.x + max.x,
		min.y + max.y,
		min.z + max.z
	) * 0.5f;
}

bool AABB::contains(const Point &vec) const {
	return isValid() &&
		(vec.x >= min.x && vec.x <= max.x) &&
		(vec.y >= min.y && vec.y <= max.y) &&
		(vec.z >= min.z && vec.z <= max.z);
}

void AABB::expandBy(const Point &vec) {
	min.x = std::min(min.x, vec.x);
	min.y = std::min(min.y, vec.y);
	min.z = std::min(min.z, vec.z);
	max.x = std::max(max.x, vec.x);
	max.y = std::max(max.y, vec.y);
	max.z = std::max(max.z, vec.z);
}

void AABB::expandBy(const AABB &aabb) {
	min.x = std::min(min.x, aabb.min.x);
	min.y = std::min(min.y, aabb.min.y);
	min.z = std::min(min.z, aabb.min.z);
	max.x = std::max(max.x, aabb.max.x);
	max.y = std::max(max.y, aabb.max.y);
	max.z = std::max(max.z, aabb.max.z);
}

Float AABB::distanceTo(const Point &p) const {
	Float result = 0;
	for (int i=0; i<3; ++i) {
		Float value = 0;
		if (p[i] < min[i])
			value = min[i] - p[i];
		if (p[i] > max[i])
			value = p[i] - max[i];
		result += value*value;
	}
	return std::sqrt(result);
}

std::string AABB::toString() const {
	std::ostringstream oss;
	oss << "AABB[";
	if (!isValid()) {
		oss << "invalid";
	} else {
		oss << "min=" << min.toString()
			<< ", max=" << max.toString();
	}
	oss	<< "]";
	return oss.str();
}

MTS_NAMESPACE_END
