/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#include <mitsuba/render/common.h>

MTS_NAMESPACE_BEGIN

std::ostream &operator<<(std::ostream &os, const ETransportMode &quantity) {
	switch (quantity) {
		case EImportance: os << "importance"; break;
		case ERadiance:   os << "radiance"; break;
		default:          os << "invalid"; break;
	};
	return os;
}

std::ostream &operator<<(std::ostream &os, const EMeasure &measure) {
	switch (measure) {
		case ESolidAngle: os << "solidAngle"; break;
		case EInterval:   os << "interval"; break;
		case EDiscrete:   os << "discrete"; break;
		default:          os << "invalid"; break;
	};
	return os;
}


MTS_NAMESPACE_END
