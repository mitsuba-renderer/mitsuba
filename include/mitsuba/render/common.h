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

#if !defined(__RENDER_COMMON_H)
#define __RENDER_COMMON_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * Specifies the transported quantity when sampling / evaluating a BSDF
 */
enum ETransportQuantity {
	ERadiance = 1,
	EImportance = 2
};

extern MTS_EXPORT_RENDER std::ostream &operator<<(std::ostream &os, const ETransportQuantity &quantity);

MTS_NAMESPACE_END

#endif /* __RENDER_COMMON_H */
