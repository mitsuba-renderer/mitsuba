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

#if !defined(__RENDER_COMMON_H)
#define __RENDER_COMMON_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Specifies the transported quantity when 
 * sampling or evaluating a scattering function
 * \ingroup librender
 */
enum ETransportMode {
	/* Note to self: do not change these enumeration
	   values, some code depends on them. */

	/// Radiance transport
	ERadiance = 0,
	/// Importance transport
	EImportance = 1,
	/// Specifies the number of supported transport modes
	ETransportModes = 2
};
/**
 * \brief Specifies the measure associated with 
 * a scattering function 
 * \ingroup librender
 */
enum EMeasure {
	ESolidAngle = 1,
	EInterval,
	EDiscrete
};


/// \cond
extern MTS_EXPORT_RENDER std::ostream &operator<<(std::ostream &os, const ETransportMode &quantity);
extern MTS_EXPORT_RENDER std::ostream &operator<<(std::ostream &os, const EMeasure &measure);
/// \endcond

MTS_NAMESPACE_END

#endif /* __RENDER_COMMON_H */
