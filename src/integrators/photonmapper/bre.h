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

#if !defined(__BEAM_RADIANCE_ESTIMATE_H)
#define __BEAM_RADIANCE_ESTIMATE_H

#include <mitsuba/render/photonmap.h>

MTS_NAMESPACE_BEGIN

/**
 * Implements the beam radiance estimate described in
 * "The Beam Radiance Estimate for Volumetric Photon Mapping"
 * by Wojciech Jarosz, Matthias Zwicker, and Henrik Wann Jensen.
 */

class BeamRadianceEstimate {
public:
	/**
	 * \brief Create a BRE acceleration data structure from
	 * an existing volumetric photon map
	 */
	BeamRadianceEstimate(PhotonMap *pmap);

	MTS_DECLARE_CLASS()
private:
};

MTS_NAMESPACE_END

#endif /* __BEAM_RADIANCE_ESTIMATE_H */
