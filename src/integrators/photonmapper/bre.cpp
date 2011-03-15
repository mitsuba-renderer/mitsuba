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

#include "bre.h"

MTS_NAMESPACE_BEGIN

BeamRadianceEstimate::BeamRadianceEstimate(PhotonMap *pmap) {
	int n = 100;
	
	PhotonMap::search_result *results = 
		new PhotonMap::search_result[n+1];

	Float *radii = new Float[pmap->getPhotonCount()];

	Log(EInfo, "Computing photon radii ..");
	for (size_t i=0; i<pmap->getPhotonCount(); ++i) {
		const Photon &photon = pmap->getPhoton(i);

		Float searchRadiusSqr = std::numeric_limits<Float>::infinity();
		pmap->nnSearch(photon.getPosition(), searchRadiusSqr, n, results);
		radii[i] = std::sqrt(searchRadiusSqr);
	}

	delete[] radii;
	delete[] results;
}

MTS_IMPLEMENT_CLASS(BeamRadianceEstimate, false, Object)
MTS_NAMESPACE_END
