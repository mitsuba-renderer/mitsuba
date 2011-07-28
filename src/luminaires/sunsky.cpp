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

#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>

#define SAMPLE_UNIFORMLY 1

MTS_NAMESPACE_BEGIN

/*!\plugin{sky}{Sun and sky luminaire}
 * \parameters{
 *     \parameter{turbidity}{\Float}{
 *         This parameter determines the amount of scattering particles (or
 *         `haze') in the atmosphere. Smaller values ($\sim 2$) produce a 
 *         clear blue sky, larger values ($\sim 8$) lead to an overcast sky, 
 *         and a very high values ($\sim 20$) cause a color shift towards 
 *         orange and red. \default{3}
 *     }
 *     \parameter{day}{\Integer}{Solar day used to compute the sun's position. 
 *       Must be in the range between 1 and 365. \default{180}}
 *     \parameter{time}{\Float}{Fractional time used to compute the sun's
 *       position. A time of 4:15 PM corresponds to 16.25. \default{15.00}}
 *     \parameter{latitude, longitude}{\Float}{
 *       These two parameters specify the oberver's latitude and longitude 
 *       in degrees, which are required to compute the sun's position.
 *       \default{35.6894, 139.6917 --- Tokyo, Japan}
 *     }
 *     \parameter{standardMeridian}{\Integer}{Denotes the
 *       standard meridian of the time zone for finding
 *       the sun's position \default{135 --- Japan standard time}
 *     }
 *     \parameter{sunDirection}{\Vector}{Allows to manually 
 *       override the sun direction in world space. When this value
 *       is provided, parameters pertaining to the computation 
 *       of the sun direction (\code{day, time, latitude, longitude,} 
 *       and \code{standardMeridian}) are unnecessary. \default{none}
 *     }
 *     \parameter{extend}{\Boolean}{
 *         Extend luminaire below the horizon? \default{\code{false}}
 *     }
 *     \parameter{resolution}{\Integer}{Specifies the resolution of the precomputed
 *         image that is used to represent the sky environment map
 *         \default{256}}
 *     \parameter{scale}{\Float}{
 *         This parameter can be used to scale the the amount of illumination
 *         emitted by the sky luminaire, for instance to change its units. To
 *         switch from photometric ($\nicefrac{W}{m^2\cdot sr}$) 
 *         to arbitrary but convenient units in the $[0,1]$ range, set 
 *         this parameter to \code{1e-5}.\default{1}.
 *     }
 * }
 *
 * This plugin implements the physically-based skylight model proposed by 
 * Preetham et al. \cite{Preetham1999Practical}. It can be used for realistic 
 * daylight renderings of scenes under clear and overcast skies, assuming
 * that the sky is observed from a position either on or close to the surface 
 * of the earth. 
 *
 * This is a convenience plugin, which has the sole purpose of instantiating 
 * \pluginref{sun} and \pluginref{sky} at the same time. Please refer to these
 * plugins individually for more detail
 */
class SunSkyLuminaire : public Luminaire {
public:
	SunSkyLuminaire(const Properties &_props)
		: Luminaire(_props) {
		Properties props(_props);
		props.setPluginName("sky");
		m_sky = static_cast<Luminaire *>(
			PluginManager::getInstance()->createObject(
			MTS_CLASS(Luminaire), props));

		/* Avoid unused parameter warnings */
		std::vector<std::string> propNames;
		props.putPropertyNames(propNames);
		for (size_t i=0; i<propNames.size(); ++i)
			if (props.wasQueried(propNames[i]))
				_props.markQueried(propNames[i]);
	}

	SunSkyLuminaire(Stream *stream, InstanceManager *manager)
		: Luminaire(stream, manager) {
		m_sky = static_cast<Luminaire *>(manager->getInstance(stream));
	}

	void serialize(Stream *stream, InstanceManager *manager) {
		Luminaire::serialize(stream, manager);
		manager->serialize(stream, m_sky.get());
	}

	void configure() {
		Luminaire::configure();
		m_sky->configure();
	}

	bool isCompound() const {
		return true;
	}

	Luminaire *getElement(int i) {
		if (i == 0)
			return m_sky;
		else
			return NULL;
	}

	MTS_DECLARE_CLASS()
private:
	Properties m_props;
	ref<Luminaire> m_sky;
};

MTS_IMPLEMENT_CLASS_S(SunSkyLuminaire, false, Luminaire)
MTS_EXPORT_PLUGIN(SunSkyLuminaire, "Preetham sun sky luminaire");
MTS_NAMESPACE_END

