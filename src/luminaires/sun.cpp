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
#include <mitsuba/core/bitmap.h>
#include <mitsuba/core/plugin.h>
#include "sun.h"

MTS_NAMESPACE_BEGIN

/*!\plugin{sun}{Sun luminaire}
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
 *     \parameter{resolution}{\Integer}{Specifies the resolution of the precomputed
 *         image that is used to represent the sun environment map
 *         \default{256}}
 *     \parameter{sunScale}{\Float}{
 *         This parameter can be used to scale the the amount of illumination
 *         emitted by the sun luminaire, for instance to change its units. To
 *         switch from photometric ($\nicefrac{W}{m^2\cdot sr}$) 
 *         to arbitrary but convenient units in the $[0,1]$ range, set 
 *         this parameter to \code{1e-5}. \default{1}
 *     }
 * }
 */
class SunLuminaire : public Luminaire {
public:
	SunLuminaire(const Properties &props)
			: Luminaire(props) {
		m_scale = props.getFloat("sunScale", 1.0f);
		m_resolution = props.getInteger("resolution", 512);
		Point2 sunPos = configureSunPosition(props);
		m_sunTheta = sunPos.x;
		m_sunDir = toSphere(sunPos.x, sunPos.y);
		m_sunDiskScale = props.getFloat("sunDiskScale", 15.0f);
		m_turbidity = props.getFloat("turbidity", 3.0f);
	}

	SunLuminaire(Stream *stream, InstanceManager *manager) 
		    : Luminaire(stream, manager) {
		m_scale = stream->readFloat();
		m_sunDiskScale = stream->readFloat();
		m_sunTheta = stream->readFloat();
		m_turbidity = stream->readFloat();
		m_resolution = stream->readInt();
		m_sunDir = Vector(stream);
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Luminaire::serialize(stream, manager);
		stream->writeFloat(m_scale);
		stream->writeFloat(m_sunDiskScale);
		stream->writeFloat(m_sunTheta);
		stream->writeFloat(m_turbidity);
		stream->writeInt(m_resolution);
		m_sunDir.serialize(stream);
	}

	void configure() {
		m_sunDiskTheta = 0.0046333f;
		Float sunDiskSA = 2*M_PI * (1-std::cos(m_sunDiskTheta));
		m_sunDiskTheta *= m_sunDiskScale;
		Float sunDiskSAScaled = 2*M_PI * (1-std::cos(m_sunDiskTheta));

		m_intensity = computeSunRadiance(m_sunTheta, m_turbidity) 
			* m_scale * sunDiskSA/sunDiskSAScaled;
	}

	bool isCompound() const {
		return true;
	}

	Luminaire *getElement(int i) {
		if (i != 0)
			return NULL;
		int thetaBins = m_resolution, phiBins = m_resolution*2;

		ref<Bitmap> bitmap = new Bitmap(phiBins, thetaBins, 128);
		Point2 factor(M_PI / thetaBins, (2*M_PI) / phiBins);
		float *target = bitmap->getFloatData();
		for (int i=0; i<thetaBins; ++i) {
			Float theta = (i+.5f)*factor.x;
			for (int j=0; j<phiBins; ++j) {
				Float phi = (j+.5f)*factor.y;
				Spectrum s = Le(Ray(Point(0.0f), toSphere(theta, phi), 0))
						* m_scale;
				Float r, g, b;
				s.toLinearRGB(r, g, b);
				*target++ = r; *target++ = g;
				*target++ = b; *target++ = 1;
			}
		}

		/* Instantiate a nested envmap plugin */
		Properties props("envmap");
		Properties::Data bitmapData;
		bitmapData.ptr = (uint8_t *) bitmap.get();
		bitmapData.size = sizeof(Bitmap);
		props.setData("bitmap", bitmapData);
		props.setTransform("toWorld", m_luminaireToWorld);
		props.setFloat("samplingWeight", m_samplingWeight);
		Luminaire *luminaire = static_cast<Luminaire *>(
			PluginManager::getInstance()->createObject(
			MTS_CLASS(Luminaire), props));
		luminaire->configure();
		return luminaire;
	}

	Spectrum Le(const Ray &ray) const {
		Float theta = unitAngle(ray.d, m_sunDir);

		if (theta < m_sunDiskTheta) 
			return m_intensity * smoothStep(m_sunDiskTheta, 0, theta);
		else
			return Spectrum(0.0f);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SunLuminaire[" << endl
			<< "  sunDir = " << m_sunDir.toString() << "," << endl 
			<< "  sunDiskScale = " << m_sunDiskScale << "," << endl 
			<< "  sunScale = " << m_scale << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	/* Environment map resolution */
	int m_resolution;
	/* Constant scale factor applied to the model */
	Float m_scale;
	/* Direction of the sun in luminaire-space */
	Vector m_sunDir;
	/* Angle cutoff for the sun disk */
	Float m_sunDiskTheta;
	/* Sun disk region scale factors */
	Float m_sunDiskScale;
	/* The turbidity of the sky ranges normally from 1 to 30.
	   For clear skies values in range [2,6] are useful. */
	Float m_turbidity;
	/* Intensity of the sun disk */
	Spectrum m_intensity;
	/* Elevation in radians */
	Float m_sunTheta;
};

MTS_IMPLEMENT_CLASS_S(SunLuminaire, false, Luminaire)
MTS_EXPORT_PLUGIN(SunLuminaire, "Preetham sky luminaire");
MTS_NAMESPACE_END

