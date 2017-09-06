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

#include <mitsuba/render/volume.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{constvolume}{Constant-valued volume data source}
 * \parameters{
 *     \parameter{value}{\Float\Or\Spectrum\Or\Vector}{
 *       Specifies the value of the volume
 *     }
 * }
 *
 * This plugin provides a volume data source that is
 * constant throughout its domain. Depending on how it is used,
 * its value can either be a scalar, a color spectrum,
 * or a 3D vector.\vspace{4mm}
 *
 * \begin{xml}[caption={Definition of a heterogeneous medium with homogeneous contents}]
 * <medium type="heterogeneous">
 *     <volume type="constvolume" name="density">
 *         <float name="value" value="1"/>
 *     </volume>
 *     <volume type="constvolume" name="albedo">
 *         <rgb name="value" value="0.9 0.9 0.7"/>
 *     </volume>
 *     <volume type="constvolume" name="orientation">
 *         <vector name="value" x="0" y="1" z="0"/>
 *     </volume>
 *
 *     <!-- .... remaining parameters for
 *          the 'heterogeneous' plugin .... -->
 * </medium>
 * \end{xml}
 */
class ConstantDataSource : public VolumeDataSource {
public:
    ConstantDataSource(const Properties &props)
        : VolumeDataSource(props) {
        m_type = props.getType("value");

        if (m_type == Properties::EFloat)
            m_float = props.getFloat("value");
        else if (m_type == Properties::EVector)
            m_vector = props.getVector("value");
        else if (m_type == Properties::ESpectrum)
            m_spectrum = props.getSpectrum("value");
        else
            Log(EError, "The value of a 'constvolume' must have "
                "one of the following types: float, vector, spectrum");
    }

    ConstantDataSource(Stream *stream, InstanceManager *manager)
        : VolumeDataSource(stream, manager) {
        m_type = stream->readInt();
        if (m_type == Properties::EFloat)
            m_float = stream->readFloat();
        else if (m_type == Properties::EVector)
            m_vector = Vector(stream);
        else if (m_type == Properties::ESpectrum)
            m_spectrum = Spectrum(stream);
        else
            Log(EError, "Internal error - unknown data type");
    }

    virtual ~ConstantDataSource() {
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        VolumeDataSource::serialize(stream, manager);
        stream->writeInt(m_type);
        if (m_type == Properties::ESpectrum)
            m_spectrum.serialize(stream);
        else if (m_type == Properties::EFloat)
            stream->writeFloat(m_float);
        else if (m_type == Properties::EVector)
            m_vector.serialize(stream);
        else
            Log(EError, "Internal error - unknown data type");
    }

    Float lookupFloat(const Point &p) const {
        return m_float;
    }

    Spectrum lookupSpectrum(const Point &p) const {
        return m_spectrum;
    }

    Vector lookupVector(const Point &p) const {
        return m_vector;
    }

    bool supportsFloatLookups() const {
        return m_type == Properties::EFloat;
    }

    bool supportsSpectrumLookups() const {
        return m_type == Properties::ESpectrum;
    }

    bool supportsVectorLookups() const {
        return m_type == Properties::EVector;
    }

    Float getStepSize() const {
        return std::numeric_limits<Float>::infinity();
    }

    Float getMaximumFloatValue() const {
        return m_float;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "ConstantDataSource[value=";
        if (m_type == Properties::EFloat)
            oss << m_float;
        else if (m_type == Properties::EVector)
            oss << m_vector.toString();
        else if (m_type == Properties::ESpectrum)
            oss << m_spectrum.toString();
        else
            Log(EError, "Invalid volume data type!");
        oss << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
protected:
    int m_type;
    Float m_float;
    Vector m_vector;
    Spectrum m_spectrum;
};

MTS_IMPLEMENT_CLASS_S(ConstantDataSource, false, VolumeDataSource);
MTS_EXPORT_PLUGIN(ConstantDataSource, "Constant data source");
MTS_NAMESPACE_END
