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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/renderproc.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{field}{Field extraction integrator}
 * \order{17}
 * \parameters{
 *     \parameter{field}{\String}{Denotes the name of the field that should be extracted.
 *        The following choices are possible:
 *        \begin{itemize}
 *            \setlength{\itemsep}{1pt}
 *            \setlength{\parskip}{1pt}
 *            \item \code{position}: 3D position in world space
 *            \item \code{relPosition}: 3D position in camera space
 *            \item \code{distance}: Ray distance to the shading point
 *            \item \code{geoNormal}: Geometric surface normal
 *            \item \code{shNormal}: Shading surface normal
 *            \item \code{uv}: UV coordinate value
 *            \item \code{albedo}: Albedo value of the BSDF
 *            \item \code{shapeIndex}: Integer index of the high-level shape
 *            \item \code{primIndex}: Integer shape primitive index
 *        \end{itemize}
 *     }
 *     \parameter{undefined}{\Spectrum\Or\Float}{Value that should be returned when
 *                           there is no intersection \default{0}}
 * }
 *
 * This integrator extracts a requested field of from the intersection records of shading
 * points and converts the resulting data into color values. It is meant to be used in conjunction with
 * \pluginref{multichannel} to dump auxiliary information (such as depth or surface normals
 * of surfaces seen by the camera) into extra channels of a rendered image, for instance to
 * create benchmark data for computer vision applications.
 * Please refer to the documentation of \pluginref{multichannel} for an example.
 */

class FieldIntegrator : public SamplingIntegrator {
public:
    enum EField {
        EPosition,
        ERelativePosition,
        EDistance,
        EGeometricNormal,
        EShadingNormal,
        EUV,
        EAlbedo,
        EShapeIndex,
        EPrimIndex
    };

    FieldIntegrator(const Properties &props) : SamplingIntegrator(props) {
        std::string field = props.getString("field");

        if (field == "position") {
            m_field = EPosition;
        } else if (field == "relPosition") {
            m_field = ERelativePosition;
        } else if (field == "distance") {
            m_field = EDistance;
        } else if (field == "geoNormal") {
            m_field = EGeometricNormal;
        } else if (field == "shNormal") {
            m_field = EShadingNormal;
        } else if (field == "uv") {
            m_field = EUV;
        } else if (field == "albedo") {
            m_field = EAlbedo;
        } else if (field == "shapeIndex") {
            m_field = EShapeIndex;
        } else if (field == "primIndex") {
            m_field = EPrimIndex;
        } else {
            Log(EError, "Invalid 'field' parameter. Must be one of 'position', "
                "'relPosition', 'distance', 'geoNormal', 'shNormal', "
                "'primIndex', 'shapeIndex', or 'uv'!");
        }

        if (props.hasProperty("undefined")) {
            if (props.getType("undefined") == Properties::EFloat)
                m_undefined = Spectrum(props.getFloat("undefined"));
            else
                m_undefined = props.getSpectrum("undefined", Spectrum(0.0f));
        } else {
            m_undefined = Spectrum(0.0f);
        }

        if (SPECTRUM_SAMPLES != 3 && (m_field == EUV || m_field == EShadingNormal || m_field == EGeometricNormal
                || m_field == ERelativePosition || m_field == EPosition)) {
            Log(EError, "The field integrator implementation requires renderings to be done in RGB when "
                    "extracting positional data or surface normals / UV coordinates.");
        }
    }

    FieldIntegrator(Stream *stream, InstanceManager *manager)
     : SamplingIntegrator(stream, manager) {
         m_field = (EField) stream->readInt();
         m_undefined = Spectrum(stream);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        SamplingIntegrator::serialize(stream, manager);
        stream->writeInt((int) m_field);
        m_undefined.serialize(stream);
    }

    Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
        Spectrum result(m_undefined);

        if (!rRec.rayIntersect(ray))
            return result;

        Intersection &its = rRec.its;

        switch (m_field) {
            case EPosition:
                result.fromLinearRGB(its.p.x, its.p.y, its.p.z);
                break;
            case ERelativePosition: {
                    const Sensor *sensor = rRec.scene->getSensor();
                    const Transform &t = sensor->getWorldTransform()->eval(its.t).inverse();
                    Point p = t(its.p);
                    result.fromLinearRGB(p.x, p.y, p.z);
                }
                break;
            case EDistance:
                result = Spectrum(its.t);
                break;
            case EGeometricNormal:
                result.fromLinearRGB(its.geoFrame.n.x, its.geoFrame.n.y, its.geoFrame.n.z);
                break;
            case EShadingNormal:
                result.fromLinearRGB(its.shFrame.n.x, its.shFrame.n.y, its.shFrame.n.z);
                break;
            case EUV:
                result.fromLinearRGB(its.uv.x, its.uv.y, 0);
                break;
            case EAlbedo:
                result = its.shape->getBSDF()->getDiffuseReflectance(its);
                break;
            case EShapeIndex: {
                    const ref_vector<Shape> &shapes = rRec.scene->getShapes();
                    result = Spectrum((Float) -1);
                    for (size_t i=0; i<shapes.size(); ++i) {
                        if (shapes[i] == its.shape) {
                            result = Spectrum((Float) i);
                            break;
                        }
                    }
                }
                break;
            case EPrimIndex:
                result = Spectrum((Float) its.primIndex);
                break;
            default:
                Log(EError, "Internal error!");
        }

        return result;
    }

    std::string toString() const {
        return "FieldIntegrator[]";
    }

    MTS_DECLARE_CLASS()
private:
    EField m_field;
    Spectrum m_undefined;
};

MTS_IMPLEMENT_CLASS_S(FieldIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(FieldIntegrator, "Field extraction integrator");
MTS_NAMESPACE_END
