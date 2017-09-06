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

#include <mitsuba/render/trimesh.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/timer.h>
#include <ply/ply_parser.hpp>
#include <functional>

MTS_NAMESPACE_BEGIN

/*!\plugin{ply}{PLY (Stanford Triangle Format) mesh loader}
 * \order{6}
 * \parameters{
 *     \parameter{filename}{\String}{
 *       Filename of the PLY file that should be loaded
 *     }
 *     \parameter{faceNormals}{\Boolean}{
 *       When set to \code{true}, Mitsuba will use face normals when rendering
 *       the object, which will give it a faceted appearance. \default{\code{false}}
 *     }
 *     \parameter{maxSmoothAngle}{\Float}{
 *       When specified, Mitsuba will discard all vertex normals in the input mesh and rebuild
 *       them in a way that is sensitive to the presence of creases and corners. For more
 *       details on this parameter, see page~\pageref{sec:maxSmoothAngle}. Disabled by default.
 *     }
 *     \parameter{flipNormals}{\Boolean}{
 *       Optional flag to flip all normals. \default{\code{false}, i.e.
 *       the normals are left unchanged}.
 *     }
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional linear object-to-world transformation.
 *        \default{none (i.e. object space $=$ world space)}
 *     }
 *     \parameter{srgb}{\Boolean}{
 *       When set to \code{true}, any vertex colors will be interpreted as sRGB,
 *       instead of linear RGB \default{\code{true}}.
 *     }
 * }
 * \renderings{
 *     \rendering{The PLY plugin is useful for loading large geometry. (Dragon
 *         statue courtesy of XYZ RGB)}{shape_ply_dragon}
 *     \rendering{The Stanford bunny loaded with \code{faceNormals=true}. Note
 *         the faceted appearance.}{shape_ply_bunny}
 * }
 * This plugin implements a fast loader for the Stanford PLY format (both
 * the ASCII and binary format). It is based on the \code{libply} library by
 * Ares Lagae (\url{http://people.cs.kuleuven.be/~ares.lagae/libply}).
 * The current plugin implementation supports triangle meshes with optional
 * UV coordinates, vertex normals, and vertex colors.
 *
 * When loading meshes that contain vertex colors, note that they need to be
 * explicitly referenced in a BSDF using a special texture named
 * \pluginref{vertexcolors}.
 */
class PLYLoader : public TriMesh {
public:
    PLYLoader(const Properties &props) : TriMesh(props) {
        fs::path filePath = Thread::getThread()->getFileResolver()->resolve(
            props.getString("filename"));
        m_name = filePath.stem().string();

        /* Determines whether vertex colors should be
           treated as linear RGB or sRGB. */
        m_sRGB = props.getBoolean("srgb", true);

        /* Object-space -> World-space transformation */
        m_objectToWorld = props.getTransform("toWorld", Transform());

        /* Load the geometry */
        Log(EInfo, "Loading geometry from \"%s\" ..", filePath.filename().string().c_str());
        if (!fs::exists(filePath))
            Log(EError, "PLY file \"%s\" could not be found!", filePath.string().c_str());

        m_triangleCount = m_vertexCount = 0;
        m_vertexCtr = m_faceCount = m_faceCtr = m_indexCtr = 0;
        m_normal = Normal(0.0f);
        m_uv = Point2(0.0f);
        m_hasNormals = false;
        m_hasTexCoords = false;
        memset(&m_face, 0, sizeof(uint32_t)*4);
        loadPLY(filePath);

        if (m_triangleCount == 0 || m_vertexCount == 0)
            Log(EError, "Unable to load \"%s\" (no triangles or vertices found)!");

        Assert(m_faceCtr   == m_faceCount);
        Assert(m_vertexCtr == m_vertexCount);

        if (props.hasProperty("maxSmoothAngle")) {
            if (m_faceNormals)
                Log(EError, "The properties 'maxSmoothAngle' and 'faceNormals' "
                "can't be specified at the same time!");
            rebuildTopology(props.getFloat("maxSmoothAngle"));
        }

        if (m_triangleCount < m_faceCount * 2) {
            /* Needed less memory than the earlier conservative estimate -- free it! */
            Triangle *temp = new Triangle[m_triangleCount];
            memcpy(temp, m_triangles, sizeof(Triangle) * m_triangleCount);
            delete[] m_triangles;
            m_triangles = temp;
        }
    }


    PLYLoader(Stream *stream, InstanceManager *manager) : TriMesh(stream, manager) { }

    void loadPLY(const fs::path &path);

    void info_callback(const std::string& filename, std::size_t line_number,
            const std::string& message) {
        Log(EInfo, "\"%s\" [line %i] info: %s", filename.c_str(), line_number,
            message.c_str());
    }

    void warning_callback(const std::string& filename, std::size_t line_number,
            const std::string& message) {
        Log(EWarn, "\"%s\" [line %i] warning: %s", filename.c_str(), line_number,
            message.c_str());
    }

    void error_callback(const std::string& filename, std::size_t line_number,
            const std::string& message) {
        Log(EError, "\"%s\" [line %i] error: %s", filename.c_str(), line_number,
            message.c_str());
    }

    template<typename ValueType> std::function <void (ValueType)>
        scalar_property_definition_callback(const std::string& element_name,
        const std::string& property_name);

    template<typename SizeType, typename IndexType> std::tuple<std::function<void (SizeType)>,
        std::function<void (IndexType)>, std::function<void ()> >
        list_property_definition_callback(const std::string& element_name,
        const std::string& property_name);

    std::tuple<std::function<void()>, std::function<void()> >
        element_definition_callback(const std::string& element_name, std::size_t count) {
        if (element_name == "vertex") {
            m_vertexCount = count;
            m_positions = new Point[m_vertexCount];
            return std::tuple<std::function<void()>,
                std::function<void()> >(
                    std::bind(&PLYLoader::vertex_begin_callback, this),
                    std::bind(&PLYLoader::vertex_end_callback, this)
            );
        } else if (element_name == "face") {
            m_faceCount = count;
            m_triangles = new Triangle[m_faceCount*2];
            return std::tuple<std::function<void()>,
                std::function<void()> >(
                std::bind(&PLYLoader::face_begin_callback, this),
                std::bind(&PLYLoader::face_end_callback, this)
            );
        } else {
            return
                std::tuple<std::function<void()>,
                std::function<void()> >(nullptr, nullptr);
        }
    }

    void vertex_x_callback(ply::float32 x) { m_position.x = x; }
    void vertex_y_callback(ply::float32 y) { m_position.y = y; }
    void vertex_z_callback(ply::float32 z) { m_position.z = z; }
    void normal_x_callback(ply::float32 x) {
        if (!m_normals)
            m_normals = new Normal[m_vertexCount];
        m_normal.x = x;
    }
    void normal_y_callback(ply::float32 y) { m_normal.y = y; }
    void normal_z_callback(ply::float32 z) { m_normal.z = z; }
    void texcoord_u_callback(ply::float32 x) {
        if (!m_texcoords)
            m_texcoords = new Point2[m_vertexCount];
        m_uv.x = x;
    }
    void texcoord_v_callback(ply::float32 y) { m_uv.y = y; }

    inline Float fromSRGBComponent(Float value) {
        if (value <= (Float) 0.04045)
            return value / (Float) 12.92;
        return std::pow((value + (Float) 0.055)
            / (Float) (1.0 + 0.055), (Float) 2.4);
    }

    void vertex_begin_callback() { }
    void vertex_end_callback() {
        Point p = m_objectToWorld(m_position);
        m_aabb.expandBy(p);
        m_positions[m_vertexCtr] = p;
        if (m_normals)
            m_normals[m_vertexCtr] = normalize(m_objectToWorld(m_normal));
        if (m_texcoords)
            m_texcoords[m_vertexCtr] = m_uv;
        if (m_colors) {
            if (m_sRGB)
                m_colors[m_vertexCtr] = Color3(
                    fromSRGBComponent(m_red),
                    fromSRGBComponent(m_green),
                    fromSRGBComponent(m_blue));
            else
                m_colors[m_vertexCtr] = Color3(m_red, m_green, m_blue);
        }
        m_vertexCtr++;
    }

    void red_callback_uint8(ply::uint8 r) {
        if (!m_colors)
            m_colors = new Color3[m_vertexCount];
        m_red = r / 255.0f;
    }
    void green_callback_uint8(ply::uint8 g) { m_green = g / 255.0f; }
    void blue_callback_uint8(ply::uint8 b) { m_blue = b / 255.0f; }
    void red_callback(ply::float32 r) {
        if (!m_colors)
            m_colors = new Color3[m_vertexCount];
        m_red = r;
    }
    void green_callback(ply::float32 g) { m_green = g; }
    void blue_callback(ply::float32 b) { m_blue = b; }

    /* Face colors are unsupported */
    void face_red_callback_uint8(ply::uint8 r) { }
    void face_green_callback_uint8(ply::uint8 g) { }
    void face_blue_callback_uint8(ply::uint8 b) { }
    void face_red_callback(ply::float32 r) { }
    void face_green_callback(ply::float32 g) { }
    void face_blue_callback(ply::float32 b) { }

    void face_begin_callback() { }
    void face_end_callback() { }

    void face_vertex_indices_begin_uint8(ply::uint8 size) {
        if (size != 3 && size != 4)
            Log(EError, "Encountered a face with %i vertices! "
                "Only triangle and quad-based PLY meshes are supported for now.", size);
        m_indexCtr = 0;
    }

    void face_vertex_indices_begin_uint32(ply::uint32 size) {
        if (size != 3 && size != 4)
            Log(EError, "Only triangle and quad-based PLY meshes are supported for now.");
        m_indexCtr = 0;
    }

    void face_vertex_indices_element_int32(ply::int32 element) {
        Assert(m_indexCtr < 4);
        Assert((size_t) element < m_vertexCount);
        m_face[m_indexCtr++] = element;
    }

    void face_vertex_indices_element_uint32(ply::uint32 element) {
        Assert(m_indexCtr < 4);
        Assert((size_t) element < m_vertexCount);
        m_face[m_indexCtr++] = element;
    }

    void face_vertex_indices_end() {
        Assert(m_indexCtr == 3 || m_indexCtr == 4);

        Triangle t;
        t.idx[0] = m_face[0]; t.idx[1] = m_face[1]; t.idx[2] = m_face[2];
        m_triangles[m_triangleCount++] = t;

        if (m_indexCtr == 4) {
            t.idx[0] = m_face[3]; t.idx[1] = m_face[0]; t.idx[2] = m_face[2];
            m_triangles[m_triangleCount++] = t;
        }

        m_faceCtr++;
    }

    MTS_DECLARE_CLASS()
private:
    Point m_position;
    Normal m_normal;
    Float m_red, m_green, m_blue;
    Transform m_objectToWorld;
    size_t m_faceCount, m_vertexCtr;
    size_t m_faceCtr, m_indexCtr;
    uint32_t m_face[4];
    bool m_hasNormals, m_hasTexCoords;
    Point2 m_uv;
    bool m_sRGB;
};

template<> std::function <void (ply::float32)>
    PLYLoader::scalar_property_definition_callback(const std::string& element_name,
    const std::string& property_name) {
    if (element_name == "vertex") {
        if (property_name == "x") {
            return std::bind(&PLYLoader::vertex_x_callback, this,  _1);
        } else if (property_name == "y") {
            return std::bind(&PLYLoader::vertex_y_callback, this,  _1);
        } else if (property_name == "z") {
            return std::bind(&PLYLoader::vertex_z_callback, this,  _1);
        } else if (property_name == "nx") {
            m_hasNormals = true;
            return std::bind(&PLYLoader::normal_x_callback, this,  _1);
        } else if (property_name == "ny") {
            return std::bind(&PLYLoader::normal_y_callback, this,  _1);
        } else if (property_name == "nz") {
            return std::bind(&PLYLoader::normal_z_callback, this,  _1);
        } else if (property_name == "u" || property_name == "texture_u" || property_name == "s") {
            m_hasTexCoords = true;
            return std::bind(&PLYLoader::texcoord_u_callback, this,  _1);
        } else if (property_name == "v" || property_name == "texture_v" || property_name == "t") {
            return std::bind(&PLYLoader::texcoord_v_callback, this,  _1);
        } else if (property_name == "diffuse_red" || property_name == "red") {
            return std::bind(&PLYLoader::red_callback, this,  _1);
        } else if (property_name == "diffuse_green" || property_name == "green") {
            return std::bind(&PLYLoader::green_callback, this,  _1);
        } else if (property_name == "diffuse_blue" || property_name == "blue") {
            return std::bind(&PLYLoader::blue_callback, this,  _1);
        }
    } else if (element_name == "face") {
        if (property_name == "diffuse_red" || property_name == "red") {
            return std::bind(&PLYLoader::face_red_callback, this,  _1);
        } else if (property_name == "diffuse_green" || property_name == "green") {
            return std::bind(&PLYLoader::face_green_callback, this,  _1);
        } else if (property_name == "diffuse_blue" || property_name == "blue") {
            return std::bind(&PLYLoader::face_blue_callback, this,  _1);
        }
    }
    return nullptr;
}

template<> std::function <void (ply::uint8)>
    PLYLoader::scalar_property_definition_callback(const std::string& element_name,
    const std::string& property_name) {
    if (element_name == "vertex") {
        if (property_name == "diffuse_red" || property_name == "red") {
            return std::bind(&PLYLoader::red_callback_uint8, this,  _1);
        } else if (property_name == "diffuse_green" || property_name == "green") {
            return std::bind(&PLYLoader::green_callback_uint8, this,  _1);
        } else if (property_name == "diffuse_blue" || property_name == "blue") {
            return std::bind(&PLYLoader::blue_callback_uint8, this,  _1);
        }
    } else if (element_name == "face") {
        if (property_name == "diffuse_red" || property_name == "red") {
            return std::bind(&PLYLoader::face_red_callback_uint8, this,  _1);
        } else if (property_name == "diffuse_green" || property_name == "green") {
            return std::bind(&PLYLoader::face_green_callback_uint8, this,  _1);
        } else if (property_name == "diffuse_blue" || property_name == "blue") {
            return std::bind(&PLYLoader::face_blue_callback_uint8, this,  _1);
        }
    }
    return nullptr;
}

template<> std::tuple<std::function<void (ply::uint8)>,
    std::function<void (ply::int32)>, std::function<void ()> >
    PLYLoader::list_property_definition_callback(const std::string& element_name,
    const std::string& property_name) {
    if ((element_name == "face") && (property_name == "vertex_indices" || property_name == "vertex_index")) {
        return std::tuple<std::function<void (ply::uint8)>,
            std::function<void (ply::int32)>, std::function<void ()> >(
            std::bind(&PLYLoader::face_vertex_indices_begin_uint8, this, _1),
            std::bind(&PLYLoader::face_vertex_indices_element_int32, this, _1),
            std::bind(&PLYLoader::face_vertex_indices_end, this)
        );
    } else {
        return std::tuple<std::function<void (ply::uint8)>,
            std::function<void (ply::int32)>,
            std::function<void ()> >(nullptr, nullptr, nullptr);
    }
}

template<> std::tuple<std::function<void (ply::uint32)>,
    std::function<void (ply::int32)>, std::function<void ()> >
    PLYLoader::list_property_definition_callback(const std::string& element_name,
    const std::string& property_name) {
    if ((element_name == "face") && (property_name == "vertex_indices" || property_name == "vertex_index")) {
        return std::tuple<std::function<void (ply::uint32)>,
            std::function<void (ply::int32)>, std::function<void ()> >(
            std::bind(&PLYLoader::face_vertex_indices_begin_uint32, this, _1),
            std::bind(&PLYLoader::face_vertex_indices_element_int32, this, _1),
            std::bind(&PLYLoader::face_vertex_indices_end, this)
        );
    } else {
        return std::tuple<std::function<void (ply::uint32)>,
            std::function<void (ply::int32)>,
            std::function<void ()> >(0, 0, 0);
    }
}

template<> std::tuple<std::function<void (ply::uint8)>,
    std::function<void (ply::uint32)>, std::function<void ()> >
    PLYLoader::list_property_definition_callback(const std::string& element_name,
    const std::string& property_name) {
    if ((element_name == "face") && (property_name == "vertex_indices" || property_name == "vertex_index")) {
        return std::tuple<std::function<void (ply::uint8)>,
            std::function<void (ply::uint32)>, std::function<void ()> >(
            std::bind(&PLYLoader::face_vertex_indices_begin_uint8, this, _1),
            std::bind(&PLYLoader::face_vertex_indices_element_uint32, this, _1),
            std::bind(&PLYLoader::face_vertex_indices_end, this)
        );
    } else {
        return std::tuple<std::function<void (ply::uint8)>,
            std::function<void (ply::uint32)>,
            std::function<void ()> >(0, 0, 0);
    }
}

template<> std::tuple<std::function<void (ply::uint32)>,
    std::function<void (ply::uint32)>, std::function<void ()> >
    PLYLoader::list_property_definition_callback(const std::string& element_name,
    const std::string& property_name) {
    if ((element_name == "face") && (property_name == "vertex_indices" || property_name == "vertex_index")) {
        return std::tuple<std::function<void (ply::uint32)>,
            std::function<void (ply::uint32)>, std::function<void ()> >(
            std::bind(&PLYLoader::face_vertex_indices_begin_uint32, this, _1),
            std::bind(&PLYLoader::face_vertex_indices_element_uint32, this, _1),
            std::bind(&PLYLoader::face_vertex_indices_end, this)
        );
    } else {
        return std::tuple<std::function<void (ply::uint32)>,
            std::function<void (ply::uint32)>,
            std::function<void ()> >(0, 0, 0);
    }
}


void PLYLoader::loadPLY(const fs::path &path) {
    ply::ply_parser ply_parser;
    ply_parser.info_callback(std::bind(&PLYLoader::info_callback,
        this, std::ref(m_name), _1, _2));
    ply_parser.warning_callback(std::bind(&PLYLoader::warning_callback,
        this, std::ref(m_name), _1, _2));
    ply_parser.error_callback(std::bind(&PLYLoader::error_callback,
        this, std::ref(m_name), _1, _2));

    ply_parser.element_definition_callback(std::bind(&PLYLoader::element_definition_callback,
        this, _1, _2));

    ply::ply_parser::scalar_property_definition_callbacks_type scalar_property_definition_callbacks;
    ply::ply_parser::list_property_definition_callbacks_type list_property_definition_callbacks;

    ply::at<ply::float32>(scalar_property_definition_callbacks) = std::bind(
        &PLYLoader::scalar_property_definition_callback<ply::float32>, this, _1, _2);

    ply::at<ply::uint8>(scalar_property_definition_callbacks) = std::bind(
        &PLYLoader::scalar_property_definition_callback<ply::uint8>, this, _1, _2);

    ply::at<ply::uint8, ply::int32>(list_property_definition_callbacks) = std::bind(
        &PLYLoader::list_property_definition_callback<ply::uint8, ply::int32>, this, _1, _2);

    ply::at<ply::uint32, ply::int32>(list_property_definition_callbacks) = std::bind(
        &PLYLoader::list_property_definition_callback<ply::uint32, ply::int32>, this, _1, _2);

    ply::at<ply::uint8, ply::uint32>(list_property_definition_callbacks) = std::bind(
        &PLYLoader::list_property_definition_callback<ply::uint8, ply::uint32>, this, _1, _2);

    ply::at<ply::uint32, ply::uint32>(list_property_definition_callbacks) = std::bind(
        &PLYLoader::list_property_definition_callback<ply::uint32, ply::uint32>, this, _1, _2);

    ply_parser.scalar_property_definition_callbacks(scalar_property_definition_callbacks);
    ply_parser.list_property_definition_callbacks(list_property_definition_callbacks);

    ref<Timer> timer = new Timer();
    ply_parser.parse(path.string());

    size_t vertexSize = sizeof(Point);
    if (m_normals)
        vertexSize += sizeof(Normal);
    if (m_colors)
        vertexSize += sizeof(Spectrum);
    if (m_texcoords)
        vertexSize += sizeof(Point2);

    Log(EInfo, "\"%s\": Loaded " SIZE_T_FMT " triangles, " SIZE_T_FMT
            " vertices (%s in %i ms).", m_name.c_str(), m_triangleCount, m_vertexCount,
            memString(sizeof(uint32_t) * m_triangleCount * 3 + vertexSize * m_vertexCount).c_str(),
            timer->getMilliseconds());
}


MTS_IMPLEMENT_CLASS_S(PLYLoader, false, TriMesh)
MTS_EXPORT_PLUGIN(PLYLoader, "PLY mesh loader");
MTS_NAMESPACE_END
