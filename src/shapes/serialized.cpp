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
#include <mitsuba/core/properties.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/lrucache.h>

#include <boost/make_shared.hpp>

/// How many files to keep open in the cache, per thread
#define MTS_SERIALIZED_CACHE_SIZE 4

MTS_NAMESPACE_BEGIN

/* Avoid having to include scenehandler.h */
extern MTS_EXPORT_RENDER void pushSceneCleanupHandler(void (*cleanup)());

/*!\plugin{serialized}{Serialized mesh loader}
 * \order{7}
 * \parameters{
 *     \parameter{filename}{\String}{
 *       Filename of the geometry file that should be loaded
 *     }
 *     \parameter{shapeIndex}{\Integer}{
 *         A \code{.serialized} file may contain several separate meshes. This parameter
 *         specifies which one should be loaded. \default{\code{0}, i.e. the first one}
 *     }
 *     \parameter{faceNormals}{\Boolean}{
 *       When set to \code{true}, any existing or computed vertex normals are
 *       discarded and \emph{face normals} will instead be used during rendering.
 *       This gives the rendered object a faceted appearance.\default{\code{false}}
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
 * }
 * The serialized mesh format represents the most space and time-efficient way
 * of getting geometry information into Mitsuba. It stores indexed triangle meshes
 * in a lossless gzip-based encoding that (after decompression) nicely matches up
 * with the internally used data structures. Loading such files is considerably
 * faster than the \pluginref{ply} plugin and orders of magnitude faster than
 * the \pluginref{obj} plugin. \vspace{-3mm}
 *
 * \paragraph{Format description:}
 * The \code{serialized} file format uses the little endian encoding, hence
 * all fields below should be interpreted accordingly. The contents are structured as
 * follows:
 * \begin{center}
 * \begin{longtable}{>{\bfseries}p{2cm}p{11cm}}
 * \toprule
 * Type & Content\\
 * \midrule
 * \code{uint16}&   File format identifier: \ \  \code{0x041C}\\
 * \code{uint16}&   File version identifier. Currently set to \ \  \code{0x0004}\\
 * \midrule
 * \multicolumn{2}{|c|}{\emph{From this point on, the stream is
 * compressed by the \code{DEFLATE} algorithm.}}\\
 * \multicolumn{2}{|c|}{\emph{The used encoding is that of
 * the \code{zlib} library.}}\\
 * \midrule
 * \code{uint32}&An 32-bit integer whose bits can be used
 * to specify the following flags:\\[-4mm]
 * &\begin{description}
 *  \item[\code{0x0001}] The mesh data includes per-vertex normals\vspace{-2mm}
 *  \item[\code{0x0002}] The mesh data includes texture coordinates\vspace{-2mm}
 *  \item[\code{0x0008}] The mesh data includes vertex colors\vspace{-2mm}
 *  \item[\code{0x0010}] Use face normals instead of smothly interpolated vertex
 *  normals. Equivalent to specifying \code{faceNormals=true} to the plugin. \vspace{-6mm}
 *  \item[\code{0x1000}] The subsequent content is represented in single precision\vspace{-2mm}
 *  \item[\code{0x2000}] The subsequent content is represented in double precision
 * \end{description}\\[-4mm]
 * \code{string}&A null-terminated string (utf-8), which denotes the name of the shape.\\
 * \code{uint64}&Number of vertices in the mesh\\
 * \code{uint64}&Number of triangles in the mesh\\
 * \code{array}&Array of all vertex positions (X, Y, Z, X, Y, Z, ...)
 * specified in binary single or double precision format (as denoted by the
 * flags)\\
 * \code{array}&Array of all vertex normal directions (X, Y, Z, X, Y, Z, ...)
 * specified in binary single or double precision format. When the mesh has
 * no vertex normals, this field is omitted.\\
 * \code{array}&Array of all vertex texture coordinates (U, V, U, V, ...)
 * specified in binary single or double precision format. When the mesh has
 * no texture coordinates, this field is omitted.\\
 * \code{array}&Array of all vertex colors (R, G, B, R, G, B, ...)
 * specified in binary single or double precision format. When the mesh has
 * no vertex colors, this field is omitted.\\
 * \code{array}&Indexed triangle data (\code{[i1, i2, i3]}, \code{[i1, i2, i3]}, ..)
 * specified in \code{uint32} or in \code{uint64} format (the latter
 * is used when the number of vertices exceeds \code{0xFFFFFFFF}).\\
 * \bottomrule
 * \end{longtable}
 * \end{center}
 * \paragraph{Multiple shapes:}
 * It is possible to store multiple meshes in a single \code{.serialized}
 * file. This is done by simply concatenating their data streams,
 * where every one is structured according to the above description.
 * Hence, after each mesh, the stream briefly reverts back to an
 * uncompressed format, followed by an uncompressed header, and so on.
 * This is neccessary for efficient read access to arbitrary sub-meshes.
 *
 * \paragraph{End-of-file dictionary:} In addition to the previous table,
 * a \code{.serialized} file also concludes with a brief summary at the end of
 * the file, which specifies the starting position of each sub-mesh:
 * \begin{center}
 * \begin{longtable}{>{\bfseries}p{2cm}p{11cm}}
 * \toprule
 * Type & Content\\
 * \midrule
 * \code{uint64}& File offset of the first mesh (in bytes)---this is always zero.\\
 * \code{uint64}& File offset of the second mesh\\
 * $\cdots$ & $\cdots$\\
 * \code{uint64}& File offset of the last sub-shape\\
 * \code{uint32}& Total number of meshes in the \code{.serialized} file\\
 * \bottomrule
 * \end{longtable}
 * \end{center}
 */
class SerializedMesh : public TriMesh {
public:
    SerializedMesh(const Properties &props) : TriMesh(props) {
        fs::path filePath = Thread::getThread()->getFileResolver()->resolve(
            props.getString("filename"));

        /* Object-space -> World-space transformation */
        Transform objectToWorld = props.getTransform("toWorld", Transform());

        /// When the file contains multiple meshes, this index specifies which one to load
        int shapeIndex = props.getInteger("shapeIndex", 0);
        AssertEx(shapeIndex >= 0, "Shape index must be nonnegative!");

        std::string name = (props.getID() != "unnamed") ? props.getID()
            : formatString("%s@%i", filePath.stem().string().c_str(), shapeIndex);

        /* Load the geometry */
        Log(EInfo, "Loading shape %i from \"%s\" ..", shapeIndex, filePath.filename().string().c_str());
        ref<Timer> timer = new Timer();
        loadCompressed(filePath, shapeIndex);
        Log(EDebug, "Done (" SIZE_T_FMT " triangles, " SIZE_T_FMT " vertices, %i ms)",
            m_triangleCount, m_vertexCount, timer->getMilliseconds());

        if (m_name.empty())
            m_name = name;

        /* By default, any existing normals will be used for
           rendering. If no normals are found, Mitsuba will
           automatically generate smooth vertex normals.
           Setting the 'faceNormals' parameter instead forces
           the use of face normals, which will result in a faceted
           appearance.
        */
        m_faceNormals = props.getBoolean("faceNormals", false);

        /* Causes all normals to be flipped */
        m_flipNormals = props.getBoolean("flipNormals", false);

        if (!objectToWorld.isIdentity()) {
            m_aabb.reset();
            for (size_t i=0; i<m_vertexCount; ++i) {
                Point p = objectToWorld(m_positions[i]);
                m_positions[i] = p;
                m_aabb.expandBy(p);
            }
            if (m_normals) {
                for (size_t i=0; i<m_vertexCount; ++i)
                    m_normals[i] = normalize(objectToWorld(m_normals[i]));
            }
        }

        if (objectToWorld.det3x3() < 0) {
            for (size_t i=0; i<m_triangleCount; ++i) {
                Triangle &t = m_triangles[i];
                std::swap(t.idx[0], t.idx[1]);
            }
        }

        if (props.hasProperty("maxSmoothAngle")) {
            if (m_faceNormals)
                Log(EError, "The properties 'maxSmoothAngle' and 'faceNormals' "
                "can't be specified at the same time!");
            rebuildTopology(props.getFloat("maxSmoothAngle"));
        }
    }

    SerializedMesh(Stream *stream, InstanceManager *manager)
        : TriMesh(stream, manager) { }

    MTS_DECLARE_CLASS()

private:

    /**
     * Helper class for loading serialized meshes from the same file
     * repeatedly: it is common for scene to load multiple meshes from the same
     * file, most times even in ascending order. This class loads the mesh
     * offsets dictionary only once and keeps the stream open.
     *
     * Instances of this class are not thread safe.
     */
    class MeshLoader {
    public:
        MeshLoader(const fs::path& filePath) {
            m_fstream = new FileStream(filePath, FileStream::EReadOnly);
            m_fstream->setByteOrder(Stream::ELittleEndian);
            const short version = SerializedMesh::readHeader(m_fstream);
            if (SerializedMesh::readOffsetDictionary(m_fstream,
                version, m_offsets) < 0) {
                // Assume there is a single mesh in the file at offset 0
                m_offsets.resize(1, 0);
            }
        }

        /**
         * Positions the stream at the location for the given shape index.
         * Returns the modified stream.
         */
        inline FileStream* seekStream(size_t shapeIndex) {
            if (shapeIndex > m_offsets.size()) {
                SLog(EError, "Unable to unserialize mesh, "
                    "shape index is out of range! (requested %i out of 0..%i)",
                    shapeIndex, (int) (m_offsets.size()-1));
            }
            const size_t pos = m_offsets[shapeIndex];
            m_fstream->seek(pos);
            return m_fstream;
        }

    private:
        std::vector<size_t> m_offsets;
        ref<FileStream> m_fstream;
    };

    typedef LRUCache<fs::path, std::less<fs::path>,
        boost::shared_ptr<MeshLoader> > MeshLoaderCache;

    class FileStreamCache : MeshLoaderCache {
    public:
        inline boost::shared_ptr<MeshLoader> get(const fs::path& path) {
            bool dummy;
            return MeshLoaderCache::get(path, dummy);
        }

        FileStreamCache() : MeshLoaderCache(MTS_SERIALIZED_CACHE_SIZE,
            &FileStreamCache::create) { }

    private:
        inline static boost::shared_ptr<MeshLoader> create(const fs::path &path) {
            return boost::make_shared<MeshLoader>(path);
        }
    };

    /// Release all currently held offset caches / file streams
    static void flushCache() {
        m_cache.set(NULL);
    }

    /// Loads the mesh from the thread-local file stream cache
    void loadCompressed(const fs::path& filePath, const int idx) {
        if (EXPECT_NOT_TAKEN(idx < 0)) {
            Log(EError, "Unable to unserialize mesh, "
                "shape index is negative! (requested %i out of 0..%i)", idx);
        }

        // Get the thread local cache; create it if this is the first time
        FileStreamCache* cache = m_cache.get();
        if (EXPECT_NOT_TAKEN(cache == NULL)) {
            cache = new FileStreamCache();
            m_cache.set(cache);
            mitsuba::pushSceneCleanupHandler(&SerializedMesh::flushCache);
        }

        boost::shared_ptr<MeshLoader> meshLoader = cache->get(filePath);
        Assert(meshLoader != NULL);
        TriMesh::loadCompressed(meshLoader->seekStream((size_t) idx));
    }

    static ThreadLocal<FileStreamCache> m_cache;
};

ThreadLocal<SerializedMesh::FileStreamCache> SerializedMesh::m_cache;

MTS_IMPLEMENT_CLASS_S(SerializedMesh, false, TriMesh)
MTS_EXPORT_PLUGIN(SerializedMesh, "Serialized mesh loader");
MTS_NAMESPACE_END
