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

MTS_NAMESPACE_BEGIN

static Float CubeData_vertexPositions[][3] = {{1, -1, -1}, {1, -1, 1}, {-1, -1, 1}, {-1, -1, -1}, {1, 1, -1}, {-1, 1, -1}, {-1, 1, 1}, {1, 1, 1}, {1, -1, -1}, {1, 1, -1}, {1, 1, 1}, {1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {-1, -1, 1}, {-1, 1, 1}, {-1, 1, -1}, {-1, -1, -1}, {1, 1, -1}, {1, -1, -1}, {-1, -1, -1}, {-1, 1, -1}};

static Float CubeData_vertexNormals[][3] = {{0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {1, 0, 0}, {1, 0, 0}, {1, 0, 0}, {1, 0, 0}, {0, 0, 1}, {0, 0, 1}, {0, 0, 1}, {0, 0, 1}, {-1, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}};

static Float CubeData_texcoords[][2] = {{0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}};

static uint32_t CubeData_triangles[][3] = {{0, 1, 2}, {3, 0, 2}, {4, 5, 6}, {7, 4, 6}, {8, 9, 10}, {11, 8, 10}, {12, 13, 14}, {15, 12, 14}, {16, 17, 18}, {19, 16, 18}, {20, 21, 22}, {23, 20, 22}};

/*!\plugin{cube}{Cube intersection primitive}
 * \order{0}
 * \parameters{
 *     \parameter{toWorld}{\Transform\Or\Animation}{
 *        Specifies an optional linear object-to-world transformation.
 *        \default{none (i.e. object space $=$ world space)}
 *     }
 *     \parameter{flipNormals}{\Boolean}{
 *        Is the cube inverted, i.e. should the normal vectors
 *        be flipped? \default{\code{false}, i.e. the normals point outside}
 *     }
 * }
 *
 * \renderings{
 *     \rendering{Basic example}
 *         {shape_cube_basic}
 *     \rendering{A textured and stretched cube with the default parameterization
 *      (Listing~\ref{lst:cube-example})}
 *         {shape_cube_parameterization}
 * }
 *
 * This shape plugin describes a simple cube/cuboid intersection primitive. By
 * default, it creates a cube between the world-space positions $(-1, -1, -1)$ and $(1,1,1)$.
 * However, an arbitrary linear transformation may be specified to translate, rotate, scale
 * or skew it as desired. The parameterization of this shape maps every face onto the
 * rectangle $[0,1]^2$ in $uv$ space.
 * \vspace{5mm}
 * \begin{xml}[caption={Example of a textured and stretched cube}, label=lst:cube-example]
 * <shape type="cube">
 *   <transform name="toWorld">
 *     <scale z="2"/>
 *   </transform>
 *
 *   <bsdf type="diffuse">
 *     <texture type="checkerboard" name="reflectance">
 *       <float name="uvscale" value="6"/>
 *     </texture>
 *   </bsdf>
 * </shape>
 * \end{xml}
 */
class Cube : public TriMesh {
public:
    Cube(const Properties &props) : TriMesh(props) {
        m_triangleCount = 12;
        m_vertexCount = 24;
        m_positions = new Point[m_vertexCount];
        m_texcoords = new Point2[m_vertexCount];
        m_normals = new Normal[m_vertexCount];
        m_triangles = new Triangle[m_triangleCount];

        Transform toWorld = props.getTransform("toWorld", Transform());
        for (uint32_t i=0; i<m_vertexCount; ++i) {
            Normal n;
            Point p;
            Point2 uv;

            for (int j=0; j<3; ++j) {
                p[j] = CubeData_vertexPositions[i][j];
                n[j] = CubeData_vertexNormals[i][j];
            }
            for (int j=0; j<2; ++j)
                uv[j] = CubeData_texcoords[i][j];

            m_positions[i] = toWorld(p);
            m_normals[i] = normalize(toWorld(n));
            m_texcoords[i] = uv;
        }

        for (uint32_t i=0; i<m_triangleCount; ++i)
            for (int j=0; j<3; ++j)
                m_triangles[i].idx[j] = CubeData_triangles[i][j];
    }

    Cube(Stream *stream, InstanceManager *manager)
        : TriMesh(stream, manager) { }

    MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS_S(Cube, false, TriMesh)
MTS_EXPORT_PLUGIN(Cube, "Cube intersection primitive");
MTS_NAMESPACE_END
