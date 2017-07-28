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

#pragma once
#if !defined(__MITSUBA_HW_VPL_H_)
#define __MITSUBA_HW_VPL_H_

#include <mitsuba/render/vpl.h>
#include <mitsuba/hw/shadow.h>
#include <mitsuba/hw/gpugeometry.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief This class is responsible for the on-demand creation of
 * GPU shaders to render shapes that are illuminated by virtual
 * point light sources.
 *
 * This is used to drive the \c vpl integrator as well as the
 * interactive preview in the Mitsuba GUI.
 *
 * For each encountered BSDF-VPL pair, a custom piece of code
 * describing the characteristic light transport between them is
 * created and cached. The implementation carefully looks at
 * the tree of shader dependencies and creates code that can be
 * shared with other materials that have a similar configuration.
 * This is necessary to avoid generating a potentially huge O(N^2)
 * number of very similar programs and brings it down to O(n^2)
 * where n << N.
 *
 * \ingroup libhw
 */
class MTS_EXPORT_HW VPLShaderManager : public Object {
public:
    /// Create a new shader manager
    VPLShaderManager(Renderer *renderer);

    /// Upload requisite shaders to the GPU. To be called once
    void init();

    /// Issue a final cleanup (before destroying the shader manager)
    void cleanup();

    /**
     * \brief Associate the shader manager with a new scene
     *
     * This function uploads all relevant triangle geometry
     * to the GPU, and releases any resources held for the
     * previous scene. It is valid to call the function with
     * a \c NULL argument.
     */
    void setScene(const Scene *scene);

    /// Return the currently bound scene
    inline const Scene *getScene() const { return m_scene.get(); }

    /**
     * \brief Prepare the shader manager for rendering
     * with a new VPL
     *
     * Must be called after \ref setScene() and before
     * \ref bind(). This function creates a suitable
     * shadow map for the VPL.
     */
    void setVPL(const VPL &vpl);

    /**
     * \brief Bind a shader for rendering a certain
     * VPL/BSDF/Emitter triplet
     *
     * \param vpl
     *    Record describing the virtual point light source
     *
     * \param bsdf
     *    Material of the object to be rendered
     *
     * \param sensor
     *    Sensor from which it is viewed
     *
     * \param instanceTransform
     *    An additional object-to-world transformation that
     *    is applied to the rendered geometry. Used for instancing.
     *
     * \param faceNormals
     *    When set to \c true, a special shader is used to
     *    create face normals for the geometry
     */
    void bind(const VPL &vpl, const BSDF *bsdf,
        const Sensor *sensor, const Emitter *emitter,
        const Matrix4x4 &instanceTransform, bool faceNormals);

    /**
     * \brief Release the currently bound shader and
     * any resources (textures,..) that it references
     */
    void unbind();

    /**
     * \brief Convenience function for rendering all geometry
     * for a given VPL
     *
     * This function issues the necessary calls to \ref bind()
     * \ref unbind(), etc., and schedules draw calls for all
     * of the scene geometry.
     */
    void drawAllGeometryForVPL(const VPL &vpl, const Sensor *sensor);

    /// Draw the background if there is an environment emitter
    void drawBackground(const Sensor *sensor,
            const Transform &projectionTransform, Float scaleFactor);

    /// Set the clamping distance
    inline void setClamping(Float clamping) { m_clamping = clamping; }

    /// Return the clamping distance
    inline Float getClamping() const { return m_clamping; }

    /// Set the current shadow map resolution
    inline void setShadowMapResolution(int resolution) { m_shadowMapResolution = resolution; }

    /// Return the current shadow map resolution
    inline int getShadowMapResolution() const { return m_shadowMapResolution; }

    /// Set whether or not non-diffuse VPLs are used
    inline void setDiffuseSources(bool diffuseSources) { m_diffuseSources = diffuseSources; }

    /// Return whether or not non-diffuse VPLs are used
    inline bool getDiffuseSources() const { return m_diffuseSources; }

    /// Set whether or not surfaces are drawn assumed to be diffuse
    inline void setDiffuseReceivers(bool diffuseReceivers) { m_diffuseReceivers = diffuseReceivers; }

    /// Return whether or not surfaces are assumed to be diffuse
    inline bool getDiffuseReceivers() const { return m_diffuseReceivers; }

    /// Reset the VPL counter
    inline void resetCounter() { m_vplIndex = 0; }
protected:
    /// Virtual destructor
    virtual ~VPLShaderManager();

    /**
     * \brief This helper class stores a reference to a \ref Shader and all
     * sub-shaders that are part of its evaluation
     *
     * It also allows to generate GLSL code for the entire chain, and to
     * bind any uniform parameters to the target parameters
     */
    struct DependencyNode {
        const Shader *shader;
        std::vector<DependencyNode> children;
        std::vector<int> parameterIDs;

        /// Create from a \ref Shader object
        inline DependencyNode(Shader *shader = NULL) : shader(shader) {
            if (!shader)
                return;
            std::vector<Shader *> deps;
            shader->putDependencies(deps);
            for (std::vector<Shader *>::iterator it = deps.begin();
                it != deps.end(); ++it)
                children.push_back(DependencyNode(*it));
        }

        /// Copy constructor
        inline DependencyNode(const DependencyNode &node)
            : shader(node.shader), children(node.children),
              parameterIDs(node.parameterIDs) { }

        /// Generate GLSL code for the entire shader chain
        inline std::string generateCode(std::ostringstream &oss, int &id) const {
            std::vector<std::string> depNames;
            for (size_t i=0; i<children.size(); ++i)
                depNames.push_back(children[i].generateCode(oss, id));
            std::string evalName = formatString("shader_%i", id++);
            shader->generateCode(oss, evalName, depNames);
            oss << endl;
            return evalName;
        }

        /// Resolve all parameters of the shader chain
        inline void resolve(GPUProgram *program, int &id) {
            std::vector<std::string> depNames;
            for (size_t i=0; i<children.size(); ++i)
                children[i].resolve(program, id);

            std::string evalName = formatString("shader_%i", id++);
            shader->resolve(program, evalName, parameterIDs);
        }

        /// Bind all referenced resources (textures etc)
        inline void bind(GPUProgram *program, const DependencyNode &targetNode, int &textureUnitOffset) {
            if (!shader)
                return;
            for (size_t i=0; i<children.size(); ++i)
                children[i].bind(program, targetNode.children[i], textureUnitOffset);
            shader->bind(program, targetNode.parameterIDs, textureUnitOffset);
        }

        /// Release resources that were bound by \ref bind()
        inline void unbind() {
            if (!shader)
                return;
            shader->unbind();
            for (size_t i=0; i<children.size(); ++i)
                children[i].unbind();
        }

        /// Generate a textual summary of the entire shader chain
        inline void toString(std::ostringstream &oss) const {
            if (!shader)
                return;
            oss << shader->getClass()->getName();
            if (children.size() > 0) {
                oss << '[';
                for (size_t i=0; i<children.size(); ++i) {
                    children[i].toString(oss);
                    if (i+1<children.size())
                        oss << ", ";
                }
                oss << "]";
            }
        }

        inline std::string toString() const {
            std::ostringstream oss;
            toString(oss);
            return oss.str();
        }
    };

    /**
     * \brief Describes the configuration of a (vpl, bsdf, emitter)
     * shader chain triplet
     */
    struct VPLConfiguration {
        DependencyNode vpl, bsdf, emitter;
        bool faceNormals;
        GPUProgram *program;

        /* GLSL program paramter IDs */
        int param_instanceTransform, param_vplTransform;
        int param_vplPosition, param_vplDirection;
        int param_camPosition, param_camDirection;
        int param_shadowMap, param_depthRange;
        int param_emitterScale, param_vplPower;
        int param_vplFrame, param_vplUV, param_vplWi;
        int param_minDistSqr;

        /// Dummy constructor
        inline VPLConfiguration() : program(NULL) { }

        /// Create a new configuration for the given (vpl, bsdf, emitter) triplet
        inline VPLConfiguration(Shader *vpl, Shader *bsdf, Shader *emitter, bool faceNormals)
            : vpl(vpl), bsdf(bsdf), emitter(emitter), faceNormals(faceNormals), program(NULL) { }

        /// Generate GLSL code for the entire shader chain
        inline void generateCode(std::ostringstream &oss, std::string &vplEvalName,
                std::string &bsdfEvalName, std::string &emitterEvalName) const {
            int id = 0;
            vplEvalName = vpl.generateCode(oss, id);
            bsdfEvalName = bsdf.generateCode(oss, id);
            if (emitter.shader)
                emitterEvalName = emitter.generateCode(oss, id);
        }

        /// Resolve all parameters of the shader chain
        inline void resolve(GPUProgram *program) {
            int id = 0;
            vpl.resolve(program, id);
            bsdf.resolve(program, id);
            if (emitter.shader)
                emitter.resolve(program, id);
        }

        /// Bind all referenced resources (textures etc)
        inline void bind(const VPLConfiguration &targetConf, int textureUnitOffset) {
            vpl.bind(targetConf.program, targetConf.vpl, textureUnitOffset);
            bsdf.bind(targetConf.program, targetConf.bsdf, textureUnitOffset);
            if (emitter.shader)
                emitter.bind(targetConf.program, targetConf.emitter, textureUnitOffset);
        }

        /// Release resources that were bound by \ref bind()
        inline void unbind() {
            vpl.unbind();
            bsdf.unbind();
            if (emitter.shader)
                emitter.unbind();
        }

        /// Generate a textual summary of the entire shader chain
        inline std::string toString() const {
            std::ostringstream oss;
            oss << "vpl=";
            vpl.toString(oss);
            oss << ", bsdf=";
            bsdf.toString(oss);
            if (emitter.shader) {
                oss << ", emitter=";
                emitter.toString(oss);
            }
            if (faceNormals)
                oss << ", faceNormals";
            return oss.str();
        }
    };

    /**
     * \brief Order materials so that they can be drawn with the least
     * number of GPU pipeline flushes. Draw transparent objects last.
     */
    struct MaterialOrder {
        const std::vector<Renderer::TransformedGPUGeometry> &geo;

        MaterialOrder(const std::vector<Renderer::TransformedGPUGeometry> &geo)
            : geo(geo) { }

        inline bool operator()(size_t idx1, size_t idx2) const {
            const Shader *shader1 = geo[idx1].first->getShader();
            const Shader *shader2 = geo[idx2].first->getShader();

            if (shader1 && (shader1->getFlags() & Shader::ETransparent))
                shader1 = NULL;
            if (shader2 && (shader2->getFlags() & Shader::ETransparent))
                shader2 = NULL;

            return shader1 < shader2;
        }
    };

    /// Helper data structure to keep track of shapes that are undergoing keyframe animations
    struct AnimatedGeometryRecord {
        const AnimatedTransform *trafo;
        ssize_t geometryIndex;
        ssize_t opaqueGeometryIndex;

        AnimatedGeometryRecord(const AnimatedTransform *trafo,
            ssize_t geometryIndex, ssize_t opaqueGeometryIndex) :
            trafo(trafo), geometryIndex(geometryIndex),
            opaqueGeometryIndex(opaqueGeometryIndex) { }
    };

    MTS_DECLARE_CLASS()
private:
    ref<Renderer> m_renderer;
    ref<const Scene> m_scene;

    /* On-GPU geometry references */
    std::vector<Renderer::TransformedGPUGeometry> m_geometry;
    std::vector<Renderer::TransformedGPUGeometry> m_opaqueGeometry;
    std::vector<AnimatedGeometryRecord> m_animatedGeometry;

    /* Shader & dependency management */
    std::map<std::string, VPLConfiguration> m_configurations;
    VPLConfiguration m_targetConfiguration;
    VPLConfiguration m_currentProgram;

    /* Background rendering - related */
    ref<GPUProgram> m_backgroundProgram;
    DependencyNode m_backgroundDependencies;
    int m_backgroundParam_camPosition;
    int m_backgroundParam_camDirection;
    int m_backgroundParam_clipToWorld;
    int m_backgroundParam_emitterScale;

    /* Shadow map - related */
    ShadowMapGenerator::EShadowMapType m_shadowMapType;
    ref<ShadowMapGenerator> m_shadowGen;
    ref<GPUTexture> m_shadowMapCube;
    ref<GPUTexture> m_shadowMap2D;
    GPUTexture *m_shadowMap;
    Transform m_shadowMapTransform;
    Float m_nearClip, m_farClip;

    /* Other rendering parameters */
    bool m_diffuseSources, m_diffuseReceivers;
    int m_shadowMapResolution;
    uint32_t m_vplIndex;
    Float m_clamping, m_alpha;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_HW_VPL_H_ */
