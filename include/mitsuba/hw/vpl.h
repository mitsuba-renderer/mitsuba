/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

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

#if !defined(__VPL_HW_H)
#define __VPL_HW_H

#include <mitsuba/render/vpl.h>
#include <mitsuba/hw/renderer.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief This class is responsible for the on-demand creation of
 * GPU shaders to render meshes that are illuminated by a virtual 
 * point light source. 
 *
 * For each encountered BSDF-VPL pair, a custom piece of code 
 * describing the characteristic light transport between them is created 
 * and cached. To avoid generating a potentially huge (N squared) 
 * number of very similar programs, the implementation passes some 
 * properties using uniforms, in which case already existing code 
 * can be reused and we get something more like lower case n squared
 * (where lower n is the number of material types).
 *
 * \ingroup libhw
 */
class MTS_EXPORT_HW VPLShaderManager : public Object {
public:
	VPLShaderManager(const Scene *scene, Renderer *renderer);

	/// To be called once before use
	void init();

	/// Generate the shadow map for a particular VPL
	void setVPL(const VPL &vpl);

	/// Prepare for rendering a material with BSDF 'bsdf' illuminated by VPL 'vpl'.
	void configure(const VPL &vpl, const BSDF *bsdf, 
		const Luminaire *luminaire, const Point &camPos, bool faceNormals);

	/// Draw the background if there is an environment luminaire
	void drawBackground(const Transform &clipToWorld, const Point &camPos, Float scaleFactor);

	/// Release bound resources
	void unbind();

	/// Return all bound opaque triangle meshes
	inline const std::vector<std::pair<const TriMesh *, Transform> > &getMeshes() const { return m_meshes; }

	/// Return all bound transparent triangle meshes
	inline const std::vector<std::pair<const TriMesh *, Transform> > &getTransparentMeshes() const { return m_transparentMeshes; }

	/// Return the shadow cube map for debugging purposes
	inline GPUTexture *getShadowMap() { return m_shadowMap; }

	/// Should the shadow map generation be done in a single pass? (requires geometry shader support)
	inline void setSinglePass(bool singlePass) { m_singlePass = singlePass; }
	
	/// Is shadow map generation generation performed in a single pass?
	inline bool isSinglePass() const { return m_singlePass; }

	/// Set whether or not non-diffuse VPLs are used
	inline void setDiffuseSources(bool diffuseSources) { m_diffuseSources = diffuseSources; }

	/// Return whether or not non-diffuse VPLs are used
	inline bool getDiffuseSources() const { return m_diffuseSources; }
	
	/// Set whether or not surfaces are drawn assumed to be diffuse
	inline void setDiffuseReceivers(bool diffuseReceivers) { m_diffuseReceivers = diffuseReceivers; }

	/// Return whether or not surfaces are assumed to be diffuse
	inline bool getDiffuseReceivers() const { return m_diffuseReceivers; }

	/// Set the current shadow map resolution
	inline void setShadowMapResolution(int resolution) { m_shadowMapResolution = resolution; }

	/// Return the current shadow map resolution
	inline int getShadowMapResolution() const { return m_shadowMapResolution; }

	/// Set the max. shadow map far plane distance
	inline void setMaxClipDist(Float maxClipDist) { m_maxClipDist = maxClipDist; }
	
	/// Return the max. shadow map far plane distance
	inline Float getMaxClipDist() const { return m_maxClipDist; }

	/// Set the clamping distance
	inline void setClamping(Float clamping) { m_clamping = clamping; }

	/// Return the clamping distance
	inline Float getClamping() const { return m_clamping; }

	/// Return the associated scene
	inline const Scene *getScene() const { return m_scene.get(); }

	/// To be called once before destruction
	void cleanup();

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~VPLShaderManager();
private:
	struct VPLDependencyNode {
		Shader *shader;
		std::vector<VPLDependencyNode> children;
		std::vector<int> parameterIDs;

		inline VPLDependencyNode(const VPLDependencyNode &node) 
			: shader(node.shader), children(node.children), parameterIDs(node.parameterIDs) {
		}

		inline VPLDependencyNode(Shader *shader = NULL) : shader(shader) {
			if (shader == NULL)
				return;
			std::vector<Shader *> deps;
			shader->putDependencies(deps);
			for (std::vector<Shader *>::iterator it = deps.begin();
				it != deps.end(); ++it)
				children.push_back(VPLDependencyNode(*it));
		}
		
		std::string recursiveGenerateCode(std::ostringstream &oss, int &id) const {
			std::vector<std::string> depNames;
			for (size_t i=0; i<children.size(); ++i)
				depNames.push_back(children[i].recursiveGenerateCode(oss, id));
			std::string evalName = formatString("shader_%i", id++);
			shader->generateCode(oss, evalName, depNames);
			oss << endl;
			return evalName;
		}

		void recursiveResolve(GPUProgram *program, int &id) {
			std::vector<std::string> depNames;
			for (size_t i=0; i<children.size(); ++i)
				children[i].recursiveResolve(program, id);

			std::string evalName = formatString("shader_%i", id++);
			shader->resolve(program, evalName, parameterIDs);
		}

		void recursiveBind(GPUProgram *program, const VPLDependencyNode &targetNode, int &textureUnitOffset) {
			for (size_t i=0; i<children.size(); ++i)
				children[i].recursiveBind(program, targetNode.children[i], textureUnitOffset);
			shader->bind(program, targetNode.parameterIDs, textureUnitOffset);
		}

		void recursiveUnbind() {
			shader->unbind();
			for (size_t i=0; i<children.size(); ++i)
				children[i].recursiveUnbind();
		}

		inline void toString(std::ostringstream &oss) const {
			oss << shader->getClass()->getName();
			if (children.size() > 0) {
				oss << '{';
				for (size_t i=0; i<children.size(); ++i) {
					children[i].toString(oss);
					if (i+1<children.size())
						oss << ',';
				}
				oss << "}";
			}
		}
	};

	struct VPLProgramConfiguration {
		VPLDependencyNode vpl, bsdf, luminaire;
		bool hasLuminaire, faceNormals;
		int param_shadowMap, param_vplPos, param_camPos, param_vplPower;
		int param_vplN, param_vplS, param_vplT, param_vplWi, param_vplUV;
		int param_nearClip, param_invClipRange, param_minDist;
		int param_diffuseSources, param_diffuseReceivers, param_alpha;

		inline VPLProgramConfiguration() { }

		inline VPLProgramConfiguration(Shader *vpl, Shader *bsdf, Shader *luminaire, bool faceNormals)
			: vpl(vpl), bsdf(bsdf), luminaire(luminaire), faceNormals(faceNormals) {
			hasLuminaire = (luminaire != NULL);
		}

		void generateCode(std::ostringstream &oss, std::string &vplEvalName,
				std::string &bsdfEvalName, std::string &luminaireEvalName) const {
			int id = 0;
			vplEvalName = vpl.recursiveGenerateCode(oss, id);
			bsdfEvalName = bsdf.recursiveGenerateCode(oss, id);
			if (hasLuminaire)
				luminaireEvalName = luminaire.recursiveGenerateCode(oss, id);
		}

		void resolve(GPUProgram *program) {
			int id = 0;
			vpl.recursiveResolve(program, id);
			bsdf.recursiveResolve(program, id);
			if (hasLuminaire)
				luminaire.recursiveResolve(program, id);
		}
	
		inline void bind(GPUProgram *program, const VPLProgramConfiguration &targetConf, int &textureUnitOffset) {
			vpl.recursiveBind(program, targetConf.vpl, textureUnitOffset);
			bsdf.recursiveBind(program, targetConf.bsdf, textureUnitOffset);
			if (hasLuminaire)
				luminaire.recursiveBind(program, targetConf.luminaire, textureUnitOffset);
		}
		
		inline void unbind() {
			vpl.recursiveUnbind();
			bsdf.recursiveUnbind();
			if (hasLuminaire)
				luminaire.recursiveUnbind();
		}

		inline void toString(std::ostringstream &oss) const {
			oss << "vpl=";
			vpl.toString(oss);
			oss << ", bsdf=";
			bsdf.toString(oss);
			if (hasLuminaire) {
				oss << ", luminaire=";
				luminaire.toString(oss);
			}
			if (faceNormals)
				oss << ", faceNormals";
		}
	};

	struct ProgramAndConfiguration {
		GPUProgram *program;
		VPLProgramConfiguration config;

		inline ProgramAndConfiguration() : program(NULL) {
		}

		inline ProgramAndConfiguration(GPUProgram *program, 
			const VPLProgramConfiguration &config)
			: program(program), config(config) {
		}
	};

	/* General */
	ref<const Scene> m_scene;
	ref<Renderer> m_renderer;
	Float m_clamping, m_minDist;
	Float m_maxClipDist;
	bool m_initialized;
	
	/* Shadow mapping related */
	ref<GPUProgram> m_shadowProgram;
	ref<GPUProgram> m_altShadowProgram;
	int m_shadowProgramParam_cubeMapTransform[6];
	int m_shadowProgramParam_depthVec[6];
	int m_altShadowProgramParam_cubeMapTransform;
	int m_altShadowProgramParam_depthVec;
	ref<GPUTexture> m_shadowMap;
	Float m_nearClip, m_invClipRange;
	int m_shadowMapResolution;
	bool m_singlePass;
	bool m_diffuseSources, m_diffuseReceivers;

	/* Rendering related */
	std::map<std::string, ProgramAndConfiguration> m_programs;
	ProgramAndConfiguration m_current;
	VPLProgramConfiguration m_targetConfig;
	ref<GPUProgram> m_backgroundProgram;
	VPLDependencyNode m_backgroundDependencies;
	std::vector<std::pair<const TriMesh *, Transform> > m_meshes;
	std::vector<std::pair<const TriMesh *, Transform> > m_transparentMeshes;
	std::vector<std::pair<const GPUGeometry *, Transform> > m_drawList;
};

MTS_NAMESPACE_END

#endif /* __VPL_HW_H */
