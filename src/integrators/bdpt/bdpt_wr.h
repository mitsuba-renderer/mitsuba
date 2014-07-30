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

#if !defined(__BDPT_WR_H)
#define __BDPT_WR_H

#include <mitsuba/render/imageblock.h>
#include <mitsuba/core/fresolver.h>
#include "bdpt.h"

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                             Work result                              */
/* ==================================================================== */

/**
   Bidirectional path tracing needs its own WorkResult implementation,
   since each rendering thread simultaneously renders to a small 'camera
   image' block and potentially a full-resolution 'light image'.
*/
class BDPTWorkResult : public WorkResult {
public:
	BDPTWorkResult(const BDPTConfiguration &conf, const ReconstructionFilter *filter,
			Vector2i blockSize = Vector2i(-1, -1));

	// Clear the contents of the work result
	void clear();

	/// Fill the work result with content acquired from a binary data stream
	virtual void load(Stream *stream);

	/// Serialize a work result to a binary data stream
	virtual void save(Stream *stream) const;

	/// Aaccumulate another work result into this one
	void put(const BDPTWorkResult *workResult);

#if BDPT_DEBUG == 1
	/* In debug mode, this function allows to dump the contributions of
	   the individual sampling strategies to a series of images */
	void dump(const BDPTConfiguration &conf,
			const fs::path &prefix, const fs::path &stem) const;

	inline void putDebugSample(int s, int t, const Point2 &sample,
			const Spectrum &spec) {
		m_debugBlocks[strategyIndex(s, t)]->put(sample, (const Float *) &spec);
	}
#endif

	inline void putSample(const Point2 &sample, const Spectrum &spec) {
		m_block->put(sample, spec, 1.0f);
	}

	inline void putLightSample(const Point2 &sample, const Spectrum &spec) {
		m_lightImage->put(sample, spec, 1.0f);
	}

	inline const ImageBlock *getImageBlock() const {
		return m_block.get();
	}

	inline const ImageBlock *getLightImage() const {
		return m_lightImage.get();
	}

	inline void setSize(const Vector2i &size) {
		m_block->setSize(size);
	}

	inline void setOffset(const Point2i &offset) {
		m_block->setOffset(offset);
	}

	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~BDPTWorkResult();

	inline int strategyIndex(int s, int t) const {
		int above = s+t-2;
		return s + above*(5+above)/2;
	}
protected:
#if BDPT_DEBUG == 1
	ref_vector<ImageBlock> m_debugBlocks;
#endif
	ref<ImageBlock> m_block, m_lightImage;
};

MTS_NAMESPACE_END

#endif /* __BDPT_WR_H */
