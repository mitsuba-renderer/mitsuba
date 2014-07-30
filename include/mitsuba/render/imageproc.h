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
#if !defined(__MITSUBA_RENDER_IMAGEPROC_H_)
#define __MITSUBA_RENDER_IMAGEPROC_H_

#include <mitsuba/core/sched.h>

MTS_NAMESPACE_BEGIN

/**
 * Abstract parallel process, which performs a certain task (to be defined by
 * the subclass) on the pixels of an image where work on adjacent pixels
 * is independent. For preview purposes, a spiraling pattern of square
 * pixel blocks is generated.
 *
 * \ingroup librender
 */
class MTS_EXPORT_RENDER BlockedImageProcess : public ParallelProcess {
public:
	// ======================================================================
	//! @{ \name Implementation of the ParallelProcess interface
	// ======================================================================

	virtual EStatus generateWork(WorkUnit *unit, int worker);

	//! @}
	// ======================================================================

	MTS_DECLARE_CLASS()
protected:
	/**
	 * Initialize the image process
	 *
	 * \param offset
	 *    Integer offset of the image region to be processed
	 * \param size
	 *    Size of the image region to be processed
	 * \param blockSize
	 *    Size of the generated square pixel blocks
	 */
	void init(const Point2i &offset, const Vector2i &size, uint32_t blockSize);

	/// Protected constructor
	inline BlockedImageProcess() { }
	/// Virtual destructor
	virtual ~BlockedImageProcess() { }
protected:
	enum EDirection {
		ERight = 0,
		EDown,
		ELeft,
		EUp
	};

	Point2i m_offset;
	Vector2i m_size, m_numBlocks;
	Point2i m_curBlock;
	int m_direction, m_numSteps;
	int m_stepsLeft, m_numBlocksTotal;
	int m_numBlocksGenerated;
	int m_blockSize;
};

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_IMAGEPROC_H_ */
