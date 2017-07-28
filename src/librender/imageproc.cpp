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

#include <mitsuba/render/imageproc.h>
#include <mitsuba/render/rectwu.h>

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                          BlockedImageProcess                         */
/* ==================================================================== */

void BlockedImageProcess::init(const Point2i &offset, const Vector2i &size, uint32_t blockSize) {
    m_offset = offset;
    m_size = size;
    m_blockSize = (int) blockSize;
    m_direction = ERight;
    m_numBlocks = Vector2i(
        (int) std::ceil((Float) size.x / (Float) blockSize),
        (int) std::ceil((Float) size.y / (Float) blockSize));
    m_numBlocksTotal = m_numBlocks.x * m_numBlocks.y;
    m_numBlocksGenerated = 0;
    m_curBlock = Point2i(m_numBlocks / 2);
    m_stepsLeft = 1;
    m_numSteps = 1;
}

ParallelProcess::EStatus BlockedImageProcess::generateWork(WorkUnit *unit, int worker) {
    /* Reimplementation of the spiraling block generator by Adam Arbree */
    RectangularWorkUnit &rect = *static_cast<RectangularWorkUnit *>(unit);

    if (m_numBlocksTotal == m_numBlocksGenerated)
        return EFailure;

    Point2i pos = m_curBlock * m_blockSize;
    rect.setOffset(pos + m_offset);
    rect.setSize(Vector2i(
        std::min(m_size.x-pos.x, m_blockSize),
        std::min(m_size.y-pos.y, m_blockSize)));

    if (++m_numBlocksGenerated == m_numBlocksTotal)
        return ESuccess;

    do {
        switch (m_direction) {
            case ERight: ++m_curBlock.x; break;
            case EDown:  ++m_curBlock.y; break;
            case ELeft:  --m_curBlock.x; break;
            case EUp:    --m_curBlock.y; break;
        }

        if (--m_stepsLeft == 0) {
            m_direction = (m_direction + 1) % 4;
            if (m_direction == ELeft || m_direction == ERight)
                ++m_numSteps;
            m_stepsLeft = m_numSteps;
        }
    } while (m_curBlock.x < 0 || m_curBlock.y < 0
        || m_curBlock.x >= m_numBlocks.x
        || m_curBlock.y >= m_numBlocks.y);

    return ESuccess;
}

MTS_IMPLEMENT_CLASS(BlockedImageProcess, true, ParallelProcess)
MTS_NAMESPACE_END
