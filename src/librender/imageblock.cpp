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

#include <mitsuba/render/imageblock.h>

MTS_NAMESPACE_BEGIN

/* ==================================================================== */
/*                              ImageBlock                              */
/* ==================================================================== */

ImageBlock::ImageBlock(const Vector2i &maxBlockSize, int borderSize, 
	bool supportWeights, bool supportAlpha, bool supportSnapshot, 
	bool supportStatistics) : border(borderSize), 
	  maxBlockSize(maxBlockSize), alpha(NULL), weights(NULL), 
	  variances(NULL), nSamples(NULL), pixelSnapshot(NULL), 
	  weightSnapshot(NULL), alphaSnapshot(NULL), extra(0) {

	int maxArraySize = (maxBlockSize.x + 2*border)
		*(maxBlockSize.y + 2*border);
	pixels = (Spectrum *) allocAligned(sizeof(Spectrum)*maxArraySize);

	if (supportAlpha)
		alpha = (Float *) allocAligned(sizeof(Float)*maxArraySize);

	if (supportWeights)
		weights = (Float *) allocAligned(sizeof(Float)*maxArraySize);

	if (supportSnapshot) {
		pixelSnapshot = (Spectrum *) allocAligned(sizeof(Spectrum)*(2*border+1)*(2*border+1));
		alphaSnapshot = (Float*) allocAligned(sizeof(Float)*(2*border+1)*(2*border+1));
		weightSnapshot = (Float *) allocAligned(sizeof(Float)*(2*border+1)*(2*border+1));
	}

	if (supportStatistics) {
		variances = (Spectrum *) allocAligned(sizeof(Spectrum)*maxArraySize);
		nSamples = (uint32_t *) allocAligned(sizeof(uint32_t)*maxArraySize);
	}
}

ImageBlock::~ImageBlock() {
	freeAligned(pixels);
	if (alpha)
		freeAligned(alpha);
	if (weights)
		freeAligned(weights);
	if (pixelSnapshot) {
		freeAligned(pixelSnapshot);
		freeAligned(alphaSnapshot);
		freeAligned(weightSnapshot);
	}
	if (variances) {
		freeAligned(variances);
		freeAligned(nSamples);
	}
}
	
void ImageBlock::clear() {
	int numEntries = fullSize.x*fullSize.y;

	memset(pixels, 0, sizeof(Spectrum) * numEntries);

	if (alpha)
		memset(alpha, 0, sizeof(Float) * numEntries);
	if (weights)
		memset(weights, 0, sizeof(Float) * numEntries);
	if (variances) {
		memset(variances, 0, sizeof(Spectrum) * numEntries);
		memset(nSamples, 0, sizeof(int) * numEntries);
	}
	extra = 0;
}

void ImageBlock::load(Stream *stream) {
	Assert(sizeof(Spectrum) == sizeof(Float)*SPECTRUM_SAMPLES);
	offset = Point2i(stream);
	size = Vector2i(stream);
	fullSize.x = size.x + 2*border;
	fullSize.y = size.y + 2*border;
	size_t nEntries = fullSize.x * fullSize.y;
	stream->readFloatArray(reinterpret_cast<Float *>(pixels), nEntries*SPECTRUM_SAMPLES);
	if (alpha)
		stream->readFloatArray(alpha, nEntries);
	if (weights)
		stream->readFloatArray(weights, nEntries);
	if (variances) {
		stream->readFloatArray(reinterpret_cast<Float *>(variances), nEntries*SPECTRUM_SAMPLES);
		stream->readUIntArray(nSamples, nEntries);
	}
	extra = stream->readInt();
}

void ImageBlock::save(Stream *stream) const {
	Assert(sizeof(Spectrum) == sizeof(Float)*SPECTRUM_SAMPLES);
	offset.serialize(stream);
	size.serialize(stream);
	size_t nEntries = fullSize.x * fullSize.y;
	stream->writeFloatArray(reinterpret_cast<Float *>(pixels), nEntries*SPECTRUM_SAMPLES);
	if (alpha)
		stream->writeFloatArray(alpha, nEntries);
	if (weights)
		stream->writeFloatArray(weights, nEntries);
	if (variances) {
		stream->writeFloatArray(reinterpret_cast<Float *>(variances), nEntries*SPECTRUM_SAMPLES);
		stream->writeUIntArray(nSamples, nEntries);
	}
	stream->writeInt(extra);
}

void ImageBlock::add(const ImageBlock *block) {
	int entry=0, imageY = block->offset.y-offset.y-1;

	Point2i topLeft = offset - Vector2i(border, border);
	Point2i bottomRight = offset + fullSize;

	for (int y=0; y<block->fullSize.y; ++y) {
		if (++imageY < topLeft.y || imageY >= bottomRight.y) {
			/// Skip a row if it is outside of the crop region
			entry += block->fullSize.x;
			continue;
		}

		int imageX = block->offset.x - offset.x - 1;
		for (int x=0; x<block->fullSize.x; ++x) {
			if (++imageX < topLeft.x || imageX >= bottomRight.x) {
				++entry;
				continue;
			}

			size_t idx = imageY*fullSize.x + imageX;

			pixels[idx] += block->pixels[entry];
			if (alpha != NULL)
				alpha[idx] += block->alpha[entry];
			if (weights != NULL)
				weights[idx] += block->weights[entry];
			entry++;
		}
	}
}

std::string ImageBlock::toString() const {
	std::ostringstream oss;
	oss << "ImageBlock[" << endl
		<< "\tborder = " << border << "," << endl
		<< "\toffset = " << offset.toString() << "," << endl
		<< "\tsize = " << size.toString() << "," << endl
		<< "\tfullSize = " << fullSize.toString() << "," << endl
		<< "\thasVariances = " << (variances != NULL) << "," << endl
		<< "\textra = " << extra << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS(ImageBlock, false, WorkResult)
MTS_NAMESPACE_END
