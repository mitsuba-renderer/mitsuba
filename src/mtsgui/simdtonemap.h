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
#if !defined(__SIMDTONEMAP_H)
#define __SIMDTONEMAP_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/object.h>

using mitsuba::Class;

// SIMD CPU tone mapper, adapted from HDRI Tools
class TonemapCPU : public mitsuba::Object {

public:
	struct Params
	{
		float invGamma;
		float invWhitePoint;
		float multiplier;

		// key / averageLogLuminance
		float scale;

		bool isSRGB;

		float avgLogLum;
		float maxLum;

		Params() :
		invGamma(1.0f/2.2f), invWhitePoint(1.0f), multiplier(1.0f), scale(1.0f),
		isSRGB(true), avgLogLum(0.18f), maxLum(1.0f)
		{}
	};

	inline mitsuba::Float logAvgLuminance() const {
		return m_params.avgLogLum;
	}

	inline mitsuba::Float maxLuminance() const {
		return m_params.maxLum;
	}

	inline mitsuba::Float multiplier() const {
		return m_params.multiplier;
	}

	inline void setInvWhitePoint(mitsuba::Float invWhitePoint) {
		m_params.invWhitePoint = static_cast<float>(invWhitePoint);
	}

	inline void setInvGamma(mitsuba::Float invGamma) {
		m_params.invGamma = static_cast<float>(invGamma);
	}

	inline void setScale(mitsuba::Float scale) {
		m_params.scale = static_cast<float>(scale);
	}

	inline void setMultiplier(mitsuba::Float multiplier) {
		m_params.multiplier = static_cast<float>(multiplier);
	}

	inline void setSRGB(bool srgb) {
		m_params.isSRGB = srgb;
	}

	// Source: RGBA32F, Target: RGBA8
	bool gammaTonemap(const mitsuba::Bitmap* source, mitsuba::Bitmap* target) const;
	bool reinhardTonemap(const mitsuba::Bitmap* source, mitsuba::Bitmap* target) const;

	// Source: RGBA32F
	bool setLuminanceInfo(const mitsuba::Bitmap* source,
		mitsuba::Float multiplier = 1);

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~TonemapCPU() {}
private:
	Params m_params;
};

#endif // __SIMDTONEMAP_H
