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
#if !defined(__MITSUBA_MITSUBA_H_)
#define __MITSUBA_MITSUBA_H_

#include <mitsuba/core/platform.h>
#include <boost/version.hpp>
#include <sstream>
#include <string>
#include <map>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <limits>

using std::cout;
using std::cerr;
using std::endl;

/**
 * Include a basic subset of the core classes
 */
#include <mitsuba/core/constants.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/core/math.h>
#include <mitsuba/core/object.h>
#include <mitsuba/core/ref.h>
#include <mitsuba/core/tls.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/thread.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/core/point.h>
#include <mitsuba/core/normal.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/util.h>

#endif /* __MITSUBA_MITSUBA_H_ */
