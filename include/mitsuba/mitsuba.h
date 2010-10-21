/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__MITSUBA_H)
#define __MITSUBA_H

#include <mitsuba/core/platform.h>

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

/// Current release of Mitsuba
#define MTS_VERSION "0.2.0"

/// Year of this release
#define MTS_YEAR "2010"

/// Default port of <tt>mtssrv</tt>
#define MTS_DEFAULT_PORT 7554

#if defined(__LINUX__)
/// Default resource directory for packaged releases on Linux
#define MTS_RESOURCE_DIR "/usr/share/mitsuba"
#endif

using std::cout;
using std::cerr;
using std::endl;

/**
 * Include a basic subset of the core classes
 */
#include <mitsuba/core/constants.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/core/stl.h>
#include <mitsuba/core/object.h>
#include <mitsuba/core/ref.h>
#include <mitsuba/core/sse.h>
#include <mitsuba/core/tls.h>
#include <mitsuba/core/logger.h>
#include <mitsuba/core/thread.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/core/point.h>
#include <mitsuba/core/normal.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/util.h>

#endif /* __MITSUBA_H */
