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

#if !defined(__IOR_DATA_H)
#define __IOR_DATA_H

#include <mitsuba/mitsuba.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

struct IOREntry {
	const char *name;
	Float value;
};

/**
 * Many values are taken from Hecht, Optics,
 * Fourth Edition.
 *
 * The IOR values are from measurements between
 * 0 and 20 degrees celsius at ~589 nm.
 */
static IOREntry iorData[] = {
	{ "vacuum",                1.0f      },
	{ "helium",                1.000036f },
	{ "hydrogen",              1.000132f },
	{ "air",                   1.000277f },
	{ "carbon dioxide",        1.00045f  },
	//////////////////////////////////////
	{ "water",                 1.3330f   },
	{ "acetone",               1.36f     },
	{ "ethanol",               1.361f    },
	{ "carbon tetrachloride",  1.461f    },
	{ "glycerol",              1.4729f   },
	{ "benzene",               1.501f    },
	{ "silicone oil",          1.52045f  },
	{ "bromine",               1.661f    },
	//////////////////////////////////////
	{ "water ice",             1.31f     },
	{ "fused quartz",          1.458f    },
	{ "pyrex",                 1.470f    },
	{ "acrylic glass",         1.49f     },
	{ "polypropylene",         1.49f     },
	{ "bk7",                   1.5046f   },
	{ "sodium chloride",       1.544f    },
	{ "amber",                 1.55f     },
	{ "pet",                   1.5750f   },
	{ "diamond",               2.419f    },

	{ NULL,                    0.0f      }
};

static Float lookupIOR(const std::string &name) {
	std::string lowerCase = boost::to_lower_copy(name);
	IOREntry *ior = iorData;

	while (ior->name) {
		if (lowerCase == ior->name)
			return ior->value;
		++ior;
	}

	std::ostringstream oss;
	oss << "Unable to find an IOR value for \"" << lowerCase
		<< "\"! Valid choices are:";

	/* Unable to find the IOR value by name -- print an error
	   message that lists all possible options */
	for (ior = iorData; ior->name != NULL; ++ior) {
		oss << ior->name;
		if ((ior+1)->name)
			oss << ", ";
	}

	SLog(EError, "%s", oss.str().c_str());
	return 0.0f;
}

inline Float lookupIOR(const Properties &props, const std::string &paramName, const std::string &defaultValue) {
	if (props.hasProperty(paramName) && props.getType(paramName) == Properties::EFloat)
		return props.getFloat(paramName);
	else
		return lookupIOR(props.getString(paramName, defaultValue));
}

inline Float lookupIOR(const Properties &props, const std::string &paramName, Float defaultValue) {
	if (props.hasProperty(paramName)) {
		if (props.getType(paramName) == Properties::EFloat)
			return props.getFloat(paramName);
		else
			return lookupIOR(props.getString(paramName));
	} else {
		return defaultValue;
	}
}

MTS_NAMESPACE_END

#endif /* __IOR_DATA_H */
