/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

	This particular file is based on code by Piti Irawan, which is
	redistributed with permission.

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

#if !defined(__IRAWAN_H)
#define __IRAWAN_H

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/home/qi/numeric/real.hpp>
#include <boost/version.hpp>

#if BOOST_VERSION >= 105800
#include <boost/phoenix/bind/bind_member_variable.hpp>
#include <boost/phoenix/bind/bind_member_function.hpp>
#include <boost/phoenix/statement/if.hpp>
#else
#include <boost/spirit/home/phoenix/bind/bind_member_variable.hpp>
#include <boost/spirit/home/phoenix/bind/bind_member_function.hpp>
#include <boost/spirit/home/phoenix/statement/if.hpp>
#endif

MTS_NAMESPACE_BEGIN

namespace spirit = boost::spirit;
namespace qi = boost::spirit::qi;
namespace ph = boost::phoenix;

/// Data structure describing the properties of a single yarn
struct Yarn {
	enum EYarnType {
		EWarp = 0,
		EWeft = 1
	};

	/// Type of yarn (warp or weft)
	EYarnType type;
	/// Fiber twist angle
	Float psi;
	// Maximum inclination angle
	Float umax;
	/// Spine curvature
	Float kappa;
	/// Width of segment rectangle
	Float width;
	/// Length of segment rectangle
	Float length;
	/*! u coordinate of the yarn segment center, 
	 * assumes that the tile covers 0 <= u, v <= 1.
	 * (0, 0) is lower left corner of the weave pattern
	 */
	Float centerU;
	/// v coordinate of the yarn segment center
	Float centerV;
	/// Diffuse color
	Spectrum kd;
	/// Specular color
	Spectrum ks;

	Yarn() : type(EWarp),
		psi(0), umax(0), kappa(0), width(0), length(0),
		centerU(0), centerV(0), kd(0.0f), ks(0.0f) { }

	Yarn(Stream *stream) {
		type = (EYarnType) stream->readInt();
		psi = stream->readFloat();
		umax = stream->readFloat();
		kappa = stream->readFloat();
		width = stream->readFloat();
		length = stream->readFloat();
		centerU = stream->readFloat();
		centerV = stream->readFloat();
		kd = Spectrum(stream);
		ks = Spectrum(stream);
	}

	void serialize(Stream *stream) const {
		stream->writeInt(type);
		stream->writeFloat(psi);
		stream->writeFloat(umax);
		stream->writeFloat(kappa);
		stream->writeFloat(width);
		stream->writeFloat(length);
		stream->writeFloat(centerU);
		stream->writeFloat(centerV);
		kd.serialize(stream);
		ks.serialize(stream);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "yarn {" << endl 
			<< "  type = " << ((type == EWarp) ? "warp" : "weft") << "," << endl;
		if (psi != 0)
			oss << "  /* Fiber twist angle */" << endl
				<< "  psi = " << psi * 180 / M_PI << "," << endl;
		oss << "  /* Maximum inclination angle */" << endl
			<< "  umax = " << umax * 180 / M_PI << "," << endl;
		if (kappa != 0)
			oss << "  /* Spine curvature */" << endl
				<< "  kappa = " << kappa << "," << endl;
		oss << "  /* Width and length of the segment rectangle */" << endl
			<< "  width = " << width << "," << endl
			<< "  length = " << length << "," << endl
			<< "  /* Yarn segment center in tile space */" << endl
			<< "  centerU = " << centerU << "," << endl
			<< "  centerV = " << centerV << "," << endl
			<< "  /* Diffuse and specular color */" << endl;
		Float r, g, b;
		kd.toLinearRGB(r, g, b);
		oss << "  kd = {" << r << ", " << g << ", " << b << "}," << endl;
		ks.toLinearRGB(r, g, b);
		oss << "  ks = {" << r << ", " << g << ", " << b << "}," << endl
			<< "}";
		return oss.str();
	}
};

struct WeavePattern {
	/// Name of the weave pattern
	std::string name;
	/// Uniform scattering parameter
	Float alpha;
	/// Forward scattering parameter
	Float beta;
	/// Filament smoothing
	Float ss;
	/// Highlight width
	Float hWidth;
	/// Combined area taken up by the warp & weft
	Float warpArea, weftArea;

	/// Size of the weave pattern
	uint32_t tileWidth, tileHeight;

	/* Noise-related parameters */
	Float dWarpUmaxOverDWarp;
	Float dWarpUmaxOverDWeft;
	Float dWeftUmaxOverDWarp;
	Float dWeftUmaxOverDWeft;
	Float fineness, period;

	/// Detailed weave pattern
	std::vector<uint32_t> pattern;

	/// List of all yarns referenced in \c pattern
	std::vector<Yarn> yarns;

	inline WeavePattern() : name(""),
		alpha(0), beta(0), ss(0), hWidth(0),
		warpArea(0), weftArea(0), tileWidth(0), tileHeight(0),
		dWarpUmaxOverDWarp(0), dWarpUmaxOverDWeft(0),
		dWeftUmaxOverDWarp(0), dWeftUmaxOverDWeft(0),
		fineness(0), period(0) { }

	WeavePattern(Stream *stream) {
		name = stream->readString();
		alpha = stream->readFloat();
		beta = stream->readFloat();
		ss = stream->readFloat();
		hWidth = stream->readFloat();
		warpArea = stream->readFloat();
		weftArea = stream->readFloat();
		tileWidth = stream->readUInt();
		tileHeight = stream->readUInt();
		dWarpUmaxOverDWarp = stream->readFloat();
		dWarpUmaxOverDWeft = stream->readFloat();
		dWeftUmaxOverDWarp = stream->readFloat();
		dWeftUmaxOverDWeft = stream->readFloat();
		fineness = stream->readFloat();
		period = stream->readFloat();
		pattern.resize(tileWidth * tileHeight);
		stream->readUIntArray(&pattern[0], pattern.size());
		size_t yarnCount = stream->readSize();
		yarns.resize(yarnCount);
		for (size_t i=0; i<yarnCount; ++i)
			yarns[i] = Yarn(stream);
	}

	void serialize(Stream *stream) const {
		stream->writeString(name);
		stream->writeFloat(alpha);
		stream->writeFloat(beta);
		stream->writeFloat(ss);
		stream->writeFloat(hWidth);
		stream->writeFloat(warpArea);
		stream->writeFloat(weftArea);
		stream->writeUInt(tileWidth);
		stream->writeUInt(tileHeight);
		stream->writeFloat(dWarpUmaxOverDWarp);
		stream->writeFloat(dWarpUmaxOverDWeft);
		stream->writeFloat(dWeftUmaxOverDWarp);
		stream->writeFloat(dWeftUmaxOverDWeft);
		stream->writeFloat(fineness);
		stream->writeFloat(period);
		stream->writeUIntArray(&pattern[0], pattern.size());
		stream->writeSize(yarns.size());
		for (size_t i=0; i<yarns.size(); ++i)
			yarns[i].serialize(stream);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "weave {" << endl
			<< "  name = \"" << name << "\"," << endl << endl
			<< "  /* Tile size of the weave pattern */" << endl
			<< "  tileWidth = " << tileWidth << "," << endl
			<< "  tileHeight = " << tileHeight << "," << endl << endl
			<< "  /* Uniform and forward scattering parameters */" << endl
			<< "  alpha = " << alpha << "," << endl
			<< "  beta = " << beta << "," << endl << endl
			<< "  /* Filament smoothing */" << endl
			<< "  ss = " << ss << "," << endl << endl
			<< "  /* Highlight width */" << endl
			<< "  hWidth = " << hWidth << "," << endl << endl
			<< "  /* Combined warp/weft size */" << endl
			<< "  warpArea = " << warpArea << "," << endl
			<< "  weftArea = " << weftArea << "," << endl << endl
			<< "  /* Noise-related parameters */" << endl;
		if (dWarpUmaxOverDWarp != 0)
			oss << "  dWarpUmaxOverDWarp = " << dWarpUmaxOverDWarp * 180 / M_PI << "," << endl;
		if (dWeftUmaxOverDWeft != 0)
			oss << "  dWarpUmaxOverDWeft = " << dWarpUmaxOverDWeft * 180 / M_PI << "," << endl;
		if (dWarpUmaxOverDWarp != 0)
			oss << "  dWeftUmaxOverDWarp = " << dWeftUmaxOverDWarp * 180 / M_PI << "," << endl;
		if (dWeftUmaxOverDWeft != 0)
			oss << "  dWeftUmaxOverDWeft = " << dWeftUmaxOverDWeft * 180 / M_PI << "," << endl;
		if (fineness != 0)
			oss << "  fineness = " << fineness << "," << endl;
		if (period != 0)
			oss << "  period = " << period << "," << endl;
		oss << endl
			<< "  /* Weave pattern description */" << endl
			<< "  pattern {" << endl
			<< "    ";
		for (size_t i=0; i<pattern.size(); ++i) {
			oss << (int) pattern[i];
			if (i+1<pattern.size())
				oss << ", ";
		}
		oss << endl
			<< "  }," << endl
			<< endl
			<< "  /* Listing of all used yarns */" << endl;
		for (size_t i=0; i<yarns.size(); ++i) {
			oss << "  " << indent(yarns[i].toString());
			if (i+1<yarns.size())
				oss << "," << endl;
			oss << endl;
		}

		oss << "}";
		return oss.str();
	}
};

template <typename Iterator> struct SkipGrammar : qi::grammar<Iterator> {
	SkipGrammar () : SkipGrammar::base_type(start) {
		using qi::char_;
		using qi::space;
		using qi::eol;
		start = space | ("/*" >> *(char_ - "*/") >> "*/");
	}
	qi::rule<Iterator> start;
};

#if defined(SINGLE_PRECISION)
	#define Float_ qi::float_
#else
	#define Float_ qi::double_
#endif

template <typename Iterator> struct YarnGrammar : qi::grammar<Iterator, Yarn(), SkipGrammar<Iterator> > {
	YarnGrammar(const Properties &props) 
			: YarnGrammar::base_type(start), props(props) {
		using namespace qi::labels;
		using qi::float_;
		using qi::double_;
		using qi::char_;
		using qi::lit;
		using qi::_val;
		using qi::_1;
		using qi::_2;
		using qi::_3;

		type = (qi::string("warp") | qi::string("weft"))
			[ ph::if_else(_1 == "warp", _val = Yarn::EWarp, _val = Yarn::EWeft ) ];

		identifier = qi::lexeme[ lit('$') >> (qi::alpha | char_('_')) 
			>> *(qi::alnum | char_('_')) ];

		spec = ((lit("{") >> Float_ >> lit(",") >> Float_ >> lit(",") >> Float_ >> lit("}")) 
					[ ph::bind(&Spectrum::fromLinearRGB, _val, _1, _2, _3, Spectrum::EReflectance) ])
		     | (identifier [ _val = ph::bind(&Properties::getSpectrum, ph::ref(props), _1)]);

		flt = (Float_ [ _val = _1 ])
		    | (identifier [ _val = ph::bind(&Properties::getFloat, ph::ref(props), _1)]);

		start = lit("yarn")
			>> lit("{")
			>> (
				 (lit("type")     >> lit("=") >> type  [ ph::bind(&Yarn::type,    _val) = _1 ])
			   | (lit("psi")      >> lit("=") >> flt   [ ph::bind(&Yarn::psi,     _val) = _1 * M_PI / 180 ])
			   | (lit("umax")     >> lit("=") >> flt   [ ph::bind(&Yarn::umax,    _val) = _1 * M_PI / 180 ])
			   | (lit("kappa")    >> lit("=") >> flt   [ ph::bind(&Yarn::kappa,   _val) = _1 ])
			   | (lit("width")    >> lit("=") >> flt   [ ph::bind(&Yarn::width,   _val) = _1 ])
			   | (lit("length")   >> lit("=") >> flt   [ ph::bind(&Yarn::length,  _val) = _1 ])
			   | (lit("centerU")  >> lit("=") >> flt   [ ph::bind(&Yarn::centerU, _val) = _1 ])
			   | (lit("centerV")  >> lit("=") >> flt   [ ph::bind(&Yarn::centerV, _val) = _1 ])
			   | (lit("kd")       >> lit("=") >> spec  [ ph::bind(&Yarn::kd,      _val) = _1 ])
			   | (lit("ks")       >> lit("=") >> spec  [ ph::bind(&Yarn::ks,      _val) = _1 ])
			) % ','
			>> lit("}");
	}

	qi::rule<Iterator, Yarn::EYarnType(), SkipGrammar<Iterator> > type;
	qi::rule<Iterator, Yarn(), SkipGrammar<Iterator> > start;
	qi::rule<Iterator, Spectrum(), SkipGrammar<Iterator> > spec;
	qi::rule<Iterator, Float(), SkipGrammar<Iterator> > flt;
	qi::rule<Iterator, std::string()> identifier;
	const Properties &props;
};

template <typename Iterator> struct WeavePatternGrammar : qi::grammar<Iterator, WeavePattern(), SkipGrammar<Iterator> > {
	WeavePatternGrammar(const Properties &props) 
			: WeavePatternGrammar::base_type(start), yarn(props), props(props) {
		using namespace qi::labels;
		using qi::float_;
		using qi::double_;
		using qi::uint_;
		using qi::char_;
		using qi::lit;
		using qi::_val;
		using qi::_1;
		using ph::push_back;

		pattern = lit("pattern") >> lit("{")
				>> uint_ [ push_back(_val, _1) ] % ','
				>> lit("}");
		
		name = qi::lexeme [ lit("\"") >> *(char_ - "\"") >> lit("\"") ];

		identifier = qi::lexeme[ lit('$') >> (qi::alpha | char_('_')) 
			>> *(qi::alnum | char_('_')) ];

		flt = (Float_ [ _val = _1 ])
		      | (identifier [ _val = ph::bind(&Properties::getFloat, ph::ref(props), _1)]);

		start = lit("weave") >> lit("{") >> (
			  lit("name")               >> lit("=") >> name   [ ph::bind(&WeavePattern::name,               _val) = _1  ]
			| lit("tileWidth")          >> lit("=") >> uint_  [ ph::bind(&WeavePattern::tileWidth,          _val) = _1  ]
			| lit("tileHeight")         >> lit("=") >> uint_  [ ph::bind(&WeavePattern::tileHeight,         _val) = _1  ]
			| lit("ss")                 >> lit("=") >> flt    [ ph::bind(&WeavePattern::ss,                 _val) = _1  ]
			| lit("alpha")              >> lit("=") >> flt    [ ph::bind(&WeavePattern::alpha,              _val) = _1  ]
			| lit("beta")               >> lit("=") >> flt    [ ph::bind(&WeavePattern::beta,               _val) = _1  ]
			| lit("warpArea")           >> lit("=") >> flt    [ ph::bind(&WeavePattern::warpArea,           _val) = _1  ]
			| lit("weftArea")           >> lit("=") >> flt    [ ph::bind(&WeavePattern::weftArea,           _val) = _1  ]
			| lit("hWidth")             >> lit("=") >> flt    [ ph::bind(&WeavePattern::hWidth,             _val) = _1  ]
			| lit("dWarpUmaxOverDWarp") >> lit("=") >> flt    [ ph::bind(&WeavePattern::dWarpUmaxOverDWarp, _val) = _1 * M_PI / 180 ]
			| lit("dWarpUmaxOverDWeft") >> lit("=") >> flt    [ ph::bind(&WeavePattern::dWarpUmaxOverDWeft, _val) = _1 * M_PI / 180 ]
			| lit("dWeftUmaxOverDWarp") >> lit("=") >> flt    [ ph::bind(&WeavePattern::dWeftUmaxOverDWarp, _val) = _1 * M_PI / 180 ]
			| lit("dWeftUmaxOverDWeft") >> lit("=") >> flt    [ ph::bind(&WeavePattern::dWeftUmaxOverDWeft, _val) = _1 * M_PI / 180 ]
			| lit("fineness")           >> lit("=") >> flt    [ ph::bind(&WeavePattern::fineness,           _val) = _1  ]
			| lit("period")             >> lit("=") >> flt    [ ph::bind(&WeavePattern::period,             _val) = _1  ]
			| pattern                                         [ ph::bind(&WeavePattern::pattern,            _val) = _1  ] 
			| yarn                                            [ push_back(ph::bind(&WeavePattern::yarns,    _val),  _1) ] 
		) % ','
		>> lit("}") >> qi::eoi;
	}

	qi::rule<Iterator, WeavePattern(), SkipGrammar<Iterator> > start;
	qi::rule<Iterator, std::vector<uint32_t>(), SkipGrammar<Iterator> > pattern;
	qi::rule<Iterator, std::string(), SkipGrammar<Iterator> > name;
	qi::rule<Iterator, Float(), SkipGrammar<Iterator> > flt;
	qi::rule<Iterator, std::string()> identifier;
	YarnGrammar<Iterator> yarn;
	const Properties &props;
};


MTS_NAMESPACE_END

#endif /* __IRAWAN_H */
