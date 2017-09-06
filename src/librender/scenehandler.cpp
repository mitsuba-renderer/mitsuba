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

#include <mitsuba/core/platform.h>

// Mitsuba's "Assert" macro conflicts with Xerces' XSerializeEngine::Assert(...).
// This becomes a problem when using a PCH which contains mitsuba/core/logger.h
#if defined(Assert)
# undef Assert
#endif

#include <xercesc/parsers/SAXParser.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/sax/Locator.hpp>
#include <mitsuba/render/scenehandler.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/scene.h>
#include <boost/algorithm/string.hpp>
#include <boost/unordered_set.hpp>

MTS_NAMESPACE_BEGIN
XERCES_CPP_NAMESPACE_USE

#define TRANSCODE_BLOCKSIZE 2048

#define XMLLog(level, fmt, ...) Thread::getThread()->getLogger()->log(\
    level, NULL, __FILE__, __LINE__, "In file \"%s\" (near line %i): " fmt, \
    m_locator ? transcode(m_locator->getSystemId()).c_str() : "<unknown>", \
    m_locator ? m_locator->getLineNumber() : -1, \
    ## __VA_ARGS__)

typedef void (*CleanupFun) ();
typedef boost::unordered_set<CleanupFun> CleanupSet;
static PrimitiveThreadLocal<CleanupSet> __cleanup_tls;

SceneHandler::SceneHandler(const ParameterMap &params,
    NamedObjectMap *namedObjects, bool isIncludedFile) : m_params(params),
        m_namedObjects(namedObjects), m_isIncludedFile(isIncludedFile) {
    m_pluginManager = PluginManager::getInstance();
    m_locator = NULL;

    if (m_isIncludedFile) {
        SAssert(namedObjects != NULL);
    } else {
        SAssert(namedObjects == NULL);
        m_namedObjects = new NamedObjectMap();
    }

#if !defined(WIN32)
    setlocale(LC_NUMERIC, "C");
#endif

    /* Create a mapping from tag names to tag IDs and associated classes */
    m_tags["scene"]      = TagEntry(EScene,      MTS_CLASS(Scene));
    m_tags["shape"]      = TagEntry(EShape,      MTS_CLASS(Shape));
    m_tags["sampler"]    = TagEntry(ESampler,    MTS_CLASS(Sampler));
    m_tags["film"]       = TagEntry(EFilm,       MTS_CLASS(Film));
    m_tags["integrator"] = TagEntry(EIntegrator, MTS_CLASS(Integrator));
    m_tags["texture"]    = TagEntry(ETexture,    MTS_CLASS(Texture));
    m_tags["sensor"]     = TagEntry(ESensor,     MTS_CLASS(Sensor));
    m_tags["emitter"]    = TagEntry(EEmitter,    MTS_CLASS(Emitter));
    m_tags["subsurface"] = TagEntry(ESubsurface, MTS_CLASS(Subsurface));
    m_tags["medium"]     = TagEntry(EMedium,     MTS_CLASS(Medium));
    m_tags["volume"]     = TagEntry(EVolume,     MTS_CLASS(VolumeDataSource));
    m_tags["phase"]      = TagEntry(EPhase,      MTS_CLASS(PhaseFunction));
    m_tags["bsdf"]       = TagEntry(EBSDF,       MTS_CLASS(BSDF));
    m_tags["rfilter"]    = TagEntry(ERFilter,    MTS_CLASS(ReconstructionFilter));
    m_tags["null"]       = TagEntry(ENull,       (Class *) NULL);
    m_tags["ref"]        = TagEntry(EReference,  (Class *) NULL);
    m_tags["integer"]    = TagEntry(EInteger,    (Class *) NULL);
    m_tags["float"]      = TagEntry(EFloat,      (Class *) NULL);
    m_tags["boolean"]    = TagEntry(EBoolean,    (Class *) NULL);
    m_tags["string"]     = TagEntry(EString,     (Class *) NULL);
    m_tags["translate"]  = TagEntry(ETranslate,  (Class *) NULL);
    m_tags["rotate"]     = TagEntry(ERotate,     (Class *) NULL);
    m_tags["lookat"]     = TagEntry(ELookAt,     (Class *) NULL);
    m_tags["lookAt"]     = TagEntry(ELookAt,     (Class *) NULL);
    m_tags["scale"]      = TagEntry(EScale,      (Class *) NULL);
    m_tags["matrix"]     = TagEntry(EMatrix,     (Class *) NULL);
    m_tags["point"]      = TagEntry(EPoint,      (Class *) NULL);
    m_tags["vector"]     = TagEntry(EVector,     (Class *) NULL);
    m_tags["rgb"]        = TagEntry(ERGB,        (Class *) NULL);
    m_tags["srgb"]       = TagEntry(ESRGB,       (Class *) NULL);
    m_tags["blackbody"]  = TagEntry(EBlackBody,  (Class *) NULL);
    m_tags["spectrum"]   = TagEntry(ESpectrum,   (Class *) NULL);
    m_tags["transform"]  = TagEntry(ETransform,  (Class *) NULL);
    m_tags["animation"]  = TagEntry(EAnimation,  (Class *) NULL);
    m_tags["include"]    = TagEntry(EInclude,    (Class *) NULL);
    m_tags["alias"]      = TagEntry(EAlias,      (Class *) NULL);
    m_tags["default"]    = TagEntry(EDefault,    (Class *) NULL);

    XMLTransService::Codes failReason;
    m_transcoder = XMLPlatformUtils::fgTransService->makeNewTranscoderFor(
            "UTF-8", failReason, TRANSCODE_BLOCKSIZE);
}

SceneHandler::~SceneHandler() {
    delete m_transcoder;
    clear();
    if (!m_isIncludedFile)
        delete m_namedObjects;
}

void SceneHandler::setDocumentLocator(const xercesc::Locator* const locator) {
    m_locator = locator;
}

void SceneHandler::clear() {
    if (!m_isIncludedFile) {
        for (NamedObjectMap::iterator it = m_namedObjects->begin();
                it != m_namedObjects->end(); ++it)
            if (it->second)
                it->second->decRef();
        m_namedObjects->clear();
    }
}

std::string SceneHandler::transcode(const XMLCh * input) const {
    XMLSize_t charsToBeConsumed = XMLString::stringLen(input);
    char output[TRANSCODE_BLOCKSIZE + 4];
    XMLSize_t totalCharsConsumed = 0;
    std::string result;

    while (totalCharsConsumed < charsToBeConsumed) {
        XMLSize_t charsConsumed = 0;
        XMLSize_t charsProduced = m_transcoder->transcodeTo(input,
            std::min((XMLSize_t) 2048, charsToBeConsumed - totalCharsConsumed),
            (XMLByte *) output, TRANSCODE_BLOCKSIZE, charsConsumed,
            XMLTranscoder::UnRep_RepChar);

        totalCharsConsumed += charsConsumed;
        output[charsProduced] = '\0';
        input += charsConsumed;
        result += output;
    }

    return result;
}

// -----------------------------------------------------------------------
//  Implementation of the SAX DocumentHandler interface
// -----------------------------------------------------------------------

void SceneHandler::startDocument() {
    clear();
}

void SceneHandler::endDocument() {
    SAssert(m_scene != NULL);

    /* Call cleanup handlers */
    CleanupSet &cleanup = __cleanup_tls.get();
    for (CleanupSet::iterator it = cleanup.begin();
            it != cleanup.end(); ++it)
        (*it)();
    cleanup.clear();
}

void SceneHandler::characters(const XMLCh* const name,
        const XMLSize_t length) {
    std::string value = trim(transcode(name));
    if (value != "")
        XMLLog(EWarn, "Unexpected character data: %s", value.c_str());
}

Float SceneHandler::parseFloat(const std::string &name,
        const std::string &str, Float defVal) const {
    char *end_ptr = NULL;
    if (str.empty()) {
        if (defVal == -1)
            XMLLog(EError, "Missing floating point value (in <%s>)", name.c_str());
        return defVal;
    }

    Float result = (Float) std::strtod(str.c_str(), &end_ptr);
    if (*end_ptr != '\0')
        XMLLog(EError, "Invalid floating point value specified (in <%s>)", name.c_str());
    return result;
}

void SceneHandler::startElement(const XMLCh* const xmlName,
    AttributeList &xmlAttributes) {
    std::string name = transcode(xmlName);
    TagMap::const_iterator it = m_tags.find(name);

    if (it == m_tags.end())
        XMLLog(EError, "Unhandled tag \"%s\" encountered!", name.c_str());

    const TagEntry &tag = it->second;
    ParseContext context((name == "scene") ? NULL : &m_context.top(), tag.first);

    for (size_t i=0; i<xmlAttributes.getLength(); i++) {
        std::string attrValue = transcode(xmlAttributes.getValue(i));
        if (attrValue.length() > 0 && attrValue.find('$') != attrValue.npos) {
            for (ParameterMap::const_reverse_iterator it = m_params.rbegin(); it != m_params.rend(); ++it) {
                std::string::size_type pos = 0;
                std::string searchString = "$" + it->first;
                while ((pos = attrValue.find(searchString, pos)) != std::string::npos) {
                    attrValue.replace(pos, searchString.size(), it->second);
                    ++pos;
                }
            }
            if (attrValue.find('$') != attrValue.npos && attrValue.find('[') == attrValue.npos)
                XMLLog(EError, "The scene referenced an undefined parameter: \"%s\"", attrValue.c_str());
        }

        context.attributes[transcode(xmlAttributes.getName(i))] = attrValue;
    }

    switch (tag.first) {
        case EScene: {
                std::string versionString = context.attributes["version"];
                if (versionString == "")
                    throw VersionException(formatString("The requested scene cannot be loaded, since it "
                        "is missing version information! Since Mitsuba 0.3.0, it is "
                        "mandatory that scene XML files specify the version of Mitsuba "
                        "that was used at the time of their creation.\nThis makes it clear "
                        "how to interpret them in the presence of a changing file format. "
                        "The version should be specified within the 'scene' tag, "
                        "e.g.\n\t<scene version=\"" MTS_VERSION "\">\n"
                        "Please update your scene file with the right version number and try reloading it."),
                        Version());
                Version fileVersion(versionString), currentVersion(MTS_VERSION);
                if (!fileVersion.isCompatible(currentVersion)) {
                    if (fileVersion < currentVersion) {
                        throw VersionException(formatString("The requested scene is from an older version of Mitsuba "
                            "(file version: %s, current version: %s), hence the loading process was stopped. "
                            "Please open the scene from within Mitsuba's graphical user interface (mtsgui) -- "
                            "it will then be upgraded to the current format.",
                            fileVersion.toString().c_str(), MTS_VERSION), fileVersion);
                    } else {
                        XMLLog(EError, "The requested scene is from an incompatible future version of Mitsuba "
                            "(file version: %s, current version: %s). Giving up.",
                            fileVersion.toString().c_str(), MTS_VERSION);
                    }
                }
            }
            break;
        case ETransform:
            m_transform = Transform();
            break;
        case EAnimation: {
                m_animatedTransform = new AnimatedTransform();
            }
            break;
        default:
            break;
    }

    m_context.push(context);
}

void pushSceneCleanupHandler(void (*cleanup)()) {
    __cleanup_tls.get().insert(cleanup);
}

void SceneHandler::endElement(const XMLCh* const xmlName) {
    std::string name = transcode(xmlName);
    ParseContext &context = m_context.top();
    std::string type = boost::to_lower_copy(context.attributes["type"]);
    context.properties.setPluginName(type);
    if (context.attributes.find("id") != context.attributes.end())
        context.properties.setID(context.attributes["id"]);

    ref<ConfigurableObject> object;

    TagMap::const_iterator it = m_tags.find(name);
    if (it == m_tags.end())
        XMLLog(EError, "Unhandled tag \"%s\" encountered!", name.c_str());

    const TagEntry &tag = it->second;

    switch (tag.first) {
        case EScene:
            object = m_scene = new Scene(context.properties);
            break;

        case ENull:
            object = NULL;
            break;

        case EReference: {
                std::string id = context.attributes["id"];
                if (m_namedObjects->find(id) == m_namedObjects->end())
                    XMLLog(EError, "Referenced object '%s' not found!", id.c_str());
                object = (*m_namedObjects)[id];
            }
            break;

        case EInteger: {
                char *end_ptr = NULL;
#ifdef WIN32
                int64_t i = _strtoi64(context.attributes["value"].c_str(), &end_ptr, 10);
#else
                int64_t i = strtoll(context.attributes["value"].c_str(), &end_ptr, 10);
#endif
                if (*end_ptr != '\0')
                    XMLLog(EError, "Invalid integer value specified (in <%s>)",
                        context.attributes["name"].c_str());
                context.parent->properties.setLong(context.attributes["name"], i);
            }
            break;

        case EFloat: {
                Float value = parseFloat(name, context.attributes["value"]);
                context.parent->properties.setFloat(context.attributes["name"], value);
            }
            break;

        case EBoolean: {
                const std::string &str = context.attributes["value"];
                bool value;
                if (str == "true") {
                    value = true;
                } else if (str == "false") {
                    value = false;
                } else {
                    XMLLog(EError, "Unsupported boolean constant '%s' -- must be "
                        "'true' or 'false'!", str.c_str());
                    return;
                }
                context.parent->properties.setBoolean(context.attributes["name"], value);
            }
            break;

        case EString: {
            context.parent->properties.setString(context.attributes["name"],
                context.attributes["value"]);
            }
            break;

        case ETranslate: {
                Float x = parseFloat(name, context.attributes["x"], 0);
                Float y = parseFloat(name, context.attributes["y"], 0);
                Float z = parseFloat(name, context.attributes["z"], 0);
                m_transform = Transform::translate(Vector(x, y, z)) * m_transform;
            }
            break;

        case ERotate: {
                Float x = parseFloat(name, context.attributes["x"], 0);
                Float y = parseFloat(name, context.attributes["y"], 0);
                Float z = parseFloat(name, context.attributes["z"], 0);
                Float angle = parseFloat(name, context.attributes["angle"]);
                m_transform = Transform::rotate(Vector(x, y, z), angle) * m_transform;
            }
            break;

        case ELookAt: {
                std::vector<std::string> tokens = tokenize(context.attributes["origin"], ", ");
                if (tokens.size() != 3)
                    XMLLog(EError, "<lookat>: invalid 'origin' argument");
                Point o(
                    parseFloat(name, tokens[0]),
                    parseFloat(name, tokens[1]),
                    parseFloat(name, tokens[2]));
                tokens = tokenize(context.attributes["target"], ", ");
                if (tokens.size() != 3)
                    XMLLog(EError, "<lookat>: invalid 'target' argument");
                Point t(
                    parseFloat(name, tokens[0]),
                    parseFloat(name, tokens[1]),
                    parseFloat(name, tokens[2]));
                Vector u(0.0f);
                tokens = tokenize(context.attributes["up"], ", ");
                if (tokens.size() == 3)
                    u = Vector(
                        parseFloat(name, tokens[0]),
                        parseFloat(name, tokens[1]),
                        parseFloat(name, tokens[2]));
                else if (tokens.size() == 0)
                    ;
                else
                    XMLLog(EError, "<lookat>: invalid 'up' argument");

                if (u.lengthSquared() == 0) {
                    /* If 'up' was not specified, use an arbitrary axis */
                    Vector unused;
                    coordinateSystem(normalize(t-o), u, unused);
                }

                m_transform =  Transform::lookAt(o, t, u) * m_transform;
            }
            break;

        case EScale: {
                bool hasXYZ =
                    context.attributes["x"] != "" ||
                    context.attributes["y"] != "" ||
                    context.attributes["z"] != "";
                bool hasValue =
                    context.attributes["value"] != "";
                Float x=0, y=0, z=0;

                if (hasXYZ && hasValue) {
                    XMLLog(EError, "<scale>: provided both xyz and value arguments!");
                } else if (hasXYZ) {
                    x = parseFloat(name, context.attributes["x"], 1);
                    y = parseFloat(name, context.attributes["y"], 1);
                    z = parseFloat(name, context.attributes["z"], 1);
                } else if (hasValue) {
                    x = y = z = parseFloat(name, context.attributes["value"]);
                } else {
                    XMLLog(EError, "<scale>: provided neither xyz nor value arguments!");
                }

                m_transform = Transform::scale(Vector(x, y, z)) * m_transform;
            }
            break;

        case EMatrix: {
                std::vector<std::string> tokens = tokenize(
                    context.attributes["value"], ", ");
                if (tokens.size() != 16)
                    XMLLog(EError, "Invalid matrix specified");
                int index = 0;
                Matrix4x4 mtx;

                for (int i=0; i<4; ++i)
                    for (int j=0; j<4; ++j)
                        mtx.m[i][j] = parseFloat(name, tokens[index++]);

                m_transform = Transform(mtx) * m_transform;
            }
            break;

        case EPoint: {
                Float x = parseFloat(name, context.attributes["x"]);
                Float y = parseFloat(name, context.attributes["y"]);
                Float z = parseFloat(name, context.attributes["z"]);

                context.parent->properties.setPoint(context.attributes["name"], Point(x, y, z));
            }
            break;

        case EVector: {
                Float x = parseFloat(name, context.attributes["x"]);
                Float y = parseFloat(name, context.attributes["y"]);
                Float z = parseFloat(name, context.attributes["z"]);

                context.parent->properties.setVector(context.attributes["name"], Vector(x, y, z));
            }
            break;

        case ERGB: {
                Spectrum::EConversionIntent intent = Spectrum::EReflectance;
                if (context.parent->tag == EEmitter)
                    intent = Spectrum::EIlluminant;

                if (context.attributes.find("intent") != context.attributes.end()) {
                    std::string intentString = boost::to_lower_copy(context.attributes["intent"]);
                    if (intentString == "reflectance")
                        intent = Spectrum::EReflectance;
                    else if (intentString == "illuminant")
                        intent = Spectrum::EIlluminant;
                    else
                        XMLLog(EError, "Invalid intent \"%s\", must be "
                            "\"reflectance\" or \"illuminant\"", intentString.c_str());
                }

                std::string valueStr = context.attributes["value"];
                std::vector<std::string> tokens = tokenize(valueStr, ", ");
                Float value[3];
                if (tokens.size() == 1 && tokens[0].length() == 7 && tokens[0][0] == '#') {
                    char *end_ptr = NULL;
                    /* Parse HTML-style hexadecimal colors */
                    int encoded = strtol(tokens[0].c_str()+1, &end_ptr, 16);
                    if (*end_ptr != '\0')
                        XMLLog(EError, "Invalid rgb value specified (in <%s>)", context.attributes["name"].c_str());
                    value[0] = ((encoded & 0xFF0000) >> 16) / 255.0f;
                    value[1] = ((encoded & 0x00FF00) >> 8) / 255.0f;
                    value[2] =  (encoded & 0x0000FF) / 255.0f;
                } else if (tokens.size() == 1) {
                    value[0] = value[1] = value[2] = parseFloat(name, tokens[0]);
                } else if (tokens.size() == 3) {
                    for (int i=0; i<3; i++)
                        value[i] = parseFloat(name, tokens[i]);
                } else {
                    value[0] = value[1] = value[2] = 0; // avoid warning
                    XMLLog(EError, "Invalid RGB value specified");
                }
                Spectrum specValue;
                specValue.fromLinearRGB(value[0], value[1], value[2], intent);
                context.parent->properties.setSpectrum(context.attributes["name"],
                    specValue);
            }
            break;

        case ESRGB: {
                std::string valueStr = context.attributes["value"];
                std::vector<std::string> tokens = tokenize(valueStr, ", ");
                Float value[3];
                if (tokens.size() == 1 && tokens[0].length() == 7 && tokens[0][0] == '#') {
                    char *end_ptr = NULL;
                    /* Parse HTML-style hexadecimal colors */
                    int encoded = strtol(tokens[0].c_str()+1, &end_ptr, 16);
                    if (*end_ptr != '\0')
                        XMLLog(EError, "Invalid sRGB value specified (in <%s>)", context.attributes["name"].c_str());
                    value[0] = ((encoded & 0xFF0000) >> 16) / 255.0f;
                    value[1] = ((encoded & 0x00FF00) >> 8) / 255.0f;
                    value[2] =  (encoded & 0x0000FF) / 255.0f;
                } else if (tokens.size() == 1) {
                    value[0] = value[1] = value[2] = parseFloat(name, tokens[0]);
                } else if (tokens.size() == 3) {
                    for (int i=0; i<3; i++)
                        value[i] = parseFloat(name, tokens[i]);
                } else {
                    value[0] = value[1] = value[2] = 0; // avoid warning
                    XMLLog(EError, "Invalid sRGB value specified");
                }
                Spectrum specValue;
                specValue.fromSRGB(value[0], value[1], value[2]);
                context.parent->properties.setSpectrum(context.attributes["name"],
                    specValue);
            }
            break;

        case EBlackBody: {
                std::string temperature = trim(context.attributes["temperature"]);
                if (temperature.length() > 0 && std::toupper(temperature[temperature.length()-1]) == 'K')
                    temperature = temperature.substr(0, temperature.length()-1);
                Float temperatureValue = parseFloat(name, temperature);
                Float scale = 1;
                if (context.attributes.find("scale") != context.attributes.end())
                    scale = parseFloat(name, context.attributes["scale"]);
                BlackBodySpectrum bb(temperatureValue);
                Spectrum discrete;
                discrete.fromContinuousSpectrum(bb);
                discrete.clampNegative();
                context.parent->properties.setSpectrum(context.attributes["name"], discrete * scale);
            }
            break;

        case ESpectrum: {
                bool hasValue = context.attributes.find("value") != context.attributes.end();
                bool hasFilename = context.attributes.find("filename") != context.attributes.end();
                bool hasIntent = context.attributes.find("intent") != context.attributes.end();

                if (hasValue == hasFilename) {
                    XMLLog(EError, "<spectrum>: please provide one of 'value' or 'filename'");
                } else if (hasFilename) {
                    if (hasIntent)
                        XMLLog(EError, "<spectrum>: 'intent' and 'filename' cannot be specified at the same time!");
                    FileResolver *resolver = Thread::getThread()->getFileResolver();
                    fs::path path = resolver->resolve(context.attributes["filename"]);
                    InterpolatedSpectrum interp(path);
                    interp.zeroExtend();
                    Spectrum discrete;
                    discrete.fromContinuousSpectrum(interp);
                    discrete.clampNegative();
                    context.parent->properties.setSpectrum(context.attributes["name"], discrete);
                } else if (hasValue) {
                    std::vector<std::string> tokens = tokenize(
                        context.attributes["value"], ", ");
                    Float value[SPECTRUM_SAMPLES];
                    if (tokens.size() == 1 && tokens[0].find(':') == std::string::npos) {
                        Spectrum::EConversionIntent intent = Spectrum::EReflectance;
                        if (context.parent->tag == EEmitter)
                            intent = Spectrum::EIlluminant;

                        if (hasIntent) {
                            std::string intentString = boost::to_lower_copy(context.attributes["intent"]);
                            if (intentString == "reflectance")
                                intent = Spectrum::EReflectance;
                            else if (intentString == "illuminant")
                                intent = Spectrum::EIlluminant;
                            else
                                XMLLog(EError, "Invalid intent \"%s\", must be "
                                    "\"reflectance\" or \"illuminant\"", intentString.c_str());
                        }
                        value[0] = parseFloat(name, tokens[0]);
                        Spectrum spec;
                        if (intent == Spectrum::EReflectance)
                            spec = Spectrum(value[0]);
                        else
                            spec = Spectrum::getD65() * value[0];
                        context.parent->properties.setSpectrum(context.attributes["name"], spec);
                    } else {
                        if (hasIntent)
                            XMLLog(EError, "<spectrum>: 'intent' can only be specified when given a single-valued argument.");
                        if (tokens[0].find(':') != std::string::npos) {
                            InterpolatedSpectrum interp(tokens.size());
                            /* Wavelength -> Value mapping */
                            for (size_t i=0; i<tokens.size(); i++) {
                                std::vector<std::string> tokens2 = tokenize(tokens[i], ":");
                                if (tokens2.size() != 2)
                                    XMLLog(EError, "Invalid spectrum->value mapping specified");
                                Float wavelength = parseFloat(name, tokens2[0]);
                                Float value = parseFloat(name, tokens2[1]);
                                interp.append(wavelength, value);
                            }
                            interp.zeroExtend();
                            Spectrum discrete;
                            discrete.fromContinuousSpectrum(interp);
                            discrete.clampNegative();
                            context.parent->properties.setSpectrum(context.attributes["name"],
                                discrete);
                        } else {
                            if (tokens.size() != SPECTRUM_SAMPLES)
                                XMLLog(EError, "Invalid spectrum value specified (length does not match the current spectral discretization!)");
                            for (int i=0; i<SPECTRUM_SAMPLES; i++)
                                value[i] = parseFloat(name, tokens[i]);
                            context.parent->properties.setSpectrum(context.attributes["name"],
                                Spectrum(value));
                        }
                    }
                }
            }
            break;

        case EAnimation: {
                m_animatedTransform->sortAndSimplify();
                context.parent->properties.setAnimatedTransform(
                    context.attributes["name"], m_animatedTransform);
                m_animatedTransform = NULL;
            }
            break;

        case ETransform: {
                if (!m_animatedTransform.get()) {
                    context.parent->properties.setTransform(
                        context.attributes["name"], m_transform);
                } else {
                    Float time = parseFloat("time", context.attributes["time"]);
                    m_animatedTransform->appendTransform(time, m_transform);
                }
            }
            break;

        case EAlias: {
                std::string id = context.attributes["id"], as = context.attributes["as"];
                if (m_namedObjects->find(id) == m_namedObjects->end())
                    XMLLog(EError, "Referenced object '%s' not found!", id.c_str());
                ConfigurableObject *obj = (*m_namedObjects)[id];
                if (m_namedObjects->find(as) != m_namedObjects->end())
                    XMLLog(EError, "Duplicate ID '%s' used in scene description!", id.c_str());
                obj->incRef();
                (*m_namedObjects)[as] = obj;
            }
            break;

        case EInclude: {
                SAXParser* parser = new SAXParser();
                FileResolver *resolver = Thread::getThread()->getFileResolver();
                fs::path schemaPath = resolver->resolveAbsolute("data/schema/scene.xsd");

                /* Check against the 'scene.xsd' XML Schema */
                parser->setDoSchema(true);
                parser->setValidationSchemaFullChecking(true);
                parser->setValidationScheme(SAXParser::Val_Always);
                parser->setExternalNoNamespaceSchemaLocation(schemaPath.c_str());

                /* Set the handler and start parsing */
                SceneHandler *handler = new SceneHandler(m_params, m_namedObjects, true);
                parser->setDoNamespaces(true);
                parser->setDocumentHandler(handler);
                parser->setErrorHandler(handler);
                fs::path path = resolver->resolve(context.attributes["filename"]);
                XMLLog(EInfo, "Parsing included file \"%s\" ..", path.filename().string().c_str());
                parser->parse(path.c_str());

                object = handler->getScene();
                delete parser;
                delete handler;
            }
            break;

        case EDefault: {
                if (m_params.find(context.attributes["name"]) == m_params.end())
                    m_params[context.attributes["name"]] = context.attributes["value"];
            }
            break;

        default: {
                if (tag.second == NULL)
                    XMLLog(EError, "Internal error: could not instantiate an object "
                        "corresponding to the tag '%s'", name.c_str());

                Properties &props = context.properties;

                /* Convenience hack: allow passing animated transforms to arbitrary shapes
                   and then internally rewrite this into a shape group + animated instance */
                if (tag.second == MTS_CLASS(Shape)
                    && props.hasProperty("toWorld")
                    && props.getType("toWorld") == Properties::EAnimatedTransform
                    && (props.getPluginName() != "instance" && props.getPluginName() != "disk")) {
                    /* (The 'disk' plugin also directly supports animated transformations, so
                        the instancing trick isn't required for it) */

                    ref<const AnimatedTransform> trafo = props.getAnimatedTransform("toWorld");
                    props.removeProperty("toWorld");

                    if (trafo->isStatic())
                        props.setTransform("toWorld", trafo->eval(0));

                    object = m_pluginManager->createObject(tag.second, props);

                    if (!trafo->isStatic()) {
                        object = m_pluginManager->createObject(tag.second, props);
                        /* If the object has children, append them */
                        for (std::vector<std::pair<std::string, ConfigurableObject *> >
                                ::iterator it = context.children.begin();
                                it != context.children.end(); ++it) {
                            if (it->second != NULL) {
                                object->addChild(it->first, it->second);
                                it->second->setParent(object);
                                it->second->decRef();
                            }
                        }
                        context.children.clear();

                        object->configure();

                        ref<Shape> shapeGroup = static_cast<Shape *> (
                            m_pluginManager->createObject(MTS_CLASS(Shape), Properties("shapegroup")));
                        shapeGroup->addChild(object);
                        shapeGroup->configure();

                        Properties instanceProps("instance");
                        instanceProps.setAnimatedTransform("toWorld", trafo);
                        object = m_pluginManager->createObject(instanceProps);
                        object->addChild(shapeGroup);

                    }
                } else {
                    try {
                        object = m_pluginManager->createObject(tag.second, props);
                    } catch (const std::exception &ex) {
                        XMLLog(EError, "Error while creating object: %s", ex.what());
                    }
                }
            }
            break;
    }

    if (object != NULL || name == "null") {
        std::string id = context.attributes["id"];
        std::string nodeName = context.attributes["name"];

        if (object) {
            /* If the object has a parent, add it to the parent's children list */
            if (context.parent != NULL) {
                object->incRef();
                context.parent->children.push_back(
                    std::pair<std::string, ConfigurableObject *>(nodeName, object));
            }

            /* If the object has children, append them */
            for (std::vector<std::pair<std::string, ConfigurableObject *> >
                    ::iterator it = context.children.begin();
                    it != context.children.end(); ++it) {
                if (it->second != NULL) {
                    object->addChild(it->first, it->second);
                    it->second->setParent(object);
                    it->second->decRef();
                }
            }

            /* Don't configure a scene object if it is from an included file */
            if (name != "include" && (!m_isIncludedFile || !object->getClass()->derivesFrom(MTS_CLASS(Scene))))
                object->configure();

            if (object->getClass()->derivesFrom(MTS_CLASS(Texture)))
                object = static_cast<Texture *>(object.get())->expand();
        }

        if (id != "" && name != "ref") {
            if (m_namedObjects->find(id) != m_namedObjects->end())
                XMLLog(EError, "Duplicate ID '%s' used in scene description!", id.c_str());
            (*m_namedObjects)[id] = object;
            if (object)
                object->incRef();
        }
    }

    /* Warn about unqueried properties */
    std::vector<std::string> unq = context.properties.getUnqueried();
    for (unsigned int i=0; i<unq.size(); ++i)
        XMLLog(EWarn, "Unqueried attribute \"%s\" in element \"%s\"", unq[i].c_str(), name.c_str());

    m_context.pop();
}

// -----------------------------------------------------------------------
//  Implementation of the SAX ErrorHandler interface
// -----------------------------------------------------------------------

void SceneHandler::warning(const SAXParseException& e) {
    SLog(EWarn, "Warning in file \"%s\" (line %i): %s",
        transcode(e.getSystemId()).c_str(), e.getLineNumber(),
        transcode(e.getMessage()).c_str());
}

void SceneHandler::error(const SAXParseException& e) {
    SLog(EError, "Error in file \"%s\" (line %i): %s",
        transcode(e.getSystemId()).c_str(), e.getLineNumber(),
        transcode(e.getMessage()).c_str());
}

void SceneHandler::fatalError(const SAXParseException& e) {
    SLog(EError, "Fatal error in file \"%s\" (line %i): %s",
        transcode(e.getSystemId()).c_str(), e.getLineNumber(),
        transcode(e.getMessage()).c_str());
}

// -----------------------------------------------------------------------

ref<Scene> SceneHandler::loadScene(const fs::path &filename, const ParameterMap &params) {
    /* Prepare for parsing scene descriptions */
    FileResolver *resolver = Thread::getThread()->getFileResolver();
    SAXParser* parser = new SAXParser();
    fs::path schemaPath = resolver->resolveAbsolute("data/schema/scene.xsd");
    SLog(EDebug, "Loading scene \"%s\" ..", filename.string().c_str());

    /* Check against the 'scene.xsd' XML Schema */
    parser->setDoSchema(true);
    parser->setValidationSchemaFullChecking(true);
    parser->setValidationScheme(SAXParser::Val_Always);
    parser->setExternalNoNamespaceSchemaLocation(schemaPath.c_str());

    SceneHandler *handler = new SceneHandler(params);
    parser->setDoNamespaces(true);
    parser->setDocumentHandler(handler);
    parser->setErrorHandler(handler);

    parser->parse(filename.c_str());
    ref<Scene> scene = handler->getScene();

    delete parser;
    delete handler;

    return scene;
}

ref<Scene> SceneHandler::loadSceneFromString(const std::string &content, const ParameterMap &params) {
    /* Prepare for parsing scene descriptions */
    FileResolver *resolver = Thread::getThread()->getFileResolver();
    SAXParser* parser = new SAXParser();
    fs::path schemaPath = resolver->resolveAbsolute("data/schema/scene.xsd");

    /* Check against the 'scene.xsd' XML Schema */
    parser->setDoSchema(true);
    parser->setValidationSchemaFullChecking(true);
    parser->setValidationScheme(SAXParser::Val_Always);
    parser->setExternalNoNamespaceSchemaLocation(schemaPath.c_str());

    SceneHandler *handler = new SceneHandler(params);
    parser->setDoNamespaces(true);
    parser->setDocumentHandler(handler);
    parser->setErrorHandler(handler);

    XMLCh *inputName = XMLString::transcode("<string input>");

    MemBufInputSource input((const XMLByte *) content.c_str(),
            content.length(), inputName);
    parser->parse(input);
    ref<Scene> scene = handler->getScene();
    XMLString::release(&inputName);

    delete parser;
    delete handler;

    return scene;
}


void SceneHandler::staticInitialization() {
    /* Initialize Xerces-C */
    try {
        XMLPlatformUtils::Initialize();
    } catch (const XMLException &toCatch) {
        SLog(EError, "Error during Xerces initialization: %s",
            XMLString::transcode(toCatch.getMessage()));
    }
}

void SceneHandler::staticShutdown() {
    XMLPlatformUtils::Terminate();
}

VersionException::~VersionException() throw () {}

MTS_NAMESPACE_END
