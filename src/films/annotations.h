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
#if !defined(__ANNOTATIONS_H)
#define __ANNOTATIONS_H

#include <mitsuba/hw/font.h>
#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

/**
 * This function implements a parser for the 'label[]' and 'metadata[]'
 * annotations supported by the ldrfilm and hdrfilm plugins
 */
void annotate(const Scene *scene, const Properties &properties,
		Bitmap *bitmap, Float renderTime, Float gamma) {
	/* Attach the custom annotations */
	Properties &metadata = bitmap->getMetadata();
	std::vector<std::string> keys = properties.getPropertyNames();
	ref<Font> font;

	for (size_t i=0; i<keys.size(); ++i) {
		std::string key = keys[i];
		key.erase(std::remove_if(key.begin(), key.end(), ::isspace), key.end());
		std::string lkey = boost::to_lower_copy(key);
		Point2i offset(0, 0);
		bool labelAnnotation = false;

		if (boost::starts_with(lkey, "metadata['") && boost::ends_with(lkey, "']")) {
			key = key.substr(10, key.length()-12);
		} else if (boost::starts_with(lkey, "label[") && boost::ends_with(lkey, "]")) {
			std::vector<std::string> args = tokenize(key.substr(6, key.length()-7), " ,");
			if (args.size() != 2)
				SLog(EError, "Label command '%s' has an invalid number of arguments!", key.c_str());
			char *end_ptr = NULL;
			offset.x = strtol(args[0].c_str(), &end_ptr, 10);
			if (*end_ptr != '\0')
				SLog(EError, "Label command '%s' has an invalid position argument!", key.c_str());
			offset.y = strtol(args[1].c_str(), &end_ptr, 10);
			if (*end_ptr != '\0')
				SLog(EError, "Label command '%s' has an invalid position argument!", key.c_str());
			labelAnnotation = true;

			if (font == NULL) {
				font = new Font(Font::EBitstreamVeraMono14);
				font->convert(bitmap->getPixelFormat(), bitmap->getComponentFormat(), gamma);
			}
		} else {
			continue;
		}

		Properties::EPropertyType type = properties.getType(keys[i]);
		if (type == Properties::EString) {
			std::string value = properties.getString(keys[i]);

			while (true) {
				char *strt;
				if (!(strt = strchr((char *) value.c_str(), '$')))
					break;

				char *par1, *par2;
				if (!(par1 = strchr(strt, '[')))
					break;
				if (!(par2 = strchr(par1, ']')))
					break;

				std::string propSource   = value.substr(strt-value.c_str()+1, par1-strt-1);
				std::string propKey = value.substr(par1-value.c_str()+1, par2-par1-1);
				propSource.erase(std::remove_if(propSource.begin(), propSource.end(), ::isspace), propSource.end());
				propKey.erase(std::remove_if(propKey.begin(), propKey.end(), ::isspace), propKey.end());

				if (!boost::starts_with(propKey, "'") || !boost::ends_with(propKey, "'"))
					SLog(EError, "Encountered invalid key '%s'", propKey.c_str());

				propKey = propKey.substr(1, propKey.length()-2);

				const ConfigurableObject *source = NULL;
				if (propSource == "scene")
					source = scene;
				else if (propSource == "film")
					source = scene->getFilm();
				else if (propSource == "sampler")
					source = scene->getSampler();
				else if (propSource == "sensor")
					source = scene->getSensor();
				else if (propSource == "integrator")
					source = scene->getIntegrator();
				else
					SLog(EError, "Unknown data source '%s' (must be film/sampler/sensor/integrator)", propSource.c_str());

				std::string replacement;

				if (source == scene) {
					if (propKey == "renderTime") {
						replacement = timeString(renderTime);
					} else if (propKey == "renderTimePrecise") {
						replacement = timeString(renderTime, true);
					} else if (propKey == "memUsage") {
						replacement = memString(getPrivateMemoryUsage());
					} else if (propKey == "memUsagePrecise") {
						replacement = memString(getPrivateMemoryUsage(), true);
					} else if (propKey == "coreCount") {
						replacement = formatString("%i", Scheduler::getInstance()->getCoreCount());
					} else if (propKey == "blockSize") {
						replacement = formatString("%i", scene->getBlockSize());
					} else if (propKey == "sourceFile") {
						replacement = scene->getSourceFile().string();
					} else if (propKey == "destFile") {
						replacement = scene->getDestinationFile().string();
					}
				}

				if (replacement.empty()) {
					if (propKey == "type")
						replacement = source->getProperties().getPluginName();
					else
						replacement = source->getProperties().getAsString(propKey);
				}

				value.replace(strt-value.c_str(), par2-strt+1, replacement);
			}

			if (labelAnnotation) {
				Vector2i size = font->getSize(value);
				bitmap->fillRect(offset-Vector2i(4, 4), size + Vector2i(8, 8), Spectrum(0.0f));
				font->drawText(bitmap, offset, value);
			} else {
				metadata.setString(key, value);
			}
		} else {
			if (labelAnnotation)
				SLog(EError, "Only string-valued label annotations are supported!");
			metadata.copyAttribute(properties, keys[i], key);
		}
	}
}

MTS_NAMESPACE_END

#endif /* __ANNOTATIONS_H */
