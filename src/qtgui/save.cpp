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

#include "save.h"
#include <QtXml>

static QDomElement findUniqueChild(QDomElement element, const char *tagName) {
	QDomElement result;
	QDomNode n = element.firstChild();

	bool found = false;
	while (!n.isNull()) {
		QDomElement e = n.toElement();
		if (!e.isNull()) {
			if (e.tagName() == tagName) {
				if (found)
					SLog(EError, "Unexpectedly found multiple tags named \"%s\"", tagName);
				found = true;
				result = e;
			}
		}
		n = n.nextSibling();
	}
	return result;
}

static QDomElement findProperty(QDomElement element, const QString &propertyName) {
	QDomElement result;
	QDomNode n = element.firstChild();

	bool found = false;
	while (!n.isNull()) {
		QDomElement e = n.toElement();
		if (!e.isNull()) {
			if (e.attribute("name") == propertyName) {
				if (found)
					SLog(EError, "Unexpectedly found multiple properties named \"%s\"", qPrintable(propertyName));
				found = true;
				result = e;
			}
		}
		n = n.nextSibling();
	}
	return result;
}


static void removeChildren(QDomElement element) {
	while (!element.lastChild().isNull())
		element.removeChild(element.lastChild());
}

static void setProperties(QDomDocument &doc, QDomElement &element, 
		const Properties &props) {
	element.setAttribute("type", props.getPluginName().c_str());

	std::vector<std::string> propertyNames;
	props.putPropertyNames(propertyNames);
	for (std::vector<std::string>::const_iterator it = propertyNames.begin();
		it != propertyNames.end(); ++it) {
		QDomElement property;
		switch (props.getType(*it)) {
			case Properties::EBoolean:
				property = doc.createElement("boolean");
				property.setAttribute("value", props.getBoolean(*it) ? "true" : "false");
				break;
			case Properties::EInteger:
				property = doc.createElement("integer");
				property.setAttribute("value", QString::number(props.getInteger(*it)));
				break;
			case Properties::EFloat:
				property = doc.createElement("float");
				property.setAttribute("value", QString::number(props.getFloat(*it)));
				break;
			case Properties::EString:
				property = doc.createElement("string");
				property.setAttribute("value", props.getString(*it).c_str());
				break;
			default:
				SLog(EError, "setProperties(): \"%s\": Unable to handle elements of type %i", 
					(*it).c_str(), props.getType(*it));
		}
		property.setAttribute("name", (*it).c_str());
		element.appendChild(property);
	}
}

void saveScene(QWidget *parent, SceneContext *ctx, const QString &targetFile) {
	QDomDocument doc("MitsubaScene");
	QFile file(ctx->fileName);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		QMessageBox::critical(parent, parent->tr("Unable to save"),
			parent->tr("Unable to save changes: could not open the original scene file <b>%1</b>. "
			"Was this file moved or deleted in the meantime?").arg(ctx->fileName), 
			QMessageBox::Ok);
		return;
	}
	if (!doc.setContent(&file)) {
		QMessageBox::critical(parent, parent->tr("Unable to save"),
			parent->tr("Unable to save changes: could not parse the original scene file <b>%1</b>. "
			"Was this file modified in the meantime?").arg(ctx->fileName), 
			QMessageBox::Ok);
		file.close();
		return;
	}
	file.close();

	QDomElement root = doc.documentElement();

	// ====================================================================
	//   Serialize the camera configuration
	// ====================================================================

	QDomElement camera = findUniqueChild(root, "camera");
	if (camera.isNull()) {
		camera = doc.createElement("camera");
		camera.setAttribute("type", "perspective");
		root.insertAfter(camera, QDomNode());
	}
	const PerspectiveCamera *sceneCamera = static_cast<const PerspectiveCamera *>(ctx->scene->getCamera());
	
	QDomElement fovProperty = findProperty(camera, "fov");
	if (fovProperty.isNull()) {
		fovProperty = doc.createElement("float");
		fovProperty.setAttribute("name", "fov");
		camera.insertBefore(fovProperty, QDomNode());
	}
	fovProperty.setAttribute("value", QString::number(sceneCamera->getFov()));

	QDomElement cameraTransform = findUniqueChild(camera, "transform");
	if (cameraTransform.isNull()) {
		cameraTransform = doc.createElement("transform");
		cameraTransform.setAttribute("name", "toWorld");
		camera.insertBefore(cameraTransform, QDomNode());
	}

	removeChildren(cameraTransform);
	QDomElement lookAt = doc.createElement("lookAt");
	cameraTransform.insertAfter(lookAt, QDomNode());

	Vector direction = sceneCamera->getInverseViewTransform()(Vector(0,0,1)),
		   u = ctx->up;
	Point t, p = sceneCamera->getInverseViewTransform()(Point(0,0,0));

	if (sceneCamera->getViewTransform().det3x3() > 0) {
		QDomElement scale = doc.createElement("scale");
		scale.setAttribute("x", "-1");
		cameraTransform.insertBefore(scale, lookAt);
	}
	t = p + direction;

	lookAt.setAttribute("ox", QString::number(p.x));
	lookAt.setAttribute("oy", QString::number(p.y));
	lookAt.setAttribute("oz", QString::number(p.z));
	lookAt.setAttribute("tx", QString::number(t.x));
	lookAt.setAttribute("ty", QString::number(t.y));
	lookAt.setAttribute("tz", QString::number(t.z));
	lookAt.setAttribute("ux", QString::number(u.x));
	lookAt.setAttribute("uy", QString::number(u.y));
	lookAt.setAttribute("uz", QString::number(u.z));

	// ====================================================================
	//   Serialize the sampler configuration
	// ====================================================================

	QDomElement newSampler, sampler = findUniqueChild(camera, "sampler");
	if (sampler.isNull()) {
		newSampler = doc.createElement("sampler");
		camera.appendChild(newSampler);
	} else {
		newSampler = doc.createElement("sampler");
		camera.insertAfter(newSampler, sampler);
		camera.removeChild(sampler);
	}

	setProperties(doc, newSampler, ctx->scene->getSampler()->getProperties());

	// ====================================================================
	//   Serialize the film configuration
	// ====================================================================

	QDomElement film = findUniqueChild(camera, "film");
	if (film.isNull()) {
		film = doc.createElement("film");
		film.setAttribute("type", "exrfilm");
		camera.insertAfter(film, QDomNode());
	}

	QDomElement widthProperty = findProperty(film, "width");
	QDomElement heightProperty = findProperty(film, "height");

	if (widthProperty.isNull()) {
		widthProperty = doc.createElement("integer");
		widthProperty.setAttribute("name", "width");
		film.insertBefore(widthProperty, QDomNode());
	}

	if (heightProperty.isNull()) {
		heightProperty = doc.createElement("integer");
		heightProperty.setAttribute("name", "height");
		film.insertBefore(heightProperty, QDomNode());
	}

	if (film.attribute("type") == "pngfilm") {
		/* Set tonemapping attributes */
		QDomElement method = findProperty(film, "toneMappingMethod");
		QDomElement reinhardBurn = findProperty(film, "reinhardBurn");
		QDomElement reinhardKey = findProperty(film, "reinhardKey");
		QDomElement exposure = findProperty(film, "exposure");
		QDomElement gamma = findProperty(film, "gamma");

		if (method.isNull()) {
			method = doc.createElement("string");
			method.setAttribute("name", "toneMappingMethod");
			film.insertBefore(method, QDomNode());
		}
		method.setAttribute("value", ctx->toneMappingMethod == EReinhard ? "reinhard" : "gamma");
		
		if (gamma.isNull()) {
			gamma = doc.createElement("float");
			gamma.setAttribute("name", "gamma");
			film.insertBefore(gamma, QDomNode());
		}
		gamma.setAttribute("value", QString::number(ctx->srgb ? (Float) -1 : ctx->gamma));

		if (ctx->toneMappingMethod == EGamma) {
			if (exposure.isNull()) {
				exposure = doc.createElement("float");
				exposure.setAttribute("name", "exposure");
				film.insertBefore(exposure, QDomNode());
			}
			exposure.setAttribute("value", QString::number(ctx->exposure));
			if (!reinhardKey.isNull())
				film.removeChild(reinhardKey);
			if (!reinhardBurn.isNull())
				film.removeChild(reinhardBurn);
		} else {
			if (reinhardKey.isNull()) {
				reinhardKey = doc.createElement("float");
				reinhardKey.setAttribute("name", "reinhardKey");
				film.insertBefore(reinhardKey, QDomNode());
			}
			if (reinhardBurn.isNull()) {
				reinhardBurn = doc.createElement("float");
				reinhardBurn.setAttribute("name", "reinhardBurn");
				film.insertBefore(reinhardBurn, QDomNode());
			}
			reinhardKey.setAttribute("value", QString::number(ctx->reinhardKey));
			reinhardBurn.setAttribute("value", QString::number(ctx->reinhardBurn));
			if (!exposure.isNull())
				film.removeChild(exposure);
		}
	}

	Vector2i filmSize = sceneCamera->getFilm()->getSize();
	widthProperty.setAttribute("value", QString::number(filmSize.x));
	heightProperty.setAttribute("value", QString::number(filmSize.y));

	// ====================================================================
	//   Serialize the reconstruction filter configuration
	// ====================================================================

	QDomElement newRFilter, rfilter = findUniqueChild(film, "rfilter");
	if (rfilter.isNull()) {
		newRFilter = doc.createElement("rfilter");
		film.appendChild(rfilter);
	} else {
		newRFilter = doc.createElement("rfilter");
		film.insertAfter(newRFilter, rfilter);
		film.removeChild(rfilter);
	}

	setProperties(doc, newRFilter, ctx->scene->getCamera()->
		getFilm()->getReconstructionFilter()->getProperties());

	// ====================================================================
	//   Serialize the integrator configuration
	// ====================================================================

	QDomElement oldIntegratorNode = findUniqueChild(root, "integrator");
	QDomElement newIntegratorNode = doc.createElement("integrator");
	if (oldIntegratorNode.isNull()) {
		film.appendChild(oldIntegratorNode);
	} else {
		film.insertAfter(newIntegratorNode, oldIntegratorNode);
		film.removeChild(oldIntegratorNode);
	}

	const Integrator *integrator = ctx->scene->getIntegrator();
	setProperties(doc, newIntegratorNode, integrator->getProperties());
	QDomElement currentIntegratorNode = newIntegratorNode;

	while (integrator->getSubIntegrator() != NULL) {
		integrator = integrator->getSubIntegrator();
		QDomElement childIntegratorNode = doc.createElement("integrator");
		setProperties(doc, childIntegratorNode, integrator->getProperties());
		currentIntegratorNode.appendChild(childIntegratorNode);
		currentIntegratorNode = childIntegratorNode;
	}

	root.insertAfter(newIntegratorNode, oldIntegratorNode);
	root.removeChild(oldIntegratorNode);

	file.setFileName(targetFile);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text)) {
		QMessageBox::critical(parent, parent->tr("Unable to save"),
			parent->tr("Unable to save changes: could not open the destination file <b>%1</b>.").arg(targetFile),
			QMessageBox::Ok);
		return;
	}

	/* Clean up the XML output generated by Qt so that it
	   is more suitable for human consumption ..
	   Beware: the code below is tailored to Qt's
	   output and won't work on arbitrary XML files */
	QString textContent = doc.toString();
	QTextStream input(&textContent);
	QTextStream output(&file);
	QRegExp 
		filenameRegExp("filename=\"[^\"]*\""),
		nameRegExp("name=\"[^\"]*\""),
		tagRegExp("^\\s*<([a-zA-Z]+) "),
		leadingSpaces("^ *"),
		closeTag("^\\s*</");
	bool inComment = false, hasContents = false;

	while (!input.atEnd()) {
		QString line = input.readLine();
		bool startComment = line.contains("<!--");
		bool endComment = line.contains("-->");

		if (startComment)
			inComment = true;
		if (endComment)
			inComment = false;

		if (inComment || endComment) {
			if (startComment) {
				/* Turn leading spaces into tabs */
				if (leadingSpaces.indexIn(line) == 0)
					line = QString('\t').repeated(leadingSpaces.matchedLength()) +
						line.mid(leadingSpaces.matchedLength());
			}
			output << line << endl;
			continue;
		}

		/* Make sure that the name attribute is always the first one */
		int tagMatch = tagRegExp.indexIn(line),
			tagLength = tagRegExp.matchedLength();
		int nameMatch = nameRegExp.indexIn(line),
			filenameMatch = filenameRegExp.indexIn(line),
			nameLength = nameRegExp.matchedLength();
		if (tagMatch != -1 && nameMatch != -1 && filenameMatch == -1) {
			line = line.left(tagLength) + line.mid(nameMatch, nameLength) + " "
				+ line.mid(tagMatch+tagLength, nameMatch-(tagMatch+tagLength))
				+ line.mid(nameMatch+nameLength);
		}

		/* Add an extra newline if this is an object tag, and if there
		   have been siblings before it */
		if (tagMatch != -1) {
			const QString &el = tagRegExp.cap(1);
			bool isObject = true;

			isObject &= (el != "string");
			isObject &= (el != "integer");
			isObject &= (el != "float");
			isObject &= (el != "boolean");
			isObject &= (el != "vector");
			isObject &= (el != "point");
			isObject &= (el != "transform");
			isObject &= (el != "spectrum");
			isObject &= (el != "rgb");
			isObject &= (el != "scale");
			isObject &= (el != "translate");
			isObject &= (el != "rotate");
			isObject &= (el != "lookAt");
			isObject &= (el != "matrix");

			if (isObject && hasContents) {
				output << endl;
				hasContents = false;
			}

			if (!isObject)
				hasContents = true;
		}

		/* Turn leading spaces into tabs */
		if (leadingSpaces.indexIn(line) == 0)
			line = QString('\t').repeated(leadingSpaces.matchedLength()) +
				line.mid(leadingSpaces.matchedLength());

		/* Remove ugly spaces */
		if (line.endsWith("  />")) {
			line = line.left(line.size()-4) + QString("/>");
			hasContents = true;
		} else if (line.endsWith(" />")) {
			line = line.left(line.size()-3) + QString("/>");
			hasContents = true;
		} else if (line.endsWith(" >")) {
			line = line.left(line.size()-2) + QString(">");
		}

		if (closeTag.indexIn(line) == 0)
			hasContents = true;

		output << line << endl;
	}
	file.close();
}

