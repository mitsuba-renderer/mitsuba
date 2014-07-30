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

#include "save.h"

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

static QList<QDomElement> findAllChildren(QDomElement element, const char *tagName) {
	QList<QDomElement> result;
	QDomNode n = element.firstChild();

	while (!n.isNull()) {
		QDomElement e = n.toElement();
		if (!e.isNull() && e.tagName() == tagName)
			result.append(e);
		n = n.nextSibling();
	}
	return result;
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
			case Properties::EAnimatedTransform: {
					const AnimatedTransform *trafo = props.getAnimatedTransform(*it);

					std::set<Float> times;
					trafo->collectKeyframes(times);

					property = doc.createElement("animation");

					for (std::set<Float>::iterator it2 = times.begin(); it2 != times.end(); ++it2) {
						const Matrix4x4 &matrix = trafo->eval(*it2).getMatrix();
						QDomElement trafoTag = doc.createElement("transform");
						QDomElement matrixTag = doc.createElement("matrix");

						QString value;
						for (int i=0; i<4; ++i)
							for (int j=0; j<4; ++j)
								value += QString("%1 ").arg(matrix(i, j));

						matrixTag.setAttribute("value", value);
						trafoTag.setAttribute("time", *it2);
						trafoTag.appendChild(matrixTag);
						property.appendChild(trafoTag);
					}
				}
				break;
			case Properties::ETransform: {
					/* Captures the subset of transformations that are used by
					   Mitsuba's perspective and orthographic camera classes */
					property = doc.createElement("transform");
					Transform trafo = props.getTransform(*it);
					if (trafo.hasScale()) {
						QDomElement scale = doc.createElement("scale");
						Float valueX = trafo(Vector(1, 0, 0)).length();
						Float valueY = trafo(Vector(0, 1, 0)).length();
						Float valueZ = trafo(Vector(0, 0, 1)).length();
						if (std::abs(1-valueX) < 1e-3f)
							valueX = 1.0f;
						if (std::abs(1-valueY) < 1e-3f)
							valueY = 1.0f;
						if (std::abs(1-valueZ) < 1e-3f)
							valueZ = 1.0f;
						scale.setAttribute("x", QString("%1").arg(valueX));
						scale.setAttribute("y", QString("%1").arg(valueY));
						scale.setAttribute("z", QString("%1").arg(valueZ));
						property.appendChild(scale);
						trafo = trafo * Transform::scale(Vector(1.0f/valueX, 1.0f/valueY, 1.0f/valueZ));
					}
					QDomElement lookAt = doc.createElement("lookat");
					property.appendChild(lookAt);
					if (trafo.det3x3() < 0) {
						QDomElement scale = doc.createElement("scale");
						scale.setAttribute("x", "-1");
						property.insertBefore(scale, lookAt);
					}
					Point p = trafo(Point(0,0,0));
					Point t = p + trafo(Vector(0,0,1));
					Vector u = trafo(Vector(0,1,0));

					lookAt.setAttribute("origin", QString("%1, %2, %3").arg(p.x).arg(p.y).arg(p.z));
					lookAt.setAttribute("up",     QString("%1, %2, %3").arg(u.x).arg(u.y).arg(u.z));
					lookAt.setAttribute("target", QString("%1, %2, %3").arg(t.x).arg(t.y).arg(t.z));
				}
				break;
			default:
				SLog(EError, "setProperties(): \"%s\": Unable to handle elements of type %i",
					(*it).c_str(), props.getType(*it));
		}
		property.setAttribute("name", (*it).c_str());
		element.appendChild(property);
	}
}

void processSubIntegrators(QDomDocument &doc, const Integrator *integrator, QDomElement integratorNode) {
	int idx = 0;
	while (integrator->getSubIntegrator(idx) != NULL) {
		const Integrator *childIntegrator = integrator->getSubIntegrator(idx);
		QDomElement childIntegratorNode = doc.createElement("integrator");
		setProperties(doc, childIntegratorNode, childIntegrator->getProperties());
		integratorNode.appendChild(childIntegratorNode);
		processSubIntegrators(doc, childIntegrator, childIntegratorNode);
		idx++;
	}
}

void saveScene(QWidget *parent, SceneContext *ctx, const QString &targetFile) {
	QDomElement root = ctx->doc.documentElement();

	// ====================================================================
	//   Serialize the sensor configuration
	// ====================================================================

	QList<QDomElement> oldSensors = findAllChildren(root, "sensor");

	const ref_vector<Sensor> sensors = ctx->scene->getSensors();
	ref_vector<Sensor>::const_iterator it = std::find(sensors.begin(),
			sensors.end(), ctx->scene->getSensor());
	if (it == sensors.end())
		SLog(EError, "Number of sensors did not match between loaded scene and XML file!");

	QDomElement sensor, oldSensor;

	if (oldSensors.size() == 0)
		; // ok -- scene did not contain a sensor before
	else if ((size_t) oldSensors.size() != ctx->scene->getSensors().size())
		SLog(EError, "Number of sensors did not match between loaded scene and XML file!");
	else
		oldSensor = oldSensors[it-sensors.begin()];

	if (oldSensor.isNull()) {
		sensor = ctx->doc.createElement("sensor");
		root.insertAfter(sensor, QDomNode());
	} else {
		sensor = ctx->doc.createElement("sensor");
		root.insertAfter(sensor, oldSensor);
		root.removeChild(oldSensor);
	}

	setProperties(ctx->doc, sensor, ctx->scene->getSensor()->getProperties());

	// ====================================================================
	//   Serialize the sampler configuration
	// ====================================================================

	QDomElement sampler = ctx->doc.createElement("sampler");
	sensor.appendChild(sampler);

	setProperties(ctx->doc, sampler,
		ctx->scene->getSampler()->getProperties());

	// ====================================================================
	//   Serialize the film configuration
	// ====================================================================

	QDomElement film = ctx->doc.createElement("film");
	sensor.appendChild(film);

	Properties filmProps(ctx->scene->getFilm()->getProperties());
	if (filmProps.getPluginName() == "ldrfilm") {
		/* Also export the tonemapper settings */
		if (ctx->toneMappingMethod == EGamma) {
			filmProps.setString("tonemapMethod", "gamma", false);
			filmProps.setFloat("exposure", ctx->exposure, false);
			filmProps.removeProperty("key");
			filmProps.removeProperty("burn");
		} else {
			filmProps.setString("tonemapMethod", "reinhard", false);
			filmProps.setFloat("key", ctx->reinhardKey, false);
			filmProps.setFloat("burn", (ctx->reinhardBurn + 10) / 20.0f, false);
			filmProps.removeProperty("exposure");
		}
		filmProps.setFloat("gamma", ctx->srgb ? -1 : ctx->gamma, false);
	}

	setProperties(ctx->doc, film, filmProps);

	// ====================================================================
	//   Serialize the reconstruction filter configuration
	// ====================================================================

	QDomElement rfilter = ctx->doc.createElement("rfilter");
	film.appendChild(rfilter);

	setProperties(ctx->doc, rfilter,
		ctx->scene->getFilm()->getReconstructionFilter()->getProperties());

	// ====================================================================
	//   Serialize medium references of the sensor
	// ====================================================================

	QList<QDomElement> oldSensorReferences = findAllChildren(oldSensor, "ref");
	oldSensorReferences.append(findAllChildren(oldSensor, "medium"));

	for (int i=0; i<oldSensorReferences.size(); ++i)
		sensor.appendChild(ctx->doc.importNode(oldSensorReferences[i], true));

	// ====================================================================
	//   Serialize the integrator configuration
	// ====================================================================

	QDomElement oldIntegratorNode = findUniqueChild(root, "integrator");
	QDomElement newIntegratorNode = ctx->doc.createElement("integrator");

	const Integrator *integrator = ctx->scene->getIntegrator();
	setProperties(ctx->doc, newIntegratorNode, integrator->getProperties());
	processSubIntegrators(ctx->doc, integrator, newIntegratorNode);

	root.insertBefore(newIntegratorNode, oldIntegratorNode);
	root.removeChild(oldIntegratorNode);

	QFile file;
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
	QString textContent = ctx->doc.toString();
	QTextStream input(&textContent);
	QTextStream output(&file);
	cleanupXML(input, output);
	file.close();
}

void cleanupXML(QTextStream &input, QTextStream &output) {
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
			QString a = line.mid(tagLength, nameMatch-tagLength).trimmed(),
				b = line.mid(nameMatch+nameLength).trimmed();
			line = line.left(tagLength) + line.mid(nameMatch, nameLength);
			if (a.length() > 0)
				line += " " + a;
			if (b.length() > 0)
				line += " " + b;
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
		} else if (line.endsWith("/>")) {
			hasContents = true;
		} else if (line.endsWith(" >")) {
			line = line.left(line.size()-2) + QString(">");
		} else if (line.endsWith("?>")) {
			hasContents = true;
		}

		if (closeTag.indexIn(line) == 0)
			hasContents = true;

		output << line << endl;
	}
}
