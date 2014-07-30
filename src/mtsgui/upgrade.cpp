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

#include "upgrade.h"
#include "save.h"
#include <mitsuba/core/fresolver.h>
#include <boost/algorithm/string.hpp>
#include <QtXmlPatterns/QtXmlPatterns>

struct VersionComparator {
	inline bool operator()(const std::pair<Version, fs::path> &s1,
		const std::pair<Version, fs::path> &s2) const {
		return s1.first < s2.first;
	}
};

class XSLTMessageHandler : public QAbstractMessageHandler {
public:
	XSLTMessageHandler() : m_fatalError(false) { }

	void handleMessage(QtMsgType type, const QString &_descr, const QUrl &id, const QSourceLocation &loc) {
		QString descr(_descr);
		int paragraphStart = descr.indexOf("<p>");
		int paragraphEnd = descr.indexOf("</p>");
		if (paragraphStart != -1 && paragraphEnd != -1)
			descr = descr.mid(paragraphStart+3, paragraphEnd-paragraphStart-3);

		if (loc.isNull())
			SLog(EWarn, "%s", qPrintable(descr));
		else
			SLog(EWarn, "%s (line %i, column %i, url=%s)", qPrintable(descr),
				loc.line(), loc.column(), qPrintable(id.toString()));

		if (type == QtFatalMsg)
			m_fatalError = true;
	}

	inline bool fatalError() const { return m_fatalError; }
private:
	bool m_fatalError;
};

UpgradeManager::UpgradeManager(const FileResolver *resolver) : m_resolver(resolver){
	fs::path transformationPath =
		resolver->resolveAbsolute("data/schema/scene.xsd").parent_path();

	fs::directory_iterator it(transformationPath), end;
	SLog(EInfo, "Searching for transformations..");

	for (; it != end; ++it) {
		fs::path file = *it;
		std::string extension = file.extension().string(),
			filename = file.filename().string();
		if (boost::to_lower_copy(extension) != ".xsl" ||
           !boost::starts_with(filename, "upgrade_"))
			continue;
		Version version(filename.substr(8, filename.length()-12));
		m_transformations.push_back(std::make_pair(version, file));
	}

	std::sort(m_transformations.begin(), m_transformations.end(), VersionComparator());

	for (size_t i=0; i<m_transformations.size(); ++i)
		SLog(EInfo, "  - registered transformation \"%s\", which updates to version %s",
			m_transformations[i].second.filename().string().c_str(),
			m_transformations[i].first.toString().c_str());
}

void UpgradeManager::performUpgrade(const QString &filename, const Version &version) {
	QString backupFilename = filename;
	if (backupFilename.endsWith(".xml", Qt::CaseInsensitive))
		backupFilename.replace(backupFilename.length()-4, 4, ".bak");
	SLog(EInfo, "Saving a backup copy of \"%s\" as \"%s\"",
		QStringToUTF8(filename), QStringToUTF8(backupFilename));

	QFile file(filename), backupFile(backupFilename);
	if (backupFile.exists())
		backupFile.remove();
	if (!file.copy(backupFilename))
		SLog(EError, "Could not create a backup copy -- "
			"stopping the upgrade operation!");

	if (!file.open(QIODevice::ReadOnly))
		SLog(EError, "Unable to open \"%s\" with read access -- stopping "
			"the upgrade operation.", QStringToUTF8(filename));
	QByteArray inputArray = file.readAll(), outputArray;
	file.close();

	QBuffer inputBuffer(&inputArray);
	QBuffer outputBuffer(&outputArray);

	XSLTMessageHandler handler;

	int nTransformations=0;
	for (size_t i=0; i<m_transformations.size(); ++i) {
		if (m_transformations[i].first <= version)
			continue;
		inputBuffer.open(QIODevice::ReadOnly);
		outputBuffer.open(QIODevice::WriteOnly | QIODevice::Truncate);

		SLog(EInfo, "Applying transformation \"%s\" ..",
			m_transformations[i].second.filename().string().c_str());
		QString trafoFilename = fromFsPath(m_transformations[i].second);
		QFile trafoFile(trafoFilename);
		if (!trafoFile.open(QIODevice::ReadOnly | QIODevice::Text))
			SLog(EError, "Unable to open the stylesheet \"%s\" -- stopping "
				"the upgrade operation.", QStringToUTF8(trafoFilename));

		QXmlQuery query(QXmlQuery::XSLT20);
		query.setMessageHandler(&handler);
		query.setFocus(&inputBuffer);
		query.setQuery(&trafoFile);
		if (!query.isValid())
			SLog(EError, "Unable to parse the stylesheet \"%s\" -- stopping "
				"the upgrade operation.", QStringToUTF8(trafoFilename));

		SLog(EInfo, "Transformation is ready, running it now..");
		QXmlFormatter formatter(query, &outputBuffer);
		formatter.setIndentationDepth(1);
		query.evaluateTo(&formatter);

		if (handler.fatalError())
			SLog(EError, "A fatal error was encountered -- stopping "
				"the upgrade operation.");

		/* Swap the outputs */
		inputBuffer.close();
		outputBuffer.close();
		inputArray = outputArray;
		++nTransformations;
	}

	if (nTransformations == 0) {
		SLog(EError, "Strange -- no transformations were applied? "
			"Stopping the upgrade operation.");
	}

	/* Done, write back to disk */
	SLog(EInfo, "Successfully applied %i transformations, writing the result "
		"to \"%s\" ..", nTransformations, QStringToUTF8(filename));

	if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
		SLog(EError, "Unable to open \"%s\" with write access -- stopping "
			"the upgrade operation.", QStringToUTF8(filename));
	}

	int line, column;
	QString errorMsg;
	QDomDocument doc;
	if (!doc.setContent(inputArray, &errorMsg, &line, &column))
		SLog(EError, "Unable to parse file: error %s at line %i, colum %i",
			QStringToUTF8(errorMsg), line, column);
	QDomElement root = doc.documentElement();

	/* Search for include nodes */
	QDomNode n = root.firstChild();
	while (!n.isNull()) {
		QDomElement e = n.toElement();
		if (!e.isNull()) {
			if (e.tagName() == "include") {
				fs::path path = m_resolver->resolve(toFsPath(e.attribute("filename")));
				SLog(EInfo, "Recursively upgrading include file \"%s\" ..",
						path.string().c_str());
				performUpgrade(fromFsPath(path), version);
			}
		}
		n = n.nextSibling();
	}

	QTextStream input(inputArray);
	QTextStream output(&file);
	output << "<?xml version=\"1.0\" encoding=\"utf-8\"?>"
		<< endl << endl;
	cleanupXML(input, output);
	file.close();

	QDomDocument tooltipDoc;

}
