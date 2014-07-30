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

#if !defined(__XMLTREEMODEL_H)
#define __XMLTREEMODEL_H

#include "common.h"
#include <QtXml/QtXml>

/**
 * Simplistic tree item used in the XML widget -- used to
 * store (key, value) pairs and categories.
 */
class TreeItem {
public:
	TreeItem(const QString &name, const QString &readableName,
			const QVariant &data = QVariant(),
			const QVariant &defaultValue = QVariant(),
			TreeItem *parent = 0,
			int importance = 10);
	~TreeItem();

	inline const QString &getName() const { return m_itemName; }
	inline int getImportance() const { return m_importance; }
	inline void setValue(const QVariant &value) { m_itemValue = value; }

	void appendChild(TreeItem *child);
	void removeChild(TreeItem *child);
	TreeItem *child(int row);
	int childCount() const;
	int columnCount() const;
	QVariant data(int column) const;
	int row() const;
	TreeItem *parent();
	bool isCategory() const;
	const QString &toolTip() const;
	void setToolTip(const QString &str);
	bool setData(int column, const QVariant &value);
	void setProperties(const Properties &props);
	void setProperty(const std::string &name, const Properties &props);
	void putProperties(Properties &props) const;
private:
	QList<TreeItem*> m_childItems;
	QString m_itemName, m_readableName;
	QVariant m_itemValue, m_itemDefault;
	QString m_toolTip;
	TreeItem *m_parentItem;
	int m_importance;
};

/**
 * Allows to edit the properties associated with a ConfigurableObject
 * using a tree widget. Metadata describing the flavor/name/description
 * of property entries is fetched from a supplied XML tree
 */
class XMLTreeModel : public QAbstractItemModel {
	Q_OBJECT

public:
	XMLTreeModel(QDomElement docRoot, const QPalette &pal, QObject *parent = 0);
	~XMLTreeModel();

	TreeItem *registerClass(const QString &className, const QString &readableName);
	TreeItem *updateClass(TreeItem *prev, const QString &_className,
				const QString &readableName);
	void setProperties(TreeItem *item, const Properties &props);

	QVariant data(const QModelIndex &index, int role) const;
	Qt::ItemFlags flags(const QModelIndex &index) const;
	QVariant headerData(int section, Qt::Orientation orientation,
						int role = Qt::DisplayRole) const;
	QModelIndex index(int row, int column,
					  const QModelIndex &parent = QModelIndex()) const;
	QModelIndex parent(const QModelIndex &index) const;
	int rowCount(const QModelIndex &parent = QModelIndex()) const;
	int columnCount(const QModelIndex &parent = QModelIndex()) const;
	bool setData(const QModelIndex &index, const QVariant &value, int role);
protected:
	void populate(const QString &className, TreeItem *parent);
private:
	TreeItem *m_rootItem;
	QDomElement m_docRoot;
	QVariant m_unimportantColor;
};

#endif // __XMLTREEMODEL_H
