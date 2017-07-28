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

#include "xmltreemodel.h"

// ====================================================================
//    TreeItem implementation
// ====================================================================

TreeItem::TreeItem(const QString &name, const QString &readableName,
                   const QVariant &data, const QVariant &defaultValue,
                   TreeItem *parent, int importance) {
    m_parentItem = parent;
    m_itemName = name;
    m_readableName = readableName;
    m_itemValue = data;
    m_itemDefault = defaultValue;
    m_importance = importance;
}

TreeItem::~TreeItem() {
    qDeleteAll(m_childItems);
}

void TreeItem::appendChild(TreeItem *item) {
    m_childItems.append(item);
}

void TreeItem::removeChild(TreeItem *item) {
    m_childItems.removeOne(item);
}

TreeItem *TreeItem::child(int row) {
    return m_childItems.value(row);
}

int TreeItem::childCount() const {
    return m_childItems.count();
}

int TreeItem::columnCount() const {
    if (isCategory())
        return 1;
    else
        return 2;
}

QVariant TreeItem::data(int column) const {
    if (column == 0)
        return QVariant(m_readableName);
    else if (column == 1 && !isCategory())
        return m_itemValue;
    else return QVariant();
}

int TreeItem::row() const {
    if (m_parentItem)
         return m_parentItem->m_childItems.indexOf(const_cast<TreeItem*>(this));
     return 0;
}

TreeItem *TreeItem::parent() {
    return m_parentItem;
}

void TreeItem::setToolTip(const QString &str) {
    m_toolTip = str;
}

const QString &TreeItem::toolTip() const {
    return m_toolTip;
}

bool TreeItem::isCategory() const {
    return m_itemValue.isNull();
}

bool TreeItem::setData(int column, const QVariant &value) {
    if (isCategory() || column != 1)
        return false;

    m_itemValue = value;
    return true;
}

void TreeItem::setProperty(const std::string &name, const Properties &props) {
    QVariant data;
    switch (props.getType(name)) {
        case Properties::EBoolean:
            data = QVariant(props.getBoolean(name));
            break;
        case Properties::EInteger:
            data = QVariant(props.getInteger(name));
            break;
        case Properties::EFloat:
            data = QVariant((double) props.getFloat(name));
            break;
        case Properties::EString:
            data = QVariant(props.getString(name).c_str());
            break;
        default:
            SLog(EError, "TreeItem::getProperties(): \"%s\": Unable to handle elements of type %i",
                name.c_str(), props.getType(name));
    }

    bool found = false;
    for (int i=0; i<m_childItems.size(); ++i) {
        if (m_childItems[i]->m_itemName == name.c_str()) {
            m_childItems[i]->setData(1, data);
            found = true;
            break;
        }
    }
    if (!found)
        SLog(EWarn, "TreeItem::getProperties(): \"%s\": Unable to find element in tree!", name.c_str());
}

void TreeItem::setProperties(const Properties &props) {
    std::vector<std::string> propertyNames;
    props.putPropertyNames(propertyNames);

    for (std::vector<std::string>::const_iterator it = propertyNames.begin();
        it != propertyNames.end(); ++it)
        setProperty(*it, props);
}

void TreeItem::putProperties(Properties &props) const {
    for (int i=0; i<m_childItems.size(); ++i) {
        const TreeItem *item = m_childItems[i];
        const std::string name(item->m_itemName.toStdString());
        const QVariant &value(item->m_itemValue);

        if (value == item->m_itemDefault)
            continue;

        switch (value.type()) {
            case QVariant::Int:
                props.setInteger(name, value.toInt());
                break;
            case QVariant::Double:
                props.setFloat(name, (Float) value.toDouble());
                break;
            case QVariant::String:
                props.setString(name, value.toString().toStdString());
                break;
            case QVariant::Bool:
                props.setBoolean(name, value.toBool());
                break;
            default:
                SLog(EError, "TreeItem::putProperties(): \"%s\": Unable to handle elements of type %i",
                    name.c_str(), value.type());
        }
    }
}

// ====================================================================
//    XMLTreeModel implementation
// ====================================================================


XMLTreeModel::XMLTreeModel(QDomElement docRoot, const QPalette &palette, QObject *parent)
    : QAbstractItemModel(parent), m_docRoot(docRoot) {
    m_rootItem = new TreeItem("Root", "Root");
    m_unimportantColor = QVariant(palette.color(QPalette::Disabled, QPalette::Text));
}

TreeItem *XMLTreeModel::registerClass(const QString &_className, const QString &readableName) {
    QString className(_className);

    if (className == "")
        return NULL;

    TreeItem *parent = new TreeItem(className, readableName,
        QVariant(), QVariant(), m_rootItem);
    populate(className, parent);

    if (parent->childCount() == 0) {
        delete parent;
        return NULL;
    }

    beginInsertRows(QModelIndex(), m_rootItem->childCount(), m_rootItem->childCount());
    m_rootItem->appendChild(parent);
    endInsertRows();

    emit layoutChanged();
    return parent;
}

TreeItem *XMLTreeModel::updateClass(TreeItem *prev, const QString &className, const QString &readableName) {
    if (prev != NULL) {
        beginRemoveRows(QModelIndex(), prev->row(), prev->row());
        m_rootItem->removeChild(prev);
        endRemoveRows();
        if (className == prev->getName()) {
            /* Just re-insert */
            m_rootItem->removeChild(prev);
            beginInsertRows(QModelIndex(), m_rootItem->childCount(), m_rootItem->childCount());
            m_rootItem->appendChild(prev);
            endInsertRows();
            emit layoutChanged();
            return prev;
        }
    }
    return registerClass(className, readableName);
}

void XMLTreeModel::setProperties(TreeItem *item, const Properties &props) {
    if (item == NULL)
        return;
    item->setProperties(props);
    int childCount = item->childCount();
    if (childCount > 0) {
        emit dataChanged(createIndex(0, 0, item->child(0)),
            createIndex(item->row(), 1, item->child(childCount-1)));
    }
}

void XMLTreeModel::populate(const QString &className, TreeItem *parent) {
    QDomElement plugin;
    for (QDomElement e = m_docRoot.firstChildElement("plugin"); !e.isNull();
         e = e.nextSiblingElement("plugin")) {
        if (e.attribute("className") == className) {
            plugin = e;
            break;
        }
    }

    if (plugin.hasAttribute("extends"))
        populate(plugin.attribute("extends"), parent);

    if (plugin.isNull())
        return;

    for (QDomElement e = plugin.firstChildElement("param"); !e.isNull();
        e = e.nextSiblingElement("param")) {
        QDomDocument tooltipDoc;
        QDomElement root = tooltipDoc.createElement("p");
        tooltipDoc.appendChild(root);

        for (QDomNode child = e.firstChild(); !child.isNull();
            child = child.nextSibling()) {
            root.appendChild(tooltipDoc.importNode(child, true));
        }
        QString toolTip = tooltipDoc.toString(),
            type = e.attribute("type"),
            value = e.attribute("default"),
            name = e.attribute("name"),
            readableName = e.attribute("readableName"),
            importanceValue = e.attribute("importance");

        QVariant variantValue, defaultValue;
        bool ok = true;
        if (type == "string")
            variantValue = QVariant(value);
        else if (type == "integer" || type == "long") {
            variantValue = QVariant((int) value.toInt(&ok));
        } else if (type == "float")
            variantValue = QVariant(value.toDouble(&ok));
        else if (type == "boolean") {
            value = value.toLower();
            if (value == "true")
                variantValue = QVariant(true);
            else if (value == "false")
                variantValue = QVariant(false);
            else
                ok = false;
        } else {
            SLog(EError, "Unexpected property type!");
        }

        if (!ok)
            SLog(EError, "Error in number conversion!");

        int importance = 10;
        if (importanceValue != "")
            importance = importanceValue.toInt();

        TreeItem *item = new TreeItem(name, readableName,
            variantValue, variantValue, parent, importance);

        item->setToolTip(toolTip);
        parent->appendChild(item);
    }
}

XMLTreeModel::~XMLTreeModel() {
    delete m_rootItem;
}

int XMLTreeModel::columnCount(const QModelIndex &) const {
    return 2;
}

QVariant XMLTreeModel::data(const QModelIndex &index, int role) const {
    if (!index.isValid())
        return QVariant();

    TreeItem *item = static_cast<TreeItem*>(index.internalPointer());

    if (role == Qt::DisplayRole || role == Qt::EditRole)
        return item->data(index.column());
    else if (role == Qt::BackgroundColorRole && item->isCategory())
        return QVariant(QColor(qRgb(0x90, 0x90, 0x90)));
    else if (role == Qt::TextColorRole && item->isCategory())
        return QVariant(QColor(qRgb(0xff, 0xff, 0xff)));
    else if (role == Qt::ToolTipRole)
        return QVariant(item->toolTip());
    else if (role == Qt::ForegroundRole && index.column() == 0)
        return item->getImportance() < 5 ?
            m_unimportantColor : QVariant();
    else
        return QVariant();
}

bool XMLTreeModel::setData(const QModelIndex &index, const QVariant &value,
                        int role) {
    if (role != Qt::EditRole)
        return false;

    TreeItem *item = static_cast<TreeItem*>(index.internalPointer());
    bool result = item->setData(index.column(), value);

    if (result)
        emit dataChanged(index, index);

    return result;
}

Qt::ItemFlags XMLTreeModel::flags(const QModelIndex &index) const {
    if (!index.isValid())
        return 0;

    TreeItem *item = static_cast<TreeItem*>(index.internalPointer());
    if (item->isCategory())
        return Qt::ItemIsEnabled;
    else if (index.column() == 0)
        return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
    else
        return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable;

}


QVariant XMLTreeModel::headerData(int section, Qt::Orientation orientation,
                               int role) const {
    if (orientation == Qt::Horizontal && role == Qt::DisplayRole) {
        if (section == 0)
            return QVariant(tr("Property"));
        else
            return QVariant(tr("Value"));
    }

    return QVariant();
}

QModelIndex XMLTreeModel::index(int row, int column, const QModelIndex &parent)
            const {
    if (!hasIndex(row, column, parent))
        return QModelIndex();

    TreeItem *parentItem;

    if (!parent.isValid())
        parentItem = m_rootItem;
    else
        parentItem = static_cast<TreeItem*>(parent.internalPointer());

    TreeItem *childItem = parentItem->child(row);
    if (childItem)
        return createIndex(row, column, childItem);
    else
        return QModelIndex();
}

QModelIndex XMLTreeModel::parent(const QModelIndex &index) const {
    if (!index.isValid())
        return QModelIndex();

    TreeItem *childItem = static_cast<TreeItem*>(index.internalPointer());
    TreeItem *parentItem = childItem->parent();

    if (parentItem == m_rootItem)
        return QModelIndex();

    return createIndex(parentItem->row(), 0, parentItem);
}

int XMLTreeModel::rowCount(const QModelIndex &parent) const {
    TreeItem *parentItem;
    if (parent.column() > 0)
        return 0;

    if (!parent.isValid())
        parentItem = m_rootItem;
    else
        parentItem = static_cast<TreeItem*>(parent.internalPointer());

    return parentItem->childCount();
}
