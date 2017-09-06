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

#if !defined(__RENDERSETTINGSDLG_H)
#define __RENDERSETTINGSDLG_H

#include "xmltreemodel.h"

namespace Ui {
    class RenderSettingsDialog;
}

/* Custom delegate for rendering & editing property data */
class PropertyDelegate : public QStyledItemDelegate {
    Q_OBJECT
public:
    PropertyDelegate(QObject *parent = NULL);
    virtual ~PropertyDelegate();

    QString displayText(const QVariant &value, const QLocale &locale) const;
    QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
            const QModelIndex &index) const;
    void setEditorData(QWidget *editor, const QModelIndex &index) const;
    void setModelData(QWidget *editor, QAbstractItemModel *model,
        const QModelIndex &index) const;
    void updateEditorGeometry(QWidget *editor,
        const QStyleOptionViewItem &option, const QModelIndex &index) const;
protected slots:
    void updateWidgetData();
};

class RenderSettingsDialog : public QDialog {
    Q_OBJECT
public:
    RenderSettingsDialog(QWidget *parent = 0);
    virtual ~RenderSettingsDialog();

    void load(const SceneContext *scene);
    void apply(SceneContext *scene);
    bool resolutionHasChanged() const;
protected slots:
    void onTreeSelectionChange(const QItemSelection &selected, const QItemSelection &deselected);
    void cbHighlighted(int index);
    void chkBoxPressed();
    void update();
    void refresh();
    void dataChanged();
protected:
    void changeEvent(QEvent *e);
    void setDocumentation(const QString &text);
    QStringList validateConfiguration() const;
private:
    Ui::RenderSettingsDialog *ui;
    QDomDocument m_document;
    XMLTreeModel *m_model;
    TreeItem *m_integratorNode, *m_samplerNode;
    TreeItem *m_rFilterNode, *m_icNode;
    TreeItem *m_aiNode;
    QString m_originalResolution;
    QString m_currentDocumentation;
    QStringList m_statusMessages;
};

#endif // __RENDERSETTINGSDLG_H
