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

	QString	displayText(const QVariant &value, const QLocale &locale) const;
	QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &option,
			const QModelIndex &index) const;
	void setEditorData(QWidget *editor, const QModelIndex &index) const;
	void setModelData(QWidget *editor, QAbstractItemModel *model,
		const QModelIndex &index) const;
	void updateEditorGeometry(QWidget *editor,
		const QStyleOptionViewItem &option, const QModelIndex &index) const;
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
