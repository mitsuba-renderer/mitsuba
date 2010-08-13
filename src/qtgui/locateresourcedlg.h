#if !defined(__LOCATERESOURCEDLG_H)
#define __LOCATERESOURCEDLG_H

#include "common.h"

namespace Ui {
	class LocateResourceDialog;
}

class LocateResourceDialog : public QDialog {
    Q_OBJECT
public:
	LocateResourceDialog(QWidget *parent, const QString &resourceName);
	~LocateResourceDialog();

	inline const QString &getFilename() const {
		return m_filename;
	}
protected slots:
	void on_pathBrowse_clicked();
protected:
    void changeEvent(QEvent *e);
private:
	Ui::LocateResourceDialog *ui;
	QString m_filename;
};

#endif // __LOCATERESOURCEDLG_H
