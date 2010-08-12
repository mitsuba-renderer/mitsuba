#if !defined(__IMPORTDLG_H)
#define __IMPORTDLG_H

#include "common.h"

namespace Ui {
	class ImportDialog;
}

class ImportDialog : public QDialog {
    Q_OBJECT
public:
	ImportDialog(QWidget *parent);
	~ImportDialog();
public slots:
	void accept();
protected slots:
	void on_inputBrowse_clicked(bool checked);
	void on_directoryBrowse_clicked(bool checked);
	void on_adjustmentBrowse_clicked(bool checked);
	void refresh();
protected:
    void changeEvent(QEvent *e);
private:
	Ui::ImportDialog *ui;
};

#endif // __IMPORTDLG_H
