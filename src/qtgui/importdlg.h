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
protected:
    void changeEvent(QEvent *e);
private:
	Ui::ImportDialog *ui;
};

#endif // __IMPORTDLG_H
