#if !defined(__ABOUTDLG_H)
#define __ABOUTDLG_H

#include "common.h"

namespace Ui {
	class AboutDialog;
}

class AboutDialog : public QDialog {
    Q_OBJECT
public:
	AboutDialog(QWidget *parent);
	~AboutDialog();

protected slots:
	void onCredits();
protected:
    void changeEvent(QEvent *e);
private:
	Ui::AboutDialog *ui;
};

#endif // __ABOUTDLG_H
