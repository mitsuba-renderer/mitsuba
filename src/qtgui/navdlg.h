#if !defined(__NAVIGATIONDLG_H)
#define __NAVIGATIONDLG_H

#include "common.h"

namespace Ui {
	class NavigationDialog;
}

class NavigationDialog : public QDialog {
    Q_OBJECT
public:
	NavigationDialog(QWidget *parent);
	virtual ~NavigationDialog();
protected:
    void changeEvent(QEvent *e);
private:
	Ui::NavigationDialog *ui;
};

#endif // __NAVIGATIONDLG_H
