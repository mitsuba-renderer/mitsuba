#if !defined(__ACKNOWLEDGEMENTDLG_H)
#define __ACKNOWLEDGEMENTDLG_H

#include "common.h"

namespace Ui {
	class AcknowledgmentDialog;
}

class AcknowledgmentDialog : public QDialog {
    Q_OBJECT
public:
	AcknowledgmentDialog(QWidget *parent);
	~AcknowledgmentDialog();
protected:
    void changeEvent(QEvent *e);
private:
	Ui::AcknowledgmentDialog *ui;
};

#endif // __ACKNOWLEDGEMENTDLG_H
