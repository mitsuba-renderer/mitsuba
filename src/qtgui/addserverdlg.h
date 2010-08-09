#if !defined(__ADDSERVERDLG_H)
#define __ADDSERVERDLG_H

#include "common.h"

namespace Ui {
	class AddServerDialog;
}

class AddServerDialog : public QDialog {
    Q_OBJECT
public:
	AddServerDialog(QWidget *parent);
	~AddServerDialog();

	ServerConnection getConnection() const;
protected slots:
	void on_sshConnection_toggled();
	void on_directConnection_toggled();
	void validate();
protected:
    void changeEvent(QEvent *e);
private:
	Ui::AddServerDialog *ui;
};

#endif // __ADDSERVERDLG_H
