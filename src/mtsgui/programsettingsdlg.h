#if !defined(__PROGRAMSETTINGSDLG_H)
#define __PROGRAMSETTINGSDLG_H

#include "common.h"
#include "ui_programsettingsdlg.h"
#include <mitsuba/mitsuba.h>

class ProgramSettingsDialog : public QDialog {
    Q_OBJECT
public:
    ProgramSettingsDialog(QWidget *parent);
    ~ProgramSettingsDialog();

    inline QList<ServerConnection> &getConnections() {
        return m_connections;
    }

    void setConnections(QList<ServerConnection> &connections);

    inline QString getNodeName() const {
        return ui->nodeName->text();
    }

    inline void setNodeName(const QString &name) const {
        ui->nodeName->setText(name);
    }

    inline int getListenPort() const {
        return ui->listenPort->text().toInt();
    }

    inline void setListenPort(int port) const {
        ui->listenPort->setText(QString("%1").arg(port));
    }

    inline bool getInvertMouse() const {
        return ui->invertMouseBox->checkState() == Qt::Checked;
    }

    inline ENavigationMode getNavigationMode() const {
        return (ENavigationMode) ui->navigationModeBox->currentIndex();
    }

    inline void setNavigationMode(ENavigationMode mode) const {
        ui->navigationModeBox->setCurrentIndex(mode);
    }

    inline void setInvertMouse(bool value) {
        ui->invertMouseBox->setCheckState(value ? Qt::Checked : Qt::Unchecked);
    }

    inline int getMouseSensitivity() const {
        return ui->sensitivitySlider->value();
    }

    inline void setMouseSensitivity(int value) {
        ui->sensitivitySlider->setValue(value);
    }

    inline bool getCheckForUpdates() const {
        return ui->checkForUpdatesBox->checkState() == Qt::Checked;
    }

    inline void setCheckForUpdates(bool value) {
        ui->checkForUpdatesBox->setCheckState(value ? Qt::Checked : Qt::Unchecked);
    }

    inline void setSearchPaths(const QStringList &paths) {
        ui->searchPathList->clear();
        ui->searchPathList->addItems(paths);
    }

    inline void setLocalWorkerCount(size_t count) {
        ui->localWorkerBox->setValue((int) count);
    }

    inline int getLocalWorkerCount() {
        return ui->localWorkerBox->value();
    }

    inline QStringList getSearchPaths() const {
        QStringList result;
        for (int i=0; i<ui->searchPathList->count(); ++i)
            result << ui->searchPathList->item(i)->text();
        return result;
    }

    inline void setLogLevel(int logLevel) {
        switch (logLevel) {
            case ETrace: ui->logVerbosityBox->setCurrentIndex(0); break;
            case EDebug: ui->logVerbosityBox->setCurrentIndex(1); break;
            case EInfo: ui->logVerbosityBox->setCurrentIndex(2); break;
            case EWarn: ui->logVerbosityBox->setCurrentIndex(3); break;
            case EError: ui->logVerbosityBox->setCurrentIndex(4); break;
            default:
                SLog(EError, "Unknown verbosity level!");
        }
    }

    inline ELogLevel getLogLevel() const {
        switch (ui->logVerbosityBox->currentIndex()) {
            case 0: return ETrace;
            case 1: return EDebug;
            case 2: return EInfo;
            case 3: return EWarn;
            case 4: return EError;
            default:
                SLog(EError, "Unknown verbosity level!");
                return EDebug;
        }
    }

    inline void setWorkerPriority(Thread::EThreadPriority logLevel) {
        switch (logLevel) {
            case Thread::EIdlePriority: ui->workerPriorityBox->setCurrentIndex(0); break;
            case Thread::ELowestPriority: ui->workerPriorityBox->setCurrentIndex(1); break;
            case Thread::ELowPriority: ui->workerPriorityBox->setCurrentIndex(2); break;
            case Thread::ENormalPriority: ui->workerPriorityBox->setCurrentIndex(3); break;
            case Thread::EHighPriority: ui->workerPriorityBox->setCurrentIndex(4); break;
            case Thread::EHighestPriority: ui->workerPriorityBox->setCurrentIndex(5); break;
            default:
                SLog(EError, "Unknown worker priority!");
        }
    }

    inline Thread::EThreadPriority getWorkerPriority() const {
        switch (ui->workerPriorityBox->currentIndex()) {
            case 0: return Thread::EIdlePriority;
            case 1: return Thread::ELowestPriority;
            case 2: return Thread::ELowPriority;
            case 3: return Thread::ENormalPriority;
            case 4: return Thread::EHighPriority;
            case 5: return Thread::EHighestPriority;
            default:
                SLog(EError, "Unknown worker priority!");
                return Thread::ENormalPriority;
        }
    }

    void setBlockSize(int size) {
        switch (size) {
            case 2: ui->blockSizeBox->setCurrentIndex(0); break;
            case 4: ui->blockSizeBox->setCurrentIndex(1); break;
            case 8: ui->blockSizeBox->setCurrentIndex(2); break;
            case 16: ui->blockSizeBox->setCurrentIndex(3); break;
            case 32: ui->blockSizeBox->setCurrentIndex(4); break;
            case 64: ui->blockSizeBox->setCurrentIndex(5); break;
            case 128: ui->blockSizeBox->setCurrentIndex(6); break;
            default:
                SLog(EError, "Unknown block size!");
        }
    }

    inline int getBlockSize() const {
        switch (ui->blockSizeBox->currentIndex()) {
            case 0: return 2;
            case 1: return 4;
            case 2: return 8;
            case 3: return 16;
            case 4: return 32;
            case 5: return 64;
            case 6: return 128;
            default:
                SLog(EError, "Unknown block size!");
                return 0;
        }
    }
protected:
    void changeEvent(QEvent *e);

public slots:
    void refresh();

protected slots:
    void on_addPathButton_clicked();
    void on_removePathButton_clicked();
    void on_addConnectionButton_clicked();
    void on_removeConnectionButton_clicked();
    void on_searchPathList_currentItemChanged(QListWidgetItem *cur, QListWidgetItem *prev);
    void on_connectionList_currentItemChanged(QListWidgetItem *cur, QListWidgetItem *prev);

private:
    Ui::ProgramSettingsDialog *ui;
    QList<ServerConnection> m_connections;
};

#endif // __PROGRAMSETTINGSDLG_H
