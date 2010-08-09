#if !defined(__TABBAR_H)
#define __TABBAR_H

#include <QtGui>

class CustomTabBar : public QTabBar {
    Q_OBJECT
public:
	CustomTabBar(QWidget *parent);
	virtual ~CustomTabBar();
};

#endif // __TABBAR_H
