#include "mainwindow.h"
#include "tabbar.h"

#if defined(__OSX__)
/** Fix a bug in Qt/Cocoa where a tab is not drawn focus after a sub-dialog receives focus **/
class TabBarProxyStyle : public QProxyStyle {
public:
	TabBarProxyStyle(MainWindow *parent) : QProxyStyle(NULL), m_parent(parent) {
	}

	virtual void drawControl(ControlElement ce, const QStyleOption * option, QPainter *painter, const QWidget *widget = 0) const {
		if (ce == CE_TabBarTabShape) {
			if (const QStyleOptionTabV3 *tabOptV3 = qstyleoption_cast<const QStyleOptionTabV3 *>(option)) {
				QStyleOptionTabV3 so(*tabOptV3);
				if (m_parent->isActive())
					so.state |= QStyle::State_Active;
				
				QProxyStyle::drawControl(ce, &so, painter, widget);
				return;
			}
		}
		QProxyStyle::drawControl(ce, option, painter, widget);
	}
private:
	MainWindow *m_parent;
};
#endif

CustomTabBar::CustomTabBar(QWidget *parent) : QTabBar(parent) {
#if defined(__OSX__)
	qApp->setStyle(new TabBarProxyStyle((MainWindow *) parent->window()));
#endif
}

CustomTabBar::~CustomTabBar() {
}
