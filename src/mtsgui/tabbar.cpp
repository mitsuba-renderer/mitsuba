/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

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
