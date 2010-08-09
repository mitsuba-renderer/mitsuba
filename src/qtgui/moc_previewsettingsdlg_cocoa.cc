/****************************************************************************
** Meta object code from reading C++ file 'previewsettingsdlg_cocoa.h'
**
** Created: Sun Jul 18 22:08:20 2010
**      by: The Qt Meta Object Compiler version 61 (Qt 4.5.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "previewsettingsdlg_cocoa.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'previewsettingsdlg_cocoa.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 61
#error "This file was generated using the moc from 4.5.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_PreviewSettingsDlg[] = {

 // content:
       2,       // revision
       0,       // classname
       0,    0, // classinfo
      10,   12, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors

 // signals: signature, parameters, type, tag, flags
      27,   20,   19,   19, 0x05,
      61,   50,   19,   19, 0x05,
     102,   93,   19,   19, 0x05,
     134,  125,   19,   19, 0x05,
     168,  157,   19,   19, 0x05,
     200,  193,   19,   19, 0x05,
     237,  193,   19,   19, 0x05,
     286,  282,   19,   19, 0x05,
     317,  312,   19,   19, 0x05,
     344,   19,   19,   19, 0x05,

       0        // eod
};

static const char qt_meta_stringdata_PreviewSettingsDlg[] = {
    "PreviewSettingsDlg\0\0length\0"
    "pathLengthChanged(int)\0resolution\0"
    "shadowMapResolutionChanged(int)\0"
    "clamping\0clampingChanged(Float)\0"
    "exposure\0exposureChanged(Float)\0"
    "srgb,gamma\0gammaChanged(bool,Float)\0"
    "method\0previewMethodChanged(EPreviewMethod)\0"
    "toneMappingMethodChanged(EToneMappingMethod)\0"
    "key\0reinhardKeyChanged(Float)\0burn\0"
    "reinhardBurnChanged(Float)\0close()\0"
};

const QMetaObject PreviewSettingsDlg::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_PreviewSettingsDlg,
      qt_meta_data_PreviewSettingsDlg, 0 }
};

const QMetaObject *PreviewSettingsDlg::metaObject() const
{
    return &staticMetaObject;
}

void *PreviewSettingsDlg::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_PreviewSettingsDlg))
        return static_cast<void*>(const_cast< PreviewSettingsDlg*>(this));
    return QObject::qt_metacast(_clname);
}

int PreviewSettingsDlg::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: pathLengthChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: shadowMapResolutionChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: clampingChanged((*reinterpret_cast< Float(*)>(_a[1]))); break;
        case 3: exposureChanged((*reinterpret_cast< Float(*)>(_a[1]))); break;
        case 4: gammaChanged((*reinterpret_cast< bool(*)>(_a[1])),(*reinterpret_cast< Float(*)>(_a[2]))); break;
        case 5: previewMethodChanged((*reinterpret_cast< EPreviewMethod(*)>(_a[1]))); break;
        case 6: toneMappingMethodChanged((*reinterpret_cast< EToneMappingMethod(*)>(_a[1]))); break;
        case 7: reinhardKeyChanged((*reinterpret_cast< Float(*)>(_a[1]))); break;
        case 8: reinhardBurnChanged((*reinterpret_cast< Float(*)>(_a[1]))); break;
        case 9: close(); break;
        default: ;
        }
        _id -= 10;
    }
    return _id;
}

// SIGNAL 0
void PreviewSettingsDlg::pathLengthChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void PreviewSettingsDlg::shadowMapResolutionChanged(int _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void PreviewSettingsDlg::clampingChanged(Float _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void PreviewSettingsDlg::exposureChanged(Float _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void PreviewSettingsDlg::gammaChanged(bool _t1, Float _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void PreviewSettingsDlg::previewMethodChanged(EPreviewMethod _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}

// SIGNAL 6
void PreviewSettingsDlg::toneMappingMethodChanged(EToneMappingMethod _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}

// SIGNAL 7
void PreviewSettingsDlg::reinhardKeyChanged(Float _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 7, _a);
}

// SIGNAL 8
void PreviewSettingsDlg::reinhardBurnChanged(Float _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 8, _a);
}

// SIGNAL 9
void PreviewSettingsDlg::close()
{
    QMetaObject::activate(this, &staticMetaObject, 9, 0);
}
QT_END_MOC_NAMESPACE
