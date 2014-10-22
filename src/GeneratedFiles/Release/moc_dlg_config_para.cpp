/****************************************************************************
** Meta object code from reading C++ file 'dlg_config_para.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.1.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../UI/dlg_config_para.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'dlg_config_para.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.1.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_ConfigParaDlg_t {
    QByteArrayData data[10];
    char stringdata[103];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    offsetof(qt_meta_stringdata_ConfigParaDlg_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData) \
    )
static const qt_meta_stringdata_ConfigParaDlg_t qt_meta_stringdata_ConfigParaDlg = {
    {
QT_MOC_LITERAL(0, 0, 13),
QT_MOC_LITERAL(1, 14, 16),
QT_MOC_LITERAL(2, 31, 0),
QT_MOC_LITERAL(3, 32, 11),
QT_MOC_LITERAL(4, 44, 7),
QT_MOC_LITERAL(5, 52, 4),
QT_MOC_LITERAL(6, 57, 8),
QT_MOC_LITERAL(7, 66, 9),
QT_MOC_LITERAL(8, 76, 12),
QT_MOC_LITERAL(9, 89, 12)
    },
    "ConfigParaDlg\0parameterChanged\0\0"
    "initWidgets\0getData\0_val\0getRigid\0"
    "getSmooth\0getTolerance\0applyChanges\0"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_ConfigParaDlg[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   49,    2, 0x05,

 // slots: name, argc, parameters, tag, flags
       3,    0,   50,    2, 0x0a,
       4,    1,   51,    2, 0x08,
       6,    1,   54,    2, 0x08,
       7,    1,   57,    2, 0x08,
       8,    1,   60,    2, 0x08,
       9,    0,   63,    2, 0x08,

 // signals: parameters
    QMetaType::Void,

 // slots: parameters
    QMetaType::Bool,
    QMetaType::Void, QMetaType::Double,    5,
    QMetaType::Void, QMetaType::Double,    5,
    QMetaType::Void, QMetaType::Double,    5,
    QMetaType::Void, QMetaType::Double,    5,
    QMetaType::Void,

       0        // eod
};

void ConfigParaDlg::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        ConfigParaDlg *_t = static_cast<ConfigParaDlg *>(_o);
        switch (_id) {
        case 0: _t->parameterChanged(); break;
        case 1: { bool _r = _t->initWidgets();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 2: _t->getData((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 3: _t->getRigid((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 4: _t->getSmooth((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 5: _t->getTolerance((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 6: _t->applyChanges(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (ConfigParaDlg::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&ConfigParaDlg::parameterChanged)) {
                *result = 0;
            }
        }
    }
}

const QMetaObject ConfigParaDlg::staticMetaObject = {
    { &QFrame::staticMetaObject, qt_meta_stringdata_ConfigParaDlg.data,
      qt_meta_data_ConfigParaDlg,  qt_static_metacall, 0, 0}
};


const QMetaObject *ConfigParaDlg::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *ConfigParaDlg::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ConfigParaDlg.stringdata))
        return static_cast<void*>(const_cast< ConfigParaDlg*>(this));
    return QFrame::qt_metacast(_clname);
}

int ConfigParaDlg::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QFrame::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 7;
    }
    return _id;
}

// SIGNAL 0
void ConfigParaDlg::parameterChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}
QT_END_MOC_NAMESPACE
