/****************************************************************************
** Meta object code from reading C++ file 'RateStateSimWindow.h'
**
** Created: Mon Apr 23 13:51:39 2012
**      by: The Qt Meta Object Compiler version 62 (Qt 4.7.4)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "RateStateSimWindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'RateStateSimWindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.7.4. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_RateStateSimWindow[] = {

 // content:
       5,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      28,   20,   19,   19, 0x0a,
      52,   20,   19,   19, 0x0a,
      77,   20,   19,   19, 0x0a,
     101,   20,   19,   19, 0x0a,
     125,   19,   19,   19, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_RateStateSimWindow[] = {
    "RateStateSimWindow\0\0new_val\0"
    "a_param_changed(double)\0"
    "ab_ratio_changed(double)\0"
    "k_param_changed(double)\0r_param_changed(double)\0"
    "recalc()\0"
};

const QMetaObject RateStateSimWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_RateStateSimWindow,
      qt_meta_data_RateStateSimWindow, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &RateStateSimWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *RateStateSimWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *RateStateSimWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_RateStateSimWindow))
        return static_cast<void*>(const_cast< RateStateSimWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int RateStateSimWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: a_param_changed((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 1: ab_ratio_changed((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 2: k_param_changed((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 3: r_param_changed((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 4: recalc(); break;
        default: ;
        }
        _id -= 5;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
