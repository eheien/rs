/****************************************************************************
** Meta object code from reading C++ file 'RateStateSimWindow.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "RateStateSimWindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'RateStateSimWindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_RateStateSimWindow[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      19,   43,   51,   51, 0x0a,
      52,   43,   51,   51, 0x0a,
      77,   43,   51,   51, 0x0a,
     101,   43,   51,   51, 0x0a,
     125,   51,   51,   51, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_RateStateSimWindow[] = {
    "RateStateSimWindow\0a_param_changed(double)\0"
    "new_val\0\0ab_ratio_changed(double)\0"
    "k_param_changed(double)\0r_param_changed(double)\0"
    "recalc()\0"
};

void RateStateSimWindow::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        RateStateSimWindow *_t = static_cast<RateStateSimWindow *>(_o);
        switch (_id) {
        case 0: _t->a_param_changed((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 1: _t->ab_ratio_changed((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 2: _t->k_param_changed((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 3: _t->r_param_changed((*reinterpret_cast< const double(*)>(_a[1]))); break;
        case 4: _t->recalc(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData RateStateSimWindow::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject RateStateSimWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_RateStateSimWindow,
      qt_meta_data_RateStateSimWindow, &staticMetaObjectExtraData }
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
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
