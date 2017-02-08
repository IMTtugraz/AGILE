/****************************************************************************
** Meta object code from reading C++ file 'pathsetting.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../../apps/imt_irgn/gui/widgets/pathsetting.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'pathsetting.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_PathSetting[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
       6,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x08,
      32,   12,   12,   12, 0x08,
      51,   12,   12,   12, 0x08,
      71,   12,   12,   12, 0x08,
      82,   12,   12,   12, 0x08,
      97,   12,   12,   12, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_PathSetting[] = {
    "PathSetting\0\0browsedcmrawpath()\0"
    "browseoutputpath()\0browsearchivepath()\0"
    "okbutton()\0cancelbutton()\0archiveactive()\0"
};

void PathSetting::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        PathSetting *_t = static_cast<PathSetting *>(_o);
        switch (_id) {
        case 0: _t->browsedcmrawpath(); break;
        case 1: _t->browseoutputpath(); break;
        case 2: _t->browsearchivepath(); break;
        case 3: _t->okbutton(); break;
        case 4: _t->cancelbutton(); break;
        case 5: _t->archiveactive(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObjectExtraData PathSetting::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject PathSetting::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_PathSetting,
      qt_meta_data_PathSetting, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &PathSetting::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *PathSetting::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *PathSetting::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_PathSetting))
        return static_cast<void*>(const_cast< PathSetting*>(this));
    return QDialog::qt_metacast(_clname);
}

int PathSetting::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 6)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 6;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
