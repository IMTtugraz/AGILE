/****************************************************************************
** Meta object code from reading C++ file 'imt_mrgui.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../apps/imt_irgn/gui/imt_mrgui.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'imt_mrgui.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_IMT_MRGui[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      16,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: signature, parameters, type, tag, flags
      29,   11,   10,   10, 0x05,
      64,   52,   10,   10, 0x05,
     101,   85,   10,   10, 0x05,

 // slots: signature, parameters, type, tag, flags
     131,   10,   10,   10, 0x08,
     166,   10,   10,   10, 0x08,
     201,   10,   10,   10, 0x08,
     225,   10,   10,   10, 0x08,
     256,   11,   10,   10, 0x08,
     279,   52,   10,   10, 0x08,
     300,   85,   10,   10, 0x08,
     330,   10,   10,   10, 0x08,
     375,  357,   10,   10, 0x08,
     410,   10,   10,   10, 0x08,
     436,   10,   10,   10, 0x08,
     473,   10,   10,   10, 0x08,
     507,   10,   10,   10, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_IMT_MRGui[] = {
    "IMT_MRGui\0\0infotext,progress\0"
    "sendInfo(QString,bool)\0warningtext\0"
    "sendWarning(QString)\0title,extension\0"
    "sendSaveFile(QString,QString)\0"
    "on_actionRead_DicomRaw_triggered()\0"
    "on_actionOpen_meas_dat_triggered()\0"
    "on_pushButton_clicked()\0"
    "on_actionParameter_triggered()\0"
    "set_Info(QString,bool)\0set_Warning(QString)\0"
    "set_SaveFile(QString,QString)\0"
    "on_pb_outputpath_clicked()\0filepath,filename\0"
    "startauto(QStringList,QStringList)\0"
    "on_pb_automatic_clicked()\0"
    "on_actionOpen_matlab_bin_triggered()\0"
    "on_actionAuto_Archive_triggered()\0"
    "on_actionSpecial_triggered()\0"
};

void IMT_MRGui::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        IMT_MRGui *_t = static_cast<IMT_MRGui *>(_o);
        switch (_id) {
        case 0: _t->sendInfo((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 1: _t->sendWarning((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 2: _t->sendSaveFile((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2]))); break;
        case 3: _t->on_actionRead_DicomRaw_triggered(); break;
        case 4: _t->on_actionOpen_meas_dat_triggered(); break;
        case 5: _t->on_pushButton_clicked(); break;
        case 6: _t->on_actionParameter_triggered(); break;
        case 7: _t->set_Info((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< bool(*)>(_a[2]))); break;
        case 8: _t->set_Warning((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 9: _t->set_SaveFile((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2]))); break;
        case 10: _t->on_pb_outputpath_clicked(); break;
        case 11: _t->startauto((*reinterpret_cast< QStringList(*)>(_a[1])),(*reinterpret_cast< QStringList(*)>(_a[2]))); break;
        case 12: _t->on_pb_automatic_clicked(); break;
        case 13: _t->on_actionOpen_matlab_bin_triggered(); break;
        case 14: _t->on_actionAuto_Archive_triggered(); break;
        case 15: _t->on_actionSpecial_triggered(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData IMT_MRGui::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject IMT_MRGui::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_IMT_MRGui,
      qt_meta_data_IMT_MRGui, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &IMT_MRGui::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *IMT_MRGui::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *IMT_MRGui::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_IMT_MRGui))
        return static_cast<void*>(const_cast< IMT_MRGui*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int IMT_MRGui::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 16)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 16;
    }
    return _id;
}

// SIGNAL 0
void IMT_MRGui::sendInfo(QString _t1, bool _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void IMT_MRGui::sendWarning(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void IMT_MRGui::sendSaveFile(QString _t1, QString _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}
QT_END_MOC_NAMESPACE
