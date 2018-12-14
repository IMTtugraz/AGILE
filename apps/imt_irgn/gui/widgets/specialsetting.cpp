#include <Qt>
#include "specialsetting.h"

//! [0]
SpecialSetting::SpecialSetting(QWidget *parent)
    : QDialog(parent)
{
  readSettings();

  QString str_autoloadtimeout;
  str_autoloadtimeout.setNum(_autoload_timeout);

  _le_autoload_timeout = createLineEdit(str_autoloadtimeout);

  _lb_autoload_timeout = new QLabel("Auto Load Timeout: ");

  _spacer_line = new QSpacerItem( 100, 20, QSizePolicy::Minimum,
                                 QSizePolicy::Expanding );


  QPushButton *pb_ok = createButton("Ok", SLOT(okbutton()));
  QPushButton *pb_cancle = createButton("Cancel", SLOT(cancelbutton()));

  QHBoxLayout* buttonsLayout = new QHBoxLayout;
  //buttonsLayout->addStretch();
  buttonsLayout->addItem(_spacer_line);
  buttonsLayout->addWidget(pb_ok);
  buttonsLayout->addWidget(pb_cancle);


    QGridLayout *pathLayout = new QGridLayout;
    pathLayout->addWidget(_lb_autoload_timeout, 0, 0);
    pathLayout->addWidget(_le_autoload_timeout, 0, 1);

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addLayout(pathLayout);
    layout->addLayout(buttonsLayout);

    setLayout(layout);

    setWindowTitle(tr("Special Settings"));
    //resize(150,75);
}

QPushButton *SpecialSetting::createButton(const QString &text, const char *member)
{
    QPushButton *button = new QPushButton(text);
    connect(button, SIGNAL(clicked()), this, member);
    return button;
}

QLineEdit *SpecialSetting::createLineEdit(const QString &text)
{
    QLineEdit *LineEdit = new QLineEdit;
    LineEdit->setText(text);
    LineEdit->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    return LineEdit;
}

void SpecialSetting::readSettings()
{
    QSettings settings("Cuda IRGN -GUI", "IRGN Reconstruction - Read Dicom Raw / meas.dat - Write Dicom by Heigl");
    _autoload_timeout = settings.value("autoload_timeout", 20).toUInt();

    _autoload_timeout_start = _autoload_timeout;
}

void SpecialSetting::writeSettings()
{
    QSettings settings("Cuda IRGN -GUI", "IRGN Reconstruction - Read Dicom Raw / meas.dat - Write Dicom by Heigl");
    settings.setValue("autoload_timeout", _autoload_timeout);
}

bool SpecialSetting::pathChanged()
{
  //qDebug()<<"\nasdfköj: "<<_archiveactive_start<<"  "<<_archiveactive;
  if ((_autoload_timeout_start == _autoload_timeout))
     return false;    //nothing changed
  else
    return true;

}

void SpecialSetting::cancelbutton()
{
  QString str_autoloadtimeout;
  str_autoloadtimeout.setNum(_autoload_timeout);

  _le_autoload_timeout->setText(str_autoloadtimeout);

  this->close();
}


void SpecialSetting::okbutton()
{
  _autoload_timeout = _le_autoload_timeout->text().toUInt();

  this->close();
}
