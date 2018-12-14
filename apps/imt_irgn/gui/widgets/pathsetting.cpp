#include <Qt>
#include "pathsetting.h"

//! [0]
PathSetting::PathSetting(QWidget *parent)
    : QDialog(parent)
{
  readSettings();

  _dcmrawpathLineEdit = createLineEdit(_dcmrawpath);
  _outputpathLineEdit = createLineEdit(_outputpath);
  _archivepathLineEdit = createLineEdit(_archivepath);

  _dcmrawpathbrowseButton = createButton(tr("&Browse..."), SLOT(browsedcmrawpath()));
  _outputpathbrowseButton = createButton(tr("&Browse..."), SLOT(browseoutputpath()));
  _archivepathbrowseButton = createButton(tr("&Browse..."), SLOT(browsearchivepath()));

  _dcmrawpathLabel = new QLabel("Dicom-RAW Path:");
  _outputpathLabel = new QLabel("Output Path:   ");

  _cb_archiveactive = new QCheckBox("Auto Archive:");
  _cb_archiveactive->setChecked(_archiveactive);
  archiveactive();
  connect(_cb_archiveactive, SIGNAL(stateChanged(int)), this, SLOT(archiveactive()));

  _spacer_line = new QSpacerItem( 450, 20, QSizePolicy::Minimum,
                                 QSizePolicy::Expanding );


  QPushButton *pb_ok = createButton("Ok", SLOT(okbutton()));
  QPushButton *pb_cancle = createButton("Cancel", SLOT(cancelbutton()));

  QHBoxLayout* buttonsLayout = new QHBoxLayout;
  //buttonsLayout->addStretch();
  buttonsLayout->addItem(_spacer_line);
  buttonsLayout->addWidget(pb_ok);
  buttonsLayout->addWidget(pb_cancle);


    QGridLayout *pathLayout = new QGridLayout;
    pathLayout->addWidget(_dcmrawpathLabel, 0, 0);
    pathLayout->addWidget(_dcmrawpathLineEdit, 0, 1);
    pathLayout->addWidget(_dcmrawpathbrowseButton, 0, 2);

    pathLayout->addWidget(_outputpathLabel, 1, 0);
    pathLayout->addWidget(_outputpathLineEdit, 1, 1);
    pathLayout->addWidget(_outputpathbrowseButton, 1, 2);

    pathLayout->addWidget(_cb_archiveactive, 3, 0);
    pathLayout->addWidget(_archivepathLineEdit, 3, 1);
    pathLayout->addWidget(_archivepathbrowseButton, 3, 2);

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addLayout(pathLayout);
    layout->addLayout(buttonsLayout);

    setLayout(layout);

    setWindowTitle(tr("Path Settings"));
    //resize(600, 175);
}

//! [2]
void PathSetting::browsedcmrawpath()
{
    QString directory = QFileDialog::getExistingDirectory(this, tr("Set Dicom RAW Path"),_dcmrawpath,
                                                          QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

    if (!directory.isEmpty()) {
      _dcmrawpathLineEdit->setText(directory);
    }
}

void PathSetting::browseoutputpath()
{
    QString directory = QFileDialog::getExistingDirectory(this, tr("Set Output Path"),_outputpath,
                                                          QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

    if (!directory.isEmpty()) {
        _outputpathLineEdit->setText(directory);
    }
}


void PathSetting::browsearchivepath()
{
    QString directory = QFileDialog::getExistingDirectory(this, tr("Set Archive Path"),_archivepath,
                                                          QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

    if (!directory.isEmpty()) {
       _archivepathLineEdit->setText(directory);
    }
}


QPushButton *PathSetting::createButton(const QString &text, const char *member)
{
    QPushButton *button = new QPushButton(text);
    connect(button, SIGNAL(clicked()), this, member);
    return button;
}

QLineEdit *PathSetting::createLineEdit(const QString &text)
{
    QLineEdit *LineEdit = new QLineEdit;
    LineEdit->setText(text);
    LineEdit->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    return LineEdit;
}

void PathSetting::readSettings()
{
    QSettings settings("Cuda IRGN -GUI", "IRGN Reconstruction - Read Dicom Raw / meas.dat - Write Dicom by Heigl");
    _dcmrawpath = settings.value("dicom-raw-path", QDir::currentPath()).toString();
    _measdatpath = settings.value("measdat-path", QDir::currentPath()).toString();
    _matlabpath = settings.value("matlab-path", QDir::currentPath()).toString();
    _outputpath = settings.value("output-path", QDir::currentPath()).toString();
    _archivepath = settings.value("archive-path", QDir::currentPath()).toString();
    _archiveactive = settings.value("archive-active", true).toBool();

    _dcmrawpath_start = _dcmrawpath;
    _measdatpath_start= _measdatpath;
    _outputpath_start= _outputpath;
    _matlabpath_start= _matlabpath;
    _archivepath_start= _archivepath;
    _archiveactive_start= _archiveactive;
}

void PathSetting::writeSettings()
{
    QSettings settings("Cuda IRGN -GUI", "IRGN Reconstruction - Read Dicom Raw / meas.dat - Write Dicom by Heigl");
    settings.setValue("dicom-raw-path", _dcmrawpath);
    settings.setValue("measdat-path", _measdatpath);
    settings.setValue("matlab-path", _matlabpath);
    settings.setValue("output-path", _outputpath);
    settings.setValue("archive-path", _archivepath);
    settings.setValue("archive-active", _archiveactive);
}

bool PathSetting::pathChanged()
{
  //qDebug()<<"\nasdfköj: "<<_archiveactive_start<<"  "<<_archiveactive;
  if ((_dcmrawpath_start == _dcmrawpath) && (_measdatpath_start == _measdatpath)
      && (_outputpath_start == _outputpath) && (_matlabpath_start == _matlabpath)
      && (_archivepath_start == _archivepath) && (_archiveactive_start == _archiveactive))
     return false;    //nothing changed
  else
    return true;

}


void PathSetting::cancelbutton()
{
  _dcmrawpathLineEdit->setText(_dcmrawpath);
  _outputpathLineEdit->setText(_outputpath);
  _archivepathLineEdit->setText(_archivepath);

  this->close();
}


void PathSetting::okbutton()
{
  _dcmrawpath = _dcmrawpathLineEdit->text();
  _outputpath = _outputpathLineEdit->text();
  _archivepath = _archivepathLineEdit->text();

  this->close();
}

void PathSetting::archiveactive()
{
  _archiveactive = _cb_archiveactive->isChecked();
  if(_archiveactive)
  {
    _archivepathLineEdit->setEnabled(true);
    _archivepathbrowseButton->setEnabled(true);
  }
  else
  {
    _archivepathLineEdit->setEnabled(false);
    _archivepathbrowseButton->setEnabled(false);
  }
}
