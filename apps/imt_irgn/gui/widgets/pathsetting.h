#ifndef PATHSETTING_H
#define PATHSETTING_H

#include <QDialog>
#include <QDir>
#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QString>
#include <QFileDialog>
#include <QSettings>
#include <QDebug>
#include <QSpacerItem>
#include <QCheckBox>

QT_BEGIN_NAMESPACE
class QLineEdit;
class QLabel;
class QPushButton;
class QTableWidget;
class QTableWidgetItem;
QT_END_NAMESPACE

class PathSetting : public QDialog
{
    Q_OBJECT

public:
    PathSetting(QWidget *parent = 0);

    void readSettings();
    void writeSettings();
    bool pathChanged();

    QString get_dcmrawpath()
    {
      return _dcmrawpath;
    }

    QString get_outputpath()
    {
      return _outputpath;
    }

    QString get_archivepath()
    {
      return _archivepath;
    }

    QString get_measdatpath()
    {
      return _measdatpath;
    }

    QString get_matlabpath()
    {
      return _matlabpath;
    }

    bool get_archiveactive()
    {
      return _archiveactive;
    }


    void set_measdatpath(QString measdatpath)
    {
      _measdatpath=measdatpath;
    }

    void set_matlabpath(QString matlabpath)
    {
      _matlabpath=matlabpath;
    }

private slots:
    void browsedcmrawpath();
    void browseoutputpath();
    void browsearchivepath();
    void okbutton();
    void cancelbutton();
    void archiveactive();

private:
    QPushButton *createButton(const QString &text, const char *member);
    QLineEdit *createLineEdit(const QString &text = QString());

    QLineEdit *_dcmrawpathLineEdit;
    QLineEdit *_outputpathLineEdit;
    QLineEdit *_archivepathLineEdit;
    QLabel *_dcmrawpathLabel;
    QLabel *_outputpathLabel;
    QPushButton *_dcmrawpathbrowseButton;
    QPushButton *_outputpathbrowseButton;
    QPushButton *_archivepathbrowseButton;

    QCheckBox* _cb_archiveactive;
    bool _archiveactive;
    bool _archiveactive_start;

    QSpacerItem *_spacer_line;

    QString _dcmrawpath;
    QString _measdatpath;
    QString _outputpath;
    QString _matlabpath;
    QString _archivepath;
    QString _dcmrawpath_start;
    QString _measdatpath_start;
    QString _outputpath_start;
    QString _matlabpath_start;
    QString _archivepath_start;
};

#endif  //PATHSETTING_H
