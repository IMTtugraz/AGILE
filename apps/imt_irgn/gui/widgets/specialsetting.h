#ifndef SPECIALSETTING_H
#define SPECIALSETTING_H

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

class SpecialSetting : public QDialog
{
    Q_OBJECT

public:
    SpecialSetting(QWidget *parent = 0);

    void readSettings();
    void writeSettings();
    bool pathChanged();

    unsigned get_autoload_timeout()
    {
      return _autoload_timeout;
    }

private slots:
    void okbutton();
    void cancelbutton();

private:
    QPushButton *createButton(const QString &text, const char *member);
    QLineEdit *createLineEdit(const QString &text = QString());

    QLineEdit *_le_autoload_timeout;
    QLabel *_lb_autoload_timeout;

    QSpacerItem *_spacer_line;

    unsigned _autoload_timeout;
    unsigned _autoload_timeout_start;
};

#endif  //SPECIALSETTING_H
