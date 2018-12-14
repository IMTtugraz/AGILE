#ifndef CYCLEREAD_H
#define CYCLEREAD_H


#include <QMutex>
#include <QSize>
#include <QThread>
#include <QWaitCondition>
#include <QStringList>

#include "agile/io/readSiemensVD11.hpp"
#include "agile/io/dicom.hpp"
#include "fileblock.h"
#include <time.h>

//const unsigned MAX_ACQ_TIME = 20;   //in sec

class CycleRead : public QThread
{
    Q_OBJECT

public:
    CycleRead(QObject *parent = 0);
    ~CycleRead();

    void cycle_read(QString dcmrawpath);
    void cycle_stop();

    void set_max_acq_time(unsigned max_acq_time)
    {
      _max_acq_time = max_acq_time;
    }

signals:
    void send_filenames(QStringList filepath, QStringList filename);

protected:
    void run();

private:

    void read_dcmrawdatapath(QString dcmrawpath);

    QMutex mutex;
    QWaitCondition condition;

    QList<FileBlock> _fileblocks;
    QStringList _all_filenames;

    bool restart;
    bool abort;

    QString _dcmrawpath;
    agile::ReadSiemensVD11* read_siemens;

    unsigned _max_acq_time;
};

#endif // CYCLEREAD_H
