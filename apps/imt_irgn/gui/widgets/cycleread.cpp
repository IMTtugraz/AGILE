#include <QtGui>

#include "cycleread.h"

CycleRead::CycleRead(QObject *parent)
    : QThread(parent)
{
    restart = false;
    abort = false;

    read_siemens = new agile::ReadSiemensVD11();

    _all_filenames.clear();
}

CycleRead::~CycleRead()
{
    mutex.lock();
    abort = true;
    condition.wakeOne();
    mutex.unlock();

    delete read_siemens;
    read_siemens = 0;

    wait();
}

void CycleRead::cycle_read(QString dcmrawpath)
{
  _dcmrawpath = dcmrawpath;

    QMutexLocker locker(&mutex);

    if (!isRunning()) {
        abort = false;
        start(LowPriority);
    } else {
        abort = false;
        condition.wakeOne();
    }
}

void CycleRead::cycle_stop()
{
  abort = true;
}

void CycleRead::run()
{
  time_t seconds;


  forever {
      mutex.lock();
      QString dcmrawpath = this->_dcmrawpath;
      mutex.unlock();

      read_dcmrawdatapath(dcmrawpath);

//      qDebug()<<"\n _fileblocks.size(): "<< _fileblocks.size();
      for(int i=0; i<_fileblocks.size();++i)
      {
        if(_fileblocks[i].get_lastfile() == false)
        { // set file status to lastfile if max_acq_time is reached
          seconds = time (NULL);
          if(_fileblocks[i].get_cnt_lastfile() == 0)  //first time here
          {
            _fileblocks[i].set_cnt_lastfile(unsigned(seconds));
          }
          if((_fileblocks[i].get_cnt_lastfile()+_max_acq_time) <= unsigned(seconds))
          {
            _fileblocks[i].set_lastfile(true);
          }
        }

//        qDebug()<<"\n fileblocks: "<<_fileblocks[i];

        if(_fileblocks[i].get_lastfile() == true)
        { // start calculation

          emit send_filenames(_fileblocks[i].get_filepath(), _fileblocks[i].get_filename());

          _fileblocks.removeAt(i);
          break;

        }

      }


      //emit send_filenames(_p_file_pool);

      QThread::sleep(int(1));   //sleep for 1 sec
      if (abort)
        return;

  }
}


void CycleRead::read_dcmrawdatapath(QString dcmrawpath)
{
  QFileInfo file;
  QDir dcmrawdir;

  dcmrawdir.setPath(dcmrawpath);


  QFileInfoList fList = dcmrawdir.entryInfoList(QDir::NoDotAndDotDot | QDir::Files , QDir::Name);

  if (!fList.isEmpty())
//    QMessageBox::warning(this, tr("Warning"), "No Files in Directory!\n"+file.absoluteFilePath(),
//                               QMessageBox::Ok);

  {
    for(int i=0;i < fList.size(); ++i)
    {
      file = fList.at(i);

      if(("dcm" == file.completeSuffix()) && (!_all_filenames.contains(file.filePath())))
      {
        _all_filenames.push_back(file.filePath());

        if(_fileblocks.isEmpty())
        { //if its the fist file - and no fileblocks exist - create first
          FileBlock fileblock;
          fileblock.set_file(file);
          _fileblocks.push_back(fileblock);
        }
        else
        {
          if(!_fileblocks.last().set_file(file))
          { //if actual file is not for actual block
            FileBlock fileblock;
            fileblock.set_file(file);
            _fileblocks.push_back(fileblock);  //create new block at the back
          }
        }
      }
    }
  }
}
