#ifndef FILEBLOCK_H
#define FILEBLOCK_H

#include <QObject>
#include <QFileInfo>
#include <QString>
#include <QStringList>

#include "agile/io/readSiemensVD11.hpp"
#include "agile/io/dicom.hpp"

class FileInfo// : public QObject
{
  //Q_OBJECT

  friend QDebug operator<<(QDebug dbg, const FileInfo &fileinfo);

  public:
    FileInfo()//QObject *parent = 0)
    {}
   ~FileInfo()
    {}

    bool operator== ( const FileInfo & rhs ) const;
    //friend std::ostream &operator<<(std::ostream &, const FileInfo &);

    void set_filepath(QString filepath)
    {
      this->_filepath = filepath;
    }

    void set_filename(QString filename)
    {
      this->_filename = filename;
    }

    void set_filesize(unsigned filesize)
    {
      this->_filesize = filesize;
    }

    QString get_filepath()
    {
      return _filepath;
    }

    QString get_filename()
    {
      return _filename;
    }

    unsigned get_filesize()
    {
      return _filesize;
    }
  private:

    QString _filepath;
    QString _filename;
    unsigned _filesize;
};


class FileBlock// : public QObject
{
  //Q_OBJECT

  friend QDebug operator<<(QDebug dbg, const FileBlock &fileblock);

  public:
     FileBlock();//QObject *parent = 0);
    ~FileBlock();

     bool set_file(QFileInfo file);

     bool get_lastfile()
     {
       return _lastfile;
     }

     void set_lastfile(bool lastfile)
     {
      _lastfile = lastfile;
     }

     unsigned get_cnt_lastfile()
     {
        return _cnt_lastfile;
     }

     void set_cnt_lastfile(unsigned cnt_lastfile)
     {
        _cnt_lastfile = cnt_lastfile;
     }

     QStringList get_filepath()
     {
       QStringList filepath;

       for(int i=0; i<_files.size();++i)
       {
         filepath.push_back(_files[i].get_filepath());
       }
       return filepath;
     }

     QStringList get_filename()
     {
       QStringList filename;

       for(int i=0; i<_files.size();++i)
       {
         filename.push_back(_files[i].get_filename());
       }
       return filename;
     }

     //std::ostream &operator<<(std::ostream &, const FileBlock &);



  private:

    QString _studyinstuid;
    QString _seriesinstuid;
    QList<FileInfo> _files;

    bool _lastfile;
    unsigned _cnt_lastfile;

};

#endif // FILEBLOCK_H
