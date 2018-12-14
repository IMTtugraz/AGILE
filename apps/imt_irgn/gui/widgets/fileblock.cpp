
#include <QtGui>

#include "fileblock.h"

FileBlock::FileBlock()//QObject *parent)
  //: QObject(parent)
{
  _studyinstuid ="";
  _seriesinstuid ="";
  _lastfile=false;
  _cnt_lastfile=0;
}


FileBlock::~FileBlock()
{

}

bool FileInfo::operator==(const FileInfo &rhs) const
{
  if( this->_filename != rhs._filename)   return false;
  if( this->_filepath != rhs._filepath)   return false;

  return true;
}

//std::ostream &FileInfo::operator<<(std::ostream &output, const FileInfo &fileinfo)
//{
//  output << 'filename: ' << fileinfo.filename << 'filepath: ' << fileinfo.filepath
//            << 'filesize: ' << fileinfo.filesize;
//  return output;
//}
QDebug operator<<(QDebug dbg, const FileInfo &fileinfo)
 {
     dbg << "\n| filename: " << fileinfo._filename << "  filepath: " << fileinfo._filepath
         << "  filesize: "<< fileinfo._filesize;
     return dbg.maybeSpace();
 }

//std::ostream &FileBlock::operator<<(std::ostream &output, const FileBlock &fileblock)
//{
//  output << fileblock._seriesinstuid << '\n ' << fileblock._studyinstuid << '\n ' << fileblock._files
//         << '\n _lastfile:' <<fileblock._lastfile;
//  return output;
//}

QDebug operator<<(QDebug dbg, const FileBlock &fileblock)
 {
     dbg << fileblock._seriesinstuid << "\n " << fileblock._studyinstuid << "\n  "<< fileblock._files
         << "\n _lastfile:" <<fileblock._lastfile;
     return dbg.maybeSpace();
 }

bool FileBlock::set_file(QFileInfo file)
{
  FileInfo fileinfo;

  agile::DICOM dicomfile;
  dicomfile.readdicom_info(file.filePath().toStdString());
  agile::DicomInfo dcminfo=dicomfile.get_dicominfo();

  if((_studyinstuid == "")&&(_seriesinstuid == ""))
  {
    // fist file for this block
    _studyinstuid = dcminfo.str_studyinstuid.c_str();
    _seriesinstuid = dcminfo.str_seriesinstuid.c_str();
  }

  if((_studyinstuid == dcminfo.str_studyinstuid.c_str())&&(_seriesinstuid == dcminfo.str_seriesinstuid.c_str()))
  {
    // new file for this block
    fileinfo.set_filename(file.fileName());
    fileinfo.set_filepath(file.filePath());
    fileinfo.set_filesize(unsigned(file.size()));

    if(!_files.contains(fileinfo))  //file not in list
    {
      _files.push_back(fileinfo);
      this->set_cnt_lastfile(0);
    }
  }
  else
    return false; //if file not for this block

  for(int i=1;i<_files.size();++i)
  {
    if(_files[i].get_filesize() < _files[i-1].get_filesize()) //if appended file is smaller than previous
      _lastfile=true;
  }

  return true; //if file was added to this block
}
