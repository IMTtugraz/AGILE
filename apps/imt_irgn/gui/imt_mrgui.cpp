#include "imt_mrgui.h"
#include "ui_imt_mrgui.h"
#include <QtGui>

// -- Agile - includes
#include "agile/io/file.hpp"
//#include "agile/io/readSiemensVD11.hpp"
#include "agile/io/dicom.hpp"
#include "agile/matrixhelper.h"
#include "agile/calc/irgn.hpp"


IMT_MRGui::IMT_MRGui(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::IMT_MRGui)
{
  // Initialize the first GPU
  GPU0_.allocateGPU(0);
  GPU0_.printInformation(std::cout);

  size_t free;
  size_t total;
  cuMemGetInfo(&free, &total);
  std::cout << "free memory: " << free / 1024 / 1024 << "mb, total memory: " << total / 1024 / 1024 << "mb" << std::endl;

  ui->setupUi(this);

  read_siemens = 0;
  fftobj_ = new agile::FFT<TType>();
  kspacefovobj_ = new agile::KSpaceFOV<TType>();

  open_irgnpara_window = new IRGN_para(); // Be sure to destroy you window somewhere
  open_irgnpara_window->hide();

  _pathsetting = new PathSetting(); // Be sure to destroy you window somewhere
  _pathsetting->hide();

  _specialsetting = new SpecialSetting(); // Be sure to destroy you window somewhere
  _specialsetting->hide();


  setStatusBar();
  setTitleText("");

  future = new QFuture<void>;
  watcher = new QFutureWatcher<void>;
  _cycleread_thread = new CycleRead;
  _cycleread_thread->set_max_acq_time(_specialsetting->get_autoload_timeout());


  _infotext="";

  QObject::connect(this,SIGNAL(sendInfo(QString, bool)),this,SLOT(set_Info(QString, bool)));
  QObject::connect(this,SIGNAL(sendWarning(QString)),this,SLOT(set_Warning(QString)));
  QObject::connect(this,SIGNAL(sendSaveFile(QString,QString)),this,SLOT(set_SaveFile(QString,QString)),Qt::BlockingQueuedConnection);

  QObject::connect(_cycleread_thread, SIGNAL(send_filenames(QStringList, QStringList)),
          this, SLOT(startauto(QStringList, QStringList)));

  _act_index = 0;
  _automatic_on = false;
  _write_file_index=0;

  ui->cb_savingtype->setCurrentIndex(2);
  ui->pb_automatic->setText("Auto Off");

  _postprocess = new agile::PostProcess<TType, TType_real>();


  QLabel* imageLabel = new QLabel(this);
  imageLabel->setBackgroundRole(QPalette::Base);
  imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
  imageLabel->setScaledContents(true);
  QImage image(":/images/wappen.png");
  imageLabel->setPixmap(QPixmap::fromImage(image));
  imageLabel->resize(imageLabel->pixmap()->size());
  imageLabel->setGeometry((this->width() - imageLabel->sizeHint().width()), 40, imageLabel->sizeHint().width(), imageLabel->sizeHint().height());
  imageLabel->show();
}


IMT_MRGui::~IMT_MRGui()
{
  delete fftobj_;
  fftobj_ = 0;
  delete kspacefovobj_;
  kspacefovobj_ = 0;
  delete read_siemens;
  read_siemens = 0;

  delete _cycleread_thread;

  delete open_irgnpara_window;
  delete _pathsetting;
  delete _specialsetting;
  delete ui;

  delete future;
  delete watcher;

  delete _postprocess;

}

void IMT_MRGui::set_Info(QString infotext, bool progress)
{
  //str_prog.setNum(progress);
  //qDebug()<<"\ninfotext: "<<infotext;
  //ui->text_info->setText(infotext+"  "+str_prog);
  ui->text_info->setText(infotext);

  QTextCursor text_cursor(ui->text_info->textCursor());
  text_cursor.movePosition(QTextCursor::End);
  ui->text_info->setTextCursor(text_cursor);

  if((progress == true))// || (_automatic_on == true))  //calculation is runnuing
  {
    ui->pushButton->setEnabled(false);
    ui->actionOpen_meas_dat->setEnabled(false);
    ui->actionRead_DicomRaw->setEnabled(false);
    ui->actionOpen_matlab_bin->setEnabled(false);
    ui->pb_outputpath->setEnabled(false);
    //progressbar to Wait
    ui->progressBar->setMaximum(0);
    ui->progressBar->setMinimum(0);
    ui->progressBar->repaint();
  }
  else
  {
    ui->pushButton->setEnabled(true);
    ui->actionOpen_meas_dat->setEnabled(true);
    ui->actionRead_DicomRaw->setEnabled(true);
    ui->actionOpen_matlab_bin->setEnabled(true);
    ui->pb_outputpath->setEnabled(true);
    ui->pb_automatic->setEnabled(true);
    //progressbar to Wait
    ui->progressBar->setMaximum(100);
    ui->progressBar->setMinimum(0);
    ui->progressBar->setValue(0);
    ui->progressBar->repaint();
  }

}



void IMT_MRGui::set_Warning(QString warningtext)
{
  QMessageBox::warning(this, tr("Warning"), "Load Data first!\n",
                             QMessageBox::Ok);
}


void IMT_MRGui::setTitleText(QString text)
{
    setWindowTitle("Cuda Powered IRGN - Gui  "+text);
}

void IMT_MRGui::setStatusBar()
{
    _statusBarLabelleft = new QLabel(this);
    //_statusBarLabelright = new QLabel(this);

    //_statusBarLabelleft->setText(tr("Dicom RawPath: ")+_pathsetting->get_dcmrawpath());
    _statusBarLabelleft->setText(tr("Output Path: ")+_pathsetting->get_outputpath());

    _statusBarLabelleft->setAlignment(Qt::AlignLeft);
    //_statusBarLabelright->setAlignment(Qt::AlignRight);

    statusBar()->addWidget(_statusBarLabelleft);
    //statusBar()->addWidget(_statusBarLabelright,1);
}

void IMT_MRGui::closeEvent(QCloseEvent *event)
{
  _cycleread_thread->cycle_stop();
  if (maybeSave()) {
      open_irgnpara_window->close();
      _pathsetting->close();
      event->accept();
  } else {
      event->ignore();
  }
}

bool IMT_MRGui::maybeSave()
{
   QMessageBox::StandardButton ret;

   if((!_pathsetting->pathChanged()) && (!_specialsetting->pathChanged()))
      return true;    //nothing changed

   ret = QMessageBox::warning(this, tr("Application"), tr("Do you want to save new Settings?\n"),
                               QMessageBox::Save | QMessageBox::Close | QMessageBox::Cancel);
   if (ret == QMessageBox::Cancel)
     return false;

   if (ret == QMessageBox::Save)
   {
     _pathsetting->writeSettings();
     _specialsetting->writeSettings();
     return true;
   }

    return true;
}

void IMT_MRGui::on_actionRead_DicomRaw_triggered()
{

    _readmatlab = false;

  //001_000004_000001.dcm
  //001_000004_000002.dcm

  _dicom=true;

  bool b_ok = false;
  int cntfile=0;

  QStringList filenames;
  QFileInfo file, fileoutname;
  QDir dcmrawdir;
  dcmrawdir.setPath(_pathsetting->get_dcmrawpath());

  QFileInfoList fList = dcmrawdir.entryInfoList(QDir::NoDotAndDotDot | QDir::Files , QDir::Name);

  if (fList.size() == 0)
    QMessageBox::warning(this, tr("Warning"), "No Files in Directory!\n"+file.absoluteFilePath(),
                               QMessageBox::Ok);


  for(int i=0;i < fList.size(); ++i)
  {
    file = fList.at(i);

    if( "dcm" == file.completeSuffix())
    {
      b_ok = true;
      filenames.push_back(file.absoluteFilePath());
      _filenames_to_archive.push_back(file.fileName());
      //qDebug()<<"\n filen: "<<file.fileName();

      if(cntfile == 0)    //only first dcm-file
        fileoutname = file;
      cntfile++;
    }
  }

  _outfilename = fileoutname.baseName()+"-"+QString::number(cntfile);   //set filename for outputfile

  setTitleText(_outfilename);

  if(b_ok == false)
    QMessageBox::warning(this, tr("Warning"), "No DicomRaw Files in Directory!\n"+file.absolutePath(),
                             QMessageBox::Ok);
  else
    get_dcmrawdata(filenames);

}

void IMT_MRGui::get_dcmrawdata(QStringList filenames)
{

  Uint8* pData;
  unsigned long length = 0;
  int error=0;

  std::string dcm_std_filename = _pathsetting->get_dcmrawpath().toStdString()+"/dcm-meas.dat";

  std::ofstream dcmmeas_file((char*)&dcm_std_filename[0], std::ofstream::binary);
  if (!dcmmeas_file.is_open()) {
    std::cerr << "not able to write to file: " << "dcm-meas.dat" << std::endl;
  }

  for(int i=0; i<filenames.size();++i)
  {
    error = _in_dicomfile.readdicom(filenames.at(i).toStdString(),dcmmeas_file,length, pData);
/*    qDebug()<<"\n length: "<<length;

    std::string str_studyinstuid;
    std::string str_seriesinstuid;
    agile::readdicom_info(filenames.at(i).toStdString(),str_studyinstuid, str_seriesinstuid);
    QString qstr_studyinstuid(str_studyinstuid.c_str());
    QString qstr_seriesinstuid(str_seriesinstuid.c_str());
    qDebug()<<"\n str_studyinstuid: "<<qstr_studyinstuid;
    qDebug()<<"\n str_seriesinstuid: "<<qstr_seriesinstuid;
*/
  }
  error = _in_dicomfile.readdicom_info(filenames.front().toStdString());

  dcmmeas_file.close();

  QString qstr(dcm_std_filename.c_str());
  if(error == 0)
    gui_readSiemens(qstr, true);
  else
    std::cerr<<"\n Error Reading Dicom Raw File(s)";

}

void IMT_MRGui::startauto(QStringList filepath, QStringList filename)
{
  _dicom=true;

//  qDebug()<<"\n filepath->size(): "<<filepath.size();
  if(!filepath.empty())
  {
    _cycleread_thread->cycle_stop();
    //qDebug()<<"\n filepath: "<<filepath;
    _filenames_to_archive = filename;

    _outfilename = filename.first()+"-"+QString::number(filename.size());   //set filename for outputfile
    setTitleText(_outfilename);

    get_dcmrawdata(filepath);
    createstartcalcthread();
  }
}


void IMT_MRGui::on_actionOpen_meas_dat_triggered()
{

  _readmatlab = false;

  QString path( QFileDialog::getOpenFileName( this, tr("Open File"),
                                             _pathsetting->get_measdatpath(),
                                             tr("*.dat") ) );
  QFileInfo pathinfo;
  QString filename;

  if ( ! path.isNull() ) {
    pathinfo.setFile(path);
    filename = pathinfo.fileName();

    _pathsetting->set_measdatpath(pathinfo.absolutePath());

    _outfilename = pathinfo.baseName();   //set filename for outputfile

    setTitleText(_outfilename);
    gui_readSiemens(path, false);
  }

  _dicom=false;

}

void IMT_MRGui::gui_readSiemens(QString path, bool dicom)
{
  // ========= ReadSiemens VD11 =========
  unsigned short int acq, sli, par, echo, pha, rep, set, seg;
  unsigned short int acq1=0, sli1=0, par1=0, echo1=0, pha1=0, rep1=0, set1=0, seg1=0;

  QString numrows;
  QString numcolumns;
  QString numcoils;
  QString str_acq, str_sli, str_par, str_echo, str_pha, str_rep, str_set, str_seg;

  if(read_siemens != 0)
  {
    delete read_siemens;
    read_siemens = 0;
  }
  read_siemens = new agile::ReadSiemensVD11(path.toStdString());

  int error = read_siemens->readfile(dicom);    //false for meas.dat - file structur
                                                //true for dicom - file structur

  error=0;  //siemens dicom-rawdata has no acq-end flag set ??!!!  ->ignore error
  if (error == 0)
  {
    _infotext ="";      //reset info message

    read_siemens->getRawdata_info(_num_rows_data, _num_columns_data, _num_coils_data, acq, sli, par, echo, pha, rep, set, seg);

    //qDebug() <<" "<<acq<<" "<< sli<<" "<< par<<" "<< echo<<" "<< pha<<" "<< rep<<" "<< set<<" "<< seg;

    numrows.setNum(_num_rows_data);
    numcolumns.setNum(_num_columns_data);
    numcoils.setNum(_num_coils_data);

    if(acq>0)
      acq1=acq+1;
    if(sli>0)
      sli1=sli+1;
    if(par>0)
      par1=par+1;
    if(echo>0)
      echo1=echo+1;
    if(pha>0)
      pha1=pha+1;
    if(rep>0)
      rep1=rep+1;
    if(set>0)
      set1=set+1;
    if(seg>0)
      seg1=seg+1;

    str_acq.setNum(acq1);
    str_sli.setNum(sli1);
    str_par.setNum(par1);
    str_echo.setNum(echo1);
    str_pha.setNum(pha1);
    str_rep.setNum(rep1);
    str_set.setNum(set1);
    str_seg.setNum(seg1);

    ui->label_coilsval->setText(numcoils);
    ui->label_rowsval->setText(numrows);
    ui->label_columnsval->setText(numcolumns);

    ui->label_acqval->setText(str_acq);
    ui->label_slival->setText(str_sli);
    ui->label_parval->setText(str_par);
    ui->label_echoval->setText(str_echo);
    ui->label_phaval->setText(str_pha);
    ui->label_repval->setText(str_rep);
    ui->label_setval->setText(str_set);
    ui->label_segval->setText(str_seg);

    QString str_int;
    ui->combo_sli_select->clear();
    ui->combo_acq_select->clear();
    ui->combo_echo_select->clear();
    ui->combo_par_select->clear();
    ui->combo_pha_select->clear();
    ui->combo_rep_select->clear();
    ui->combo_seg_select->clear();
    ui->combo_set_select->clear();

    if(sli > 0)
    {
      ui->combo_sli_select->addItem("All");
      for(int i=0;i<=sli;++i)
        ui->combo_sli_select->addItem(str_int.setNum(i));
    }
    if(acq > 0)
    {
      ui->combo_acq_select->addItem("All");
      for(int i=0;i<=acq;++i)
        ui->combo_acq_select->addItem(str_int.setNum(i));
    }
    if(par > 0)
    {
      ui->combo_par_select->addItem("All");
      for(int i=0;i<=par;++i)
        ui->combo_par_select->addItem(str_int.setNum(i));
    }
    if(echo > 0)
    {
      ui->combo_echo_select->addItem("All");
      for(int i=0;i<=echo;++i)
        ui->combo_echo_select->addItem(str_int.setNum(i));
    }
    if(pha > 0)
    {
      ui->combo_pha_select->addItem("All");
      for(int i=0;i<=pha;++i)
        ui->combo_pha_select->addItem(str_int.setNum(i));
    }
    if(rep > 0)
    {
      ui->combo_rep_select->addItem("All");
      for(int i=0;i<=rep;++i)
        ui->combo_rep_select->addItem(str_int.setNum(i));
    }
    if(set > 0)
    {
      ui->combo_set_select->addItem("All");
      for(int i=0;i<=set;++i)
        ui->combo_set_select->addItem(str_int.setNum(i));
    }
    if(seg > 0)
    {
      ui->combo_seg_select->addItem("All");
      for(int i=0;i<=seg;++i)
        ui->combo_seg_select->addItem(str_int.setNum(i));
    }

    statusBar()->showMessage(tr(" num: ")+numrows+tr("  ")+numcolumns+tr("  ")+numcoils);
  }
  else
  {
    std::cerr<<"\n Error: reading rawdata-file - error-number: "<<error;
  }
}

//void IMT_MRGui::on_actionSet_DicomRaw_Path_triggered()
//{
//  QString path;
//  path = QFileDialog::getExistingDirectory(this, tr("Set DicomRaw Path"),_dcmrawpath,
//                                           QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

//  if ( ! path.isNull() )
//    _dcmrawpath = path;

//  _statusBarLabelleft->setText(tr("Dicom RawPath: ")+_dcmrawpath);
//  statusBar()->repaint();

//}

void IMT_MRGui::on_pushButton_clicked()
{
  ui->pb_automatic->setEnabled(false);
  createstartcalcthread();
}

void IMT_MRGui::createstartcalcthread()
{
  emit sendInfo("Reading...",true);

  *future = QtConcurrent::run(this,&IMT_MRGui::readcalcwrite_thread);
  watcher->setFuture(*future);

  //watcher->waitForFinished();
}


void IMT_MRGui::generate_idx(agile::ReadSiemensVD11* read_siemens, std::vector<INDEX_select>& index_vector)
{
  INDEX_select index;
  unsigned short int acq_read, sli_read, par_read, echo_read, pha_read, rep_read, set_read, seg_read;
  unsigned int num_rows = 0;
  unsigned int num_columns = 0;
  unsigned int num_coils = 0;
  unsigned short int acq, sli, par, echo, pha, rep, set, seg;
  unsigned short int val_acq, val_sli, val_par, val_echo, val_pha, val_rep, val_set, val_seg;
  unsigned short int acq_select, sli_select, par_select, echo_select, pha_select, rep_select, set_select, seg_select;

  read_siemens->getRawdata_info(num_rows, num_columns, num_coils, acq_read, sli_read, par_read, echo_read, pha_read, rep_read, set_read, seg_read);

  sli_select = ( ui->combo_sli_select->currentIndex() == int(-1) ) ? 0 : ui->combo_sli_select->currentIndex()-1; //-1 for index All
  acq_select = ( ui->combo_acq_select->currentIndex() == int(-1) ) ? 0 : ui->combo_acq_select->currentIndex()-1; //-1 for index All
  par_select = ( ui->combo_par_select->currentIndex() == int(-1) ) ? 0 : ui->combo_par_select->currentIndex()-1; //-1 for index All
  echo_select = ( ui->combo_echo_select->currentIndex() == int(-1) ) ? 0 : ui->combo_echo_select->currentIndex()-1; //-1 for index All
  pha_select = ( ui->combo_pha_select->currentIndex() == int(-1) ) ? 0 : ui->combo_pha_select->currentIndex()-1; //-1 for index All
  rep_select = ( ui->combo_rep_select->currentIndex() == int(-1) ) ? 0 : ui->combo_rep_select->currentIndex()-1; //-1 for index All
  set_select = ( ui->combo_set_select->currentIndex() == int(-1) ) ? 0 : ui->combo_set_select->currentIndex()-1; //-1 for index All
  seg_select = ( ui->combo_seg_select->currentIndex() == int(-1) ) ? 0 : ui->combo_seg_select->currentIndex()-1; //-1 for index All

  if (ui->combo_sli_select->currentText() == "All")
    sli = sli_read + 1;
  else
    sli = 1;

  if (ui->combo_acq_select->currentText() == "All")
    acq = acq_read + 1;
  else
    acq = 1;

  if (ui->combo_par_select->currentText() == "All")
    par = par_read + 1;
  else
    par = 1;

  if (ui->combo_echo_select->currentText() == "All")
    echo = echo_read + 1;
  else
    echo = 1;

  if (ui->combo_pha_select->currentText() == "All")
    pha = pha_read + 1;
  else
    pha = 1;

  if (ui->combo_rep_select->currentText() == "All")
    rep = rep_read + 1;
  else
    rep = 1;

  if (ui->combo_set_select->currentText() == "All")
    set = set_read + 1;
  else
    set = 1;

  //unsigned all_select = acq_select * sli_select * par_select * echo_select * pha_select * rep_select * set_select * set_select;

  for(int i_sli=0; i_sli< sli; ++i_sli)
    for(int i_acq=0; i_acq< acq; ++i_acq)
      for(int i_par=0; i_par< par; ++i_par)
          for(int i_echo=0; i_echo< echo; ++i_echo)
              for(int i_pha=0; i_pha< pha; ++i_pha)
                  for(int i_rep=0; i_rep< rep; ++i_rep)
                      for(int i_set=0; i_set< set; ++i_set)
                          {
                            val_sli = (ui->combo_sli_select->currentText() == "All") ? i_sli : sli_select;
                            val_acq = (ui->combo_acq_select->currentText() == "All") ? i_acq : acq_select;
                            val_par = (ui->combo_par_select->currentText() == "All") ? i_par : par_select;
                            val_echo = (ui->combo_echo_select->currentText() == "All") ? i_echo : echo_select;
                            val_pha = (ui->combo_pha_select->currentText() == "All") ? i_pha : pha_select;
                            val_rep = (ui->combo_rep_select->currentText() == "All") ? i_rep : rep_select;
                            val_set = (ui->combo_set_select->currentText() == "All") ? i_set : set_select;
                            val_seg = (ui->combo_seg_select->currentText() == "All") ? seg_read+1 : seg_select;

                            index.sli = val_sli;
                            index.acq = val_acq;
                            index.par = val_par;
                            index.echo = val_echo;
                            index.pha = val_pha;
                            index.rep = val_rep;
                            index.set = val_set;
                            index.seg = val_seg;

                            index_vector.push_back(index);
                          }

}

void IMT_MRGui::readcalcwrite_thread()
{
  unsigned int num_rows = 0;
  unsigned int num_columns = 0;
  unsigned int num_coils = 0;
  std::vector<TType_data> matrix_read;
  std::vector<TType_image> image_data;
  QString name_index;
  unsigned short int acq, sli, par, echo, pha, rep, set, seg;
  unsigned short int acq_info, sli_info, par_info, echo_info, pha_info, rep_info, set_info, seg_info;
  std::vector<INDEX_select> index_vector;

  int error = 0;

  matrix_read.clear();        //clear all elements
  image_data.clear();         //clear all elements

  _seconds = time (NULL);

 if( ((read_siemens == 0) && !_readmatlab) || (_readmatlab && (_matlab_num_rows == 0)))
  {
    emit sendInfo("Load Data first!", false);
    emit sendWarning("Load Data first!");
  }
  else
  {
    if(!_readmatlab)
      read_siemens->getRawdata_info(num_rows, num_columns, num_coils, acq_info, sli_info, par_info, echo_info,
                                  pha_info, rep_info, set_info, seg_info);
    else
    {
      num_rows = _matlab_num_rows;
      num_columns = _matlab_num_columns;
      num_coils = _matlab_num_coils;
    }

    if(ui->cb_cropFOV->isChecked() == true)
    {
      _num_rows_calc = num_rows;
      _num_columns_calc = num_rows;
    }
    else
    {
      _num_rows_calc = num_rows;
      _num_columns_calc = num_columns;
    }

    if(!_readmatlab)
      generate_idx(read_siemens, index_vector);

    for(int i=0; i<_num_columns_calc*_num_rows_calc;  ++i)
      zero_vec_cpu_.push_back(TType_real(0));

    matrix_read.reserve(2*_num_columns_data*_num_rows_data*_num_coils_data);
    image_data.reserve(_num_columns_calc*_num_rows_calc);

    if(!_readmatlab)
    {
      for(int i=0; i<index_vector.size();++i)     //get threw all user-selected idx
      {
  //      endselect = setcheckindex(read_siemens, acq, sli, par, echo, pha, rep, set, seg );
        acq = index_vector[i].acq;
        sli = index_vector[i].sli;
        par = index_vector[i].par;
        echo = index_vector[i].echo;
        pha = index_vector[i].pha;
        rep = index_vector[i].rep;
        set = index_vector[i].set;
        seg = index_vector[i].seg;

        qDebug() <<"\n\nactual selection: "<<sli<<" "<< acq<<" "<< rep<<" "<< par<<" "<< echo<<" "<< pha<<" "<< set<<" "<< seg;

        if (ui->combo_seg_select->currentText() == "All")
        {
          for(int i=0; i<= seg_info;++i)
          {
            read_siemens->getRawdata(matrix_read,acq, sli, par, echo, pha, rep, set, i);
          }
        }
        else
        {
          read_siemens->getRawdata(matrix_read,acq, sli, par, echo, pha, rep, set, seg);
        }


        // TEST READ
        //agile::readMatrixFile3D("brain_151x256x6.bin", num_rows, num_columns, num_coils, matrix_read);
        //agile::readMatrixFile3D("braindata.bin", num_rows, num_columns, num_coils, matrix_read);

        QString str_acq, str_sli, str_par, str_echo, str_pha, str_rep, str_set, str_seg;
        str_acq.setNum(acq);
        str_sli.setNum(sli);
        str_par.setNum(par);
        str_echo.setNum(echo);
        str_pha.setNum(pha);
        str_rep.setNum(rep);
        str_set.setNum(set);
        str_seg.setNum(seg);

        _infotext.append("\n===================\n");
        //_strselection = "===================\n";
        _infotext.append("| "+str_sli+" - "+str_acq+" - "+str_par+" - "+str_echo+" - "+str_pha+" - "+str_rep+" - "
                             +str_set+" - "+str_seg+" |");
        _infotext.append("\n-------------------------------------");


        error = startcalculation(_num_rows_calc, _num_columns_calc, _num_coils_data, matrix_read, image_data);

        if(error == 0)
        {
          time_t timet=time (NULL);
          unsigned time = unsigned(timet) - unsigned(_seconds);
          QString str_time;
          str_time.setNum(time);
          _infotext.append("\nFinished in "+str_time+"sec\n");
        }

        if(error != 0)
          std::cerr<<"\n Calculation Error  - nb: "<<error;
        else
        {
          unsigned dummy=0;

          std::ostringstream filename;
          filename <<std::setfill ('0') << std::setw (4)<<_write_file_index
                   <<std::setfill ('0') << std::setw (2)<<dummy
                   <<std::setfill ('0') << std::setw (2)<<ui->cb_savingtype->currentIndex()+1
                   <<std::setfill ('0') << std::setw (2)<<sli
                   <<std::setfill ('0') << std::setw (2)<<echo
                   <<std::setfill ('0') << std::setw (2)<<dummy
                   <<std::setfill ('0') << std::setw (2)<<acq;

  //        filename << "_"<<std::setfill ('0') << std::setw (4)<<sli
  //                 <<std::setfill ('0') << std::setw (4)<<acq
  //                 <<std::setfill ('0') << std::setw (4)<<par
  //                 <<std::setfill ('0') << std::setw (4)<<echo
  //                 <<std::setfill ('0') << std::setw (4)<<pha
  //                 <<std::setfill ('0') << std::setw (4)<<rep
  //                 <<std::setfill ('0') << std::setw (4)<<set
  //                 <<std::setfill ('0') << std::setw (4)<<seg;
          name_index = (filename.str().c_str());

          writeoutfile(_num_rows_calc, _num_columns_calc, _num_coils_data, image_data, matrix_read,name_index);

          _write_file_index++;
        }
      }
    }
    else
    {
      QString path = _pathsetting->get_matlabpath() +"/"+ _outfilename+".bin";
      agile::readMatrixFile3D(path.toStdString().c_str(), num_rows, num_columns, num_coils, matrix_read);

            //qDebug()<<"\n path: "<<path<<"  size: "<<matrix_read.size();
      _infotext.append("\n===================");

      error = startcalculation(_num_rows_calc, _num_columns_calc, _num_coils_data, matrix_read, image_data);

      if(error == 0)
      {
        time_t timet=time (NULL);
        unsigned time = unsigned(timet) - unsigned(_seconds);
        QString str_time;
        str_time.setNum(time);
        _infotext.append("\nFinished in "+str_time+"sec\n");
      }

      if(error != 0)
        std::cerr<<"\n Calculation Error  - nb: "<<error;
      else
      {
        unsigned dummy=0;

        std::ostringstream filename;
        filename <<std::setfill ('0') << std::setw (4)<<_write_file_index
                 <<std::setfill ('0') << std::setw (2)<<dummy
                 <<std::setfill ('0') << std::setw (2)<<ui->cb_savingtype->currentIndex()+1
                 <<std::setfill ('0') << std::setw (2)<<sli
                 <<std::setfill ('0') << std::setw (2)<<echo
                 <<std::setfill ('0') << std::setw (2)<<dummy
                 <<std::setfill ('0') << std::setw (2)<<acq;

//        filename << "_"<<std::setfill ('0') << std::setw (4)<<sli
//                 <<std::setfill ('0') << std::setw (4)<<acq
//                 <<std::setfill ('0') << std::setw (4)<<par
//                 <<std::setfill ('0') << std::setw (4)<<echo
//                 <<std::setfill ('0') << std::setw (4)<<pha
//                 <<std::setfill ('0') << std::setw (4)<<rep
//                 <<std::setfill ('0') << std::setw (4)<<set
//                 <<std::setfill ('0') << std::setw (4)<<seg;
        name_index = (filename.str().c_str());

        writeoutfile(_num_rows_calc, _num_columns_calc, _num_coils_data, image_data, matrix_read,name_index);

        _write_file_index++;
      }

    }

    if(error==0)
    {
      if(_pathsetting->get_archiveactive())
        copytoarchive();

      emit sendInfo(_infotext,false);

      _act_index++;
    }
    else
    {
      _infotext.append("\nError\n");
      emit sendInfo(_infotext,false);
    }

    if(_automatic_on == true)
      _cycleread_thread->cycle_read(_pathsetting->get_dcmrawpath());

  }
}

void IMT_MRGui::set_SaveFile(QString title, QString extension)
{
  _str_outfilename="";

  _str_outfilename = QFileDialog::getSaveFileName(this,
           title, "",
           extension);
}

int IMT_MRGui::writeoutfile(unsigned num_rows, unsigned num_columns, unsigned num_coils, const std::vector<TType_image>& image_data,
                            const std::vector<TType_data>& matrix_read, QString name_index)
{
  QFileInfo file;

  QDateTime datum = QDateTime::currentDateTime() ;
  QString dat_str;
  dat_str = datum.toString("hh-mm-ss");

  agile::DICOM dicomfile;
  agile::GPUMatrix<TType_data> Coil_assigned;
  agile::GPUMatrix<TType_data>* Coil_fov = new agile::GPUMatrix<TType_data> [num_coils];
  std::vector<TType_data > matrix_read_fov, matrix_read_fov_help;
  agile::KSpaceFOV<TType_data>* kspacefovobj = new agile::KSpaceFOV<TType_data>();

  matrix_read_fov.resize(num_rows * num_columns * num_coils);
  matrix_read_fov.clear();

  if(ui->radioButton_rawdata->isChecked())    // no calculation - only RAWdata will be written
  {
    if(ui->cb_autofilename->isChecked())
      _str_outfilename = _pathsetting->get_outputpath()+"/"+name_index+".bin";
    else
      emit sendSaveFile(tr("Save file as"),tr("Rawdata (*.bin)"));

//    qDebug()<<"\n _str_outfilename: "<<_str_outfilename;
    if(ui->cb_cropFOV->isChecked())
    {
      for(unsigned i=0; i<num_coils;++i)
      {
        Coil_assigned.assignFromHost(_num_rows_data, _num_columns_data, &matrix_read[_num_rows_data*_num_columns_data*i]);
        kspacefovobj->genkspace_fov(Coil_assigned,Coil_fov[i]);
        Coil_fov[i].copyToHost(matrix_read_fov_help);

//        for(unsigned lauf=0;lauf<matrix_read_fov_help.size();++lauf)
//        {
//          matrix_read_fov.push_back(matrix_read_fov_help[lauf]);
//        }

        matrix_read_fov.insert(matrix_read_fov.end(),matrix_read_fov_help.begin() , matrix_read_fov_help.end());
      }
    agile::writeMatrixFile3D(_str_outfilename.toStdString().c_str(),num_rows,num_columns,num_coils,matrix_read_fov);   //write RAWdata
    }
    else
      agile::writeMatrixFile3D(_str_outfilename.toStdString().c_str(),num_rows,num_columns,num_coils,matrix_read);   //write RAWdata
  }
  else    // image_data will be written to dicom file
  {
    if(ui->cb_autofilename->isChecked())
      _str_outfilename = _pathsetting->get_outputpath()+"/"+name_index+".dcm";
    else
      emit sendSaveFile(tr("Save file as"),tr("Dicom (*.dcm);;Image (*.img)"));

    file.setFile(_str_outfilename);

    if( "dcm" == file.completeSuffix())
    {
      dicomfile.set_dicominfo(_in_dicomfile.get_dicominfo());
      dicomfile.gendicom(_str_outfilename.toStdString().c_str(), image_data, num_rows, num_columns);   //generate Dicom file
    }
    else //if *.img
      agile::writeVectorFile(_str_outfilename.toStdString().c_str(),image_data);
  }

  delete kspacefovobj;
  kspacefovobj = 0;
}

void IMT_MRGui::copytoarchive()
{
  QDateTime datum = QDateTime::currentDateTime() ;
  QString dat_str;
  dat_str = datum.toString("dd-MM-yyyy_hh-mm-ss");

  QString archive_path =_pathsetting->get_archivepath()+"/"+dat_str;
  QDir archive_dir;
  QFile filetocopy;

  if(_dicom)
  {

    archive_dir.mkpath(archive_path);

    for(int i=0; i<_filenames_to_archive.size();++i)
    {
      filetocopy.setFileName(_pathsetting->get_dcmrawpath()+"/"+_filenames_to_archive.at(i));
      filetocopy.copy(archive_path+"/"+_filenames_to_archive.at(i));
      filetocopy.remove();
    }
  }
}

//start calculation
int IMT_MRGui::startcalculation(unsigned num_rows,unsigned num_columns,unsigned num_coils,
                                const std::vector<TType_data>& matrix_read, std::vector<TType_image>& image_data)
{
  std::vector<TType_real> image_cpu_vec;
  //qDebug()<<"\n numrows: "<<num_rows<<"  "<<"num_columns: "<<num_columns<<"   matrix_read.size(): "<<matrix_read.size();
  //           <<"image_data_.size(): "<<image_data_.size();
  //qDebug() << "Hello from thread " << QThread::currentThread();

  std::vector<TType> matrix_calc(matrix_read.begin(), matrix_read.end());

  /*qDebug()<<"matrix_read.size(): "<<matrix_read.size();
  qDebug()<<"matrix_calc.size(): "<<matrix_calc.size();
  output("matrix_read[i]: ",1,10,matrix_read);
  output("matrix_calc[i]: ",1,10,matrix_read);*/


  if(matrix_calc.empty())
    return -1;               // return 1 - if matrix is empty

  agile::GPUMatrix<TType> Coil_assigned;
  agile::GPUMatrix<TType>* Coil;
  Coil = new agile::GPUMatrix<TType> [num_coils];
  agile::GPUMatrix<TType>* Coil_erg;
  Coil_erg = new agile::GPUMatrix<TType> [num_coils];

  _postprocess->set_size(num_rows,num_columns,num_coils);

  irgnpara = open_irgnpara_window->getIRGN_params();

  //------------------------------------ Generate FOV -------------------------------------

  for(unsigned i=0; i<num_coils;++i)
  {
    if(ui->cb_cropFOV->isChecked()) //generate FOV - yes:
    {
      Coil_assigned.assignFromHost(_num_rows_data, _num_columns_data, &matrix_calc[_num_rows_data*_num_columns_data*i]);
      kspacefovobj_->genkspace_fov(Coil_assigned,Coil[i]);
    }
    else                            //generate FOV - no:
    {
      Coil[i].assignFromHost(_num_rows_data, _num_columns_data, &matrix_calc[_num_rows_data*_num_columns_data*i]);
    }
    Coil_erg[i].resize(Coil[i].getNumRows(),Coil[i].getNumColumns());
  }

  //qDebug()<<"\n Coil[i].getNumRows(): "<<Coil[0].getNumRows()<<"    - Coil[i].getNumColumns():"<<Coil[0].getNumColumns();

  //------------------------------------ Undersampling Simulation--------------------------

  //bool simu_active = true;
  IRGN_Simul irgn_simul = open_irgnpara_window->getIRGN_simul();
  if(irgn_simul.active)
  {
    _factor_x = irgn_simul.factor_x;
    _factor_y = irgn_simul.factor_y;
    _ref_lines = irgn_simul.ref_lines;

    simu_undersampling(Coil, num_coils);
  }

  //--------------------------- Start IRGN, IFFT, RAW Calculation  -------------------------
  if(ui->radioButton_irgntgv->isChecked())
  {
    _infotext.append("\n IRGN TGV selected \n Calculating...\n");
    emit sendInfo(_infotext,true);

    IRGNobj = new agile::TGVSolve<TType,TType_real>(Coil,num_coils, irgnpara);
    irgncalc(IRGNobj, image_data);

    delete IRGNobj;
  }
  else if(ui->radioButton_irgntv->isChecked())
  {
    _infotext.append("\n IRGN TV selected \n Calculating...\n");
    emit sendInfo(_infotext,true);
    IRGNobj = new agile::TVSolve<TType,TType_real>(Coil,num_coils, irgnpara);
    irgncalc(IRGNobj, image_data);

    delete IRGNobj;
  }
  else if(ui->radioButton_irgnl2->isChecked())
  {
    _infotext.append("\n IRGN L2 selected \n Calculating...\n");
    emit sendInfo(_infotext,true);
    IRGNobj = new agile::L2Solve<TType,TType_real>(Coil,num_coils, irgnpara);
    irgncalc(IRGNobj, image_data);

    delete IRGNobj;
  }
  else if(ui->radioButton_ifft->isChecked())
  {
    _infotext.append("\nifft selected\n");
    agile::GPUMatrix<TType_real> erg_val(num_rows,num_columns,&zero_vec_cpu_[0]);

    int error = 0;

    fftobj_->setfftplan(num_rows,num_columns);

    for(unsigned i=0;i<num_coils;++i)
    {
      fftobj_->calc_pattern(Coil[0]);
      error = fftobj_->CenterdIFFT(Coil[i],Coil_erg[i]);       //ifft each Coil Matrix

      if(error != 0)
      {
        qDebug()<<"\n!err: "<<error;
        return error;
      }
    }

    calc_postprocess(&Coil_erg[0],erg_val);   
    erg_val.copyToHost(image_cpu_vec);                  //copy solution back to host
    std::vector<TType_image> img(image_cpu_vec.begin(), image_cpu_vec.end());

    image_data.assign(img.begin(), img.end());
  }
  else if(ui->radioButton_rawdata->isChecked())
  {
    _infotext.append("\nRawdata selected \n");
    //qDebug()<<"\n Rawdata slected";
  }

  _infotext.prepend(_strselection);

  emit sendInfo(_infotext,true);

  delete[] Coil;
  delete[] Coil_erg;

  Coil=0;
  Coil_erg=0;

  return 0;
}

//------------- Undersampling Simulation -------------
void IMT_MRGui::simu_undersampling(agile::GPUMatrix<TType>* Coil, unsigned num_coils)
{
  unsigned num_rows = Coil->getNumRows();
  unsigned num_columns = Coil->getNumColumns();

  agile::GPUMatrix<TType_real> Simul_matrix;

  std::vector<TType_real> simul_vector;
  std::vector<TType_real> simul_vector_2;

  simul_vector.clear();
  for(int r=0; r<num_rows; ++r)
  {
    if( ((r % _factor_y) == 0) || ((r >= (num_rows-_ref_lines)/2) && (r < (num_rows+_ref_lines)/2)) )
      for(int c=0; c<num_columns; ++c)
      {
        if( ((c % _factor_x) == 0) || ((c >= (num_columns-_ref_lines)/2) && (c < (num_columns+_ref_lines)/2)) )
          simul_vector.push_back(1);
        else
          simul_vector.push_back(0);
      }
    else
      for(int c=0; c<num_columns; ++c)
      {
        simul_vector.push_back(0);
      }
  }

  unsigned num_rows_2;
  unsigned num_columns_2;
  bool fileexist = true;

  QString path="simul_pattern.bin";

  std::ifstream mat_file(path.toStdString().c_str(), std::ifstream::binary);
  if (!mat_file.is_open())  //check if file exists
    fileexist=false;
  else
    agile::readMatrixFile(path.toStdString().c_str(), num_rows_2, num_columns_2, simul_vector_2);

  if(fileexist && ((num_rows_2 != num_rows)||(num_columns_2 != num_columns)))
  {
    Simul_matrix.assignFromHost(num_rows, num_columns, &zero_vec_cpu_[0]);
    std::cerr<<"\n Simulation Pattern Matrix Dimension problem - check simul_pattern.bin";
  }
  else if(fileexist)
  {
    Simul_matrix.assignFromHost(num_rows, num_columns, &simul_vector_2[0]);
    qDebug()<<"\n Simulation Pattern Matrix - simul_pattern.bin - assigned";
  }
  else // assign calculated simul_vector
    Simul_matrix.assignFromHost(num_rows, num_columns, &simul_vector[0]);

  for(int i=0; i<num_coils; ++i)
  {
    agile::multiplyElementwise(Simul_matrix, Coil[i], Coil[i]);
  }


}


void IMT_MRGui::irgncalc(agile::IRGN<TType,TType_real>* IRGNobj, std::vector<TType_image>& image_data)
{
  agile::GPUMatrix<TType>* us_matrix;
  agile::GPUMatrix<TType_real> erg_val;
  std::vector<TType_real> image_cpu_vec;

  QString resnormtext;
  QString breaktext;
  QString str_iter;
  QString str_iterall="";
  int runto= 0;

  IRGNobj->HighFreqPenalty();
  IRGNobj->Normalize();
  IRGNobj->Iteration();
  IRGNobj->Postprocess();


  us_matrix = IRGNobj->get_us_mat();

  _postprocess->set_size(IRGNobj->get_coil()->getNumRows(),IRGNobj->get_coil()->getNumColumns());
  calc_postprocess(&us_matrix[0],erg_val);
  erg_val.copyToHost(image_cpu_vec);
  image_data.assign(image_cpu_vec.begin(), image_cpu_vec.end());

//      qDebug()<<"\n image_data - size: "<<image_data_.size();
//      agile::DICOM dicomfile;
//      dicomfile.set_dicominfo(_in_dicomfile.get_dicominfo());
//      dicomfile.gendicom("test.dcm", image_data_, 256, 256);   //generate Dicom file


  str_iterall.append("\n Residual Norm: ");

  if(IRGNobj->get_nr_k().back() == TType_real(-1))   //if break happend
    runto = IRGNobj->get_nr_k().size()-1;
  else
    runto = IRGNobj->get_nr_k().size();

  for(int i=1;i<runto;++i)
  {
    str_iter.setNum(i);
    resnormtext.setNum(IRGNobj->get_nr_k()[i]);
    str_iterall.append("\n Iteration "+str_iter+": "+resnormtext);
  }
  breaktext.setNum(IRGNobj->get_nr_k().size()-2);

  if(IRGNobj->get_nr_k().back() == TType_real(-1))
    str_iterall.append("\n  Break at Iteration: "+breaktext);

  _infotext.append(str_iterall);


}



//void IMT_MRGui::on_actionSet_Output_Path_triggered()
//{
//  QString path;
//  path = QFileDialog::getExistingDirectory(this, tr("Set Output Path"),_outputpath,
//                                                   QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

//  if ( ! path.isNull() )
//    _outputpath = path;

//  _statusBarLabelright->setText(tr("Output Path: ")+_outputpath);
//  statusBar()->repaint();

//}

void IMT_MRGui::on_actionParameter_triggered()
{
  open_irgnpara_window->show();
}

void IMT_MRGui::on_pb_outputpath_clicked()
{
  QString path;
  path = QFileDialog::getOpenFileName(this, tr("Open Selected File"), _pathsetting->get_outputpath(), "Data (*.dcm *.raw *.img)");
  if ( ! path.isNull() )
    QDesktopServices::openUrl(QUrl::fromLocalFile(path));
}

void IMT_MRGui::on_pb_automatic_clicked()
{
  _readmatlab = false;
  if(_automatic_on)
  {   //stop automatic
    _automatic_on =false;
    ui->pb_automatic->setStyleSheet("* { background-color: rgb(200,200,200) }");
    ui->pb_automatic->setText("Auto Off");

    ui->pushButton->setEnabled(true);

    _cycleread_thread->cycle_stop();
  }
  else
  {   //start automatic
    _automatic_on =true;
    ui->pb_automatic->setStyleSheet("* { background-color: rgb(125,255,100) }");
    ui->pb_automatic->setText("Auto On");

    ui->pushButton->setEnabled(false);

    _cycleread_thread->set_max_acq_time(_specialsetting->get_autoload_timeout());
    _cycleread_thread->cycle_read(_pathsetting->get_dcmrawpath());

  }

}


void IMT_MRGui::calc_postprocess(const agile::GPUMatrix<TType>* in_mat, agile::GPUMatrix<TType_real>& out_mat)
{
  switch(ui->cb_savingtype->currentIndex())
  {
  //real:
   case 0: _postprocess->calc_real(in_mat,out_mat);
           break;
  //imag:
   case 1: _postprocess->calc_imag(in_mat,out_mat);
           break;
  //abs:
   case 2: _postprocess->calc_abs(in_mat,out_mat);
           break;
  //phase:
   case 3: _postprocess->calc_phase(in_mat,out_mat);
           break;
  }
}

void IMT_MRGui::on_actionOpen_matlab_bin_triggered()
{

  QString path( QFileDialog::getOpenFileName( this, tr("Open File"),
                                             _pathsetting->get_matlabpath(),
                                             tr("*.bin") ) );
  QFileInfo pathinfo;
  QString filename;

  if ( ! path.isNull() ) {
    pathinfo.setFile(path);
    filename = pathinfo.fileName();

    _pathsetting->set_matlabpath(pathinfo.absolutePath());

    _outfilename = pathinfo.baseName();   //set filename for outputfile

    setTitleText(_outfilename);
    gui_readMatlab(path);
  }

  _dicom=false;

}


void IMT_MRGui::gui_readMatlab(QString path)
{
  QString numrows;
  QString numcolumns;
  QString numcoils;

  std::vector<TType_data> matrix_read;

  bool ok = agile::readMatrixFile3D(path.toStdString().c_str(), _num_rows_data, _num_columns_data, _num_coils_data, matrix_read);

  if(ok)
  {
    ui->label_acqval->setText("0");
    ui->label_slival->setText("0");
    ui->label_parval->setText("0");
    ui->label_echoval->setText("0");
    ui->label_phaval->setText("0");
    ui->label_repval->setText("0");
    ui->label_setval->setText("0");
    ui->label_segval->setText("0");

    ui->combo_sli_select->clear();
    ui->combo_acq_select->clear();
    ui->combo_echo_select->clear();
    ui->combo_par_select->clear();
    ui->combo_pha_select->clear();
    ui->combo_rep_select->clear();
    ui->combo_seg_select->clear();
    ui->combo_set_select->clear();

    numrows.setNum(_num_rows_data);
    numcolumns.setNum(_num_columns_data);
    numcoils.setNum(_num_coils_data);

    _matlab_num_rows = _num_rows_data;
    _matlab_num_columns = _num_columns_data;
    _matlab_num_coils = _num_coils_data;


    ui->label_coilsval->setText(numcoils);
    ui->label_rowsval->setText(numrows);
    ui->label_columnsval->setText(numcolumns);

    statusBar()->showMessage(tr(" num: ")+numrows+tr("  ")+numcolumns+tr("  ")+numcoils);

    _readmatlab = true;

  }
  else
  {
    std::cerr<<"\n Error: reading rawdata-file - error-number: "<<ok;
  }
}

void IMT_MRGui::on_actionAuto_Archive_triggered()
{
  _pathsetting->show();
}

void IMT_MRGui::on_actionSpecial_triggered()
{
  _specialsetting->show();
}
