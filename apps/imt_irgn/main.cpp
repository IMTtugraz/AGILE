// -- QtGui - includes
#include <QtGui/QApplication>
#include "gui/imt_mrgui.h"

int main(int argc, char *argv[])
{
  Q_INIT_RESOURCE(images);

  // start GUI Application
  QApplication a(argc, argv);
  IMT_MRGui w;
  //w.resize(400,640);
  w.setFixedSize(705,570);
  w.show();

  QIcon icon(":/images/icon24x24.png");
  w.setWindowIcon(icon);

  a.exec();


  // ===============================================================================
/*
  unsigned int num_rows = 0;
  unsigned int num_columns = 0;
  unsigned int num_coils = 0;
  std::vector<TType> matrix_read;

  Uint8* pData;
  unsigned long length = 0;

  std::ofstream dcmmeas_file("dcm-meas.dat", std::ofstream::binary);
  if (!dcmmeas_file.is_open()) {
    std::cerr << "not able to write to file: " << "dcm-meas.dat" << std::endl;
    return false;
  }

  readdicom("001_000004_000001.dcm",pData, length, dcmmeas_file);
  readdicom("001_000004_000002.dcm",pData, length, dcmmeas_file);
  readdicom("001_000005_000001.dcm",pData, length, dcmmeas_file);
  readdicom("001_000005_000002.dcm",pData, length, dcmmeas_file);
  readdicom("001_000005_000003.dcm",pData, length, dcmmeas_file);
  readdicom("001_000005_000004.dcm",pData, length, dcmmeas_file);


  dcmmeas_file.close();


  // ========= ReadSiemens =========

  unsigned short int acq, sli, par, echo, pha, rep, set, seg;

  ReadSiemensVD11* read_siemens;
  //read_siemens = new ReadSiemens("meas_MID253_t2_tse_384_singleSlice_triple_FID36961.dat");
  //read_siemens = new ReadSiemensVD11("meas_MID131_localizer_FID7889.dat");

  //read_siemens = new ReadSiemensVD11("meas_MID93_t1_fl2d_tra_FID11734.dat");
  //read_siemens = new ReadSiemensVD11("meas_MID94_t2_tse_tra_FID11735.dat");

  //read_siemens = new ReadSiemensVD11("meas_MID95_t1_fl2d_tra_2averages_FID11736.dat");
  //read_siemens = new ReadSiemensVD11("meas_MID131_localizer_FID7889.dat");
  read_siemens = new ReadSiemensVD11("dcm-meas.dat");

  int error = read_siemens->readfile(true);   //true for dicom-read
  if (error == 0)
  {
    read_siemens->getRawdata_info(acq, sli, par, echo, pha, rep, set, seg);

    std::cout<<"\n acq: "<<acq<<"  sli: "<<sli<<"  par: "<<par<<"  echo: "<<echo<<"  pha: "<<pha<<"  rep: "<<rep<<"  set: "<<set<<"  seg: "<<seg;

    read_siemens->getRawdata(num_rows, num_columns, num_coils, matrix_read,0,1,0,0,0,0,0,0);

    std::cout<<"\n num: "<<num_rows<<"  "<<num_columns<<"  "<<num_coils;

    agile::writeMatrixFile3D("measrawimage.dat",num_rows,num_columns,num_coils,matrix_read);


  }
  //output("\n rawdata: ",10,10,matrix_read);

  //agile::readMatrixFile3D("braindata.bin", num_rows, num_columns, num_coils, matrix_read);

  delete read_siemens;
  read_siemens = 0;


*/
  return 0;
}
