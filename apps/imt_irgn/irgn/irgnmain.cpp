// Copyright (C) 2010-2011 Institute of Medical Engineering,
// Graz University of Technology
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses>.

// $Id: gpu_matrix_test.cpp 476 2011-06-16 08:54:14Z freiberger $

// We have to include the headers for the environment, for the GPU vector and
// for the GPU matrix.


//#include "agile/operator/irgn_all.hpp"
#include "irgn_all.hpp"
//#include "irgn.hpp"
//#include "l2solve.hpp"
//#include "tvsolve.hpp"
//#include "tgvsolve.hpp"


#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_matrix.hpp"
#include "agile/gpu_vector.hpp"

#include "agile/io/file.hpp"
#include "agile/gpu_timer.hpp"

#include <iostream>
#include <iomanip>
#include <cufft.h>
#include <limits>
#include <cstdlib>
#include <algorithm>

#include "matrixhelper.h"
#include "readSiemensVD11.hpp"
#include "gendicom.hpp"
#include "readdicom.hpp"

/*
#include "gdcmReader.h"
#include "gdcmImageReader.h"
#include "gdcmImageWriter.h"
#include "gdcmWriter.h"
#include "gdcmDataSet.h"
#include "gdcmAttribute.h"
*/




//typedef float TType;
//typedef double TType;
typedef std::complex<float> TType;
//typedef std::complex<double> TType;


//Here starts the main program
int main()
{
  // Initialize the first GPU and print information as done in the first step.
  agile::GPUEnvironment GPU0;
  GPU0.allocateGPU(0);
  GPU0.printInformation(std::cout);

  std::fstream myfile2;
  //init_matrixlog(myfile2,"logmain.txt");

  std::cout << std::endl;

  // Create a vector first on the CPU and transfer it to the GPU.
  std::vector<TType> matrix_read;
  std::vector<TType> matrix_read1;


  std::vector<TType> A_host;
  std::vector<float> float_host;


  unsigned int num_rows = 10;
  unsigned int num_columns = 10;
  unsigned int num_coils = 6;

  // Now we create a dense matrix. The elements are stored in row-major order.
  for (unsigned coil = 0; coil < num_coils; ++coil)
    for (unsigned column = 0; column < num_columns; ++column)
      for (unsigned row = 0; row < num_rows; ++row)
      {
        float_host.push_back(float(row*column+1));

      }

/*
  agile::GPUMatrix<agile::to_real_type<TType>::type>* Coil_save;
  Coil_save = new agile::GPUMatrix<agile::to_real_type<TType>::type> [num_coils];

  agile::GPUMatrix<agile::to_real_type<TType>::type>* Mat_2;
  Mat_2 = new agile::GPUMatrix<agile::to_real_type<TType>::type> [num_coils];
*/


  std::cout<<"\n\n\n test 1: ";

  agile::GPUMatrix<float>* float_test;
  float_test = new agile::GPUMatrix<float> [num_coils];

  agile::GPUMatrix<float> Ergebnis(num_rows,num_columns,NULL);
  agile::GPUMatrix<float> Ergebnis2(num_rows,num_columns,NULL);



  agile::GPUMatrix<TType>* Complex_test;
  Complex_test = new agile::GPUMatrix<TType> [num_coils];


  for(unsigned i=0; i<num_coils;i++)
  {
    float_test[i].assignFromHost(num_rows, num_columns, &float_host[num_rows*num_columns*i]);
    Complex_test[i].assignFromHost(num_rows, num_columns,NULL);
  }

/*
  std::cout<<"\noutput:";
  num_rows=700;
  std::vector<typename agile::to_real_type<TType>::type> data_host;
  std::fstream myfile;
  init_matrixlog(myfile,"matrixlog.txt");

  agile::GPUVector<typename agile::to_real_type<TType>::type> linspace_x_vec(num_rows,0);
  agile::GPUVector<typename agile::to_real_type<TType>::type> linspace_y_vec(num_rows,0);
  agile::GPUMatrix<typename agile::to_real_type<TType>::type> xi_mat(num_rows,num_rows,NULL);
  agile::GPUMatrix<typename agile::to_real_type<TType>::type> eta_mat(num_rows,num_rows,NULL);

  agile::linspace(linspace_x_vec, 700 , 1);
  //agile::linspace(linspace_y_vec, -0.5f , 0.5f);

  linspace_x_vec.copyToHost(data_host);
  std::cout<<"\n linspace_x_vec: ";
  matrixlog(myfile,"linspace_x_vec ",700,1, data_host);
  output("  ",10,10, data_host);

  close_matrixlog(myfile);


  agile::meshgrid(xi_mat,eta_mat,linspace_x_vec,linspace_y_vec);

      xi_mat.copyToHost(data_host);
      std::cout<<"\n Mx: ";
      output("  ",10,10, data_host);

      eta_mat.copyToHost(data_host);
      std::cout<<"\n My: ";
      output("  ",10,10, data_host);
*/



  //"001_000004_000001.dcm"
  Uint8* pData;
  unsigned long length = 0;

  std::ofstream dcmmeas_file("dcm-meas.dat", std::ofstream::binary);
  if (!dcmmeas_file.is_open()) {
    std::cerr << "not able to write to file: " << "dcm-meas.dat" << std::endl;
    return false;
  }

  readdicom("001_000004_000001.dcm",pData, length, dcmmeas_file);
  readdicom("001_000004_000002.dcm",pData, length, dcmmeas_file);
  /*readdicom("001_000005_000001.dcm",pData, length, dcmmeas_file);
  readdicom("001_000005_000002.dcm",pData, length, dcmmeas_file);
  readdicom("001_000005_000003.dcm",pData, length, dcmmeas_file);
  readdicom("001_000005_000004.dcm",pData, length, dcmmeas_file);
  */

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

  int error = read_siemens->readfile(true);
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


  //num_columns = num_coils;

  // ========= IRGN =========

    //return 0;

  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  std::string filepath;
  std::ostringstream filename;
  std::string date(asctime(timeinfo),4,6);

  size_t free;
  size_t total;
  cuMemGetInfo(&free, &total);
  if ((free/1024/1024) <= 150) //break if free memory is lower then 150MB
    return -1;
  std::cout << "\nfree memory: " << free / 1024 / 1024 << "mb, total memory: " << total / 1024 / 1024 << "mb" << std::endl;


  agile::GPUTimer Timer;

  double timer_value;
  Timer.start();


  //agile::readMatrixFile3D("10x10x5imag.bin", num_rows, num_columns, num_coils, matrix_read);
  //agile::readMatrixFile3D("10x10x6imag.bin", num_rows, num_columns, num_coils, matrix_read);
  agile::readMatrixFile3D("braindata.bin", num_rows, num_columns, num_coils, matrix_read);
  //agile::readMatrixFile3D("128x130x6_brain.dat", num_rows, num_columns, num_coils, matrix_read);

  //agile::readMatrixFile3D("newbrain.bin", num_rows, num_columns, num_coils, matrix_read);

  std::cout<<"\n num: "<<num_rows<<"  "<<num_columns<<"  "<<num_coils;

  agile::GPUMatrix<TType>* Coil;
  Coil = new agile::GPUMatrix<TType> [num_coils];

  agile::GPUMatrix<agile::to_real_type<TType>::type>* image;
  std::vector<agile::to_real_type<TType>::type> image_host;
  IRGN_Params irgnpara;

  for(unsigned i=0; i<num_coils;i++)
  {
    Coil[i].assignFromHost(num_rows, num_columns, &matrix_read[num_rows*num_columns*i]);
  } 

  //TVSolve:
  irgnpara.alpha0 = 1;
  irgnpara.alpha_min = std::numeric_limits<float>::min();   //1.17549e-38
  irgnpara.alpha_q = 0.1;
  irgnpara.beta0 = 1;
  irgnpara.beta_min = std::numeric_limits<float>::min();
  irgnpara.beta_q = 0.2;

  irgnpara.tvmax = 1000;
  irgnpara.tvits = 20;
  irgnpara.maxit = 5;

  IRGN<TType>* irgn_tvsolve;
  irgn_tvsolve = new TGVSolve<TType>(Coil,num_coils, irgnpara);

  irgn_tvsolve->HighFreqPenalty();
  irgn_tvsolve->Normalize();
  irgn_tvsolve->Iteration();
  irgn_tvsolve->Postprocess();

  image = irgn_tvsolve->get_image();
  (*image).copyToHost(image_host);

  //agile::writeVectorFile("/home/hheigl/Desktop/IRGN_Pics/img_irgn_TGVsolve_alphaq-0.00025_betaq-0.0005_iter-2.dat",image_host);
  //agile::writeVectorFile("/home/hheigl/Desktop/IRGN_Pics/img_irgn_TGVsolve_alphaq-0.1_betaq-0.2_iter-2.dat",image_host);


  agile::writeVectorFile("bild.dat",image_host);


/*
  //filepath = "/home/hheigl/Desktop/IRGN_Pics/"+date+"/";
  filepath = "./"+date+"/";
  filename << "imgirgn TGV resnorm-"<<std::setprecision(4)<<(irgn_tvsolve->get_nr_k()).back()<<".dat";
  filepath = filepath+filename.str().c_str();

  agile::writeVectorFile(&filepath[0],image_host);
*/

  timer_value = Timer.stop();
  timer_value = timer_value/1000;
  std::cout << "\n\nTime-IRGN-Calc: " << std::setprecision(5)<< timer_value << "[s]  "<< timer_value/60 << "[min]";


  delete irgn_tvsolve;
  irgn_tvsolve = 0;


  gendicom("bild.dcm", image_host, num_rows, num_columns);

  return 0;

}



// End of $Id: gpu_matrix_test.cpp 476 2011-06-16 08:54:14Z freiberger $.
