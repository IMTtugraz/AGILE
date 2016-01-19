
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


#include "agile/calc/l2solve.hpp"
#include "agile/calc/tvsolve.hpp"
#include "agile/calc/tgvsolve.hpp"
//#include "irgn_solve.hpp"

#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_matrix.hpp"
#include "agile/gpu_vector.hpp"

#include "agile/io/file.hpp"
//#include "agile/gpu_complex.hpp"

#include <iostream>
#include <iomanip>
#include <cufft.h>
#include "agile/gpu_timer.hpp"


#include "matrixhelper.h"



//typedef float TType;
//typedef double TType;
typedef std::complex<float> TType;
typedef float TType2;
//typedef std::complex<double> TType;


//Here starts the main program
int main()
{
  // Initialize the first GPU and print information as done in the first step.
  agile::GPUEnvironment GPU0;
  GPU0.allocateGPU(0);
  GPU0.printInformation(std::cout);

  std::fstream myfile2;
  init_matrixlog(myfile2,"logmain.txt");

  std::cout << std::endl;

  // Create a vector first on the CPU and transfer it to the GPU.
  std::vector<std::complex<float> > matrix_read;
  std::vector<std::complex<float> > matrix_read1;


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



  agile::GPUMatrix<TType>* Complex_test;
  Complex_test = new agile::GPUMatrix<TType> [num_coils];


  for(unsigned i=0; i<num_coils;i++)
  {
    float_test[i].assignFromHost(num_rows, num_columns, &float_host[num_rows*num_columns*i]);
    Complex_test[i].assignFromHost(num_rows, num_columns,NULL);
  }





  std::cout<<"\noutput:";





  // ========= IRGN =========


  size_t free;
  size_t total;
  cuMemGetInfo(&free, &total);
  if ((free/1024/1024) <= 50) //break if free memory is lower then 50MB
    return -1;
  std::cout << "\nfree memory: " << free / 1024 / 1024 << "mb, total memory: " << total / 1024 / 1024 << "mb" << std::endl;


  agile::GPUTimer Timer;

  double timer_value;
  Timer.start();


  //agile::readMatrixFile3D("10x10x5imag.bin", num_rows, num_columns, num_coils, matrix_read);
  //agile::readMatrixFile3D("10x10x6imag.bin", num_rows, num_columns, num_coils, matrix_read);
  agile::readMatrixFile3D("braindata.bin", num_rows, num_columns, num_coils, matrix_read);

  agile::GPUMatrix<TType>* Coil;
  Coil = new agile::GPUMatrix<TType> [num_coils];

  agile::GPUMatrix<agile::to_real_type<TType>::type>* image;
  std::vector<agile::to_real_type<TType>::type> image_host;
  agile::IRGN_Params irgnpara;


  for(unsigned i=0; i<num_coils;i++)
  {
    Coil[i].assignFromHost(num_rows, num_columns, &matrix_read[num_rows*num_columns*i]);
  }

  agile::FFT<TType> fftobj_;
  fftobj_.setfftplan(num_rows, num_columns);

  std::vector<TType> data_host;

  Coil[0].copyToHost(data_host);
  matrixlog(myfile2,"orig-data Coil ",1,10, data_host);

/*
  irgn_tvsolve->CenterdIFFTpattern(Coil,Complex_test,6);
  Complex_test[0].copyToHost(data_host);
  matrixlog(myfile2,"CenterdIFFTpattern  ",1,100, data_host);

  irgn_tvsolve->CenterdFFTpattern(Coil,Complex_test,6);
  Complex_test[0].copyToHost(data_host);
  matrixlog(myfile2,"CenterdFFTpattern = Coil ??? ",1,100, data_host);
*/

  cufftResult_t cufftResult;

  std::cout<<"\nsize matrix_read: "<<matrix_read.size();

  //agile::writeMatrixFile3D("braindata.bin",num_rows,num_columns,num_coils,matrix_read);


  float sum=0;
  for(int i=0; i<num_rows*num_columns ;++i)
  {
    std::complex<float> val;
    val = matrix_read[i];

    sum=sum+(val.real()*val.real() + val.imag()*val.imag());
  }
  std::cout<<"\nNorm CPU: "<<std::sqrt(sum);
  float dscale = 0;
  dscale = agile::norm2(Coil[0]);
  std::cout<<"\n Normierung vom Rohdatensatz: "<<dscale;


  //=================================================================================================
    std::cout<<"\n\n\n\n VARIANTE 0\n\n";
    //  ----------no-pattern:___________

        agile::GPUMatrix<TType>* Complex_indata;
        Complex_indata = new agile::GPUMatrix<TType> [num_coils];
        agile::GPUMatrix<TType> Complex_outdata(num_rows,num_columns,NULL);
        agile::readMatrixFile3D("brainraw.bin", num_rows, num_columns, num_coils, matrix_read);
        std::vector<TType> matrix_calc(matrix_read.begin(),matrix_read.end());
        Complex_indata[0].assignFromHost(num_rows, num_columns, &matrix_calc[0]);


        fftobj_.CenterdFFT(Complex_indata[0],Complex_indata[0]);
        Complex_indata[0].copyToHost(data_host);
        matrixlog(myfile2,"CenterdFFT   ",1,10, data_host);

        dscale = 0;
        dscale = agile::norm2(Complex_indata[0]);
        std::cout<<"\n\n Normierung nach CenterdFFT: "<<dscale;

        fftobj_.CenterdIFFT(Complex_indata[0],Complex_outdata);
        Complex_outdata.copyToHost(data_host);
        matrixlog(myfile2,"CenterdIFFT  ",1,10, data_host);

        dscale = 0;
        dscale = agile::norm2(Complex_outdata);
        std::cout<<"\n Normierung nach CenterdIFFT: "<<dscale;

    //  -----------pattern:___________

        agile::readMatrixFile3D("braindata.bin", num_rows, num_columns, num_coils, matrix_read);
        Complex_indata[0].assignFromHost(num_rows, num_columns, &matrix_read[0]);
        fftobj_.calc_pattern(Complex_indata[0]);

        fftobj_.CenterdFFTpattern(Complex_indata[0],Complex_indata[0]);
        Complex_indata[0].copyToHost(data_host);
        matrixlog(myfile2,"CenterdFFTpattern   ",1,10, data_host);

        dscale = 0;
        dscale = agile::norm2(Complex_indata[0]);
        std::cout<<"\n\n Normierung nach CenterdFFTpattern: "<<dscale;

        fftobj_.CenterdIFFTpattern(Complex_indata[0],Complex_outdata);
        Complex_outdata.copyToHost(data_host);
        matrixlog(myfile2,"CenterdIFFTpattern  ",1,10, data_host);

        dscale = 0;
        dscale = agile::norm2(Complex_outdata);
        std::cout<<"\n Normierung nach CenterdIFFTpattern: "<<dscale;


//=================================================================================================
  std::cout<<"\n\n\n\n VARIANTE 1\n\n";

//  agile::GPUMatrix<TType>* Complex_indata;
//  Complex_indata = new agile::GPUMatrix<TType> [num_coils];
//  agile::GPUMatrix<TType> Complex_outdata(num_rows,num_columns,NULL);


  agile::readMatrixFile3D("braindata.bin", num_rows, num_columns, num_coils, matrix_read);
  Complex_indata[0].assignFromHost(num_rows, num_columns, &matrix_read[0]);

  cufftHandle fftplan_;
  cufftPlan2d(&fftplan_, num_rows, num_columns, CUFFT_C2C);

  //agile::ifftshift(Complex_indata[0]);
  cufftResult = cufftExecC2C(fftplan_,
               (agile::substitute_gpu_complex<TType>::cublas_type*)(Complex_indata[0].data()),
               (agile::substitute_gpu_complex<TType>::cublas_type*)(Complex_indata[0].data()),
               CUFFT_FORWARD);
  //agile::fftshift(Complex_indata[0]);
  float val_sqrt;
  //val_sqrt = float(1)/std::sqrt(float(num_rows * num_columns));
  //agile::scale(val_sqrt, Complex_indata[0], Complex_indata[0]);


  Complex_indata[0].copyToHost(data_host);
  matrixlog(myfile2,"CenterdFFT   ",1,10, data_host);

  dscale = 0;
  dscale = agile::norm2(Complex_indata[0]);
  std::cout<<"\n\n Normierung nach FFT: "<<dscale;

  if(cufftResult != 0)
    std::cout<<"\n\nfft-ERROR\n\n";

  //irgn_tvsolve->CenterdIFFT(Complex_test,Complex_test,1);
  //agile::fftshift(Complex_indata[0]);
  cufftResult = cufftExecC2C(fftplan_,
               (cufftComplex*)Complex_indata[0].data(),
               (cufftComplex*)Complex_outdata.data(),
               CUFFT_INVERSE);
  //agile::ifftshift(Complex_outdata);
  val_sqrt = TType::value_type(1)/(Complex_outdata.getNumRows() * Complex_outdata.getNumColumns());
  agile::scale(val_sqrt, Complex_outdata, Complex_outdata);
  Complex_outdata.copyToHost(data_host);

  matrixlog(myfile2,"CenterdIFFT  ",1,10, data_host);

  dscale = 0;
  dscale = agile::norm2(Complex_outdata);
  std::cout<<"\n Normierung nach IFFT: "<<dscale;


  if(cufftResult != 0)
    std::cout<<"\n\nfft-ERROR\n\n";




//=================================================================================================
  std::cout<<"\n\n\n\n VARIANTE 2\n\n";

  int NX1=256;
  int NY1=256;

  unsigned byte_size1 = NX1 * NY1 * sizeof(TType);
  matrix_read1.resize((byte_size1 + sizeof(TType) - 1)
                       / sizeof(TType));

  cufftHandle plan1;

  TType *ttype_idata, *ttype_odata;

  agile::readMatrixFile3D("braindata.bin", num_rows, num_columns, num_coils, matrix_read);

  cudaMalloc((void**)&ttype_idata, sizeof(TType)*NX1*NY1);
  cudaMalloc((void**)&ttype_odata, sizeof(TType)*NX1*NY1);

  CUBLAS_SAFE_CALL(cublasSetMatrix(NX1, NY1, sizeof(TType),&matrix_read[0], NX1, ttype_idata, NX1));

  {
    CUBLAS_SAFE_CALL(cublasGetMatrix (NX1, NY1, sizeof(TType),ttype_idata, NX1, &matrix_read1[0], NX1));

    agile::GPUMatrix<TType> odata_mat(num_rows, num_columns, &matrix_read1[0]);

    float dscale = 0;
    dscale = agile::norm2(odata_mat);
    std::cout<<"\n Vor FFT  2:  2: "<<dscale;

    std::vector<TType> data_host;
    odata_mat.copyToHost(data_host);
    std::cout<<"\n Vor FFT  2: ";
    output("  ",1,10, data_host);

    matrixlog(myfile2," Vor FFT  2: ",1,10, data_host);
  }


  /* Create a 2D FFT plan. */
  cufftPlan2d(&plan1, NX1, NY1, CUFFT_C2C);
  /* Use the CUFFT plan to transform the signal out of place. */
  cufftExecC2C(plan1,
               (agile::substitute_gpu_complex<TType>::cublas_type*) ttype_idata,
               (agile::substitute_gpu_complex<TType>::cublas_type*) ttype_odata,
               CUFFT_FORWARD);


  {
    cublasGetMatrix (NX1, NY1, sizeof(TType),ttype_odata, NX1, &matrix_read1[0], NX1);

    agile::GPUMatrix<TType> odata_mat(num_rows, num_columns, &matrix_read1[0]);

    float dscale = 0;
    dscale = agile::norm2(odata_mat);
    std::cout<<"\n Normierung nach FFT  2: "<<dscale;

    std::vector<TType> data_host;
    odata_mat.copyToHost(data_host);
    std::cout<<"\n Normierung nach FFT  2: ";
    output("  ",1,10, data_host);

    matrixlog(myfile2," Normierung nach FFT  2: ",1,10, data_host);
  }


  /* Note: idata != odata indicates an out‐of‐place transformation
           to CUFFT at execution time. */
  /* Inverse transform the signal in place */
  cufftExecC2C(plan1, (agile::substitute_gpu_complex<TType>::cublas_type*)ttype_odata,
               (agile::substitute_gpu_complex<TType>::cublas_type*)ttype_odata, CUFFT_INVERSE);

  {
    cublasGetMatrix (NX1, NY1, sizeof(TType),ttype_odata, NX1, &matrix_read1[0], NX1);
    agile::GPUMatrix<TType> odata_mat(num_rows, num_columns, &matrix_read1[0]);

    float val_sqrt = TType::value_type(1)/(odata_mat.getNumRows() * odata_mat.getNumColumns());
    agile::scale(val_sqrt, odata_mat, odata_mat);


    float dscale = 0;
    dscale = agile::norm2(odata_mat);
    std::cout<<"\n Normierung nach IFFT  2: "<<dscale;

    std::vector<TType> data_host;
    odata_mat.copyToHost(data_host);
    std::cout<<"\n Normierung nach IFFT  2: ";
    output("  ",1,10, data_host);

    matrixlog(myfile2," Normierung nach IFFT  2: ",1,10, data_host);
  }


  /* Destroy the CUFFT plan. */
  cufftDestroy(plan1);
  cudaFree(ttype_idata); cudaFree(ttype_odata);


//=================================================================================================




//=================================================================================================
  std::cout<<"\n\n\n\n VARIANTE 3\n\n";

  int NX=256;
  int NY=256;

  unsigned byte_size = NX * NY * sizeof(TType);
  matrix_read1.resize((byte_size + sizeof(TType) - 1)
                       / sizeof(TType));

  cufftHandle plan;
  cufftComplex *idata, *odata;

  agile::readMatrixFile3D("braindata.bin", num_rows, num_columns, num_coils, matrix_read);

  cudaMalloc((void**)&idata, sizeof(TType)*NX*NY);
  cudaMalloc((void**)&odata, sizeof(TType)*NX*NY);

  CUBLAS_SAFE_CALL(cublasSetMatrix(NX, NY, sizeof(TType),&matrix_read[0], NX, idata, NX));

  {
    CUBLAS_SAFE_CALL(cublasGetMatrix (NX, NY, sizeof(TType),idata, NX, &matrix_read1[0], NX));

    agile::GPUMatrix<TType> odata_mat(num_rows, num_columns, &matrix_read1[0]);

    float dscale = 0;
    dscale = agile::norm2(odata_mat);
    std::cout<<"\n Vor FFT  2:  2: "<<dscale;

    std::vector<TType> data_host;
    odata_mat.copyToHost(data_host);
    std::cout<<"\n Vor FFT  2: ";
    output("  ",1,10, data_host);

    matrixlog(myfile2," Vor FFT  2: ",1,10, data_host);
  }


  /* Create a 2D FFT plan. */
  cufftPlan2d(&plan, NX, NY, CUFFT_C2C);
  /* Use the CUFFT plan to transform the signal out of place. */
  cufftExecC2C(plan, idata, odata, CUFFT_FORWARD);


  {
    cublasGetMatrix (NX, NY, sizeof(TType),odata, NX, &matrix_read1[0], NX);

    agile::GPUMatrix<TType> odata_mat(num_rows, num_columns, &matrix_read1[0]);

    float dscale = 0;
    dscale = agile::norm2(odata_mat);
    std::cout<<"\n Normierung nach FFT  2: "<<dscale;

    std::vector<TType> data_host;
    odata_mat.copyToHost(data_host);
    std::cout<<"\n Normierung nach FFT  2: ";
    output("  ",1,10, data_host);

    matrixlog(myfile2," Normierung nach FFT  2: ",1,10, data_host);
  }


  /* Note: idata != odata indicates an out‐of‐place transformation
           to CUFFT at execution time. */
  /* Inverse transform the signal in place */
  cufftExecC2C(plan, odata, odata, CUFFT_INVERSE);

  {
    cublasGetMatrix (NX, NY, sizeof(TType),odata, NX, &matrix_read1[0], NX);
    agile::GPUMatrix<TType> odata_mat(num_rows, num_columns, &matrix_read1[0]);

    float val_sqrt = TType::value_type(1)/(odata_mat.getNumRows() * odata_mat.getNumColumns());
    agile::scale(val_sqrt, odata_mat, odata_mat);


    float dscale = 0;
    dscale = agile::norm2(odata_mat);
    std::cout<<"\n Normierung nach IFFT  2: "<<dscale;

    std::vector<TType> data_host;
    odata_mat.copyToHost(data_host);
    std::cout<<"\n Normierung nach IFFT  2: ";
    output("  ",1,10, data_host);

    matrixlog(myfile2," Normierung nach IFFT  2: ",1,10, data_host);
  }


  /* Destroy the CUFFT plan. */
  cufftDestroy(plan);
  cudaFree(idata); cudaFree(odata);


/*
  std::vector<float> matrix_read1;
  int num_r = 4;
  int num_c = 4;

  for (unsigned row = 0; row < num_r; ++row)
    for (unsigned column = 0; column < num_c; ++column)
    {
        matrix_read1.push_back(float(column*row));

    }

  std::cout<<"\nMatrix Trans Mul:";
  output("  ",num_r,num_c, matrix_read1);
  agile::GPUMatrix<float> test_mat(num_r,num_c,&matrix_read1[0]);
  float erg;
  agile::dotProduct(test_mat,test_mat,erg);
  std::cout<<"\nerg: "<<erg;
*/

  //  irgn.HighFreqPenalty();

  //Coil = irgn.get_coil();

/*  for(unsigned i=0; i<irgn.get_numcoils();i++)
  {

      Coil[i].copyToHost(A_host);
      std::cout<<"\nCoil["<<i<<"]";
      output("1  ",Coil[i].getNumRows(),Coil[i].getNumColumns(), A_host);
  }

*/
  /*for(unsigned i=0; i<num_coils;i++)
  {

    fftshift(Coil[i]);
    cufftExecC2C(plan,
                 (cufftComplex*)Coil[i].data(),
                 (cufftComplex*)Coil[i].data(),
                 CUFFT_INVERSE);
    //cufftExecC2C(plan, odata, odata, CUFFT_INVERSE);
    ifftshift(Coil[i]);



    //agile::absMatrix(irgn.get_coil()[i],Coil_save[i]);
    agile::absMatrix(Coil[i],Coil_save[i]);
    agile::multiplyElementwise(Coil_save[i],Coil_save[i],Coil_save[i]);
    std::cout<<"\nnum_row: "<<Coil_save[i].getNumRows()<<"  num_columns: "<<Coil_save[i].getNumColumns();

    if(i>0)
      agile::addMatrix(Coil_save[i-1], Coil_save[i], Coil_save[i]);

  }
  agile::sqrt(Coil_save[num_coils-1],Coil_save[0]);
  cufftDestroy(plan);



  std::vector<agile::to_real_type<TType>::type> Z_host;
  Coil_save[0].copyToHost(Z_host);
  std::cout<<"\nPicture1";
  output("  \n",num_rows, num_columns, Z_host);
  */
  /*
  agile::writeVectorFile("sos_pic.dat",Z_host);
*/

/*
  for(unsigned i=0; i<num_coils;i++)
  {
      Coil_ifft[i].copyToHost(A_host);
      std::cout<<"\nCoil["<<i<<"]";
      output("  ",Coil_ifft[i].getNumRows(),Coil_ifft[i].getNumColumns(), A_host);
  }
*/


//pattern = abs(rawdata(:,:,1))>0;
/*
  agile::GPUMatrix<agile::to_real_type<TType>::type> Z(num_rows, num_columns, NULL);   //Result
  agile::GPUMatrix<TType> Z1(num_rows, num_columns, NULL);   //Result
  agile::GPUMatrix<float> absdata(num_rows, num_columns, NULL);   //Result

  agile::GPUVector<float> Linspace_x(num_columns,NULL);
  agile::GPUVector<float> Linspace_y(num_columns,NULL);
  agile::GPUMatrix<float> Mesh_x(1,1,NULL);
  agile::GPUMatrix<float> Mesh_y(1,1,NULL);

  agile::linspace(Linspace_x, -0.5 , 0.5);
  agile::linspace(Linspace_y, -0.5 , 0.5);
  agile::meshgrid(Mesh_x,Mesh_y,Linspace_x,Linspace_y);



  agile::linspace(Linspace_x, 1 , num_columns);
  agile::linspace(Linspace_y, 1 , num_columns);

  agile::meshgrid(Mesh_x,Mesh_y,Linspace_x,Linspace_y);

  agile::absMatrix(Coil[0],absdata);

  agile::conjMatrix(Coil[2],Z1);


  std::vector<float> absdata_host;
  std::vector<TType> Z_host;
  Coil[0].copyToHost(Z_host);
  output("Coil[0]:  ",absdata.getNumRows(),absdata.getNumColumns(), Z_host);
  absdata.copyToHost(absdata_host);
  output("absdata:  ",absdata.getNumRows(),absdata.getNumColumns(), absdata_host);




  float h=8;
  agile::pow(h,Z,Z);


    if (agile::is_complex<TType>::value)
        std::cout<<"hallo\n";


  std::vector<TType> Z1_host;

  std::vector<float> ALL_host;
*/


  //std::cout << "\nColumns: " << A.getNumColumns() << "    Rows: " << A.getNumRows() << "    Coils: " << num_coils << " \n\n";

  //FFTShift of Matrix A
/*  fftshift(Coil1);
  fftshift(Coil2);
  fftshift(Coil3);
  fftshift(Coil4);
  fftshift(Coil5);
  fftshift(Coil6);
*/

//  agile::multiplyElementwise(C,D,Z);
//  Z.copyToHost(Z_host);

  // Transfer the result back to the host and print it.
  //A.copyToHost(A_host);
  //std::cout<<std::endl;
  //output("A = ",A.getNumRows(),A.getNumColumns(), A_host);
  //B.copyToHost(B_host);
  //output("fftshift(B) = ",B.getNumRows(),B.getNumColumns(), B_host);
  /*output("C*D=Z:  ",Z.getNumRows(),Z.getNumColumns(), Z_host);

  fftshift(Z);
  */

//  Z.copyToHost(Z_host);
//  output("Z:  ",Z.getNumRows(),Z.getNumColumns(), Z_host);
/*  Z1.copyToHost(Z1_host);
  output("Z1:  ",Z1.getNumRows(),Z1.getNumColumns(), Z1_host);

  Linspace_x.copyToHost(ALL_host);
  output("Linspace_x:  ", ALL_host);

  Linspace_y.copyToHost(ALL_host);
  output("Linspace_y:  ", ALL_host);

  Mesh_x.copyToHost(ALL_host);
  output("Mesh_x:  ", Mesh_x.getNumRows(),Mesh_x.getNumColumns(), ALL_host);

  Mesh_y.copyToHost(ALL_host);
  output("Mesh_y:  ", Mesh_y.getNumRows(),Mesh_y.getNumColumns(), ALL_host);
*/

/*
  ALL.copyToHost(ALL_host);
  output("ALL:  ",ALL.getNumRows(),ALL.getNumColumns(), ALL_host);

  fftshift(ALL);
  ALL.copyToHost(ALL_host);
  output("ALL:  ",ALL.getNumRows(),ALL.getNumColumns(), ALL_host);

*/

  close_matrixlog(myfile2);

  return 0;
}



// End of $Id: gpu_matrix_test.cpp 476 2011-06-16 08:54:14Z freiberger $.
