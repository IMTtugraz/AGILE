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

// $Id: gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp 504 2011-07-27 09:41:14Z freiberger $

// We have to include the headers for the environment, for the GPU vector and
// for the GPU matrix.
#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_matrix.hpp"
#include "agile/gpu_matrix_pitched.hpp"
#include "agile/gpu_timer.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>


void multiplyMatrixVector();
void multiplyMatrixVectorHermitian();
void multiplyMatrixScalar();

typedef double TType;

// Function to print a vector to \p std::cout.
template <typename TType>
void output(const char* string, const std::vector<TType>& x)
{
  std::cout << string;
  for (unsigned counter = 0; counter < x.size(); ++counter)
    std::cout << x[counter] <<std::setprecision(20) << " ";
  std::cout << std::endl;
};

// Function to print a vector to \p std::cout.
template <typename TType>
void output(std::vector<TType> x, std::vector<TType> y)
{
  for (unsigned counter = 0; counter < x.size(); ++counter)
    {
        std::cout << x[counter] << " ";
        std::cout << y[counter] << " ";
        if(x[counter] != y[counter]) {std::cout << " FALSE ";}
        std::cout << std::endl;
    }
  std::cout << std::endl;
};

//Function to print a matrix to \p std::cout.
template <typename TType>
void output(const char* string, unsigned num_rows, unsigned num_columns,
            const std::vector<TType>& data)
{
  typename std::vector<TType>::const_iterator iter = data.begin();
  std::cout << string << std::endl;
  for (unsigned row = 0; row < num_rows; ++row)
  {
    std::cout << "  ";
    for (unsigned column = 0; column < num_columns; ++column)
      std::cout <<std::setprecision(20) << std::setw(4) << *iter++;
    std::cout << std::endl;
  }
  std::cout << std::endl;
}


int main()
{
    // Initialize the GPU and print information as done in the first step.
    agile::GPUEnvironment::allocateGPU(0);
    agile::GPUEnvironment::printInformation(std::cout);
    std::cout << std::endl;

    unsigned int test_nr = 3;

    switch(test_nr)
    {
    case 1:
        multiplyMatrixVector();
    break;

    case 2:
        multiplyMatrixVectorHermitian();
    break;

    case 3:
        multiplyMatrixScalar();
    }

}

//Here starts the main program when testing MultiplyMatrixVector
void multiplyMatrixVector()
{

   agile::GPUTimer Timer;

   double timer_value;

   unsigned int num_rows;
   unsigned int num_columns;

   //Initialisation of the Output Text File
   std::fstream myfile;
   const char * file_name = "../test/gpu_matrix_vs_gpu_matrix_pitched/Timer Results 1404/MultiplyMatrixScalar/R2C_CFloat.xls";
   myfile.open (file_name, std::fstream::out);

   if (!myfile.is_open())
   {
       std::cerr << "File not found: " << file_name << std::endl;
   }

  // Create a vector first on the CPU and transfer it to the GPU.
  std::vector<TType> vector_host;
  std::vector<TType> result_vector_CUBLAS_host;
  std::vector<TType> result_vector_CUDA_host;
  // Create a dense matrix on the CPU. The elements are stored in row-major order.
  std::vector<TType> matrix_host;

  //GPU Vector x
  agile::GPUVector<TType> x;
  //GPU Vector y
  agile::GPUVector<TType> y;
  //GPU Matrix A (CUBLAS)
  agile::GPUMatrix<TType> Matrix_CUBLAS;
  //pitched GPU_Matrix B (CUDA)
  agile::GPUMatrixPitched<TType> Matrix_CUDA;

  //Delete Vector Contents and set size back to 0
  vector_host.clear();
  matrix_host.clear();

  for(unsigned test_ctr = 1; test_ctr <= 1000; test_ctr++)
  {

      num_rows = test_ctr*2;
      num_columns = test_ctr;

      //create Vectors
      for(unsigned length = 0; length < num_columns; length++)
      {

        if((length % 2) == 1)
        {
            vector_host.push_back((TType((length+1) % 30)));
        }
        else
        {
            vector_host.push_back((TType((length+1) % 30))*(TType(-1.0)));
        }
      }


      x.assignFromHost(vector_host.begin(), vector_host.end());
      y.assign(num_rows, 0);

      //create Matrix
      for (unsigned row = 0; row < num_rows; ++row)
        for (unsigned column = 0; column < num_columns; ++column)
          matrix_host.push_back(TType((column + 1) % 30));


      //Transfer the Matrix data from the Host to the GPUMatrix
      Matrix_CUBLAS.assignFromHost(num_rows, num_columns, &matrix_host[0]);
      Matrix_CUDA.assignFromHost(num_rows, num_columns, &matrix_host[0]);


      matrix_host.clear();
      Matrix_CUBLAS.copyToHost(matrix_host);
      //output("Matrix: ", num_rows, num_columns, matrix_host);

      //Multiply Matrix Vector CUDA
      //start Timer
      Timer.start();
      for(unsigned i=0; i<2000; i++)
      {
          // Now we can use our matrix. Perform the matrix-vector product
          multiply(Matrix_CUDA, x, y);
      }
      //stop Timer and write Timer Values to std::cout and file
      timer_value = Timer.stop();
      myfile << std::endl << num_rows << "\t" << num_columns << "\t" << std::setprecision(8)<< timer_value;
      std::cout << std::endl << "NumRows: " << num_rows << "   NumColumns: " << num_columns << "   ------>    CUDA: " << std::setprecision(8)<< timer_value << "[ms]";
      //End Multiply Matrix Vector CUDA

      y.copyToHost(result_vector_CUDA_host);
      y.clear();
      y.assign(num_rows, 0);
      timer_value = 0;

      //Multiply Matrix Vector CUBLAS
      //start Timer
      Timer.start();
      for(unsigned i=0; i<2000; i++)
      {
          // Now we can use our matrix. Perform the matrix-vector product
          // \f$ y \leftarrow Ax \f$.
          multiply(Matrix_CUBLAS, x, y);
      }
      //stop Timer and write Timer Values to std::cout and file
      timer_value = Timer.stop();
      myfile  << "\t" << std::setprecision(8)<< timer_value;
      std::cout << "      CUBLAS: " << std::setprecision(8)<< timer_value << "[ms]";
      //End Multiply Matrix Vector CUBLAS

      y.copyToHost(result_vector_CUBLAS_host);


      x.clear();
      y.clear();
      matrix_host.clear();
      vector_host.clear();

      //Check if the results are equal
      if(result_vector_CUBLAS_host == result_vector_CUDA_host)
      {
          std::cout << "     TRUE";
          //output(result_vector_CUDA_host, result_vector_CUBLAS_host);
          //myfile  << "\t" << "TRUE";
      }
      else
      {
          std::cout << "     FALSE";
          myfile  << "\t" << "FALSE";
          //output(result_vector_CUDA_host, result_vector_CUBLAS_host);
      }

  }


  std::cout << std::endl;
  //Delete Vector Contents and set size back to 0
  vector_host.clear();
  matrix_host.clear();

 Matrix_CUBLAS.~GPUMatrix();
 Matrix_CUDA.~GPUMatrixPitched();

  myfile.close();

}



//Here starts the main program when testing MultiplyMatrixVectorHermitian
void multiplyMatrixVectorHermitian()
{

    agile::GPUTimer Timer;

    double timer_value;

    unsigned int num_rows;
    unsigned int num_columns;

    //Initialisation of the Output Text File
    std::fstream myfile;
    const char * file_name = "../test/gpu_matrix_vs_gpu_matrix_pitched/Timer Results 1404/MultiplyMatrixScalar/R2C_CFloat.xls";
    myfile.open (file_name, std::fstream::out);

    if (!myfile.is_open())
    {
        std::cerr << "File not found: " << file_name << std::endl;
    }

   // Create a vector first on the CPU and transfer it to the GPU.
   std::vector<TType> vector_host;
   std::vector<TType> result_vector_CUBLAS_host;
   std::vector<TType> result_vector_CUDA_host;
   // Create a dense matrix on the CPU. The elements are stored in row-major order.
   std::vector<TType> matrix_host;

   //GPU Vector x
   agile::GPUVector<TType> x;
   //GPU Vector y
   agile::GPUVector<TType> y;
   //GPU Matrix A (CUBLAS)
   agile::GPUMatrix<TType> Matrix_CUBLAS;
   //pitched GPU_Matrix B (CUDA)
   agile::GPUMatrixPitched<TType> Matrix_CUDA;

   //Delete Vector Contents and set size back to 0
   vector_host.clear();
   matrix_host.clear();

   for(unsigned test_ctr = 1; test_ctr <= 1000; test_ctr++)
   {

       num_rows = test_ctr*2;
       num_columns = test_ctr;

       //create Vectors
       for(unsigned length = 0; length < num_rows; length++)
       {

         if((length % 2) == 1)
         {
             vector_host.push_back((TType((length+1) % 30)));
         }
         else
         {
             vector_host.push_back((TType((length+1) % 30))*(TType(-1.0)));
         }
       }


       x.assignFromHost(vector_host.begin(), vector_host.end());
       y.assign(num_rows, 0);

       //create Matrix
       for (unsigned row = 0; row < num_rows; ++row)
         for (unsigned column = 0; column < num_columns; ++column)
           matrix_host.push_back(TType((column + 1) % 30));


       //Transfer the Matrix data from the Host to the GPUMatrix
       Matrix_CUBLAS.assignFromHost(num_rows, num_columns, &matrix_host[0]);
       Matrix_CUDA.assignFromHost(num_rows, num_columns, &matrix_host[0]);


       matrix_host.clear();
       Matrix_CUBLAS.copyToHost(matrix_host);
       //output("Matrix: ", num_rows, num_columns, matrix_host);

       //Multiply Matrix Vector Hermitian CUDA
       //start Timer
       Timer.start();
       for(unsigned i=0; i<2000; i++)
       {
           // Now we can use our matrix. Perform the matrix-vector product
           multiply( x, Matrix_CUDA, y);
       }
       //stop Timer and write Timer Values to std::cout and file
       timer_value = Timer.stop();
       myfile << std::endl << num_rows << "\t" << num_columns << "\t" << std::setprecision(8)<< timer_value;
       std::cout << std::endl << "NumRows: " << num_rows << "   NumColumns: " << num_columns << "   ------>    CUDA: " << std::setprecision(8)<< timer_value << "[ms]";
       //End Multiply Matrix Vector CUDA

       y.copyToHost(result_vector_CUDA_host);
       y.clear();
       y.assign(num_rows, 0);
       timer_value = 0;

       //Multiply Matrix Vector Hermitian CUBLAS
       //start Timer
       Timer.start();
       for(unsigned i=0; i<2000; i++)
       {
           // Now we can use our matrix. Perform the matrix-vector product
           // \f$ y \leftarrow Ax \f$.
           multiply( x, Matrix_CUBLAS, y);
       }
       //stop Timer and write Timer Values to std::cout and file
       timer_value = Timer.stop();
       myfile  << "\t" << std::setprecision(8)<< timer_value;
       std::cout << "      CUBLAS: " << std::setprecision(8)<< timer_value << "[ms]";
       //End Multiply Matrix Vector CUBLAS

       y.copyToHost(result_vector_CUBLAS_host);


       x.clear();
       y.clear();
       matrix_host.clear();
       vector_host.clear();

       //Check if the results are equal
       if(result_vector_CUBLAS_host == result_vector_CUDA_host)
       {
           std::cout << "     TRUE";
           //output(result_vector_CUDA_host, result_vector_CUBLAS_host);
           //myfile  << "\t" << "TRUE";
       }
       else
       {
           std::cout << "     FALSE";
           myfile  << "\t" << "FALSE";
           //output(result_vector_CUDA_host, result_vector_CUBLAS_host);
       }

   }


   std::cout << std::endl;
   //Delete Vector Contents and set size back to 0
   vector_host.clear();
   matrix_host.clear();

  Matrix_CUBLAS.~GPUMatrix();
  Matrix_CUDA.~GPUMatrixPitched();

   myfile.close();
}




//Here starts the main program when testing Multiply Matrix_Scalar
void multiplyMatrixScalar()
{

   agile::GPUTimer Timer;

   double timer_value;


   //Initialisation of the Output Text File
   std::fstream myfile;
   const char * file_name = "../test/gpu_matrix_vs_gpu_matrix_pitched/Timer Results 1404/MultiplyMatrixScalar/R2C_CFloat.xls";
   myfile.open (file_name, std::fstream::out);

   if (!myfile.is_open())
   {
       std::cerr << "File not found: " << file_name << std::endl;
   }

  // Create a vector first on the CPU and transfer it to the GPU.
  std::vector<TType> result_vector_CUBLAS_host;
  std::vector<TType> result_vector_CUDA_host;
  // Create a dense matrix on the CPU. The elements are stored in row-major order.
  std::vector<TType> matrix_host;


  //GPU Matrix A (CUBLAS)
  agile::GPUMatrix<TType> Matrix_CUBLAS;
  //pitched GPU_Matrix B (CUDA)
  agile::GPUMatrixPitched<TType> Matrix_CUDA;

  //GPU Matrix for Result  (CUBLAS)
  agile::GPUMatrix<TType> ResultMatrix_CUBLAS;
  //pitched GPU_Matrix for Result (CUDA)
  agile::GPUMatrixPitched<TType> ResultMatrix_CUDA;

  //Delete Vector Contents and set size back to 0
  matrix_host.clear();
  result_vector_CUDA_host.clear();
  result_vector_CUDA_host.clear();

  TType alpha = 10.123456789;

  unsigned int num_rows;
  unsigned int num_columns;

  //myfile << "\n\tMultiply Matrix Scalar CUBLAS";
  for(unsigned test_ctr = 1; test_ctr <= 1000; test_ctr++)
  {

      num_rows = test_ctr*2;
      num_columns = test_ctr;

      //create Matrix
      for (unsigned row = 0; row < num_rows; ++row)
        for (unsigned column = 0; column < num_columns; ++column)
          {
            if((column % 2) == 1)
            {
                matrix_host.push_back((row % 30 + 1)*(column % 30 + 1));
            }
            else
            {
               matrix_host.push_back((row % 30 + 1)*(column % 30 + 1)*(-1));
            }


          }

      //Transfer the Matrix data from the Host to the GPUMatrix
      Matrix_CUBLAS.assignFromHost(num_rows, num_columns, &matrix_host[0]);
      Matrix_CUDA.assignFromHost(num_rows, num_columns, &matrix_host[0]);

      ResultMatrix_CUBLAS.assignFromHost(num_rows, num_columns, NULL);
      ResultMatrix_CUDA.assignFromHost(num_rows, num_columns, NULL);

      matrix_host.clear();
      result_vector_CUDA_host.clear();
      result_vector_CUDA_host.clear();

      //Multiply Matrix Scalar CUDA
      //start Timer
      Timer.start();
      for(unsigned i=0; i<2000; i++)
      {
          // Now we can use our matrix. Perform the matrix-scalar product
          // \f$ y \leftarrow Ax \f$.
          scale(alpha, Matrix_CUDA, ResultMatrix_CUDA);
      }
      //stop Timer and write Timer Values to std::cout and file
      timer_value = Timer.stop();
      myfile << std::endl << num_rows << "\t" << num_columns << "\t" << std::setprecision(8)<< timer_value;
      std::cout << std::endl << "NumRows: " << num_rows << "   NumColumns: " << num_columns << "   ------>    CUDA: " << std::setprecision(8)<< timer_value << "[ms]";
      //End Multiply Matrix Vector CUDA

      ResultMatrix_CUDA.copyToHost(result_vector_CUDA_host);
      timer_value = 0;

      //Multiply Matrix Scalar CUBLAS
      //start Timer
      Timer.start();
      for(unsigned i=0; i<2000; i++)
      {
          // Now we can use our matrix. Perform the matrix-scalar product
          // \f$ y \leftarrow Ax \f$.
          scale(alpha, Matrix_CUBLAS, ResultMatrix_CUBLAS);
      }
      //stop Timer and write Timer Values to std::cout and file
      timer_value = Timer.stop();
      myfile  << "\t" << std::setprecision(8)<< timer_value;
      std::cout << "      CUBLAS: " << std::setprecision(8)<< timer_value << "[ms]";
      //End Multiply Matrix Vector CUBLAS

      ResultMatrix_CUBLAS.copyToHost(result_vector_CUBLAS_host);



      //Check if the results are equal
      if(result_vector_CUBLAS_host == result_vector_CUDA_host)
      {
          std::cout << "     TRUE";
          //output(result_vector_CUDA_host, result_vector_CUBLAS_host);
          //myfile  << "\t" << "TRUE";
      }
      else
      {
          std::cout << "     FALSE";
          myfile  << "\t" << "FALSE";
          //output(result_vector_CUDA_host, result_vector_CUBLAS_host);
      }

  }

  Matrix_CUDA.~GPUMatrixPitched();
  Matrix_CUBLAS.~GPUMatrix();
  ResultMatrix_CUDA.~GPUMatrixPitched();
  ResultMatrix_CUBLAS.~GPUMatrix();

  std::cout << std::endl;
  //Delete Vector Contents and set size back to 0
  matrix_host.clear();

  myfile.close();

}


// End of $Id: gpu_matrix_vs_gpu_matrix_pitched_double_test.cpp 504 2011-07-27 09:41:14Z freiberger $.
