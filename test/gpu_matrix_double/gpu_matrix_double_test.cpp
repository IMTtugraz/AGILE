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

// $Id: gpu_matrix_double_test.cpp 504 2011-07-27 09:41:14Z freiberger $

// We have to include the headers for the environment, for the GPU vector and
// for the GPU matrix.
#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_matrix.hpp"
#include "agile/gpu_matrix_pitched.hpp"

#include <iostream>
#include <iomanip>



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




//Here starts the main program
int main()
{
  // Initialize the first GPU and print information as done in the first step.
  agile::GPUEnvironment::allocateGPU(0);
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;

  // Create a vector first on the CPU and transfer it to the GPU.
  std::vector<TType> x_host;
  std::vector<TType> matrix_data;

  unsigned int num_rows = 4;
  unsigned int num_columns = 4;

  for (unsigned length = 0;length < num_columns; ++length)
  {
    x_host.push_back(length % 100 + 1);
  }
  agile::GPUVector<TType> x;
  x.assignFromHost(x_host.begin(), x_host.end());

  x.copyToHost(x_host);
  output<TType>("x: ", x_host);


  // We need another vector to store the result of our matrix vector
  // multiplications. NOTE: This matrix has to have the correct dimensions
  // because the library is too lazy to check this!
  agile::GPUVector<TType> y(num_rows);


  // Now we create a dense matrix. The elements are stored in row-major order.
  for (unsigned row = 0; row < num_rows; ++row)
    for (unsigned column = 0; column < num_columns; ++column)
    {
        matrix_data.push_back(TType(column % 30 + 1));
    }

  // We transfer the matrix to the GPU. This can be done using the constructor
  // of \p GPUMatrix, which takes the number of rows, the number of columns
  // and a pointer to an array of size (rows * columns) holding the matrix
  // elements.
  agile::GPUMatrix<TType> A(num_rows, num_columns, &matrix_data[0]);

  // Transfer the result back to the host and print it.
  std::vector<TType> A_host;
  A.copyToHost(A_host);
  output<TType>("A: ", num_rows, num_columns,  A_host);
  std::cout << "Columns: " << A.getNumColumns() << "    Rows: " << A.getNumRows() << " \n\n";

  //Create a Scalar
  //in case of real datatype
  TType scalar=2;
  //in case of complex datatype
  //TType scalar (2,2);
  std::cout << "alpha: " << scalar << std::endl;


  // Now we can use our matrix. Perform the matrix-vector product
  // \f$ y \leftarrow Ax \f$.
  multiply(A, x, y);



  // Transfer the result back to the host and print it.
  std::vector<TType> y_host;
  y.copyToHost(y_host);
  output("A * x: ", y_host);


  // Also the multiplication with the hermitian matrix \f$ A^H = \bar A^T \f$
  // is implemented. It can be evaluated by changing the order of arguments
  // to the \p multiply function (i.e. vector-in, matrix-in, vector-out).
  multiply(x, A, y);

  // output_TType this result, too.
  y.copyToHost(y_host);
  output("A^H * x: ", y_host);



  // Create a Matrix B
  agile::GPUMatrix<TType> B(num_rows, num_columns, NULL);


  // Perform the matrix-scalar product
  scale(scalar, A, B);

  // Transfer the result back to the host and print it.
  std::vector<TType> B_host;
  B.copyToHost(B_host);
  output("B = A * alpha: ",B.getNumRows(), B.getNumColumns(), B_host);



  // Create a Matrix C
  agile::GPUMatrix<TType> C(x.size(), x.size(), NULL);

  //Elementwise Multiplication of 2 matrices
  multiplyElementwise(A, B, C);

  // Transfer the result back to the host and print it.
  std::vector<TType> C_host;
  C.copyToHost(C_host);
  output("C = A(elem) * B(elem): ",C.getNumRows(), C.getNumColumns(), C_host);

//  //FFTShift of Matrix A
//  fftshift(A);

//  // Transfer the result back to the host and print it.
//  A.copyToHost(A_host);
//  output("fftshift(A) = ",A.getNumRows(),A.getNumColumns(), A_host);


  return 0;
}



// End of $Id: gpu_matrix_double_test.cpp 504 2011-07-27 09:41:14Z freiberger $.
