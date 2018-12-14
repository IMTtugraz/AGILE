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

// $Id: gpu_crs_matrix_test.cpp 476 2011-06-16 08:54:14Z freiberger $

#include <iostream>

#include "agile/gpu_environment.hpp"
#include "agile/gpu_cs_matrix.hpp"
#include "../test_defines.h"

#define MAT_ROWS 6
#define MAT_COLS 5


template <typename TType>
void makeCRS(const TType* matrix, unsigned num_rows, unsigned num_columns,
             std::vector<unsigned>& row_nnz,
             std::vector<unsigned>& column_index,
             std::vector<TType>& data)
{
  const TType* ptr = matrix;
  row_nnz.assign(num_rows, 0);
  column_index.clear();
  data.clear();
  for (unsigned row_counter = 0; row_counter < num_rows; ++row_counter)
    for (unsigned column_counter = 0; column_counter < num_columns;
         ++column_counter)
    {
      if (*ptr != 0)
      {
        ++row_nnz[row_counter];
        column_index.push_back(column_counter);
        data.push_back(*ptr);
      }
      ++ptr;
    }
}

int main()
{
  // acquire the device and print some specs
  agile::GPUEnvironment::allocateGPU(0);
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;

  // create the CRS matrix
  float matrixA[MAT_ROWS][MAT_COLS] = {{1, 0, 0, 4, 7},
                                       {8, 2, 0, 0, 0},
                                       {0, 0, 0, 7, 2},
                                       {0, 0, 0, 0, 0},
                                       {0, 0, 0, 1, 0},
                                       {1, 9, 0, 0, 0}};
                         
  // print out testmatrix (DEBUG)
  //-------------------------------------------------------------------
  std::cout << "full test matrix A:" << std::endl;
  for (unsigned row = 0; row < MAT_ROWS; ++row)
  {
    for (unsigned column = 0; column < MAT_COLS; ++column)
      std::cout << matrixA[row][column] << " ";
    std::cout << std::endl;
  }
  //-------------------------------------------------------------------

  // make CRS matrix
  std::vector<unsigned> row_nnz;
  std::vector<unsigned> column_index;
  std::vector<float> data;
  makeCRS((float*)matrixA, MAT_ROWS, MAT_COLS, row_nnz, column_index, data);
  agile::GPUCSMatrix<float> A(row_nnz, column_index, data);
  

  PRINT_SECTION("real matrix, real vector multiplication");
  //------------------------------------------------------------------
  PRINT_VEC("row_nnz [A]", row_nnz);
  PRINT_VEC("column_index [A]", column_index);
  PRINT_VEC("data [A]", data);
  std::vector<float> x_host(MAT_COLS);
  for (unsigned counter = 0; counter < x_host.size(); ++counter)
    x_host[counter] = (counter + 1) * (counter + 1);
  agile::GPUVector<float> x;
  x.assignFromHost(x_host.begin(), x_host.end());
  PRINT_VEC("x", x_host);
  agile::GPUVector<float> y(MAT_ROWS);
  agile::multiply(A, x, y);
  // copy the result to the host
  std::vector<float> y_host;
  y.copyToHost(y_host);
  PRINT_VEC("y = Ax [real A, real x]", y_host);


  PRINT_SECTION("complex matrix, real vector multiplication");
  //------------------------------------------------------------------
  std::vector<std::complex<float> > complex_data(data.size());
  for (unsigned counter = 0; counter < complex_data.size(); ++counter)
    if (data[counter] != 0)
      complex_data[counter]
        = std::complex<float>(data[counter],
                              ((counter % 2) ? 1 : -1) * data[counter]);
  PRINT_VEC("row_nnz [B]", row_nnz);
  PRINT_VEC("column_index [B]", column_index);
  PRINT_VEC("complex_data [B]", complex_data);
  agile::GPUCSMatrix<std::complex<float> > B(row_nnz, column_index,
                                               complex_data);

  agile::GPUVector<std::complex<float> > z(MAT_ROWS);
  agile::multiply(B, x, z);
  // copy the result to the host
  std::vector<std::complex<float> > z_host;
  z.copyToHost(z_host);
  PRINT_VEC("z = Bx [complex B, real x]", z_host);


  PRINT_SECTION("real matrix, complex vector multiplication");
  //------------------------------------------------------------------
  agile::GPUVector<std::complex<float> > u(MAT_COLS);
  std::vector<std::complex<float> > u_host(u.size());
  for (unsigned counter = 0; counter < u_host.size(); ++counter)
    u_host[counter] = std::complex<float>(
                        float(counter) * 3 - 4, -float(counter) * counter);
  u.assignFromHost(u_host.begin(), u_host.end());
  PRINT_VEC("u", u_host);
  agile::multiply(A, u, z);
  // copy the result to the host
  z.copyToHost(z_host);
  PRINT_VEC("z = Au [real A, complex u]", z_host);
  
  
  PRINT_SECTION("complex matrix, complex vector multiplication");
  //------------------------------------------------------------------
  PRINT_VEC("u", u_host);
  agile::multiply(B, u, z);
  // copy the result to the host
  z.copyToHost(z_host);
  PRINT_VEC("z = Bu [complex B, complex u]", z_host);

  return 0;
}

// End of $Id: gpu_crs_matrix_test.cpp 476 2011-06-16 08:54:14Z freiberger $.

