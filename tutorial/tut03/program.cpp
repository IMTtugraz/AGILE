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

// $Id: program.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"

// Instead of including the header for the dense matrices (which was called
// \p gpu_matrix.hpp), we include one for a matrix using compressed row/column
// storage (CRS/CCS). As both are contained in one single class, the class is
// found in the header \p gpu_cs_matrix.hpp
#include "agile/gpu_cs_matrix.hpp"

#include <iostream>
#include <iomanip>
#include <algorithm>  // for \p std::max_element

// \subsection printing Printing functions

// Define a small function that prints a vector to \p std::cout.
void output(const char* string, const std::vector<float>& x)
{
  std::cout << string;
  for (unsigned counter = 0; counter < x.size(); ++counter)
    std::cout << x[counter] << " ";
  std::cout << std::endl;
}

// Another function to print a CRS matrix.
void output(const char* string, const std::vector<unsigned>& row_nnz,
            const std::vector<unsigned>& column_index,
            const std::vector<float>& data)
{
  unsigned num_columns = 1 + *std::max_element(column_index.begin(),
                                               column_index.end());
  // An iterator to the current entry in the current row:
  std::vector<unsigned>::const_iterator row_iter = column_index.begin();
  std::vector<float>::const_iterator data_iter = data.begin();
  std::cout << string << std::endl;
  for (unsigned row = 0; row < row_nnz.size(); ++row)
  {
    // An iterator past the last entry of this row
    std::vector<unsigned>::const_iterator row_end = row_iter + row_nnz[row];

    std::cout << "  ";
    for (unsigned column = 0; column < num_columns; ++column)
      if (column < *row_iter || row_iter == row_end)
        std::cout << std::setw(4) << 0;
      else
      {
        std::cout << std::setw(4) << *data_iter;
        ++data_iter;
        ++row_iter;
      }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

// \subsection main Main program
int main()
{
  // Initialize the first GPU and print information as done in the first step.
  agile::GPUEnvironment::allocateGPU(0);
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;
  
  // Create a vector first on the CPU and transfer it to the GPU.
  std::vector<float> x_host;
  for (unsigned counter = 0; counter < 10; ++counter)
    x_host.push_back(counter * 2 + 1);
  output("x: ", x_host);
  agile::GPUVector<float> x;
  x.assignFromHost(x_host.begin(), x_host.end());

  // \subsubsection compressed Creation from compressed row storage

  // Compressed row storage uses three vectors to describe a matrix.
  // \p row_nnz is a vector containing the number of non-zero (nnz) entries
  // per row (the \p i-th element holds the number of non-zeros in the \p i-th
  // row). This vector also defines the amount of rows the matrix has by its
  // length.
  // \p data is a vector containing the value of the non-zero entries. They
  // are sorted by the row index: First all elements of row 0 are stored,
  // followed by the elements of row 1 and so on. It is not necessary to sort
  // the elements by the column indices. However, this is recommended because
  // mulitplications should be faster then as the number of cache misses is
  // reduced.
  // Finally, \p column_index is a vector containing the column index of each
  // non-zero entry which has the same layout as \p data.

  // We set up these three vectors for a "checkerboard" matrix.
  std::vector<unsigned> row_nnz(x.size(), 0);
  std::vector<unsigned> column_index;
  std::vector<float> data;
  for (unsigned row = 0; row < x.size(); ++row)
    for (unsigned column = 0; column < x.size(); ++column)
      if ((row + column) % 2 == 0)
      {
        ++row_nnz[row];
        column_index.push_back(column);
        data.push_back(float(row + 1) * float(2 * column + 1));
      }

  // Print the matrix to \p std::cout.
  output("A: ", row_nnz, column_index, data);

  // The vectors have to be aligned correctly such that access is efficient
  // for the GPU. This alignment happens automatically when transfering the
  // matrix to the GPU by the constructor of \p GPUCSMatrix. Note: The second
  // template parameter of this class defaults to \p false, which constructs
  // a CRS matrix. If you use \p GPUCSMatrix<float, true>, you will get a
  // compressed column storage (CCS) matrix of type \p float.
  agile::GPUCSMatrix<float> A(row_nnz, column_index, data);

  // A vector for the result of the product.
  agile::GPUVector<float> y(x.size());

  // No we can use the matrix in the same fashion as the dense matrix.
  // For example, computing the product \f$ y \leftarrow Ax \f$ is done with
  // \p multiply.
  multiply(A, x, y);

  std::vector<float> y_host;
  y.copyToHost(y_host);
  output("A * x: ", y_host);

  // Also the multiplication with the hermitian matrix \f$ A^H = \bar A^T \f$
  // is implemented. However, CRS is optimized for traversal along rows, so
  // this operation requires the usage of temporary memory and is usually highly
  // inefficient.
  multiply(x, A, y);

  // Output this result, too.
  y.copyToHost(y_host);
  output("A^H * x: ", y_host);

  // \subsubsection coord Creation from coordinate storage

  // Sometimes it can be practical to generate a CRS matrix from coordinate
  // storage which uses two vectors for the row and column indices and a
  // third one for the data. The same matrix as above in coordinate storage
  // is created now. Note that the order of elements in the coordinate storage
  // will be different due to the switch of the inner and outer loop. This is
  // intentional to test if the CRS matrix will still be generated correctly.
  std::vector<unsigned> coordinate_row;
  std::vector<unsigned> coordinate_column;
  std::vector<float> coordinate_data;
  for (unsigned column = 0; column < x.size(); ++column)
    for (unsigned row = 0; row < x.size(); ++row)
      if ((row + column) % 2 == 0)
      {
        coordinate_row.push_back(row);
        coordinate_column.push_back(column);
        coordinate_data.push_back(float(row + 1) * float(2 * column + 1));
      }

  // This storage format can be converted into a CRS matrix using the
  // method \p createFromCoordinateStorage(). There is no requirement on the
  // layout of the three vectors except that they are of the same size. The
  // function will sort the elements and align them for the GPU automatically.
  agile::GPUCSMatrix<float> A2;
  A2.createFromCoordinateStorage(coordinate_row, coordinate_column,
                                 coordinate_data);

  // We perform the same operations as above to make sure, the coordinate
  // storage works.
  multiply(A2, x, y);
  y.copyToHost(y_host);
  output("A2 * x: ", y_host);

  multiply(x, A2, y);
  y.copyToHost(y_host);
  output("A2^H * x: ", y_host);

  return 0;
}

// End of $Id: program.cpp 476 2011-06-16 08:54:14Z freiberger $.
