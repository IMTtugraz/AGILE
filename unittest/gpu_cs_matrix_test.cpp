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

// $Id: gpu_cs_matrix_test.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/gpu_environment.hpp"
#include "agile/gpu_cs_matrix.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// A fixture that creates a CRS matrix and vectors for the multiplication.
struct GPUCRSMatrixFixture
{
  GPUCRSMatrixFixture()
  {
    agile::GPUEnvironment::allocateGPU(0);

    // create the matrix
    const unsigned NUM_ROWS = 6;
    const unsigned NUM_COLUMNS = 5;
    float matrix[NUM_ROWS][NUM_COLUMNS] = {{1, 0, 0, 4, 7},
                                           {8, 2, 0, 0, 0},
                                           {0, 0, 0, 7, 2},
                                           {0, 0, 0, 0, 0},
                                           {0, 0, 0, 1, 0},
                                           {1, 9, 0, 0, 0}};
    std::vector<unsigned> row_nnz(NUM_ROWS, 0);
    std::vector<unsigned> column_index;
    std::vector<float> data;
    for (unsigned row = 0; row < NUM_ROWS; ++row)
      for (unsigned column = 0; column < NUM_COLUMNS; ++column)
        if (matrix[row][column] != 0)
        {
          ++row_nnz[row];
          column_index.push_back(column);
          data.push_back(matrix[row][column]);
        }
    A = new agile::GPUCSMatrix<float>(row_nnz, column_index, data);

    std::vector<float> x_host, y_host;
    for (unsigned i = 0;
         i < (NUM_ROWS > NUM_COLUMNS ? NUM_ROWS : NUM_COLUMNS); ++i)
    {
      if (i < NUM_COLUMNS)
        x_host.push_back(2 * i + 1);
      if (i < NUM_ROWS)
        y_host.push_back(2 * i + 2);
    }
    x.assignFromHost(x_host.begin(), x_host.end());
    y.assignFromHost(y_host.begin(), y_host.end());
  }

  ~GPUCRSMatrixFixture()
  {
    delete A;
  }

  agile::GPUVector<float> x, y;
  agile::GPUCSMatrix<float>* A;
};

// A fixture that creates a CCS matrix and vectors for the multiplication.
struct GPUCCSMatrixFixture
{
  GPUCCSMatrixFixture()
  {
    agile::GPUEnvironment::allocateGPU(0);

    // create the matrix
    const unsigned NUM_ROWS = 6;
    const unsigned NUM_COLUMNS = 5;
    float matrix[NUM_ROWS][NUM_COLUMNS] = {{1, 0, 0, 4, 7},
                                           {8, 2, 0, 0, 0},
                                           {0, 0, 0, 7, 2},
                                           {0, 0, 0, 0, 0},
                                           {0, 0, 0, 1, 0},
                                           {1, 9, 0, 0, 0}};
    std::vector<unsigned> column_nnz(NUM_COLUMNS, 0);
    std::vector<unsigned> row_index;
    std::vector<float> data;
    for (unsigned column = 0; column < NUM_COLUMNS; ++column)
      for (unsigned row = 0; row < NUM_ROWS; ++row)
        if (matrix[row][column] != 0)
        {
          ++column_nnz[column];
          row_index.push_back(row);
          data.push_back(matrix[row][column]);
        }
    A = new agile::GPUCSMatrix<float, true>(column_nnz, row_index, data);

    std::vector<float> x_host, y_host;
    for (unsigned i = 0;
         i < (NUM_ROWS > NUM_COLUMNS ? NUM_ROWS : NUM_COLUMNS); ++i)
    {
      if (i < NUM_COLUMNS)
        x_host.push_back(2 * i + 1);
      if (i < NUM_ROWS)
        y_host.push_back(2 * i + 2);
    }
    x.assignFromHost(x_host.begin(), x_host.end());
    y.assignFromHost(y_host.begin(), y_host.end());
  }

  ~GPUCCSMatrixFixture()
  {
    delete A;
  }

  agile::GPUVector<float> x, y;
  agile::GPUCSMatrix<float, true>* A;
};

BOOST_AUTO_TEST_SUITE( GPUCSMatrix_test_suite )


BOOST_FIXTURE_TEST_CASE( crs_multiplication_test, GPUCRSMatrixFixture )
{
  agile::GPUVector<float> z(A->getNumRows());
  std::vector<float> z_host;
  multiply(*A, x, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), 6 );
  BOOST_CHECK_EQUAL( z_host[0], 92 );
  BOOST_CHECK_EQUAL( z_host[1], 14 );
  BOOST_CHECK_EQUAL( z_host[2], 67 );
  BOOST_CHECK_EQUAL( z_host[3],  0 );
  BOOST_CHECK_EQUAL( z_host[4],  7 );
  BOOST_CHECK_EQUAL( z_host[5], 28 );
}

BOOST_FIXTURE_TEST_CASE( crs_adjoint_multiplication_test, GPUCRSMatrixFixture )
{
  agile::GPUVector<float> z(A->getNumColumns());
  std::vector<float> z_host;
  multiply(y, *A, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), 5 );
  BOOST_CHECK_EQUAL( z_host[0],  46 );
  BOOST_CHECK_EQUAL( z_host[1], 116 );
  BOOST_CHECK_EQUAL( z_host[2],   0 );
  BOOST_CHECK_EQUAL( z_host[3],  60 );
  BOOST_CHECK_EQUAL( z_host[4],  26 );
}

BOOST_FIXTURE_TEST_CASE( ccs_multiplication_test, GPUCCSMatrixFixture )
{
  agile::GPUVector<float> z(A->getNumRows());
  std::vector<float> z_host;
  multiply(*A, x, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), 6 );
  BOOST_CHECK_EQUAL( z_host[0], 92 );
  BOOST_CHECK_EQUAL( z_host[1], 14 );
  BOOST_CHECK_EQUAL( z_host[2], 67 );
  BOOST_CHECK_EQUAL( z_host[3],  0 );
  BOOST_CHECK_EQUAL( z_host[4],  7 );
  BOOST_CHECK_EQUAL( z_host[5], 28 );
}

BOOST_FIXTURE_TEST_CASE( ccs_adjoint_multiplication_test, GPUCCSMatrixFixture )
{
  agile::GPUVector<float> z(A->getNumColumns());
  std::vector<float> z_host;
  multiply(y, *A, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), 5 );
  BOOST_CHECK_EQUAL( z_host[0],  46 );
  BOOST_CHECK_EQUAL( z_host[1], 116 );
  BOOST_CHECK_EQUAL( z_host[2],   0 );
  BOOST_CHECK_EQUAL( z_host[3],  60 );
  BOOST_CHECK_EQUAL( z_host[4],  26 );
}

BOOST_AUTO_TEST_SUITE_END()

// End of $Id: gpu_cs_matrix_test.cpp 476 2011-06-16 08:54:14Z freiberger $.
