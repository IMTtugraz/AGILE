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

#include "agile/gpu_environment.hpp"
#include "agile/gpu_matrix_pitched.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

struct GPUMatrixFixture
{
  GPUMatrixFixture()
  {
    agile::GPUEnvironment::allocateGPU(0);

    // create the matrix
    float data[4][5] = {{1, 0, 3, 4, 7},
                        {8, 2, 0, 0, 0},
                        {0, 0, 0, 7, 2},
                        {0, 0, 0, 1, 0}};
    A = new agile::GPUMatrixPitched<float>(4, 5, (float*)data);
    std::vector<float> x_host, y_host;
    for (unsigned i = 0; i < 5; ++i)
    {
      x_host.push_back(2 * i + 1);
      if (i < 4)
        y_host.push_back(2 * i + 2);
    }
    x.assignFromHost(x_host.begin(), x_host.end());
    y.assignFromHost(y_host.begin(), y_host.end());
  }

  ~GPUMatrixFixture()
  {
    delete A;
  }

  agile::GPUVector<float> x, y;
  agile::GPUMatrixPitched<float>* A;
};

struct GPUMatrixFixtureComplex
{
  GPUMatrixFixtureComplex()
  {
    agile::GPUEnvironment::allocateGPU(0);

    // create the matrix
    float real_data[4][5] = {{ 1,  0,  3,  4,  7},
                             { 8,  2,  0,  0,  0},
                             { 0,  0,  0,  7,  2},
                             { 0,  0,  0,  1,  0}};
    float imag_data[4][5] = {{ 9, -1,  3,  0,  2},
                             { 7,  0,  0,  2, -5},
                             { 0,  4,  0,  0,  1},
                             {-8,  2,  0,  1,  0}};

    std::vector<std::complex<float> > data;
    for (unsigned row = 0; row < 4; ++row)
      for (unsigned column = 0; column < 5; ++column)
        data.push_back(std::complex<float>(real_data[row][column],
                                           imag_data[row][column]));

    A = new agile::GPUMatrixPitched<std::complex<float> >(4, 5, &data[0]);

    std::vector<std::complex<float> > x_host, y_host;
    for (unsigned i = 0; i < 5; ++i)
    {
      x_host.push_back(std::complex<float>(2 * i + 1, 3 * i + 2));
      if (i < 4)
        y_host.push_back(std::complex<float>(2 * i + 2,
                                             3.0 * i + 2 * ((i % 2) ? 1 : -1)));
    }
    x.assignFromHost(x_host.begin(), x_host.end());
    y.assignFromHost(y_host.begin(), y_host.end());
  }

  ~GPUMatrixFixtureComplex()
  {
    delete A;
  }

  agile::GPUVector<std::complex<float> > x, y;
  agile::GPUMatrixPitched<std::complex<float> >* A;
};


BOOST_AUTO_TEST_SUITE( GPUMatrix_test_suite )

BOOST_AUTO_TEST_CASE( constructor_test )
{
  agile::GPUEnvironment::allocateGPU(0);

  // try to create a 10x10 matrix
  agile::GPUMatrixPitched<float> m10x10(10, 10, 0);
  BOOST_CHECK_EQUAL( m10x10.getNumRows(), 10 );
  BOOST_CHECK_EQUAL( m10x10.getNumColumns(), 10 );
  BOOST_CHECK( m10x10.getPitchElements() >= m10x10.getNumColumns() );
}

BOOST_FIXTURE_TEST_CASE( multiplication_test, GPUMatrixFixture )
{
  agile::GPUVector<float> z(A->getNumRows());
  std::vector<float> z_host;
  multiply(*A, x, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), 4 );
  BOOST_CHECK_EQUAL( z_host[0], 107 );
  BOOST_CHECK_EQUAL( z_host[1],  14 );
  BOOST_CHECK_EQUAL( z_host[2],  67 );
  BOOST_CHECK_EQUAL( z_host[3],   7 );
}

BOOST_FIXTURE_TEST_CASE( adjoint_multiplication_test, GPUMatrixFixture )
{
  agile::GPUVector<float> z(A->getNumColumns());
  std::vector<float> z_host;
  multiply(y, *A, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), 5 );
  BOOST_CHECK_EQUAL( z_host[0], 34 );
  BOOST_CHECK_EQUAL( z_host[1],  8 );
  BOOST_CHECK_EQUAL( z_host[2],  6 );
  BOOST_CHECK_EQUAL( z_host[3], 58 );
  BOOST_CHECK_EQUAL( z_host[4], 26 );
}

BOOST_AUTO_TEST_CASE( complex_constructor_test )
{
  agile::GPUEnvironment::allocateGPU(0);

  // try to create a 10x10 matrix
  agile::GPUMatrixPitched<std::complex<float> > m10x10(10, 10, 0);
  BOOST_CHECK_EQUAL( m10x10.getNumRows(), 10 );
  BOOST_CHECK_EQUAL( m10x10.getNumColumns(), 10 );
  BOOST_CHECK( m10x10.getPitchElements() >= m10x10.getNumColumns() );
}

BOOST_FIXTURE_TEST_CASE( complex_multiplication_test, GPUMatrixFixtureComplex )
{
  agile::GPUVector<std::complex<float> > z(A->getNumRows());
  std::vector<std::complex<float> > z_host;
  multiply(*A, x, z);
  z.copyToHost(z_host);

  BOOST_CHECK_EQUAL( z_host.size(), 4 );
  BOOST_CHECK_EQUAL( real(z_host[0]),  42 );
  BOOST_CHECK_EQUAL( real(z_host[1]),  48 );
  BOOST_CHECK_EQUAL( real(z_host[2]),  33 );
  BOOST_CHECK_EQUAL( real(z_host[3]),   2 );
  BOOST_CHECK_EQUAL( imag(z_host[0]), 207 );
  BOOST_CHECK_EQUAL( imag(z_host[1]),   2 );
  BOOST_CHECK_EQUAL( imag(z_host[2]), 126 );
  BOOST_CHECK_EQUAL( imag(z_host[3]),  16 );
}

BOOST_FIXTURE_TEST_CASE( complex_adjoint_multiplication_test,
                         GPUMatrixFixtureComplex )
{
  agile::GPUVector<std::complex<float> > z(A->getNumColumns());
  std::vector<std::complex<float> > z_host;
  multiply(y, *A, z);
  z.copyToHost(z_host);

  BOOST_CHECK_EQUAL( z_host.size(), 5 );
  BOOST_CHECK_EQUAL( real(z_host[0]), -37 );
  BOOST_CHECK_EQUAL( real(z_host[1]),  48 );
  BOOST_CHECK_EQUAL( real(z_host[2]),   0 );
  BOOST_CHECK_EQUAL( real(z_host[3]),  79 );
  BOOST_CHECK_EQUAL( real(z_host[4]),   1 );
  BOOST_CHECK_EQUAL( imag(z_host[0]),  56 );
  BOOST_CHECK_EQUAL( imag(z_host[1]), -28 );
  BOOST_CHECK_EQUAL( imag(z_host[2]), -12 );
  BOOST_CHECK_EQUAL( imag(z_host[3]),  15 );
  BOOST_CHECK_EQUAL( imag(z_host[4]),   4 );
}

BOOST_AUTO_TEST_SUITE_END()

// End of $Id: gpu_matrix_test.cpp 476 2011-06-16 08:54:14Z freiberger $.
