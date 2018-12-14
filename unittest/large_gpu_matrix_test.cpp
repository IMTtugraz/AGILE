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

// $Id: large_gpu_matrix_test.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/gpu_environment.hpp"
#include "agile/gpu_matrix_pitched.hpp"
#include "agile/cpu_matrix.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// a size; this is intentially no power of 2 to make sure the operations are
// not dependent on that
const unsigned NUM_ROWS = 5001;
const unsigned NUM_COLUMNS = 2987;
// tolerance for floating point comparison (absolute)
const float TOLERANCE = 1e-2;

// A fixture that initializes the GPU and sets everything.
struct LargeGPUMatrixFixture
{
  LargeGPUMatrixFixture()
  {
    agile::GPUEnvironment::allocateGPU(0);

    x_host.resize(NUM_COLUMNS);
    y_host.resize(NUM_ROWS);
    std::vector<float> matrix_data(NUM_ROWS * NUM_COLUMNS);
    for (unsigned counter = 0; counter < NUM_COLUMNS; ++counter)
      x_host[counter] = float(rand()) / RAND_MAX * 10 - 5;
    for (unsigned counter = 0; counter < NUM_ROWS; ++counter)
      y_host[counter] = float(rand()) / RAND_MAX * 10 - 5;
    for (unsigned row = 0; row < NUM_ROWS; ++row)
      for (unsigned column = 0; column < NUM_COLUMNS; ++column)
        matrix_data[row * NUM_COLUMNS + column]
          = float(rand()) / RAND_MAX * 10 - 5;

    x.assignFromHost(x_host.begin(), x_host.end());
    y.assignFromHost(y_host.begin(), y_host.end());
    A_host = new agile::CPUMatrix<float>(NUM_ROWS, NUM_COLUMNS,
                                           &matrix_data[0]);
    A = new agile::GPUMatrixPitched<float>(NUM_ROWS, NUM_COLUMNS, &matrix_data[0]);
  }

  ~LargeGPUMatrixFixture()
  {
    delete A_host;
    delete A;
  }

  std::vector<float> x_host, y_host;
  agile::GPUVector<float> x, y;
  agile::CPUMatrix<float>* A_host;
  agile::GPUMatrixPitched<float>* A;
};

BOOST_FIXTURE_TEST_SUITE( Large_GPUMatrix_test_suite, LargeGPUMatrixFixture )

BOOST_AUTO_TEST_CASE( setup_test )
{
  BOOST_CHECK_EQUAL( x_host.size(), NUM_COLUMNS );
  BOOST_CHECK_EQUAL( y_host.size(), NUM_ROWS );
  BOOST_CHECK_EQUAL( x.size(), NUM_COLUMNS );
  BOOST_CHECK_EQUAL( y.size(), NUM_ROWS );
  BOOST_CHECK_EQUAL( A_host->getNumRows(), NUM_ROWS );
  BOOST_CHECK_EQUAL( A_host->getNumColumns(), NUM_COLUMNS );
  BOOST_CHECK_EQUAL( A->getNumRows(), NUM_ROWS );
  BOOST_CHECK_EQUAL( A->getNumColumns(), NUM_COLUMNS );
}

BOOST_AUTO_TEST_CASE( multiplication_test )
{
  agile::GPUVector<float> z(NUM_ROWS);
  std::vector<float> z_host;
  std::vector<float> z_reference(NUM_ROWS);
  multiply(*A, x, z);
  agile::multiply(*A_host, x_host, z_reference);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), z_reference.size() );
  BOOST_CHECK_EQUAL( z_host.size(), NUM_ROWS );

  for (unsigned counter = 0; counter < z_host.size(); ++counter)
    BOOST_CHECK(std::abs(z_host[counter] - z_reference[counter]) < TOLERANCE);
}

BOOST_AUTO_TEST_CASE( adjoint_multiplication_test )
{
  agile::GPUVector<float> z(NUM_COLUMNS);
  std::vector<float> z_host;
  std::vector<float> z_reference(NUM_COLUMNS);
  multiply(y, *A, z);
  agile::multiply(y_host, *A_host, z_reference);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), z_reference.size() );
  BOOST_CHECK_EQUAL( z_host.size(), NUM_COLUMNS );

  for (unsigned counter = 0; counter < z_host.size(); ++counter)
    BOOST_CHECK(std::abs(z_host[counter] - z_reference[counter]) < TOLERANCE);
}

BOOST_AUTO_TEST_SUITE_END()
#if 0
// A fixture that initializes the GPU and sets everything.
struct LargeGPUMatrixFixtureComplex
{
  LargeGPUMatrixFixtureComplex()
  {
std::cout << "LargeGPUMatrixFixtureComplex" << std::endl;
    agile::GPUEnvironment::allocateGPU(0);

    x_host.resize(NUM_COLUMNS);
    y_host.resize(NUM_ROWS);
    std::vector<std::complex<float> > matrix_data(NUM_ROWS * NUM_COLUMNS);
    for (unsigned counter = 0; counter < NUM_COLUMNS; ++counter)
      x_host[counter] = std::complex<float>(float(rand()) / RAND_MAX * 10 - 5,
                                            float(rand()) / RAND_MAX * 10 - 5);
    for (unsigned counter = 0; counter < NUM_ROWS; ++counter)
      y_host[counter] = std::complex<float>(float(rand()) / RAND_MAX * 10 - 5,
                                            float(rand()) / RAND_MAX * 10 - 5);
    for (unsigned row = 0; row < NUM_ROWS; ++row)
      for (unsigned column = 0; column < NUM_COLUMNS; ++column)
        matrix_data[row * NUM_COLUMNS + column]
          = std::complex<float>(float(rand()) / RAND_MAX * 10 - 5,
                                float(rand()) / RAND_MAX * 10 - 5);

    x.assignFromHost(x_host.begin(), x_host.end());
    y.assignFromHost(y_host.begin(), y_host.end());
    A_host = new agile::CPUMatrix<std::complex<float> >(NUM_ROWS, NUM_COLUMNS,
                                                          &matrix_data[0]);
    A = new agile::GPUMatrixPitched<std::complex<float> >(NUM_ROWS, NUM_COLUMNS,
                                                     &matrix_data[0]);
std::cout << "LargeGPUMatrixFixtureComplex done" << std::endl;
  }

  ~LargeGPUMatrixFixtureComplex()
  {
    delete A_host;
    delete A;
  }

  std::vector<std::complex<float> > x_host, y_host;
  agile::GPUVector<std::complex<float> > x, y;
  agile::CPUMatrix<std::complex<float> >* A_host;
  agile::GPUMatrixPitched<std::complex<float> >* A;
};

BOOST_FIXTURE_TEST_SUITE( Large_GPUMatrix_test_suite_complex,
                          LargeGPUMatrixFixtureComplex )

BOOST_AUTO_TEST_CASE( setup_test )
{
std::cout << "setup_test started" << std::endl;
  BOOST_CHECK_EQUAL( x_host.size(), NUM_COLUMNS );
  BOOST_CHECK_EQUAL( y_host.size(), NUM_ROWS );
  BOOST_CHECK_EQUAL( x.size(), NUM_COLUMNS );
  BOOST_CHECK_EQUAL( y.size(), NUM_ROWS );
  BOOST_CHECK_EQUAL( A_host->getNumRows(), NUM_ROWS );
  BOOST_CHECK_EQUAL( A_host->getNumColumns(), NUM_COLUMNS );
  BOOST_CHECK_EQUAL( A->getNumRows(), NUM_ROWS );
  BOOST_CHECK_EQUAL( A->getNumColumns(), NUM_COLUMNS );
std::cout << "setup_test done" << std::endl;
}

BOOST_AUTO_TEST_CASE( multiplication_test )
{
std::cout << "multiplication test started" << std::endl;
  agile::GPUVector<std::complex<float> > z(NUM_ROWS);
  std::vector<std::complex<float> > z_host;
  std::vector<std::complex<float> > z_reference(NUM_ROWS);
  multiply(*A, x, z);
  agile::multiply(*A_host, x_host, z_reference);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), z_reference.size() );
  BOOST_CHECK_EQUAL( z_host.size(), NUM_ROWS );

  for (unsigned counter = 0; counter < z_host.size(); ++counter)
    BOOST_CHECK(std::abs(z_host[counter] - z_reference[counter]) < TOLERANCE);
std::cout << "multiplication test done" << std::endl;
}

BOOST_AUTO_TEST_CASE( adjoint_multiplication_test )
{
std::cout << "adjoint multiplication test started" << std::endl;
  agile::GPUVector<std::complex<float> > z(NUM_COLUMNS);
  std::vector<std::complex<float> > z_host;
  std::vector<std::complex<float> > z_reference(NUM_COLUMNS);
  multiply(y, *A, z);
  agile::multiply(y_host, *A_host, z_reference);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), z_reference.size() );
  BOOST_CHECK_EQUAL( z_host.size(), NUM_COLUMNS );

  for (unsigned counter = 0; counter < z_host.size(); ++counter)
    BOOST_CHECK(std::abs(z_host[counter] - z_reference[counter]) < TOLERANCE);
std::cout << "adjoint multiplication test done" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
#endif
// End of $Id: large_gpu_matrix_test.cpp 476 2011-06-16 08:54:14Z freiberger $.
