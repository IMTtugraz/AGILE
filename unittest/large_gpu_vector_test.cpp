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

// $Id: large_gpu_vector_test.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/cpu_vector.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// a size; this is intentially no power of 2 to make sure the operations are
// not dependent on that
const unsigned SIZE = 1024*1024+1;
// tolerance for floating point comparison (in percentage units)
const float TOLERANCE = 1e-3;
const float SCALAR_PRODUCT_TOLERANCE = 1e-3;
// absolute tolerance for subScaledVector
const float SUB_TOLERANCE = 1e-2;

// A fixture that initializes the GPU and sets up two GPUVectors.
struct LargeGPUVectorFixture
{
  LargeGPUVectorFixture()
  {
    agile::GPUEnvironment::allocateGPU(0);

    x_host.resize(SIZE);
    y_host.resize(SIZE);
    for (unsigned counter = 0; counter < SIZE; ++counter)
    {
      x_host[counter] = float(rand()) / RAND_MAX * 10 - 5;
      y_host[counter] = float(rand()) / RAND_MAX * 10 - 5;
    }
    x.assignFromHost(x_host.begin(), x_host.end());
    y.assignFromHost(y_host.begin(), y_host.end());
  }

  ~LargeGPUVectorFixture()
  {
  }

  std::vector<float> x_host, y_host;
  agile::GPUVector<float> x, y;
};

BOOST_FIXTURE_TEST_SUITE( Large_GPUVector_test_suite, LargeGPUVectorFixture )

BOOST_AUTO_TEST_CASE( setup_test )
{
  BOOST_CHECK_EQUAL( x_host.size(), SIZE );
  BOOST_CHECK_EQUAL( y_host.size(), SIZE );
  BOOST_CHECK_EQUAL( x.size(), SIZE );
  BOOST_CHECK_EQUAL( y.size(), SIZE );
}

BOOST_AUTO_TEST_CASE( addition_test )
{
  agile::GPUVector<float> z(SIZE);
  std::vector<float> z_host;
  std::vector<float> z_reference(SIZE);
  addVector(x, y, z);
  agile::addVector(x_host, y_host, z_reference);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), z_reference.size() );
  BOOST_CHECK_EQUAL( z_host.size(), SIZE );

  for (unsigned counter = 0; counter < SIZE; ++counter)
    BOOST_CHECK_CLOSE(z_host[counter], z_reference[counter], TOLERANCE);
}

BOOST_AUTO_TEST_CASE( subtraction_test )
{
  agile::GPUVector<float> z(SIZE);
  std::vector<float> z_host;
  std::vector<float> z_reference(SIZE);
  subVector(x, y, z);
  agile::subVector(x_host, y_host, z_reference);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), z_reference.size() );
  BOOST_CHECK_EQUAL( z_host.size(), SIZE );

  for (unsigned counter = 0; counter < SIZE; ++counter)
    BOOST_CHECK_CLOSE(z_host[counter], z_reference[counter], TOLERANCE);
}

BOOST_AUTO_TEST_CASE( scaling_test )
{
  agile::GPUVector<float> z(SIZE);
  std::vector<float> z_host;
  std::vector<float> z_reference(SIZE);
  scale(float(3), x, z);
  agile::scale(float(3), x_host, z_reference);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), z_reference.size() );
  BOOST_CHECK_EQUAL( z_host.size(), SIZE );

  for (unsigned counter = 0; counter < SIZE; ++counter)
    BOOST_CHECK_CLOSE(z_host[counter], z_reference[counter], TOLERANCE);
}

BOOST_AUTO_TEST_CASE( scalarproduct_test )
{
  float sp = getScalarProduct(x, y);
  float sp_reference = agile::getScalarProduct(x_host, y_host);
  BOOST_CHECK_CLOSE(sp, sp_reference, SCALAR_PRODUCT_TOLERANCE);
}

BOOST_AUTO_TEST_CASE( scale_and_add_test )
{
  agile::GPUVector<float> z(SIZE);
  std::vector<float> z_host;
  std::vector<float> z_reference(SIZE);
  addScaledVector(x, float(4), y, z);
  agile::addScaledVector(x_host, float(4), y_host, z_reference);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), z_reference.size() );
  BOOST_CHECK_EQUAL( z_host.size(), SIZE );

  for (unsigned counter = 0; counter < SIZE; ++counter)
    BOOST_CHECK_CLOSE(z_host[counter], z_reference[counter], TOLERANCE);
}

BOOST_AUTO_TEST_CASE( scale_and_sub_test )
{
  agile::GPUVector<float> z(SIZE);
  std::vector<float> z_host;
  std::vector<float> z_reference(SIZE);
  subScaledVector(x, float(5), y, z);
  agile::subScaledVector(x_host, float(5), y_host, z_reference);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), z_reference.size() );
  BOOST_CHECK_EQUAL( z_host.size(), SIZE );

  for (unsigned counter = 0; counter < SIZE; ++counter)
    BOOST_CHECK(std::abs(z_host[counter]-z_reference[counter]) <SUB_TOLERANCE);
}

BOOST_AUTO_TEST_CASE( elementwise_multiplication_test )
{
  agile::GPUVector<float> z(SIZE);
  std::vector<float> z_host;
  std::vector<float> z_reference(SIZE);
  multiplyElementwise(x, y, z);
  agile::multiplyElementwise(x_host, y_host, z_reference);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host.size(), z_reference.size() );
  BOOST_CHECK_EQUAL( z_host.size(), SIZE );

  for (unsigned counter = 0; counter < SIZE; ++counter)
    BOOST_CHECK_CLOSE(z_host[counter], z_reference[counter], TOLERANCE);
}

BOOST_AUTO_TEST_SUITE_END()

// End of $Id: large_gpu_vector_test.cpp 476 2011-06-16 08:54:14Z freiberger $.
