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

// $Id: gpu_vector_test.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// A fixture that initializes the GPU and sets up two GPUVectors.
struct GPUVectorFixture
{
  GPUVectorFixture()
  {
    agile::GPUEnvironment::allocateGPU(0);

    std::vector<float> x_host, y_host;
    for (unsigned i = 0; i < 3; ++i)
    {
      x_host.push_back(2 * i + 1);
      y_host.push_back(2 * i + 2);
    }
    x.assignFromHost(x_host.begin(), x_host.end());
    y.assignFromHost(y_host.begin(), y_host.end());
  }

  ~GPUVectorFixture()
  {
  }

  agile::GPUVector<float> x, y;
};

BOOST_AUTO_TEST_SUITE( GPUVector_test_suite )

BOOST_AUTO_TEST_CASE( constructor_test )
{
  agile::GPUEnvironment::allocateGPU(0);

  // test for an empty vector
  agile::GPUVector<float> v0;
  BOOST_CHECK_EQUAL( v0.size(), 0 );
  BOOST_CHECK( v0.empty() );

  // test a vector of given size
  agile::GPUVector<float> v10(10);
  BOOST_CHECK_EQUAL( v10.size(), 10);
  BOOST_CHECK( !v10.empty() );
}

BOOST_FIXTURE_TEST_CASE( assignment_test, GPUVectorFixture )
{
  // copy back to the host and check values
  std::vector<float> x_host, y_host;
  x.copyToHost(x_host);
  y.copyToHost(y_host);
  BOOST_CHECK_EQUAL( x_host.size(), 3 );
  BOOST_CHECK_EQUAL( y_host.size(), 3 );
  BOOST_CHECK_EQUAL( x_host[0], 1 );
  BOOST_CHECK_EQUAL( x_host[1], 3 );
  BOOST_CHECK_EQUAL( x_host[2], 5 );
  BOOST_CHECK_EQUAL( y_host[0], 2 );
  BOOST_CHECK_EQUAL( y_host[1], 4 );
  BOOST_CHECK_EQUAL( y_host[2], 6 );
}

BOOST_FIXTURE_TEST_CASE( addition_test, GPUVectorFixture )
{
  agile::GPUVector<float> z(3);
  std::vector<float> z_host;
  addVector(x, y, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host[0],  3 );
  BOOST_CHECK_EQUAL( z_host[1],  7 );
  BOOST_CHECK_EQUAL( z_host[2], 11 );
  // change order of operators
  z_host.clear();
  addVector(y, x, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host[0],  3 );
  BOOST_CHECK_EQUAL( z_host[1],  7 );
  BOOST_CHECK_EQUAL( z_host[2], 11 );
}

BOOST_FIXTURE_TEST_CASE( substraction_test, GPUVectorFixture )
{
  agile::GPUVector<float> z(3);
  std::vector<float> z_host;
  subVector(x, y, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host[0], -1 );
  BOOST_CHECK_EQUAL( z_host[1], -1 );
  BOOST_CHECK_EQUAL( z_host[2], -1 );
  // change order of operators
  z_host.clear();
  subVector(y, x, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host[0], 1 );
  BOOST_CHECK_EQUAL( z_host[1], 1 );
  BOOST_CHECK_EQUAL( z_host[2], 1 );
}

BOOST_FIXTURE_TEST_CASE( scaling_test, GPUVectorFixture )
{
  agile::GPUVector<float> z(3);
  std::vector<float> z_host;
  scale(float(2), x, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host[0],  2 );
  BOOST_CHECK_EQUAL( z_host[1],  6 );
  BOOST_CHECK_EQUAL( z_host[2], 10 );
}

BOOST_FIXTURE_TEST_CASE( scalarproduct_test, GPUVectorFixture )
{
  float sp = getScalarProduct(x, y);
  BOOST_CHECK_EQUAL( sp, 44 );
  // change order of operators
  sp = getScalarProduct(y, x);
  BOOST_CHECK_EQUAL( sp, 44 );
}

BOOST_FIXTURE_TEST_CASE( scale_and_add_test, GPUVectorFixture )
{
  agile::GPUVector<float> z(3);
  std::vector<float> z_host;
  addScaledVector(x, float(2), y, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host[0],  5 );
  BOOST_CHECK_EQUAL( z_host[1], 11 );
  BOOST_CHECK_EQUAL( z_host[2], 17 );
}

BOOST_FIXTURE_TEST_CASE( scale_and_sub_test, GPUVectorFixture )
{
  agile::GPUVector<float> z(3);
  std::vector<float> z_host;
  subScaledVector(x, float(3), y, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host[0],  -5 );
  BOOST_CHECK_EQUAL( z_host[1],  -9 );
  BOOST_CHECK_EQUAL( z_host[2], -13 );
}

BOOST_FIXTURE_TEST_CASE( elementwise_multiplication_test, GPUVectorFixture )
{
  agile::GPUVector<float> z(3);
  std::vector<float> z_host;
  multiplyElementwise(x, y, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host[0],  2 );
  BOOST_CHECK_EQUAL( z_host[1], 12 );
  BOOST_CHECK_EQUAL( z_host[2], 30 );
  // change order or operators
  z_host.clear();
  multiplyElementwise(y, x, z);
  z.copyToHost(z_host);
  BOOST_CHECK_EQUAL( z_host[0],  2 );
  BOOST_CHECK_EQUAL( z_host[1], 12 );
  BOOST_CHECK_EQUAL( z_host[2], 30 );
}

BOOST_AUTO_TEST_SUITE_END()

// End of $Id: gpu_vector_test.cpp 476 2011-06-16 08:54:14Z freiberger $.
