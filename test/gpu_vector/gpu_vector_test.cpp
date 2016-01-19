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

#include <iostream>

#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "../test_defines.h"

int main()
{
  // acquire the device and print some specs
  agile::GPUEnvironment::allocateGPU(0);
  agile::GPUEnvironment::printInformation(std::cout);
  std::cout << std::endl;

  const unsigned SIZE = 5;

  // initialize the vectors
  std::cout << "Initializing vectors" << std::endl;
  PRINT_SECTION("Initializing vectors");
  std::cout << "  x[i] = 2" << std::endl;
  agile::GPUVector<float> x(SIZE, 2.0);
  std::cout << "  y[i] = 5" << std::endl;
  agile::GPUVector<float> y(SIZE, 5.0);
  std::cout << "  u[i] = (3,4)" << std::endl;
  agile::GPUVector<std::complex<float> > u(SIZE, std::complex<float>(3, 4));

  PRINT_SECTION("Testing addition z <= x + y");
  agile::GPUVector<float> z(SIZE);
  // add the two vectors
  agile::addVector(x, y, z);
  // copy the result to the host
  std::vector<float> z_host;
  z.copyToHost(z_host);
  // output the first index
  std::cout << "  z[0] = " << z_host[0] << std::endl;

  PRINT_SECTION("Testing addition v <= u + x");
  agile::GPUVector<std::complex<float> > v(SIZE);
  agile::addVector(u, x, v);
  // copy the result to the host
  std::vector<std::complex<float> > v_host;
  v.copyToHost(v_host);
  // output the first index
  std::cout << "  v[0] = " << v_host[0] << std::endl;

  PRINT_SECTION("Testing addition w <= y + u");
  agile::GPUVector<std::complex<float> > w(SIZE);
  agile::addVector(y, u, w);
  // copy the result to the host
  std::vector<std::complex<float> > w_host;
  w.copyToHost(w_host);
  // output the first index
  std::cout << "  w[0] = " << w_host[0] << std::endl;

  PRINT_SECTION("Testing addition w <= u + v");
  agile::addVector(u, v, w);
  // copy the result to the host
  w.copyToHost(w_host);
  // output the first index
  std::cout << "  w[0] = " << w_host[0] << std::endl;

  PRINT_SECTION("Testing subtraction z <= x - y");
  // subtract the two vectors
  agile::subVector(x, y, z);
  // copy the result to the host
  z.copyToHost(z_host);
  // output the first index
  std::cout << "  z[0] = " << z_host[0] << std::endl;

  PRINT_SECTION("Testing subtraction v <= u - x");
  agile::subVector(u, x, v);
  // copy the result to the host
  v.copyToHost(v_host);
  // output the first index
  std::cout << "  v[0] = " << v_host[0] << std::endl;

  PRINT_SECTION("Testing subtraction w <= y - u");
  agile::subVector(y, u, w);
  // copy the result to the host
  w.copyToHost(w_host);
  // output the first index
  std::cout << "  w[0] = " << w_host[0] << std::endl;

  PRINT_SECTION("Testing subtraction w <= u - v");
  agile::subVector(u, v, w);
  // copy the result to the host
  w.copyToHost(w_host);
  // output the first index
  std::cout << "  w[0] = " << w_host[0] << std::endl;

  PRINT_SECTION("Scalar product");
  std::cout << "(x, y) = " << agile::getScalarProduct(x, y) << std::endl;
  std::cout << "(y, x) = " << agile::getScalarProduct(y, x) << std::endl;
  std::cout << "(u, x) = " << agile::getScalarProduct(u, x) << std::endl;
  std::cout << "(x, u) = " << agile::getScalarProduct(x, u) << std::endl;
  std::cout << "(u, y) = " << agile::getScalarProduct(u, y) << std::endl;
  std::cout << "(y, u) = " << agile::getScalarProduct(y, u) << std::endl;
  std::cout << "(u, v) = " << agile::getScalarProduct(u, v) << std::endl;
  std::cout << "(v, u) = " << agile::getScalarProduct(v, u) << std::endl;

  PRINT_SECTION("2D interpolation");
  //typedef std::complex<float> InterpDataType;
  typedef float InterpDataType;
  std::vector<InterpDataType> iSrcHost(9);
  for (unsigned i=iSrcHost.size(); i--; ) iSrcHost[i] = InterpDataType(i);
  agile::GPUVector<InterpDataType> iSrc;
  iSrc.assignFromHost(iSrcHost.begin(), iSrcHost.end());
  PRINT_VEC("Source vector", iSrcHost);

  std::vector<std::complex<float> > iPosHost;
  iPosHost.resize(9);
  for (unsigned x=3; x--; )
    for (unsigned y=3; y--; )
      iPosHost[y*3+x] = std::complex<float>(x, y);

  agile::GPUVector<std::complex<float> > iPos;
  iPos.assignFromHost(iPosHost.begin(), iPosHost.end());
  PRINT_VEC("Positions", iPosHost);
  agile::GPUVector<InterpDataType> iRes(iPosHost.size());

  agile::interpolate2d(iSrc, 3, 3, false, iPos, iRes);

  // copy the result to the host
  std::vector<InterpDataType> iResHost;
  iRes.copyToHost(iResHost);
  PRINT_VEC("Result (column major)", iResHost);

  agile::interpolate2d(iSrc, 3, 3, true, iPos, iRes);

  // copy the result to the host
  iRes.copyToHost(iResHost);
  PRINT_VEC("Result (row major)", iResHost);

  return 0;
}

// End of $Id: gpu_vector_test.cpp 476 2011-06-16 08:54:14Z freiberger $.
