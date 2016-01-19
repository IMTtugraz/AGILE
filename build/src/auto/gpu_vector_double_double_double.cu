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

// $Id: gpu_vector.cu.in 452 2011-05-31 12:00:18Z freiberger $

/* This file was generated automatically by CMake. You have to modify '/home2/GIT/AGILE/src/gpu_vector.cu.in' if you want to make changes. */

#define TType2EqualTType3 1
#define TType1IsComplex 0
#define TType2IsComplex 0

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

// Unfortunately, textures have to have file scope, which is why we have to
// use this ugly preprocessor definitions.
#define AGILE_TEXTURE agile_matrix_texture_doubledouble
texture<agile::to_tuple_type<double >::texture_type> AGILE_TEXTURE;

#define AGILE_TEXTURE_2D agile_matrix_texture_2d_doubledouble
texture<agile::to_tuple_type<double >::texture_type, 2> AGILE_TEXTURE_2D;


#include "gpu_vector.ipp"

namespace agile
{
  
  template 
  void copy<double >(const GPUVector<double >& x, GPUVector<double >& y);
  
  template 
  void maxElement<double >(const GPUVector<double >& x, int* maxVal);

namespace lowlevel
{
  // **************************************************************************
  // functions that depend on one type only
  // **************************************************************************


#if TType2EqualTType3
#if !TType2IsComplex


  template void interpolate2d<double, double >(
    const double* src, unsigned numColumns, unsigned numRows,
    bool reshapeRowMajor, const std::complex<double>* pos,
    double* res, unsigned size);


  template void fftshift<double >(double* x, unsigned size1,
                                    unsigned size2);

  template void ifftshift<double >(double* x, unsigned size1,
                                    unsigned size2);


  template void absVector<double >(
      const double* x,
      typename to_real_type<double >::type* y, unsigned size);

  template void meshgrid<double >(
      double* mesh_x, double* mesh_y,
      const double* x, unsigned x_size, const double* y, unsigned y_size);


  template void imag<double >(
    const double* x,
    typename to_real_type<double >::type* y, unsigned size);

  template typename to_real_type<double >::type norm1(
    const double* x, unsigned size);
  
  template typename to_real_type<double >::type norm2(
    const double* x, unsigned size);

  template void real<double >(
    const double* x,
    typename to_real_type<double >::type* y, unsigned size);

  template void setVectorConstant<double >(
    const double& value, double* x, unsigned size);

  template void pattern<double >(
      const double* x, typename to_real_type<double >::type* z, unsigned size);

  template void diff<double >(
    const unsigned dim, const unsigned x_size, const double* x, double* y, unsigned size);
  
  template void difftrans<double >(
    const unsigned dim, const unsigned x_size, const double* x, double* y, unsigned size);

  template void diff3<double >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const double* x, double* y, unsigned size, bool borderWrap);
  
  template void diff3trans<double >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const double* x, double* y, unsigned size, bool borderWrap);
  
  template void bdiff3<double >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const double* x, double* y, unsigned size, bool borderWrap);
  
  template void bdiff3trans<double >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const double* x, double* y, unsigned size, bool borderWrap);

  template void sqrt<double >(const double* x, double* y, unsigned size);

  template void expand_rowdim<double >(const double* x_data, const double* delta_o, const double* delta_u,
                                        unsigned rows, unsigned cols, unsigned row_o, unsigned row_u,
                                        double* z);

  template void expand_coldim<double >(const double* x_data, const double* delta_o, const double* delta_u,
                                        unsigned rows, unsigned cols, unsigned col_o, unsigned col_u,
                                        double* z);

  template void get_content<double >(const double* x_data, unsigned rows, unsigned cols,
                     unsigned row_offset, unsigned col_offset, double* z, unsigned z_rows, unsigned z_cols);

#if !TType1IsComplex

  template void linspace<double >(double* x, unsigned size,
                                    float a, float b);

  template void pow<double, double>(const double& alpha,
                              const double* x,
                              double* y, unsigned size);

#endif  // !TType1IsComplex
#endif  // !TType2IsComplex
#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on two types
  // **************************************************************************
#if TType2EqualTType3

  template void addVector<double, double >(
    const double* x, const double* y,
    typename promote<double, double >::type* z, unsigned size);

  template void divideVector<double, double >(
    const double& alpha, const double* x,
    typename promote<double, double >::type* y, unsigned size);

  template typename promote<double, double >::type getScalarProduct(
    const double* x, const double* y, unsigned size);

  template void multiplyConjElementwise<double, double >(
    const double* x, const double* y,
    typename promote<double, double >::type* z, unsigned size);

  template void multiplyElementwise<double, double >(
    const double* x, const double* y,
    typename promote<double, double >::type* z, unsigned size);

  template void divideElementwise<double, double >(
    const double* x, const double* y,
    typename promote<double, double >::type* z, unsigned size);

  template void scale<double, double >(
    const double& alpha, const double* x,
    typename promote<double, double >::type* y, unsigned size);
    
  template void subVector<double, double >(
    const double* x, const double* y,
    typename promote<double, double >::type* z, unsigned size);

  template void conjVector<double >(
      const double* x,  double* z, unsigned size);
  
  template void expVector<double >(
      const double* x,  double* z, unsigned size);

  template void max<double, double >(
      const double* x1, const double* x2, typename promote<double, double >::type* y, unsigned size);

  template void max<double, double >(
      const double* x1, const double & x2, typename promote<double, double >::type* y, unsigned size);

#if TType1IsComplex
#if !TType2IsComplex
template void phaseVector<double, double >(
    const double* x,
    double* y, unsigned size);
#endif  // TType1IsComplex
#endif  // !TType2IsComplex


#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on three types
  // **************************************************************************

  template void addScaledVector<double, double, double >(
    const double* x, const double& scale, const double* y,
    typename promote<typename promote<double, double >::type,
                     double >::type* z,
    unsigned size);

  template void subScaledVector<double, double, double >(
    const double* x, const double& scale, const double* y,
    typename promote<typename promote<double, double >::type,
                     double >::type* z,
    unsigned size);

} // namespace lowlevel
} // namespace agile

// End of $Id: gpu_vector.cu.in 452 2011-05-31 12:00:18Z freiberger $.
