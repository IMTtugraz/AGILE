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

/* This file was generated automatically by CMake. You have to modify '/home/dieheart/workspace/AGILE/src/gpu_vector.cu.in' if you want to make changes. */

#define TType2EqualTType3 0
#define TType1IsComplex 0
#define TType2IsComplex 0

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

// Unfortunately, textures have to have file scope, which is why we have to
// use this ugly preprocessor definitions.
#define AGILE_TEXTURE agile_matrix_texture_floatfloat
texture<agile::to_tuple_type<float >::texture_type> AGILE_TEXTURE;

#define AGILE_TEXTURE_2D agile_matrix_texture_2d_floatfloat
texture<agile::to_tuple_type<float >::texture_type, 2> AGILE_TEXTURE_2D;


#include "gpu_vector.ipp"

namespace agile
{
  
  template 
  void copy<float >(const GPUVector<float >& x, GPUVector<float >& y);
  
  template 
  void maxElement<float >(const GPUVector<float >& x, int* maxVal);

namespace lowlevel
{
  // **************************************************************************
  // functions that depend on one type only
  // **************************************************************************


#if TType2EqualTType3
#if !TType2IsComplex


  template void interpolate2d<float, float >(
    const float* src, unsigned numColumns, unsigned numRows,
    bool reshapeRowMajor, const std::complex<float>* pos,
    float* res, unsigned size);


  template void fftshift<float >(float* x, unsigned size1,
                                    unsigned size2);

  template void ifftshift<float >(float* x, unsigned size1,
                                    unsigned size2);


  template void absVector<float >(
      const float* x,
      typename to_real_type<float >::type* y, unsigned size);

  template void meshgrid<float >(
      float* mesh_x, float* mesh_y,
      const float* x, unsigned x_size, const float* y, unsigned y_size);


  template void imag<float >(
    const float* x,
    typename to_real_type<float >::type* y, unsigned size);

  template typename to_real_type<float >::type norm1(
    const float* x, unsigned size);
  
  template typename to_real_type<float >::type norm2(
    const float* x, unsigned size);

  template void real<float >(
    const float* x,
    typename to_real_type<float >::type* y, unsigned size);

  template void setVectorConstant<float >(
    const float& value, float* x, unsigned size);

  template void pattern<float >(
      const float* x, typename to_real_type<float >::type* z, unsigned size);

  template void diff<float >(
    const unsigned dim, const unsigned x_size, const float* x, float* y, unsigned size);
  
  template void difftrans<float >(
    const unsigned dim, const unsigned x_size, const float* x, float* y, unsigned size);

  template void diff3<float >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const float* x, float* y, unsigned size, bool borderWrap);
  
  template void diff3sym<float >(
      const unsigned dim, const unsigned x_size, const unsigned y_size, const float* x, float* y, unsigned size);
  
  template void diff3trans<float >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const float* x, float* y, unsigned size, bool borderWrap);
  
  template void bdiff3<float >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const float* x, float* y, unsigned size, bool borderWrap);
  
  template void bdiff3sym<float >(
      const unsigned dim, const unsigned x_size, const unsigned y_size, const float* x, float* y, unsigned size);

  template void bdiff3trans<float >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const float* x, float* y, unsigned size, bool borderWrap);

  template void sqrt<float >(const float* x, float* y, unsigned size);

  template void expand_rowdim<float >(const float* x_data, const float* delta_o, const float* delta_u,
                                        unsigned rows, unsigned cols, unsigned row_o, unsigned row_u,
                                        float* z);

  template void expand_coldim<float >(const float* x_data, const float* delta_o, const float* delta_u,
                                        unsigned rows, unsigned cols, unsigned col_o, unsigned col_u,
                                        float* z);

  template void get_content<float >(const float* x_data, unsigned rows, unsigned cols,
                     unsigned row_offset, unsigned col_offset, float* z, unsigned z_rows, unsigned z_cols);

#if !TType1IsComplex

  template void linspace<float >(float* x, unsigned size,
                                    float a, float b);

  template void pow<float, float>(const float& alpha,
                              const float* x,
                              float* y, unsigned size);

#endif  // !TType1IsComplex
#endif  // !TType2IsComplex
#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on two types
  // **************************************************************************
#if TType2EqualTType3

  template void addVector<float, float >(
    const float* x, const float* y,
    typename promote<float, float >::type* z, unsigned size);

  template void divideVector<float, float >(
    const float& alpha, const float* x,
    typename promote<float, float >::type* y, unsigned size);

  template typename promote<float, float >::type getScalarProduct(
    const float* x, const float* y, unsigned size);

  template void multiplyConjElementwise<float, float >(
    const float* x, const float* y,
    typename promote<float, float >::type* z, unsigned size);

  template void multiplyElementwise<float, float >(
    const float* x, const float* y,
    typename promote<float, float >::type* z, unsigned size);

  template void divideElementwise<float, float >(
    const float* x, const float* y,
    typename promote<float, float >::type* z, unsigned size);

  template void scale<float, float >(
    const float& alpha, const float* x,
    typename promote<float, float >::type* y, unsigned size);
    
  template void subVector<float, float >(
    const float* x, const float* y,
    typename promote<float, float >::type* z, unsigned size);

  template void conjVector<float >(
      const float* x,  float* z, unsigned size);
  
  template void expVector<float >(
      const float* x,  float* z, unsigned size);

  template void max<float, float >(
      const float* x1, const float* x2, typename promote<float, float >::type* y, unsigned size);

  template void max<float, float >(
      const float* x1, const float & x2, typename promote<float, float >::type* y, unsigned size);

#if TType1IsComplex
#if !TType2IsComplex
template void phaseVector<float, float >(
    const float* x,
    float* y, unsigned size);
#endif  // TType1IsComplex
#endif  // !TType2IsComplex


#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on three types
  // **************************************************************************

  template void addScaledVector<float, float, std::complex<float> >(
    const float* x, const float& scale, const std::complex<float>* y,
    typename promote<typename promote<float, float >::type,
                     std::complex<float> >::type* z,
    unsigned size);

  template void subScaledVector<float, float, std::complex<float> >(
    const float* x, const float& scale, const std::complex<float>* y,
    typename promote<typename promote<float, float >::type,
                     std::complex<float> >::type* z,
    unsigned size);

} // namespace lowlevel
} // namespace agile

// End of $Id: gpu_vector.cu.in 452 2011-05-31 12:00:18Z freiberger $.
