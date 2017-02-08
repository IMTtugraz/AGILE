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

#define TType2EqualTType3 1
#define TType1IsComplex 1
#define TType2IsComplex 0

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

// Unfortunately, textures have to have file scope, which is why we have to
// use this ugly preprocessor definitions.
#define AGILE_TEXTURE agile_matrix_texture_stdcomplexfloatfloat
texture<agile::to_tuple_type<std::complex<float> >::texture_type> AGILE_TEXTURE;

#define AGILE_TEXTURE_2D agile_matrix_texture_2d_stdcomplexfloatfloat
texture<agile::to_tuple_type<std::complex<float> >::texture_type, 2> AGILE_TEXTURE_2D;


#include "gpu_vector.ipp"

namespace agile
{
  
  template 
  void copy<std::complex<float> >(const GPUVector<std::complex<float> >& x, GPUVector<std::complex<float> >& y);
  
  template 
  void maxElement<std::complex<float> >(const GPUVector<std::complex<float> >& x, int* maxVal);

namespace lowlevel
{
  // **************************************************************************
  // functions that depend on one type only
  // **************************************************************************


#if TType2EqualTType3
#if !TType2IsComplex


  template void interpolate2d<std::complex<float>, float >(
    const std::complex<float>* src, unsigned numColumns, unsigned numRows,
    bool reshapeRowMajor, const std::complex<float>* pos,
    std::complex<float>* res, unsigned size);


  template void fftshift<std::complex<float> >(std::complex<float>* x, unsigned size1,
                                    unsigned size2);

  template void ifftshift<std::complex<float> >(std::complex<float>* x, unsigned size1,
                                    unsigned size2);


  template void absVector<std::complex<float> >(
      const std::complex<float>* x,
      typename to_real_type<std::complex<float> >::type* y, unsigned size);

  template void meshgrid<std::complex<float> >(
      std::complex<float>* mesh_x, std::complex<float>* mesh_y,
      const std::complex<float>* x, unsigned x_size, const std::complex<float>* y, unsigned y_size);


  template void imag<std::complex<float> >(
    const std::complex<float>* x,
    typename to_real_type<std::complex<float> >::type* y, unsigned size);

  template typename to_real_type<std::complex<float> >::type norm1(
    const std::complex<float>* x, unsigned size);
  
  template typename to_real_type<std::complex<float> >::type norm2(
    const std::complex<float>* x, unsigned size);

  template void real<std::complex<float> >(
    const std::complex<float>* x,
    typename to_real_type<std::complex<float> >::type* y, unsigned size);

  template void setVectorConstant<std::complex<float> >(
    const std::complex<float>& value, std::complex<float>* x, unsigned size);

  template void pattern<std::complex<float> >(
      const std::complex<float>* x, typename to_real_type<std::complex<float> >::type* z, unsigned size);

  template void diff<std::complex<float> >(
    const unsigned dim, const unsigned x_size, const std::complex<float>* x, std::complex<float>* y, unsigned size);
  
  template void difftrans<std::complex<float> >(
    const unsigned dim, const unsigned x_size, const std::complex<float>* x, std::complex<float>* y, unsigned size);

  template void diff3<std::complex<float> >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const std::complex<float>* x, std::complex<float>* y, unsigned size, bool borderWrap);
  
  template void diff3sym<std::complex<float> >(
      const unsigned dim, const unsigned x_size, const unsigned y_size, const std::complex<float>* x, std::complex<float>* y, unsigned size);
  
  template void diff3trans<std::complex<float> >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const std::complex<float>* x, std::complex<float>* y, unsigned size, bool borderWrap);
  
  template void bdiff3<std::complex<float> >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const std::complex<float>* x, std::complex<float>* y, unsigned size, bool borderWrap);
  
  template void bdiff3sym<std::complex<float> >(
      const unsigned dim, const unsigned x_size, const unsigned y_size, const std::complex<float>* x, std::complex<float>* y, unsigned size);

  template void bdiff3trans<std::complex<float> >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const std::complex<float>* x, std::complex<float>* y, unsigned size, bool borderWrap);

  template void sqrt<std::complex<float> >(const std::complex<float>* x, std::complex<float>* y, unsigned size);

  template void expand_rowdim<std::complex<float> >(const std::complex<float>* x_data, const std::complex<float>* delta_o, const std::complex<float>* delta_u,
                                        unsigned rows, unsigned cols, unsigned row_o, unsigned row_u,
                                        std::complex<float>* z);

  template void expand_coldim<std::complex<float> >(const std::complex<float>* x_data, const std::complex<float>* delta_o, const std::complex<float>* delta_u,
                                        unsigned rows, unsigned cols, unsigned col_o, unsigned col_u,
                                        std::complex<float>* z);

  template void get_content<std::complex<float> >(const std::complex<float>* x_data, unsigned rows, unsigned cols,
                     unsigned row_offset, unsigned col_offset, std::complex<float>* z, unsigned z_rows, unsigned z_cols);

#if !TType1IsComplex

  template void linspace<std::complex<float> >(std::complex<float>* x, unsigned size,
                                    float a, float b);

  template void pow<std::complex<float>, float>(const std::complex<float>& alpha,
                              const float* x,
                              float* y, unsigned size);

#endif  // !TType1IsComplex
#endif  // !TType2IsComplex
#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on two types
  // **************************************************************************
#if TType2EqualTType3

  template void addVector<std::complex<float>, float >(
    const std::complex<float>* x, const float* y,
    typename promote<std::complex<float>, float >::type* z, unsigned size);

  template void divideVector<std::complex<float>, float >(
    const std::complex<float>& alpha, const float* x,
    typename promote<std::complex<float>, float >::type* y, unsigned size);

  template typename promote<std::complex<float>, float >::type getScalarProduct(
    const std::complex<float>* x, const float* y, unsigned size);

  template void multiplyConjElementwise<std::complex<float>, float >(
    const std::complex<float>* x, const float* y,
    typename promote<std::complex<float>, float >::type* z, unsigned size);

  template void multiplyElementwise<std::complex<float>, float >(
    const std::complex<float>* x, const float* y,
    typename promote<std::complex<float>, float >::type* z, unsigned size);

  template void divideElementwise<std::complex<float>, float >(
    const std::complex<float>* x, const float* y,
    typename promote<std::complex<float>, float >::type* z, unsigned size);

  template void scale<std::complex<float>, float >(
    const std::complex<float>& alpha, const float* x,
    typename promote<std::complex<float>, float >::type* y, unsigned size);
    
  template void subVector<std::complex<float>, float >(
    const std::complex<float>* x, const float* y,
    typename promote<std::complex<float>, float >::type* z, unsigned size);

  template void conjVector<std::complex<float> >(
      const std::complex<float>* x,  std::complex<float>* z, unsigned size);
  
  template void expVector<std::complex<float> >(
      const std::complex<float>* x,  std::complex<float>* z, unsigned size);

  template void max<std::complex<float>, float >(
      const std::complex<float>* x1, const float* x2, typename promote<std::complex<float>, float >::type* y, unsigned size);

  template void max<std::complex<float>, float >(
      const std::complex<float>* x1, const float & x2, typename promote<std::complex<float>, float >::type* y, unsigned size);

#if TType1IsComplex
#if !TType2IsComplex
template void phaseVector<std::complex<float>, float >(
    const std::complex<float>* x,
    float* y, unsigned size);
#endif  // TType1IsComplex
#endif  // !TType2IsComplex


#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on three types
  // **************************************************************************

  template void addScaledVector<std::complex<float>, float, float >(
    const std::complex<float>* x, const float& scale, const float* y,
    typename promote<typename promote<std::complex<float>, float >::type,
                     float >::type* z,
    unsigned size);

  template void subScaledVector<std::complex<float>, float, float >(
    const std::complex<float>* x, const float& scale, const float* y,
    typename promote<typename promote<std::complex<float>, float >::type,
                     float >::type* z,
    unsigned size);

} // namespace lowlevel
} // namespace agile

// End of $Id: gpu_vector.cu.in 452 2011-05-31 12:00:18Z freiberger $.
