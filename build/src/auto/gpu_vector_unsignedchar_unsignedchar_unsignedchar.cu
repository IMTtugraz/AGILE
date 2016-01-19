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
#define AGILE_TEXTURE agile_matrix_texture_unsignedcharunsignedchar
texture<agile::to_tuple_type<unsigned char >::texture_type> AGILE_TEXTURE;

#define AGILE_TEXTURE_2D agile_matrix_texture_2d_unsignedcharunsignedchar
texture<agile::to_tuple_type<unsigned char >::texture_type, 2> AGILE_TEXTURE_2D;


#include "gpu_vector.ipp"

namespace agile
{
  
  template 
  void copy<unsigned char >(const GPUVector<unsigned char >& x, GPUVector<unsigned char >& y);
  
  template 
  void maxElement<unsigned char >(const GPUVector<unsigned char >& x, int* maxVal);

namespace lowlevel
{
  // **************************************************************************
  // functions that depend on one type only
  // **************************************************************************


#if TType2EqualTType3
#if !TType2IsComplex


  template void interpolate2d<unsigned char, unsigned char >(
    const unsigned char* src, unsigned numColumns, unsigned numRows,
    bool reshapeRowMajor, const std::complex<unsigned char>* pos,
    unsigned char* res, unsigned size);


  template void fftshift<unsigned char >(unsigned char* x, unsigned size1,
                                    unsigned size2);

  template void ifftshift<unsigned char >(unsigned char* x, unsigned size1,
                                    unsigned size2);


  template void absVector<unsigned char >(
      const unsigned char* x,
      typename to_real_type<unsigned char >::type* y, unsigned size);

  template void meshgrid<unsigned char >(
      unsigned char* mesh_x, unsigned char* mesh_y,
      const unsigned char* x, unsigned x_size, const unsigned char* y, unsigned y_size);


  template void imag<unsigned char >(
    const unsigned char* x,
    typename to_real_type<unsigned char >::type* y, unsigned size);

  template typename to_real_type<unsigned char >::type norm1(
    const unsigned char* x, unsigned size);
  
  template typename to_real_type<unsigned char >::type norm2(
    const unsigned char* x, unsigned size);

  template void real<unsigned char >(
    const unsigned char* x,
    typename to_real_type<unsigned char >::type* y, unsigned size);

  template void setVectorConstant<unsigned char >(
    const unsigned char& value, unsigned char* x, unsigned size);

  template void pattern<unsigned char >(
      const unsigned char* x, typename to_real_type<unsigned char >::type* z, unsigned size);

  template void diff<unsigned char >(
    const unsigned dim, const unsigned x_size, const unsigned char* x, unsigned char* y, unsigned size);
  
  template void difftrans<unsigned char >(
    const unsigned dim, const unsigned x_size, const unsigned char* x, unsigned char* y, unsigned size);

  template void diff3<unsigned char >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned char* x, unsigned char* y, unsigned size, bool borderWrap);
  
  template void diff3trans<unsigned char >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned char* x, unsigned char* y, unsigned size, bool borderWrap);
  
  template void bdiff3<unsigned char >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned char* x, unsigned char* y, unsigned size, bool borderWrap);
  
  template void bdiff3trans<unsigned char >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned char* x, unsigned char* y, unsigned size, bool borderWrap);

  template void sqrt<unsigned char >(const unsigned char* x, unsigned char* y, unsigned size);

  template void expand_rowdim<unsigned char >(const unsigned char* x_data, const unsigned char* delta_o, const unsigned char* delta_u,
                                        unsigned rows, unsigned cols, unsigned row_o, unsigned row_u,
                                        unsigned char* z);

  template void expand_coldim<unsigned char >(const unsigned char* x_data, const unsigned char* delta_o, const unsigned char* delta_u,
                                        unsigned rows, unsigned cols, unsigned col_o, unsigned col_u,
                                        unsigned char* z);

  template void get_content<unsigned char >(const unsigned char* x_data, unsigned rows, unsigned cols,
                     unsigned row_offset, unsigned col_offset, unsigned char* z, unsigned z_rows, unsigned z_cols);

#if !TType1IsComplex

  template void linspace<unsigned char >(unsigned char* x, unsigned size,
                                    float a, float b);

  template void pow<unsigned char, unsigned char>(const unsigned char& alpha,
                              const unsigned char* x,
                              unsigned char* y, unsigned size);

#endif  // !TType1IsComplex
#endif  // !TType2IsComplex
#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on two types
  // **************************************************************************
#if TType2EqualTType3

  template void addVector<unsigned char, unsigned char >(
    const unsigned char* x, const unsigned char* y,
    typename promote<unsigned char, unsigned char >::type* z, unsigned size);

  template void divideVector<unsigned char, unsigned char >(
    const unsigned char& alpha, const unsigned char* x,
    typename promote<unsigned char, unsigned char >::type* y, unsigned size);

  template typename promote<unsigned char, unsigned char >::type getScalarProduct(
    const unsigned char* x, const unsigned char* y, unsigned size);

  template void multiplyConjElementwise<unsigned char, unsigned char >(
    const unsigned char* x, const unsigned char* y,
    typename promote<unsigned char, unsigned char >::type* z, unsigned size);

  template void multiplyElementwise<unsigned char, unsigned char >(
    const unsigned char* x, const unsigned char* y,
    typename promote<unsigned char, unsigned char >::type* z, unsigned size);

  template void divideElementwise<unsigned char, unsigned char >(
    const unsigned char* x, const unsigned char* y,
    typename promote<unsigned char, unsigned char >::type* z, unsigned size);

  template void scale<unsigned char, unsigned char >(
    const unsigned char& alpha, const unsigned char* x,
    typename promote<unsigned char, unsigned char >::type* y, unsigned size);
    
  template void subVector<unsigned char, unsigned char >(
    const unsigned char* x, const unsigned char* y,
    typename promote<unsigned char, unsigned char >::type* z, unsigned size);

  template void conjVector<unsigned char >(
      const unsigned char* x,  unsigned char* z, unsigned size);
  
  template void expVector<unsigned char >(
      const unsigned char* x,  unsigned char* z, unsigned size);

  template void max<unsigned char, unsigned char >(
      const unsigned char* x1, const unsigned char* x2, typename promote<unsigned char, unsigned char >::type* y, unsigned size);

  template void max<unsigned char, unsigned char >(
      const unsigned char* x1, const unsigned char & x2, typename promote<unsigned char, unsigned char >::type* y, unsigned size);

#if TType1IsComplex
#if !TType2IsComplex
template void phaseVector<unsigned char, unsigned char >(
    const unsigned char* x,
    unsigned char* y, unsigned size);
#endif  // TType1IsComplex
#endif  // !TType2IsComplex


#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on three types
  // **************************************************************************

  template void addScaledVector<unsigned char, unsigned char, unsigned char >(
    const unsigned char* x, const unsigned char& scale, const unsigned char* y,
    typename promote<typename promote<unsigned char, unsigned char >::type,
                     unsigned char >::type* z,
    unsigned size);

  template void subScaledVector<unsigned char, unsigned char, unsigned char >(
    const unsigned char* x, const unsigned char& scale, const unsigned char* y,
    typename promote<typename promote<unsigned char, unsigned char >::type,
                     unsigned char >::type* z,
    unsigned size);

} // namespace lowlevel
} // namespace agile

// End of $Id: gpu_vector.cu.in 452 2011-05-31 12:00:18Z freiberger $.
