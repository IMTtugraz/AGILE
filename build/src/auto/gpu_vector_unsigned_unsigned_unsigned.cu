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
#define AGILE_TEXTURE agile_matrix_texture_unsignedunsigned
texture<agile::to_tuple_type<unsigned >::texture_type> AGILE_TEXTURE;

#define AGILE_TEXTURE_2D agile_matrix_texture_2d_unsignedunsigned
texture<agile::to_tuple_type<unsigned >::texture_type, 2> AGILE_TEXTURE_2D;


#include "gpu_vector.ipp"

namespace agile
{
  
  template 
  void copy<unsigned >(const GPUVector<unsigned >& x, GPUVector<unsigned >& y);
  
  template 
  void maxElement<unsigned >(const GPUVector<unsigned >& x, int* maxVal);

namespace lowlevel
{
  // **************************************************************************
  // functions that depend on one type only
  // **************************************************************************


#if TType2EqualTType3
#if !TType2IsComplex


  template void interpolate2d<unsigned, unsigned >(
    const unsigned* src, unsigned numColumns, unsigned numRows,
    bool reshapeRowMajor, const std::complex<unsigned>* pos,
    unsigned* res, unsigned size);


  template void fftshift<unsigned >(unsigned* x, unsigned size1,
                                    unsigned size2);

  template void ifftshift<unsigned >(unsigned* x, unsigned size1,
                                    unsigned size2);


  template void absVector<unsigned >(
      const unsigned* x,
      typename to_real_type<unsigned >::type* y, unsigned size);

  template void meshgrid<unsigned >(
      unsigned* mesh_x, unsigned* mesh_y,
      const unsigned* x, unsigned x_size, const unsigned* y, unsigned y_size);


  template void imag<unsigned >(
    const unsigned* x,
    typename to_real_type<unsigned >::type* y, unsigned size);

  template typename to_real_type<unsigned >::type norm1(
    const unsigned* x, unsigned size);
  
  template typename to_real_type<unsigned >::type norm2(
    const unsigned* x, unsigned size);

  template void real<unsigned >(
    const unsigned* x,
    typename to_real_type<unsigned >::type* y, unsigned size);

  template void setVectorConstant<unsigned >(
    const unsigned& value, unsigned* x, unsigned size);

  template void pattern<unsigned >(
      const unsigned* x, typename to_real_type<unsigned >::type* z, unsigned size);

  template void diff<unsigned >(
    const unsigned dim, const unsigned x_size, const unsigned* x, unsigned* y, unsigned size);
  
  template void difftrans<unsigned >(
    const unsigned dim, const unsigned x_size, const unsigned* x, unsigned* y, unsigned size);

  template void diff3<unsigned >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned* x, unsigned* y, unsigned size, bool borderWrap);
  
  template void diff3trans<unsigned >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned* x, unsigned* y, unsigned size, bool borderWrap);
  
  template void bdiff3<unsigned >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned* x, unsigned* y, unsigned size, bool borderWrap);
  
  template void bdiff3trans<unsigned >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned* x, unsigned* y, unsigned size, bool borderWrap);

  template void sqrt<unsigned >(const unsigned* x, unsigned* y, unsigned size);

  template void expand_rowdim<unsigned >(const unsigned* x_data, const unsigned* delta_o, const unsigned* delta_u,
                                        unsigned rows, unsigned cols, unsigned row_o, unsigned row_u,
                                        unsigned* z);

  template void expand_coldim<unsigned >(const unsigned* x_data, const unsigned* delta_o, const unsigned* delta_u,
                                        unsigned rows, unsigned cols, unsigned col_o, unsigned col_u,
                                        unsigned* z);

  template void get_content<unsigned >(const unsigned* x_data, unsigned rows, unsigned cols,
                     unsigned row_offset, unsigned col_offset, unsigned* z, unsigned z_rows, unsigned z_cols);

#if !TType1IsComplex

  template void linspace<unsigned >(unsigned* x, unsigned size,
                                    float a, float b);

  template void pow<unsigned, unsigned>(const unsigned& alpha,
                              const unsigned* x,
                              unsigned* y, unsigned size);

#endif  // !TType1IsComplex
#endif  // !TType2IsComplex
#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on two types
  // **************************************************************************
#if TType2EqualTType3

  template void addVector<unsigned, unsigned >(
    const unsigned* x, const unsigned* y,
    typename promote<unsigned, unsigned >::type* z, unsigned size);

  template void divideVector<unsigned, unsigned >(
    const unsigned& alpha, const unsigned* x,
    typename promote<unsigned, unsigned >::type* y, unsigned size);

  template typename promote<unsigned, unsigned >::type getScalarProduct(
    const unsigned* x, const unsigned* y, unsigned size);

  template void multiplyConjElementwise<unsigned, unsigned >(
    const unsigned* x, const unsigned* y,
    typename promote<unsigned, unsigned >::type* z, unsigned size);

  template void multiplyElementwise<unsigned, unsigned >(
    const unsigned* x, const unsigned* y,
    typename promote<unsigned, unsigned >::type* z, unsigned size);

  template void divideElementwise<unsigned, unsigned >(
    const unsigned* x, const unsigned* y,
    typename promote<unsigned, unsigned >::type* z, unsigned size);

  template void scale<unsigned, unsigned >(
    const unsigned& alpha, const unsigned* x,
    typename promote<unsigned, unsigned >::type* y, unsigned size);
    
  template void subVector<unsigned, unsigned >(
    const unsigned* x, const unsigned* y,
    typename promote<unsigned, unsigned >::type* z, unsigned size);

  template void conjVector<unsigned >(
      const unsigned* x,  unsigned* z, unsigned size);
  
  template void expVector<unsigned >(
      const unsigned* x,  unsigned* z, unsigned size);

  template void max<unsigned, unsigned >(
      const unsigned* x1, const unsigned* x2, typename promote<unsigned, unsigned >::type* y, unsigned size);

  template void max<unsigned, unsigned >(
      const unsigned* x1, const unsigned & x2, typename promote<unsigned, unsigned >::type* y, unsigned size);

#if TType1IsComplex
#if !TType2IsComplex
template void phaseVector<unsigned, unsigned >(
    const unsigned* x,
    unsigned* y, unsigned size);
#endif  // TType1IsComplex
#endif  // !TType2IsComplex


#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on three types
  // **************************************************************************

  template void addScaledVector<unsigned, unsigned, unsigned >(
    const unsigned* x, const unsigned& scale, const unsigned* y,
    typename promote<typename promote<unsigned, unsigned >::type,
                     unsigned >::type* z,
    unsigned size);

  template void subScaledVector<unsigned, unsigned, unsigned >(
    const unsigned* x, const unsigned& scale, const unsigned* y,
    typename promote<typename promote<unsigned, unsigned >::type,
                     unsigned >::type* z,
    unsigned size);

} // namespace lowlevel
} // namespace agile

// End of $Id: gpu_vector.cu.in 452 2011-05-31 12:00:18Z freiberger $.
