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
#define AGILE_TEXTURE agile_matrix_texture_intint
texture<agile::to_tuple_type<int >::texture_type> AGILE_TEXTURE;

#define AGILE_TEXTURE_2D agile_matrix_texture_2d_intint
texture<agile::to_tuple_type<int >::texture_type, 2> AGILE_TEXTURE_2D;


#include "gpu_vector.ipp"

namespace agile
{
  
  template 
  void copy<int >(const GPUVector<int >& x, GPUVector<int >& y);
  
  template 
  void maxElement<int >(const GPUVector<int >& x, int* maxVal);

namespace lowlevel
{
  // **************************************************************************
  // functions that depend on one type only
  // **************************************************************************


#if TType2EqualTType3
#if !TType2IsComplex


  template void interpolate2d<int, int >(
    const int* src, unsigned numColumns, unsigned numRows,
    bool reshapeRowMajor, const std::complex<int>* pos,
    int* res, unsigned size);


  template void fftshift<int >(int* x, unsigned size1,
                                    unsigned size2);

  template void ifftshift<int >(int* x, unsigned size1,
                                    unsigned size2);


  template void absVector<int >(
      const int* x,
      typename to_real_type<int >::type* y, unsigned size);

  template void meshgrid<int >(
      int* mesh_x, int* mesh_y,
      const int* x, unsigned x_size, const int* y, unsigned y_size);


  template void imag<int >(
    const int* x,
    typename to_real_type<int >::type* y, unsigned size);

  template typename to_real_type<int >::type norm1(
    const int* x, unsigned size);
  
  template typename to_real_type<int >::type norm2(
    const int* x, unsigned size);

  template void real<int >(
    const int* x,
    typename to_real_type<int >::type* y, unsigned size);

  template void setVectorConstant<int >(
    const int& value, int* x, unsigned size);

  template void pattern<int >(
      const int* x, typename to_real_type<int >::type* z, unsigned size);

  template void diff<int >(
    const unsigned dim, const unsigned x_size, const int* x, int* y, unsigned size);
  
  template void difftrans<int >(
    const unsigned dim, const unsigned x_size, const int* x, int* y, unsigned size);

  template void diff3<int >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const int* x, int* y, unsigned size, bool borderWrap);
  
  template void diff3trans<int >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const int* x, int* y, unsigned size, bool borderWrap);
  
  template void bdiff3<int >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const int* x, int* y, unsigned size, bool borderWrap);
  
  template void bdiff3trans<int >(
    const unsigned dim, const unsigned x_size, const unsigned y_size, const int* x, int* y, unsigned size, bool borderWrap);

  template void sqrt<int >(const int* x, int* y, unsigned size);

  template void expand_rowdim<int >(const int* x_data, const int* delta_o, const int* delta_u,
                                        unsigned rows, unsigned cols, unsigned row_o, unsigned row_u,
                                        int* z);

  template void expand_coldim<int >(const int* x_data, const int* delta_o, const int* delta_u,
                                        unsigned rows, unsigned cols, unsigned col_o, unsigned col_u,
                                        int* z);

  template void get_content<int >(const int* x_data, unsigned rows, unsigned cols,
                     unsigned row_offset, unsigned col_offset, int* z, unsigned z_rows, unsigned z_cols);

#if !TType1IsComplex

  template void linspace<int >(int* x, unsigned size,
                                    float a, float b);

  template void pow<int, int>(const int& alpha,
                              const int* x,
                              int* y, unsigned size);

#endif  // !TType1IsComplex
#endif  // !TType2IsComplex
#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on two types
  // **************************************************************************
#if TType2EqualTType3

  template void addVector<int, int >(
    const int* x, const int* y,
    typename promote<int, int >::type* z, unsigned size);

  template void divideVector<int, int >(
    const int& alpha, const int* x,
    typename promote<int, int >::type* y, unsigned size);

  template typename promote<int, int >::type getScalarProduct(
    const int* x, const int* y, unsigned size);

  template void multiplyConjElementwise<int, int >(
    const int* x, const int* y,
    typename promote<int, int >::type* z, unsigned size);

  template void multiplyElementwise<int, int >(
    const int* x, const int* y,
    typename promote<int, int >::type* z, unsigned size);

  template void divideElementwise<int, int >(
    const int* x, const int* y,
    typename promote<int, int >::type* z, unsigned size);

  template void scale<int, int >(
    const int& alpha, const int* x,
    typename promote<int, int >::type* y, unsigned size);
    
  template void subVector<int, int >(
    const int* x, const int* y,
    typename promote<int, int >::type* z, unsigned size);

  template void conjVector<int >(
      const int* x,  int* z, unsigned size);
  
  template void expVector<int >(
      const int* x,  int* z, unsigned size);

  template void max<int, int >(
      const int* x1, const int* x2, typename promote<int, int >::type* y, unsigned size);

  template void max<int, int >(
      const int* x1, const int & x2, typename promote<int, int >::type* y, unsigned size);

#if TType1IsComplex
#if !TType2IsComplex
template void phaseVector<int, int >(
    const int* x,
    int* y, unsigned size);
#endif  // TType1IsComplex
#endif  // !TType2IsComplex


#endif  // TType2EqualTType3


  // **************************************************************************
  // functions that depend on three types
  // **************************************************************************

  template void addScaledVector<int, int, int >(
    const int* x, const int& scale, const int* y,
    typename promote<typename promote<int, int >::type,
                     int >::type* z,
    unsigned size);

  template void subScaledVector<int, int, int >(
    const int* x, const int& scale, const int* y,
    typename promote<typename promote<int, int >::type,
                     int >::type* z,
    unsigned size);

} // namespace lowlevel
} // namespace agile

// End of $Id: gpu_vector.cu.in 452 2011-05-31 12:00:18Z freiberger $.
