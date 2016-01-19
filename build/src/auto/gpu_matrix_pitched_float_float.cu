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

// $Id$

/* This file was generated automatically by CMake. You have to modify '/home2/GIT/AGILE/src/gpu_matrix_pitched.cu.in' if you want to make changes. */

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

#include "gpu_matrix_pitched.ipp"

namespace agile
{
  template void multiply<float, float >(
    const GPUMatrixPitched<float >& A, const GPUVector<float >& x,
    GPUVector<typename promote<float, float >::type>& y);

  template void multiply<float, float >(
    const GPUVector<float >& x, const GPUMatrixPitched<float >& A,
    GPUVector<typename promote<float, float >::type>& y);

  template void multiplyElementwise<float, float >(
    const GPUMatrixPitched<float >& A, const GPUMatrixPitched<float >& B,
    GPUMatrixPitched<typename promote<float, float >::type>& Z);

  template void scale<float, float >(
    const float& alpha, const GPUMatrixPitched<float >& A,
    GPUMatrixPitched<typename promote<float, float >::type>& B);

#if !TType2IsComplex

  // second type may not be complex because, complex is already in signature
  template void interp2d<float, float >(
    const GPUMatrixPitched<float >& M,
    const GPUVector<std::complex<float > >& pos,
    GPUVector<float >& res);

#endif  // TType2IsComplex

} // namespace agile

// End of $Id: gpu_matrix.cu.in 376 2010-02-11 13:16:45Z freiberger $.
