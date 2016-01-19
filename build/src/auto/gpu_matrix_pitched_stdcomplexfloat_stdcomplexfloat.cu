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

#define TType1IsComplex 1
#define TType2IsComplex 1

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

// Unfortunately, textures have to have file scope, which is why we have to
// use this ugly preprocessor definitions.
#define AGILE_TEXTURE agile_matrix_texture_stdcomplexfloatstdcomplexfloat
texture<agile::to_tuple_type<std::complex<float> >::texture_type> AGILE_TEXTURE;

#define AGILE_TEXTURE_2D agile_matrix_texture_2d_stdcomplexfloatstdcomplexfloat
texture<agile::to_tuple_type<std::complex<float> >::texture_type, 2> AGILE_TEXTURE_2D;

#include "gpu_matrix_pitched.ipp"

namespace agile
{
  template void multiply<std::complex<float>, std::complex<float> >(
    const GPUMatrixPitched<std::complex<float> >& A, const GPUVector<std::complex<float> >& x,
    GPUVector<typename promote<std::complex<float>, std::complex<float> >::type>& y);

  template void multiply<std::complex<float>, std::complex<float> >(
    const GPUVector<std::complex<float> >& x, const GPUMatrixPitched<std::complex<float> >& A,
    GPUVector<typename promote<std::complex<float>, std::complex<float> >::type>& y);

  template void multiplyElementwise<std::complex<float>, std::complex<float> >(
    const GPUMatrixPitched<std::complex<float> >& A, const GPUMatrixPitched<std::complex<float> >& B,
    GPUMatrixPitched<typename promote<std::complex<float>, std::complex<float> >::type>& Z);

  template void scale<std::complex<float>, std::complex<float> >(
    const std::complex<float>& alpha, const GPUMatrixPitched<std::complex<float> >& A,
    GPUMatrixPitched<typename promote<std::complex<float>, std::complex<float> >::type>& B);

#if !TType2IsComplex

  // second type may not be complex because, complex is already in signature
  template void interp2d<std::complex<float>, std::complex<float> >(
    const GPUMatrixPitched<std::complex<float> >& M,
    const GPUVector<std::complex<std::complex<float> > >& pos,
    GPUVector<std::complex<float> >& res);

#endif  // TType2IsComplex

} // namespace agile

// End of $Id: gpu_matrix.cu.in 376 2010-02-11 13:16:45Z freiberger $.
