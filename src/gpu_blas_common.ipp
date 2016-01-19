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

// $Id: gpu_blas_common.ipp 476 2011-06-16 08:54:14Z freiberger $

// ================ common gpu functions for matrix and vector =================

//! \brief Multiply a GPU vector with a scalar (GPU function)
template <typename TType1, typename TType2, typename TType3>
__global__ void multiplyScalar_GPU(
  const TType1 alpha, const TType2* x, TType3* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  if (thread_id < size)
    agile::BasicBinaryScalarOperation<TType1, TType2>::multiply(
      alpha, x[thread_id], y[thread_id]);
}

// End of $Id: gpu_blas_common.ipp 476 2011-06-16 08:54:14Z freiberger $.
