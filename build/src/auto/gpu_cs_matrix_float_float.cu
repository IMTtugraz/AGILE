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

// $Id: gpu_cs_matrix.cu.in 452 2011-05-31 12:00:18Z freiberger $

/* This file was generated automatically by CMake. You have to modify '/home2/GIT/AGILE/src/gpu_cs_matrix.cu.in' if you want to make changes. */

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

// Unfortunately, textures have to have file scope, which is why we have to
// use this ugly preprocessor definition.
#define AGILE_TEXTURE agile_cs_matrix_texture_floatfloat
texture<agile::to_tuple_type<float >::type> AGILE_TEXTURE ;

#include "gpu_cs_matrix.ipp"

namespace agile
{
  namespace detail
  {
    // ==================== explicit instantiations ====================
    template void multiplyCRS<float, float >(
      const GPUVector<unsigned>& row_nnz, const GPUVector<unsigned>& row_offset,
      const GPUVector<unsigned>& column_index,
      const GPUVector<float >& matrix_data, const GPUVector<float >& x,
      GPUVector<typename promote<float, float >::type>& y);

    template void multiplyCRS<float, float >(
      const GPUVector<float >& x,
      const GPUVector<unsigned>& row_nnz, const GPUVector<unsigned>& row_offset,
      const GPUVector<unsigned>& column_index,
      const GPUVector<float >& matrix_data,
      GPUVector<typename promote<float, float >::type>& y);

  } // namespace detail
} // namespace agile

// End of $Id: gpu_cs_matrix.cu.in 452 2011-05-31 12:00:18Z freiberger $.
