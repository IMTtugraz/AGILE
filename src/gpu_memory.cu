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

// $Id: gpu_memory.cu 452 2011-05-31 12:00:18Z freiberger $

#include "agile/gpu_config.hpp"
#include "agile/gpu_complex.hpp"
#include "agile/gpu_environment.hpp"
#include "agile/gpu_memory.hpp"
#include "agile/gpu_type_traits.hpp"

#include <complex>
#include <cuda.h>

//! \brief Initialize memory on the GPU with a constant (GPU function).
template <typename TType>
__global__ void initializeGPUMemory_GPU(TType* data, unsigned size,
                                        TType value)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  if (thread_id < size)
    data[thread_id] = value;
}

//! \brief Add the columns of a rectangular memory area (GPU function).
//!
//! The reason why this function is located in \p gpu_memory.cu
//! (vs. \p gpu_matrix.cu or \p gpu_vector.cu, where it would fit, too) is
//! that it requires only one template argument and currently it is quite
//! difficult to automatically instantiate templates with a different number
//! of template parameters.
//! Call this kernel with as many threads as there are elements in \p input
//! which is equal to (\p num_rows * \p num_columns ).
//! \param[in] input Pointer to the memory whose columns are to be summed.
//! \param[in] num_rows The number of rows in \p input and \p output,
//! respectively.
//! \param[in] num_columns The number of columns in \p input.
//! \param[out] output The sum of the input vectors. Has to be of size
//! \p num_rows.
template <typename TType>
__global__ void addMemoryColumns_GPU(
  TType* input, unsigned num_rows, unsigned num_columns, TType* output)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  // calculate on which row we have to operate on
  unsigned row = thread_id / num_columns;
  if (row < num_rows)
  {
    // calculate which column we have to deal with
    unsigned column = thread_id % num_columns;
  
    for (unsigned counter = num_columns >> 1; counter != 0; counter >>= 1)
    {
      if (column < counter)
      {
        // input[row * num_columns + column]
        //   += input[row * num_columns + column + counter]
        input[thread_id] += input[thread_id + counter];
      }
      __syncthreads();
    }
  
    // assign the result to the final output vector
    if (column == 0)
      output[row] = input[thread_id];
  }
}

namespace agile
{
  //! \brief Set GPU memory to a constant (host function).
  //!
  //! \param[in] data Pointer to the GPU to initialize.
  //! \param[in] size The amount of elements to initialize.
  //! \param[in] value The value for the elements.
  template <typename TType>
  void initializeGPUMemory(TType* data, unsigned size, const TType& value)
  {
    unsigned grid_size = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
                         / GPUEnvironment::getMaxNumThreadsPerBlock();
    initializeGPUMemory_GPU<<<grid_size,
                              GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (typename substitute_gpu_complex<TType>::type*)data,
      size, (const typename substitute_gpu_complex<TType>::type&)value);
  }

  //! \brief Add the columns of a rectangular memory area (host function).
  template <typename TType>
  void addMemoryColumns(TType* input, unsigned num_rows, unsigned num_columns,
                        TType* output)
  {
    unsigned grid_size = (num_rows * num_columns
                          + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
                         / GPUEnvironment::getMaxNumThreadsPerBlock();
    addMemoryColumns_GPU<<<grid_size,
                           GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (typename substitute_gpu_complex<TType>::type*)input,
      num_rows, num_columns,
      (typename substitute_gpu_complex<TType>::type*)output);
  }

  // ========================= explicit instantiation =========================
  template void initializeGPUMemory<unsigned char>(
    unsigned char*, unsigned, const unsigned char&);
  template void initializeGPUMemory<int>(
    int*, unsigned, const int&);
  template void initializeGPUMemory<unsigned>(
    unsigned*, unsigned, const unsigned&);
  template void initializeGPUMemory<float>(
    float*, unsigned, const float&);
  template void initializeGPUMemory<std::complex<float> >(
    std::complex<float>*, unsigned, const std::complex<float>&);

  template void addMemoryColumns<float>(
    float*, unsigned, unsigned, float*);
  template void addMemoryColumns<std::complex<float> >(
    std::complex<float>*, unsigned, unsigned, std::complex<float>*);

} // namespace agile

// End of $Id: gpu_memory.cu 452 2011-05-31 12:00:18Z freiberger $.
