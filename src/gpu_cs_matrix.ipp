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

// $Id: gpu_cs_matrix.ipp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/exception.hpp"
#include "agile/gpu_config.hpp"
#include "agile/gpu_complex.hpp"
#include "agile/gpu_cs_matrix.hpp"
#include "agile/gpu_environment.hpp"
#include "agile/gpu_memory.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

// Unfortunately, textures have to have file scope, which is why we have to
// use this ugly preprocessor definition.
#ifndef AGILE_TEXTURE
  #error CUDA textures have to have file scope. Thus, you have to define AGILE_TEXTURE.
#endif

//! \brief A*x = y with a CRS matrix A and a GPU vector x (GPU function).
//!
//! Every thread computes the inner product of one matrix row with the vector,
//! that is the \p i-th thread computes \f$ y(i) \leftarrow A(i, :) * x \f$.
//! \p TType2 is used to deduce the type of the vector \p x.
//! \param[in] row_nnz Amount of non-zero entries in each row.
//! \param[in] row_offset Vector containing the offset in \p column_index and
//!            \p data for each row.
//! \param[in] column_index The column index for every entry.
//! \param[in] data The matrix entries.
//! \param[in] num_rows The amount of rows.
//! \param[in] stride The stride used to when aligning the matrix entries for
//!            the GPU.
//! \param[out] y The multiplication result.
template <typename TType1, typename TType2, typename TType3>
__global__ void multiplyCRSMatrixVector_GPU(
  const unsigned* row_nnz, const unsigned* row_offset,
  const unsigned* column_index, const TType1* data, unsigned num_rows,
  unsigned stride, TType3* y)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  if (thread_id < num_rows)
  {
    // calculate the first and the last index to multiply
    unsigned element_index_start = row_offset[thread_id];
    unsigned element_index_end = row_offset[thread_id]
                                 + stride * row_nnz[thread_id];

    TType3 sum(0);
    for (unsigned index = element_index_start;
         index < element_index_end; index += stride)
    {
      unsigned crs_index = column_index[index];
      TType2 x_value = tex1Dfetch(AGILE_TEXTURE, crs_index);
      sum += data[index] * x_value;
    }
    y[thread_id] = sum;
  }
}

//! \brief A^H*x = y with a CRS matrix A and a GPU vector x (GPU function).
//!
//! The transposed/hermitian matrix-vector product is somehow problematic when
//! the matrix is stored in CRS format. The reason is that with CRS it is
//! efficient to traverse along a row but it is rather complicated to go along
//! a column.
//! The compromise is the following: The \p i-th thread performs the
//! operation \f$ A(i, :)*x(i) + A(i+N, :)*x(i+N) + A(i+2N, :)*x(i+2N) + ...\f$,
//! where \f$ N \f$ is a natural number larger than zero, and stores the output
//! in the vector \f$ y_i \f$. Afterwards these vectors have to be summed to
//! get the final result.
//! When calling this function, the total number of threads has to be equal
//! to the row stride \f$ N \f$ and the output vector \p y has to be of size
//! (\p num_columns x \p row_stride ).
//! \param[in] x The vector to be multiplied with \f$ A^H \f$.
//! \param[in] row_nnz Amount of non-zero entries in each row of \p A.
//! \param[in] row_offset Vector containing the offset in \p column_index and
//!            \p data for each row.
//! \param[in] column_index The column index for every entry.
//! \param[in] data The matrix entries.
//! \param[in] num_rows The amount of rows.
//! \param[in] element_stride The stride used to when aligning the matrix
//!            entries for the GPU.
//! \param[in] row_stride The row stride \f$ N \f$.
//! \param[out] y The multiplication result.
template <typename TType1, typename TType2, typename TType3>
__global__ void multiplyHermCRSMatrixVector_GPU(
  const TType1* x,
  const unsigned* row_nnz, const unsigned* row_offset,
  const unsigned* column_index, const TType2* data, unsigned num_rows,
  unsigned element_stride, unsigned row_stride, TType3* y)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  for (unsigned row_counter = thread_id; row_counter < num_rows;
       row_counter += row_stride)
  {
    // get the vector element
    TType1 x_value = x[row_counter];
    
    // calculate the first and the last index to multiply
    unsigned element_index_start = row_offset[row_counter];
    unsigned element_index_end = row_offset[row_counter]
                                 + element_stride * row_nnz[row_counter];

    // add the product A(i+row_counter*N, :) * x[row_counter] to the temporary
    // vector y[i]
    for (unsigned index = element_index_start;
         index < element_index_end; index += element_stride)
    {
      unsigned crs_index = column_index[index];
      y[crs_index * row_stride + thread_id]
        += agile::conj(data[index]) * x_value;
    }
  }
}

namespace agile
{
  namespace detail
  {
    using namespace StandardException;

    //! \brief Multiply a GPU CRS matrix with a GPU vector (host function).
    template <typename TType1, typename TType2>
    void multiplyCRS(const GPUVector<unsigned>& row_nnz,
                     const GPUVector<unsigned>& row_offset,
                     const GPUVector<unsigned>& column_index,
                     const GPUVector<TType1>& matrix_data,
                     const GPUVector<TType2>& x,
                     GPUVector<typename promote<TType1, TType2>::type>& y)
    {
      AGILE_ASSERT(row_nnz.size() == row_offset.size(),
                    ExceptionSizeMismatch("row_nnz", "row_offset",
                                          row_nnz.size(), row_offset.size()));
      AGILE_ASSERT(column_index.size() == matrix_data.size(),
                    ExceptionSizeMismatch(
                      "column_index", "matrix_data",
                      column_index.size(), matrix_data.size()));
      AGILE_ASSERT(row_nnz.size() == y.size(),
                    ExceptionSizeMismatch("row_nnz", "y",
                                          row_nnz.size(), y.size()));

      unsigned num_rows = row_nnz.size();
      // bind a texture to the x-vector as we need random access to it
      cudaBindTexture(
        0, AGILE_TEXTURE, x.data(),
        x.size() * sizeof(typename substitute_gpu_complex<TType2>::type));
      // calculate the grid size
      unsigned grid_size
        = (num_rows + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
          / GPUEnvironment::getMaxNumThreadsPerBlock();
      // call the kernel
      multiplyCRSMatrixVector_GPU<
        typename substitute_gpu_complex<TType1>::type,
        typename substitute_gpu_complex<TType2>::type,
        typename substitute_gpu_complex<
           typename promote<TType1, TType2>::type>::type>
      <<<grid_size, GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
        row_nnz.data(), row_offset.data(), column_index.data(),
        (const typename substitute_gpu_complex<TType1>::type*)
          matrix_data.data(),
        num_rows, CRS_BLOCK_SIZE,
        (typename substitute_gpu_complex<
           typename promote<TType1, TType2>::type>::type*)y.data());

      // free the texture again
      cudaUnbindTexture(AGILE_TEXTURE);
    }

    //! \brief Hermitian matrix-vector product with CRS matrix (host fn).
    template <typename TType1, typename TType2>
    void multiplyCRS(const GPUVector<TType1>& x,
                     const GPUVector<unsigned>& row_nnz,
                     const GPUVector<unsigned>& row_offset,
                     const GPUVector<unsigned>& column_index,
                     const GPUVector<TType2>& matrix_data,
                     GPUVector<typename promote<TType1, TType2>::type>& y)
    {
      AGILE_ASSERT(row_nnz.size() == row_offset.size(),
                    ExceptionSizeMismatch("row_nnz", "row_offset",
                                          row_nnz.size(), row_offset.size()));
      AGILE_ASSERT(column_index.size() == matrix_data.size(),
                    ExceptionSizeMismatch(
                      "column_index", "matrix_data",
                      column_index.size(), matrix_data.size()));
      AGILE_ASSERT(row_nnz.size() == x.size(),
                    ExceptionSizeMismatch("row_nnz", "x",
                                          row_nnz.size(), x.size()));

      // create a temporary vector
//TODO: with csr matrix of size 1024x1024 the following line results in
//      cudaMalloc of 512 MB!!! (far too much for my gpu, gerald buchgraber)
      GPUVector<typename promote<TType1, TType2>::type> y_temp(
        y.size() * CRS_BLOCK_SIZE);
      unsigned num_rows = row_nnz.size();
      unsigned grid_size
        = (CRS_BLOCK_SIZE + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
          / GPUEnvironment::getMaxNumThreadsPerBlock();
      multiplyHermCRSMatrixVector_GPU<
        typename substitute_gpu_complex<TType1>::type,
        typename substitute_gpu_complex<TType2>::type,
        typename substitute_gpu_complex<
           typename promote<TType1, TType2>::type>::type>
      <<<grid_size, CRS_BLOCK_SIZE>>>(
        (const typename substitute_gpu_complex<TType1>::type*)x.data(),
        row_nnz.data(), row_offset.data(), column_index.data(),
        (const typename substitute_gpu_complex<TType2>::type*)
          matrix_data.data(),
        num_rows, CRS_BLOCK_SIZE, CRS_BLOCK_SIZE,
        (typename substitute_gpu_complex<
           typename promote<TType1, TType2>::type>::type*)y_temp.data());

      // sum the partial results
      addMemoryColumns(y_temp.data(), y.size(), CRS_BLOCK_SIZE, y.data());
    }

  } // namespace detail
} // namespace agile

// End of $Id: gpu_cs_matrix.ipp 476 2011-06-16 08:54:14Z freiberger $.
