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

// $Id: gpu_matrix_pitched.ipp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/exception.hpp"
#include "agile/gpu_config.hpp"
#include "agile/gpu_environment.hpp"
#include "agile/gpu_matrix_pitched.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

#ifndef TType2IsComplex
  #error TType2IsComplex is not defined.
#endif

// Unfortunately, textures have to have file scope, which is why we have to
// use this ugly preprocessor definition.
#ifndef AGILE_TEXTURE
  #error CUDA textures have to have file scope. Thus, you have to define AGILE_TEXTURE.
#endif

// ================================== products =================================
//! \brief Multiply a GPU matrix with a scalar (GPU function)
template <typename TType1, typename TType2, typename TType3>
__global__ void multiplyMatrixScalar_GPU(
  const TType1 alpha, const TType2* x, TType3* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  if (thread_id < size)
    y[thread_id] = alpha * x[thread_id];
}

//! \brief Multiply a dense GPU matrix with a GPU vector (GPU function).
//!
//! Every block calculates the product of one matrix row with the vector \p x.
//! \param[in] num_columns The amount of columns of the matrix.
//! \param[in] aligned_row_size The amount of entries per row after alignment.
//! \param[in] matrix Vector containing the matrix elements.
//! \param[in] x The vector which has to be multiplied with the matrix.
//! \param[out] y The multiplication result.
template <typename TType1, typename TType2, typename TType3>
__global__ void multiplyMatrixVector_GPU(
  unsigned num_columns, unsigned aligned_row_size, const TType1* matrix,
  const TType2* x, TType3* y)
{
  // hold the partial sums
  __shared__ TType3 partial_row_sum[agile::MAX_NUM_THREADS_PER_BLOCK];

  unsigned row_offset = aligned_row_size * blockIdx.x;
  // calculate the product of the current row with the vector x with a stride
  // equal to the amount of threads
  TType3 sum(0);
  for (unsigned index = threadIdx.x; index < num_columns; index += blockDim.x)
    sum += matrix[row_offset + index] * x[index];

  partial_row_sum[threadIdx.x] = sum;
  // reduce the partial sums
  for (unsigned counter = blockDim.x >> 1; counter != 0; counter >>= 1)
  {
    __syncthreads();
    if (threadIdx.x < counter)
      partial_row_sum[threadIdx.x] += partial_row_sum[threadIdx.x + counter];
  }
  // write the result into the output vector
  if (threadIdx.x == 0)
    y[blockIdx.x] = partial_row_sum[0];
}

//! \brief Hermitian matrix-vector product (GPU function).
//!
//! Implements the product \f$ A^H x \f$. Every thread calculates the inner
//! product over one column of \f$ A \f$.
//! \param[in] num_rows The amount of rows of the matrix.
//! \param[in] aligned_row_size The amount of entries per row after alignment.
//! \param[in] num_columns The amount of columns of the matrix.
//! \param[in] x The vector to be multiplied with \f$ A^H \f$.
//! \param[in] matrix Vector containing the matrix elements.
//! \param[out] y The multiplication result.
template <typename TType1, typename TType2, typename TType3>
__global__ void multiplyHermMatrixVector_GPU(
   unsigned num_rows, unsigned aligned_row_size, unsigned num_columns,
   const TType1* x, const TType2* matrix, TType3* y)
{
  __shared__ TType1 local_x[agile::MAX_NUM_THREADS_PER_BLOCK];
  // compute the index of the column which the first thread in this block
  // has to process
  unsigned column_offset = blockDim.x * blockIdx.x;
  while (column_offset < num_columns)
  {
    // compute the index of the column which is currently processed by this thread
    unsigned column_index = column_offset + threadIdx.x;
    TType3 sum(0);
    // loop over the rows in blocks
    for (unsigned row_offset = 0; row_offset < num_rows; row_offset += blockDim.x)
    {
      // cache x
      __syncthreads();
      if (row_offset + threadIdx.x < num_rows)
        local_x[threadIdx.x] = x[row_offset + threadIdx.x];
      __syncthreads();

      if (column_index < num_columns)
      {
        unsigned row_end = row_offset + blockDim.x;
        if (row_end > num_rows)
          row_end = num_rows;
        for (unsigned row = row_offset, x_counter = 0; row < row_end;
             ++row, ++x_counter)
        {
          sum += agile::conj(matrix[row * aligned_row_size + column_index])
                 * local_x[x_counter];
        }
      }
    }

    if (column_index < num_columns)
      y[column_index] = sum;

    column_offset += blockDim.x * gridDim.x;
  }

}

//! \brief Element-wise matrix-matrix product (GPU function).
template <typename TType1, typename TType2, typename TType3>
__global__ void multiplyMatrixMatrixElementwise_GPU(
  const TType1* matrix1, const TType2* matrix2,
  TType3* matrix3, unsigned num_elements)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  if (thread_id < num_elements)
    matrix3[thread_id] = matrix1[thread_id] * matrix2[thread_id];
}

// =============================== interpolation ===============================
//! \brief Bilinear interpolation (GPU function).
template <typename TType1, typename TType2>
__global__ void interp2d_GPU
(const TType2* pos, TType1* res, unsigned num_elements)
{
  unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < num_elements)
  {
      typename agile::to_tuple_type<TType1>::texture_type res_helper
               = tex2D(AGILE_TEXTURE_2D, real(pos[idx])+0.5, imag(pos[idx])+0.5);
      res[idx] = agile::to_tuple_type<TType1>::texture2type(res_helper);
  }
}


namespace agile
{
  using namespace StandardException;

  //! \brief Multiply a dense GPU matrix with a GPU vector (host function).
  template <typename TType1, typename TType2>
  void multiply(const GPUMatrixPitched<TType1>& A, const GPUVector<TType2>& x,
                GPUVector<typename promote<TType1, TType2>::type>& y)
  {
    AGILE_ASSERT(A.getNumColumns() == x.size(),
                  ExceptionSizeMismatch("columns of A", "x",
                                        A.getNumColumns(), x.size()));
    AGILE_ASSERT(A.getNumRows() == y.size(),
                  ExceptionSizeMismatch("rows of A", "y",
                                        A.getNumRows(), y.size()));

    multiplyMatrixVector_GPU<<<A.getNumRows(),
                               GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      A.getNumColumns(), A.getPitchElements(),
      (const typename substitute_gpu_complex<TType1>::type*)A.data(),
      (const typename substitute_gpu_complex<TType2>::type*)x.data(),
      (typename substitute_gpu_complex<
         typename promote<TType1, TType2>::type>::type*)y.data());
  }

  //! \brief \f$ A^H x \f$ with dense GPU matrix A and GPU vector x (host fn).
  template <typename TType1, typename TType2>
  void multiply(const GPUVector<TType1>& x, const GPUMatrixPitched<TType2>& A,
                GPUVector<typename promote<TType1, TType2>::type>& y)
  {
    AGILE_ASSERT(A.getNumColumns() == y.size(),
                  ExceptionSizeMismatch("columns of A", "y",
                                        A.getNumColumns(), y.size()));
    AGILE_ASSERT(A.getNumRows() == x.size(),
                  ExceptionSizeMismatch("rows of A", "x",
                                        A.getNumRows(), x.size()));

    unsigned num_columns = A.getNumColumns();

    unsigned grid_size
      = (num_columns + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    if (grid_size > 30)
        grid_size = 30;
    multiplyHermMatrixVector_GPU<
      typename substitute_gpu_complex<TType1>::type,
      typename substitute_gpu_complex<TType2>::type,
      typename substitute_gpu_complex<
        typename promote<TType1, TType2>::type>::type>
      <<<grid_size, GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      A.getNumRows(), A.getPitchElements(), A.getNumColumns(),
      (const typename substitute_gpu_complex<TType1>::type*)x.data(),
      (const typename substitute_gpu_complex<TType2>::type*)A.data(),
      (typename substitute_gpu_complex<
         typename promote<TType1, TType2>::type>::type*)y.data());
  }

  //! \brief Element-wise matrix-matrix product (host function).
  template <typename TType1, typename TType2>
  void multiplyElementwise(const GPUMatrixPitched<TType1>& A,
                           const GPUMatrixPitched<TType2>& B,
                           GPUMatrixPitched<typename promote<TType1, TType2>::type>& Z)
  {
    AGILE_ASSERT(A.getNumRows() == B.getNumRows(),
                  ExceptionSizeMismatch("rows of A", "rows of B",
                                        A.getNumRows(), B.getNumRows()));
    AGILE_ASSERT(A.getNumColumns() == B.getNumColumns(),
                  ExceptionSizeMismatch("columns of A", "columns of B",
                                        A.getNumColumns(), B.getNumColumns()));
    AGILE_ASSERT(A.getNumRows() == Z.getNumRows(),
                  ExceptionSizeMismatch("rows of A", "rows of Z",
                                        A.getNumRows(), Z.getNumRows()));
    AGILE_ASSERT(A.getNumColumns() == Z.getNumColumns(),
                  ExceptionSizeMismatch("columns of A", "columns of Z",
                                        A.getNumColumns(), Z.getNumColumns()));

    unsigned a_size = A.getPitchElements() * A.getNumRows();
    unsigned grid_size
      = (a_size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    multiplyMatrixMatrixElementwise_GPU
      <<<grid_size, GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)A.data(),
      (const typename substitute_gpu_complex<TType2>::type*)B.data(),
      (typename substitute_gpu_complex<
         typename promote<TType1, TType2>::type>::type*)Z.data(),
      a_size);
  }

  //! \brief Multiply a GPU matrix with a scalar (host function).
  template <typename TType1, typename TType2>
  void scale(const TType1& alpha, const GPUMatrixPitched<TType2>& A,
             GPUMatrixPitched<typename promote<TType1, TType2>::type>& B)
  {
    AGILE_ASSERT(A.getNumColumns() == B.getNumColumns(),
                  ExceptionSizeMismatch("A", "B",
                                        A.getNumColumns(), B.getNumColumns()));

    AGILE_ASSERT(A.getNumRows() == B.getNumRows(),
                  ExceptionSizeMismatch("A", "B",
                                        A.getNumRows(), B.getNumRows()));
    unsigned num_elements = A.getPitchElements() * A.getNumRows();
    unsigned block_dim = GPUEnvironment::getMaxNumThreadsPerBlock();
    unsigned grid_dim = (num_elements + block_dim - 1) / block_dim;

    multiplyMatrixScalar_GPU<<<grid_dim, block_dim>>>(
      *(const typename substitute_gpu_complex<TType1>::type*)(&alpha),
      (const typename substitute_gpu_complex<TType2>::type*)A.data(),
      (typename substitute_gpu_complex<
         typename promote<TType1, TType2>::type>::type*)B.data(),
      num_elements);
  }


#if !TType2IsComplex

  //! \brief Bilinear interpolation (host function).
  template <typename TType1, typename TType2>
  void interp2d(const GPUMatrixPitched<TType1>& M,
    const GPUVector<std::complex<TType2> >& pos, GPUVector<TType1>& res)
  {
    AGILE_ASSERT(pos.size() == res.size(),
                  ExceptionSizeMismatch("pos", "res",
                                        pos.size(), res.size()));

    AGILE_TEXTURE_2D.filterMode = cudaFilterModeLinear;

#if CUDA_VERSION < 2020
    // in cuda version lower than 2.2 cuda arrays
    // are necessary for interpolation
    cudaArray *dev_source;
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<
      typename substitute_gpu_complex<TType1>::type >();

    CUDA_SAFE_CALL(cudaMallocArray(&dev_source, &channelDesc,
      M.getNumColumns(), M.getNumRows()));

    CUDA_SAFE_CALL(cudaMemcpy2DToArray(dev_source, 0, 0,
      (const typename substitute_gpu_complex<TType1>::type*)M.data(),
      M.getPitch(), 
      M.getNumColumns() * sizeof(typename substitute_gpu_complex<TType1>::type),
      M.getNumRows(), cudaMemcpyDeviceToDevice));

    // bind a 2d texture to the matrix as we need random access to it
    CUDA_SAFE_CALL(cudaBindTextureToArray(AGILE_TEXTURE_2D, dev_source));

#else

    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<
      typename substitute_gpu_complex<TType1>::type >();
    CUDA_SAFE_CALL(cudaBindTexture2D(0, &AGILE_TEXTURE_2D,
      (const typename substitute_gpu_complex<TType1>::type*)M.data(),
      &channelDesc, M.getNumColumns(), M.getNumRows(), M.getPitch()));

#endif

    unsigned block_dim = GPUEnvironment::getMaxNumThreadsPerBlock();
    unsigned grid_dim = (res.size() + block_dim - 1) / block_dim;

    interp2d_GPU<<<grid_dim, block_dim>>>(
      (const typename substitute_gpu_complex<std::complex<TType2> >::type*)
        pos.data(),
      (typename substitute_gpu_complex<TType1>::type*)res.data(),
      res.size());

    // unbind texture and reset changes to global texture
    CUDA_SAFE_CALL(cudaUnbindTexture(AGILE_TEXTURE_2D));
    AGILE_TEXTURE_2D.filterMode = cudaFilterModePoint;

#if CUDA_VERSION < 2020
    CUDA_SAFE_CALL(cudaFreeArray(dev_source));
#endif
  }

#endif  // TType2IsComplex

} // namespace agile

// End of $Id: gpu_matrix_pitched.ipp 476 2011-06-16 08:54:14Z freiberger $.
