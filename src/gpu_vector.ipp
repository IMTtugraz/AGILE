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

// $Id: gpu_vector.ipp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/gpu_config.hpp"
#include "agile/gpu_complex.hpp"
#include "agile/gpu_environment.hpp"
#include "agile/gpu_host_vector.hpp"
#include "agile/gpu_type_traits.hpp"
#include "agile/gpu_vector.hpp"
#include <cuda.h>
#include <stdio.h>
#include <math.h>


//! \brief Add a scaled GPU vector to another GPU vector (GPU function).
//!
//! This function multiplies the vector \p y with a constant \p scale and adds
//! the vector \p x. The result is stored in vector \p z:
//! z <- x + scale * y
template <typename TType1, typename TType2, typename TType3, typename TType4>
__global__ void addScaledVector_GPU(
  const TType1* x, const TType2 scale, const TType3* y, TType4* z,
  unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  while (thread_id < size)
  {
    z[thread_id] = x[thread_id] + scale * y[thread_id];
    thread_id += blockDim.x*gridDim.x; 
  }
}

//! \brief Add two vectors (GPU function).
//!
//! z <- x + y
template <typename TType1, typename TType2, typename TType3>
__global__ void addVector_GPU(const TType1* x, const TType2* y, TType3* z,
                              unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    z[thread_id] = x[thread_id] + y[thread_id];
    thread_id += blockDim.x*gridDim.x; 
  }
}

//! \brief Divide a scalar by a GPU vector (elementwise; GPU function).
//!
//! y[i] <- alpha / x[i]
template <typename TType1, typename TType2, typename TType3>
__global__ void divideVector_GPU(
  const TType1 alpha, const TType2* x, TType3* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    y[thread_id] = alpha / x[thread_id];
    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Calculate the biinear form of two vectors (GPU function).
//!
//! The bilinear form is calculated only partially on the GPU. In the first
//! step the i-th thread sums the products of the elements i, i + S, i + 2*S...,
//! where the stride S is equal to blockDim * gridDim (i.e. the total number
//! of threads).
//! In the second step each thread block reduces its partial sums in a parallel
//! manner and stores the result to the vector \p w.
//! \param[in] x First vector.
//! \param[in] y Second vector.
//! \param[out] w Partial bilinear form which has to be reduced by the CPU.
//!               The vector has to be of length \p grid_size.
//! \param[in] size Size of the vectors.
template <typename TType1, typename TType2, typename TType3>
__global__ void getBilinearForm_GPU(const TType1* x, const TType2* y,
                                    TType3* partial_result, unsigned size)
{
  __shared__ TType3 partial_sum[agile::MAX_NUM_THREADS_PER_BLOCK];

  // sum the products of the vector elements with a stride
  // threads_per_block * grid_size
  TType3 sum(0);
  for (unsigned counter = blockDim.x * blockIdx.x + threadIdx.x;
       counter < size; counter += blockDim.x * gridDim.x)
    sum += x[counter] * y[counter];

  // store this partial sum within this thread block
  partial_sum[threadIdx.x] = sum;

  // Reduce the partial sums within this thread block by parallel sums. In the
  // first iteration the threads 0 .. blockDim.x/2-1 add the contribution
  // of the threads blockDim.x/2 .. blockDim.x-1. In the second iteration
  // the threads 0 .. blockDim.x/4-1 add the value of
  // blockDim.x/4..blockDim.x/2 and so on until in the last iteration
  // thread 0 adds the value of thread 1.
  for (unsigned counter = blockDim.x >> 1; counter != 0; counter >>= 1)
  {
    __syncthreads();
    if (threadIdx.x < counter)
      partial_sum[threadIdx.x] += partial_sum[threadIdx.x + counter];
  }
  // thread 0 writes the output to the vector partial_result
  if (threadIdx.x == 0)
    partial_result[blockIdx.x] = partial_sum[0];
}

//! \brief Calculate the scalar product of two vectors (GPU function).
//!
//! The scalar product is calculated only partially on the GPU. In the first
//! step the i-th thread sums the products of the elements i, i + S, i + 2*S...,
//! where the stride S is equal to blockDim * gridDim (i.e. the total number
//! of threads).
//! In the second step each thread block reduces its partial sums in a parallel
//! manner and stores the result to the vector \p w.
//! \param[in] x First vector.
//! \param[in] y Second vector.
//! \param[out] w Partial scalar product which has to be reduced by the CPU.
//!               The vector has to be of length \p grid_size.
//! \param[in] size Size of the vectors.
template <typename TType1, typename TType2, typename TType3>
__global__ void getScalarProduct_GPU(const TType1* x, const TType2* y,
                                     TType3* partial_result, unsigned size)
{
  __shared__ TType3 partial_sum[agile::MAX_NUM_THREADS_PER_BLOCK];

  // sum the products of the vector elements with a stride
  // threads_per_block * grid_size
  TType3 sum(0);
  for (unsigned counter = blockDim.x * blockIdx.x + threadIdx.x;
       counter < size; counter += blockDim.x * gridDim.x)
    sum += agile::conj(x[counter]) * y[counter];

  // store this partial sum within this thread block
  partial_sum[threadIdx.x] = sum;

  // Reduce the partial sums within this thread block by parallel sums. In the
  // first iteration the threads 0 .. blockDim.x/2-1 add the contribution
  // of the threads blockDim.x/2 .. blockDim.x-1. In the second iteration
  // the threads 0 .. blockDim.x/4-1 add the value of
  // blockDim.x/4..blockDim.x/2 and so on until in the last iteration
  // thread 0 adds the value of thread 1.
  for (unsigned counter = blockDim.x >> 1; counter != 0; counter >>= 1)
  {
    __syncthreads();
    if (threadIdx.x < counter)
      partial_sum[threadIdx.x] += partial_sum[threadIdx.x + counter];
  }
  // thread 0 writes the output to the vector partial_result
  if (threadIdx.x == 0)
    partial_result[blockIdx.x] = partial_sum[0];
}

//! \brief Extract the imaginary part of a GPU vector (GPU function).
//!
//! y <- imag(x)
template <typename TType1, typename TType2>
__global__ void imag_GPU(
  const TType1* x, TType2* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    y[thread_id] = agile::imag(x[thread_id]);
    thread_id += blockDim.x*gridDim.x;
  }
}


//////////////////////////////////////////////////////////////////////////
//! \brief Bilinear interpolation (GPU function).
template <typename TType1, typename TType2>
__global__ void interpolate2DVectorColumnMajor_GPU
(const TType2* pos, TType1* res, unsigned num_elements)
{
  //typename agile::to_tuple_type<TType1>::texture_type res_helper = 0;
  unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;

  if (idx < num_elements) {
     typename agile::to_tuple_type<TType1>::texture_type res_helper
       = tex2D(AGILE_TEXTURE_2D, imag(pos[idx])+0.5, real(pos[idx])+0.5);
     res[idx] = agile::to_tuple_type<TType1>::texture2type(res_helper);
  }
}


//! \brief Bilinear interpolation (GPU function).
template <typename TType1, typename TType2>
__global__ void interpolate2DVectorRowMajor_GPU
(const TType2* pos, TType1* res, unsigned num_elements)
{
  unsigned idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < num_elements) {
    typename agile::to_tuple_type<TType1>::texture_type res_helper
     = tex2D(AGILE_TEXTURE_2D, real(pos[idx])+0.5, imag(pos[idx])+0.5);
    res[idx] = agile::to_tuple_type<TType1>::texture2type(res_helper);
  }
}



//! \brief Multiply conj(x) and y element-wise and store result in z (GPU function).
//!
//! z[i] <- conj(x[i]) * y[i]
template <typename TType1, typename TType2, typename TType3>
__global__ void multiplyConjElementwise_GPU(
  const TType1* x, const TType2* y, TType3* z, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  while (thread_id < size)
  {
    z[thread_id] = agile::conj(x[thread_id]) * y[thread_id];
    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Multiply x and y element-wise and store result in z (GPU function).
//!
//! z[i] <- x[i] * y[i]
template <typename TType1, typename TType2, typename TType3>
__global__ void multiplyElementwise_GPU(
  const TType1* x, const TType2* y, TType3* z, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  while (thread_id < size)
  {
    z[thread_id] = x[thread_id] * y[thread_id];
    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Divide x and y element-wise and store result in z (GPU function).
//!
//! z[i] <- x[i] * y[i]
template <typename TType1, typename TType2, typename TType3>
__global__ void divideElementwise_GPU(
  const TType1* x, const TType2* y, TType3* z, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  while (thread_id < size)
  { 
    z[thread_id] = x[thread_id] / y[thread_id];
    thread_id += blockDim.x*gridDim.x;
  }
}


//! \brief Extract the real part of a GPU vector (GPU function).
//!
//! y <- real(x)
template <typename TType1, typename TType2>
__global__ void real_GPU(
  const TType1* x, TType2* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    y[thread_id] = agile::real(x[thread_id]);
    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Multiply a GPU vector with a scalar (GPU function).
//!
//! y <- alpha * x
template <typename TType1, typename TType2, typename TType3>
__global__ void scale_GPU(
  const TType1 alpha, const TType2* x, TType3* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    y[thread_id] = alpha * x[thread_id];
    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Set every vector element to a constant (GPU function).
//!
//! x[i] <- value
template <typename TType>
__global__ void setVectorConstant_GPU(
  const TType value, TType* x, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    x[thread_id] = value;
    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Compute the square root of a GPU vector (GPU function).
//!
//! \p TType2 is the type cast in order to call the sqrt function.
//! y[i] <- sqrt(x[i])
//template <typename TType1, typename TType2>
template <typename TType1>
__global__ void sqrt_GPU(
  const TType1* x, TType1* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {

    y[thread_id] = agile::cusqrt(x[thread_id]);
    //y[thread_id] = (x[thread_id]);

    //printf("\n complexsqrt: thread_id:%d - x[thread_id]: %f - y[thread_id]: %f",thread_id,x[thread_id],y[thread_id]);
    //printf("\n sqrt(2): %f",agile::cusqrt(TType1(2.0)));
    thread_id += blockDim.x*gridDim.x;
  }
   // y[thread_id] = sqrtf(x[thread_id]);
//    y[thread_id] = sqrt(TType2(x[thread_id]));
}

/*
//! \brief Compute Average over last dimension
//!
//! z = mean(x,dim)
template <typename TType1, typename TType2>
__global__ void averVector_GPU(const TType1* x,
                                  TType2* z, unsigned size, unsigned aver)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned img_id = size / aver;
  unsigned idx = 0;
  float tmp = 0.0;
  if (thread_id < img_id)
  {
    for(int idx = 0; idx < aver-1; idx++)
    {
      tmp += x[img_id + idx*img_id]; 
    }
    z[thread_id] = tmp / aver;
    thread_id += blockDim.x*gridDim.x;
  }
}
*/

//! \brief Subtract two vectors (GPU function).
//!
//! z <- x - y
template <typename TType1, typename TType2, typename TType3>
__global__ void subVector_GPU(const TType1* x, const TType2* y,
                                  TType3* z, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  while (thread_id < size)
  {
    z[thread_id] = x[thread_id] - y[thread_id];
    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Subtract a scaled GPU vector from another GPU vector (GPU function).
//!
//! This function multiplies the vector \p y with a constant \p scale and
//! subtracts this scaled vector from the vector \p x. The result is stored in
//! vector \p z:
//! z <- x - scale * y
template <typename TType1, typename TType2, typename TType3, typename TType4>
__global__ void subScaledVector_GPU(
  const TType1* x, const TType2 scale, const TType3* y, TType4* z,
  unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  while (thread_id < size)
  {
    z[thread_id] = x[thread_id] - scale * y[thread_id];
    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Compute the absolute sum of GPU vector elements (GPU function).
//!
//! \return \f$ \sum_i |x_i| \f$
template <typename TType1, typename TType2>
__global__ void sumabs_GPU(const TType1* x, TType2* partial_result,
                                 unsigned size)
{
  __shared__ TType2 partial_sum[agile::MAX_NUM_THREADS_PER_BLOCK];

  // sum the products of the vector elements with a stride
  // threads_per_block * grid_size
  TType1 sum(0);
  for (unsigned counter = blockDim.x * blockIdx.x + threadIdx.x;
       counter < size; counter += blockDim.x * gridDim.x)
    sum += agile::cusqrt(agile::norm(x[counter]));

  // store this partial sum within this thread block
  partial_sum[threadIdx.x] = agile::real(sum);

  // Reduce the partial sums within this thread block by parallel sums. In the
  // first iteration the threads 0 .. blockDim.x/2-1 add the contribution
  // of the threads blockDim.x/2 .. blockDim.x-1. In the second iteration
  // the threads 0 .. blockDim.x/4-1 add the value of
  // blockDim.x/4..blockDim.x/2 and so on until in the last iteration
  // thread 0 adds the value of thread 1.
  for (unsigned counter = blockDim.x >> 1; counter != 0; counter >>= 1)
  {
    __syncthreads();
    if (threadIdx.x < counter)
      partial_sum[threadIdx.x] += partial_sum[threadIdx.x + counter];
  }
  // thread 0 writes the output to the vector partial_result
  if (threadIdx.x == 0)
    partial_result[blockIdx.x] = partial_sum[0];
}


//! \brief Compute the sum of squares of GPU vector elements (GPU function).
//!
//! \return \f$ \sum_i x_i^2 \f$
template <typename TType1, typename TType2>
__global__ void sumofsquares_GPU(const TType1* x, TType2* partial_result,
                                 unsigned size)
{
  __shared__ TType2 partial_sum[agile::MAX_NUM_THREADS_PER_BLOCK];

  // sum the products of the vector elements with a stride
  // threads_per_block * grid_size
  TType1 sum(0);
  for (unsigned counter = blockDim.x * blockIdx.x + threadIdx.x;
       counter < size; counter += blockDim.x * gridDim.x)
    sum += agile::norm(x[counter]);

  // store this partial sum within this thread block
  partial_sum[threadIdx.x] = agile::real(sum);

  // Reduce the partial sums within this thread block by parallel sums. In the
  // first iteration the threads 0 .. blockDim.x/2-1 add the contribution
  // of the threads blockDim.x/2 .. blockDim.x-1. In the second iteration
  // the threads 0 .. blockDim.x/4-1 add the value of
  // blockDim.x/4..blockDim.x/2 and so on until in the last iteration
  // thread 0 adds the value of thread 1.
  for (unsigned counter = blockDim.x >> 1; counter != 0; counter >>= 1)
  {
    __syncthreads();
    if (threadIdx.x < counter)
      partial_sum[threadIdx.x] += partial_sum[threadIdx.x + counter];
  }
  // thread 0 writes the output to the vector partial_result
  if (threadIdx.x == 0)
    partial_result[blockIdx.x] = partial_sum[0];
}

//! \brief  Conjugate GPU vector (GPU function).
//!
//! z[i] <- conj(x[i])
template <typename TType1>
__global__ void conjVector_GPU(const TType1* x, TType1* z, unsigned size)
{
   unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
   while (thread_id < size)
   {
     z[thread_id] = agile::conj(x[thread_id]);
     thread_id += blockDim.x*gridDim.x;
   }
}


//! \brief  absolute-value of GPU vector (GPU function).
//!
//! y[i] <- abs(x[i])
template <typename TTypeFloat, typename TType1, typename TType2>
__global__ void abs_GPU(
  const TType1* x, TType2* y, unsigned size)
{
   unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

   while (thread_id < size)
   {
     y[thread_id] = TType2(sqrt(TTypeFloat(agile::norm(x[thread_id]))));
     thread_id += blockDim.x*gridDim.x;
   }
}


//! \brief  phase-value of GPU vector (GPU function).
//!
//! y[i] <- arctan(imag/real)
template <typename TTypeFloat, typename TType1, typename TType2>
__global__ void phase_GPU(
  const TType1* x, TType2* y, unsigned size)
{
   unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

   while (thread_id < size)
   {
     y[thread_id] = TType2(atan2( x[thread_id].imag(), x[thread_id].real()));
     thread_id += blockDim.x*gridDim.x;
   }
}

//! \brief exponential of GPU vector (GPU function).
//!
//! y[i] <- exp(y[i])
template <typename TType1>
__global__ void exp_GPU(
  const TType1* x, TType1* y, unsigned size)
{
   unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

   while (thread_id < size)
   {
     y[thread_id] = agile::exp(x[thread_id]);
     thread_id += blockDim.x*gridDim.x;
   }
}


//! \brief Generate linearly spaced vector between a and b
//!
//! x = linspace(a,b)
template <typename TType1>
__global__ void linspace_GPU(
  TType1* x, float a, float b, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    x[thread_id] = a + thread_id * ((b-a) / (size-1));
    thread_id += blockDim.x*gridDim.x;
  }
}


template <typename TType1>
__global__ void meshgrid_GPU(TType1* m_y, const TType1* y, unsigned x_size, unsigned y_size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  unsigned y_idx = thread_id/x_size;

  extern __shared__ TType1 shValues_y[];
  unsigned cnt = threadIdx.x;
  while (cnt < y_size)
  {
    shValues_y[cnt] = y[cnt];
    cnt += blockDim.x;
  }
    __syncthreads();

  if (thread_id < y_size*x_size)
    m_y[thread_id] = shValues_y[y_idx];

}

//! \brief Compute the power of alpha for every element of a GPU Vector (GPU function).
//!
//! \p TType2 is the type cast in order to call the pow function.
//! y[i] <- pow(x[i],alpha)
template <typename TType1, typename TType2>
__global__ void pow_GPU(
  const TType1 alpha, const TType2* x, TType2* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    y[thread_id] = pow(x[thread_id] , alpha );
    thread_id += blockDim.x*gridDim.x;
  } 
}


//! \brief generate pattern of a vector.    (GPU function)
//! \brief z = abs(x)>0;
template <typename TTypeFloat, typename TType1, typename TType2>
__global__ void pattern_GPU(const TType1* x, TType2* z, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

 // TTypeFloat safe;
  while (thread_id < size)
  {
    //safe = TTypeFloat(sqrt(TTypeFloat(agile::norm(x[thread_id]))));
    z[thread_id] = (TTypeFloat(sqrt(TTypeFloat(agile::norm(x[thread_id])))) > 0) ? TType2(1) : TType2(0);
    thread_id += blockDim.x*gridDim.x;
  }
}


//! \brief Compute the difference between each value in a GPU vector (GPU function)
//! \brief last (x_size) value with first value
//!
//! y = diff(x)
template <typename TType1>
__global__ void diff_GPU(const unsigned dim, const unsigned x_size, const TType1* x, TType1* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    int val = thread_id-size+x_size;
    if (dim == 1)
      if ( ((thread_id+1) % x_size ) != 0) //last values not reached
        y[thread_id] = x[thread_id+1] - x[thread_id];
      else
        y[thread_id] = x[thread_id+1-x_size] - x[thread_id];

    if (dim == 2)
    {
      if ((val) < 0) //last values not reached
        y[thread_id] = x[thread_id+x_size] - x[thread_id];
      else
        //y[thread_id] = x[thread_id-val] - x[thread_id];
        y[thread_id] = x[thread_id-size+x_size] - x[thread_id];
    }
    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Compute the difference between each value in a GPU vector (GPU function)
//! \brief considering a 3D data layout.
//!
//! y = diff(x)
template <typename TType1>
__global__ void diff3_GPU(const unsigned dim, const unsigned x_size, const unsigned y_size, const TType1* x, TType1* y, unsigned size, bool borderWrap)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    int val = thread_id-size+x_size;
    int N = x_size * y_size;
    if (dim == 1)
      if ( ((thread_id+1) % x_size ) != 0) //last values not reached
        y[thread_id] = x[thread_id+1] - x[thread_id];
      else
        y[thread_id] = borderWrap ? x[thread_id+1-x_size] - x[thread_id] : 0;

    if (dim == 2)
    {
      int y_pos = (int)((thread_id - ((int)(thread_id / N))*N)/x_size)+1;
      if (y_pos % y_size != 0) //last values not reached
        y[thread_id] = x[thread_id+x_size] - x[thread_id];
      else
        y[thread_id] = borderWrap ? x[thread_id-N+x_size] - x[thread_id] : 0;
    }

    if (dim == 3)
    {
      val = thread_id - size + N;
      if ((val) < 0) 
        y[thread_id] = x[thread_id+N] - x[thread_id];
      else
        y[thread_id] = borderWrap ? x[thread_id-size+N] - x[thread_id] : 0;
    }

    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Compute the difference (sym. boundary cond.) between each value in a GPU vector (GPU function)
//! \brief considering a 3D data layout.
//!
//! y = diff(x)
template <typename TType1>
__global__ void diff3sym_GPU(const unsigned dim, const unsigned x_size, const unsigned y_size, const TType1* x, TType1* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    int val = thread_id-size+x_size;
    int N = x_size * y_size;
    if (dim == 1)
      if ( ((thread_id+1) % x_size ) != 0) //last values not reached
        y[thread_id] = x[thread_id+1] - x[thread_id];
      else
        y[thread_id] = x[thread_id-1] - x[thread_id];

    if (dim == 2)
    {
      int y_pos = (int)((thread_id - ((int)(thread_id / N))*N)/x_size)+1;
      if (y_pos % y_size != 0) //last values not reached
        y[thread_id] = x[thread_id+x_size] - x[thread_id];
      else
        y[thread_id] = x[thread_id-x_size] - x[thread_id];
    }

    if (dim == 3)
    {
      val = thread_id - size + N;
      if ((val) < 0)
        y[thread_id] = x[thread_id+N] - x[thread_id];
      else
        y[thread_id] = x[thread_id-N] - x[thread_id];
    }

    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Compute the backward difference between each value in a GPU vector (GPU function)
//! \brief considering a 3D data layout.
//!
//! y = bdiff(x)
//!
//! basically bdiff3(x) = -diff3_trans(x)
template <typename TType1>
__global__ void bdiff3_GPU(const unsigned dim, const unsigned x_size, const unsigned y_size, const TType1* x, TType1* y, unsigned size, bool borderWrap)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  int N = x_size * y_size;

  while (thread_id < size)
  {
    if (dim == 1)
    {
      int xPos = ((thread_id+1) % x_size);
      if (xPos == 0) // right border
        y[thread_id] = (borderWrap ? x[thread_id]: 0) - x[thread_id-1];
      else if (xPos == 1) // left border
        y[thread_id] = x[thread_id] - (borderWrap ? x[thread_id+x_size-1]: 0);
      else
        y[thread_id] = x[thread_id] - x[thread_id-1];
    }

    if (dim == 2)
    {
      int yPos = (int)((thread_id - ((int)(thread_id / N))*N)/x_size)+1;
      if (yPos % y_size == 0) // bottom border
        y[thread_id] = (borderWrap ? x[thread_id]: 0) - x[thread_id-x_size];
      else if (yPos % y_size == 1) //top border
        y[thread_id] = x[thread_id] - (borderWrap ? x[thread_id+(x_size*(y_size-1))]: 0);
      else
        y[thread_id] = x[thread_id] - x[thread_id-x_size];
    }

    if (dim == 3)
    {
      int zPos = (int)(thread_id/N);
      int zSize = size / N; 
      
      if (zSize == 1)
        y[thread_id] = 0;
      else if (zPos == 0) // first slice
        y[thread_id] = x[thread_id] - (borderWrap ? x[thread_id+size-N] : 0); 
      else if (zPos == (zSize-1)) // last slice
        y[thread_id] = (borderWrap ? x[thread_id] : 0) - x[thread_id-N];
      else
        y[thread_id] = x[thread_id] - x[thread_id-N]; 
    }

    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Compute the backward difference (sym. boundary cond.) between each value in a GPU vector (GPU function)
//! \brief considering a 3D data layout.
//!
//! y = bdiff(x)
//!
//! basically bdiff3(x) = -diff3_trans(x)
template <typename TType1>
__global__ void bdiff3sym_GPU(const unsigned dim, const unsigned x_size, const unsigned y_size, const TType1* x, TType1* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  int N = x_size * y_size;

  while (thread_id < size)
  {
    if (dim == 1)
    {
      int xPos = ((thread_id+1) % x_size);
      if (xPos == 0) // right border
        y[thread_id] = x[thread_id-1] - x[thread_id];
      else if (xPos == 1) // left border
        y[thread_id] = -2.0 * x[thread_id];
      else
        y[thread_id] = x[thread_id] - x[thread_id-1];
    }

    if (dim == 2)
    {
      int yPos = (int)((thread_id - ((int)(thread_id / N))*N)/x_size)+1;
      if (yPos % y_size == 0) // bottom border
        y[thread_id] = x[thread_id-x_size] - x[thread_id];
      else if (yPos % y_size == 1) //top border
        y[thread_id] = -2.0 * x[thread_id];
      else
        y[thread_id] = x[thread_id] - x[thread_id-x_size];
    }

    if (dim == 3)
    {
      int zPos = (int)(thread_id/N);
      int zSize = size / N;

      if (zSize == 1)
        y[thread_id] = 0;
      else if (zPos == 0) // first slice
        y[thread_id] = -2.0 * x[thread_id];
      else if (zPos == (zSize-1)) // last slice
        y[thread_id] = x[thread_id-N] - x[thread_id];
      else
        y[thread_id] = x[thread_id] - x[thread_id-N];
    }

    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Compute the difference between each value in a GPU vector (GPU function)
//! \brief last (x_size) value with first value
//!
//! y = difftrans(x)
template <typename TType1>
__global__ void difftrans_GPU(const unsigned dim, const unsigned x_size, const TType1* x, TType1* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    if (dim == 1)
    {
      if (((thread_id+1) % x_size ) != 0) //last values not reached
        y[thread_id+1] = x[thread_id] - x[thread_id+1];
      else
        y[thread_id+1-x_size] = x[thread_id] - x[thread_id+1-x_size];
    }

    if (dim == 2)
    {
      if (thread_id < size-x_size) //last values not reached
        y[thread_id+x_size] = x[thread_id] - x[thread_id+x_size];
      else
        y[thread_id-size+x_size] = x[thread_id] - x[thread_id-size+x_size];
    }

    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Compute the difference between each value in a GPU vector (GPU function)
//! \brief considering a 3D data layout.
//!
//! y = diff(x)
template <typename TType1>
__global__ void diff3trans_GPU(const unsigned dim, const unsigned x_size, const unsigned y_size, const TType1* x, TType1* y, unsigned size, bool borderWrap)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  int N = x_size * y_size;

  while (thread_id < size)
  {
    if (dim == 1)
    {
      int xPos = ((thread_id+1) % x_size);
      if (xPos == 0) // right border
        y[thread_id] = x[thread_id-1] - (borderWrap ? x[thread_id]: 0);
      else if (xPos == 1) // left border
        y[thread_id] = (borderWrap ? x[thread_id+x_size-1]: 0) - x[thread_id];
      else
        y[thread_id] = x[thread_id-1] - x[thread_id];
    }

    if (dim == 2)
    {
      int yPos = (int)((thread_id - ((int)(thread_id / N))*N)/x_size)+1;
      if (yPos % y_size == 0) // bottom border
        y[thread_id] = x[thread_id-x_size] - (borderWrap ? x[thread_id]: 0);
      else if (yPos % y_size == 1) //top border
        y[thread_id] = (borderWrap ? x[thread_id+(x_size*(y_size-1))]: 0) - x[thread_id];
      else
        y[thread_id] = x[thread_id-x_size]-x[thread_id];
    }

    if (dim == 3)
    {
      int zPos = (int)(thread_id/N);
      int zSize = size / N; 
      
      if (zSize == 1)
        y[thread_id] = 0;
      else if (zPos == 0) // first slice
        y[thread_id] = (borderWrap ? x[thread_id+size-N] : 0) - x[thread_id]; 
      else if (zPos == (zSize-1)) // last slice
        y[thread_id] = x[thread_id-N] - (borderWrap ? x[thread_id] : 0);
      else
        y[thread_id] = x[thread_id-N] - x[thread_id]; 
    }

    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Compute the backward difference between each value in a GPU vector (GPU function)
//! \brief considering a 3D data layout.
//!
//! y = diff(x)
template <typename TType1>
__global__ void bdiff3trans_GPU(const unsigned dim, const unsigned x_size, const unsigned y_size, const TType1* x, TType1* y, unsigned size, bool borderWrap)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    int val = thread_id-size+x_size;
    int N = x_size * y_size;
    if (dim == 1)
      if ( ((thread_id+1) % x_size ) != 0) //last values not reached
        y[thread_id] = x[thread_id] - x[thread_id+1];
      else
        y[thread_id] = borderWrap ? x[thread_id] - x[thread_id+1-x_size] : 0;

    if (dim == 2)
    {
      int y_pos = (int)((thread_id - ((int)(thread_id / N))*N)/x_size)+1;
      if (y_pos % y_size != 0) //last values not reached
        y[thread_id] = x[thread_id] - x[thread_id+x_size];
      else
        y[thread_id] = borderWrap ? x[thread_id] - x[thread_id-N+x_size] : 0;
    }

    if (dim == 3)
    {
      val = thread_id - size + N;
      if ((val) < 0) 
        y[thread_id] = x[thread_id] - x[thread_id+N];
      else
        y[thread_id] = borderWrap ? x[thread_id] - x[thread_id-size+N] : 0;
    }

    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief Compute the difference between each value in a GPU vector (GPU function)
//! \brief considering a 4D data layout
template <typename TType1>
__global__ void diff4_GPU(const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned z_size, const TType1* x, TType1* y, unsigned size, bool borderWrap)
{
  //blockDim.x number of threads in a block
  //blockIdx.x number of blocks in a grid
  //blockIdx.x index of threads in a block
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  while (thread_id < size)
  {
    int val = thread_id-size+x_size;
    int Nxy = x_size * y_size;
    int Nxyz = x_size * y_size * z_size;
    if (dim == 1)
      if ( ((thread_id+1) % x_size ) != 0) //last values not reached
        y[thread_id] = x[thread_id+1] - x[thread_id];
      else
        y[thread_id] = borderWrap ? x[thread_id+1-x_size] - x[thread_id] : 0;

    if (dim == 2)
    {
      int y_pos = (int)((thread_id - ((int)(thread_id / Nxy))*Nxy)/x_size)+1;
      if (y_pos % y_size != 0) //last values not reached
        y[thread_id] = x[thread_id+x_size] - x[thread_id];
      else
        y[thread_id] = borderWrap ? x[thread_id-Nxy+x_size] - x[thread_id] : 0;
    }

    if (dim == 3)
    {
      int z_pos = (int)((thread_id - ((int)(thread_id / Nxyz))*Nxyz)/Nxy)+1;
      if (z_pos % z_size != 0)
        y[thread_id] = x[thread_id+Nxy] - x[thread_id];
      else
        y[thread_id] = borderWrap ? x[thread_id-Nxyz+Nxy] - x[thread_id] : 0;
    }

    if (dim == 4)
    {
      val = thread_id - size + Nxyz;
      if ((val) < 0)
        y[thread_id] = x[thread_id+Nxyz] - x[thread_id];
      else
        y[thread_id] = borderWrap ? x[thread_id-size+Nxyz] - x[thread_id] : 0;
    }


    thread_id += blockDim.x*gridDim.x; // gives the number of threads in a grid
  }

}

//! \brief Compute the backward difference between each value in a GPU vector (GPU function)
//! \brief considering a 4D data layout.
//!
//! y = bdiff(x)
//!
//! basically bdiff4(x) = -diff4_trans(x)
//! 1, width, height, depth, data_gpu[0].data(), gradient[0].data(), N, false);
template <typename TType1>
__global__ void bdiff4_GPU(const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned z_size, const TType1* x, TType1* y, unsigned size, bool borderWrap)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;
  int Nxy = x_size * y_size;
  int Nxyz = x_size * y_size * z_size;

  if (dim == 1)
  {
    int xPos = ((thread_id+1) % x_size);
    if (xPos == 0) // right border
      y[thread_id] = (borderWrap ? x[thread_id]: 0) - x[thread_id-1];
    else if (xPos == 1) // left border
      y[thread_id] = x[thread_id] - (borderWrap ? x[thread_id+x_size-1]: 0);
    else
      y[thread_id] = x[thread_id] - x[thread_id-1];
  }
  if (dim == 2)
  {
    int yPos = (int)((thread_id - ((int)(thread_id / Nxy))*Nxy)/x_size)+1;
    if (yPos % y_size == 0) // bottom border
      y[thread_id] = (borderWrap ? x[thread_id]: 0) - x[thread_id-x_size];
    else if (yPos % y_size == 1) //top border
      y[thread_id] = x[thread_id] - (borderWrap ? x[thread_id+(x_size*(y_size-1))]: 0);
    else
      y[thread_id] = x[thread_id] - x[thread_id-x_size];
  }
  if (dim == 3)
  {
    int zPos = (int)((thread_id - ((int)(thread_id / Nxyz))*Nxyz)/Nxy)+1;
    if (zPos % z_size == 0) // last slice
      y[thread_id] = (borderWrap ? x[thread_id]: 0) - x[thread_id-Nxy];
    else if (zPos % z_size == 1) // first slice
      y[thread_id] = x[thread_id] - (borderWrap ? x[thread_id+(Nxy*(z_size-1))]: 0);
    else
      y[thread_id] = x[thread_id] - x[thread_id - Nxy];
  }
  if (dim == 4)
  {
    int tPos = (int)(thread_id/Nxyz);
    int tSize = size / Nxyz;

    if (tSize == 1) // only one dimesion no gradient calculation possible
      y[thread_id] = 0;
    else if (tPos == 0) // first timeframe
      y[thread_id] = x[thread_id] - (borderWrap ? x[thread_id + size - Nxyz] : 0);
    else if (tPos == (tSize-1)) // last timeframe
      y[thread_id] = (borderWrap ? x[thread_id] : 0) - x[thread_id - Nxyz];
    else
      y[thread_id] = x[thread_id] - x[thread_id - Nxyz];
  }
}

//! \brief generate max-value vector of two vector (elementwise).    (GPU function)
//! \brief z = abs(x1)>abs(x1 ? x1 : x2;
template <typename TTypeFloat, typename TType1, typename TType2, typename TType3>
__global__ void max_GPU(const TType1* x1,const TType2* x2, TType3* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  TTypeFloat abs1, abs2;

  while (thread_id < size)
  {
    abs1 = TTypeFloat(sqrt(TTypeFloat(agile::norm(x1[thread_id]))));
    abs2 = TTypeFloat(sqrt(TTypeFloat(agile::norm(x2[thread_id]))));
    y[thread_id] = (abs1 > abs2) ? x1[thread_id] : x2[thread_id];

    thread_id += blockDim.x*gridDim.x;
  }
}

//! \brief generate max-value vector of two vector (elementwise).    (GPU function)
//! \brief z = abs(x1)>abs(x1 ? x1 : x2;
template <typename TTypeFloat, typename TType1, typename TType2, typename TType3>
__global__ void max_GPU(const TType1* x1,const TType2 x2, TType3* y, unsigned size)
{
  unsigned thread_id = blockDim.x * blockIdx.x + threadIdx.x;

  TTypeFloat abs1, abs2;

  while (thread_id < size)
  {
    abs1 = TTypeFloat(sqrt(TTypeFloat(agile::norm(x1[thread_id]))));
    abs2 = TTypeFloat(sqrt(TTypeFloat(agile::norm(x2))));
    y[thread_id] = (abs1 > abs2) ? x1[thread_id] : x2;

    thread_id += blockDim.x*gridDim.x;
  }
}


// *****************************************************************************
// ************************* CPU functions from here on ************************
// *****************************************************************************

namespace agile
{

namespace lowlevel
{
  //! \brief Add a scaled GPU vector to another vector (host function).
  template <typename TType1, typename TType2, typename TType3>
  void addScaledVector(const TType1* x, const TType2& scale, const TType3* y,
                       typename promote<typename promote<TType1, TType2>::type,
                                        TType3>::type* z,
                       unsigned size)
  {
    unsigned grid_size
      = std::min((size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock(),65520u);
    addScaledVector_GPU<<<grid_size,
                          GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      *(const typename substitute_gpu_complex<TType2>::type*)(&scale),
      (const typename substitute_gpu_complex<TType3>::type*)y,
      (typename substitute_gpu_complex<
         typename promote<typename promote<TType1, TType2>::type,
                          TType3>::type>::type*)z,
      size);
  }

  //! \brief Add two vectors (host function).
  template <typename TType1, typename TType2>
  void addVector(const TType1* x, const TType2* y,
                 typename promote<TType1, TType2>::type* z, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    addVector_GPU<<<std::min(grid_size,65520u),
                    GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (const typename substitute_gpu_complex<TType2>::type*)y,
      (typename substitute_gpu_complex<
         typename promote<TType1, TType2>::type>::type*)z,
      size);
  }
  
  //! \brief Divide a scalar by a vector (elementwise; host function).
  template <typename TType1, typename TType2>
  void divideVector(const TType1& alpha, const TType2* x,
                    typename promote<TType1, TType2>::type* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    divideVector_GPU<<<std::min(grid_size,65520u),
                       GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      *(const typename substitute_gpu_complex<TType1>::type*)(&alpha),
      (const typename substitute_gpu_complex<TType2>::type*)x,
      (typename substitute_gpu_complex<
         typename promote<TType1, TType2>::type>::type*)y,
      size);
  }

  //! \brief Compute the bilinear form of two vectors (host function).
  template <typename TType1, typename TType2>
  typename promote<TType1, TType2>::type getBilinearForm(
    const TType1* x, const TType2* y, unsigned size)
  {
    typedef typename promote<TType1, TType2>::type promote_type;

    static GPUVector<promote_type> gpu_partial_result;
    static GPUHostVector<promote_type> host_partial_result;

    // resize the temporary storage for the partial results
    gpu_partial_result.resize(GPUEnvironment::getNumMultiprocessors());
    host_partial_result.resize(GPUEnvironment::getNumMultiprocessors());

    // get the scalar product partially
    getBilinearForm_GPU<<<GPUEnvironment::getNumMultiprocessors(),
                          GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (const typename substitute_gpu_complex<TType2>::type*)y,
      (typename substitute_gpu_complex<promote_type>::type*)
        gpu_partial_result.data(),
      size);
    // transfer the data to the host
    gpu_partial_result.copyToHost(host_partial_result);
    // do the last reduction on the CPU
    promote_type sum(0);
    for (unsigned counter = 0;
         counter < GPUEnvironment::getNumMultiprocessors(); ++counter)
      sum += host_partial_result[counter];
    return sum;
  }

  //! \brief Compute the scalar product of two vectors (host function).
  template <typename TType1, typename TType2>
  typename promote<TType1, TType2>::type getScalarProduct(
    const TType1* x, const TType2* y, unsigned size)
  {
    typedef typename promote<TType1, TType2>::type promote_type;

    static GPUVector<promote_type> gpu_partial_result;
    static GPUHostVector<promote_type> host_partial_result;

    // resize the temporary storage for the partial results
    gpu_partial_result.resize(GPUEnvironment::getNumMultiprocessors());
    host_partial_result.resize(GPUEnvironment::getNumMultiprocessors());

    // get the scalar product partially
    getScalarProduct_GPU<<<GPUEnvironment::getNumMultiprocessors(),
                           GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (const typename substitute_gpu_complex<TType2>::type*)y,
      (typename substitute_gpu_complex<promote_type>::type*)
        gpu_partial_result.data(),
      size);
    // transfer the data to the host
    gpu_partial_result.copyToHost(host_partial_result);
    // do the last reduction on the CPU
    promote_type sum(0);
    for (unsigned counter = 0;
         counter < GPUEnvironment::getNumMultiprocessors(); ++counter)
      sum += host_partial_result[counter];
    return sum;
  }

  //! Used to shift between centered, and not centered Fast Fourier Transform
  template <typename TType>
  void fftshift(TType* x, unsigned rows, unsigned cols)
  { 
    unsigned numRow_first;
    unsigned numCol_first;

    //calc the middle:
    if( (rows)%2 == 0)
      numRow_first = rows/2;    // if even
    else
      numRow_first = rows/2+1;  // if odd

    if(cols%2 == 0)
      numCol_first = cols/2;    // if even
    else
      numCol_first = cols/2+1;  // if odd

    unsigned numRow_second = rows-numRow_first;
    unsigned numCol_second = cols-numCol_first;

    //printf("\n fftshift begin");

    GPUVector<TType> s(rows*cols);             //gpu_partial_result

    unsigned A_offset_2 = numCol_first;
    unsigned A_offset_3 = cols * numRow_first;
    unsigned A_offset_4 = cols * numRow_first + numCol_first;

    unsigned S_offset_2 = numCol_second;
    unsigned S_offset_3 = cols * numRow_second;
    unsigned S_offset_4 = cols * numRow_second + numCol_second;
/*
        {                                       //INSIDE x
            std::vector<TType> iResHost2;
            iResHost2.resize(s.size());

            CUDA_SAFE_CALL(cudaMemcpy(
                &iResHost2[0],x,s.size() * sizeof(TType),cudaMemcpyDeviceToHost)
                );
            std::cout <<std::endl<<"inside x: ";
            for (unsigned i=0; i<iResHost2.size(); ++i)
              std::cout << iResHost2[i] << " ";
        }                                       //INSIDE x


        {                                       //INSIDE S
            std::vector<TType> iResHost;
            s.copyToHost(iResHost);
            std::cout <<"inside s: ";
            for (unsigned i=0; i<iResHost.size(); ++i)
              std::cout << iResHost[i] << " ";
        }                                       //INSIDE S
*/

   // printf("\n fftshift cudamemcopy2D");

    //A_1 -> S_4
    CUDA_SAFE_CALL(cudaMemcpy2D(
        s.data() + S_offset_4, cols * sizeof(TType),
        x, cols * sizeof(TType),
        numCol_first * sizeof(TType), numRow_first,
        cudaMemcpyDeviceToDevice)
        );

    //A_2 -> S_3
    CUDA_SAFE_CALL(cudaMemcpy2D(
        s.data() + S_offset_3, cols * sizeof(TType),
        x + A_offset_2, cols * sizeof(TType),
        numCol_second * sizeof(TType), numRow_first,
        cudaMemcpyDeviceToDevice)
        );

    //A_3 -> S_2
    CUDA_SAFE_CALL(cudaMemcpy2D(
        s.data() + S_offset_2, cols * sizeof(TType),
        x + A_offset_3, cols * sizeof(TType),
        numCol_first * sizeof(TType), numRow_second,
        cudaMemcpyDeviceToDevice)
        );

    //A_4 -> S_1
    CUDA_SAFE_CALL(cudaMemcpy2D(
        s.data(), cols * sizeof(TType),
        x + A_offset_4, cols * sizeof(TType),
        numCol_second * sizeof(TType), numRow_second,
        cudaMemcpyDeviceToDevice)
        );

    CUDA_SAFE_CALL(cudaMemcpy(
        x,s.data(),s.size() * sizeof(TType),cudaMemcpyDeviceToDevice)
        );

  //  printf("\n fftshift end");

//<- added HEIGL
  }


  //! inverse fftshift
  // needed for matrix with row != col
  template <typename TType>
  void ifftshift(TType* x, unsigned rows, unsigned cols)
  {
    unsigned numRow_second;
    unsigned numCol_second;

    //calc the middle:
    if(rows%2 == 0)
      numRow_second = rows/2;    // if even
    else
      numRow_second = rows/2+1;  // if odd

    if(cols%2 == 0)
      numCol_second = cols/2;    // if even
    else
      numCol_second = cols/2+1;  // if odd

    unsigned numRow_first = rows-numRow_second;
    unsigned numCol_first = cols-numCol_second;

    GPUVector<TType> s(rows*cols);             //gpu_partial_result

    unsigned A_offset_2 = numCol_first;
    unsigned A_offset_3 = cols * numRow_first;
    unsigned A_offset_4 = cols * numRow_first + numCol_first;

    unsigned S_offset_2 = numCol_second;
    unsigned S_offset_3 = cols * numRow_second;
    unsigned S_offset_4 = cols * numRow_second + numCol_second;

    //A_1 -> S_4
    CUDA_SAFE_CALL(cudaMemcpy2D(
        s.data() + S_offset_4, cols * sizeof(TType),
        x, cols * sizeof(TType),
        numCol_first * sizeof(TType), numRow_first,
        cudaMemcpyDeviceToDevice)
        );

    //A_2 -> S_3
    CUDA_SAFE_CALL(cudaMemcpy2D(
        s.data() + S_offset_3, cols * sizeof(TType),
        x + A_offset_2, cols * sizeof(TType),
        numCol_second * sizeof(TType), numRow_first,
        cudaMemcpyDeviceToDevice)
        );

    //A_3 -> S_2
    CUDA_SAFE_CALL(cudaMemcpy2D(
        s.data() + S_offset_2, cols * sizeof(TType),
        x + A_offset_3, cols * sizeof(TType),
        numCol_first * sizeof(TType), numRow_second,
        cudaMemcpyDeviceToDevice)
        );

    //A_4 -> S_1
    CUDA_SAFE_CALL(cudaMemcpy2D(
        s.data(), cols * sizeof(TType),
        x + A_offset_4, cols * sizeof(TType),
        numCol_second * sizeof(TType), numRow_second,
        cudaMemcpyDeviceToDevice)
        );

    CUDA_SAFE_CALL(cudaMemcpy(
        x,s.data(),s.size() * sizeof(TType),cudaMemcpyDeviceToDevice)
        );

  }


  //! \brief Extract the imaginary part of a GPU vector (host function).
  template <typename TType1>
  void imag(const TType1* x,
            typename to_real_type<TType1>::type* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    imag_GPU<<<std::min(grid_size,65520u),
               GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (typename to_real_type<TType1>::type*)y,
      size);
  }

  //! \brief Bilinear interpolation with reshape of 1d vector (host function).
  template <typename TType1, typename TType2>
  void interpolate2d(const TType1* src, unsigned numColumns,
                     unsigned numRows, bool reshapeRowMajor,
                     const std::complex<TType2>* pos, TType1* res,
                     unsigned size)
  {
    AGILE_TEXTURE_2D.filterMode = cudaFilterModeLinear;

    cudaArray *dev_source;
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<
      typename substitute_gpu_complex<TType1>::type >();

    CUDA_SAFE_CALL(cudaMallocArray(&dev_source, &channelDesc,
      numColumns, numRows));

    CUDA_SAFE_CALL(cudaMemcpy2DToArray(dev_source, 0, 0,
      (const typename substitute_gpu_complex<TType1>::type*)src,
      numColumns * sizeof(typename substitute_gpu_complex<TType1>::type),
      numColumns * sizeof(typename substitute_gpu_complex<TType1>::type),
      numRows, cudaMemcpyDeviceToDevice));

    // bind a 2d texture to the matrix as we need random access to it
    CUDA_SAFE_CALL(cudaBindTextureToArray(AGILE_TEXTURE_2D, dev_source));

    unsigned block_dim = GPUEnvironment::getMaxNumThreadsPerBlock();
    unsigned grid_dim = (size + block_dim - 1) / block_dim;

    if (reshapeRowMajor) {
      interpolate2DVectorRowMajor_GPU<<<grid_dim, block_dim>>>(
        (const typename substitute_gpu_complex<std::complex<TType2> >::type*)
        pos, (typename substitute_gpu_complex<TType1>::type*)res,
        size);
    } else {
      interpolate2DVectorColumnMajor_GPU<<<grid_dim, block_dim>>>(
        (const typename substitute_gpu_complex<std::complex<TType2> >::type*)
        pos, (typename substitute_gpu_complex<TType1>::type*)res,
        size);
    }

    // unbind texture and reset changes to global texture
    CUDA_SAFE_CALL(cudaUnbindTexture(AGILE_TEXTURE_2D));
    AGILE_TEXTURE_2D.filterMode = cudaFilterModePoint;
    CUDA_SAFE_CALL(cudaFreeArray(dev_source));
  }

  //! \brief Multiply a conjugate GPU vector with another one element-wise (host function).
  template <typename TType1, typename TType2>
  void multiplyConjElementwise(const TType1* x, const TType2* y,
                               typename promote<TType1, TType2>::type* z,
                               unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    multiplyConjElementwise_GPU<<<std::min(grid_size,65520u),
                                  GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (const typename substitute_gpu_complex<TType2>::type*)y,
      (typename substitute_gpu_complex<
         typename promote<TType1, TType2>::type>::type*)z,
      size);
  }

  //! \brief Multiply two GPU vectors element-wise (host function).
  template <typename TType1, typename TType2>
  void multiplyElementwise(const TType1* x, const TType2* y,
                           typename promote<TType1, TType2>::type* z,
                           unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    multiplyElementwise_GPU<<<std::min(grid_size,65520u),
                              GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (const typename substitute_gpu_complex<TType2>::type*)y,
      (typename substitute_gpu_complex<
         typename promote<TType1, TType2>::type>::type*)z,
      size);
  }

  //! \brief Divide two GPU vectors element-wise (host function).
  template <typename TType1, typename TType2>
  void divideElementwise(const TType1* x, const TType2* y,
                           typename promote<TType1, TType2>::type* z,
                           unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    divideElementwise_GPU<<<std::min(grid_size,65520u),
                              GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (const typename substitute_gpu_complex<TType2>::type*)y,
      (typename substitute_gpu_complex<
         typename promote<TType1, TType2>::type>::type*)z,
      size);
  }

  //! \brief Compute the l1-norm of a GPU vector (host function).
  template <typename TType1>
  typename to_real_type<TType1>::type norm1(const TType1* x, unsigned size)
  {
    typedef typename to_real_type<TType1>::type real_type;

    static GPUVector<real_type> gpu_partial_result;
    static GPUHostVector<real_type> host_partial_result;

    // resize the temporary storage for the partial results
    gpu_partial_result.resize(GPUEnvironment::getNumMultiprocessors());
    host_partial_result.resize(GPUEnvironment::getNumMultiprocessors());

    // get the sum of squares partially
    sumabs_GPU<<<GPUEnvironment::getNumMultiprocessors(),
                       GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (typename substitute_gpu_complex<real_type>::type*)
        gpu_partial_result.data(),
      size);
    // transfer the data to the host
    gpu_partial_result.copyToHost(host_partial_result);
    // do the last reduction on the CPU
    real_type sum(0);
    for (unsigned counter = 0;
         counter < GPUEnvironment::getNumMultiprocessors(); ++counter)
      sum += host_partial_result[counter];
    return sum;
  }

  //! \brief Compute the l2-norm of a GPU vector (host function).
  template <typename TType1>
  typename to_real_type<TType1>::type norm2(const TType1* x, unsigned size)
  {
    typedef typename to_real_type<TType1>::type real_type;

    static GPUVector<real_type> gpu_partial_result;
    static GPUHostVector<real_type> host_partial_result;

    // resize the temporary storage for the partial results
    gpu_partial_result.resize(GPUEnvironment::getNumMultiprocessors());
    host_partial_result.resize(GPUEnvironment::getNumMultiprocessors());

    // get the sum of squares partially
    sumofsquares_GPU<<<GPUEnvironment::getNumMultiprocessors(),
                       GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (typename substitute_gpu_complex<real_type>::type*)
        gpu_partial_result.data(),
      size);
    // transfer the data to the host
    gpu_partial_result.copyToHost(host_partial_result);
    // do the last reduction on the CPU
    real_type sum(0);
    for (unsigned counter = 0;
         counter < GPUEnvironment::getNumMultiprocessors(); ++counter)
      sum += host_partial_result[counter];
    // return the square root of the sum of squares
    return std::sqrt(sum);
  }

  //! \brief Extract the real part of a GPU vector (host function).
  template <typename TType1>
  void real(const TType1* x,
            typename to_real_type<TType1>::type* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    real_GPU<<<std::min(grid_size,65520u),
               GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (typename to_real_type<TType1>::type*)y,
      size);
  }

  //! \brief Multiply a GPU vector with a scalar (host function).
  template <typename TType1, typename TType2>
  void scale(const TType1& alpha, const TType2* x,
             typename promote<TType1, TType2>::type* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    scale_GPU<<<std::min(grid_size,65520u),
                GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      *(const typename substitute_gpu_complex<TType1>::type*)(&alpha),
      (const typename substitute_gpu_complex<TType2>::type*)x,
      (typename substitute_gpu_complex<
         typename promote<TType1, TType2>::type>::type*)y,
      size);
  }
  
  //! \brief Set every vector element to a constant (host function).
  template <typename TType>
  void setVectorConstant(const TType& value, TType* x, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    setVectorConstant_GPU<<<std::min(grid_size,65520u),
                            GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      *(const typename substitute_gpu_complex<TType>::type*)(&value),
      (typename substitute_gpu_complex<TType>::type*)x,
      size);
  }

  //! \brief Compute the square root of a GPU vector (host function).
  template <typename TType1>
  void sqrt(const TType1* x, TType1* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    //sqrt_GPU<TType1, typename to_float_type<TType1>::type>
    sqrt_GPU
      <<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
        (typename substitute_gpu_complex<TType1>::type*)x,
        (typename substitute_gpu_complex<TType1>::type*)y, size);
  }

  /*
  //! \brief compute average over last dimension (host function).
  template <typename TType1>
  void averVector(const TType1* x,
                 typename promote<TType1>::type* z, unsigned size, unsigned aver)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    averVector_GPU<<<std::min(grid_size,65520u),
                    GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (typename substitute_gpu_complex<
         typename promote<TType1>::type>::type*)z,
      size, aver);
  }
*/

  //! \brief Subtract two vectors (host function).
  template <typename TType1, typename TType2>
  void subVector(const TType1* x, const TType2* y,
                 typename promote<TType1, TType2>::type* z, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    subVector_GPU<<<std::min(grid_size,65520u),
                    GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      (const typename substitute_gpu_complex<TType2>::type*)y,
      (typename substitute_gpu_complex<
         typename promote<TType1, TType2>::type>::type*)z,
      size);
  }

  //! \brief Subtract a scaled GPU vector from another vector (host function).
  template <typename TType1, typename TType2, typename TType3>
  void subScaledVector(const TType1* x, const TType2& scale, const TType3* y,
                       typename promote<typename promote<TType1, TType2>::type,
                                        TType3>::type* z,
                       unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    subScaledVector_GPU<<<std::min(grid_size,65520u),
                          GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (const typename substitute_gpu_complex<TType1>::type*)x,
      *(const typename substitute_gpu_complex<TType2>::type*)(&scale),
      (const typename substitute_gpu_complex<TType3>::type*)y,
      (typename substitute_gpu_complex<
         typename promote<typename promote<TType1, TType2>::type,
                          TType3>::type>::type*)z,
      size);
  }


  //! \brief Conjugate GPU vector (host function).
  template <typename TType1>
  void conjVector(const TType1* x, TType1* z, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    conjVector_GPU<<<std::min(grid_size,65520u),
                                  GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (typename substitute_gpu_complex<TType1>::type*)x,
      (typename substitute_gpu_complex<TType1>::type*)z,
      size);

  }
  
  //! \brief Exponential of GPU vector (host function).
  template <typename TType1>
  void expVector(const TType1* x, TType1* z, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    exp_GPU<<<std::min(grid_size,65520u),
                                  GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
      (typename substitute_gpu_complex<TType1>::type*)x,
      (typename substitute_gpu_complex<TType1>::type*)z,
      size);

  }


  //! \brief returns abs-value in y (elementwise; host function).
  template <typename TType1>
  void absVector(const TType1* x,
                    typename to_real_type<TType1>::type* y, unsigned size)
  {
    unsigned grid_size = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
                         / GPUEnvironment::getMaxNumThreadsPerBlock();
    abs_GPU<typename to_float_type<TType1>::type>
      <<<std::min(grid_size,65520u),
                  GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
        (const typename substitute_gpu_complex<TType1>::type*)x,
        (typename to_real_type<TType1>::type*)y , size);
  }

  //! \brief returns phase-value in y (elementwise; host function).
  template <typename TType1, typename TType2>
  void phaseVector(const TType1* x,
                    TType2* y, unsigned size)
  {
    unsigned grid_size = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
                         / GPUEnvironment::getMaxNumThreadsPerBlock();
    phase_GPU<typename to_float_type<TType1>::type>
      <<<std::min(grid_size,65520u),
                  GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
        (const typename substitute_gpu_complex<TType1>::type*)x,
        (TType2*)y , size);
  }

  //! \brief Generate linearly spaced vector between a and b with n numbers
  //! x = linspace(a,b,n)               (host function)
  template <typename TType1>
  void linspace(TType1* x, unsigned size, float a, float b)
  {
    unsigned grid_size = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
                       / GPUEnvironment::getMaxNumThreadsPerBlock();
    linspace_GPU
      <<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>
        ( (typename to_float_type<TType1>::type*)x,
          a , b,
          size);
  }

  //! \brief Generates meshgrid Vectors for input Vector x and y
  //! meshgrid(x,xsize,y,ysize)               (host function)
  template <typename TType1>
  void meshgrid(TType1* mesh_x, TType1* mesh_y,
                const TType1* x, unsigned x_size, const TType1* y, unsigned y_size)
  {
    for(unsigned c=0;c<y_size;c++)
    {
        CUDA_SAFE_CALL(cudaMemcpy(
            mesh_x + x_size * c,
            x, x_size * sizeof(TType1),
            cudaMemcpyDeviceToDevice)
            );
    }

    unsigned grid_size = (x_size*y_size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
                       / GPUEnvironment::getMaxNumThreadsPerBlock();

    unsigned smSize=y_size*sizeof(TType1);
//    unsigned smSize=GPUEnvironment::getMaxNumThreadsPerBlock()*sizeof(TType1);

    meshgrid_GPU<<<std::min(grid_size,65520u), GPUEnvironment::getMaxNumThreadsPerBlock(), smSize>>>(
          (typename substitute_gpu_complex<TType1>::type*)mesh_y,
          (typename substitute_gpu_complex<TType1>::type*)y, x_size, y_size);

  }

  //! \brief Compute the power of alpha for every element of a GPU Vector (host function).
  template <typename TType1, typename TType2>
  void pow(const TType1& alpha,
           const TType2* x,
                 TType2* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    pow_GPU<<<std::min(grid_size,65520u),
         GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
                    *(const typename to_float_type<TType1>::type*)(&alpha),
                    (const typename to_float_type<TType2>::type*) x,
                    (typename to_float_type<TType2>::type*) y, size);
  }

  //! \brief generate pattern of a vector. (host function).
  template <typename TType1>
  void pattern(const TType1* x, typename to_real_type<TType1 >::type* z, unsigned size)
  {

    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    pattern_GPU<typename to_float_type<TType1>::type>
      <<<std::min(grid_size,65520u), GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
          (const typename substitute_gpu_complex<TType1>::type*)x,
          (typename to_real_type<TType1>::type*) z, size);

  }

  //! \brief Compute the difference between each value in a GPU vector  (host function).
  template <typename TType1>
  void diff(const unsigned dim, const unsigned x_size,
           const TType1* x,
                 TType1* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    diff_GPU<<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
                  dim, x_size,
                  (const typename substitute_gpu_complex<TType1>::type*) x,
                  (typename substitute_gpu_complex<TType1>::type*) y, size);
  }

  //! \brief Compute the difference between each value in a GPU vector transposed (host function).
  template <typename TType1>
  void difftrans(const unsigned dim, const unsigned x_size,
           const TType1* x,
                 TType1* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    difftrans_GPU<<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
                  dim, x_size,
                  (const typename substitute_gpu_complex<TType1>::type*) x,
                  (typename substitute_gpu_complex<TType1>::type*) y, size);
  }

  //! \brief Compute the difference between each value in a GPU Vector considering a 3d data layout (host function)
  template <typename TType1>
  void diff3(const unsigned dim, const unsigned x_size, const unsigned y_size,
           const TType1* x,
                 TType1* y, unsigned size, bool borderWrap)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    diff3_GPU<<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
                  dim, x_size, y_size,
                  (const typename substitute_gpu_complex<TType1>::type*) x,
                  (typename substitute_gpu_complex<TType1>::type*) y, size, borderWrap);
  }

  //! \brief Compute the difference (sym. boundary cond.) between each value in a GPU Vector considering a 3d data layout (host function)
  template <typename TType1>
  void diff3sym(const unsigned dim, const unsigned x_size, const unsigned y_size,
           const TType1* x,
                 TType1* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    diff3sym_GPU<<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
                  dim, x_size, y_size,
                  (const typename substitute_gpu_complex<TType1>::type*) x,
                  (typename substitute_gpu_complex<TType1>::type*) y, size);
  }


  //! \brief Compute the difference between each value in a GPU Vector considering a transposed 3d data layout (host function)
  template <typename TType1>
  void diff3trans(const unsigned dim, const unsigned x_size, const unsigned y_size,
           const TType1* x,
                 TType1* y, unsigned size, bool borderWrap)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    diff3trans_GPU<<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
                  dim, x_size, y_size,
                  (const typename substitute_gpu_complex<TType1>::type*) x,
                  (typename substitute_gpu_complex<TType1>::type*) y, size, borderWrap);
  }

  //! \brief Compute the backward difference between each value in a GPU Vector considering a 3d data layout (host function)
  template <typename TType1>
  void bdiff3(const unsigned dim, const unsigned x_size, const unsigned y_size,
           const TType1* x,
                 TType1* y, unsigned size, bool borderWrap)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    bdiff3_GPU<<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
                  dim, x_size, y_size,
                  (const typename substitute_gpu_complex<TType1>::type*) x,
                  (typename substitute_gpu_complex<TType1>::type*) y, size, borderWrap);
  }

  //! \brief Compute the backward difference between each value in a GPU Vector considering a 3d data layout (host function)
  template <typename TType1>
  void bdiff3sym(const unsigned dim, const unsigned x_size, const unsigned y_size,
                 const TType1* x,
                       TType1* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
            / GPUEnvironment::getMaxNumThreadsPerBlock();
    bdiff3sym_GPU<<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
                         dim, x_size, y_size,
                         (const typename substitute_gpu_complex<TType1>::type*) x,
                         (typename substitute_gpu_complex<TType1>::type*) y, size);
  }

  //! \brief Compute the backward difference between each value in a GPU Vector considering a transposed 3d data layout (host function)
  template <typename TType1>
  void bdiff3trans(const unsigned dim, const unsigned x_size, const unsigned y_size,
           const TType1* x,
                 TType1* y, unsigned size, bool borderWrap)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    bdiff3trans_GPU<<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
                  dim, x_size, y_size,
                  (const typename substitute_gpu_complex<TType1>::type*) x,
                  (typename substitute_gpu_complex<TType1>::type*) y, size, borderWrap);
  }

  //! \brief Compute the difference between each value in a GPU Vector considering a 4d data layout (host function)
  template <typename TType1>
  void diff4(const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned z_size,
                 const TType1* x, TType1* y, unsigned size, bool borderWrap)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    diff4_GPU<<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
                  dim, x_size, y_size, z_size,
                  (const typename substitute_gpu_complex<TType1>::type*) x,
                  (typename substitute_gpu_complex<TType1>::type*) y, size, borderWrap);
  }

  //! \brief Compute the backward difference between each value in a GPU Vector considering a 4d data layout (host function)
  template <typename TType1>
  void bdiff4(const unsigned dim, const unsigned x_size, const unsigned y_size, const unsigned z_size,
                  const TType1* x, TType1* y, unsigned size, bool borderWrap)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    bdiff4_GPU<<<std::min(grid_size,65520u),GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
                  dim, x_size, y_size, z_size,
                  (const typename substitute_gpu_complex<TType1>::type*) x,
                  (typename substitute_gpu_complex<TType1>::type*) y, size, borderWrap);
  }

  //! \brief generate max-value vector of two vector (elementwise) (host function).
  template <typename TType1, typename TType2>
  void max(const TType1* x1, const TType2* x2,
            typename promote<TType1, TType2>::type* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    max_GPU<typename to_float_type<TType1>::type>
      <<<std::min(grid_size,65520u), GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
          (const typename substitute_gpu_complex<TType1>::type*)x1,
          (const typename substitute_gpu_complex<TType2>::type*)x2,
          (typename substitute_gpu_complex<
             typename promote<TType1, TType2>::type>::type*) y, size);

  }

  //! \brief generate max-value vector of vector (elementwise) and scalar (host function).
  template <typename TType1, typename TType2>
  void max(const TType1* x1, const TType2& x2,
            typename promote<TType1, TType2>::type* y, unsigned size)
  {
    unsigned grid_size
      = (size + GPUEnvironment::getMaxNumThreadsPerBlock() - 1)
        / GPUEnvironment::getMaxNumThreadsPerBlock();
    max_GPU<typename to_float_type<TType1>::type>
      <<<std::min(grid_size,65520u), GPUEnvironment::getMaxNumThreadsPerBlock()>>>(
          (const typename substitute_gpu_complex<TType1>::type*)x1,
          *(const typename substitute_gpu_complex<TType2>::type*)&x2,
          (typename substitute_gpu_complex<
             typename promote<TType1, TType2>::type>::type*) y, size);
  }


  //! Generate data vector with value-padding in col-dimension
  template <typename TType>
  void expand_coldim(const TType* x_data, const TType* delta_o, const TType* delta_u,
                      unsigned rows, unsigned cols, unsigned col_o, unsigned col_u, TType* z)
  {
    //GPUVector<TType> delta_o(rows*col_o, 0);    //left zero vector
    //GPUVector<TType> delta_u(rows*col_u, 0);    //right zero vector


    //  |-----------|---------|-----------|
    //  |  delta_o  |    x    |  delta_u  |
    //  |-----------|---------|-----------|


    //delta_o
    CUDA_SAFE_CALL(cudaMemcpy2D(
        z, (col_o+col_u+cols)* sizeof(TType),
        delta_o, col_o * sizeof(TType),
        col_o * sizeof(TType), rows,
        cudaMemcpyDeviceToDevice)
        );

    //x_data
    CUDA_SAFE_CALL(cudaMemcpy2D(
        z+col_o, (col_o+col_u+cols)* sizeof(TType),
        x_data, cols * sizeof(TType),
        cols * sizeof(TType), rows,
        cudaMemcpyDeviceToDevice)
        );

    //delta_u
    CUDA_SAFE_CALL(cudaMemcpy2D(
        z+col_o+cols, (col_o+col_u+cols)* sizeof(TType),
        delta_u, col_u * sizeof(TType),
        col_u * sizeof(TType), rows,
        cudaMemcpyDeviceToDevice)
        );
  }


  //! Generate data vector with value-padding in row-dimension
  template <typename TType>
  void expand_rowdim(const TType* x_data, const TType* delta_o, const TType* delta_u,
                      unsigned rows, unsigned cols, unsigned row_o, unsigned row_u, TType* z)
  {
    //GPUVector<TType> delta_o(rows*row_o, 0);    //upper zero vector
    //GPUVector<TType> delta_u(rows*row_u, 0);    //bottom zero vector

    //  -------------
    //  |  delta_o  |
    //  -------------
    //  |     x     |
    //  -------------
    //  |  delta_u  |
    //  -------------


    //delta_o
    CUDA_SAFE_CALL(cudaMemcpy2D(
        z, cols * sizeof(TType),
        delta_o, cols * sizeof(TType),
        cols * sizeof(TType), row_o,
        cudaMemcpyDeviceToDevice)
        );

    //x_data
    CUDA_SAFE_CALL(cudaMemcpy2D(
        z+cols*row_o, cols * sizeof(TType),
        x_data, cols * sizeof(TType),
        cols * sizeof(TType), rows,
        cudaMemcpyDeviceToDevice)
        );

    //delta_u
    CUDA_SAFE_CALL(cudaMemcpy2D(
        z+cols*row_o+cols*rows, cols * sizeof(TType),
        delta_u, cols * sizeof(TType),
        cols * sizeof(TType), row_u,
        cudaMemcpyDeviceToDevice)
        );
  }

  //!copy data from a given input
  template <typename TType>
  void get_content(const TType* x_data, unsigned rows, unsigned cols,
                   unsigned row_offset, unsigned col_offset, TType* z, unsigned z_rows, unsigned z_cols)
  {

    //x_data -> z
    CUDA_SAFE_CALL(cudaMemcpy2D(
        z, z_cols * sizeof(TType),
        x_data + row_offset*cols +col_offset, cols * sizeof(TType),
        z_cols * sizeof(TType), z_rows,
        cudaMemcpyDeviceToDevice)
        );
  }
} // namespace lowlevel

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //Max Vector (implemented in CUBLAS)
  template <typename TType>
  struct max_Helper
  {
      static void max(unsigned int numElements, const TType* vec_in, int* max_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
        AGILE_ASSERT(false, ExceptionInvalidProgramFlow());
      }
  };

  template <>
  struct max_Helper<float>
  {
      static void max(unsigned int numElements, const float* vec_in, int* max_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasIsamax (handle, numElements, vec_in, unsigned(1), max_out));
      }
  };

  template <>
  struct max_Helper<double>
  {
      static void max(unsigned int numElements, const double* vec_in, int* max_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasIdamax(handle, numElements, vec_in, unsigned(1), max_out));
      }
  };

  template <>
  struct max_Helper<cuComplex>
  {
      static void max(unsigned int numElements, const cuComplex* vec_in, int* max_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasIcamax (handle, numElements, vec_in, unsigned(1), max_out));
      }
  };

  template <>
  struct max_Helper<cuDoubleComplex>
  {
      static void max(unsigned int numElements, const cuDoubleComplex* vec_in, int* max_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasIzamax(handle, numElements, vec_in, unsigned(1), max_out));
      }
  };
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //Copy Vector (implemented in CUBLAS)
  template <typename TType>
  struct copy_Helper
  {
      static void copy(unsigned int numElements, const TType* vec_in,
                        TType* vec_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
        AGILE_ASSERT(false, ExceptionInvalidProgramFlow());
      }
  };

  template <>
  struct copy_Helper<float>
  {
      static void copy(unsigned int numElements, const float* vec_in,
                        float* vec_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasScopy (handle, numElements, vec_in, unsigned(1), vec_out, unsigned(1)));
      }
  };

  template <>
  struct copy_Helper<double>
  {
      static void copy(unsigned int numElements, const double* vec_in,
                        double* vec_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasDcopy (handle, numElements, vec_in, unsigned(1), vec_out, unsigned(1)));
      }
  };

  template <>
  struct copy_Helper<cuComplex>
  {
      static void copy(unsigned int numElements, const cuComplex* vec_in,
                         cuComplex* vec_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasCcopy (handle, numElements, vec_in, unsigned(1), vec_out, unsigned(1)));
      }
  };

  template <>
  struct copy_Helper<cuDoubleComplex>
  {
      static void copy(unsigned int numElements, const cuDoubleComplex* vec_in,
                         cuDoubleComplex* vec_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasZcopy (handle, numElements, vec_in, unsigned(1), vec_out, unsigned(1)));
      }
  };

  //! \brief copy a GPU vector x to y (host function).
  template <typename TType>
  void copy(const GPUVector<TType>& x, GPUVector<TType>& y)
  {
      AGILE_ASSERT(x.size() == y.size(),
                    ExceptionSizeMismatch("x", "y",
                                          x.size(), y.size()));

      copy_Helper<typename substitute_gpu_complex<TType>::cublas_type>::copy(
                         x.size(), 
                         (const typename substitute_gpu_complex<TType>::cublas_type*) x.data(),
                         (typename substitute_gpu_complex<TType>::cublas_type*) y.data());
  };
  
  //! \brief derive maximum element of a GPU vector x (host function).
  template <typename TType>
  void maxElement(const GPUVector<TType>& x, int* maxVal)
  {
      max_Helper<typename substitute_gpu_complex<TType>::cublas_type>::max(
                         x.size(), 
                         (const typename substitute_gpu_complex<TType>::cublas_type*) x.data(),
                         maxVal);
  };
} // namespace agile

// End of $Id: gpu_vector.ipp 476 2011-06-16 08:54:14Z freiberger $.
