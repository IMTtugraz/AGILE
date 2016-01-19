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

// $Id: example.cu 452 2011-05-31 12:00:18Z freiberger $

#include <iostream>
#include <vector>
#include <cuda.h>

//! \brief Set all elements of vector x to the value a.
//!
//! Each thread sets one element to the given value which means we need at
//! least one thread per vector element. Therefore, it is required that
//! threads_per_block * grid_size >= n.
//! \param[in] x Vector to be set.
//! \param[in] a Value for each element.
//! \param[in] n Size of the vector.
__global__ void vector_set(float *x, float a, int n)
{
  // calculate the vector index to be processed by this thread
  int index = blockDim.x * blockIdx.x + threadIdx.x;
  // set the element to the value given
  if (index < n)
    x[index] = a;
}

//! \brief Add to vectors.
//!
//! The i-th thread calculates the sum x[i] + y[i] and stores the result in
//! z[i]. Thus, we need at least \p n threads which means that
//! threads_per_block * grid_size >= n.
//! \param[in] x First vector.
//! \param[in] y Second vector.
//! \param[out] z Sum of \p x and \p y.
//! \param[in] n Size of the vectors.
__global__ void vector_add(float *x, float *y, float *z, int n)
{
  int index = blockDim.x * blockIdx.x + threadIdx.x;
  // each thread sums two elements
  if (index < n)
    z[index] = x[index] + y[index];
}

//! \brief Calculate the scalar product of two vectors.
//!
//! The scalar product is calculated only partially on the GPU. In the first
//! step the i-th thread sums the products of the elements i, i + S, i + 2*S...,
//! where the stride S is equal to threads_per_block * grid_size.
//! In the second step each thread block reduces its partial sums in a parallel
//! manner and stores the result to the vector \p w.
//! \param[in] x First vector.
//! \param[in] y Second vector.
//! \param[out] w Partial scalar product which has to be reduced by the CPU.
//!               The vector has to be of length \p grid_size.
//! \param[in] n Size of the vectors.
__global__ void scalar_product(float *x, float *y, float *w, int n)
{
  __shared__ float partial_sum[512];

  // sum the products of the vector elements with a stride
  // threads_per_block * grid_size
  float s = 0.0f;
  for (int j = blockDim.x*blockIdx.x + threadIdx.x; j < n;
       j += blockDim.x * gridDim.x)
    s += x[j] * y[j];
  // store this partial sum within this thread block
  partial_sum[threadIdx.x] = s;

  // Reduce the partial sums within this thread block by parallel sums. In the
  // first iteration the threads 0 .. blockDim.x/2-1 add the contribution
  // of the threads blockDim.x/2 .. blockDim.x-1. In the second iteration
  // the threads 0 .. blockDim.x/4-1 add the value of
  // blockDim.x/4..blockDim.x/2 and so on until in the last iteration
  // thread 0 adds the value of thread 1.
  for (int i = blockDim.x >> 1; i > 0; i >>= 1)
  {
    __syncthreads();
    if (threadIdx.x < i)
      partial_sum[threadIdx.x] += partial_sum[threadIdx.x + i];
  }
  // thread 0 writes the output to the vector w
  if (threadIdx.x == 0)
    w[blockIdx.x] = partial_sum[0];
}

//! Vector maximum norm.
//!
//! Also the maximum norm is calculated partially only. The algorithm is very
//! similar to the scalar product.
//! \param[in] x Vector to apply the maximum norm.
//! \param[out] w Vector of partial maximum elements which has to be reduced by
//!               the CPU. This vector has to be of length \p grid_size.
//! \param[in] n Size of the vector.
__global__ void vector_max_norm(float *x, float *w, int n)
{
  // vector to store the maximum elements within this thread block
  __shared__ float max_elements[512];

  // look for maximum elements with stride threads_per_block * grid_size
  float s = 0.0f;
  for (int j = blockDim.x * blockIdx.x + threadIdx.x; j < n;
       j += blockDim.x * gridDim.x)
  {
    float t = fabs(x[j]);
    if (s < t)
      s = t;
  }
  max_elements[threadIdx.x] = s;

  for(int i = blockDim.x >> 1; i > 0; i >>= 1)
  {
    __syncthreads();
    if (threadIdx.x < i)
    {
      if (max_elements[threadIdx.x] < max_elements[threadIdx.x + i])
        max_elements[threadIdx.x] = max_elements[threadIdx.x + i];
    }
  }
  if (threadIdx.x == 0)
    w[blockIdx.x] = max_elements[0];
}

int main (int argc, char** argv)
{
  // get the amount of GPUs
  int num_gpus = 0;
  cudaGetDeviceCount(&num_gpus);
  std::cout << "Found " << num_gpus << " GPU(s)" << std::endl;

  // assign a device to the current thread
  int device = 0;
  cudaSetDevice(device);

  // get the device properties and print them to the terminal
  cudaDeviceProp properties;
  cudaGetDeviceProperties(&properties, device);
  std::cout << "Device name:                  " << properties.name << std::endl;
  std::cout << "Compute capability:           " << properties.major
            << "." << properties.minor << std::endl;
  std::cout << "Clock frequency (MHz):        " << properties.clockRate / 1000
            << std::endl;
  std::cout << "32-bit registers per block:   " << properties.regsPerBlock
            << std::endl;
  std::cout << "Total global memory (MB):     "
            << properties.totalGlobalMem / (1024 * 1024) << std::endl;
  std::cout << "Shared memory per block (kB): "
            << properties.sharedMemPerBlock / 1024 << std::endl;
  std::cout << "Total const memory (kB):      "
            << properties.totalConstMem / 1024 << std::endl;
  std::cout << "Number of multiprocessors:    "
            << properties.multiProcessorCount << std::endl;
  std::cout << "Max threads per block:        "
            << properties.maxThreadsPerBlock << std::endl;
  std::cout << std::endl;

  // choose a vector size
  int n = 1024 * 1024;
  // set the amount of threads per block; this choice is somethat arbitrary
  const int THREADS_PER_BLOCK = 256;
  // however, the grid dimension has to be set such that every vector gets
  // one corresponding thread
  const int GRID_SIZE = (n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

  float* dev_x;
  float* dev_y;
  float* dev_z;
  // allocate memory on the GPU
  cudaMalloc((void**)&dev_x, n * sizeof(float));
  cudaMalloc((void**)&dev_y, n * sizeof(float));
  cudaMalloc((void**)&dev_z, n * sizeof(float));

  // initialize the vectors
  std::cout << "Initializing vectors" << std::endl;
  std::cout << "  x[i] = 2" << std::endl;
  vector_set<<<GRID_SIZE, THREADS_PER_BLOCK>>>(dev_x, 2.0, n);
  std::cout << "  y[i] = 5" << std::endl;
  vector_set<<<GRID_SIZE, THREADS_PER_BLOCK>>>(dev_y, 5.0, n);

  // add two vectors
  std::cout << "Testing addition z <= x + y" << std::endl;
  vector_add<<<GRID_SIZE, THREADS_PER_BLOCK>>>(dev_x, dev_y, dev_z, n);
  // copy the vector z back to the CPU
  std::vector<float> host_z(n, 0);
  cudaMemcpy(&host_z[0], dev_z, n * sizeof(float), cudaMemcpyDeviceToHost);
  // print the first element of z
  std::cout << "  z[0] = " << host_z[0] << std::endl;

  // each multiprocessor calculates one partial scalar product
  std::cout << "Scalar product" << std::endl;
  float* dev_w;
  cudaMalloc((void**)&dev_w, properties.multiProcessorCount * sizeof(float));
  scalar_product<<<properties.multiProcessorCount, THREADS_PER_BLOCK>>>(
    dev_x, dev_y, dev_w, n);
  // copy the result back to the CPU
  std::vector<float> host_w(properties.multiProcessorCount, 0);
  cudaMemcpy(&host_w[0], dev_w, properties.multiProcessorCount * sizeof(float),
             cudaMemcpyDeviceToHost);
  // reduce on the CPU
  float result = 0;
  for (int i = 0; i < host_w.size(); ++i)
    result += host_w[i];
  std::cout << "  (x, y) = " << result << std::endl;

  // test the maximum norm
  std::cout << "Maximum norm" << std::endl;
  vector_max_norm<<<properties.multiProcessorCount, THREADS_PER_BLOCK>>>(
    dev_x, dev_w, n);
  cudaMemcpy(&host_w[0], dev_w, properties.multiProcessorCount * sizeof(float),
             cudaMemcpyDeviceToHost);
  result = 0;
  for (int i = 0; i < host_w.size(); ++i)
    if (host_w[i] > result)
      result = host_w[i];
  std::cout << "  ||x||_oo = " << result << std::endl;
  vector_max_norm<<<properties.multiProcessorCount, THREADS_PER_BLOCK>>>(
    dev_y, dev_w, n);
  cudaMemcpy(&host_w[0], dev_w, properties.multiProcessorCount * sizeof(float),
             cudaMemcpyDeviceToHost);
  result = 0;
  for (int i = 0; i < host_w.size(); ++i)
    if (host_w[i] > result)
      result = host_w[i];
  std::cout << "  ||y||_oo = " << result << std::endl;
  vector_max_norm<<<properties.multiProcessorCount, THREADS_PER_BLOCK>>>(
    dev_z, dev_w, n);
  cudaMemcpy(&host_w[0], dev_w, properties.multiProcessorCount * sizeof(float),
             cudaMemcpyDeviceToHost);
  result = 0;
  for (int i = 0; i < host_w.size(); ++i)
    if (host_w[i] > result)
      result = host_w[i];
  std::cout << "  ||z||_oo = " << result << std::endl;

  // free what was allocated before
  cudaFree(dev_x);
  cudaFree(dev_y);
  cudaFree(dev_z);
  cudaFree(dev_w);

  return 0;
}

// End of $Id: example.cu 452 2011-05-31 12:00:18Z freiberger $.
