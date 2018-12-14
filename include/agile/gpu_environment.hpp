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

// $Id: gpu_environment.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_ENVIRONMENT_HPP
#define AGILE_GPU_ENVIRONMENT_HPP

#include "agile/gpu_config.hpp"
#include "agile/exception.hpp"
#include <ostream>
#include <cuda_runtime.h>
#include <cublas_v2.h>

namespace agile
{
  class GPUEnvironment
  {
    public:
      //! \brief Get the number of GPUs in this PC.
      //!
      //! This method determines the number of GPUs in this PC and returns it.
      //! \return The number of GPUs.
      static unsigned getNumGPUs()
      {
        unsigned num_gpus;
        CUDA_SAFE_CALL(cudaGetDeviceCount((int*)&num_gpus));
        return num_gpus;
      }

      //! \brief Allocate a GPU.
      //!
      //! This method assigns a GPU to the current process. If the current
      //! process already has a device allocated, nothing is done.
      //! \param[in] device_id The id of the GPU to allocate This number can
      //! be in the range from 0 to \p getNumGPUs() - 1.
      static void allocateGPU(unsigned device_id)
      {
        // if a device is already allocated, we must not allocate it again
        if (m_initialized)
          return;

        CUDA_SAFE_CALL(cudaSetDevice(device_id));
        CUDA_SAFE_CALL(cudaGetDeviceProperties(&m_property, device_id));
        m_device_id = device_id;
        m_initialized = true;
        // constrain the maximum amount of threads
        if (m_property.maxThreadsPerBlock > int(MAX_NUM_THREADS_PER_BLOCK))
          m_property.maxThreadsPerBlock = MAX_NUM_THREADS_PER_BLOCK;
        m_max_threads_per_block_dim2 = 1;
        while (m_max_threads_per_block_dim2 < MAX_NUM_THREADS_PER_BLOCK_DIM2
               && m_max_threads_per_block_dim2 * m_max_threads_per_block_dim2
                  <= unsigned(m_property.maxThreadsPerBlock))
          m_max_threads_per_block_dim2 *= 2;
        m_max_threads_per_block_dim2 /= 2;

        CUBLAS_SAFE_CALL(cublasCreate (&cublas_handle_ )) ;
      }

      //! \brief Return the id of the allocated device.
      //!
      //! This method returns the id of the allocated GPU. This is the
      //! same number as was supplied to \p allocateGPU for the allocation.
      //! \return The id of the allocated GPU.
      static unsigned getDeviceID()
      {
        AGILE_ASSERT(m_initialized, StandardException::ExceptionMessage(
                                       "GPUEnvironment is not initialized"));
        return m_device_id;
      }

      //! \brief Determine if a GPU has been allocated.
      //!
      //! This method returns true if a GPU has already been initialized from
      //! this process.
      //! \return True, if a GPU has been allocated.
      static bool initialized()
      {
        return m_initialized;
      }

      //! \brief Get the maximum amount of threads per block.
      //!
      //! This function returns the maximum amount of threads in a single
      //! block. Note that this value is constrained by the constant
      //! \p MAX_NUM_THREADS_PER_BLOCK defined in \p gpu_config.hpp. Thus this
      //! value is actually the minimum of \p MAX_NUM_THREADS_PER_BLOCK and
      //! the maximum value specified by the graphics card.
      //! \return The maximum number of threads in a single block the GPU
      //! allows or \p MAX_NUM_THREADS_PER_BLOCK if this value should be
      //! smaller.
      inline
      static unsigned getMaxNumThreadsPerBlock()
      {
        AGILE_ASSERT(m_initialized, StandardException::ExceptionMessage(
                                       "GPUEnvironment is not initialized"));
        // the constraining of the value is done in allocateGPU
        return m_property.maxThreadsPerBlock;
      }

      //! \brief Get the maximum amount of threads per block for 2-dimensional grids.
      //!
      //! This function returns the maximum number of threads for a single
      //! dimension for 2-dimensional grids. It is guaranteed that
      //! the value returned is not larger than \p MAX_NUM_THREADS_PER_BLOCK_DIM2
      //! and that the square of the value returned is not larger than
      //! \p getMaxNumThreadsPerBlock().
      //! \return The maximum number of blocks fro 2-dimensional grids.
      inline
      static unsigned getMaxNumThreadsPerBlockDim2()
      {
        AGILE_ASSERT(m_initialized, StandardException::ExceptionMessage(
                                       "GPUEnvironment is not initialized"));
        // the constraining of the value is done in allocateGPU
        return m_max_threads_per_block_dim2;
      }

      //! \brief Get the number of multiprocessors.
      //!
      //! \return The number of multiprocessors.
      inline
      static unsigned getNumMultiprocessors()
      {
        AGILE_ASSERT(m_initialized, StandardException::ExceptionMessage(
                                       "GPUEnvironment is not initialized"));
        return m_property.multiProcessorCount;
      }

      //! \brief Output some information.
      //!
      //! This method prints some information about the allocated GPU to
      //! the given stream.
      //! \param stream The stream to print the information to.
      static void printInformation(std::ostream& stream);

      //! \brief Output memory usage information.
      //!
      //! This method prints the usage (available memory/total memory)
      //! information about the allocated GPU to
      //! the given stream.
      //! \param stream The stream to print the information to.
      static void printUsage(std::ostream& stream);

      //! Destructor to destroy the cublasHandle
      virtual ~GPUEnvironment()
      {
        CUBLAS_SAFE_CALL(cublasDestroy (cublas_handle_ )) ;
      }

      //! return the cublas Handle
      inline
      static cublasHandle_t getcublasHandle()
      {
        AGILE_ASSERT(cublas_handle_, StandardException::ExceptionMessage(
                                       "GPUEnvironment is not initialized"));
        return cublas_handle_;
      }


    private:
      //! cublas handle
      static cublasHandle_t cublas_handle_;

      //! Flag to store if a GPU is allocated.
      static bool m_initialized;

      //! The id of the allocated GPU.
      static unsigned m_device_id;

      //! The properties of the allocated GPU.
      static cudaDeviceProp m_property;

      //! Maximum amount of blocks for 2-dimensional grids.
      static unsigned m_max_threads_per_block_dim2;
  };

} // namespace agile

#endif // AGILE_GPU_ENVIRONMENT_HPP

// End of $Id: gpu_environment.hpp 476 2011-06-16 08:54:14Z freiberger $.
