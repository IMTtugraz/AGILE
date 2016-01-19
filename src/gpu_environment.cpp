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

// $Id: gpu_environment.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/gpu_environment.hpp"

namespace agile
{
  bool GPUEnvironment::m_initialized = false;
  unsigned GPUEnvironment::m_device_id;
  cudaDeviceProp GPUEnvironment::m_property;
  unsigned GPUEnvironment::m_max_threads_per_block_dim2;

  cublasHandle_t GPUEnvironment::cublas_handle_;

  void GPUEnvironment::printInformation(std::ostream& stream)
  {
    if (!m_initialized)
    {
      stream << "No device allocated, yet." << std::endl;
       return;
    }

    stream << "Device name:                  " << m_property.name << std::endl;
    stream << "Compute capability:           " << m_property.major
           << "." << m_property.minor << std::endl;
    stream << "Clock frequency (MHz):        " << m_property.clockRate / 1000
           << std::endl;
    stream << "32-bit registers per block:   " << m_property.regsPerBlock
           << std::endl;
    stream << "Total global memory (MB):     "
           << m_property.totalGlobalMem / (1024 * 1024) << std::endl;
    stream << "Shared memory per block (kB): "
           << m_property.sharedMemPerBlock / 1024 << std::endl;
    stream << "Total const memory (kB):      "
           << m_property.totalConstMem / 1024 << std::endl;
    stream << "Number of multiprocessors:    "
           << m_property.multiProcessorCount << std::endl;
    stream << "Max threads per block:        "
           << m_property.maxThreadsPerBlock << std::endl;
  }
  
  void GPUEnvironment::printUsage(std::ostream& stream)
  {
    if (!m_initialized)
    {
      stream << "No device allocated, yet." << std::endl;
       return;
    }
    size_t free_mem = 0;
    size_t total_mem = 0;
    CUDA_SAFE_CALL(cudaMemGetInfo(&free_mem, &total_mem));

    stream << "Device name:                   " << m_property.name << std::endl;
    stream << "Memory usage (used/total) (MB):" << (total_mem-free_mem)/(1024*1024) 
      << "/" << total_mem/(1024*1024)<< std::endl;
  }

} // namespace agile

// End of $Id: gpu_environment.cpp 476 2011-06-16 08:54:14Z freiberger $.
