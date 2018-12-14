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

// $Id: gpu_timer.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_TIMER_HPP
#define AGILE_GPU_TIMER_HPP

#include <iostream>
#include <iomanip>

#include "agile/gpu_config.hpp"
#include "cuda.h"
#include "cuda_runtime.h"

namespace agile
{
  //! \brief A simple GPU timer class.
  class GPUTimer
  {
    public:

      //! \brief Construct the timer.
      GPUTimer()
      {
        CUDA_SAFE_CALL( cudaEventCreate(&m_start) );
        CUDA_SAFE_CALL( cudaEventCreate(&m_end) );
        start();
      }

      //! \brief Destructor.
      virtual ~GPUTimer()
      {
        CUDA_SAFE_CALL( cudaEventDestroy(m_start) );
        CUDA_SAFE_CALL( cudaEventDestroy(m_end) );
      }

      //! \brief Start/Reset the timer.
      void start()
      {
        CUDA_SAFE_CALL( cudaEventRecord(m_start,0) );
      }

      //! \brief Get the elapsed time in seconds
      //! \return The elapsed time between start and stop in milliseconds.
      float stop()
      {
        float elapsed_time;
        CUDA_SAFE_CALL( cudaEventRecord(m_end, 0) );
        CUDA_SAFE_CALL( cudaEventSynchronize(m_end) );
        CUDA_SAFE_CALL( cudaEventElapsedTime(&elapsed_time, m_start, m_end) );
        return elapsed_time;
      }

    private:

      //! \brief Start event.
      cudaEvent_t m_start;

      //! \brief Stop event.
      cudaEvent_t m_end;
  };

} // namespace agile

#endif // AGILE_GPU_TIMER_HPP

// End of $Id: gpu_timer.hpp 476 2011-06-16 08:54:14Z freiberger $.
