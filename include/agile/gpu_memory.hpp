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

// $Id: gpu_memory.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_MEMORY_HPP
#define AGILE_GPU_MEMORY_HPP

#include "agile/gpu_config.hpp"

namespace agile
{
  //! \brief Set GPU memory to a constant.
  //!
  //! This function sets all the elements \p data[0], ... \p data[\p size -1]
  //! to \p value.
  //! \param[in] data Pointer to the GPU to initialize.
  //! \param[in] size The amount of elements to initialize.
  //! \param[in] value The value for the elements.
  //! NOTE: This function will be removed as it is redundant with
  //! \p setVectorConstant
  template <typename TType>
  void initializeGPUMemory(TType* data, unsigned size, const TType& value);

  //! \brief Add the columns of a rectangular memory area.
  //!
  //! This function takes a pointer to a memory area of size \p num_rows *
  //! \p num_columns. It adds the elements along the rows and stores the
  //! sum in the memory area given by \p output. Thus, \p output has to be
  //! of size \p num_rows. The input is destroyed upon calling this function
  //! because it is needed to store partial sums.
  //! \param[in] input Pointer to memory to be summed (elements will be
  //! destroyed).
  //! \param[in] num_rows Number of rows of input.
  //! \param[in] num_columns Number of columns of input.
  //! \param[out] output Pointer to output memory of size \p num_rows. Will
  //! contain the sum of the input columns.
  template <typename TType>
  void addMemoryColumns(TType* input, unsigned num_rows, unsigned num_columns,
                        TType* output);

} // namespace agile

#endif // AGILE_GPU_MEMORY_HPP

// End of $Id: gpu_memory.hpp 476 2011-06-16 08:54:14Z freiberger $.
