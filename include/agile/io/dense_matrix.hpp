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

// $Id: dense_matrix.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_IO_DENSE_MATRIX_HPP
#define AGILE_IO_DENSE_MATRIX_HPP

#include "agile/gpu_config.hpp"
#include "agile/exception.hpp"
#include "agile/io/io_utils.hpp"

namespace agile
{
  namespace detail
  {
    //! \brief Low-level function for writing a dense matrix.
    template <typename TType>
    void writeDenseMatrix(const std::string& file_name,
                          unsigned num_rows, unsigned num_columns,
                          const TType* data, const std::string& comment = "")
    {
      std::ofstream file(file_name.c_str());
      AGILE_ASSERT(file.is_open(),
                    StandardException::ExceptionFileIO(file_name));

      std::vector<unsigned> temp_size(2);
      temp_size[0] = num_rows;
      temp_size[1] = num_columns;
      writeHeader(file, "matrix dense " + type_to_file_string<TType>::get(),
                  temp_size, comment);
      file.write((const char*)data, num_rows * num_columns * sizeof(TType));
    }

  } // namespace detail

  //! \brief Write a dense GPU matrix to a file.
  //!
  //! \param[in] file_name The name of the ouput file.
  //! \param[in] matrix The matrix to be written to the file.
  //! \param[in] comment An additional comment to be written.
  template <typename TType>
  void writeToFile(const std::string& file_name, const GPUMatrixPitched<TType>& matrix,
                   const std::string& comment = "")
  {
    // copy the matrix data to the host
    std::vector<TType> host_data;
    matrix.copyToHost(host_data);
    detail::writeDenseMatrix(file_name, matrix.getNumRows(),
                             matrix.getNumColumns(), &host_data[0], comment);
  }

  //! \brief Write a dense CPU matrix to a file.
  //!
  //! \param[in] file_name The name of the ouput file.
  //! \param[in] matrix The matrix to be written to the file.
  //! \param[in] comment An additional comment to be written.
  template <typename TType>
  void writeToFile(const std::string& file_name, const CPUMatrix<TType>& matrix,
                   const std::string& comment = "")
  {
    detail::writeDenseMatrix(file_name, matrix.getNumRows(),
                             matrix.getNumColumns(), matrix.data(),
                             comment);
  }

} // namespace agile

#endif // AGILE_IO_DENSE_MATRIX_HPP

// End of $Id: dense_matrix.hpp 476 2011-06-16 08:54:14Z freiberger $.
