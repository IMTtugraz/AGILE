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

// $Id: dense_vector.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_IO_DENSE_VECTOR_HPP
#define AGILE_IO_DENSE_VECTOR_HPP

#include "agile/gpu_config.hpp"
#include "agile/exception.hpp"
#include "agile/io/io_utils.hpp"

namespace agile
{
  namespace detail
  {
    //! \brief Low-level function for writing a dense vector.
    template <typename TType>
    void writeDenseVector(const std::string& file_name, unsigned size,
                          const TType* data, const std::string& comment = "")
    {
      std::ofstream file(file_name.c_str());
      AGILE_ASSERT(file.is_open(),
                    StandardException::ExceptionFileIO(file_name));

      std::vector<unsigned> temp_size(1, size);
      writeHeader(file, "vector dense " + type_to_file_string<TType>::get(),
                  temp_size, comment);
      file.write((const char*)data, size * sizeof(TType));
    }

  } // namespace detail

  //! \brief Write a dense GPU vector to a file.
  //!
  //! \param[in] file_name The name of the ouput file.
  //! \param[in] vector The vector to be written to the file.
  //! \param[in] comment An additional comment to be written.
  template <typename TType>
  void writeToFile(const std::string& file_name, const GPUVector<TType>& vector,
                   const std::string& comment = "")
  {
    // copy the vector data to the host
    std::vector<TType> host_data;
    vector.copyToHost(host_data);
    detail::writeDenseVector(file_name, vector.size(), &host_data[0], comment);
  }

  //! \brief Write a dense CPU vector to a file.
  //!
  //! \param[in] file_name The name of the ouput file.
  //! \param[in] vector The vector to be written to the file.
  //! \param[in] comment An additional comment to be written.
  template <typename TType>
  void writeToFile(const std::string& file_name,
                   const std::vector<TType>& vector,
                   const std::string& comment = "")
  {
    detail::writeDenseVector(file_name, vector.size(), &vector[0], comment);
  }

  //! \brief Read a dense GPU vector from a file.
  template <typename TType>
  void readFromFile(const std::string& file_name, GPUVector<TType>&  vector)
  {
    std::vector<TType> host_data;
    readFromFile(file_name, host_data);
    vector.assignFromHost(host_data.begin(), host_data.end());
  }

  //! \brief Read a dense CPU vector from a file.
  //!
  //! The output vector will be resized such that it is large enough to hold
  //! all the bytes from the file. However, it does not make sense to load
  //! a 32-bit integer file into a 32-bit floating point vector, for example.
  //! Thus, it is your responsibility to make sure that the types match.
  //! \param[out] vector The output vector.
  template <typename TType>
  void readFromFile(const std::string& file_name, std::vector<TType>& vector)
  {
    std::ifstream file(file_name.c_str());
    AGILE_ASSERT(file.is_open(),
                  StandardException::ExceptionFileIO(file_name));

    // read the header and make sure, the type is "vector dense <something>"
    FileHeader header(file);
    header.checkType("vector dense");
    // compute how large the vector has to be such that all elements fit in
    // there
    size_t byte_size = header.size[0] * header.type_size;
    vector.resize((byte_size + sizeof(TType) - 1) / sizeof(TType));
    file.read((char*)vector[0], byte_size);
  }

} // namespace agile

#endif // AGILE_IO_DENSE_VECTOR_HPP

// End of $Id: dense_vector.hpp 476 2011-06-16 08:54:14Z freiberger $.
