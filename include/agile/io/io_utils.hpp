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

// $Id: io_utils.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_IO_IO_UTILS_HPP
#define AGILE_IO_IO_UTILS_HPP

#include "agile/gpu_config.hpp"
#include "agile/exception.hpp"
#include <complex>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>

namespace agile
{
  //! \brief Convert a type to a string.
  template <typename TType>
  struct type_to_file_string
  {
    static std::string get()
    {
      AGILE_ASSERT(false, StandardException::ExceptionMessage(
                             "Cannot convert this type"));
      return "";
    }
  };

  //! \brief Partial specialisation for complex types.
  template <typename TType>
  struct type_to_file_string<std::complex<TType> >
  {
    static std::string get()
    {
      return "complex " + type_to_file_string<TType>::get();
    }
  };

  //! \brief Specialisation for unsigned.
  template<>
  struct type_to_file_string<unsigned>
  {
    static std::string get() { return "uint32"; }
  };

  //! \brief Specialisation for float.
  template<>
  struct type_to_file_string<float>
  {
    static std::string get() { return "float32"; }
  };

  //! \brief Specialisation for double.
  template<>
  struct type_to_file_string<double>
  {
    static std::string get() { return "float64"; }
  };



  //! A struct holding the header of an AGILE file.
  struct FileHeader
  {
    //! Read a file header from a file stream. The file stream has to be
    //! open when this method is called.
    void read(std::ifstream& ifstr);

    //! Check if the file is of a certain type. The method will check the file
    //! type from left to right against the type supplied via \p check_type.
    //! The type hierarchy in \p check_type has to be separated by
    //! spaces. For example, if \p check_type is set to <tt>"vector dense"</tt>,
    //! the first type of the file has to be \p vector and the second has to
    //! be \p dense. From the third type on, no checks are performex. Thus,
    //! the vector could be of type \p uint32 or <tt>complex float 32</tt> and
    //! the method would return \p true in both cases.
    bool checkType(const std::string& check_type) const;

    //! The version of the file header. This string is only set if a header
    //! was read from a file. It should contain the value "1.0". For writing
    //! a header this variable is ignored.
    std::string version;

    //! The file type in a hierarchic manner separated by spaces. Examples
    //! are "vector dense uint32" or "matrix dense complex float32".
    std::string file_type;

    //! The size of the variable stored in the file in every dimension. The
    //! exact meaning depends on \p file_type. Typically, for vectors \p size
    //! will contain only one element (the vector's length). For dense matrices
    //! \p size[0] stores the number of rows, and \p size[1] holds the number
    //! of columns.
    std::vector<unsigned> size;

    //! An arbitrary one-line comment. This string must not contain a new-line
    //! character. The write function <b>won't perform sanity checks</b>.
    std::string comment;

    //! The file type split at the space delimiter. This value is only valid
    //! after reading a header from a file and is not used for writing the
    //! file header. In the latter case the type has to be supplied using
    //! \p file_type.
    std::vector<std::string> split_file_type;

    //! This parameter holds the size in bytes of one element and is valid only
    //! after reading a header from a file. For writing a header to a file,
    //! this variable is ignored.
    //! The size is determined from the last entry of the file type (e.g. if
    //! the last entry is uint32, the type size will be 4). If the last but
    //! one file type entry is "complex", the size will be doubled. Thus,
    //! "complex double" will have a \p type_size of 16.
    //! If the size could not be determined from the file type, this entry is
    //! set to zero.
    size_t type_size;
  };


  //! \brief Write the header to a stream.
  //!
  //! This function writes the common header to a stream. The stream has to
  //! be openend before calling this function.
  //! \param[in,out] ofstr The output stream to write to. This stream has to
  //! be open when calling this function.
  //! \param[in] type The file type to be written.
  //! \param[in] size A vector of size information.
  //! \param[in] comment An optional comment to be written.
  void writeHeader(std::ofstream& ofstr, const std::string& type,
                   const std::vector<unsigned>& size,
                   const std::string& comment = "");

} // namespace agile

#endif // AGILE_IO_IO_UTILS_HPP

// End of $Id: io_utils.hpp 476 2011-06-16 08:54:14Z freiberger $.
