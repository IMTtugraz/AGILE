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

// $Id: file.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_IO_FILE_HPP
#define AGILE_IO_FILE_HPP

#include "agile/gpu_type_traits.hpp"

#include <complex>
#include <fstream>
#include <iostream>
#include <string.h>

#define AGILE_IO_FILE_FLAG_COMPLEX 1

namespace agile
{

  //! \brief Load a REAL CRS/CCS matrix from a binary file.
  //!
  //! The binary crs/ccs matrix file format is defined as follows:
  //! 1x char ....... lowest bit set indicates a complex matrix
  //! 1x unsigned ... number of matrix rows (rowCnt)
  //! 1x unsigned ... number of matrix columns
  //! 1x unsigned ... number of nonzero values in matrix (nnzCnt)
  //! rowCnt x unsigned  ...  number of nonzero values in row
  //! nnzCnt x unsigned  ...  column index in row, corresponds to data list
  //! nnzCnt x double    ...  data, list of nonzero, real values
  //!
  //! \param[in] file_name Name of binary matrix file to be read
  //! \param[out] num_rows Number of matrix rows
  //! \param[out] num_columns Number of matrix columns
  //! \param[out] major_nnz List of nonzero count per row
  //! \param[out] minor_index List of column indices in row
  //! \param[out] data List of (nonzero) matrix values
  template <typename TSizeType, typename TValueType>
  bool readCSMatrixFile(const char* file_name,
                        unsigned& num_rows, unsigned& num_columns,
                        std::vector<TSizeType>& major_nnz,
                        std::vector<TSizeType>& minor_index,
                        std::vector<TValueType>& data)
  {
    std::ifstream mat_file(file_name, std::ifstream::binary);
    if (!mat_file.is_open()) {
      std::cerr << "file not found: " << file_name << std::endl;
      return false;
    }

    // read vector header, complex flag
    char is_complex;
    mat_file.read((char*)&is_complex, sizeof(is_complex));

    // read matrix dimensions
    unsigned num_nnz;
    mat_file.read((char*)&num_rows, sizeof(num_rows));
    mat_file.read((char*)&num_columns, sizeof(num_columns));
    mat_file.read((char*)&num_nnz, sizeof(num_nnz));

    double tempd;
    unsigned tempd_size = sizeof(tempd);
    unsigned tempu;
    unsigned tempu_size = sizeof(tempu);

    // read number of nonzeros per row (CRS) or column (CCS)
    major_nnz.resize(num_rows, TSizeType(0));
    for (unsigned counter = 0; counter < num_rows; ++counter)
    {
      mat_file.read((char*)&tempu, tempu_size);
      major_nnz[counter] = (TSizeType)tempu;
    }

    // read column (CRS) or row (CCS) indices
    minor_index.resize(num_nnz, TSizeType(0));
    for (unsigned counter = 0; counter < num_nnz; ++counter)
    {
      mat_file.read((char*)&tempu, tempu_size);
      minor_index[counter] = (TSizeType)tempu;
    }

    // read values
    data.resize(num_nnz, TValueType(0));
    for (unsigned counter = 0; counter < num_nnz; ++counter)
    {
      mat_file.read((char*)&tempd, tempd_size);
      data[counter] = (TValueType)tempd;
    }

    mat_file.close();
    return true;
  }

  //! \brief Load a COMPLEX CRS/CCS matrix from a binary file.
  //!
  //! The binary crs/ccs matrix file format is defined as follows:
  //! 1x char ....... lowest bit set indicates a complex matrix
  //! 1x unsigned ... number of matrix rows (rowCnt)
  //! 1x unsigned ... number of matrix columns
  //! 1x unsigned ... number of nonzero values in matrix (nnzCnt)
  //! rowCnt x unsigned  ...  number of nonzero values in row
  //! nnzCnt x unsigned  ...  column index in row, corresponds to data list
  //! nnzCnt x double    ...  data, list of nonzero, real or complex values
  //!
  //! \param[in] file_name Name of binary matrix file to be read
  //! \param[out] num_rows Number of matrix rows
  //! \param[out] num_columns Number of matrix columns
  //! \param[out] major_nnz List of nonzero count per row
  //! \param[out] minor_index List of column indices in row
  //! \param[out] data List of (nonzero) matrix values
  template <typename TSizeType, typename TValueType>
  bool readCSMatrixFile(const char* file_name,
                        unsigned& num_rows, unsigned& num_columns,
                        std::vector<TSizeType>& major_nnz,
                        std::vector<TSizeType>& minor_index,
                        std::vector<std::complex<TValueType> >& data)
  {
    std::ifstream mat_file(file_name, std::ifstream::binary);
    if (!mat_file.is_open()) {
      std::cerr << "file not found: " << file_name << std::endl;
      return false;
    }

    // read vector header, complex flag
    char is_complex;
    mat_file.read((char*)&is_complex, sizeof(is_complex));

    // read matrix dimensions
    unsigned num_nnz;
    mat_file.read((char*)&num_rows, sizeof(num_rows));
    mat_file.read((char*)&num_columns, sizeof(num_columns));
    mat_file.read((char*)&num_nnz, sizeof(num_nnz));

    double tempd;
    unsigned tempd_size = sizeof(tempd);
    unsigned tempu;
    unsigned tempu_size = sizeof(tempu);

    // read number of nonzeros per row (CRS) or column (CCS)
    major_nnz.resize(num_rows, TSizeType(0));
    for (unsigned counter = 0; counter < num_rows; ++counter)
    {
      mat_file.read((char*)&tempu, tempu_size);
      major_nnz[counter] = (TSizeType)tempu;
    }

    // read column (CRS) or row (CCS) indices
    minor_index.resize(num_nnz, TSizeType(0));
    for (unsigned counter = 0; counter < num_nnz; ++counter)
    {
      mat_file.read((char*)&tempu, tempu_size);
      minor_index[counter] = (TSizeType)tempu;
    }

    // read values
    data.resize(num_nnz, std::complex<TValueType>(0));
    if (is_complex & AGILE_IO_FILE_FLAG_COMPLEX)
    {
      std::vector<TValueType> tmp_data(2 * num_nnz, TValueType(0));
      for (unsigned counter = 0; counter < tmp_data.size(); ++counter)
      {
        mat_file.read((char*)&tempd, tempd_size);
        tmp_data[counter] = (TValueType)tempd;
      }
      for (unsigned counter = 0; counter < num_nnz; ++counter)
        data[counter] = std::complex<TValueType>(tmp_data[counter],
                                                 tmp_data[num_nnz + counter]);
    } else {
      for (unsigned counter = 0; counter < num_nnz; ++counter)
      {
        mat_file.read((char*)&tempd, tempd_size);
        data[counter] = (TValueType)tempd;
      }
    }

    mat_file.close();
    return true;
  }

//-> added Heigl

  //! \brief Load a COMPLEX 3D matrix from a binary file (real or complex).
  //!
  //! \param[in] file_name Name of binary vector file to be read
  //! \param[out] num_rows Number of matrix rows
  //! \param[out] num_columns Number of matrix columns
  //! \param[out] num_coils Number of matrix coils (3rd Dimension)
  //! \param[out] data List of values
  template <typename TValueType>
  bool readMatrixFile3D(const char* file_name,
                        unsigned& num_rows, unsigned& num_columns, unsigned& num_coils,
                        std::vector<std::complex<TValueType> >& data)
  {
    std::ifstream mat_file(file_name, std::ifstream::binary);
    if (!mat_file.is_open()) {
      std::cerr << "file not found: " << file_name << std::endl;
      return false;
    }

    // read vector header
    mat_file.read((char*)&num_rows, sizeof(num_rows));
    mat_file.read((char*)&num_columns, sizeof(num_columns));
    mat_file.read((char*)&num_coils, sizeof(num_coils));

    unsigned num_bytes_per_entry;
    mat_file.read((char*)&num_bytes_per_entry, sizeof(num_bytes_per_entry));

    // read vector header, complex flag
    unsigned is_complex;
    mat_file.read((char*)&is_complex, sizeof(is_complex));

    unsigned size = num_rows * num_columns * num_coils
                                * num_bytes_per_entry
                                * (is_complex ? 2 : 1);
    unsigned matrix_size = num_rows * num_columns * num_coils;

    data.resize(matrix_size, std::complex<TValueType>(0));

    mat_file.read((char*)&data[0], size);
    mat_file.close();
    return true;
  }
//<-added Heigl
  
  /*
  //! \brief Load a REAL vector from a cfl file (real or complex).
  //!
  //! \param[in] file_name Name of cfl vector file to be read
  //! \param[out] data List of values
  template <typename TValueType>
  bool readCflFile(const char* file_name, unsigned size,
                      std::vector<TValueType>& data)
  {
    std::ifstream vec_file(file_name, std::ifstream::binary);
    if (!vec_file.is_open()) {
      std::cerr << "file not found: " << file_name << std::endl;
      return false;
    }

    double tempd;
    unsigned tempd_size = sizeof(tempd);

    // read vector data
    data.resize(size, TValueType(0));
    for (unsigned counter = 0; counter < size; ++counter)
    {
      vec_file.read((char*)&tempd, tempd_size);
      data[counter] = (TValueType)tempd;
    }
    vec_file.close();
    return true;
  }
  */

  //! \brief Load a REAL vector from a binary file (real or complex).
  //!
  //! \param[in] file_name Name of binary vector file to be read
  //! \param[out] data List of values
  template <typename TValueType>
  bool readVectorFile(const char* file_name,
                      std::vector<TValueType>& data)
  {
    std::ifstream vec_file(file_name, std::ifstream::binary);
    if (!vec_file.is_open()) {
      std::cerr << "file not found: " << file_name << std::endl;
      return false;
    }

    // read vector header, complex flag
    char is_complex;
    vec_file.read((char*)&is_complex, sizeof(is_complex));

    // read vector size
    unsigned size;
    vec_file.read((char*)&size, sizeof(size));

    double tempd;
    unsigned tempd_size = sizeof(tempd);

    // read vector data
    data.resize(size, TValueType(0));
    for (unsigned counter = 0; counter < size; ++counter)
    {
      vec_file.read((char*)&tempd, tempd_size);
      data[counter] = (TValueType)tempd;
    }
    vec_file.close();
    return true;
  }

  /*
  //! \brief Load a COMPLEX vector from a cfl file (real or complex).
  //!
  //! \param[in] file_name Name of cfl vector file to be read
  //! \param[out] data List of values
  template <typename TValueType>
  bool readCflFile(const char* file_name, unsigned size,
                      std::vector<std::complex<TValueType> >& data)
  {
    std::ifstream vec_file(file_name, std::ifstream::binary);
    if (!vec_file.is_open()) {
      std::cerr << "file not found: " << file_name << std::endl;
      return false;
    }

    // read vector size
    //unsigned size;
    //vec_file.read((char*)&size, sizeof(size));

    double tempd;
    unsigned tempd_size = sizeof(tempd);

    // read vector data
    data.resize(size, std::complex<TValueType>(0));
    if (AGILE_IO_FILE_FLAG_COMPLEX)
    {
      std::vector<TValueType> tmp_data(2 * size, TValueType(0));
      for (unsigned counter = 0; counter < tmp_data.size(); ++counter)
      {
        vec_file.read((char*)&tempd, tempd_size);
        tmp_data[counter] = (TValueType)tempd;
      }
      for (unsigned counter = 0; counter < size; ++counter)
        data[counter] = std::complex<TValueType>(tmp_data[counter],
                                                 tmp_data[size + counter]);
    } else {
      for (unsigned counter = 0; counter < size; ++counter)
      {
        vec_file.read((char*)&tempd, tempd_size);
        data[counter] = (TValueType)tempd;
      }
    }
    vec_file.close();
    return true;
  }
 */

  //! \brief Load a COMPLEX vector from a binary file (real or complex).
  //!
  //! \param[in] file_name Name of binary vector file to be read
  //! \param[out] data List of values
  template <typename TValueType>
  bool readVectorFile(const char* file_name,
                      std::vector<std::complex<TValueType> >& data)
  {
    std::ifstream vec_file(file_name, std::ifstream::binary);
    if (!vec_file.is_open()) {
      std::cerr << "file not found: " << file_name << std::endl;
      return false;
    }

    // read vector header, complex flag
    char is_complex;
    vec_file.read((char*)&is_complex, sizeof(is_complex));

    // read vector size
    unsigned size;
    vec_file.read((char*)&size, sizeof(size));

    double tempd;
    unsigned tempd_size = sizeof(tempd);

    // read vector data
    data.resize(size, std::complex<TValueType>(0));
    if (is_complex & AGILE_IO_FILE_FLAG_COMPLEX)
    {
      std::vector<TValueType> tmp_data(2 * size, TValueType(0));
      for (unsigned counter = 0; counter < tmp_data.size(); ++counter)
      {
        vec_file.read((char*)&tempd, tempd_size);
        tmp_data[counter] = (TValueType)tempd;
      }
      for (unsigned counter = 0; counter < size; ++counter)
        data[counter] = std::complex<TValueType>(tmp_data[counter],
                                                 tmp_data[size + counter]);
    } else {
      for (unsigned counter = 0; counter < size; ++counter)
      {
        vec_file.read((char*)&tempd, tempd_size);
        data[counter] = (TValueType)tempd;
      }
    }
    vec_file.close();
    return true;
  }

  //! \brief Writes a REAL std::vector to a binary file.
  //!
  //! \param[in] file_name Name of binary vector file
  //! \param[out] data List of values
  template <typename TValueType>
  bool writeVectorFile(const char* file_name,
                       const std::vector<TValueType> &data)
  {
    std::ofstream vec_file(file_name, std::ofstream::binary);
    if (!vec_file.is_open()) {
      std::cerr << "not able to write to file: " << file_name << std::endl;
      return false;
    }

    // write vector header, complex flag
    unsigned char is_complex = 0;
    vec_file.write((char*)&is_complex, sizeof(is_complex));

    // write vector size
    unsigned size = data.size();
    vec_file.write((char*)&size, sizeof(size));

    double tempd;
    unsigned tempd_size = sizeof(tempd);

    // write real vector data
    for (unsigned counter = 0; counter < data.size(); ++counter)
    {
      tempd = (double)data[counter];
      vec_file.write((char*)&tempd, tempd_size);
    }

    vec_file.close();
    return true;
  }


  //! \brief Writes a COMPLEX std::vector to a binary file.
  //!
  //! \param[in] file_name Name of binary vector file
  //! \param[out] data List of values
  template <typename TValueType>
  bool writeVectorFile(const char* file_name,
                       const std::vector<std::complex<TValueType> >&data)
  {
    std::ofstream vec_file(file_name, std::ofstream::binary);
    if (!vec_file.is_open()) {
      std::cerr << "not able to write to file: " << file_name << std::endl;
      return false;
    }

    // write vector header, complex flag
    unsigned char is_complex = AGILE_IO_FILE_FLAG_COMPLEX;
    vec_file.write((char*)&is_complex, sizeof(is_complex));

    // write vector size
    unsigned size = data.size();
    vec_file.write((char*)&size, sizeof(size));

    double tempd;
    unsigned tempd_size = sizeof(tempd);

    // write real vector data
    for (unsigned counter = 0; counter < data.size(); ++counter)
    {
      tempd = (double)real(data[counter]);
      vec_file.write((char*)&tempd, tempd_size);
    }

    // write complex vector data
    for (unsigned counter = 0; counter < data.size(); ++counter)
    {
      tempd = (double)imag(data[counter]);
      vec_file.write((char*)&tempd, tempd_size);
    }

    vec_file.close();
    return true;
  }


//-> added Heigl

  //! \brief Load a COMPLEX 3D matrix from a binary file (real or complex).
  //!
  //! \param[in] file_name Name of binary vector file to be read
  //! \param[out] num_rows Number of matrix rows
  //! \param[out] num_columns Number of matrix columns
  //! \param[out] num_coils Number of matrix coils (3rd Dimension)
  //! \param[out] data List of values
  template <typename TValueType>
  bool writeMatrixFile3D(const char* file_name,
                        unsigned num_rows, unsigned num_columns, unsigned num_coils,
                        const std::vector<std::complex<TValueType> >& data)
  {
    std::ofstream mat_file(file_name, std::ofstream::binary);
    if (!mat_file.is_open()) {
      std::cerr << "not able to write to file: " << file_name << std::endl;
      return false;
    }

    // write vector header, complex flag

    mat_file.write((char*)&num_rows, sizeof(num_rows));
    mat_file.write((char*)&num_columns, sizeof(num_columns));
    mat_file.write((char*)&num_coils, sizeof(num_coils));
    unsigned num_bytes_per_entry=sizeof(num_bytes_per_entry);
    mat_file.write((char*)&num_bytes_per_entry, sizeof(num_bytes_per_entry));
    unsigned is_complex = AGILE_IO_FILE_FLAG_COMPLEX;
    mat_file.write((char*)&is_complex, sizeof(is_complex));

    float tempd;
    unsigned tempd_size = sizeof(tempd);

    // write real vector data
    for (unsigned counter = 0; counter < data.size(); ++counter)
    {
      tempd = (float)real(data[counter]);
      mat_file.write((char*)&tempd, tempd_size);
      tempd = (float)imag(data[counter]);
      mat_file.write((char*)&tempd, tempd_size);

    }

    mat_file.close();
    return true;
  }

  template <typename TValueType>
  bool readMatrixFile(const char* file_name,
                        unsigned& num_rows, unsigned& num_columns,
                        std::vector<TValueType>& data)
  {
    unsigned num_coils;
    std::ifstream mat_file(file_name, std::ifstream::binary);
    if (!mat_file.is_open()) {
      std::cerr << "file not found: " << file_name << std::endl;
      return false;
    }

    // read vector header
    mat_file.read((char*)&num_rows, sizeof(num_rows));
    mat_file.read((char*)&num_columns, sizeof(num_columns));
    mat_file.read((char*)&num_coils, sizeof(num_coils));

    unsigned num_bytes_per_entry;
    mat_file.read((char*)&num_bytes_per_entry, sizeof(num_bytes_per_entry));

    // read vector header, complex flag
    unsigned is_complex;
    mat_file.read((char*)&is_complex, sizeof(is_complex));

    unsigned size = num_rows * num_columns * num_coils
                                * num_bytes_per_entry
                                * (is_complex ? 2 : 1);
    unsigned matrix_size = num_rows * num_columns * num_coils;

    data.resize(matrix_size, TValueType(0));

    mat_file.read((char*)&data[0], size);
    mat_file.close();
    return true;
  }

  //<-added Heigl


} // namespace agile

#endif // AGILE_IO_FILE_HPP

// End of $Id: file.hpp 476 2011-06-16 08:54:14Z freiberger $.
