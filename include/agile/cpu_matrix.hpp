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

// $Id: cpu_matrix.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_CPU_MATRIX_HPP
#define AGILE_CPU_MATRIX_HPP

#include "agile/cpu_vector.hpp"
#include "agile/exception.hpp"
#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"

#include <cstring> // for memcpy
#include <vector>

namespace agile
{
  //! \brief A dense matrix on the host.
  template <typename TType>
  class CPUMatrix
  {
    public:
      //! \brief Construct a matrix.
      //!
      //! \param[in] num_rows The number of rows.
      //! \param[in] num_columns The number of columns.
      //! \param[in] data Pointer to the matrix data. If you supply a
      //! null-pointer, an empty matrix is created. Otherwise, this memory
      //! area has to be of size \p num_rows * \p num_columns. The elements
      //! have to be stored in row-major order.
      CPUMatrix(unsigned num_rows, unsigned num_columns,
                const TType* data)
        : m_num_rows(num_rows), m_num_columns(num_columns)
      {
        if (data)
        {
          m_data.resize(num_rows * num_columns);
          memcpy(&m_data[0], data, num_rows * num_columns * sizeof(TType));
        }
        else
          m_data.assign(num_rows * num_columns, 0);
      }

      //! \brief Get the amount of rows.
      unsigned getNumRows() const
      {
        return m_num_rows;
      }

      //! \brief Get the amount of columns.
      unsigned getNumColumns() const
      {
        return m_num_columns;
      }

      //! \brief Get a pointer to the data.
      //!
      //! \return A pointer to the beginning of the GPU vector's memory.
      TType* data()
      {
        return &m_data[0];
      }

      //! \brief Get a constant pointer to the data.
      const TType* data() const
      {
        return &m_data[0];
      }

    private:

      //! \brief Number of matrix rows.
      unsigned m_num_rows;

      //! \brief Number of elements per row.
      unsigned m_num_columns;

      //! \brief Vector holding the elements.
      std::vector<TType> m_data;
  };

  //! \brief Multiply a dense CPU matrix with a CPU vector.
  //!
  //! The matrix and the vectors have to have the correct dimensions when
  //! calling this function.
  //! \param[in] A Matrix to multiply with the vector.
  //! \param[in] x Vector to be multiplied with the matrix.
  //! \param[out] y The multiplication result.
  template <typename TType1, typename TType2>
  void multiply(const CPUMatrix<TType1>& A,
                const CPUVector<TType2>& x,
                CPUVector<typename agile::promote<TType1, TType2>::type>& y)
  {
    AGILE_ASSERT(A.getNumColumns() == x.size(),
                  StandardException::
                    ExceptionSizeMismatch("columns of A", "x",
                                          A.getNumColumns(), x.size()));
    AGILE_ASSERT(A.getNumRows() == y.size(),
                  StandardException::
                    ExceptionSizeMismatch("rows of A", "y",
                                          A.getNumRows(), y.size()));

    unsigned num_columns = A.getNumColumns();
    for (unsigned row_counter = 0; row_counter < A.getNumRows(); ++row_counter)
    {
      typename agile::promote<TType1, TType2>::type sum = 0;
      for (unsigned column_counter = 0; column_counter < num_columns;
           ++column_counter)
      {
        unsigned offset = row_counter * num_columns + column_counter;
        sum += A.data()[offset] * x[column_counter];
      }
      y[row_counter] = sum;
    }
  }

  //! \brief Hermitian dense CPU matrix/vector product.
  //!
  //! This function carries out the hermitian matrix vector product
  //! \f$ y \leftarrow A^H x = \bar A^T x \f$, with a dense CPU matrix \p A and
  //! a CPU vector \p x. For real matrices this operation reduces to the
  //! multiplication with a transposed matrix:
  //! \f$ y \leftarrow A^T x \f$.
  //! \param[in] x Vector for multiplication.
  //! \param[in] A Matrix for multiplication.
  //! \param[out] y The multiplication result.
  template <typename TType1, typename TType2>
  void multiply(const CPUVector<TType1>& x, const CPUMatrix<TType2>& A,
                CPUVector<typename promote<TType1, TType2>::type>& y)
  {
    AGILE_ASSERT(A.getNumColumns() == y.size(),
                  StandardException::
                    ExceptionSizeMismatch("columns of A", "y",
                                          A.getNumColumns(), y.size()));
    AGILE_ASSERT(A.getNumRows() == x.size(),
                  StandardException::
                    ExceptionSizeMismatch("rows of A", "x",
                                          A.getNumRows(), x.size()));

    unsigned num_columns = A.getNumColumns();
    for (unsigned row_counter = 0; row_counter < A.getNumRows(); ++row_counter)
    {
      typename agile::promote<TType1, TType2>::type sum = 0;
      for (unsigned column_counter = 0; column_counter < num_columns;
           ++column_counter)
      {
        unsigned offset = column_counter * num_columns + row_counter;
        sum += A.data()[offset] * x[column_counter];
      }
      y[row_counter] = sum;
    }
  }

  //! \brief Element-wise multiplication of two dense CPU matrices.
  //!
  //! \param[in] A First matrix.
  //! \param[in] B Second matrix.
  //! \param[in] Z Resultant element-wise product.
  template <typename TType1, typename TType2>
  void multiplyElementwise(
    const CPUMatrix<TType1>& A, const CPUMatrix<TType2>& B,
    CPUMatrix<typename promote<TType1, TType2>::type>& Z)
  {
    AGILE_ASSERT(A.getNumRows() == B.getNumRows(),
                  StandardException::
                    ExceptionSizeMismatch("rows of A", "rows of B",
                                          A.getNumRows(), B.getNumRows()));
    AGILE_ASSERT(A.getNumColumns() == B.getNumColumns(),
                  StandardException::
                    ExceptionSizeMismatch("columns of A", "columns of B",
                                          A.getNumColumns(), B.getNumColumns()));
    AGILE_ASSERT(A.getNumRows() == Z.getNumRows(),
                  StandardException::
                    ExceptionSizeMismatch("rows of A", "rows of Z",
                                          A.getNumRows(), Z.getNumRows()));
    AGILE_ASSERT(A.getNumColumns() == Z.getNumColumns(),
                  StandardException::
                    ExceptionSizeMismatch("columns of A", "columns of Z",
                                          A.getNumColumns(), Z.getNumColumns()));

    unsigned num_columns = A.getNumColumns();
    for (unsigned counter = 0; counter < A.getNumRows() * A.getNumColumns();
         ++counter)
      Z.data()[counter] = A.data()[counter] * B.data()[counter];
  }

} // namespace agile

#endif // AGILE_CPU_MATRIX_HPP

// End of $Id: cpu_matrix.hpp 476 2011-06-16 08:54:14Z freiberger $.
