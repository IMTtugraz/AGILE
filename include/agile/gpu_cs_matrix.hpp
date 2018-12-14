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

// $Id: gpu_cs_matrix.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_CS_MATRIX_HPP
#define AGILE_GPU_CS_MATRIX_HPP

#include "agile/gpu_config.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/radix_exchange_sort.hpp"

#include <algorithm>  // max_element

namespace agile
{
  // forward declarations
  template <typename TSizeVector, typename TValueVector>
  void alignCSforGPU(const TSizeVector& major_nnz,
                     const TSizeVector& minor_index,
                     const TValueVector& data,
                     unsigned major_block_size, unsigned stride,
                     TSizeVector& aligned_major_offset,
                     TSizeVector& aligned_minor_index,
                     TValueVector& aligned_data);

  //! \brief A compressed row/column storage matrix on the GPU.
  //!
  //! This class implements compressed row storage (CRS) and compressed
  //! column storage (CCS) matrices on the GPU. If the template argument
  //! \p TIsCCS is set to false, this will be a CRS matrix (default) otherwise
  //! it will be a CCS matrix.
  template <typename TType, bool TIsCCS = false>
  class GPUCSMatrix
  {
    public:
      //! \brief Construct an empty CRS matrix.
      GPUCSMatrix() {}

      //! \brief Construct a CRS/CCS matrix.
      //!
      //! This constructor automatically re-aligns the input vectors such that
      //! multiplications can be performed efficiently on the GPU.
      //! \param[in] major_nnz Vector holding the amount of non-zero entries
      //! in the major direction. For CRS, this is a vector containing the nnz
      //! for every row, for CCS it is a vector containing the nnz for every
      //! column. This vector also specifies the amount of rows (CRS) and
      //! columns (CCS) in the matrix by its length.
      //! \param[in] minor_index A vector holding the minor index of each
      //! element. For CRS this is the column index of every entry, for CCS it
      //! is the row index. The length has to be equal to the sum of all
      //! entries of \p major_nnz.
      //! \param[in] data A vector holding the actual matrix entries. The length
      //! has to be equal to the sum of all entries of \p major_nnz.
      //! \param[in] minor_size The size in the minor direction. For a CRS
      //! matrix this is the number of columns, for a CCS matrix it is the
      //! number of rows. The value defaults to \p unsigned(-1) which means
      //! that the size is set to the maximum value of \p minor_index + 1.
      template <typename TSizeVector, typename TValueVector>
      GPUCSMatrix(const TSizeVector& major_nnz,
                  const TSizeVector& minor_index,
                  const TValueVector& data,
                  unsigned minor_size = unsigned(-1))
        : m_major_nnz(), m_major_offset(), m_minor_index(), m_data(),
          m_minor_size(minor_size)
      {
        // create temporary vectors
        TSizeVector aligned_major_offset;
        TSizeVector aligned_minor_index;
        TValueVector aligned_data;
        alignCSforGPU(major_nnz, minor_index, data, CRS_BLOCK_SIZE,
                      CRS_BLOCK_SIZE, aligned_major_offset,
                      aligned_minor_index, aligned_data);

        // assign to GPU vectors
        m_major_nnz.assignFromHost(major_nnz.begin(), major_nnz.end());
        m_major_offset.assignFromHost(aligned_major_offset.begin(),
                                      aligned_major_offset.end());
        m_minor_index.assignFromHost(aligned_minor_index.begin(),
                                     aligned_minor_index.end());
        m_data.assignFromHost(aligned_data.begin(), aligned_data.end());

        if (m_minor_size == unsigned(-1))
          m_minor_size = 1 + *std::max_element(minor_index.begin(),
                                               minor_index.end());
      }

      //! \brief Construct a CRS/CCS matrix from coordinate storage.
      //!
      //! This method constructs a CRS or a CCS matrix using three vectors:
      //! one for the row indices, one for the column indices and one for the
      //! values. All three vectors have to have the same length (because
      //! \p row[i] is the row index of the \p i-th entry, \p column[i] is
      //! its column index and \p data[i] its value).
      //! \param[in] row Vector containing the row indices.
      //! \param[in] column Vector containing the column indices.
      //! \param[in] data Vector containing the data entries.
      //! \param[in] minor_size The size in the minor direction. For a CRS
      //! matrix this is the number of columns, for a CCS matrix it is the
      //! number of rows. The value defaults to \p unsigned(-1) which means
      //! that the size is set to the maximum value of \p column_index + 1
      //! (for a CRS matrix) and to maximum of \p row_index + 1 (for a CCS
      //! matrix), respectively.
      template <typename TSizeVector, typename TValueVector>
      void createFromCoordinateStorage(
        const TSizeVector& row, const TSizeVector& column,
        const TValueVector& data, unsigned minor_size = unsigned(-1))
      {
        AGILE_ASSERT(row.size() == column.size(),
                      ExceptionSizeMismatch("row", "column",
                                            row.size(), column.size()));
        AGILE_ASSERT(row.size() == data.size(),
                      ExceptionSizeMismatch("row", "data",
                                            row.size(), data.size()));

        // special case if one wants to create an empty matrix
        if (row.size() == 0)
        {
          m_major_nnz.clear();
          m_major_offset.clear();
          m_minor_index.clear();
          m_data.clear();
          return;
        }

        // make temporaries; sorted_major is the major for CRS and the minor
        // for CCS, respectively. sorted_minor is vice versa
        TSizeVector sorted_major(TIsCCS ? column : row);
        TSizeVector sorted_minor(TIsCCS ? row : column);
        TValueVector sorted_data(data);
        m_minor_size = minor_size;
        if (m_minor_size == unsigned(-1))
          m_minor_size = 1 + *std::max_element(sorted_minor.begin(),
                                               sorted_minor.end());

        // sort by the major indices and copy the corresponding minor indices
        // as well as the data
        RadixExchangeSortCopyCopy<
          typename TSizeVector::iterator, typename TSizeVector::iterator,
          typename TValueVector::iterator>::sort(
            sorted_major.begin(), sorted_major.end(),
            sorted_minor.begin(), sorted_minor.end(),
            sorted_data.begin(), sorted_data.end());

        // sort every major slice by the minor indices (for CRS sort every row
        // by the column indices, for CCS sort every column by the row indices)
        unsigned major_begin = 0;
        unsigned major_end = 0;
        while (major_begin < sorted_major.size())
        {
          // advance as long as both iterators point to the same major
          while (major_end < sorted_major.size()
                 && sorted_major[major_end] == sorted_major[major_begin])
            ++major_end;
          // sort this major slice by the minor indices and copy the data
          RadixExchangeSortCopy<typename TSizeVector::iterator,
                                typename TValueVector::iterator>::sort(
            sorted_minor.begin() + major_begin,
            sorted_minor.begin() + major_end,
            sorted_data.begin() + major_begin,
            sorted_data.begin() + major_end);
          // beginning of the next major slice is the end of the previous one
          major_begin = major_end;
        }

        // compute the amount of non-zeros
        TSizeVector major_nnz(sorted_major.back() + 1);
        for (typename TSizeVector::const_iterator iter = sorted_major.begin();
             iter != sorted_major.end(); ++iter)
          ++major_nnz[*iter];
        sorted_major.clear();

        // align for the GPU
        TSizeVector aligned_major_offset;
        TSizeVector aligned_minor_index;
        TValueVector aligned_data;
        alignCSforGPU(major_nnz, sorted_minor, sorted_data, CRS_BLOCK_SIZE,
                      CRS_BLOCK_SIZE, aligned_major_offset,
                      aligned_minor_index, aligned_data);

        // assign to GPU vectors
        m_major_nnz.assignFromHost(major_nnz.begin(), major_nnz.end());
        m_major_offset.assignFromHost(aligned_major_offset.begin(),
                                      aligned_major_offset.end());
        m_minor_index.assignFromHost(aligned_minor_index.begin(),
                                     aligned_minor_index.end());
        m_data.assignFromHost(aligned_data.begin(), aligned_data.end());
      }

      //! \brief Get the vector containing the amount of non-zeros per major.
      const GPUVector<unsigned>& getMajorNNZ() const
      {
        return m_major_nnz;
      }

      //! \brief Get the vector containing the first indices of all majors.
      const GPUVector<unsigned>& getMajorOffset() const
      {
        return m_major_offset;
      }

      //! \brief Get the vector containing the minor index of every entry.
      const GPUVector<unsigned>& getMinorIndex() const
      {
        return m_minor_index;
      }

      //! \brief Get the vector containing the elements.
      const GPUVector<TType>& getData() const
      {
        return m_data;
      }

      //! \brief Get the number of rows.
      //!
      //! \return The number of rows of this matrix.
      unsigned getNumRows() const
      {
        return TIsCCS ? m_minor_size : m_major_nnz.size();
      }

      //! \brief Get the number of columns.
      //!
      //! \return The number of columns of this matrix.
      unsigned getNumColumns() const
      {
        return TIsCCS ? m_major_nnz.size() : m_minor_size;
      }

    private:
      //! \brief Vector holding amount of non-zero entries per major slice.
      GPUVector<unsigned> m_major_nnz;

      //! \brief Vector holding the first index of a major slice in \p m_data.
      GPUVector<unsigned> m_major_offset;

      //! \brief Vector holding the minor index for every element.
      GPUVector<unsigned> m_minor_index;

      //! \brief Vector holding the elements.
      GPUVector<TType> m_data;

      //! \brief The size in the minor direction.
      unsigned m_minor_size;
  };

  //! \brief Align a CRS/CCS matrix for GPU.
  //!
  //! Rows of a CRS matrix are processed in blocks of size \p block_size.
  //! The function searches the maximum amount of non-zero elements per row
  //! in each block. It then places in the new array the first elements of
  //! each row in the current block, followed by the second elements of each
  //! row and so on. This is visualized in the following diagram. The rows of
  //! the first block are numbered 0, 1, 2... and those of the second block
  //! with A, B, C...
  //! +-----------------------------------------------------------+-------------
  //! |                       Block 0                             |     Block 1
  //! +------------------+------------------+--+------------------+-------------
  //! | 1. entry in row  | 2. entry in row  |..| Last entry in row|  1. entry
  //! +---+---+---+---+--+---+---+---+---+--+--+---+---+---+---+--+---+---+-----
  //! | 0 | 1 | 2 | 3 |..| 0 | 1 | 2 | 3 |..|..| 0 | 1 | 2 | 3 |..| A | B | C ..
  //! +---+---+---+---+--+---+---+---+---+--+--+---+---+---+---+--+---+---+-----
  //!   |     stride       |
  //!   |<---------------->|
  //! The distance between elements of the same row can be specified using
  //! \p stride. Obviously it has to apply that stride >= block_size.
  //! This constraint is enforced by this function.
  //! If a row has less elements than the maximum, its entries are filled up
  //! with zeros.
  //! In case this is going to be a CCS matrix, we do the same operations on
  //! the columns. Thus, the parameters of the function are called
  //! \p major_offset and \p major_nnz, for example, which is the \p row_offset
  //! and \p row_nnz for CRS and the \p column_offset and \p column_nnz for CCS.
  //!
  //! \param[in] major_block_size The amount of major slices to consider for
  //!            re-ordering the entries.
  //! \param[in] stride The distance between elements of the same major slice
  //!            in the reshaped matrix.
  //! \param[out] major_offset The index where a major slice starts in the
  //!             vectors \p minor_index and \p data.
  //! \param[out] major_nnz The amount of (structural) non-zero entries per
  //!             slice.
  //! \param[out] minor_index The minor index of all matrix elements.
  //! \param[out] data The value of all matrix elements.
  template <typename TSizeVector, typename TValueVector>
  void alignCSforGPU(const TSizeVector& major_nnz,
                     const TSizeVector& minor_index,
                     const TValueVector& data,
                     unsigned major_block_size, unsigned stride,
                     TSizeVector& gpu_major_offset,
                     TSizeVector& gpu_minor_index,
                     TValueVector& gpu_data)
  {
    // the distance between elements has to be at least the block size
    if (stride < major_block_size)
      stride = major_block_size;

    // store where a major slice starts after the alignment
    unsigned aligned_major_offset = 0;
    // store the first index of a the current block of major slices
    unsigned aligned_block_start = 0;
    // store the maximum amount of non-zeros in the current block
    unsigned max_nnz_in_block = 0;
    gpu_major_offset.resize(major_nnz.size());
    for (unsigned major_counter = 0; major_counter < major_nnz.size();
         ++major_counter)
    {
      // if this is the first line in a new block, we have to make a jump in the
      // aligned major offset so that we can store all the elements of the last
      // block in between
      if (major_counter % major_block_size == 0)
      {
        // leave enough space to insert all the non-zero elements of the
        // current major slice block
        aligned_major_offset = aligned_block_start + stride * max_nnz_in_block;
        aligned_block_start = aligned_major_offset;
        max_nnz_in_block = 0;
      }

      gpu_major_offset[major_counter] = aligned_major_offset++;

      // look for maximum
      if (max_nnz_in_block < major_nnz[major_counter])
        max_nnz_in_block = major_nnz[major_counter];
    }
    // append space for the last block
    if (max_nnz_in_block)
      aligned_major_offset += stride * (max_nnz_in_block - 1);

    gpu_minor_index.assign(aligned_major_offset, 0);
    gpu_data.assign(aligned_major_offset, 0);
    // iterators to the current elements of minor_index and data
    typename TSizeVector::const_iterator minor_iter = minor_index.begin();
    typename TValueVector::const_iterator data_iter = data.begin();
    for (unsigned major_counter = 0; major_counter < major_nnz.size();
         ++major_counter)
      for (unsigned element_counter = 0;
           element_counter < major_nnz[major_counter]; ++element_counter)
      {
        unsigned aligned_index = gpu_major_offset[major_counter]
                                 + stride * element_counter;
        gpu_minor_index[aligned_index] = *minor_iter;
        gpu_data[aligned_index] = *data_iter;
        ++minor_iter;
        ++data_iter;
      }
  }

  namespace detail
  {
    //! \brief Multiply a GPU CRS matrix with a GPU vector.
    template <typename TType1, typename TType2>
    void multiplyCRS(const GPUVector<unsigned>& row_nnz,
                     const GPUVector<unsigned>& row_offset,
                     const GPUVector<unsigned>& column_index,
                     const GPUVector<TType1>& matrix_data,
                     const GPUVector<TType2>& x,
                     GPUVector<typename promote<TType1, TType2>::type>& y);

    //! \brief Hermitian product of a GPU CRS matrix with a GPU vector.
    template <typename TType1, typename TType2>
    void multiplyCRS(const GPUVector<TType1>& x,
                     const GPUVector<unsigned>& row_nnz,
                     const GPUVector<unsigned>& row_offset,
                     const GPUVector<unsigned>& column_index,
                     const GPUVector<TType2>& matrix_data,
                     GPUVector<typename promote<TType1, TType2>::type>& y);

  } // namespace detail

  //! \brief Multiply a GPU CRS/CCS matrix with a GPU vector.
  //!
  //! The matrix and the vectors have to have the correct dimensions when
  //! calling this function.
  //! \param[in] A CRS matrix to multiply with the vector.
  //! \param[in] x Vector to be multiplied with the matrix.
  //! \param[out] y The multiplication result (\f$ y \leftarrow A x \f$).
  template <typename TType1, bool TIsCCS, typename TType2>
  void multiply(const GPUCSMatrix<TType1, TIsCCS>& A,
                const GPUVector<TType2>& x,
                GPUVector<typename agile::promote<TType1, TType2>::type>& y)
  {
    // compute the normal product for CRS and the hermitian for CCS
    if (!TIsCCS)
      detail::multiplyCRS(A.getMajorNNZ(), A.getMajorOffset(),
                          A.getMinorIndex(), A.getData(), x, y);
    else
      detail::multiplyCRS(x, A.getMajorNNZ(), A.getMajorOffset(),
                          A.getMinorIndex(), A.getData(), y);
  }

  //! \brief Hermitian product of a GPU CRS/CCS matrix with a GPU vector.
  //!
  //! The matrix and the vectors have to have the correct dimensions when
  //! calling this function.
  //! \param[in] A CRS matrix to multiply with the vector.
  //! \param[in] x Vector to be multiplied with the matrix.
  //! \param[out] y The multiplication result (\f$ y \leftarrow A^H x \f$).
  template <typename TType1, bool TIsCCS, typename TType2>
  void multiply(const GPUVector<TType1>& x,
                const GPUCSMatrix<TType2, TIsCCS>& A,
                GPUVector<typename agile::promote<TType1, TType2>::type>& y)
  {
    if (!TIsCCS)
      detail::multiplyCRS(x, A.getMajorNNZ(), A.getMajorOffset(),
                          A.getMinorIndex(), A.getData(), y);
    else
      detail::multiplyCRS(A.getMajorNNZ(), A.getMajorOffset(),
                          A.getMinorIndex(), A.getData(), x, y);
  }

} // namespace agile

#endif // AGILE_GPU_CS_MATRIX_HPP

// End of $Id: gpu_cs_matrix.hpp 476 2011-06-16 08:54:14Z freiberger $.
