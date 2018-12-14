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

// $Id: gpu_matrix_pitched.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_MATRIX_PITCHED_HPP
#define AGILE_GPU_MATRIX_PITCHED_HPP

#include "agile/gpu_config.hpp"
#include "agile/gpu_complex.hpp"
#include "agile/gpu_vector.hpp"
#include <complex>

namespace agile
{
  //! \brief A dense matrix on the GPU.
  template <typename TType>
  class GPUMatrixPitched
  {
    public:
      //! The type of the elements in this vector.
      typedef TType value_type;

      //! \brief Default constructor.
      //!
      //! The default constructor creates an empty GPU matrix.
      GPUMatrixPitched()
        : m_num_rows(0), m_num_columns(0), m_pitch(0), m_pitch_elements(0),
          m_data(0)
      {
      }

      //! \brief Construct a pitched GPU matrix.
      //!
      //! The rows of the matrix will be zero-padded such that the GPU access
      //! to the elements is efficient.
      //! \param[in] num_rows The number of rows.
      //! \param[in] num_columns The number of columns.
      //! \param[in] data Pointer to the matrix data. If you supply a
      //! null-pointer, the memory assigned to the matrix will be set to zero.
      //! Otherwise, this has to be a pointer to a memory area of size
      //! \p num_rows * \p num_columns on the host. The GPU matrix will
      //! be initialized from this memory area. The elements have to be
      //! stored in row-major order.
      GPUMatrixPitched(unsigned num_rows, unsigned num_columns,
                const TType* data)
        : m_num_rows(num_rows), m_num_columns(num_columns),
          m_pitch(0), m_pitch_elements(0), m_data(0)
      {
        if (m_num_rows * m_num_columns)
        {
          CUDA_SAFE_CALL(cudaMallocPitch((void**)&m_data, &m_pitch,
            m_num_columns * sizeof(TType), m_num_rows));
          m_pitch_elements = m_pitch / sizeof(TType);

          // the padded elements have to be zero -> zero the whole memory block
          lowlevel::setVectorConstant(TType(0), m_data,
                                      m_num_rows * m_pitch / sizeof(TType));

          if (data)
          {
            CUDA_SAFE_CALL(cudaMemcpy2D(
              m_data, m_pitch,
              data, m_num_columns * sizeof(TType),
              m_num_columns * sizeof(TType), m_num_rows,
              cudaMemcpyHostToDevice));
          }
        }
      }

      //! \brief Destructor.
      virtual ~GPUMatrixPitched()
      {
        freeMemory();
      }

      //! \brief Assign a matrix on the host to a GPU matrix.
      //!
      //! \param[in] num_rows The number of rows.
      //! \param[in] num_columns The number of columns.
      //! \param[in] data Pointer to the matrix data. If you supply a
      //! null-pointer, the memory assigned to the matrix will be set to zero.
      //! Otherwise, this has to be a pointer to a memory area of size
      //! \p num_rows * \p num_columns on the host. The GPU matrix will
      //! be initialized from this memory area. The elements have to be
      //! stored in row-major order.
      void assignFromHost(unsigned num_rows, unsigned num_columns,
                          const TType* data)
      {
        resize(num_rows, num_columns);

        if (m_num_rows * m_num_columns && data)
        {
          CUDA_SAFE_CALL(cudaMemcpy2D(
            m_data, m_pitch,
            data, m_num_columns * sizeof(TType),
            m_num_columns * sizeof(TType), m_num_rows,
            cudaMemcpyHostToDevice));
        }
      }

      //! \brief Copy the matrix from the GPU memory to the host (std::vector).
      //!
      //! This method copies the GPU memory holding the current GPU matrix
      //! back to the host. The host vector is resized such that its byte-size
      //! is large enough to hold the GPU memory allocated by the GPU matrix.
      //! This means that if you copy a GPUMatrixPitched<float> to a
      //! std::vector<double>, for example, the host vector will be resized
      //! to half of the size of the GPU matrix because sizeof(double) ==
      //! 2*sizeof(float). Then two floats from the GPU are copied into
      //! one double of the host meaning that you have to explicitly cast
      //! the host vector in order to get the content of the GPU matrix.
      //! Therefore it is recommended to perform copies between the same
      //! GPU matrix and host vector types only.
      //! \param[out] host_vector The host vector into which the GPU matrix
      //! will be copied.
      template <typename THostVector>
      void copyToHost(THostVector& host_vector) const
      {
        typedef typename THostVector::value_type host_type;
        unsigned byte_size = m_num_rows * m_num_columns * sizeof(TType);
        host_vector.resize((byte_size + sizeof(host_type) - 1)
                           / sizeof(host_type));
        CUDA_SAFE_CALL(cudaMemcpy2D(
          &host_vector[0], m_num_columns * sizeof(TType),
          m_data, m_pitch,
          m_num_columns * sizeof(TType), m_num_rows,
          cudaMemcpyDeviceToHost));
      }

      //! \brief Get a pointer to the data.
      //!
      //! \return A pointer to the beginning of the GPU matrix's memory.
      TType* data()
      {
        return m_data;
      }

      //! \brief Get a constant pointer to the data.
      //!
      //! \return A pointer to the beginning of the GPU matrix's memory.
      const TType* data() const
      {
        return m_data;
      }

      //! \brief Get the amount of columns.
      //!
      //! \return The number of columns in this matrix.
      unsigned getNumColumns() const
      {
        return m_num_columns;
      }

      //! \brief Get the amount of rows.
      //!
      //! \return The row-size of the matrix:.
      unsigned getNumRows() const
      {
        return m_num_rows;
      }

      //! \brief Get pitch.
      //!
      //! This method returns the number of bytes that one row of the GPU
      //! matrix occupies in the GPU memory. This value is usually larger
      //! than <tt>sizeof(TType) * num_columns</tt> because the rows are
      //! aligned such that access to the elements is efficient.
      //! \return The real width of one row in memory (in bytes).
      size_t getPitch() const
      {
        return m_pitch;
      }

      //! \brief Get real number of elements per row (alignment)
      //!
      //! This method returns the number of elements that fit into the
      //! memory that is occupied for a single row of the matrix. This is
      //! simply <tt>getPitch() / sizeof(TType)</tt>. Note that there will
      //! be still a gap in memory, if the pitch is not an integer multiple
      //! of <tt>sizeof(TType)</tt>.
      //! \return The real number of elements in memory per row (zero padded)
      unsigned getPitchElements() const
      {
        return m_pitch_elements;
      }

      //! \brief Resize the GPU matrix.
      //!
      //! The GPU matrix will be resized to hold the specified amount
      //! of \p rows and \p columns. All data stored currently in the
      //! matrix will be lost. The matrix' memory is not initialized.
      void resize(unsigned num_rows, unsigned num_columns)
      {
        if (num_rows != m_num_rows || num_columns != m_num_columns)
        {
          freeMemory();
          m_num_rows = num_rows;
          m_num_columns = num_columns;
          m_pitch = 0;
          m_pitch_elements = 0;

          if (m_num_rows * m_num_columns)
          {
            CUDA_SAFE_CALL(cudaMallocPitch((void**)&m_data, &m_pitch,
              m_num_columns * sizeof(TType), m_num_rows));
            m_pitch_elements = m_pitch / sizeof(TType);

            // the padded elements have to be zero -> zero the whole
            // memory block
            lowlevel::setVectorConstant(
              TType(0), m_data, m_num_rows * m_pitch / sizeof(TType));
          }
        }
      }

    private:
      //! \brief Number of matrix rows.
      unsigned m_num_rows;

      //! \brief Number of elements per row.
      unsigned m_num_columns;

      //! \brief Number of bytes that one row of the matrix occupies in memory.
      size_t m_pitch;

      //! \brief Number of elements per row (aligned).
      unsigned m_pitch_elements;

      //! \brief Pointer to GPU memory containing the matrix elements.
      TType* m_data;

      //! \brief Free GPU memory if allocated.
      void freeMemory()
      {
        if (m_data)
        {
          CUDA_SAFE_CALL(cudaFree(m_data));
          m_data = 0;
        }
      }
  };

  // ================================ functions ===============================

  //! Used to shift between centered, and not centered Fast Fourier Transform:
  //! Usually, the CUDA FFT considers the matrix element [0,0] to be the
  //! DC-part (center of kSpace). For visual display, a commonly used convention
  //! is to display the DC-part in the center. This is performed by function
  //! fftshift, which allows to toggle between the two versions.
  //! \param[in] M Centered/Not Centered kSpace matrix.
  //! \param[out] M Not Centered/Centered kSpace matrix, inplace computation.
  template <typename TType1>
  void fftshift(GPUMatrixPitched<TType1>& M)
  {
    lowlevel::fftshift(M.data(), M.getNumRows(), M.getNumColumns());
  }

  template <typename TType1>
  void ifftshift(GPUMatrixPitched<TType1>& M)
  {
    lowlevel::ifftshift(M.data(), M.getNumRows(), M.getNumColumns());
  }


  //! \brief Interpolate a dense GPU matrix with a GPU vector.
  //!
  //! The matrix and the vectors have to have the correct dimensions when
  //! calling this function.
  //! \param[in] M Interpolation source matrix.
  //! \param[in] pos Complex vector containing the interpolation positions.
  //! \param[out] res The interpolation result vector.
  template <typename TType1, typename TType2>
  void interp2d(const GPUMatrixPitched<TType1>& M,
                const GPUVector<std::complex<TType2> >& pos,
                GPUVector<TType1>& res);

  //! \brief Multiply a dense GPU matrix with a GPU vector.
  //!
  //! The matrix and the vectors have to have the correct dimensions when
  //! calling this function.
  //! \param[in] A Matrix to multiply with the vector.
  //! \param[in] x Vector to be multiplied with the matrix.
  //! \param[out] y The multiplication result.
  template <typename TType1, typename TType2>
  void multiply(const GPUMatrixPitched<TType1>& A,
                const GPUVector<TType2>& x,
                GPUVector<typename agile::promote<TType1, TType2>::type>& y);

  //! \brief Hermitian dense GPU matrix/vector product.
  //!
  //! This function carries out the hermitian matrix vector product
  //! \f$ y \leftarrow A^H x = \bar A^T x \f$, with a dense GPU matrix \p A and
  //! a GPU vector \p x. For real matrices this operation reduces to the
  //! multiplication with a transposed matrix:
  //! \f$ y \leftarrow A^T x \f$.
  //! \param[in] x Vector for multiplication.
  //! \param[in] A Matrix for multiplication.
  //! \param[out] y The multiplication result.
  template <typename TType1, typename TType2>
  void multiply(const GPUVector<TType1>& x, const GPUMatrixPitched<TType2>& A,
                GPUVector<typename promote<TType1, TType2>::type>& y);

  //! \brief Element-wise multiplication of two dense GPU matrices.
  //!
  //! \param[in] A First matrix.
  //! \param[in] B Second matrix.
  //! \param[in] Z Resultant element-wise product.
  template <typename TType1, typename TType2>
  void multiplyElementwise(
    const GPUMatrixPitched<TType1>& A, const GPUMatrixPitched<TType2>& B,
    GPUMatrixPitched<typename promote<TType1, TType2>::type>& Z);

  //! \brief Multiply a GPU matrix with a scalar.
  //!
  //! \param[in] alpha A scalar factor.
  //! \param[in] x A GPU matrix.
  //! \param[out] y The scaled GPU matrix alpha * x.
  template <typename TType1, typename TType2>
  void scale(const TType1& alpha, const GPUMatrixPitched<TType2>& A,
             GPUMatrixPitched<typename promote<TType1, TType2>::type>& B);

} // namespace agile

#endif // AGILE_GPU_MATRIX_PITCHED_HPP

// End of $Id: gpu_matrix_pitched.hpp 476 2011-06-16 08:54:14Z freiberger $.
