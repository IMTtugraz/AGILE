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

// $Id: gpu_matrix.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_MATRIX_HPP
#define AGILE_GPU_MATRIX_HPP

#include "agile/gpu_config.hpp"
#include "agile/gpu_complex.hpp"
#include "agile/gpu_vector.hpp"

#include <complex>
#include <cublas_v2.h>


#include <iostream>
#include <iomanip>


namespace agile
{
  //! \brief A dense matrix on the GPU using cuBLAS.
  template <typename TType>
  class GPUMatrix
  {
    public:
      //! The type of the elements in this vector.
      typedef TType value_type;

      //! \brief Default constructor.
      //!
      //! The default constructor creates an empty GPU matrix.
      GPUMatrix()
        : m_num_rows(0), m_num_columns(0), m_data(0)
      {

      }

      //! \brief Construct a GPU matrix.
      //!
      //! \param[in] num_rows The number of rows.
      //! \param[in] num_columns The number of columns.
      //! \param[in] data Pointer to the matrix data. If you supply a
      //! null-pointer, the memory assigned to the matrix will be set to zero.
      //! Otherwise, this has to be a pointer to a memory area of size
      //! \p num_rows * \p num_columns on the host. The GPU matrix will
      //! be initialized from this memory area. The elements have to be
      //! stored in row-major order.
      GPUMatrix(unsigned num_rows, unsigned num_columns,
                const TType* data)
        : m_num_rows(num_rows), m_num_columns(num_columns), m_data(0)
      {
       if (m_num_rows * m_num_columns)
        {

           CUDA_SAFE_CALL(cudaMalloc((void**)&m_data, m_num_columns * m_num_rows * sizeof(TType) ));

          if (data)
          {

              CUBLAS_SAFE_CALL(cublasSetMatrix(m_num_columns, m_num_rows, sizeof(TType),
                   data, m_num_columns, m_data, m_num_columns));
          }
        }
      }

      /*//! \brief Copy Constructor.
      GPUMatrix(const GPUMatrix& A)
      {
        this->m_num_columns = A.getNumColumns();
        this->m_num_rows = A.getNumRows();
        this->resize(this->m_num_rows,this->m_num_columns);

        if (A.data())
        {
          CUDA_SAFE_CALL(cudaMemcpy(
            this.data(),A.data(),
            m_num_columns*m_num_rows * sizeof(TType),cudaMemcpyDeviceToDevice)
            );
        }
      }*/


      //! \brief Destructor.
      virtual ~GPUMatrix()
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
          //std::cout<<"\ndata: "<<*data<<"  - sizeof-type"<<sizeof(TType);
            CUBLAS_SAFE_CALL(cublasSetMatrix(m_num_columns, m_num_rows,  sizeof(TType),
                 data, m_num_columns, m_data, m_num_columns));
        }
      }


      //! \brief Copy the matrix from the GPU memory to the host (std::vector).
      //!
      //! This method copies the GPU memory holding the current GPU matrix
      //! back to the host. The host vector is resized such that its byte-size
      //! is large enough to hold the GPU memory allocated by the GPU matrix.
      //! This means that if you copy a GPUMatrix<float> to a
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
        CUBLAS_SAFE_CALL(cublasGetMatrix (m_num_columns, m_num_rows, sizeof(TType),
                                          m_data, m_num_columns, &host_vector[0], m_num_columns));


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

          if (m_num_rows * m_num_columns)
          {
            CUDA_SAFE_CALL(cudaMalloc((void**)&m_data, m_num_columns * m_num_rows * sizeof(TType) ));
          }
        }
      }


    private:
      //! \brief Number of matrix rows.
      unsigned m_num_rows;

      //! \brief Number of elements per row.
      unsigned m_num_columns;

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
  //! \brief Copy Matrix A to Z
  //!
  //! \param[in] A A GPU matrix.
  //! \param[out] Z A GPU matrix.
  template <typename TType1>

  void copy(const GPUMatrix<TType1>& A,
             GPUMatrix<TType1>& Z);
/*  template <typename TType1>
  void copy(const GPUMatrix<TType1>& A,
             GPUMatrix<TType1>* B);*/

  //! \brief transposed multiplication of two dense GPU matrices.
  //!
  //! \param[in] A First matrix.
  //! \param[in] B Second matrix.
  //! \param[in] Z Resultant product. Z = X' * Y
  template <typename TType1>
  inline
  void multiplyTransMatrixMatrix(const GPUMatrix<TType1>& X, const GPUMatrix<TType1>& Y,
                                 GPUMatrix<TType1>& Z);


  //! \brief Dot Product of two dense GPU matrices (read as vectors).
  //!
  //! \param[in] A First matrix.
  //! \param[in] B Second matrix.
  //! \param[in] Z Resultant product. Z = X(:)' * Y(:)
  template <typename TType1>
  void dotProduct(const GPUMatrix<TType1>& X, const GPUMatrix<TType1>& Y,
                           TType1& val);


  //! \brief Multiply a GPU matrix with a scalar.
  //!
  //! \param[in] alpha A scalar factor.
  //! \param[in] A A GPU matrix.
  //! \param[out] B The scaled GPU matrix alpha * A.
  template <typename TType1>

  void scale(const typename to_real_type<TType1>::type& alpha, const GPUMatrix<TType1>& A,
             GPUMatrix<TType1>& B);


  //! \brief Multiply a dense GPU matrix with a GPU vector.
  //!
  //! The matrix and the vectors have to have the correct dimensions when
  //! calling this function.
  //! \param[in] A Matrix to multiply with the vector.
  //! \param[in] x Vector to be multiplied with the matrix.
  //! \param[out] y The multiplication result.
  template <typename TType1>
  inline
  void multiply(const GPUMatrix<TType1>& A,
                const GPUVector<TType1>& x,
                GPUVector<TType1>& y);


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
  template <typename TType1>
  inline
  void multiply(const GPUVector<TType1>& x, const GPUMatrix<TType1>& A,
                GPUVector<TType1>& y);

/*
  //! \brief Element-wise multiplication of two dense GPU matrices.
  //!
  //! \param[in] A First matrix.
  //! \param[in] B Second matrix.
  //! \param[in] Z Resultant element-wise product.
  template <typename TType1>
  void multiplyElementwise(const GPUMatrix<TType1>& A, const GPUMatrix<TType1>& B,
                           GPUMatrix<TType1>& Z);
*/

  //! \brief Compute the l2-norm of a GPU Matrix
  //!
  //! \param[in] A A GPU matrix.
  //! \return \f$ \left ( \sum_i x_i^2 \right)^{1/2} \f$
  template <typename TType1>
  typename to_real_type<TType1>::type norm2(const GPUMatrix<TType1>& A);


  //! Used to shift between centered, and not centered Fast Fourier Transform:
  //! Usually, the CUDA FFT considers the matrix element [0,0] to be the
  //! DC-part (center of kSpace). For visual display, a commonly used convention
  //! is to display the DC-part in the center. This is performed by function
  //! fftshift, which allows to toggle between the two versions.
  //! \param[in] M Centered/Not Centered kSpace matrix.
  //! \param[out] M Not Centered/Centered kSpace matrix, inplace computation.
  template <typename TType1>
  inline
  void fftshift(GPUMatrix<TType1>& M)
  {
    lowlevel::fftshift(M.data(), M.getNumRows(), M.getNumColumns());
  }
  template <typename TType1>
  inline
  void ifftshift(GPUMatrix<TType1>& M)
  {
    lowlevel::ifftshift(M.data(), M.getNumRows(), M.getNumColumns());
  }



  //////////////////////////////
  //! \brief Conjugation of Complex Matrix Elements.
  //!
  //! Its the user's responsibility to make sure that the sizes of the Matrix
  //! are correct before calling this function.
  //! \param[in] X Input Matrix.
  //! \param[in] z Conjugate of x \p conj(X).
  template <typename TDerived1>
  inline
  void conjMatrix(const GPUMatrix<TDerived1>& X,
                               GPUMatrix<TDerived1>& Z)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Z.getNumColumns()*Z.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Z", X.getNumColumns()*X.getNumRows(), Z.getNumColumns()*Z.getNumRows()));

    lowlevel::conjVector(X.data(), Z.data(), X.getNumColumns()*X.getNumRows());
  }


  //! \brief Conj Element-wise multiplication of GPU matrices.
  //!
  //! \param[in] X First matrix.
  //! \param[in] Y Second matrix.
  //! \param[in] Z Resultant element-wise product (Z = conj(X).*Y )
  template <typename TType1, typename TType2, typename TType3>
  inline
  void multiplyConjElementwise(const GPUMatrix<TType1>& X, const GPUMatrix<TType2>& Y,
                           GPUMatrix<TType3>& Z)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Z.getNumColumns()*Z.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Z", X.getNumColumns()*X.getNumRows(), Z.getNumColumns()*Z.getNumRows()));

    lowlevel::multiplyConjElementwise(X.data(), Y.data(), Z.data(), X.getNumColumns()*X.getNumRows());
  }



  //! Calculates Absolute-Value.
  //! \param[in] A Matrix.
  //! \param[out] Y Matrix with Absolute Values.
  template <typename TType1, typename TType2>
  inline
  void absMatrix(const GPUMatrix<TType1>& A, GPUMatrix<TType2>& Y)
  {
      AGILE_ASSERT(A.getNumColumns()*A.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                    StandardException::ExceptionSizeMismatch(
                      "A", "Y", A.getNumColumns()*A.getNumRows(), Y.getNumColumns()*Y.getNumRows()));

    lowlevel::absVector(A.data(), Y.data(), A.getNumColumns()*A.getNumRows());
  }

  //////////////////////////////
  //! \brief Generate linearly spaced vector between a and b with n numbers
  //!
  //! \param[in] a start value.
  //! \param[in] b end value.
  //! \param[in] size - amount of numbers
  //! \param[out] M Matrix to be filled
  template <typename TType1>
  inline
  void linspace(GPUMatrix<TType1>& M, float a, float b)
  {
    lowlevel::linspace(M.data(), M.getNumRows() * M.getNumColumns(), a, b);

  }

  //////////////////////////////
  //! \brief Generates meshgrid Matrix Mx and My for input Vector x,y
  //!
  //! \param[in] x Vector
  //! \param[in] y Vector
  //! \param[out] Mx and My Matrix to be filled
  template <typename TType1>
  inline
  void meshgrid(GPUMatrix<TType1>& Mx, GPUMatrix<TType1>& My,
                const GPUVector<TType1> x, const GPUVector<TType1> y)
  {
    Mx.resize(y.size(),x.size());
    My.resize(y.size(),x.size());
    lowlevel::meshgrid(Mx.data(), My.data(),
                       x.data(), x.size(), y.data(), y.size());
  }


  //! \brief Add Matrix X and Y.
  //!
  //! The Matrix have to be of equal size when calling this function as they
  //! won't be resized due to performance considerations.
  //! \param[in] X First Matrix.
  //! \param[in] Y Second Matrix.
  //! \param[out] Z Sum of the two Matrix.
  template <typename TType1, typename TType2, typename TType3>
  inline
  void addMatrix(const GPUMatrix<TType1>& X,
                 const GPUMatrix<TType2>& Y,
                 GPUMatrix<TType3>& Z)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Z.getNumColumns()*Z.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Z", X.getNumColumns()*X.getNumRows(), Z.getNumColumns()*Z.getNumRows()));

    lowlevel::addVector(X.data(), Y.data(), Z.data(), X.getNumColumns()*X.getNumRows());
  }

  //! \brief Subtract Matrix X and Y.
  //!
  //! The Matrix have to be of equal size when calling this function as they
  //! won't be resized due to performance considerations.
  //! \param[in] X First Matrix.
  //! \param[in] Y Second Matrix.
  //! \param[out] Z Difference of the two Matrix.
  template <typename TType1, typename TType2, typename TType3>
  inline
  void subMatrix(const GPUMatrix<TType1>& X,
                 const GPUMatrix<TType2>& Y,
                 GPUMatrix<TType3>& Z)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Z.getNumColumns()*Z.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Z", X.getNumColumns()*X.getNumRows(), Z.getNumColumns()*Z.getNumRows()));

    lowlevel::subVector(X.data(), Y.data(), Z.data(), X.getNumColumns()*X.getNumRows());
  }


  //! \brief Compute the square root of every element of a GPU Matrix.
  //!
  //! \param[in] X The input Matrix.
  //! \param[out] Y The Matrix Y[i] <- sqrt(X[i]).
  template <typename TDerived1>
  inline
  void sqrt(const GPUMatrix<TDerived1>& X, GPUMatrix<TDerived1>& Y)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));

    lowlevel::sqrt(X.data(), Y.data(), X.getNumColumns()*X.getNumRows());
  }

  //! \brief Compute the power of alpha for every element of a GPU Matrix.
  //!
  //! \param[in] alpha power-value.
  //! \param[in] X The input Matrix.
  //! \param[out] Y The Matrix Y[i] <- (X[i])^(alpha).
  template <typename TType1, typename TType2>
  inline
  void pow(const TType1& alpha, const GPUMatrix<TType2>& X, GPUMatrix<TType2>& Y)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));

    lowlevel::pow(alpha, X.data(), Y.data(), X.getNumColumns()*X.getNumRows());
  }

  //! \brief Generate pattern of a Matrix.
  //! \brief Z = abs(X)>0;   for complex value: 1+0i
  //!
  //! \param[in] X The input Matrix.
  //! \param[out] Z The Patterned Matrix
  template <typename TDerived1>
  inline
  void pattern(const GPUMatrix<TDerived1>& X, GPUMatrix<typename to_real_type<TDerived1>::type>& Z)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Z.getNumColumns()*Z.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Z", X.getNumColumns()*X.getNumRows(), Z.getNumColumns()*Z.getNumRows()));

    lowlevel::pattern(X.data(), Z.data(), X.getNumColumns()*X.getNumRows());
  }


  //! \brief Compute the difference between each value in a GPU Matrix.
  //! \brief (last value with first value)
  //!
  //! \param[in] X - The input Matrix.
  //! \param[in] dim - dimension to be calculated (1...rowmajor, 2...columnmajor)
  //! \param[out] Y - The Matrix Y[i] <- diff(X[i]).
  template <typename TDerived1>
  inline
  void diff(const unsigned dim, const GPUMatrix<TDerived1>& X, GPUMatrix<TDerived1>& Y)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));

    bool past=true;
    if ((dim <= 0)||(dim > 2))
      past = false;
    AGILE_ASSERT(past == true,
                 StandardException::ExceptionMessage(
               "dim - error - wrong dimension for calculation of Matrix-Value-Difference"));


    lowlevel::diff(dim, X.getNumRows(), X.data(), Y.data(), X.getNumColumns()*X.getNumRows());
  }

  //! \brief Compute the difference between each value in a GPU Matrix transposed.
  //! \brief (last value with first value)
  //!
  //! \param[in] X - The input Matrix.
  //! \param[in] dim - dimension to be calculated (1...rowmajor, 2...columnmajor)
  //! \param[out] Y - The Matrix Y[i] <- difftrans(X[i]).
  template <typename TDerived1>
  inline
  void difftrans(const unsigned dim, const GPUMatrix<TDerived1>& X, GPUMatrix<TDerived1>& Y)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));

    bool past=true;
    if ((dim <= 0)||(dim > 2))
      past = false;
    AGILE_ASSERT(past == true,
                 StandardException::ExceptionMessage(
               "dim - error - wrong dimension for calculation of Matrix-Value-Difference"));


    lowlevel::difftrans(dim, X.getNumRows(), X.data(), Y.data(), X.getNumColumns()*X.getNumRows());
  }




  //! \brief generate max-value vector of two Matrics (elementwise).
  //! \brief Y = max(X1,X2);
  //!
  //! \param[in] X1 Matrix.
  //! \param[in] X2 Matrix.
  //! \param[out] Y max Matrix.
  template <typename TDerived1, typename TDerived2>
  inline
  void max(const GPUMatrix<TDerived1>& X1, const GPUMatrix<TDerived2>& X2,
            GPUMatrix<typename promote<TDerived1, TDerived2>::type>& Y)
  {
    AGILE_ASSERT(X1.getNumColumns()*X1.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X1", "Y", X1.getNumColumns()*X1.getNumRows(), Y.getNumColumns()*Y.getNumRows()));

    AGILE_ASSERT(X1.getNumColumns()*X1.getNumRows() == X2.getNumColumns()*X2.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X1", "X2", X1.getNumColumns()*X1.getNumRows(), X2.getNumColumns()*X2.getNumRows()));


    lowlevel::max(X1.data(),X2.data(), Y.data(), X1.getNumColumns()*X1.getNumRows());
  }
  
  //! \brief generate max-value vector of Matrics (elementwise) and Scalar.
  //! \brief Y[i,j] = max(X1[i,j],x2);
  //!
  //! \param[in] X1 Matrix.
  //! \param[in] x2 Scalar.
  //! \param[out] Y max Matrix.
  template <typename TDerived1, typename TDerived2>
  inline
  void max(const GPUMatrix<TDerived1>& X1, const TDerived2& x2,
            GPUMatrix<typename promote<TDerived1, TDerived2>::type>& Y)
  {
    AGILE_ASSERT(X1.getNumColumns()*X1.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X1", "Y", X1.getNumColumns()*X1.getNumRows(), Y.getNumColumns()*Y.getNumRows()));

    lowlevel::max(X1.data(),x2, Y.data(), X1.getNumColumns()*X1.getNumRows());
  }

  //! \brief Extract the real part of a GPU Matrix.
  //!
  //! \param[in] x A GPU Matrix.
  //! \param[out] y The real part of the GPU Matrix \p x.
  template <typename TDerived1, typename TDerived2>
  inline
  void real(const GPUMatrix<TDerived1>& X, GPUMatrix<TDerived2>& Y)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));

    lowlevel::real(X.data(), Y.data(), Y.getNumColumns()*Y.getNumRows());
  }

  //! \brief Extract the imag part of a GPU Matrix.
  //!
  //! \param[in] x A GPU Matrix.
  //! \param[out] y The imag part of the GPU Matrix \p x.
  template <typename TDerived1, typename TDerived2>
  inline
  void imag(const GPUMatrix<TDerived1>& X, GPUMatrix<TDerived2>& Y)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));

    lowlevel::imag(X.data(), Y.data(), Y.getNumColumns()*Y.getNumRows());
  }


  //! Multiply Matrix Matrix Element-wise (implemented in CUDA)
  //! \brief Element-wise matrix-matrix product (host function).
  //! \param[in] x First vector.
  //! \param[in] y Second vector.
  //! \param[in] z Element-wise product of \p x and \p y.
  template <typename TDerived1, typename TDerived2, typename TDerived3>
  inline
  void multiplyElementwise(const GPUMatrix<TDerived1>& X,
                           const GPUMatrix<TDerived2>& Y,
                           GPUMatrix<TDerived3>& Z)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Z.getNumColumns()*Z.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Z", X.getNumColumns()*X.getNumRows(), Z.getNumColumns()*Z.getNumRows()));

    lowlevel::multiplyElementwise(X.data(), Y.data(), Z.data(),X.getNumColumns()*X.getNumRows());
  }

  //! Divide Matrix Matrix Element-wise (implemented in CUDA)
  //! \brief Element-wise matrix-matrix product (host function).
  //! \param[in] x First vector.
  //! \param[in] y Second vector.
  //! \param[in] z Element-wise division of \p x and \p y.
  template <typename TDerived1, typename TDerived2, typename TDerived3>
  inline
  void divideElementwise(const GPUMatrix<TDerived1>& X,
                           const GPUMatrix<TDerived2>& Y,
                           GPUMatrix<TDerived3>& Z)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Z.getNumColumns()*Z.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Z", X.getNumColumns()*X.getNumRows(), Z.getNumColumns()*Z.getNumRows()));

    lowlevel::divideElementwise(X.data(), Y.data(), Z.data(),X.getNumColumns()*X.getNumRows());
  }

  //! Generate sqare Data matrix with zero-padding in row-dimension
  //! \param[in] X input data matrix.
  //! \param[in] row_o number of rows to be added on upper side
  //! \param[in] row_u number of rows to be added at the bottom
  //! \param[in] Z generated Matrix
  template <typename TDerived1>
  inline
  void expand_rowdim(const GPUMatrix<TDerived1>& X_data,
                     const GPUMatrix<TDerived1>& Delta_o,
                     const GPUMatrix<TDerived1>& Delta_u,
                           GPUMatrix<TDerived1>& Z)
  {
    AGILE_ASSERT(Z.getNumColumns()*Z.getNumRows() ==
                 X_data.getNumColumns()*X_data.getNumRows() +
                 Delta_o.getNumColumns()*Delta_o.getNumRows() +
                 Delta_u.getNumColumns()*Delta_u.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "Z", "Delta_o + X + Delta_u", Z.getNumColumns()*Z.getNumRows(), X_data.getNumColumns()*X_data.getNumRows() +
                   Delta_o.getNumColumns()*Delta_o.getNumRows() +
                   Delta_u.getNumColumns()*Delta_u.getNumRows()));

    lowlevel::expand_rowdim(X_data.data(), Delta_o.data(), Delta_u.data(),
                            X_data.getNumRows(), X_data.getNumColumns(),
                            Delta_o.getNumRows(), Delta_u.getNumRows(), Z.data());
  }

  //! Generate sqare Data matrix with zero-padding in col-dimension
  //! \param[in] X input data matrix.
  //! \param[in] col_o number of rows to be added on right side
  //! \param[in] col_u number of rows to be added left side
  //! \param[in] Z generated Matrix
  template <typename TDerived1>
  inline
  void expand_coldim(const GPUMatrix<TDerived1>& X_data,
                     const GPUMatrix<TDerived1>& Delta_o,
                     const GPUMatrix<TDerived1>& Delta_u,
                           GPUMatrix<TDerived1>& Z)
  {
    AGILE_ASSERT(Z.getNumColumns()*Z.getNumRows() ==
                 X_data.getNumColumns()*X_data.getNumRows() +
                 Delta_o.getNumColumns()*Delta_o.getNumRows() +
                 Delta_u.getNumColumns()*Delta_u.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "Z", "Delta_o + X + Delta_u", Z.getNumColumns()*Z.getNumRows(), X_data.getNumColumns()*X_data.getNumRows() +
                   Delta_o.getNumColumns()*Delta_o.getNumRows() +
                   Delta_u.getNumColumns()*Delta_u.getNumRows()));

     lowlevel::expand_coldim(X_data.data(), Delta_o.data(), Delta_u.data(),
                             X_data.getNumRows(), X_data.getNumColumns(),
                             Delta_o.getNumColumns(), Delta_u.getNumColumns(), Z.data());
  }


  //! copy a matrix-sized data from a given input matrix
  //! \param[in] X input data matrix.
  //! \param[in] row_offset - offset in row dimension
  //! \param[in] col_offset - offset in col dimension
  //! \param[in] Z generated Matrix
  template <typename TDerived1>
  inline
  void get_content(const GPUMatrix<TDerived1>& X_data,
                     unsigned row_offset, unsigned col_offset,
                           GPUMatrix<TDerived1>& Z)
  {
    AGILE_ASSERT(Z.getNumColumns() + col_offset <= X_data.getNumColumns(),
                  StandardException::ExceptionSizeMismatch(
                    "Z-columns + col_offset", "X_data-columns", Z.getNumColumns() + col_offset,
                   X_data.getNumColumns()));

    AGILE_ASSERT(Z.getNumRows() + row_offset <= X_data.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "Z-rows + row_offset", "X_data-rows", Z.getNumRows() + row_offset,
                   X_data.getNumRows()));

     lowlevel::get_content(X_data.data(), X_data.getNumRows(), X_data.getNumColumns(),
                           row_offset, col_offset,
                           Z.data(), Z.getNumRows(), Z.getNumColumns());
  }

  //! \brief Extract the phase values of a GPU Matrix.
  //!
  //! \param[in] x A GPU Matrix.
  //! \param[out] y The phase values of the GPU Matrix \p x.
  template <typename TDerived1, typename TDerived2>
  inline
  void phase(const GPUMatrix<TDerived1>& X, GPUMatrix<TDerived2>& Y)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));

    lowlevel::phaseVector(X.data(), Y.data(), Y.getNumColumns()*Y.getNumRows());
  }


} // namespace agile

#endif // AGILE_GPU_MATRIX_HPP

// End of $Id: gpu_matrix.hpp 476 2011-06-16 08:54:14Z freiberger $.
