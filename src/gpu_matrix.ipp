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

// $Id: gpu_matrix.ipp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/exception.hpp"
#include "agile/gpu_config.hpp"
#include "agile/gpu_environment.hpp"
#include "agile/gpu_matrix.hpp"
#include "agile/gpu_type_traits.hpp"

#include <cuda.h>




namespace agile
{
  using namespace StandardException;


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //Multiply Matrix Scalar (implemented in CUBLAS)
  template <typename TType>
  struct multiplyMatrixScalar_Helper;

  template <>
  struct multiplyMatrixScalar_Helper<float>
  {
      static void scal(unsigned int numColumns, unsigned int numRows, const float* matrix,
                      const float* alpha, float* y)
      {
         cublasHandle_t handle = GPUEnvironment::getcublasHandle();

        if(matrix != y)
          CUBLAS_SAFE_CALL(cublasScopy (handle, numRows * numColumns, matrix, unsigned(1), y, unsigned(1)));

         CUBLAS_SAFE_CALL(cublasSscal (handle, numRows * numColumns, alpha, y, unsigned (1)));

      }
  };

  template <>
  struct multiplyMatrixScalar_Helper<double>
  {
      static void scal(unsigned int numColumns, unsigned int numRows, const double* matrix,
                      const double* alpha, double* y)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
        if(matrix != y)
          CUBLAS_SAFE_CALL(cublasDcopy (handle, numRows * numColumns, matrix, unsigned(1), y, unsigned(1)));

          CUBLAS_SAFE_CALL(cublasDscal (handle, numRows * numColumns, alpha, y, unsigned (1)));

      }
  };


  template <>
  struct multiplyMatrixScalar_Helper<cuComplex>
  {
      static void scal(unsigned int numColumns, unsigned int numRows, const cuComplex* matrix,
                      const float* alpha, cuComplex* y)
      {
         cublasHandle_t handle = GPUEnvironment::getcublasHandle();

         if(matrix != y)
          CUBLAS_SAFE_CALL(cublasCcopy (handle, numRows * numColumns, matrix, unsigned(1), y, unsigned(1)));

          CUBLAS_SAFE_CALL(cublasCsscal (handle, numRows * numColumns, alpha, y, unsigned (1)));
      }
  };

  template <>
  struct multiplyMatrixScalar_Helper<cuDoubleComplex>
  {
      static void scal(unsigned int numColumns, unsigned int numRows, const cuDoubleComplex* matrix,
                      const double* alpha, cuDoubleComplex* y)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
        if(matrix != y)
          CUBLAS_SAFE_CALL(cublasZcopy (handle, numRows * numColumns, matrix, unsigned(1), y, unsigned(1)));

          CUBLAS_SAFE_CALL(cublasZdscal (handle, numRows * numColumns, alpha, y, unsigned (1)));

      }
  };


  //! \brief Multiply a GPU matrix with a scalar (host function).
  template <typename TType1>
  void scale(const typename to_real_type<TType1>::type& alpha, const GPUMatrix<TType1>& A,
             GPUMatrix<TType1>& B)
  {
      AGILE_ASSERT(A.getNumColumns() == B.getNumColumns(),
                    ExceptionSizeMismatch("A", "B",
                                          A.getNumColumns(), B.getNumColumns()));

      AGILE_ASSERT(A.getNumRows() == B.getNumRows(),
                    ExceptionSizeMismatch("A", "B",
                                          A.getNumRows(), B.getNumRows()));

      multiplyMatrixScalar_Helper<typename substitute_gpu_complex<TType1>::cublas_type>::scal(
                         A.getNumColumns(), A.getNumRows(),
                         (const typename substitute_gpu_complex<TType1>::cublas_type*)A.data(),
                         (const typename to_real_type<TType1>::type*) (&alpha),
                         (typename substitute_gpu_complex<TType1>::cublas_type*)B.data());

  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////











  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //Multiply Matrix Vector (implemented in CUBLAS)
  template <typename TType>
  struct multiplyMatrixVector_Helper;

  template <>
  struct multiplyMatrixVector_Helper<float>
  {

      static void gemv(unsigned int numColumns, unsigned int numRows, const float* matrix,
                      const float* x, float* y)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
        float alpha=1;
        float beta=0;

          CUBLAS_SAFE_CALL(cublasSgemv (handle, CUBLAS_OP_T, numColumns, numRows, &alpha, matrix, numColumns, x, (unsigned int)(1),
                       &beta, y, (unsigned int)(1) ));

      }
  };


  template <>
  struct multiplyMatrixVector_Helper<double>
  {
      static void gemv(unsigned int numColumns, unsigned int numRows, const double* matrix,
                      const double* x, double* y)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
        double alpha=1;
        double beta=0;

          CUBLAS_SAFE_CALL(cublasDgemv (handle, CUBLAS_OP_T, numColumns, numRows, &alpha, matrix, numColumns, x, (unsigned int)(1),
                       &beta, y, (unsigned int)(1) ));

      }
  };


  template <>
  struct multiplyMatrixVector_Helper<cuComplex>
  {
      static void gemv(unsigned int numColumns, unsigned int numRows, const cuComplex* matrix,
                    const cuComplex* x, cuComplex* y)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
        cuComplex alpha={1,0};
        cuComplex beta={0,0};
          CUBLAS_SAFE_CALL(cublasCgemv (handle, CUBLAS_OP_T, numColumns, numRows, &alpha, matrix, numColumns, x, (unsigned int)(1),
                       &beta, y, (unsigned int)(1) ));

      }
  };



  template <>
  struct multiplyMatrixVector_Helper<cuDoubleComplex>
  {
      static void gemv(unsigned int numColumns, unsigned int numRows, const cuDoubleComplex* matrix,
                    const cuDoubleComplex* x, cuDoubleComplex* y)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
        cuDoubleComplex alpha={1,0};
        cuDoubleComplex beta={0,0};
          CUBLAS_SAFE_CALL(cublasZgemv (handle, CUBLAS_OP_T, numColumns, numRows, &alpha, matrix, numColumns, x, (unsigned int)(1),
                       &beta, y, (unsigned int)(1) ));

      }
  };


  //! \brief Multiply a dense GPU matrix with a GPU vector (host function).
  template <typename TType>
  void multiply(const GPUMatrix<TType>& A, const GPUVector<TType>& x,
                GPUVector<TType>& y)
  {
      AGILE_ASSERT(A.getNumColumns() == x.size(),
                    ExceptionSizeMismatch("columns of A", "x",
                                          A.getNumColumns(), x.size()));
      AGILE_ASSERT(A.getNumRows() == y.size(),
                    ExceptionSizeMismatch("rows of A", "y",
                                          A.getNumRows(), y.size()));

      multiplyMatrixVector_Helper<typename substitute_gpu_complex<TType>::cublas_type>::gemv(
                              A.getNumColumns(), A.getNumRows(),
                              (const typename substitute_gpu_complex<TType>::cublas_type*)A.data(),
                              (const typename substitute_gpu_complex<TType>::cublas_type*)x.data(),
                              (typename substitute_gpu_complex<TType>::cublas_type*)y.data());

  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////







  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //Multiply Matrix Vector Hermitian (implemented in CUBLAS + 2 VectorConjugations in CUDA)
  template <typename TType>
  struct multiplyHermMatrixVector_Helper;

  template <>
  struct multiplyHermMatrixVector_Helper<float>
  {

      static void gemv(unsigned int numColumns, unsigned int numRows, const float* matrix,
                      const float* x, float* y)
      {
        cublasHandle_t handle = agile::GPUEnvironment::getcublasHandle();
        float alpha=1;
        float beta=0;

          CUBLAS_SAFE_CALL(cublasSgemv (handle, CUBLAS_OP_N, numColumns, numRows, &alpha, matrix, numColumns, x, unsigned(1), &beta, y,
                       unsigned(1) ));

      }
  };


  template <>
  struct multiplyHermMatrixVector_Helper<double>
  {
      static void gemv(unsigned int numColumns, unsigned int numRows, const double* matrix,
                      const double* x, double* y)
      {
        cublasHandle_t handle = agile::GPUEnvironment::getcublasHandle();
        double alpha=1;
        double beta=0;

          CUBLAS_SAFE_CALL(cublasDgemv (handle, CUBLAS_OP_N, numColumns, numRows, &alpha, matrix, numColumns, x, unsigned(1), &beta, y,
                       unsigned(1) ));

      }
  };


  template <>
  struct multiplyHermMatrixVector_Helper<cuComplex>
  {
      static void gemv(unsigned int numColumns, unsigned int numRows, const cuComplex* matrix,
                    const cuComplex* x, cuComplex* y)
      {
        cublasHandle_t handle = agile::GPUEnvironment::getcublasHandle();
        cuComplex alpha={1,0};
        cuComplex beta={0,0};
          CUBLAS_SAFE_CALL(cublasCgemv (handle, CUBLAS_OP_N, numColumns, numRows, &alpha, matrix, numColumns, x, unsigned(1),
                       &beta, y, unsigned(1) ));

      }

  };



  template <>
  struct multiplyHermMatrixVector_Helper<cuDoubleComplex>
  {
      static void gemv(unsigned int numColumns, unsigned int numRows, const cuDoubleComplex* matrix,
                    const cuDoubleComplex* x, cuDoubleComplex* y)
      {
        cublasHandle_t handle = agile::GPUEnvironment::getcublasHandle();
        cuDoubleComplex alpha={1,0};
        cuDoubleComplex beta={0,0};

          CUBLAS_SAFE_CALL(cublasZgemv (handle, CUBLAS_OP_N, numColumns, numRows, &alpha, matrix, numColumns, x, unsigned(1),
                       &beta, y, unsigned(1) ));

      }
  };

  //! \brief \f$ A^H x \f$ with dense GPU matrix A and GPU vector x (host function).
  template <typename TType>
  void multiply(const GPUVector<TType>& x, const GPUMatrix<TType>& A,
                GPUVector<TType>& y)
  {
      AGILE_ASSERT(A.getNumColumns() == x.size(),
                    ExceptionSizeMismatch("columns of A", "x",
                                          A.getNumColumns(), x.size()));
      AGILE_ASSERT(A.getNumRows() == y.size(),
                    ExceptionSizeMismatch("rows of A", "y",
                                          A.getNumRows(), y.size()));

      //Create a help vector to conjugate x
      GPUVector<TType> x_conj(x.size());
      //Conjugation of vector x
      conjVector(x, x_conj);

      //A^H*x = y  is the same like   A^T*x_conj = y_conj
      multiplyHermMatrixVector_Helper<typename substitute_gpu_complex<TType>::cublas_type>::gemv(
              A.getNumColumns(), A.getNumRows(),
              (const typename substitute_gpu_complex<TType>::cublas_type*)A.data(),
              (const typename substitute_gpu_complex<TType>::cublas_type*)x_conj.data(),
              (typename substitute_gpu_complex<TType>::cublas_type*)y.data());

      //ReConjugation of vector y
      conjVector(y, y);
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////





  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //Multiply Matrix Matrix Element-wise (implemented in CUDA)
  //! \brief Element-wise matrix-matrix product (host function).
  template <typename TType>
  void multiplyElementwise(const GPUMatrix<TType>& A,
                           const GPUMatrix<TType>& B,
                           GPUMatrix<TType>& Z)
  {
      AGILE_ASSERT(A.getNumRows() == B.getNumRows(),
                    ExceptionSizeMismatch("rows of A", "rows of B",
                                          A.getNumRows(), B.getNumRows()));
      AGILE_ASSERT(A.getNumColumns() == B.getNumColumns(),
                    ExceptionSizeMismatch("columns of A", "columns of B",
                                          A.getNumColumns(), B.getNumColumns()));
      AGILE_ASSERT(A.getNumRows() == Z.getNumRows(),
                    ExceptionSizeMismatch("rows of A", "rows of Z",
                                          A.getNumRows(), Z.getNumRows()));
      AGILE_ASSERT(A.getNumColumns() == Z.getNumColumns(),
                    ExceptionSizeMismatch("columns of A", "columns of Z",
                                          A.getNumColumns(), Z.getNumColumns()));

      unsigned size = A.getNumColumns() * A.getNumRows();

      lowlevel::multiplyElementwise(A.data(), B.data(), Z.data(), size);
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////



  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //Multiply Matrix Scalar (implemented in CUBLAS)
  template <typename TType>
  struct norm2_Helper;

  template <>
  struct norm2_Helper<float>
  {

      static void norm2(unsigned int numColumns, unsigned int numRows, const float* matrix,
                        float* result)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasSnrm2 (handle, numRows * numColumns, matrix, unsigned(1),result));

      }
  };

  template <>
  struct norm2_Helper<double>
  {
      static void norm2(unsigned int numColumns, unsigned int numRows, const double* matrix,
                        double* result)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasDnrm2 (handle, numRows * numColumns, matrix, unsigned(1),result));

      }
  };


  template <>
  struct norm2_Helper<cuComplex>
  {
      static void norm2(unsigned int numColumns, unsigned int numRows, const cuComplex* matrix,
                        float* result)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasScnrm2 (handle, numRows * numColumns, matrix, unsigned(1),result));

      }
  };

  template <>
  struct norm2_Helper<cuDoubleComplex>
  {
      static void norm2(unsigned int numColumns, unsigned int numRows, const cuDoubleComplex* matrix,
                        double* result)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasDznrm2 (handle, numRows * numColumns, matrix, unsigned(1),result));
      }
  };


  //! \brief Multiply a GPU matrix with a scalar (host function).
  template <typename TType>
  typename to_real_type<TType>::type norm2(const GPUMatrix<TType>& A)
  {
      typename to_real_type<TType>::type result=0;

      norm2_Helper<typename substitute_gpu_complex<TType>::cublas_type>::norm2(
                         A.getNumColumns(), A.getNumRows(),
                         (const typename substitute_gpu_complex<TType>::cublas_type*) A.data(),
                         (typename to_real_type<TType>::type*) &result);

      return result;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////



  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //Multiply Matrix Scalar (implemented in CUBLAS)
  template <typename TType>
  struct copy_Helper;

  template <>
  struct copy_Helper<float>
  {

      static void copy(unsigned int numColumns, unsigned int numRows, const float* matrix_in,
                        float* matrix_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasScopy (handle, numRows * numColumns, matrix_in, unsigned(1), matrix_out, unsigned(1)));

      }
  };

  template <>
  struct copy_Helper<double>
  {
      static void copy(unsigned int numColumns, unsigned int numRows, const double* matrix_in,
                        double* matrix_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasDcopy (handle, numRows * numColumns, matrix_in, unsigned(1), matrix_out, unsigned(1)));

      }
  };


  template <>
  struct copy_Helper<cuComplex>
  {
      static void copy(unsigned int numColumns, unsigned int numRows, const cuComplex* matrix_in,
                         cuComplex* matrix_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasCcopy (handle, numRows * numColumns, matrix_in, unsigned(1), matrix_out, unsigned(1)));

      }
  };

  template <>
  struct copy_Helper<cuDoubleComplex>
  {
      static void copy(unsigned int numColumns, unsigned int numRows, const cuDoubleComplex* matrix_in,
                         cuDoubleComplex* matrix_out)
      {
        cublasHandle_t handle = GPUEnvironment::getcublasHandle();
          CUBLAS_SAFE_CALL(cublasZcopy (handle, numRows * numColumns, matrix_in, unsigned(1), matrix_out, unsigned(1)));
      }
  };


  //! \brief copy a GPU matrix A to Z (host function).
  template <typename TType>
  void copy(const GPUMatrix<TType>& A, GPUMatrix<TType>& Z)
  {
      AGILE_ASSERT(A.getNumColumns() == Z.getNumColumns(),
                    ExceptionSizeMismatch("A", "Z",
                                          A.getNumColumns(), Z.getNumColumns()));

      AGILE_ASSERT(A.getNumRows() == Z.getNumRows(),
                    ExceptionSizeMismatch("A", "Z",
                                          A.getNumRows(), Z.getNumRows()));

      copy_Helper<typename substitute_gpu_complex<TType>::cublas_type>::copy(
                         A.getNumColumns(), A.getNumRows(),
                         (const typename substitute_gpu_complex<TType>::cublas_type*) A.data(),
                         (typename substitute_gpu_complex<TType>::cublas_type*) Z.data());

  }
  /*
  //! \brief copy a GPU matrix A.data to data (host function).
  template <typename TType>
  void copy(const GPUMatrix<TType>& A, GPUMatrix<TType>* B)
  {
      copy_Helper<typename substitute_gpu_complex<TType>::cublas_type>::copy(
                         A.getNumColumns(), A.getNumRows(),
                         (const typename substitute_gpu_complex<TType>::cublas_type*) A.data(),
                         (typename substitute_gpu_complex<TType>::cublas_type*) (*B).data());

  }*/

  ///////////////////////////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //Multiply Transposed Matrix with Matrix (implemented in CUBLAS)
  template <typename TType>
  struct multiplyTransMatrixMatrix_Helper;

  template <>
  struct multiplyTransMatrixMatrix_Helper<float>
  {

    static void gemm(unsigned int numColumns, unsigned int numRows, const float* x_mat,
                     const float* y_mat, float* z_mat)
    {
      cublasHandle_t handle = GPUEnvironment::getcublasHandle();
      float alpha=1;
      float beta=0;

      CUBLAS_SAFE_CALL(cublasSgemm (handle, CUBLAS_OP_T, CUBLAS_OP_N, numRows, numColumns, numRows, &alpha, x_mat, numColumns, y_mat, numRows,
                   &beta, z_mat, numRows ));

    }
  };


  template <>
  struct multiplyTransMatrixMatrix_Helper<double>
  {
    static void gemm(unsigned int numColumns, unsigned int numRows, const double* x_mat,
                     const double* y_mat, double* z_mat)
    {
      cublasHandle_t handle = GPUEnvironment::getcublasHandle();
      double alpha=1;
      double beta=0;

      CUBLAS_SAFE_CALL(cublasDgemm (handle, CUBLAS_OP_T, CUBLAS_OP_N, numRows, numColumns, numRows, &alpha, x_mat, numColumns, y_mat, numRows,
                       &beta, z_mat, numRows ));

    }
  };


  template <>
  struct multiplyTransMatrixMatrix_Helper<cuComplex>
  {
    static void gemm(unsigned int numColumns, unsigned int numRows, const cuComplex* x_mat,
                     const cuComplex* y_mat, cuComplex* z_mat)
    {
      cublasHandle_t handle = GPUEnvironment::getcublasHandle();
      cuComplex alpha={1,0};
      cuComplex beta={0,0};
      CUBLAS_SAFE_CALL(cublasCgemm (handle, CUBLAS_OP_T, CUBLAS_OP_N, numRows, numColumns, numRows, &alpha, x_mat, numColumns, y_mat, numRows,
                       &beta, z_mat, numRows ));

    }
  };



  template <>
  struct multiplyTransMatrixMatrix_Helper<cuDoubleComplex>
  {
    static void gemm(unsigned int numColumns, unsigned int numRows, const cuDoubleComplex* x_mat,
                     const cuDoubleComplex* y_mat, cuDoubleComplex* z_mat)
    {
      cublasHandle_t handle = GPUEnvironment::getcublasHandle();
      cuDoubleComplex alpha={1,0};
      cuDoubleComplex beta={0,0};
      CUBLAS_SAFE_CALL(cublasZgemm (handle, CUBLAS_OP_T, CUBLAS_OP_N, numRows, numColumns, numRows, &alpha, x_mat, numColumns, y_mat, numRows,
                       &beta, z_mat, numRows ));

      }
  };


  //! \brief Multiply a dense GPU matrix with a GPU vector (host function).
  template <typename TType>
  void multiplyTransMatrixMatrix(const GPUMatrix<TType>& X, const GPUMatrix<TType>& Y,
                GPUMatrix<TType>& Z)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Z.getNumColumns()*Z.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Z", X.getNumColumns()*X.getNumRows(), Z.getNumColumns()*Z.getNumRows()));

      multiplyTransMatrixMatrix_Helper<typename substitute_gpu_complex<TType>::cublas_type>::gemm(
                              X.getNumColumns(), X.getNumRows(),
                              (const typename substitute_gpu_complex<TType>::cublas_type*) X.data(),
                              (const typename substitute_gpu_complex<TType>::cublas_type*) Y.data(),
                                    (typename substitute_gpu_complex<TType>::cublas_type*) Z.data());

  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //Dot Product of two dense GPU matrices (read as vectors). (implemented in CUBLAS)
  template <typename TType>
  struct dotProduct_Helper;

  template <>
  struct dotProduct_Helper<float>
  {

    static void dot(unsigned int numColumns, unsigned int numRows, const float* x_mat,
                     const float* y_mat, float* val)
    {
      cublasHandle_t handle = GPUEnvironment::getcublasHandle();

      CUBLAS_SAFE_CALL(cublasSdot (handle, numRows * numColumns, x_mat, unsigned(1), y_mat, unsigned(1), val));

    }
  };


  template <>
  struct dotProduct_Helper<double>
  {
    static void dot(unsigned int numColumns, unsigned int numRows, const double* x_mat,
                     const double* y_mat, double* val)
    {
      cublasHandle_t handle = GPUEnvironment::getcublasHandle();

      CUBLAS_SAFE_CALL(cublasDdot (handle, numRows * numColumns, x_mat, unsigned(1), y_mat, unsigned(1), val));

    }
  };


  template <>
  struct dotProduct_Helper<cuComplex>
  {
    static void dot(unsigned int numColumns, unsigned int numRows, const cuComplex* x_mat,
                     const cuComplex* y_mat, cuComplex* val)
    {
      cublasHandle_t handle = GPUEnvironment::getcublasHandle();

      CUBLAS_SAFE_CALL(cublasCdotu (handle, numRows * numColumns, x_mat, unsigned(1), y_mat, unsigned(1), val));

    }
  };



  template <>
  struct dotProduct_Helper<cuDoubleComplex>
  {
    static void dot(unsigned int numColumns, unsigned int numRows, const cuDoubleComplex* x_mat,
                     const cuDoubleComplex* y_mat, cuDoubleComplex* val)
    {
      cublasHandle_t handle = GPUEnvironment::getcublasHandle();

      CUBLAS_SAFE_CALL(cublasZdotu (handle, numRows * numColumns, x_mat, unsigned(1), y_mat, unsigned(1), val));

      }
  };


  //! \brief Dot Product of two dense GPU matrices (read as vectors). (host function).
  template <typename TType>
  void dotProduct(const GPUMatrix<TType>& X, const GPUMatrix<TType>& Y,
                  TType& val)
  {
    AGILE_ASSERT(X.getNumColumns()*X.getNumRows() == Y.getNumColumns()*Y.getNumRows(),
                  StandardException::ExceptionSizeMismatch(
                    "X", "Y", X.getNumColumns()*X.getNumRows(), Y.getNumColumns()*Y.getNumRows()));


      dotProduct_Helper<typename substitute_gpu_complex<TType>::cublas_type>::dot(
                              X.getNumColumns(), X.getNumRows(),
                              (const typename substitute_gpu_complex<TType>::cublas_type*) X.data(),
                              (const typename substitute_gpu_complex<TType>::cublas_type*) Y.data(),
                                    (typename substitute_gpu_complex<TType>::cublas_type*) (&val));
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////


} // namespace agile

// End of $Id: gpu_matrix.ipp 476 2011-06-16 08:54:14Z freiberger $.
