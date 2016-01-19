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

// $Id: gpu_vector_base.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_VECTOR_BASE_HPP
#define AGILE_GPU_VECTOR_BASE_HPP

#include "agile/exception.hpp"
#include "agile/gpu_type_traits.hpp"

namespace agile
{
  //! \brief Base class for all GPU vectors.
  //!
  //! This is the base class for all vector expressions. The minimum requirement
  //! for all derived types is that they implement:
  //! - a \p size() method returning the number of elements in the vector,
  //! - a \p data() method returning a pointer to the first vector element.
  template <typename TDerived>
  class GPUVectorBase
  {
    public:
      //! The type of the derived vector class.
      typedef TDerived derived_type;

      //! Downcast to derived type.
      inline
      const derived_type& operator() () const
      {
        return *static_cast<const derived_type*> (this);
      }

      //! Downcast to derived type.
      inline
      derived_type& operator() ()
      {
        return *static_cast<derived_type*> (this);
      }

    protected:
      //! Protected constructor.
      GPUVectorBase() {}

      //! Protected destructor.
      ~GPUVectorBase() {}
  };

  namespace lowlevel
  {
    //! \brief Add a scaled GPU vector to another vector (host function).
    //!
    //! z <- x + scale * y
    template <typename TType1, typename TType2, typename TType3>
    void addScaledVector(
      const TType1* x, const TType2& scale, const TType3* y,
      typename promote<typename promote<TType1, TType2>::type,
                       TType3>::type* z,
      unsigned size);

    //! \brief Add two vectors (host function).
    //!
    //! z <- x + y
    template <typename TType1, typename TType2>
    void addVector(const TType1* x, const TType2* y,
                   typename promote<TType1, TType2>::type* z, unsigned size);

    //! \brief Divide a constant by a vector (elementwise; host function).
    //!
    //! y[i] <- alpha / x[i]
    template <typename TType1, typename TType2>
    void divideVector(const TType1& alpha, const TType2* x,
                      typename promote<TType1, TType2>::type* y, unsigned size);

    //! Used to shift between centered, and not centered Fast Fourier Transform:
    //! Usually, the CUDA FFT considers the matrix element [0,0] to be the
    //! DC-part (center of kSpace). For visual display, a commonly used convention
    //! is to display the DC-part in the center. This is performed by function
    //! fftshift, which allows to toggle between the two versions.
    //! \param[in] M Centered/Not Centered kSpace matrix.
    //! \param[out] M Not Centered/Centered kSpace matrix, inplace computation.
    template <typename TType1>
    void fftshift(TType1* x, unsigned rows, unsigned cols);
    template <typename TType1>
    void ifftshift(TType1* x, unsigned rows, unsigned cols);


    //! \brief Compute the bilinear form of two vectors (host function).
    //!
    //! For real vectors this function will return the scalar product, for
    //! complex vector this is different from the scalar product.
    //! \return (x, y) := sum_i x[i] * y[i]
    template <typename TType1, typename TType2>
    typename promote<TType1, TType2>::type getBilinearForm(
      const TType1* x, const TType2* y, unsigned size);

    //! \brief Compute the scalar product of two vectors (host function).
    //!
    //! \return (x, y) := sum_i conj(x[i]) * y[i]
    template <typename TType1, typename TType2>
    typename promote<TType1, TType2>::type getScalarProduct(
      const TType1* x, const TType2* y, unsigned size);

    //! \brief Extract the imaginary part of a GPU vector (host function).
    //!
    //! y <- imag(x)
    template <typename TType1>
    void imag(const TType1* x,
              typename to_real_type<TType1>::type* y, unsigned size);

    //! \brief Bilinear interpolation with reshape of 1d vector (host function).
    template <typename TType1, typename TType2>
    void interpolate2d(const TType1* src, unsigned numColumns,
                       unsigned numRows, bool reshapeRowMajor,
                       const std::complex<TType2>* pos, TType1* res,
                       unsigned size);

    //! \brief Multiply conjugate GPU vector with another one element-wise (host function).
    //!
    //! z[i] <- conj(x[i]) * y[i]
    template <typename TType1, typename TType2>
    void multiplyConjElementwise(const TType1* x, const TType2* y,
                                 typename promote<TType1, TType2>::type* z,
                                 unsigned size);

    //! \brief Multiply two GPU vectors element-wise (host function).
    //!
    //! z[i] <- x[i] * y[i]
    template <typename TType1, typename TType2>
    void multiplyElementwise(const TType1* x, const TType2* y,
                             typename promote<TType1, TType2>::type* z,
                             unsigned size);

    //! \brief Divide two GPU vectors element-wise (host function).
    //!
    //! z[i] <- x[i] / y[i]
    template <typename TType1, typename TType2>
    void divideElementwise(const TType1* x, const TType2* y,
                             typename promote<TType1, TType2>::type* z,
                             unsigned size);


    //! \brief Compute the l1-norm of a GPU vector (host function).
    //!
    //! \return \f$ \sum_i \vert x_i \vert \f$
    template <typename TType1>
    typename to_real_type<TType1>::type norm1(const TType1* x, unsigned size);

    //! \brief Compute the l2-norm of a GPU vector (host function).
    //!
    //! \return \f$ \left ( \sum_i x_i^2 \right)^{1/2} \f$
    template <typename TType1>
    typename to_real_type<TType1>::type norm2(const TType1* x, unsigned size);
/*
    //! \brief Compute the maximum norm of a GPU vector (host function).
    //!
    //! \return \f$ \max_i \vert x_i \vert \f$
    template <typename TType1>
    typename to_real_type<TType1>::type normInf(const TType1* x, unsigned size);
*/
    //! \brief Extract the real part of a GPU vector (host function).
    //!
    //! y <- real(x)
    template <typename TType1>
    void real(const TType1* x,
              typename to_real_type<TType1>::type* y, unsigned size);

    //! \brief Multiply a GPU vector with a scalar (host function).
    //!
    //! y <- alpha * x
    template <typename TType1, typename TType2>
    void scale(const TType1& alpha, const TType2* x,
               typename promote<TType1, TType2>::type* y, unsigned size);

    //! \brief Set all elements of a GPU vector to a constant (host function).
    //!
    //! x[i] <- value
    template <typename TType>
    void setVectorConstant(const TType& value, TType* x, unsigned size);

    //! \brief Compute the square root of every element of a GPU vector (host function).
    //!
    //! y[i] <- sqrt(x[i])
    template <typename TType1>
    void sqrt(const TType1* x, TType1* y, unsigned size);

    //! \brief Subtract two vectors (host function).
    //!
    //! z <- x - y
    template <typename TType1, typename TType2>
    void subVector(const TType1* x, const TType2* y,
                   typename promote<TType1, TType2>::type* z, unsigned size);

    //! \brief Subtract a scaled GPU vector from another vector (host function).
    //!
    //! z <- x - scale * y
    template <typename TType1, typename TType2, typename TType3>
    void subScaledVector(
      const TType1* x, const TType2& scale, const TType3* y,
      typename promote<typename promote<TType1, TType2>::type,
                       TType3>::type* z,
      unsigned size);

    //! \brief Conjugation of one vector.
    //!
    //! \param[in] z Conjugate of x \p conj(x).
    template <typename TType1>
    void conjVector(const TType1* x, TType1* z, unsigned size);

    //! \brief Exponential of one vector.
    //!
    //! \param[in] z Exp(x).
    template <typename TType1>
    void expVector(const TType1* x, TType1* z, unsigned size);

    //! \brief Calculates the Absolute-Value of a GPU vector (host function).
    //!
    //! y <- abs(x)
    template <typename TType1>
    void absVector(const TType1* x,
                   typename to_real_type<TType1>::type* y, unsigned size);

    //! \brief Calculates the Phase-Value of a GPU vector (host function).
    //!
    //! y <- phase(x)
    template <typename TType1, typename TType2>
    void phaseVector(const TType1* x,
                   TType2* y, unsigned size);

    //////////////////////////////
    //! \brief Generate linearly spaced vector between a and b
    //!
    //! x = linspace(a,b)
    template <typename TType1>
    void linspace(TType1* x, unsigned size, float a, float b);

    //////////////////////////////
    //! \brief Generates meshgrid Vectors for input Vector x and y
    //!
    template <typename TType1>
    void meshgrid(TType1* mesh_x, TType1* mesh_y,
                  const TType1* x, unsigned x_size,const TType1* y, unsigned y_size);

    //! \brief Compute the power of alpha for every element of a GPU Vector (host function).
    //!
    //! y[i] <- (x[i])^(alpha)
    template <typename TType1, typename TType2>
    void pow(const TType1& alpha,
             const TType2* x,
                   TType2* y, unsigned size);

    //! \brief generate pattern of a vector.
    //! \brief z = abs(x)>0;   for complex value: 1+0i
    //!
    //! \param[in] x vector.
    //! \param[out] z patterned vector.
    template <typename TType1>
    void pattern(const TType1* x, typename to_real_type<TType1>::type* z, unsigned size);

    //! \brief Compute the difference between each value in a GPU vector  (host function).
    //! \brief (last value with first value)
    //!
    //! \param[in] x - The input vector.
    //! \param[in] dim - dimension to be calculated (1...rowmajor, 2...columnmajor)
    //! \param[in] x_size - x_dimension size of Matrix/Vector
    //! \param[out] y - The vector y[i] <- diff(x[i]).
    template <typename TType1>
    void diff(const unsigned dim, const unsigned x_size, const TType1* x, TType1* y, unsigned size);


    //! \brief Compute the difference between each value in a GPU vector transposed (host function).
    //! \brief (last value with first value)
    //!
    //! \param[in] x - The input vector.
    //! \param[in] dim - dimension to be calculated (1...rowmajor, 2...columnmajor)
    //! \param[in] x_size - x_dimension size of Matrix/Vector
    //! \param[out] y - The vector y[i] <- difftrans(x[i]).
    template <typename TType1>
    void difftrans(const unsigned dim, const unsigned x_size, const TType1* x, TType1* y, unsigned size);

    //! \brief Compute the difference between each value in a GPU vector 
    //! \brief considering a 3D data layout (host function).
    //!
    //! if border has to be wrapped: 
    //!
    //! dim 1: diff in first dimension: Last column with first column of the same slice
    //! dim 2: diff in second dimension: Last row with first row of the same slice
    //! dim 3: diff in third dimension: Last slice with first slice
    //!
    //! else border derivative is 0
    //!
    //! \param[in] x - The input vector.
    //! \param[in] dim - dimension to be calculated (1...row, 2...column, 3...slice)
    //! \param[in] x_size - x_dimension size of Matrix/Vector
    //! \param[in] y_size - y_dimension size of Matrix/Vector
    //! \param[in] borderWrap - true if border has to be wrapped in the required dimension 
    //! \param[out] y - The vector y[i] <- diff(x[i]).
    template <typename TType1>
    void diff3(const unsigned dim, const unsigned x_size, const unsigned y_size, const TType1* x, TType1* y, unsigned size, bool borderWrap);

    //! \brief Compute the difference between each value in a transposed GPU vector 
    //! \brief considering a 3D data layout (host function).
    //!
    //! if border has to be wrapped: 
    //!
    //! dim 1: diff in first dimension: Last column with first column of the same slice
    //! dim 2: diff in second dimension: Last row with first row of the same slice
    //! dim 3: diff in third dimension: Last slice with first slice
    //!
    //! else border derivative is 0
    //!
    //! \param[in] x - The input vector.
    //! \param[in] dim - dimension to be calculated (1...row, 2...column, 3...slice)
    //! \param[in] x_size - x_dimension size of Matrix/Vector
    //! \param[in] y_size - y_dimension size of Matrix/Vector
    //! \param[in] borderWrap - true if border has to be wrapped in the required dimension 
    //! \param[out] y - The vector y[i] <- diff(x[i]).
    template <typename TType1>
    void diff3trans(const unsigned dim, const unsigned x_size, const unsigned y_size, const TType1* x, TType1* y, unsigned size, bool borderWrap);

    //! \brief Compute the backward difference between each value in a GPU vector 
    //! \brief considering a 3D data layout (host function).
    //!
    //! if border has to be wrapped: 
    //!
    //! dim 1: diff in first dimension: Last column with first column of the same slice
    //! dim 2: diff in second dimension: Last row with first row of the same slice
    //! dim 3: diff in third dimension: Last slice with first slice
    //!
    //! else border derivative is 0
    //!
    //! \param[in] x - The input vector.
    //! \param[in] dim - dimension to be calculated (1...row, 2...column, 3...slice)
    //! \param[in] x_size - x_dimension size of Matrix/Vector
    //! \param[in] y_size - y_dimension size of Matrix/Vector
    //! \param[in] borderWrap - true if border has to be wrapped in the required dimension 
    //! \param[out] y - The vector y[i] <- diff(x[i]).
    template <typename TType1>
    void bdiff3(const unsigned dim, const unsigned x_size, const unsigned y_size, const TType1* x, TType1* y, unsigned size, bool borderWrap);

    //! \brief Compute the backward difference between each value in a transposed GPU vector 
    //! \brief considering a 3D data layout (host function).
    //!
    //! if border has to be wrapped: 
    //!
    //! dim 1: diff in first dimension: Last column with first column of the same slice
    //! dim 2: diff in second dimension: Last row with first row of the same slice
    //! dim 3: diff in third dimension: Last slice with first slice
    //!
    //! else border derivative is 0
    //!
    //! \param[in] x - The input vector.
    //! \param[in] dim - dimension to be calculated (1...row, 2...column, 3...slice)
    //! \param[in] x_size - x_dimension size of Matrix/Vector
    //! \param[in] y_size - y_dimension size of Matrix/Vector
    //! \param[in] borderWrap - true if border has to be wrapped in the required dimension 
    //! \param[out] y - The vector y[i] <- diff(x[i]).
    template <typename TType1>
    void bdiff3trans(const unsigned dim, const unsigned x_size, const unsigned y_size, const TType1* x, TType1* y, unsigned size, bool borderWrap);

    //! \brief generate max-value vector of two vector (elementwise).
    //! \brief y = max(x1,x2); checks only for real values
    //!
    //! \param[in] x1 vector.
    //! \param[in] x2 vector.
    //! \param[out] y real-value vector.
    template <typename TType1, typename TType2>
    void max(const TType1* x1, const TType2* x2,  typename promote<TType1, TType2>::type* y, unsigned size);

    //! \brief generate max-value vector of two vector (elementwise).
    //! \brief y = max(x1,x2); checks only for real values
    //!
    //! \param[in] x1 vector.
    //! \param[in] x2 scalar.
    //! \param[out] y real-value vector.
    template <typename TType1, typename TType2>
    void max(const TType1* x1, const TType2& x2,  typename promote<TType1, TType2>::type* y, unsigned size);

    //! Generate data vector with zero-padding in col-dimension
    //! \param[in] x_data input data vector.
    //! \param[in] delta_o input data vector.
    //! \param[in] delta_u input data vector.
    //! \param[in] rows x_data number of rows
    //! \param[in] cols x_data number of cols
    //! \param[in] col_o number of rows to be added on right side
    //! \param[in] col_u number of rows to be added left side
    //! \param[out] z generated vector
    template <typename TType>
    void expand_coldim(const TType* x_data, const TType* delta_o, const TType* delta_u,
                       unsigned rows, unsigned cols, unsigned col_o, unsigned col_u, TType* z);


    //! Generate data vector with zero-padding in col-dimension
    //! \param[in] x_data input data vector.
    //! \param[in] delta_o input data vector.
    //! \param[in] delta_u input data vector.
    //! \param[in] rows x_data number of rows
    //! \param[in] cols x_data number of cols
    //! \param[in] col_o number of rows to be added on right side
    //! \param[in] col_u number of rows to be added left side
    //! \param[out] z generated vector
    template <typename TType>
    void expand_rowdim(const TType* x_data, const TType* delta_o, const TType* delta_u,
                       unsigned rows, unsigned cols, unsigned row_o, unsigned row_u, TType* z);


    //! copy data from a given input
    //! \param[in] x input data.
    //! \param[in] row_offset - offset in row dimension
    //! \param[in] col_offset - offset in col dimension
    //! \param[in] z generated Vector / Matrix
    template <typename TType>
    void get_content(const TType* x_data, unsigned rows, unsigned cols,
                     unsigned row_offset, unsigned col_offset, TType* z, unsigned z_rows, unsigned z_cols);

  } // namespace lowlevel

  //! \brief Add a scaled GPU vector to another GPU vector.
  //!
  //! Multiplies the vector \p y with the scalar \p scale and adds this scaled
  //! vector to the vector \p x. The result is stored in \p z.
  //! The vectors have to have the same size when calling this function.
  //! \param[in] x First vector.
  //! \param[in] scale The scaling factor for \p y.
  //! \param[in] y Vector to be scaled and added to \p x afterwards.
  //! \param[out] z The resultant vector (z = x + scale * y).
  template <typename TDerived1, typename TScale, typename TDerived2,
            typename TDerived3>
  inline
  void addScaledVector(const GPUVectorBase<TDerived1>& x, const TScale& scale,
                       const GPUVectorBase<TDerived2>& y,
                       GPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    lowlevel::addScaledVector(x().data(), scale, y().data(), z().data(),
                              x().size());
  }

  //! \brief Add two vectors.
  //!
  //! The vectors have to be of equal size when calling this function as they
  //! won't be resized due to performance considerations.
  //! \param[in] x First vector.
  //! \param[in] y Second vector.
  //! \param[out] z Sum of the two vectors.
  template <typename TDerived1, typename TDerived2, typename TDerived3>
  inline
  void addVector(const GPUVectorBase<TDerived1>& x,
                 const GPUVectorBase<TDerived2>& y,
                 GPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    lowlevel::addVector(x().data(), y().data(), z().data(), x().size());
  }

  //! Used to shift between centered, and not centered Fast Fourier Transform:
  //! Usually, the CUDA FFT considers the matrix element [0,0] to be the
  //! DC-part (center of kSpace). For visual display, a commonly used convention
  //! is to display the DC-part in the center. This is performed by function
  //! fftshift, which allows to toggle between the two versions when col=row.
  //! \param[in] M Centered/Not Centered kSpace matrix.
  //! \param[out] M Not Centered/Centered kSpace matrix, inplace computation.
  template <typename TDerived1>
  inline
  void fftshift(GPUVectorBase<TDerived1>& x)
  {
    //lowlevel::fftshift(x().data(), x().size(), 1);
    lowlevel::fftshift(x().data(), 1, x().size());
  }

  //! Used to shift between centered, and not centered Fast Fourier Transform:
  //! Usually, the CUDA FFT considers the matrix element [0,0] to be the
  //! DC-part (center of kSpace). For visual display, a commonly used convention
  //! is to display the DC-part in the center. This is performed by function
  //! fftshift, which allows to toggle between the two versions when col=row.
  //! \param[in] M Centered/Not Centered kSpace matrix.
  //! \param[out] M Not Centered/Centered kSpace matrix, inplace computation.
  template <typename TDerived1>
  inline
  void ifftshift(GPUVectorBase<TDerived1>& x)
  {
    lowlevel::ifftshift(x().data(), 1, x().size());
  }


  //! \brief Compute the bilinear form of two vectors.
  //!
  //! This function computes x^T*y. The vectors have to be of equal size.
  //! Note that for complex vectors this function will NOT compute the
  //! scalar product.
  //! \param[in] x First vector.
  //! \param[in] y Second vector.
  //! \return The bilinear form \f$ \sum_i x[i] * y[i] \f$.
  //! \see getScalarProduct
  template <typename TDerived1, typename TDerived2>
  inline
  typename promote<typename TDerived1::value_type,
                   typename TDerived2::value_type>::type
  getBilinearForm(const GPUVectorBase<TDerived1>& x,
                  const GPUVectorBase<TDerived2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    return lowlevel::getBilinearForm(x().data(), y().data(), x().size());
  }

  //! \brief Compute the scalar product of two vectors.
  //!
  //! This function computes conj(x)*y. The vectors have to be of equal size.
  //! \param[in] x First vector.
  //! \param[in] y Second vector.
  //! \return The scalar product \f$ (x, y) := \sum_i conj(x[i]) * y[i] \f$.
  template <typename TDerived1, typename TDerived2>
  inline
  typename promote<typename TDerived1::value_type,
                   typename TDerived2::value_type>::type
  getScalarProduct(const GPUVectorBase<TDerived1>& x,
                   const GPUVectorBase<TDerived2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    return lowlevel::getScalarProduct(x().data(), y().data(), x().size());
  }

  //! \brief Extract the imaginary part of a GPU vector.
  //!
  //! \param[in] x A GPU vector.
  //! \param[out] y The imaginary part of the GPU vector \p x.
  template <typename TDerived1, typename TDerived2>
  inline
  void imag(const GPUVectorBase<TDerived1>& x, GPUVectorBase<TDerived2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    lowlevel::imag(x().data(), y().data(), x().size());
  }

  //! \brief 2D interpolation of data stored in a 1D gpu vector.
  //!
  //! This function interprets the gpu vector as a field of two dimensional data
  //! with respect to the given width \p numColumns and height \p numRows
  //! information.
  //! You have the opportunity to secify row or column major reshaping.
  //! E.g. you want to rebuild a 2D data field serialized into the 1D vector by
  //! the matlab colon ( matrix(:) ) you have to set \p reshapeRowMajor to false
  //! because Matlab uses column major ordering in that case.
  //!
  //! \param[in] x Interpolation source vector. (one dimensional)
  //! \param[in] numColumns 2D width of data stored in x (for reshaping)
  //! \param[in] numRows 2D height of data stored in x (for reshaping)
  //! \param[in] reshapeRowMajor row-major or column-major reshaping
  //! \param[in] pos Complex vector containing the interpolation positions.
  //! \param[out] res The interpolation result vector. (one dimensional)
  template <typename TDerived1, typename TDerivedComplex2>
  inline
  void interpolate2d(const GPUVectorBase<TDerived1>& x, unsigned numColumns,
                     unsigned numRows, bool reshapeRowMajor,
                     const GPUVectorBase<TDerivedComplex2>& pos,
                     GPUVectorBase<TDerived1>& res)
  {
    AGILE_ASSERT(x.size() == (numColumns * numRows),
                  StandardException::ExceptionSizeMismatch("src", "(numColumns * numRows)",
                                        x.size(), (numColumns * numRows)));
    AGILE_ASSERT(pos.size() == res.size(),
                  StandardException::ExceptionSizeMismatch("pos", "res",
                                        pos.size(), res.size()));

    lowlevel::interpolate2d(x().data(), numColumns, numRows, reshapeRowMajor,
                            pos().data(), res().data(), pos().size());
  }

  //! \brief Invert every element of a GPU vector.
  //!
  //! \param[in] x A GPU vector.
  //! \param[out] y The "inverse vector" with y[i] = 1 / x[i].
  template <typename TDerived1, typename TDerived2>
  inline
  void invertVector(const GPUVectorBase<TDerived1>& x,
                    GPUVectorBase<TDerived2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    lowlevel::divideVector(typename TDerived1::value_type(1),
                           x().data(), y().data(), x().size());
  }

  //! \brief Multiply two GPU vectors element-wise with complex conjugation of one vector.
  //!
  //! Multiplies each conjugated element of the vector \p x with the
  //! corresponding element of vector \p y and stores the result in the vector
  //! \p z, i.e.
  //! \f$ z[i] \leftarrow conj(x[i]) \cdot y[i] \f$.
  //! Its the user's responsibility to make sure that the sizes of the vectors
  //! are correct before calling this function.
  //! \param[in] x First vector.
  //! \param[in] y Second vector.
  //! \param[in] z Element-wise product of \p conj(x) and \p y.
  template <typename TDerived1, typename TDerived2, typename TDerived3>
  inline
  void multiplyConjElementwise(const GPUVectorBase<TDerived1>& x,
                               const GPUVectorBase<TDerived2>& y,
                               GPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    lowlevel::multiplyConjElementwise(x().data(), y().data(), z().data(),
                                      x().size());
  }

  //! \brief Multiply two GPU vectors element-wise.
  //!
  //! Multiplies each element of the vector \p x with the corresponding element
  //! of vector \p y and stores the result in the vector \p z, i.e.
  //! \f$ z[i] \leftarrow x[i] \cdot y[i] \f$.
  //! Its the user's responsibility to make sure that the sizes of the vectors
  //! are correct before calling this function.
  //! \param[in] x First vector.
  //! \param[in] y Second vector.
  //! \param[in] z Element-wise product of \p x and \p y.
  template <typename TDerived1, typename TDerived2, typename TDerived3>
  inline
  void multiplyElementwise(const GPUVectorBase<TDerived1>& x,
                           const GPUVectorBase<TDerived2>& y,
                           GPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    lowlevel::multiplyElementwise(x().data(), y().data(), z().data(),
                                  x().size());
  }

  //! \brief Divide two GPU vectors element-wise.
  //!
  //! Divides each element of the vector \p x with the corresponding element
  //! of vector \p y and stores the result in the vector \p z, i.e.
  //! \f$ z[i] \leftarrow x[i] / y[i] \f$.
  //! Its the user's responsibility to make sure that the sizes of the vectors
  //! are correct before calling this function.
  //! \param[in] x First vector.
  //! \param[in] y Second vector.
  //! \param[in] z Element-wise product of \p x and \p y.
  template <typename TDerived1, typename TDerived2, typename TDerived3>
  inline
  void divideElementwise(const GPUVectorBase<TDerived1>& x,
                           const GPUVectorBase<TDerived2>& y,
                           GPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    lowlevel::divideElementwise(x().data(), y().data(), z().data(),
                                  x().size());
  }

  //! \brief Compute the l1-norm of a GPU vector.
  //!
  //! \return \f$ \sum_i \vert x_i \vert \f$
  template <typename TDerived1>
  inline
  typename to_real_type<typename TDerived1::value_type>::type norm1(
    const GPUVectorBase<TDerived1>& x)
  {
    return lowlevel::norm1(x().data(), x().size());
  }

  //! \brief Compute the l1-norm of a GPU vector.
  //!
  //! \return \f$ \left ( \sum_i x_i^2 \right)^{1/2} \f$
  template <typename TDerived1>
  inline
  typename to_real_type<typename TDerived1::value_type>::type norm2(
    const GPUVectorBase<TDerived1>& x)
  {
    return lowlevel::norm2(x().data(), x().size());
  }
/*
  //! \brief Compute the maximum norm of a GPU vector.
  //!
  //! \return \f$ \max_i \vert x_i \vert \f$
  template <typename TDerived1>
  inline
  typename to_real_type<typename TDerived1::value_type>::type normInf(
    const GPUVectorBase<TDerived1>& x)
  {
    return lowlevel::normInf(x().data(), x().size());
  }
*/
  //! \brief Extract the real part of a GPU vector.
  //!
  //! \param[in] x A GPU vector.
  //! \param[out] y The real part of the GPU vector \p x.
  template <typename TDerived1, typename TDerived2>
  inline
  void real(const GPUVectorBase<TDerived1>& x, GPUVectorBase<TDerived2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    lowlevel::real(x().data(), y().data(), x().size());
  }

  //! \brief Multiply a GPU vector with a scalar.
  //!
  //! \param[in] alpha A scalar factor.
  //! \param[in] x A GPU vector.
  //! \param[out] y The scaled GPU vector alpha * x.
  template <typename TScale, typename TDerived1, typename TDerived2>
  inline
  void scale(const TScale& alpha, const GPUVectorBase<TDerived1>& x,
             GPUVectorBase<TDerived2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    lowlevel::scale(alpha, x().data(), y().data(), x().size());
  }

  //! \brief Compute the square root of every element of a GPU vector.
  //!
  //! \param[in] x The input vector.
  //! \param[out] y The vector y[i] <- sqrt(x[i]).
  template <typename TDerived1>
  inline
  void sqrt(const GPUVectorBase<TDerived1>& x, GPUVectorBase<TDerived1>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    lowlevel::sqrt(x().data(), y().data(), x().size());
  }

  //! \brief Set all elements of a GPU vector to a constant value.
  template <typename TDerived1>
  inline
  void setVectorConstant(const typename TDerived1::value_type& value,
                         GPUVectorBase<TDerived1>& x)
  {
    lowlevel::setVectorConstant(value, x().data(), x().size());
  }

  //! \brief Subtract a scaled GPU vector from another GPU vector.
  //!
  //! Multiplies the vector \p y with the scalar \p scale and subtracts this
  //! scaled vector from the vector \p x. The result is stored in \p z.
  //! The vectors have to have the same size when calling this function.
  //! \param[in] x First vector.
  //! \param[in] scale The scaling factor for \p y.
  //! \param[in] y Vector to be scaled and subtracted from \p x afterwards.
  //! \param[out] z The resultant vector (z = x - scale * y).
  template <typename TDerived1, typename TScale, typename TDerived2,
            typename TDerived3>
  inline
  void subScaledVector(const GPUVectorBase<TDerived1>& x, const TScale& scale,
                       const GPUVectorBase<TDerived2>& y,
                       GPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    lowlevel::subScaledVector(x().data(), scale, y().data(), z().data(),
                              x().size());
  }

  //! \brief Subtract two vectors.
  //!
  //! The vectors have to be of equal size when calling this function as they
  //! won't be resized due to performance considerations.
  //! \param[in] x First vector.
  //! \param[in] y Second vector.
  //! \param[out] z The difference x - y.
  template <typename TDerived1, typename TDerived2, typename TDerived3>
  inline
  void subVector(const GPUVectorBase<TDerived1>& x,
                 const GPUVectorBase<TDerived2>& y,
                 GPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    lowlevel::subVector(x().data(), y().data(), z().data(), x().size());
  }

  //////////////////////////////
  //! \brief Conjugation of one vector.
  //!
  //! Its the user's responsibility to make sure that the sizes of the vectors
  //! are correct before calling this function.
  //! \param[in] x First vector.
  //! \param[in] z Conjugate of x \p conj(x).
  template <typename TDerived1>
  inline
  void conjVector(const GPUVectorBase<TDerived1>& x,
                               GPUVectorBase<TDerived1>& z)
  {
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    lowlevel::conjVector(x().data(), z().data(), x().size());
  }

  //////////////////////////////
  //! \brief Exponential of one vector.
  //!
  //! Its the user's responsibility to make sure that the sizes of the vectors
  //! are correct before calling this function.
  //! \param[in] x First vector.
  //! \param[in] z Exp of x.
  template <typename TDerived1>
  inline
  void expVector(const GPUVectorBase<TDerived1>& x,
                               GPUVectorBase<TDerived1>& z)
  {
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    lowlevel::expVector(x().data(), z().data(), x().size());
  }

  //////////////////////////////
  //! \brief Calculates Absolute-Value for vector elements.
  //!
  //! \param[in] x First vector.
  //! \param[in] y Absolute-Value of x \p abs(x).
  template <typename TDerived1, typename TDerived2>
  inline
   void absVector(const GPUVectorBase<TDerived1>& x, GPUVectorBase<TDerived2>& y)
  {
      AGILE_ASSERT(x().size() == y().size(),
                    StandardException::ExceptionSizeMismatch(
                      "x", "y", x().size(), y().size()));

      lowlevel::absVector(x().data(), y().data() ,  x().size());

  }


  //! \brief Calculates Phase-Value for vector elements.
  //!
  //! \param[in] x First vector.
  //! \param[in] y Phase-Value of x \p phase(x).
  template <typename TDerived1, typename TDerived2>
  inline
   void phaseVector(const GPUVectorBase<TDerived1>& x, GPUVectorBase<TDerived2>& y)
  {
      AGILE_ASSERT(x().size() == y().size(),
                    StandardException::ExceptionSizeMismatch(
                      "x", "y", x().size(), y().size()));

      lowlevel::phaseVector(x().data(), y().data() ,  x().size());

  }

  //////////////////////////////
  //! \brief Generate linearly spaced vector between a and b with n numbers
  //!
  //! \param[in] a start value.
  //! \param[in] b end value.
  //! \param[in] n amount of numbers
  //! \param[out] x vector with linspace values
  template <typename TDerived1>
  inline
  void linspace(GPUVectorBase<TDerived1>& x, float a, float b)
  {
    lowlevel::linspace(x().data(), x().size(), a, b);
  }


  //! \brief Compute the power of alpha for every element of a GPU Vector.
  //!
  //! \param[in] alpha power-value.
  //! \param[in] x The input vector.
  //! \param[out] y The vector y[i] <- (x[i])^(alpha).
  template <typename TType1, typename TType2>
  inline
  void pow(const TType1& alpha, const GPUVectorBase<TType2>& x, GPUVectorBase<TType2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    lowlevel::pow(alpha, x().data(), y().data(), x().size());
  }

  //! \brief generate pattern of a vector.
  //! \brief z = abs(x)>0;
  //!
  //! \param[in] x vector.
  //! \param[out] z patterned vector
  template <typename TDerived1>
  inline
  void pattern(const GPUVectorBase<TDerived1>& x, GPUVectorBase<typename to_real_type<TDerived1>::type>& z)
  {
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    lowlevel::pattern(x().data(), z().data(), x().size());
  }

  //! \brief Compute the difference between each value in a GPU Matrix.
  //! \brief (last value with first value)
  //!
  //! \param[in] x - The input vector.
  //! \param[in] dim - dimension to be calculated (1...rowmajor, 2...columnmajor)
  //! \param[out] x - vector y[i] <- diff(x[i]).
  template <typename TDerived1>
  inline
  void diff(const unsigned dim, const GPUVectorBase<TDerived1>& x, GPUVectorBase<TDerived1>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), y().size()));


    AGILE_ASSERT(dim == 1,
                 StandardException::ExceptionMessage(
               "dim - error - dim=2 for vector not possible -"));

    lowlevel::diff(dim, x().size(), x.data(), y.data(), x().size());

  }

  //! \brief Compute the difference between each value in a GPU Matrix transposed.
  //! \brief (last value with first value)
  //!
  //! \param[in] x - The input vector.
  //! \param[in] dim - dimension to be calculated (1...rowmajor, 2...columnmajor)
  //! \param[out] y - vector y[i] <- difftrans(x[i]).
  template <typename TDerived1>
  inline
  void difftrans(const unsigned dim, const GPUVectorBase<TDerived1>& x, GPUVectorBase<TDerived1>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), y().size()));


    AGILE_ASSERT(dim == 1,
                 StandardException::ExceptionMessage(
               "dim - error - dim=2 for vector not possible -"));

    lowlevel::difftrans(dim, x().size(), x.data(), y.data(), x().size());

  }

  //! \brief generate max-value vector of two vector (elementwise).
  //! \brief y = max(x1,x2); checks abs-value
  //!
  //! \param[in] x1 vector.
  //! \param[in] x2 vector.
  //! \param[out] y real-value vector.
  template <typename TDerived1, typename TDerived2, typename TDerived3>
  inline
  void max(const GPUVectorBase<TDerived1>& x1, const GPUVectorBase<TDerived2>& x2,
              GPUVectorBase<TDerived3>& y)
  {
    AGILE_ASSERT(x1().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x1", "y", x1().size(), y().size()));
    AGILE_ASSERT(x1().size() == x2().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x1", "x2", x1().size(), x2().size()));

    lowlevel::max(x1().data(),x2().data(), y().data(), x1().size());
  }
  
  //! \brief generate max-value vector of vector (elementwise) and scalar.
  //! \brief y[k] = max(x1[k],x2); checks abs-value
  //!
  //! \param[in] x1 vector.
  //! \param[in] x2 scalar.
  //! \param[out] y real-value vector.
  template <typename TDerived1, typename TDerived2, typename TDerived3>
  inline
  void max(const GPUVectorBase<TDerived1>& x1, const TDerived2& x2,
              GPUVectorBase<TDerived3>& y)
  {
    AGILE_ASSERT(x1().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x1", "y", x1().size(), y().size()));

    lowlevel::max(x1().data(),x2, y().data(), x1().size());
  }


  //! Generate data vector with value-padding in row-dimension
  //! \param[in] x_data input data vector.
  //! \param[in] delta_o values on upper side
  //! \param[in] delta_o values on bottom side
  //! \param[in] row_o number of rows to be added on upper side
  //! \param[in] row_u number of rows to be added at the bottom
  //! \param[in] z generated vector
  template <typename TDerived1>
  inline
  void expand_rowdim(const GPUVectorBase<TDerived1>& x_data,
                     const GPUVectorBase<TDerived1>& delta_o,
                     const GPUVectorBase<TDerived1>& delta_u,
                           GPUVectorBase<TDerived1>& z)
  {
    AGILE_ASSERT(z().size() ==
                 x_data().size() +
                 delta_o().size() +
                 delta_u().size(),
                  StandardException::ExceptionSizeMismatch(
                    "z", "delta_o + x + delta_u", z().size(), x_data().size() +
                   delta_o().size() +
                   delta_u().size()));

    lowlevel::expand_rowdim(x_data().data(), delta_o().data(), delta_u().data(), 1, x_data().size(),
                            delta_o().size(), delta_u().size(), z.data());
  }

  //! Generate data vector with value-padding in col-dimension
  //! \param[in] x_data input data vector.
  //! \param[in] delta_o values on left side
  //! \param[in] delta_u values on right side
  //! \param[in] col_o number of rows to be added on left side
  //! \param[in] col_u number of rows to be added right side
  //! \param[in] z generated vector
  template <typename TDerived1>
  inline
  void expand_coldim(const GPUVectorBase<TDerived1>& x_data,
                     const GPUVectorBase<TDerived1>& delta_o,
                     const GPUVectorBase<TDerived1>& delta_u,
                           GPUVectorBase<TDerived1>& z)
  {
    AGILE_ASSERT(z().size() ==
                 x_data().size() +
                 delta_o().size() +
                 delta_u().size(),
                  StandardException::ExceptionSizeMismatch(
                    "z", "delta_o + x + delta_u", z().size(), x_data().size() +
                   delta_o().size() +
                   delta_u().size()));


    lowlevel::expand_coldim(x_data().data(), delta_o().data(), delta_u().data(), 1, x_data().size(),
                            delta_o().size(), delta_u().size(), z.data());
  }



  //! copy data from a given input vector
  //! \param[in] x input data.
  //! \param[in] row_offset - offset in row dimension
  //! \param[in] col_offset - offset in col dimension
  //! \param[in] z generated Vector
  template <typename TDerived1>
  inline
  void get_content(const GPUVectorBase<TDerived1>& x_data,
                     unsigned row_offset, unsigned col_offset,
                         GPUVectorBase<TDerived1>& z)
  {
    AGILE_ASSERT(z().size() + col_offset <= x_data().size(),
                  StandardException::ExceptionSizeMismatch(
                    "z-size + offset", "x_data size", z().size() + col_offset,
                   x_data().size()));

     lowlevel::get_content(x_data().data(), 1, x_data().size(),
                           row_offset, col_offset,
                           z().data(), 1, z().size());
  }

  ///////////////////////////////////////////////

} // namespace agile

#endif // AGILE_GPU_VECTOR_BASE_HPP

// End of $Id: gpu_vector_base.hpp 476 2011-06-16 08:54:14Z freiberger $.
