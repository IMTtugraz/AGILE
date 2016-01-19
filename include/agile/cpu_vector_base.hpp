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

// $Id: cpu_vector_base.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_CPU_VECTOR_BASE_HPP
#define AGILE_CPU_VECTOR_BASE_HPP

#include "agile/exception.hpp"
#include "agile/gpu_type_traits.hpp"
#include "agile/gpu_complex.hpp"

#include <complex>

namespace agile
{
  //! \brief Base class for all CPU vectors.
  //!
  //! This is the base class for all vector expressions. The minimum requirement
  //! for all derived types is that they implement:
  //! - a \p size() method returning the number of elements in the vector,
  //! - a \p data() method returning a pointer to the first vector element.
  template <typename TDerived>
  class CPUVectorBase
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
      CPUVectorBase() {}

      //! Protected destructor.
      ~CPUVectorBase() {}
  };

  //! \brief Add a scaled CPU vector to another CPU vector.
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
  void addScaledVector(const CPUVectorBase<TDerived1>& x, const TScale& scale,
                       const CPUVectorBase<TDerived2>& y,
                       CPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    const typename TDerived2::value_type* y_ptr = y().data();
    typename TDerived3::value_type* z_ptr = z().data();
    while (x_ptr != x_ptr_end)
      *z_ptr++ = *x_ptr++ + scale * *y_ptr++;
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
  void addVector(const CPUVectorBase<TDerived1>& x,
                 const CPUVectorBase<TDerived2>& y,
                 CPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    const typename TDerived2::value_type* y_ptr = y().data();
    typename TDerived3::value_type* z_ptr = z().data();
    while (x_ptr != x_ptr_end)
      *z_ptr++ = *x_ptr++ + *y_ptr++;
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
  getScalarProduct(const CPUVectorBase<TDerived1>& x,
                   const CPUVectorBase<TDerived2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    typename promote<typename TDerived1::value_type,
                     typename TDerived2::value_type>::type sum(0);

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    const typename TDerived2::value_type* y_ptr = y().data();
    while (x_ptr != x_ptr_end)
      sum += conj(*x_ptr++) * *y_ptr++;

    return sum;
  }

  //! \brief Extract the imaginary part of a CPU vector.
  //!
  //! \param[in] x A CPU vector.
  //! \param[out] y The imaginary part of the CPU vector \p x.
  template <typename TDerived1, typename TDerived2>
  inline
  void imag(const CPUVectorBase<TDerived1>& x, CPUVectorBase<TDerived2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    typename TDerived2::value_type* y_ptr = y().data();
    while (x_ptr != x_ptr_end)
      *y_ptr++ = agile::imag(*x_ptr++);
  }

#if 0
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
  void interpolate2d(const CPUVectorBase<TDerived1>& x, unsigned numColumns,
                     unsigned numRows, bool reshapeRowMajor,
                     const CPUVectorBase<TDerivedComplex2>& pos,
                     CPUVectorBase<TDerived1>& res)
  {
    AGILE_ASSERT(src.size() == (numColumns * numRows),
                  ExceptionSizeMismatch("src", "(numColumns * numRows)",
                                        src.size(), (numColumns * numRows)));
    AGILE_ASSERT(pos.size() == res.size(),
                  ExceptionSizeMismatch("pos", "res",
                                        pos.size(), res.size()));

    lowlevelCPU::interpolate2d(x().data(), numColumns, numRows, reshapeRowMajor,
                            pos().data(), res().data(), pos().size());
  }
#endif

  //! \brief Invert every element of a CPU vector.
  //!
  //! \param[in] x A CPU vector.
  //! \param[out] y The "inverse vector" with y[i] = 1 / x[i].
  template <typename TDerived1, typename TDerived2>
  inline
  void invertVector(const CPUVectorBase<TDerived1>& x,
                    CPUVectorBase<TDerived2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    typename TDerived2::value_type* y_ptr = y().data();
    while (x_ptr != x_ptr_end)
      *y_ptr++ = TDerived1::value_type(1) / *y_ptr++;
  }

  //! \brief Multiply two CPU vectors element-wise with complex conjugation of one vector.
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
  void multiplyConjElementwise(const CPUVectorBase<TDerived1>& x,
                               const CPUVectorBase<TDerived2>& y,
                               CPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    const typename TDerived2::value_type* y_ptr = y().data();
    typename TDerived3::value_type* z_ptr = z().data();
    while (x_ptr != x_ptr_end)
      *z_ptr++ = conj(*x_ptr++) * *y_ptr++;
  }

  //! \brief Multiply two CPU vectors element-wise.
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
  void multiplyElementwise(const CPUVectorBase<TDerived1>& x,
                           const CPUVectorBase<TDerived2>& y,
                           CPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    const typename TDerived2::value_type* y_ptr = y().data();
    typename TDerived3::value_type* z_ptr = z().data();
    while (x_ptr != x_ptr_end)
      *z_ptr++ = *x_ptr++ * *y_ptr++;
  }
/*
  //! \brief Compute the l1-norm of a CPU vector.
  //!
  //! \return \f$ \sum_i \vert x_i \vert \f$
  template <typename TDerived1>
  inline
  typename to_real_type<typename TDerived1::value_type>::type norm1(
    const CPUVectorBase<TDerived1>& x)
  {
    return lowlevelCPU::norm1(x().data(), x().size());
  }
*/
  //! \brief Compute the l1-norm of a CPU vector.
  //!
  //! \return \f$ \left ( \sum_i x_i^2 \right)^{1/2} \f$
  template <typename TDerived1>
  inline
  typename to_real_type<typename TDerived1::value_type>::type norm2(
    const CPUVectorBase<TDerived1>& x)
  {
    typename TDerived1::value_type sum_of_squares(0);

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    while (x_ptr != x_ptr_end)
      sum_of_squares += agile::norm(*x_ptr++);

    return std::sqrt(sum_of_squares);
  }
/*
  //! \brief Compute the maximum norm of a CPU vector.
  //!
  //! \return \f$ \max_i \vert x_i \vert \f$
  template <typename TDerived1>
  inline
  typename to_real_type<typename TDerived1::value_type>::type normInf(
    const CPUVectorBase<TDerived1>& x)
  {
    return lowlevelCPU::normInf(x().data(), x().size());
  }
*/
  //! \brief Extract the real part of a CPU vector.
  //!
  //! \param[in] x A CPU vector.
  //! \param[out] y The real part of the CPU vector \p x.
  template <typename TDerived1, typename TDerived2>
  inline
  void real(const CPUVectorBase<TDerived1>& x, CPUVectorBase<TDerived2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    typename TDerived2::value_type* y_ptr = y().data();
    while (x_ptr != x_ptr_end)
      *y_ptr++ = agile::real(*x_ptr++);
  }

  //! \brief Multiply a CPU vector with a scalar.
  //!
  //! \param[in] alpha A scalar factor.
  //! \param[in] x A CPU vector.
  //! \param[out] y The scaled CPU vector alpha * x.
  template <typename TScale, typename TDerived1, typename TDerived2>
  inline
  void scale(const TScale& alpha, const CPUVectorBase<TDerived1>& x,
             CPUVectorBase<TDerived2>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    typename TDerived2::value_type* y_ptr = y().data();
    while (x_ptr != x_ptr_end)
      *y_ptr++ = alpha * *x_ptr++;
  }

  //! \brief Compute the square root of every element of a CPU vector.
  //!
  //! \param[in] x The input vector.
  //! \param[out] y The vector y[i] <- sqrt(x[i]).
  template <typename TDerived1>
  inline
  void sqrt(const CPUVectorBase<TDerived1>& x, CPUVectorBase<TDerived1>& y)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    typename TDerived1::value_type* y_ptr = y().data();
    while (x_ptr != x_ptr_end)
      *y_ptr++ = std::sqrt(*x_ptr++);
  }

  //! \brief Subtract a scaled CPU vector from another CPU vector.
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
  void subScaledVector(const CPUVectorBase<TDerived1>& x, const TScale& scale,
                       const CPUVectorBase<TDerived2>& y,
                       CPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    const typename TDerived2::value_type* y_ptr = y().data();
    typename TDerived3::value_type* z_ptr = z().data();
    while (x_ptr != x_ptr_end)
      *z_ptr++ = *x_ptr++ - scale * *y_ptr++;
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
  void subVector(const CPUVectorBase<TDerived1>& x,
                 const CPUVectorBase<TDerived2>& y,
                 CPUVectorBase<TDerived3>& z)
  {
    AGILE_ASSERT(x().size() == y().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "y", x().size(), y().size()));
    AGILE_ASSERT(x().size() == z().size(),
                  StandardException::ExceptionSizeMismatch(
                    "x", "z", x().size(), z().size()));

    const typename TDerived1::value_type* x_ptr = x().data();
    const typename TDerived1::value_type* x_ptr_end = x().data() + x().size();
    const typename TDerived2::value_type* y_ptr = y().data();
    typename TDerived3::value_type* z_ptr = z().data();
    while (x_ptr != x_ptr_end)
      *z_ptr++ = *x_ptr++ - *y_ptr++;
  }

} // namespace agile

#endif // AGILE_CPU_VECTOR_BASE_HPP

// End of $Id: cpu_vector_base.hpp 476 2011-06-16 08:54:14Z freiberger $.
