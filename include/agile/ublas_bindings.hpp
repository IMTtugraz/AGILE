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

// $Id: ublas_bindings.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_UBLAS_BINDINGS_HPP
#define AGILE_UBLAS_BINDINGS_HPP

#include "boost/numeric/ublas/vector_expression.hpp"
#include "boost/numeric/ublas/matrix_expression.hpp"

namespace agile
{
  //! Short-cut for the Boost.uBLAS namespace.
  namespace ublas = boost::numeric::ublas;

  //! \brief Set a vector to a constant value.
  template <typename TConstant, typename TType>
  BOOST_UBLAS_INLINE
  void setVectorConstant(const TConstant& value, ublas::vector<TType>& x)
  {
    for (TType* x_iter = &x[0], *x_end_iter = &x[0] + x.size();
         x_iter != x_end_iter; ++x_iter)
      *x_iter = value;
  }

  //! \brief Boost.uBLAS bindings for the vector addition.
  //!
  //! \param[in] x The first vector (vector expression) to add.
  //! \param[in] y The second vector (vector expression) to add.
  //! \param[in] z The vector to store the sum of \p x and \p y.
  template <typename TExpression1, typename TExpression2, typename TVector>
  BOOST_UBLAS_INLINE
  void addVector(const ublas::vector_expression<TExpression1>& x,
                 const ublas::vector_expression<TExpression2>& y,
                 TVector& z)
  {
    z = x + y;
  }

  //! \brief Boost.uBLAS bindings for the vector subtraction.
  //!
  //! \param[in] x The first vector (vector expression).
  //! \param[in] y The second vector (vector expression) to subtract.
  //! \param[in] z The vector to store the difference \f$ x - y \f$.
  template <typename TExpression1, typename TExpression2, typename TVector>
  BOOST_UBLAS_INLINE
  void subVector(const ublas::vector_expression<TExpression1>& x,
                 const ublas::vector_expression<TExpression2>& y,
                 TVector& z)
  {
    z = x - y;
  }

  //! \brief Boost.uBLAS bindings for scaled vector addition.
  //!
  //! Multiplies the vector \p y with the scalar \p scale and adds this scaled
  //! vector to the vector \p x. The result is stored in \p z.
  //! \param[in] x First vector or vector expression.
  //! \param[in] scale The scaling factor for \p y.
  //! \param[in] y Vector (expression) to be scaled and added to \p x.
  //! \param[out] z The resultant vector (z = x + scale * y).
  template <typename TExpression1, typename TType, typename TExpression2,
            typename TVector>
  BOOST_UBLAS_INLINE
  void addScaledVector(const ublas::vector_expression<TExpression1>& x,
                       const TType& scale,
                       const ublas::vector_expression<TExpression2>& y,
                       TVector& z)
  {
    z = x + scale * y;
  }

  //! \brief Boost.uBLAS bindings for vector inversion.
  //!
  //! Compute y[i] = 1 / x[i]
  //! \param[in] x First vector or vector expression.
  //! \param[out] y The resultant vector (y = 1 / x).
  template <typename TType, typename TVector>
  BOOST_UBLAS_INLINE
  void invertVector(const ublas::vector<TType>& x, TVector& y)
  {
    for (unsigned counter = 0; counter < x.size(); ++counter)
      y[counter] = TType(1) / x[counter];
  }

  //! \brief Boost.uBLAS bindings for scaled vector subtraction.
  //!
  //! Multiplies the vector \p y with the scalar \p scale and subtracts this
  //! scaled vector from the vector \p x. The result is stored in \p z.
  //! \param[in] x First vector or vector expression.
  //! \param[in] scale The scaling factor for \p y.
  //! \param[in] y Vector (expression) to be scaled and subtracted from \p x.
  //! \param[out] z The resultant vector (z = x - scale * y).
  template <typename TExpression1, typename TType, typename TExpression2,
            typename TVector>
  BOOST_UBLAS_INLINE
  void subScaledVector(const ublas::vector_expression<TExpression1>& x,
                       const TType& scale,
                       const ublas::vector_expression<TExpression2>& y,
                       TVector& z)
  {
    z = x - scale * y;
  }

  //! \brief Boost.uBLAS bindings for multiplying a vector with a scalar.
  //!
  //! \param[in] alpha A scalar factor.
  //! \param[in] x The vector (expression) to be scaled.
  //! \param[out] y The scaled vector \f$ \alpha x \f$.
  template <typename TType, typename TExpression, typename TVector>
  BOOST_UBLAS_INLINE
  void scale(const TType& alpha,
             const ublas::vector_expression<TExpression>& x,
             TVector& y)
  {
    y = alpha * x;
  }

  //! \brief Boost.uBLAS bindings for multiplying two vectors element-wise.
  //!
  //! Multiplies each element of the vector \p x with the corresponding element
  //! of vector \p y and stores the result in the vector \p z, i.e.
  //! \f$ z[i] \leftarrow x[i] \cdot y[i] \f$.
  //! \param[in] x First vector or vector expression.
  //! \param[in] y Second vector or vector expression.
  //! \param[in] z Element-wise product of \p x and \p y.
  template <typename TExpression1, typename TExpression2, typename TVector>
  BOOST_UBLAS_INLINE
  void multiplyElementwise(const ublas::vector_expression<TExpression1>& x,
                           const ublas::vector_expression<TExpression2>& y,
                           TVector& z)
  {
    z = element_prod(x, y);
  }

  //! \brief Boost.uBLAS bindings for the bilinear form of two vectors.
  //!
  //! This function computesx \f$ sum_i x_i * y_i \f$.
  //! \param[in] x First vector or vector expression.
  //! \param[in] y Second vector or vector expression.
  //! \return The bilinear form.
  template <typename TExpression1, typename TExpression2>
  BOOST_UBLAS_INLINE
  typename ublas::vector_scalar_binary_traits<
    TExpression1, TExpression2, ublas::vector_inner_prod<
      TExpression1, TExpression2, typename ublas::promote_traits<
       typename TExpression1::value_type,
       typename TExpression2::value_type>::promote_type> >::result_type
  getBilinearForm(const ublas::vector_expression<TExpression1>& x,
                  const ublas::vector_expression<TExpression2>& y)
  {
    return inner_prod(x, y);
  }

  //! \brief Boost.uBLAS bindings for the scalar product of two vectors.
  //!
  //! This function computes \f$ sum_i conj(x_i) * y_i \f$.
  //! \param[in] x First vector or vector expression.
  //! \param[in] y Second vector or vector expression.
  //! \return The scalar product.
  template <typename TExpression1, typename TExpression2>
  BOOST_UBLAS_INLINE
  typename ublas::vector_scalar_binary_traits<
    TExpression1, TExpression2, ublas::vector_inner_prod<
      TExpression1, TExpression2, typename ublas::promote_traits<
       typename TExpression1::value_type,
       typename TExpression2::value_type>::promote_type> >::result_type
  getScalarProduct(const ublas::vector_expression<TExpression1>& x,
                   const ublas::vector_expression<TExpression2>& y)
  {
    return inner_prod(conj(x), y);
  }

  //! \brief Boost.uBLAS bindings for the matrix-vector product.
  //!
  //! \param[in] A The matrix or matrix-expression to multiply.
  //! \param[in] x The vector or vector-expression to multiply.
  //! \param[out] y The vector to store the product \f$ A x \f$.
  template <typename TExpression1, typename TExpression2, typename TVector>
  BOOST_UBLAS_INLINE
  void multiply(const ublas::matrix_expression<TExpression1>& A,
                const ublas::vector_expression<TExpression2>& x,
                TVector& y)
  {
    axpy_prod(A, x, y, true);
  }

  //! \brief Boost.uBLAS bindings for the Hermitian matrix-vector product.
  //!
  //! \param[in] x The vector or vector-expression to multiply.
  //! \param[in] A The matrix or matrix-expression to multiply.
  //! \param[out] y The vector to store the product \f$ A^H x \f$.
  template <typename TExpression1, typename TExpression2, typename TVector>
  BOOST_UBLAS_INLINE
  void multiply(const ublas::vector_expression<TExpression1>& x,
                const ublas::matrix_expression<TExpression2>& A,
                TVector& y)
  {
    axpy_prod(x, conj(A), y, true);
  }

} // namespace agile

#endif // AGILE_UBLAS_BINDINGS_HPP

// End of $Id: ublas_bindings.hpp 476 2011-06-16 08:54:14Z freiberger $.
