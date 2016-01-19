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

// $Id: binary_measure_operator.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_BINARY_MEASURE_OPERATOR_HPP
#define AGILE_OPERATOR_BINARY_MEASURE_OPERATOR_HPP

#include "agile/gpu_config.hpp"
#include "agile/operator/operator_expression.hpp"

namespace agile
{
  template <typename TExpression>
  class BinaryMeasureOperatorExpression : public OperatorExpression<TExpression>
  {
    public:
      typedef TExpression expression_type;

      //! \brief Downcast to derived type.
      inline
      const expression_type& operator() () const
      {
        return *static_cast<const expression_type*> (this);
      }

      //! \brief Downcast to derived type.
      inline
      expression_type& operator() ()
      {
        return *static_cast<expression_type*> (this);
      }

      //! \brief Apply the binary operator.
      //!
      //! \param[in] x First input vector (accumulated).
      //! \param[in] y Second input vector (distributed).
      //! \return The result of the binary operator.
      template <typename TVectorType>
      typename TVectorType::value_type operator() (
        const TVectorType& x, const TVectorType& y)
      {
        return (*this)()(y, x);
      }
  };

  //! This operator computes the scalar product of two vectors which is
  //! //! \f$ sum_i conj(x_i) * y_i \f$.
  template <typename TCommunicator>
  class ScalarProductMeasure : public BinaryMeasureOperatorExpression<
                                        ScalarProductMeasure<TCommunicator> >
  {
    public:
      ScalarProductMeasure(TCommunicator& communicator)
        : m_communicator(communicator)
      {
      }

      //! \brief Uses the scalar product as a measure.
      //!
      //! \param[in] x Accumulated vector.
      //! \param[in] y Distributed vector.
      //! \return The scalar product (x, y).
      template <typename TVectorType>
      typename TVectorType::value_type operator() (const TVectorType& x,
                                                   const TVectorType& y)
      {
        typename TVectorType::value_type partial_scalar_product
          = getScalarProduct(x, y);
        // collect the partial products from all processes
        m_communicator.collect(partial_scalar_product);

        return partial_scalar_product;
      }

    private:
      //! Reference to the communicator for parallel communication.
      TCommunicator& m_communicator;
  };

  //! This operator computes the bilinear form of two vectors which is
  //! \f$ sum_i x_i * y_i \f$. For real vectors this is equal to the scalar
  //! product. However, for complex vectors the bilinear form does not
  //! define a scalar product.
  template <typename TCommunicator>
  class BilinearFormMeasure : public BinaryMeasureOperatorExpression<
                                       BilinearFormMeasure<TCommunicator> >
  {
    public:
      BilinearFormMeasure(TCommunicator& communicator)
        : m_communicator(communicator)
      {
      }

      //! \brief Uses the bilinear form as a measure.
      //!
      //! \param[in] x Accumulated vector.
      //! \param[in] y Distributed vector.
      //! \return The bilinear form \f$ sum_i x_i * y_i \f$
      template <typename TVectorType>
      typename TVectorType::value_type operator() (const TVectorType& x,
                                                   const TVectorType& y)
      {
        typename TVectorType::value_type partial_bilinear_form
          = getBilinearForm(x, y);
        // collect the partial products from all processes
        m_communicator.collect(partial_bilinear_form);

        return partial_bilinear_form;
      }

    private:
      //! Reference to the communicator for parallel communication.
      TCommunicator& m_communicator;
  };

} // namespace agile

#endif // AGILE_OPERATOR_BINARY_MEASURE_OPERATOR_HPP

// End of $Id: binary_measure_operator.hpp 476 2011-06-16 08:54:14Z freiberger $.
