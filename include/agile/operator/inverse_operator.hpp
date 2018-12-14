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

// $Id: inverse_operator.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_INVERSE_OPERATOR_HPP
#define AGILE_OPERATOR_INVERSE_OPERATOR_HPP

#include "agile/operator/operator_expression.hpp"

namespace agile
{
  //! \brief Base class for an inverse operator expression.
  //!
  //! The inverse operator implements an <tt>operator()</tt> that takes
  //! a distributed vector as input and returns an accumulated vector.
  template <typename TExpression>
  class InverseOperatorExpression : public OperatorExpression<TExpression>
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

      //! \brief Apply the inverse operator.
      //!
      //! \param[in] y Input vector (distributed).
      //! \param[out] x Output vector (accumulated).
      template <typename TVectorType>
      void operator() (const TVectorType& y, TVectorType& x)
      {
        (*this)()(x, y);
      }
  };

  //! \brief The identity operator.
  template <typename TCommunicator>
  class InverseIdentity : public InverseOperatorExpression<
                                   InverseIdentity<TCommunicator> >
  {
    public:
      //! \brief The adjoint type is the identity itself.
      typedef InverseIdentity adjoint_type;

      //! \brief Constructor.
      InverseIdentity(TCommunicator& communicator)
        : m_communicator(communicator)
      {
      }

      //! \brief Inverse identity operator.
      //!
      //! \param[in] y Distributed vector.
      //! \param[out] x \p y as accumulated vector.
      template <typename TVectorType>
      void operator() (const TVectorType& y, TVectorType& x)
      {
        x = y;
        m_communicator.accumulate(x);
      }

      //! \brief Get the adjoint operator (which is the operator itself).
      adjoint_type getAdjoint() const { return *this; }

    private:
      //! Reference to the communicator for parallel communication.
      TCommunicator& m_communicator;
  };

  //! \brief The diagonal matrix as inverse operator.
  //!
  //! This class takes the elements of a diagonal matrix in vector form and
  //! implements <tt>operator()</tt> to take a distrubuted vector and to
  //! return an accumulated vector.
  //! \note This class does not compute the inverse of the diagonal or apply
  //! the inverse matrix or what so ever. It simply applies the diagonal
  //! matrix given to a vector. The term "inverse" in the name states that
  //! the input vector is distributed and the output vector is accumulated,
  //! i.e. this is just a sub-class of an \p InverseOperatorExpression.
  //! If you would like to apply the inverse of the diagonal, have a look
  //! at the class \p DiagonalSolver.
  template <typename TCommunicator, typename TVector, bool TIsAdjoint = false>
  class InverseDiagonalMatrix
    : public InverseOperatorExpression<
               InverseDiagonalMatrix<TCommunicator, TVector, TIsAdjoint> >
  {
    public:
      //! \brief The adjoint type.
      typedef InverseDiagonalMatrix<TCommunicator, TVector, !TIsAdjoint>
        adjoint_type;

      //! \brief Constructor.
      //!
      //! The constructor takes a reference to the \p communicator for network
      //! access and a reference to the diagonal elements of the matrix in
      //! \p vector.
      InverseDiagonalMatrix(TCommunicator& communicator, const TVector& vector)
        : m_communicator(communicator), m_vector(vector)
      {
      }

      //! \brief Apply the diagonal matrix.
      //!
      //! \param[in] x Distributed vector.
      //! \param[out] y The product D*x (accumulated vector).
      template <typename TVectorType>
      void operator() (const TVectorType& x, TVectorType& y)
      {
        if (TIsAdjoint)
          multiplyConjElementWise(x, m_vector, y);
        else
          multiplyElementWise(x, m_vector, y);
        m_communicator.accumulate(y);
      }

      //! \brief Get the adjoint operator.
      //!
      //! This method returns a new operator that applies the adjoint matrix
      //! to a vector.
      adjoint_type getAdjoint() const
      {
        return adjoint_type(m_communicator, m_vector);
      }

    private:
      //! Reference to the communicator for parallel communication.
      TCommunicator& m_communicator;

      //! Reference to the vector that holds the diagonal entries of the matrix.
      const TVector& m_vector;
  };

  // ============================= free functions =============================

  //! \brief Creator function for adjoints.
  template <typename TExpression>
  typename TExpression::adjoint_type adjoint(
    const InverseOperatorExpression<TExpression>& expression)
  {
    return expression.downcast().getAdjoint();
  }

} // namespace agile

#endif // AGILE_OPERATOR_INVERSE_OPERATOR_HPP

// End of $Id: inverse_operator.hpp 476 2011-06-16 08:54:14Z freiberger $.
