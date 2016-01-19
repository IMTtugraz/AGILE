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

// $Id: forward_operator.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_FORWARD_OPERATOR_HPP
#define AGILE_OPERATOR_FORWARD_OPERATOR_HPP

#include "agile/operator/operator_expression.hpp"

namespace agile
{
  //! \brief Base class for forward operator expressions.
  //!
  //! The parenthesis-operator takes an accumulated vector and returns a
  //! distributed vector. This can be achieved with a distributed matrix,
  //! for example.
  template <typename TExpression>
  class ForwardOperatorExpression : public OperatorExpression<TExpression>
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

      //! \brief Apply the forward operator.
      //!
      //! \param[in] x Input vector (accumulated).
      //! \param[out] y Output vector (distributed).
      template <typename TVectorType>
      void operator() (const TVectorType& x, TVectorType& y)
      {
        (*this)()(x, y);
      }
  };

  //! \brief The identity operator.
  template <typename TCommunicator>
  class ForwardIdentity : public ForwardOperatorExpression<
                                   ForwardIdentity<TCommunicator> >
  {
    public:
      //! \brief The adjoint type is the identity itself.
      typedef ForwardIdentity adjoint_type;

      //! \brief Constructor.
      ForwardIdentity(TCommunicator& communicator)
        : m_communicator(communicator)
      {
      }

      //! \brief Forward identity operator.
      //!
      //! \param[in] x Accumulated vector.
      //! \param[out] y \p x as distributed vector.
      template <typename TVectorType>
      void operator() (const TVectorType& x, TVectorType& y)
      {
        y = x;
        m_communicator.distribute(y);
      }

      //! \brief Get the adjoint operator (which is the operator itself).
      adjoint_type getAdjoint() const { return *this; }

    private:
      //! Reference to the communicator for parallel communication.
      TCommunicator& m_communicator;
  };

  //! \brief A matrix as forward expression.
  //!
  //! This class takes as template arguments:
  //! - \p TCommunicator The communicator type for parallel communication.
  //! - \p TMatrix The matrix type which shall be wrapped by this operator.
  //! - \p TIsAdjoint A boolean value. If it is \p false, the operator will
  //!   carry out the normal matrix-vector product. If the flag is \p true,
  //!   the operator will multiply the input vector with the hermitian matrix.
  //! - \p TIsDistributed A boolean specifying whether the matrix is already
  //!   distributed or not. If this flag is set to \p false, the resultant
  //!   vector after the multiplication will be distributed by hand.
  template <typename TCommunicator, typename TMatrix,
            bool TIsAdjoint = false, bool TIsDistributed = true>
  class ForwardMatrix
    : public ForwardOperatorExpression<
               ForwardMatrix<TCommunicator, TMatrix, TIsAdjoint,
               TIsDistributed> >
  {
    public:
      //! \brief The adjoint type.
      typedef ForwardMatrix<TCommunicator, TMatrix, !TIsAdjoint,
                            TIsDistributed> adjoint_type;

      //! \brief Constructor.
      ForwardMatrix(TCommunicator& communicator, const TMatrix& matrix)
        : m_communicator(communicator), m_matrix(matrix)
      {
      }

      //! \brief Copy constructor.
      ForwardMatrix(const ForwardMatrix& other)
        : m_communicator(other.m_communicator), m_matrix(other.m_matrix)
      {
      }

      //! \brief Apply the matrix.
      //!
      //! \param[in] x Accumulated vector.
      //! \param[out] y The product M*x (distributed vector).
      template <typename TVectorType>
      void operator() (const TVectorType& x, TVectorType& y)
      {
        // depending on the value of the adjoint flag, either A^H*x or A*x
        // is computed
        if (TIsAdjoint)
          multiply(x, m_matrix, y); // y <-- A^H * x
        else
          multiply(m_matrix, x, y); // y <-- A * x

        // if the matrix is not distributed, we have to distribute the vector
        // afterwards
        if (!TIsDistributed)
          m_communicator.distribute(y);
      }

      //! \brief Get the adjoint operator.
      //!
      //! This method returns a new operator that applies the adjoint matrix
      //! to a vector.
      adjoint_type getAdjoint() const
      {
        return adjoint_type(m_communicator, m_matrix);
      }

    private:
      //! Reference to the communicator for parallel communication.
      TCommunicator& m_communicator;

      //! Reference to the matrix.
      const TMatrix& m_matrix;
  };

  //! \brief A matrix and its Hermitian as forward expression.
  //!
  //! This class is an extention of the normal \p ForwardMatrixOperator. Its
  //! constructor takes references to two matrices, one in the original form
  //! and one in the Hermitian form. The two matrix types can be different
  //! so that it is possible to have an optimal storage layout for both the
  //!  normal and the Hermitian product.
  //! This class takes as template arguments:
  //! - \p TCommunicator The communicator type for parallel communication.
  //! - \p TMatrix The matrix type which shall be wrapped by this operator.
  //! - \p TAdjointMatrix The adjoint matrix type.
  //! - \p TIsAdjoint A boolean value. If it is \p false, the operator will
  //!   carry out the normal matrix-vector product. If the flag is \p true,
  //!   the operator will multiply the input vector with the Hermitian matrix.
  //! - \p TIsDistributed A boolean specifying whether both matrices are already
  //!   distributed or not. If this flag is set to \p false, the resultant
  //!   vector after the multiplication will be distributed by hand.
  template <typename TCommunicator, typename TMatrix, typename TAdjointMatrix,
            bool TIsAdjoint = false, bool TIsDistributed = true>
  class ForwardMatrixWithAdjoint
    : public ForwardOperatorExpression<
               ForwardMatrixWithAdjoint<TCommunicator, TMatrix, TAdjointMatrix,
                                        TIsAdjoint, TIsDistributed> >
  {
    public:
      //! \brief The adjoint type.
      typedef ForwardMatrixWithAdjoint<TCommunicator, TMatrix, TAdjointMatrix,
                !TIsAdjoint, TIsDistributed> adjoint_type;

      //! \brief Constructor.
      ForwardMatrixWithAdjoint(TCommunicator& communicator,
                               const TMatrix& matrix,
                               const TAdjointMatrix& adjoint_matrix)
        : m_communicator(communicator), m_matrix(matrix),
          m_adjoint_matrix(adjoint_matrix)
      {
      }

      //! \brief Apply the transposed matrix.
      //!
      //! \param[in] x Accumulated vector.
      //! \param[out] y The product M*x (distributed vector).
      template <typename TVectorType>
      void operator() (const TVectorType& x, TVectorType& y)
      {
        if (TIsAdjoint)
          multiply(x, m_adjoint_matrix, y);
        else
          multiply(m_matrix, x, y);
        // if the matrices are not distributed, we have to distribute the vector
        // afterwards
        if (!TIsDistributed)
          m_communicator.distribute(y);
      }

      //! \brief Get the adjoint operator.
      //!
      //! This method returns a new operator that applies the adjoint matrix
      //! to a vector.
      adjoint_type getAdjoint() const
      {
        return adjoint_type(m_communicator, m_matrix, m_adjoint_matrix);
      }

    private:
      //! Reference to the communicator for parallel communication.
      TCommunicator& m_communicator;

      //! Reference to the matrix.
      const TMatrix& m_matrix;

      //! Reference to the adjoint matrix.
      const TAdjointMatrix& m_adjoint_matrix;
  };

  //! \brief Forward operator defined by a diagonal matrix.
  //!
  //! This class defines a diagonal matrix as forward operator, that is the
  //! <tt>operator()</tt> takes an accumulated vector as input and returns
  //! a distributed vector. The diagonal matrix is stored as a vector (i.e.
  //! only the diagonal elements are stored).
  //! The class is templated for the communicator type and also the vector type.
  template <typename TCommunicator, typename TVector, bool TIsAdjoint = false>
  class ForwardDiagonalMatrix
    : public ForwardOperatorExpression<
               ForwardDiagonalMatrix<TCommunicator, TVector, TIsAdjoint> >
  {
    public:
      //! \brief The adjoint type.
      typedef ForwardDiagonalMatrix<TCommunicator, TVector, !TIsAdjoint>
        adjoint_type;

      //! \brief Constructor.
      //!
      //! The constructor takes a reference to the \p communicator for network
      //! access and a reference to the diagonal elements of the matrix in
      //! \p vector.
      ForwardDiagonalMatrix(TCommunicator& communicator, const TVector& vector)
        : m_communicator(communicator), m_vector(vector)
      {
      }

      //! \brief Apply the diagonal matrix.
      //!
      //! \param[in] x Accumulated vector.
      //! \param[out] y The product D*x (distributed vector).
      template <typename TVectorType>
      void operator() (const TVectorType& x, TVectorType& y)
      {
        if (TIsAdjoint)
          multiplyConjElementwise(x, m_vector, y);
        else
          multiplyElementWise(x, m_vector, y);
        m_communicator.distribute(y);
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
    const ForwardOperatorExpression<TExpression>& expression)
  {
    return expression.downcast().getAdjoint();
  }

  //! Free function to ease the creation of a forward identity operator.
  //!
  //! \param[in] com The communicator for parallel communication.
  template <typename TCommunicator>
  ForwardIdentity<TCommunicator> makeForwardIdentityOperator(
    TCommunicator& com)
  {
    return ForwardIdentity<TCommunicator>(com);
  }

  //! Free function to ease the creation of a forward matrix operator.
  //!
  //! \param[in] com The communicator for parallel communication.
  //! \param[in] matrix The matrix to be wrapped by the forward operator.
  template <typename TCommunicator, typename TMatrix>
  ForwardMatrix<TCommunicator, TMatrix>
    makeForwardMatrixOperator(TCommunicator& com, const TMatrix& matrix)
  {
    return ForwardMatrix<TCommunicator, TMatrix>(com, matrix);
  }

} // namespace agile

#endif // AGILE_OPERATOR_FORWARD_OPERATOR_HPP

// End of $Id: forward_operator.hpp 476 2011-06-16 08:54:14Z freiberger $.
