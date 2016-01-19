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

// $Id: diagonal_solver.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_DIAGONAL_SOLVER_HPP
#define AGILE_OPERATOR_DIAGONAL_SOLVER_HPP

#include "agile/gpu_vector_base.hpp"
#include "agile/operator/inverse_operator.hpp"

namespace agile
{
  //! \brief A diagonal matrix solver.
  //!
  //! This class implements a diagonal solver. It takes the elements of a
  //! diagonal matrix in vector form. The <tt>operator()</tt> method applies
  //! the inverse of the diagonal matrix to a vector.
  //! As this method sub-classes \p InverseOperatorExpression the input to
  //! the parenthesis-operator has to be a distributed vector and the output
  //! will be an accumulated vector.
  template <typename TCommunicator, typename TVector>
  class DiagonalSolver
    : public InverseOperatorExpression<
               DiagonalSolver<TCommunicator, TVector> >
  {
    public:
      //! \brief Constructor.
      //!
      //! Constructs a diagonal matrix solver.
      //! \param[in] communicator Object for parallel communication.
      //! \param[in] diagonal The diagonal of the matrix in vector format.
      //! This diagonal will be inverted and applied to the vector
      //! supplied to the <tt>operator()</tt> method.
      //! \param[in] accumulate If this flag is set, the result of the
      //! multiplication is accumulated using \p communicator before the
      //! product is returned to the caller. The flag defaults to \p true.
      DiagonalSolver(TCommunicator& communicator, const TVector& diagonal,
                     bool accumulate = true)
        : m_communicator(communicator), m_inverse_diagonal(diagonal.size()),
          m_accumulate(accumulate)
      {
        // create the inverse diagonal
        invertVector(diagonal, m_inverse_diagonal);
      }

      //! \brief Preconditioning operation.
      //!
      //! \param[in] y Distributed vector.
      //! \param[out] x \f$ \mbox{diag}(d)^{-1} * y \f$ (accumulated vector).
      template <typename TVectorType>
      void operator() (const TVectorType& y, TVectorType& x)
      {
        multiplyElementwise(y, m_inverse_diagonal, x);
        if (m_accumulate)
          m_communicator.accumulate(x);
      }

    private:
      //! Reference to the communicator for parallel communication.
      TCommunicator& m_communicator;

      //! The inverse of the diagonal.
      TVector m_inverse_diagonal;

      //! If this flag is set, the result is accumulated after multiplication
      //! with the inverse of the diagonal.
      bool m_accumulate;
  };

} // namespace agile

#endif // AGILE_OPERATOR_DIAGONAL_SOLVER_HPP

// End of $Id: diagonal_solver.hpp 476 2011-06-16 08:54:14Z freiberger $.
