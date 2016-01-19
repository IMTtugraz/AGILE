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

// $Id: jacobi.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_JACOBI_HPP
#define AGILE_OPERATOR_JACOBI_HPP

#include "agile/gpu_vector.hpp"
#include "agile/operator/inverse_operator.hpp"
#include <vector>

namespace agile
{
  //! \brief A Jacobi preconditioner.
  //!
  //! \note THIS CLASS IS OUTDATED AND SHOULD NOT BE USED ANY LONGER. USE
  //! THE MORE GENERAL \p DiagonalSolver INSTEAD!
  template <typename TCommunicator, typename TType>
  class JacobiPreconditioner
    : public InverseOperatorExpression<
               JacobiPreconditioner<TCommunicator, TType> >
  {
    public:
      //! \brief Constructor.
      //!
      //! Constructs a diagonal preconditioner given the diagonal of the forward
      //! operator.
      //! \param[in] communicator Object for parallel communication.
      //! \param[in] diagonal The diagonal of the forward operator for which
      //! the preconditioner is constructed.
      JacobiPreconditioner(TCommunicator& communicator,
                           const std::vector<TType>& diagonal,
                           bool accumulate = true)
        : m_communicator(communicator), m_inverse_diagonal(diagonal.size()),
          m_accumulate(accumulate)
      {
        // create the inverse diagonal
        std::vector<TType> inverse(diagonal.size());
        for (unsigned counter = 0; counter < diagonal.size(); ++counter)
          inverse[counter] = TType(1) / diagonal[counter];
        m_inverse_diagonal.assignFromHost(inverse.begin(), inverse.end());
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
      GPUVector<TType> m_inverse_diagonal;

      //! If this flag is set, the result is accumulated after multiplication
      //! with the inverse of the diagonal.
      bool m_accumulate;
  };

} // namespace agile

#endif // AGILE_OPERATOR_JACOBI_HPP

// End of $Id: jacobi.hpp 476 2011-06-16 08:54:14Z freiberger $.
