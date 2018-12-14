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

// $Id: cg.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_CG_HPP
#define AGILE_OPERATOR_CG_HPP

#include "agile/gpu_config.hpp"
#include "agile/operator/forward_operator.hpp"
#include "agile/operator/inverse_operator.hpp"
#include "agile/operator/binary_measure_operator.hpp"

namespace agile
{
  //! \brief Conjugate Gradient method.
  template <typename TCommunicator, typename TForward, typename TBinaryMeasure>
  class ConjugateGradient
    : public InverseOperatorExpression<
               ConjugateGradient<TCommunicator, TForward, TBinaryMeasure> >
  {
    public:
      //! \brief Constructor.
      //!
      //! \param[in] communicator Object for parallel communication.
      //! \param[in] forward The forward operator.
      //! \param[in] measure A measurement operator.
      //! \param[in] rel_tolerance The relative tolerance
      //! \f$ \varepsilon_{rel} \f$ serving as 1st stopping criterion for the
      //! algorithm.
      //! \param[in] abs_tolerance The absolute tolerance
      //! \f$ \varepsilon_{abs} \f$ (2nd stopping criterion).
      //! \param[in] max_iterations The maximum number of iterations for the
      //! inversion after which the algorithm gives up.
      ConjugateGradient(TCommunicator& communicator,
                        ForwardOperatorExpression<TForward>& forward,
                        BinaryMeasureOperatorExpression<TBinaryMeasure>&
                          measure,
                        const double& rel_tolerance,
                        const double& abs_tolerance, unsigned max_iterations)
        : m_communicator(communicator), m_forward(forward),
          m_binary_measure(measure),
          m_rel_tolerance(rel_tolerance), m_abs_tolerance(abs_tolerance),
          m_max_iterations(max_iterations), m_convergence(false),
          m_iteration(0), m_rho_0(0), m_rho_k(0)
      {
      }

      //! \brief Apply the CG method to the vector x.
      //!
      //! This operator solves the problem F(x) = y for x.
      //! The right-hand side \p y has to be a distributed vector and the
      //! solution \p x will be an accumulated one. The initial solution
      //! is supplied via \p x as well and has to be accumulated, too.
      //!
      //! The "exactness" of the solution is measured using the measurement
      //! operator \p measure provided to the constructor and is denoted
      //! as \f$ (\cdot, \cdot) \f$. The residual norm in the k-th iteration is
      //! given by \f$ \varrho_k := (r_k, r_k) \f$, where
      //! \f$ r_k = y - F(x_k) \f$ is the residual
      //! and \f$ F \f$ the forward operator.
      //! The algorithm stops in iteration \p k when
      //! - the maximum amount of iterations is reached (\f$ k >= k_{max} \f$).
      //! - the residual is less than the absolute tolerance, i.e.
      //!   \f$ \varrho_k <= \varepsilon_{abs} \f$, or
      //! - the residual has decreased significantly such that the ratio
      //!   \f$ \varrho_k / \varrho_0 <= \varepsilon_{rel} \f$.
      //!
      //! \param[in] y The right-hand side (distributed).
      //! \param[in,out] x Provides the initial guess for the CG algorithm
      //! and returns the final solution (accumulated).
      template <typename TVectorType>
      void operator() (const TVectorType& y, TVectorType& x)
      {
        typedef typename TVectorType::value_type value_type;

        // assume that the algorithm does not converge (negative thinking ;-) )
        m_convergence = false;

        // temporary vectors
        TVectorType p_k(y.size(), 0);
        TVectorType v_k(y.size(), 0);

        // the residual y - F(x_k); this vector is always distributed
        TVectorType residual_k(y);
        // compute the initial residual r_0 = y - F(x_0)
        m_forward(x, v_k);
        subVector(residual_k, v_k, residual_k);

        // v_0 = r_0 but accumulated
        v_k = residual_k;
        m_communicator.accumulate(v_k);

        // rho_0 = (v_0, r_0); v_0 is accumulated, r_0 is distributed
        value_type rho_0 = m_binary_measure(v_k, residual_k);
        value_type rho_k = rho_0;
        // we also store the initial and current residual as double values as
        // they can be checked by the calling program afterwards
        m_rho_0 = std::abs(rho_0);
        m_rho_k = std::abs(rho_k);

        // if the solution is already exact enough (e.g. because both the rhs
        // and the initial guess are zero vectors, we return immediately)
        if (m_rho_0 <= m_abs_tolerance)
        {
          m_convergence = true;
          m_iteration = 0;
          return;
        }

        // p_0 = v_0; p_0 is accumulated because v_0 is accumulated
        p_k = v_k;

        for (m_iteration = 0; m_iteration < m_max_iterations; ++m_iteration)
        {
          // v_k = F(p_k); v_k will be distributed afterwards
          m_forward(p_k, v_k);

          // sigma_k = (v_k, p_k); p_k is accumulated
          value_type sigma_k = m_binary_measure(p_k, v_k);

          // alpha_k = rho_k / sigma_k
          value_type alpha_k = rho_k / sigma_k;

          // update the solution and the residual
          // x_{k+1} = x_k + alpha_k * p_k
          addScaledVector(x, alpha_k, p_k, x);
          // r_{k+1} = r_k - alpha_k * v_k
          subScaledVector(residual_k, alpha_k, v_k, residual_k);

          // v_{k+1} = r_{k+1} but accumulated
          v_k = residual_k;
          m_communicator.accumulate(v_k);

          // rho_{k+1} = (v_{k+1}, r_{k+1}); r_k is distributed
          value_type rho_k_plus_1 = m_binary_measure(v_k, residual_k);

          // test for convergence (rho_{k+1} < epsilon_abs or
          // rho_{k+1} / rho_0 < epsilon_rel)
          m_rho_k = std::abs(rho_k_plus_1);
          if ((m_rho_k / m_rho_0 <= m_rel_tolerance)
              || (m_rho_k <= m_abs_tolerance))
          {
            m_convergence = true;
            return;
          }

          // calculate the new search direction
          // p_{k+1} = v_{k+1} + rho_{k+1} / rho_k * p_k
          addScaledVector(v_k, rho_k_plus_1 / rho_k, p_k, p_k);

          rho_k = rho_k_plus_1;
        }
      }

      //! \brief Return true, if the last CG inversion converged.
      bool convergence() const
      {
        return m_convergence;
      }

      //! \brief Returns the number of iterations needed for the last inversion.
      unsigned getIteration() const
      {
        return m_iteration;
      }

      //! \brief Get the initial residual measure.
      //!
      //! This method returns the initial residual measure defined by
      //! \f$ \vert (r_0, r_0) \vert \f$, whereby \f$ r_0 := F(x_0) \f$ is the
      //! initial residual, \p F is the forward operator, \f$ x_0 \f$
      //! is the starting guess for the inversion and
      //! \f$ (\cdot, \cdot) \f$ denotes the measurement operator.
      const double& getRho0() const
      {
        return m_rho_0;
      }

      //! \brief Get the final residual measure.
      //!
      //! In analogy to \p getRho0 this method returns the residual measure
      //! in the last iteration of the CG inversion.
      const double& getRho() const
      {
        return m_rho_k;
      }

    private:
      //! Object for parallel communication.
      TCommunicator& m_communicator;

      //! The forward operator.
      ForwardOperatorExpression<TForward>& m_forward;

      //! The object used to calculate the inner product.
      BinaryMeasureOperatorExpression<TBinaryMeasure>& m_binary_measure;

      //! Relative tolerance. If rho_{k+1} / rho_0 <= epsilon_rel, the CG
      //! inversion converged.
      double m_rel_tolerance;

      //! Absolute tolerance. If rho_{k+1} <= epsilon_abs, the inversion
      //! converged.
      double m_abs_tolerance;

      //! The maximum number of iterations before the CG is stopped
      //! without convergence.
      unsigned m_max_iterations;

      //! True, if the last CG inversion converged.
      bool m_convergence;

      //! The number of iterations needed for the last CG inversion.
      unsigned m_iteration;

      //! The initial residual measure.
      double m_rho_0;

      //! The final residual measure.
      double m_rho_k;
  };

} // namespace agile

#endif // AGILE_OPERATOR_CG_HPP

// End of $Id: cg.hpp 476 2011-06-16 08:54:14Z freiberger $.
