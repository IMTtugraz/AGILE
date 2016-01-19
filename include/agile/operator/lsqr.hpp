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

// $Id: lsqr.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_LSQR_HPP
#define AGILE_OPERATOR_LSQR_HPP

#include "agile/gpu_config.hpp"
#include "agile/operator/forward_operator.hpp"
#include "agile/operator/inverse_operator.hpp"
#include "agile/operator/binary_measure_operator.hpp"

namespace agile
{
  //! \brief LSQR method.
  template <typename TCommunicator, typename TForward, typename TBinaryMeasure>
  class LSQR
    : public InverseOperatorExpression<
               LSQR<TCommunicator, TForward, TBinaryMeasure> >
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
      LSQR(TCommunicator& communicator,
           const ForwardOperatorExpression<TForward>& forward,
           const BinaryMeasureOperatorExpression<TBinaryMeasure> &measure,
           const double& tolerance,
           unsigned max_iterations)
        : m_communicator(communicator), m_forward(forward()),
          m_binary_measure(measure()), m_tolerance(tolerance),
          m_max_iterations(max_iterations), m_convergence(false),
          m_iteration(0), m_rho_0(0), m_rho_k(0)
      {
      }

      //! \brief Apply the LSQR method to the vector x.
      //!
      //! This operator solves the problem  Ax = b for x.
      //! The right-hand side \p y has to be a distributed vector and the
      //! solution \p x will be an accumulated one. The initial solution
      //! is supplied via \p x as well and has to be accumulated, too.
      //!
      //! The "exactness" of the solution is measured using the measurement
      //! operator \p measure provided to the constructor and is denoted
      //! as \f$ (\cdot, \cdot) \f$. The residual norm in the k-th iteration is
      //! given by \f$ \varrho_k := (P^{-1}r_k, r_k) \f$, where
      //! \f$ r_k = y - F(x_k) \f$ is the residual, \f$ P \f$ is the
      //! The algorithm stops in iteration \p k when
      //! - the maximum amount of iterations is reached (\f$ k >= k_{max} \f$).
      //! - the residual is less than the absolute tolerance, i.e.
      //!   \f$ \varrho_k <= \varepsilon_{abs} \f$, or
      //! - the residual has decreased significantly such that the ratio
      //!   \f$ \varrho_k / \varrho_0 <= \varepsilon_{rel} \f$.
      //!
      //! \param[in] y The right-hand side (distributed).
      //! \param[in,out] x Provides the initial guess for the LSQR algorithm
      //! and returns the final solution (accumulated).
      template <typename TVectorType>
      void operator() (const TVectorType& y, TVectorType& x)
      {
        typedef typename TVectorType::value_type value_type;
        typedef typename to_real_type<value_type>::type real_type;

        // no convergence at the beginning
        m_convergence = false;
        m_iteration = 0;

        real_type n2y = getL2NormOfDistributedVector(y);
        real_type toly = m_tolerance * n2y;

        // the residual y - F(x_k)
        TVectorType residual_k(y.size());

        // temporary vectors
        TVectorType t_sy(y.size());
        TVectorType t_sx(x.size());

        // compute the initial residual r_0 = y - F(x_0), now t_sy is distributed
        // u = b - A*x;
        m_forward(x, t_sy);
        subVector(y, t_sy, residual_k);

        // L2 norm of residual
        // beta = norm(u);
        real_type beta = getL2NormOfDistributedVector(residual_k);
        m_rho_0 = beta;
        m_rho_k = m_rho_0;
        real_type normr = beta;

        // normalize the residual
        // u = u / beta
        if (beta != 0)
        {
          scale(real_type(1) / beta, residual_k, residual_k);
        }

        // v = A'*u;
        TVectorType v(x.size(), value_type(0));
        adjoint(m_forward)(residual_k, v);

        // normalize v
        real_type alpha = getL2NormOfDistributedVector(v);
        if (alpha != 0)
        {
          scale(real_type(1) / alpha, v, v);
        }

        real_type normar = alpha * beta;

        if (normar == 0)
        {
          // stop if norm of residum is 0
          // => x0 is the exact solution
          return;
        }

        if (n2y == 0) {
          // rhs vector y is zero, so x is also a zero vector
          x.assign(x.size(), value_type(0));
          return;
        }

        real_type norma = 0;
        real_type c = 1;
        real_type s = 0;
        real_type phibar = beta;
        TVectorType d(x.size(), value_type(0)); //d must be initialized to zero
        for (m_iteration = 0; m_iteration < m_max_iterations; ++m_iteration)
        {
          // u = A*v - alpha * u;
          m_forward(v, t_sy);
          subScaledVector(t_sy, alpha, residual_k, residual_k);
          // beta = norm(u)
          // normalize the residual
          beta = getL2NormOfDistributedVector(residual_k);
          // Normalize residual
          scale(real_type(1) / beta, residual_k, residual_k);
          //m_rho_k = beta;

          // norma = norm([norma alpha beta]);
          norma = std::sqrt(norma*norma + alpha*alpha + beta*beta);

          // thet = -s * alpha;
          real_type thet = -s * alpha;
          // rhot = c * alpha;
          real_type rhot = c * alpha;
          // sqrt(rhot^2 + beta^2);
          real_type rho = std::sqrt(rhot*rhot + beta*beta);
          // c = rhot / rho;
          c = rhot / rho;
          // s = -beta / rho;
          s = -beta / rho;
          // phi = c * phibar
          real_type phi = c * phibar;

          if (phi == 0)
          {
            // Stagnation @see matlab implementation
            break;
          }

          // phibar = s * phibar;
          phibar = s * phibar;

          // d = (v - thet * d) / rho;
          subScaledVector(v, thet, d, d);
          scale(real_type(1) / rho, d, d);

          // check for convergence in min{|b-A*x|}
          // absolute convergence
          if (normar / (norma*normr) <= m_tolerance)
          {
            m_convergence = true;
            break;
          }

          // check for convergence in A*x=b
          // relative convergence
          if (normr <= toly)
          {
            m_convergence = true;
            break;
          }

          // x = x + phi * d;
          addScaledVector(x, phi, d, x);
          // normr = abs(s) * normr;
          normr = std::abs(s) * normr;
          m_rho_k = normr;
          // vt = A' * u;
          adjoint(m_forward)(residual_k, t_sx);
          // v = vt - beta * v;
          subScaledVector(t_sx, beta, v, v);
          // alpha = norm(v);
          alpha = getL2NormOfDistributedVector(v);
          // v = v / alpha;
          scale(real_type(1) / alpha, v, v);
          // normar = alpha * abs(s * phi);
          normar = alpha * std::abs(s * phi);
        }
      }

      //! \brief Return true, if the last LSQR inversion converged.
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
      //! \f$ \vert (P^{-1}r_0, r_0) \vert \f$, whereby \f$ r_0 := F(x_0)
      //! \f$ is the initial residual, \p F is the forward operator, \f$ x_0 \f$
      //! is the starting guess for the inversion and
      //! \f$ (\cdot, \cdot) \f$ denotes the measurement operator.
      const double& getRho0() const
      {
        return m_rho_0;
      }

      //! \brief Get the final residual measure.
      //!
      //! In analogy to \p getRho0 this method returns the residual measure
      //! in the last iteration of the LSQR inversion.
      const double& getRho() const
      {
        return m_rho_k;
      }

    private:

      //! \brief Get L2 norm of the given vector.
      //!
      //! The norm is calculated using the scalar product:
      //! \f$ sqrt(\vert dot(v,v) \vert) \f$
      //! \param[in] vec The vector to compute the L2 norm of
      template <typename TVectorType>
      typename to_real_type<typename TVectorType::value_type>::type
      getL2NormOfDistributedVector(const TVectorType &vec)
      {
        //is this necessary
        // every call consumes additional gpu memory and frees it immediately
        // -> accumulate (i think this is not necessary with gpu,
        // but with finite elements)
        //
        //TVectorType vec_accum(vec);
        //m_communicator.accumulate(vec_accum);
        //return std::sqrt(std::abs(m_binary_measure(vec_accum, vec)));

        return std::sqrt(std::abs(m_binary_measure(vec, vec)));
      }

      //! Object for parallel communication.
      TCommunicator& m_communicator;

      //! The forward operator.
      TForward m_forward;

      //! The object used to calculate the inner product.
      TBinaryMeasure m_binary_measure;

      //! Absolute tolerance. If rho_{k+1} <= epsilon_abs, the inversion
      //! converged.
      double m_tolerance;

      //! The maximum number of iterations before the LSQR is stopped
      //! without convergence.
      unsigned m_max_iterations;

      //! True, if the last LSQR inversion converged.
      bool m_convergence;

      //! The number of iterations needed for the last LSQR inversion.
      unsigned m_iteration;

      //! The initial residual measure.
      double m_rho_0;

      //! The final residual measure.
      double m_rho_k;
  };

} // namespace agile

#endif // AGILE_OPERATOR_LSQR_HPP

// End of $Id: lsqr.hpp 476 2011-06-16 08:54:14Z freiberger $.
