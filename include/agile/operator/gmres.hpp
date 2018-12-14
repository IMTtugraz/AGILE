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

// $Id: gmres.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_GMERS_HPP
#define AGILE_OPERATOR_GMERS_HPP

#include "agile/gpu_config.hpp"
#include "agile/operator/forward_operator.hpp"
#include "agile/operator/inverse_operator.hpp"
#include "agile/operator/binary_measure_operator.hpp"
#include "agile/gpu_type_traits.hpp"
#include "agile/basic_math.hpp"

namespace agile
{
  //! Generate a Givens rotation.
  //!
  //! |  c  s | * | a | = | r |,    c^2 + s^2 = 1
  //! | -s  c |   | b |   | 0 |
  //!
  //! c = a / sqrt(a^2 + b^2)
  //! s = b / sqrt(a^2 + b^2)
  //! r = sqrt(a^2 + b^2)
  //!
  //! Anderson, E.: "Discontinuous Plane Rotations and the Symmetric Eigenvalue
  //! Problem". LAPACK Working Note 150, University of Tennessee, 2000.
  template <typename TType>
  inline
  void getGivensRotation(const TType& a, const TType& b,
                         TType& c, TType& s, TType& r)
  {

    if (b == TType(0))
    {
      c = (a < TType(0)) ? TType(-1) : TType(1);
      s = TType(0);
      r = std::abs(a);
    }
    else if (a == TType(0))
    {
      c = TType(0);
      s = (b < TType(0)) ? TType(-1) : TType(1);
      r = std::abs(b);
    }
    else if (std::abs(a) > std::abs(b))
    {
      const TType t = b / a;
      const TType u = std::sqrt(TType(1) + t*t)
                      * ((a < TType(0)) ? TType(-1) : TType(1));
      c = TType(1) / u;
      s = c*t;
      r = a*u;
    }
    else
    {
      const TType t = a / b;
      const TType u = std::sqrt(TType(1) + t*t)
                      * ((b < TType(0)) ? TType(-1) : TType(1));
      s = TType(1) / u;
      c = s*t;
      r = b*u;
    }
  }

  //! Generate a complex Givens rotation.
  //!
  //! |  c        s | * | f | = | r |,   c^2 + s*conj(s) = 1
  //! | -conj(s)  c |   | g |   | 0 |
  //!
  //! r = sgn(Re(f)) * sgn(f) * sqrt(|f|^2 + |g|^2)
  //! c = f / r
  //! s = conj(g) / conj(r)
  //!
  //! c is always real. The complex sgn-function is defined as
  //! sgn(x) = | 1       if x = 0
  //!          | x/|x|   otherwise
  //!
  //! Anderson, E.: "Discontinuous Plane Rotations and the Symmetric Eigenvalue
  //! Problem". LAPACK Working Note 150, University of Tennessee, 2000.
  template <typename TType>
  inline
  void getGivensRotation(const std::complex<TType>& f,
                         const std::complex<TType>& g,
                         TType& c, std::complex<TType>& s,
                         std::complex<TType>& r)
  {
    if (g == std::complex<TType>(0))
    {
      c = std::real(f) < TType(0) ? TType(-1) : TType(1);
      s = std::complex<TType>(0);
      r = f * c;
    }
    else if (f == std::complex<TType>(0))
    {
      c = 0;
      s = sgn(std::conj(g));
      r = std::abs(g);
    }
    else
    {
      TType f1 = std::abs(std::real(f)) + std::abs(std::imag(f));
      const TType g1 = std::abs(std::real(g)) + std::abs(std::imag(g));
      if (f1 >= g1)
      {
        const std::complex<TType> fs = f / f1;
        const TType f2 = std::real(fs) * std::real(fs)
                         + std::imag(fs) * std::imag(fs);
        const std::complex<TType> gs = g / f1;
        const TType g2 = std::real(gs) * std::real(gs)
                         + std::imag(gs) * std::imag(gs);
        const TType u = (std::real(f) < TType(0) ? TType(-1) : TType(1))
                        * std::sqrt(TType(1) + g2 / f2);
        c = TType(1) / u;
        s = std::conj(gs) * fs * (c / f2);
        r = f * u;
      }
      else
      {
        std::complex<TType> fs = f / g1;
        const TType f2 = std::real(fs) * std::real(fs)
                         + std::imag(fs) * std::imag(fs);
        const std::complex<TType> gs = g / g1;
        const TType g2 = std::real(gs) * std::real(gs)
                         + std::imag(gs) * std::imag(gs);
        const TType u = (std::real(f) < TType(0) ? TType(-1) : TType(1)) * g1
                        * std::sqrt(f2 + g2);
        f1 = std::abs(f);
        fs = f / f1;
        c = f1 / u;
        s = fs * (std::conj(g) / u);
        r = fs * u;
      }
    }
  }

  //! \brief Generalized minimal residual method.
  //!
  //! This code is based uppon
  //! Hanson, R.J. and Kincaid D.R.: "Notes on GMRES Algorithm Organization".
  //! Technical Report TR-05-05. University of Texas at Austin. March 2005.
  template <typename TCommunicator, typename TForward,
            typename TPreconditioner, typename TBinaryMeasure>
  class GMRES
    : public InverseOperatorExpression<
               GMRES<TCommunicator, TForward, TPreconditioner, TBinaryMeasure> >
  {
    public:
      //! \brief Constructor.
      //!
      //! \param[in] communicator Object for parallel communication.
      //! \param[in] forward The forward operator.
      //! \param[in] preconditioner The operator used for preconditioning.
      //! \param[in] measure A measurement operator.
      //! TODO: finish this doc
      GMRES(TCommunicator& communicator,
            const ForwardOperatorExpression<TForward>& forward,
            const InverseOperatorExpression<TPreconditioner>& preconditioner,
            const BinaryMeasureOperatorExpression<TBinaryMeasure>& measure,
            const double& abs_tolerance, unsigned max_inner_iterations,
            unsigned max_outer_iterations)
        : m_communicator(communicator), m_forward(forward()),
          m_preconditioner(preconditioner()), m_binary_measure(measure()),
          m_abs_tolerance(abs_tolerance),
          m_max_inner_iterations(max_inner_iterations),
          m_max_outer_iterations(max_outer_iterations), m_convergence(false),
          m_iteration(0), m_rho_0(0), m_rho_k(0)
      {
      }

      //! \param[in] y The right-hand side (distributed).
      //! \param[in,out] x Provides the initial guess for the PCG algorithm
      //! and returns the final solution (accumulated).
      template <typename TVectorType>
      void operator() (const TVectorType& y, TVectorType& x)
      {
        typedef typename TVectorType::value_type value_type;
        typedef typename to_real_type<value_type>::type real_type;

        m_convergence = false;

        // get memory for the basis vectors (we are going to construct one
        // basis vector per inner iteration)
        std::vector<TVectorType> V(m_max_inner_iterations);
        for (unsigned counter = 0; counter < V.size(); ++counter)
          V[counter].resize(x.size());

        // temporary vectors
        TVectorType p_k(y.size(), 0);
        TVectorType v_k(y.size(), 0);

        // upper Hessenberg matrix of size
        // (max_inner_iterations x max_inner_iterations)
        // rho stores the values in the first sub-diagonal from index 1 on;
        // in index 0 we store ||P^{-1}(y - F(x_0))|| for later
        std::vector<real_type> rho(m_max_inner_iterations + 1);
        // H stores the values in the upper triangle in row-major form
        std::vector<value_type> H(
          m_max_inner_iterations * (m_max_inner_iterations + 1) / 2);

        for (m_iteration = 0; m_iteration < m_max_outer_iterations;
             ++m_iteration)
        {
          // ---------- Phase 1: create basis vectors ----------

          // v_0 = x_0 (both are accumulated)
          v_k = x;
          unsigned basis_counter;
          for (basis_counter = 0;
               basis_counter <= m_max_inner_iterations; ++basis_counter)
          {
            // in the very first iteration it applies that v_k = x_0 in all
            // other iterations v_k = V_{k-1} (with v_k being accumulated)
            // p_0 = F(x_0) or p_k = F(V_{k-1}) for k > 0; p_k is distributed
            // then
            m_forward(v_k, p_k);
            // the first basis vector is
            // V_0 = r_0 / ||r_0|| with r_0 = P^{-1}(y - F(x_0))
            // while the others are orthogonalized versions of F(V_{k-1})
            if (basis_counter == 0)
              subVector(y, p_k, p_k);

            // apply the preconditioner; v_k will be accumulated
            m_preconditioner(p_k, v_k);
            // orthogonalize against the previous basis vectors; this loop is
            // skipped in the first round as there is nothing to orthogonalize
            for (unsigned counter = 0; counter < basis_counter; ++counter)
            {
              // compute the inner product with a previous basis vector; as
              // v_k is accumulated and V[counter], too, we need to distribute
              // one of them
              p_k = V[counter];
              m_communicator.distribute(p_k);
              // H[counter][basis_counter - 1] = (v_k, V[counter])
              unsigned index = counter * m_max_inner_iterations
                               - counter * (counter+1) / 2 + basis_counter - 1;
              H[index] = m_binary_measure(v_k, p_k);
              // v_k -= (v_k, V[counter]) * V[counter];
              subScaledVector(v_k, H[index], V[counter], v_k);
            }
            // at this point v_k is orthogonal to all V[i] with 0 <= i < k

            // we have to get the norm of v_k to normalize the vector;
            // unfortunately the binary measure wants one accumulated and one
            // distibuted vector but we only have an accumulated v_k until now
            p_k = v_k;
            m_communicator.distribute(p_k);
            rho[basis_counter]
              = std::sqrt(std::abs(m_binary_measure(v_k, p_k)));
            // note: H[basis_counter][basis_counter - 1] = (v_k, v_k) and this
            // is stored in rho[basis_counter]
            if (m_iteration == 0 && basis_counter == 0)
              m_rho_0 = rho[0];

            // we only need max_inner_iterations basis vectors, but calculate
            // one more inner product as residual measure; we can leave the
            // inner loop in the last round right here, therefore;
            // if the residual is smaller than the absolute tolerance, we can
            // leave, too
            if ((rho[basis_counter] <= m_abs_tolerance)
                || (basis_counter == m_max_inner_iterations))
              break;

            // create the basis vector by normalizing v_k (V[k] is accumulated)
            scale(value_type(1) / rho[basis_counter], v_k, v_k);
            V[basis_counter] = v_k;
          }
          // now we have basis_counter basis vectors
          // the part H[0:basis_counter-1][0:basis_counter-1] is filled and
          // the sub-diagonal of the Hessenberg is in rho[1:basis_counter]

          // store the current residual
          m_rho_k = rho[0];

          // ---------- Phase 2: update phase ----------

          // if basis_counter is zero, either max_inner_iterations is zero
          // too or the algorithm converged (this would be a good case!)
          // in both cases we leave the method because there is nothing to
          // update
          if (basis_counter == 0)
          {
            m_convergence = (m_rho_k <= m_abs_tolerance);
            return;
          }

          // zero the sub-diagonal of H by Givens rotations
          std::vector<value_type> z(basis_counter + 1, 0);
          z[0] = rho[0];
          // point in H to the first element of the current row
          unsigned row_offset = 0;
          for (unsigned row = 0; row < basis_counter; ++row)
          {
            // construct a Givens rotation from H[row:row+1][row]
            // to zero H[row+1][row] (which is stored in rho[row+1])
            real_type c;
            value_type s;
            getGivensRotation(H[row_offset], value_type(rho[row + 1]),
                              c, s, H[row_offset]);
            // determine the amount of elements in the current row in H
            unsigned row_length = m_max_inner_iterations - row;
            // apply the Givens rotation to H[row:row+1][row+1:end]
            for (unsigned counter = 0; counter < basis_counter - row - 1;
                 ++counter)
            {
              value_type f = H[row_offset + 1 + counter];
              value_type g = H[row_offset + row_length + counter];
              H[row_offset + 1 + counter] = c * f + s * g;
              H[row_offset + row_length + counter] = -conj(s) * f + c * g;
            }
            // apply the Givens rotation to y, too; consider that z[row+1] is
            // still zero
            value_type t = c * z[row];
            z[row + 1] = -conj(s) * z[row];
            z[row] = t;
            // go to the next row
            row_offset += row_length;
          }
          // H is now upper triangluar

          // point to the start of row (basis_counter - 1)
          row_offset = (basis_counter - 1) * (m_max_inner_iterations+1)
                       - (basis_counter - 1) * basis_counter / 2;
          // backward solve
          z[basis_counter - 1] /= H[row_offset];
          for (int row = basis_counter - 2; row >= 0; --row)
          {
            unsigned row_length = m_max_inner_iterations - row;
            // point to the beginning of the current row
            row_offset -= row_length;
            value_type sum = z[row];
            for (unsigned counter = 1; counter < row_length; ++counter)
              sum -= H[row_offset + counter] * z[row + counter];
            z[row] = sum / H[row_offset];
          }

          //update the solution
          for (unsigned counter = 0; counter < basis_counter; ++counter)
            addScaledVector(x, z[counter], V[counter], x);

#if 0
TODO
          // test stagnation
          if (std::abs(m_rho_k - std::abs(z[basis_counter]))
              <= sqrt(eps) * m_rho_0)
          {
            return;
          }
#endif
        } // outer iteration
      }

      //! \brief Return true, if the last GMRES inversion converged.
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
      const double& getRho0() const
      {
        return m_rho_0;
      }

      //! \brief Get the final residual measure.
      const double& getRho() const
      {
        return m_rho_k;
      }

    private:
      //! Object for parallel communication.
      TCommunicator& m_communicator;

      //! The forward operator.
      TForward m_forward;

      //! The preconditioner.
      TPreconditioner m_preconditioner;

      //! The object used to calculate the inner product.
      TBinaryMeasure m_binary_measure;

      //! Absolute tolerance.
      double m_abs_tolerance;

      //! The amount of inner iterations is equal to the amount of basis
      //! vectors created for the Krylov subspace.
      unsigned m_max_inner_iterations;

      //! The number of outer iterations
      unsigned m_max_outer_iterations;

      //! True, if the last PCG inversion converged.
      bool m_convergence;

      //! The number of iterations needed for the last PCG inversion.
      unsigned m_iteration;

      //! The initial residual measure.
      double m_rho_0;

      //! The final residual measure.
      double m_rho_k;
  };

} // namespace agile

#endif // AGILE_OPERATOR_GMERS_HPP

// End of $Id: gmres.hpp 476 2011-06-16 08:54:14Z freiberger $.
