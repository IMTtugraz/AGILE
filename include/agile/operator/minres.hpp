// $Id: minres.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_MINRES_HPP
#define AGILE_OPERATOR_MINRES_HPP

#include "agile/gpu_config.hpp"
#include "agile/operator/forward_operator.hpp"
#include "agile/operator/inverse_operator.hpp"
#include "agile/operator/binary_measure_operator.hpp"

namespace agile
{
  //! \brief
  template <typename TCommunicator, typename TForward,
            typename TPreconditioner, typename TBinaryMeasure>
  class MINRES
    : public InverseOperatorExpression<
               PreconditionedConjugateGradient<
                 TCommunicator, TForward, TPreconditioner, TBinaryMeasure> >
  {
    public:
      //! \brief Constructor.
      //!
      //! \param[in] communicator Object for parallel communication.
      //! \param[in] forward The forward operator.
      //! \param[in] preconditioner The operator used for preconditioning.
      //! \param[in] measure A measurement operator.
      //! \param[in] rel_tolerance The relative tolerance
      //! \f$ \varepsilon_{rel} \f$ serving as 1st stopping criterion for the
      //! algorithm.
      //! \param[in] abs_tolerance The absolute tolerance
      //! \f$ \varepsilon_{abs} \f$ (2nd stopping criterion).
      //! \param[in] max_iterations The maximum number of iterations for the
      //! inversion after which the algorithm gives up.
      MINRES(
        TCommunicator& communicator,
        const ForwardOperatorExpression<TForward>& forward,
        const InverseOperatorExpression<TPreconditioner>& preconditioner,
        const BinaryMeasureOperatorExpression<TBinaryMeasure>& measure,
        const double& rel_tolerance, const double& abs_tolerance,
        unsigned max_iterations)
        : m_communicator(communicator), m_forward(forward()),
          m_preconditioner(preconditioner()), m_binary_measure(measure()),
          m_rel_tolerance(rel_tolerance), m_abs_tolerance(abs_tolerance),
          m_max_iterations(max_iterations), m_convergence(false),
          m_iteration(0), m_rho_0(0), m_rho_k(0)
      {
      }

      template <typename TVectorType>
      void operator() (const TVectorType& y, TVectorType& x)
      {
        typedef typename TVectorType::value_type value_type;

  //setVectorConstant(x, 0);
  // r = b
  // d = r
  TVectorType d = y;
  TVectorType r = y;
  TVectorType ro;
  TVectorType Ad;
  TVectorType Ar;
  TVectorType Aro;

  // Ar = A*r ;
  m_forward(r, Ar);
  Ad = Ar;

  // res0 = r'*r;
  value_type res0 = m_binary_measure(r, r);
  m_rho_0 = std::abs(res0);
        for (m_iteration = 0; m_iteration < m_max_iterations; ++m_iteration)
        {
      ro = r;
      Aro = Ar;
      // alpha = (r'*Ar) / (Ad'*Ad);
      value_type alpha = m_binary_measure(r, Ar) / m_binary_measure(Ad, Ad);
      // x = x + alpha*d;
      addScaledVector(x, alpha, d, x);
      // r = r - alpha*Ad;
      subScaledVector(r, alpha, Ad, r);
      // res = r'*r;
      value_type res = m_binary_measure(r, r);

      m_rho_k = std::abs(res);

      if (m_rho_k < (m_rel_tolerance*m_rel_tolerance*m_rho_0))
      {
        break;
      }
      //Ar = A*r;
      m_forward(r, Ar);
      // beta = (r'*Ar) / (ro'*Aro);
      value_type beta = m_binary_measure(r, Ar) / m_binary_measure(ro, Aro);
      // d =  r + beta*d;
      addScaledVector(r, beta, d, d);
      // Ad = Ar + beta*Ad;
      addScaledVector(Ar, beta, Ad, Ad);

        }
      }

      //! \brief Return true, if the last PCG inversion converged.
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
      //! \f$ \vert (P^{-1}r_0, r_0) \vert \f$, whereby \p P is the
      //! preconditioning operator, \f$ r_0 := F(x_0) \f$ is the
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
      //! in the last iteration of the PCG inversion.
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

      //! Relative tolerance. If rho_{k+1} / rho_0 <= epsilon_rel, the PCG
      //! inversion converged.
      double m_rel_tolerance;

      //! Absolute tolerance. If rho_{k+1} <= epsilon_abs, the inversion
      //! converged.
      double m_abs_tolerance;

      //! The maximum number of iterations before the PCG is stopped
      //! without convergence.
      unsigned m_max_iterations;

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

#endif // AGILE_OPERATOR_MINRES_HPP

// End of $Id: minres.hpp 476 2011-06-16 08:54:14Z freiberger $
