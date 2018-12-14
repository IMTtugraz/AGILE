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

// $Id: irgn.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_IRGN_HPP
#define AGILE_OPERATOR_IRGN_HPP

class OperatorSamplingCartesian : public OperatorSamplingBase
{
  public:
    OperatorSampling(TV& pattern)
    : m_pattern(pattern), m_fft2()
    {
      m_fft2.createPlan(m_nx,m_nx);
    }

    //! \brief Forward sampling
    //!
    //! out = pattern*FFT(in)
    void forward(TV& in, TV& out)
    {
      m_fft2.forward(in,out);
      multiplyElementwise(m_pattern,out,out);
    }

    //! \brief Inverse sampling
    //!
    //! out = IFFT(pattern*in)
    void inverse(TV& in, TV& out)
    {
      multiplyElementwise(m_pattern,in,out);
      m_fft2.inverse(in,out);
    }

  private:
    //! Sampling pattern
    TV& m_pattern

    //! FFT operator
    agile::FFT2D m_fft2;
};

//! \brief Iteratively regularized Gauss-Newton solver.
//!
//! Minimizes \f$ \frac{1}{2} || Fx - y ||^2 + \frac{\alpha}{2} || x - x_0 ||^2 \f$
//! by solving \f$ (F'^*(x_k)F'(x_k) + \alpha I) \delta x = F'^*(x_k)(y - F(x_k)) + \alpha (x_0 - x) \f$
//!
class IRGN
{
  public:
    IRGN(TSolver& solver, const double& alpha0, const double& alpha_factor,
  unsigned max_iterations)
      : m_alpha_0(alpha0), m_alpha_factor(alpha_factor), m_alpha(m_alpha0),
  m_l(8), m_tau(2.2), m_gntol(0.75), m_max_iterations(max_iterations),
  m_iteration(0), m_solver(solver), m_x(nx*nx*(nc+1),0), m_u(m_x,0,nx*nx), m_c(mx,nx*nx,nc*nx*nx)
    {

      //TODO load data, pattern, weight

      //TODO scale data

      //TODO initial guess
      setVectorConstant(m_u,1);
      setVectorConstant(m_c,0);

      // generate sampling operator, derivative operator

      m_sampling = OperatorSamplingCartesian(m_pattern);

      //TODO transfer data to GPU

    }

    //! \brief Evaluate forward operator
    //!
    //! \f$ y = F'^*F'x+\alpha x\f$
    void forward(TV& x, TV& y)
    {
      tmp = new(x);
      derivative(x,tmp);
      adjoint_derivative(tmp,y);
      addScaledVector(y,m_alpha,x,y);
    }

    //! \brief Applies inverse Sobolev weight to sensitivities
    //!
    //! \f$ cw = FFT^{-1}(w\cdot c) \f$
    void apply_weight(TV& c,TV& cw)
    {
      multiplyElementwise(m_w,c,cw);
      m_fft2.inverse(cw,cw);
    }

    //! \brief Applies Adjoint of inverse Sobolev weighting
    //!
    //! \f$ c = \bar w \cdot FFT(cw) \f$
    void apply_weight_adjoint(TV& cw,TV& c)
    {
      m_fft2.forward(cw,c);
      multiplyConjElementwise(m_w,c,c);
    }

    //! \brief Compute the directional derivative.
    //!
    //! This method computes the directional derivative at the point \p x in
    //! the director \p dx and returns the result in \p dy. The evaluation
    //! point \p x has to be specified using \p set().
    void derivative(const TV& dx, TV& dy)
    {

      // get slices
      GPUVectorRange<typename TV::value_type> du(dx,0,m_nx*m_nx);
      tmp = new TV(m_nx,m_nx);
      for (i = 0;i<m_nc;i++)
      {
  GPUVectorRange<typename TV::value_type> dc(dx,(i+1)*m_nx*m_nx,m_nx*m_nx);
  GPUVectorRange<typename TV::value_type> c(m_cw,i*m_nx*m_nx,m_nx*m_nx);
  GPUVectorRange<typename TV::value_type> dyi(dy,i*m_nx*m_nx,m_nx*m_nx);

  apply_weight(dc,tmp);
  multiplyElementwise(m_u,tmp,dyi);
  multiplyElementwise(c,du,tmp);
  addVector(tmp,dyi,dyi);
  m_sampling.forward(dyi,dyi);
      }

    }

    //! \brief Compute the adjoint of the directional derivative.
    //!
    //! This method computes the adjoint of the directional derivative at the point \p x in
    //! the director \p dy and returns the result in \p dx. The evaluation
    //! point \p x has to be specified using \p set().
    void adjoint_derivative(const TV& dy, TV& dx)
    {

     // get slices
      GPUVectorRange<typename TV::value_type> du(dx,0,m_nx*m_nx);
      fdy = new TV(m_nx,m_nx);
      tmp = new TV(m_nx,m_nx);
      for (i = 0;i<m_nc;i++)
      {
  GPUVectorRange<typename TV::value_type> dyi(dy,i*m_nx*m_nx,m_nx*m_nx);
  GPUVectorRange<typename TV::value_type> c(m_cw,i*m_nx*m_nx,m_nx*m_nx);
  GPUVectorRange<typename TV::value_type> dc(dx,(i+1)*m_nx*m_nx,m_nx*m_nx);

  m_sampling.inverse(dyi,fdy);

  multiplyConjElementwise(c,fdy,tmp);
  addVector(tmp,du,du);
  multiplyConjElementwise(fdy,m_u,dc);
  apply_adjoint_weight(dc,dc);
      }

    }

    //! \brief Performs Gauß-Newton Iteration
    //!
    //! Solves for update \f$ \delta x = x^{(n+1)}-x^{(n)}\f$ using PCG
    void operator(TV& x, TV& y)
    {
      for (m_iteration = 0; m_iteration < m_max_iterations; ++m_iteration)
      {

  // preweight sensitivities
  apply_weight(m_c,m_cw);

  // compute residual
  m_sampling.forward(m_u,m_cw,sampled_data);
  residual = m_data - sampled_data;
  norm_res = residual.norm();

  // check for convergence
  if (norm_res < m_sigma)||(m_iteration>2)&&(norm_res > m_nr*m_gntol))
  {
    m_nr = norm_res;
    m_convergence = true;
    return;
  }
  m_nr = norm_res;
  m_pcgtol = (m_nr*m_alpha/3.)^2;

  // compute the right-hand side: rhs = F'^*(y - F(x)) + alpha* (x_0 - x)
  m_derivative_adjoint(residual, adj_residual);
  subVector(m_x0, m_x, rhs);
  addScaledVector(adj_residual, alpha, rhs, rhs);

  // tell the forward operator the current value of the regularization
  // parameter
  m_solve.getForward().setParameter(alpha);
  m_solve.setTolerance(m_pcgtol);
  m_solve(rhs, dx);

  // update iterates
  GPUVectorRange<typename TV::value_type> du(dx,0,m_nx*m_nx);
  GPUVectorRange<typename TV::value_type> dc(dx,m_nx*m_nx,nc*m_nx*m_nx);

  addVector(m_u, du, m_u);
  addVector(m_c, du, m_c);

  m_alpha = m_q*m_alpha;
      }
    }

    //! \brief Postprocesses iterate
    //!
    //! \f$ us = \frac{1}{dscale} u \cdot \sqrt{\sum_i |W(c_i)|^2}  \f$
    void postprocess()
    {
      apply_weight(c,cw);
      cscale = sqrt(multiplyConjElementwise(cw,cw,cw).sum(3));
      m_un = ScaleVector(cscale/m_dscale, u);
    }

  private:

    // Fixed parameters

    //! Initial regularization parameter.
    double m_alpha_0;

    //! Factor to reduce the regularization parameter in each iteration.
    double m_alpha_factor;

    //! Order of the Sobolev penalty for sensitivities \f$ \|\Delta^{l/2} c\|^2 \f$.
    unsigned m_l;

    //! The maximum iterations.
    unsigned m_max_iterations;

    //! Scaling for stopping criterion
    double m_tau;

    //! Fallback stopping criterion: res_k > res_{k-1}*m_gntol.
    double m_gntol;

    // Data dependent parameters

    //! Number of pixels in image.
    unsigned m_nx;

    //! Number of coils.
    unsigned m_nc;

    //! Stopping criterion: res_k < sigma
    double m_sigma;

    //! Scaling factor.
    double m_dscale;

    // Iteration dependent parameters

    //! The current iteration.
    unsigned m_iteration;

    //! Current regularization parameter.
    double m_alpha;

    //! Current residuum norm.
    double m_nr;

    //! Tolerance for PCG.
    double m_pcgtol;

    // Data

    //! Sampling pattern
    TV& m_pattern;

    //! Sobolev weight
    TV& m_w;

    //! Measured data
    TV& m_data;

    //! Initial guess
    TV& m_x0;

    // Iterates:

    //! Current iterate ((nc+1)*nx)
    TV& m_x;

    //! Current image guess (nx*nx)
    GPUVectorRange<typename TV::value_type> m_u;

    //! Current sensitivity guess (nx*nx*nc)
    GPUVectorRange<typename TV::value_type> m_c;

    //! Current sensitivity guess preweighted
    TV& m_cw;

    // Operators

    //! PCG solution operator
    TSolver& m_solver;

    //! Sampling operator
    OperatorSampling m_sampling;

};

#endif // AGILE_OPERATOR_IRGN_HPP

// End of $Id: irgn.hpp 476 2011-06-16 08:54:14Z freiberger $.

