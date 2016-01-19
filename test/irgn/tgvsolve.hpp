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
//
//==================================================================================
//
//                    TGVSOLVE.HPP
//
//==================================================================================

#ifndef AGILE_IRGN_TGVSOLVE_HPP
#define AGILE_IRGN_TGVSOLVE_HPP


#include "irgn.hpp"

using agile::GPUMatrix;
using agile::GPUVector;


template <typename TType>
class TGVSolve : public IRGN<TType>
{
  public:

    //! \brief Default constructor.
    //!
    //! The default constructor creates an empty TGVSolve Class.
    TGVSolve()
      : IRGN<TType>()
    {
    }


    //! \brief Constructor.
    //!
    //! The Constructor creates an TGVSolve.
    TGVSolve(GPUMatrix<TType>* coil, unsigned int num_coils, IRGN_Params param)
      : IRGN<TType>(coil, num_coils, param)
    {
      Init();
    }

    //! \brief Constructor.
    //!
    //! The Constructor creates an initialized TGVSolve.
    TGVSolve(GPUMatrix<TType>* coil, unsigned int num_coils)
        : IRGN<TType>(coil, num_coils)
    {
      Init();
    }

    //! \brief Destructor.
    virtual ~TGVSolve()
    {
      delete[] x2_mat_;
      x2_mat_ = 0;
      delete[] y2_mat_;
      y2_mat_ = 0;
      delete[] cb_mat_;
      cb_mat_ = 0;
      delete[] Mc_mat_;
      Mc_mat_ = 0;
      delete[] eta4_mat_;
      eta4_mat_ = 0;
      delete[] dc_old_mat_;
      dc_old_mat_ = 0;
    }

    void Solve(const GPUMatrix<TType>* u, const GPUMatrix<TType>* c,
                                   const GPUMatrix<TType>* rhs, const GPUMatrix<TType>* u0,
                                   unsigned maxits, float alpha, float beta,
                                   GPUMatrix<TType>* du, GPUMatrix<TType>* dc);


  private:

    void Init();


    GPUMatrix<TType> x1_mat_;
    GPUMatrix<TType> y1_mat_;
    GPUMatrix<TType> ub_mat_;
    GPUMatrix<TType> ukp_mat_;
    GPUMatrix<TType> Mu_mat_;
    GPUMatrix<TType> eta3_mat_;
    GPUMatrix<TType> eta1_mat_;
    GPUMatrix<TType> eta2_mat_;
    GPUMatrix<TType> y1p_mat_;
    GPUMatrix<TType> y2p_mat_;
    GPUMatrix<TType> p1_mat_;
    GPUMatrix<TType> p2_mat_;
    GPUMatrix<TType> du_old_mat_;
    GPUMatrix<TType> safe_;
    GPUMatrix<TType> safe2_;
    GPUMatrix<typename agile::to_real_type<TType>::type> safe_real_;

    GPUMatrix<TType> v1_mat_;
    GPUMatrix<TType> v2_mat_;
    GPUMatrix<TType> vb1_mat_;
    GPUMatrix<TType> vb2_mat_;
    GPUMatrix<TType> zeta1_mat_;
    GPUMatrix<TType> zeta2_mat_;
    GPUMatrix<TType> zeta3_mat_;
    GPUMatrix<TType> zeta4_mat_;
    GPUMatrix<TType> zeta5_mat_;
    GPUMatrix<TType> q1_mat_;
    GPUMatrix<TType> q2_mat_;
    GPUMatrix<TType> q3_mat_;
    GPUMatrix<TType> v1_old_mat_;
    GPUMatrix<TType> v2_old_mat_;


    GPUMatrix<TType>* x2_mat_;
    GPUMatrix<TType>* y2_mat_;
    GPUMatrix<TType>* cb_mat_;
    GPUMatrix<TType>* Mc_mat_;
    GPUMatrix<TType>* eta4_mat_;
    GPUMatrix<TType>* dc_old_mat_;

    typename agile::to_real_type<TType>::type norm_;
    typename agile::to_real_type<TType>::type one_;
    typename agile::to_real_type<TType>::type two_;

    std::complex<typename agile::to_real_type<TType>::type> l1_;
    std::complex<typename agile::to_real_type<TType>::type> l2_;
    std::complex<typename agile::to_real_type<TType>::type> l2_safe_;

};

//===========================================================================
//
//              M E T H O D E S
//
//===========================================================================


// ---------------------------------------------------------------
//! \brief Init()
//!  special initialize for class IRGN.
// ---------------------------------------------------------------
template <typename TType>
void TGVSolve<TType>::Init()
{
  norm_=0;
  one_ =1;
  two_ = 2;

  unsigned int num_rows = IRGN<TType>::num_rows_;
  unsigned int num_columns = IRGN<TType>::num_columns_;

  x2_mat_ = new GPUMatrix<TType> [IRGN<TType>::num_coils_];
  y2_mat_ = new GPUMatrix<TType> [IRGN<TType>::num_coils_];
  cb_mat_ = new GPUMatrix<TType> [IRGN<TType>::num_coils_];
  Mc_mat_ = new GPUMatrix<TType> [IRGN<TType>::num_coils_];
  eta4_mat_ = new GPUMatrix<TType> [IRGN<TType>::num_coils_];
  dc_old_mat_ = new GPUMatrix<TType> [IRGN<TType>::num_coils_];


  //initialize two-dimensional matrix
  x1_mat_.resize(num_rows, num_columns);
  y1_mat_.resize(num_rows, num_columns);
  ub_mat_.resize(num_rows, num_columns);
  ukp_mat_.resize(num_rows, num_columns);
  Mu_mat_.resize(num_rows, num_columns);
  eta1_mat_.resize(num_rows, num_columns);
  eta2_mat_.resize(num_rows, num_columns);
  y1p_mat_.resize(num_rows, num_columns);
  y2p_mat_.resize(num_rows, num_columns);
  eta3_mat_.resize(num_rows, num_columns);
  du_old_mat_.resize(num_rows, num_columns);
  safe_.resize(num_rows, num_columns);
  safe2_.resize(num_rows, num_columns);
  p1_mat_.resize(num_rows, num_columns);
  p2_mat_.resize(num_rows, num_columns);
  safe_real_.resize(num_rows, num_columns);

  v1_mat_.resize(num_rows, num_columns);
  v2_mat_.resize(num_rows, num_columns);
  vb1_mat_.resize(num_rows, num_columns);
  vb2_mat_.resize(num_rows, num_columns);

  zeta1_mat_.resize(num_rows, num_columns);
  zeta2_mat_.resize(num_rows, num_columns);
  zeta3_mat_.resize(num_rows, num_columns);
  zeta4_mat_.resize(num_rows, num_columns);
  zeta5_mat_.resize(num_rows, num_columns);
  q1_mat_.resize(num_rows, num_columns);
  q2_mat_.resize(num_rows, num_columns);
  q3_mat_.resize(num_rows, num_columns);
  v1_old_mat_.resize(num_rows, num_columns);
  v2_old_mat_.resize(num_rows, num_columns);



  //initialize multidimensional matrix
  for(unsigned i=0; i<IRGN<TType>::num_coils_;i++)
  {
    x2_mat_[i].resize(num_rows, num_columns);
    y2_mat_[i].resize(num_rows, num_columns);
    cb_mat_[i].resize(num_rows, num_columns);  //achtung - size nxnx
    Mc_mat_[i].resize(num_rows, num_columns);
    eta4_mat_[i].resize(num_rows, num_columns);
    dc_old_mat_[i].resize(num_rows, num_columns);
  }

  std::cout<<"\n\n INIT \n\n";


}


// ---------------------------------------------------------------
//! \brief TGVSolve()
//!  calculates TGVSolve
// ---------------------------------------------------------------
template <typename TType>
void TGVSolve<TType>::Solve(const GPUMatrix<TType>* u, const GPUMatrix<TType>* c,
                               const GPUMatrix<TType>* rhs, const GPUMatrix<TType>* u0,
                               unsigned maxits, float alpha, float beta,
                               GPUMatrix<TType>* du, GPUMatrix<TType>* dc)
{
  l2_.real(0);
  l2_.imag(0);

  agile::copy(IRGN<TType>::random1_mat_,x1_mat_);  //achtung - size nxnx - zero_complex_mat_nxnx_
  // initialize iterates du, dc and ub, cb
  agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,*du);  //achtung - size nxnx - zero_complex_mat_nxnx_
  agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,ub_mat_);  //achtung - size nxnx - zero_complex_mat_nxnx_
  agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,p1_mat_);  //achtung - size nxnx - zero_complex_mat_nxnx_
  agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,p2_mat_);  //achtung - size nxnx - zero_complex_mat_nxnx_
  agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,v1_mat_);  //achtung - size nxnx - zero_complex_mat_nxnx_
  agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,v2_mat_);  //achtung - size nxnx - zero_complex_mat_nxnx_
  agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,vb1_mat_);  //achtung - size nxnx - zero_complex_mat_nxnx_
  agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,vb2_mat_);  //achtung - size nxnx - zero_complex_mat_nxnx_
  agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,q1_mat_);  //achtung - size nxnx - zero_complex_mat_nxnx_
  agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,q2_mat_);  //achtung - size nxnx - zero_complex_mat_nxnx_
  agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,q3_mat_);  //achtung - size nxnx - zero_complex_mat_nxnx_
  for(unsigned i=0; i < IRGN<TType>::num_coils_;i++)
  {
    agile::copy(IRGN<TType>::random2_mat_[i],x2_mat_[i]);  //achtung - size nxnx - zero_complex_mat_nxnx_
    agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,dc[i]);  //achtung - size nxnx - zero_complex_mat_nxnx_
    agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,cb_mat_[i]);  //achtung - size nxnx - zero_complex_mat_nxnx_
  }

  std::cout<<"\n <IRGN<TType>::num_coils_: "<<IRGN<TType>::num_coils_;

  ApplyM(&y1_mat_,y2_mat_,&x1_mat_,x2_mat_);

  for (unsigned i = 0; i < 10; i++)
  {
    //x1=y1./norm(y1)
    norm_ = agile::norm2(y1_mat_);
    if (norm_ == 0)
      agile::copy(y1_mat_, x1_mat_);
    else
      agile::scale(one_/norm_,y1_mat_,x1_mat_);

    //x2=y2./norm(y2)
    norm_ = 0;
    for(unsigned i=0; i<IRGN<TType>::num_coils_;i++)
    {
      norm_ = norm_ + std::pow(agile::norm2(y2_mat_[i]),2);
    }
    norm_ = std::sqrt(norm_);

    if (norm_ == 0)
      IRGN<TType>::CopyMatrixZ(y2_mat_, x2_mat_, IRGN<TType>::num_coils_);
    else
    {
      for(unsigned i=0; i<IRGN<TType>::num_coils_;i++)
        agile::scale(one_/norm_,y2_mat_[i],x2_mat_[i]);
    }
    ApplyM(&y1_mat_,y2_mat_,&x1_mat_,x2_mat_);
  }

  agile::dotProduct(x1_mat_,y1_mat_,l1_);

  for(unsigned i=0; i<IRGN<TType>::num_coils_;i++)
  {
    agile::dotProduct(x2_mat_[i],y2_mat_[i],l2_safe_);
    l2_ = l2_ + l2_safe_;
  }

  typename agile::to_real_type<TType>::type L;
  L = 2 * ((std::abs(l1_) >= std::abs(l2_)) ? std::abs(l1_) : std::abs(l2_));
  typename agile::to_real_type<TType>::type tau;
  tau = 1/std::sqrt(8+L);

  typename agile::to_real_type<TType>::type sigma;
  sigma = 1/std::sqrt(8+L);


  for (unsigned iter = 1; iter <= maxits; iter++)
  {
/*
    % update dual (p): (eta1,eta2) = nabla(u+du-u0)
    ukp = ub + u - u0;
    eta1 = dx(ukp) - vb1;
    eta2 = dy(ukp) - vb2;
    y1 = p1 + sigma * eta1;
    y2 = p2 + sigma * eta2;
    my = proj(y1,y2);
    p1 = y1./my;
    p2 = y2./my;
*/
    //update dual (p): (eta1,eta2) = nabla(u+du-u0)

    //ukp = ub + u - u0;
    agile::addMatrix(ub_mat_,*u,safe_);


    agile::subMatrix(safe_,*u0,ukp_mat_);

    //eta1 = dx(ukp);
    agile::diff(1,ukp_mat_,eta1_mat_);
    agile::subMatrix(eta1_mat_,vb1_mat_,eta1_mat_);
    //eta2 = dy(ukp);
    agile::diff(2,ukp_mat_,eta2_mat_);
    agile::subMatrix(eta2_mat_,vb2_mat_,eta2_mat_);

    //y1 = p1 + sigma * eta1;
    agile::scale(sigma,eta1_mat_,eta1_mat_);
    agile::addMatrix(p1_mat_,eta1_mat_,y1p_mat_);
    //y2 = p2 + sigma * eta2;
    agile::scale(sigma,eta2_mat_,eta2_mat_);
    agile::addMatrix(p2_mat_,eta2_mat_,y2p_mat_);

    //proj = @(y1,y2) max(1,sqrt(y1.^2+y2.^2)/beta);
    agile::multiplyElementwise(y1p_mat_,y1p_mat_,y1p_mat_);
    agile::multiplyElementwise(y2p_mat_,y2p_mat_,y2p_mat_);
    agile::addMatrix(y1p_mat_,y2p_mat_,safe_);
    agile::real(safe_,safe_real_);
    agile::sqrt(safe_real_,safe_real_);
    if(beta == 0)
      beta=1;
    agile::scale(one_/beta,safe_real_,safe_real_);

    agile::max(this->ones_mat_nxnx_,safe_real_,safe_real_); //safe_real_==my
    //safe_real_ = 1 / safe_real_
    agile::pow(typename TType::value_type(-1),safe_real_,safe_real_);

    //p1 = y1./my;
    agile::multiplyElementwise(y1p_mat_,safe_real_,p1_mat_);
    //p2 = y2./my;
    agile::multiplyElementwise(y2p_mat_,safe_real_,p2_mat_);

    /*
    % update dual (q): (zeta1,zeta2,zeta3) = E(v); v: (zeta4,zeta5) = -p - div q
    zeta1 = dx(vb1);
    zeta2 = dy(vb2);
    zeta3 = (dy(vb1) + dx(vb2))/2;
    z1 = q1 + sigma * zeta1;
    z2 = q2 + sigma * zeta2;
    z3 = q3 + sigma * zeta3;
    mz = proj2(z1,z2,z3);
    q1 = z1./mz;
    q2 = z2./mz;
    q3 = z3./mz;
    */

    agile::diff(1,vb1_mat_,zeta1_mat_);
    agile::diff(2,vb2_mat_,zeta2_mat_);
    //dx(vb2)
    agile::diff(1,vb2_mat_,safe_);
    //dx(vb1)
    agile::diff(2,vb1_mat_,safe2_);
    agile::addMatrix(safe_,safe2_,safe_);
    //zeta3 = (dy(vb1) + dx(vb2))/2;
    agile::scale(typename TType::value_type(0.5),safe_,zeta3_mat_);

    //zx = qx + sigma * zetax;
    agile::scale(sigma,zeta1_mat_,zeta1_mat_);
    agile::scale(sigma,zeta2_mat_,zeta2_mat_);
    agile::scale(sigma,zeta3_mat_,zeta3_mat_);
    agile::addMatrix(q1_mat_,zeta1_mat_,zeta1_mat_);
    agile::addMatrix(q2_mat_,zeta2_mat_,zeta2_mat_);
    agile::addMatrix(q3_mat_,zeta3_mat_,zeta3_mat_);

    //proj2  = @(x1,x2,x3) max(1,sqrt(x1.^2+x2.^2+2*x3.^2)/gamma); % on C^2_gamma
    agile::multiplyElementwise(zeta1_mat_,zeta1_mat_,zeta1_mat_);
    agile::multiplyElementwise(zeta2_mat_,zeta2_mat_,zeta2_mat_);
    agile::multiplyElementwise(zeta3_mat_,zeta3_mat_,zeta3_mat_);
    agile::scale(typename TType::value_type(2),zeta3_mat_,zeta3_mat_);//2*x3.^2
    agile::addMatrix(zeta1_mat_,zeta2_mat_,zeta2_mat_);
    agile::addMatrix(zeta2_mat_,zeta3_mat_,zeta3_mat_);
    agile::real(zeta3_mat_,safe_real_);
    agile::sqrt(safe_real_,safe_real_);
    float gamma = (typename TType::value_type(2))*beta;   //gamma = 2*beta
    if(gamma == 0)
      gamma=1;
    agile::scale(one_/gamma,safe_real_,safe_real_);

    agile::max(this->ones_mat_nxnx_,safe_real_,safe_real_); //safe_real_==mz

    //qx = zx./mz;
    agile::multiplyElementwise(zeta1_mat_,safe_real_,q1_mat_);
    agile::multiplyElementwise(zeta2_mat_,safe_real_,q2_mat_);
    agile::multiplyElementwise(zeta3_mat_,safe_real_,q3_mat_);


    /*
    % update primal (du,dc): F'*F du  + F(u) - div(p)
    [Mu,Mc] = M(ub,cb);
    eta3 = Mu + rhsu + dxt(p1)+dyt(p2); % primal variable: image
    eta4 = Mc + rhsc + alpha*(c+cb);    % primal variable: coils
    uold = du;
    cold = dc;
    du   = du - tau * eta3;
    dc   = dc - tau * eta4;
    */

    //[Mu,Mc] = M(ub,cb);
    IRGN<TType>::ApplyM(&Mu_mat_,Mc_mat_,&ub_mat_,cb_mat_);

    //eta3 = Mu + rhsu + dxt(p1)+dyt(p2); % primal variable: image
    agile::difftrans(1,p1_mat_,safe_);
    agile::difftrans(2,p2_mat_,safe2_);
    agile::addMatrix(safe_,safe2_,safe_);
    agile::addMatrix(rhs[0],safe_,safe_);
    agile::addMatrix(Mu_mat_,safe_,eta3_mat_);
    //uold = du;
    agile::copy(*du, du_old_mat_);
    //du   = du - tau * eta3;
    agile::scale(tau,eta3_mat_,safe_);
    agile::subMatrix(*du,safe_,*du);

    //update primal leading points - ub
    agile::scale(two_,*du,safe_);
    agile::subMatrix(safe_,du_old_mat_,ub_mat_);

    //eta4=Mc+rhsc+aspha*(c+cb)           primal variable: coils

/*    if((maxits == 128) && ((iter == 17)|| (iter == 18)))     //check de 17te  - bei 18 gehts dann nimma...
    {
      std::vector<TType> data_host;
      cb_mat_[1].copyToHost(data_host);
      std::cout<<"\n cb_mat_[1]: "<<iter;
      output("  ",10,10, data_host);
    }
*/
    for(unsigned i=0; i<IRGN<TType>::num_coils_;i++)
    {
      agile::addMatrix(c[i],cb_mat_[i],safe_);

      agile::scale(alpha,safe_,safe_);
      agile::addMatrix(rhs[i+1],safe_,safe_);
      agile::addMatrix(Mc_mat_[i],safe_,eta4_mat_[i]);
/*
      if((i == i)&& (maxits == 80) && ((iter == 17)|| (iter == 18)))     //check de 17te  - bei 18 gehts dann nimma...
      {
        std::vector<TType> data_host;
        Mc_mat_[1].copyToHost(data_host);
        std::cout<<"\n drinnen - Mc_mat_[1]: "<<iter;
        output("  ",10,10, data_host);

        safe_.copyToHost(data_host);
        std::cout<<"\n drinnen - safe_: "<<iter;
        output("  ",10,10, data_host);


        eta4_mat_[1].copyToHost(data_host);
        std::cout<<"\n drinnen - eta4_mat_[1]: "<<iter;
        output("  ",10,10, data_host);

      }
*/

      agile::copy(dc[i], dc_old_mat_[i]);

      agile::scale(tau,eta4_mat_[i],safe_);
      agile::subMatrix(dc[i],safe_,dc[i]);

      //update primal leading points - cb
      agile::scale(two_,dc[i],safe_);
      agile::subMatrix(safe_,dc_old_mat_[i],cb_mat_[i]);
    }

    /*
    % update primal (v): (zeta4,zeta5) = -p - div q
    zeta4 = dxt(q1) + dyt(q3) - p1;
    zeta5 = dxt(q3) + dyt(q2) - p2;
    v1old = v1;
    v2old = v2;
    v1 = v1 - tau * zeta4;
    v2 = v2 - tau * zeta5;
    */
    //zeta4 = dxt(q1) + dyt(q3) - p1;
    agile::difftrans(1,q1_mat_,safe_);
    agile::difftrans(2,q3_mat_,safe2_);
    agile::addMatrix(safe_,safe2_,safe_);
    agile::subMatrix(safe_,p1_mat_,zeta4_mat_);
    //zeta5 = dxt(q3) + dyt(q2) - p2;
    agile::difftrans(1,q3_mat_,safe_);
    agile::difftrans(2,q2_mat_,safe2_);
    agile::addMatrix(safe_,safe2_,safe_);
    agile::subMatrix(safe_,p2_mat_,zeta5_mat_);
    //v1old = v1;   v2old = v2;
    agile::copy(v1_mat_,v1_old_mat_);
    agile::copy(v2_mat_,v2_old_mat_);
    //v1 = v1 - tau * zeta4;
    agile::scale(tau,zeta4_mat_,safe_);
    agile::subMatrix(v1_mat_,safe_,v1_mat_);
    //v2 = v2 - tau * zeta5;
    agile::scale(tau,zeta5_mat_,safe_);
    agile::subMatrix(v2_mat_,safe_,v2_mat_);
    //vb1 = 2*v1 - v1old;
    agile::scale(two_,v1_mat_,safe_);
    agile::subMatrix(safe_,v1_old_mat_,vb1_mat_);
    //vb2 = 2*v2 - v2old;
    agile::scale(two_,v2_mat_,safe_);
    agile::subMatrix(safe_,v2_old_mat_,vb2_mat_);
  }
}


#endif // AGILE_IRGN_TGVSOLVE_HPP
