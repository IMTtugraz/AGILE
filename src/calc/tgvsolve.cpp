#include "agile/calc/tgvsolve.hpp"

//===========================================================================
//
//              M E T H O D E S
//
//===========================================================================
namespace agile
{
  // ---------------------------------------------------------------
  //! \brief Init()
  //!  special initialize for class IRGN.
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void TGVSolve<TType, TType2>::Init()
  {
    norm_=0;
    one_ =1;
    two_ = 2;

    unsigned int num_rows = IRGN<TType,TType2>::num_rows_;
    unsigned int num_columns = IRGN<TType,TType2>::num_columns_;

    cb_mat_ = new GPUMatrix<TType> [IRGN<TType,TType2>::num_coils_];
    Mc_mat_ = new GPUMatrix<TType> [IRGN<TType,TType2>::num_coils_];
    eta4_mat_ = new GPUMatrix<TType> [IRGN<TType,TType2>::num_coils_];
    dc_old_mat_ = new GPUMatrix<TType> [IRGN<TType,TType2>::num_coils_];


    //initialize two-dimensional matrix
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
    safe3_.resize(num_rows, num_columns);
    p1_mat_.resize(num_rows, num_columns);
    p2_mat_.resize(num_rows, num_columns);

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
    for(unsigned i=0; i<IRGN<TType,TType2>::num_coils_;++i)
    {
      cb_mat_[i].resize(num_rows, num_columns);  //achtung - size nxnx
      Mc_mat_[i].resize(num_rows, num_columns);
      eta4_mat_[i].resize(num_rows, num_columns);
      dc_old_mat_[i].resize(num_rows, num_columns);
    }
  }


  // ---------------------------------------------------------------
  //! \brief TGVSolve()
  //!  calculates TGVSolve
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void TGVSolve<TType, TType2>::Solve(const GPUMatrix<TType>* u, const GPUMatrix<TType>* c,
                                 const GPUMatrix<TType>* rhs, const GPUMatrix<TType>* u0,
                                 unsigned maxits, typename to_real_type<TType>::type alpha,
                                 typename to_real_type<TType>::type beta,
                                 GPUMatrix<TType>* du, GPUMatrix<TType>* dc)
  {

    // initialize iterates du, dc and ub, cb
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,*du);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,ub_mat_);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,p1_mat_);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,p2_mat_);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,v1_mat_);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,v2_mat_);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,vb1_mat_);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,vb2_mat_);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,q1_mat_);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,q2_mat_);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,q3_mat_);
    for(unsigned i=0; i < IRGN<TType,TType2>::num_coils_;++i)
    {
      agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,dc[i]);
      agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,cb_mat_[i]);
    }

    typename agile::to_real_type<TType>::type L;
    L = IRGN<TType,TType2>::calcLipschitz();

    typename agile::to_real_type<TType>::type tau;
    tau = 1/std::sqrt(12+L);

    typename agile::to_real_type<TType>::type sigma;
    sigma = 1/std::sqrt(12+L);

    /*IRGN<TType>::timer_value = IRGN<TType>::Timer.stop();
    IRGN<TType>::timer_value = IRGN<TType>::timer_value/1000;
    std::cout << "\n TGVSolve bis iter: " << std::setprecision(3)<< IRGN<TType>::timer_value << "[s]  ";
    IRGN<TType>::Timer.start();*/

    for (unsigned iter = 1; iter <= maxits; ++iter)
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

      agile::sqrt(safe_,safe_);
      if(beta == 0)
        beta=std::numeric_limits<typename TType::value_type>::min();
      agile::scale(one_/beta,safe_,safe_);

      agile::max(this->ones_complex_mat_nxns_,safe_,safe_); //safe_==my

      //p1 = y1./my;
      agile::divideElementwise(y1p_mat_,safe_,p1_mat_);
      //p2 = y2./my;
      agile::divideElementwise(y2p_mat_,safe_,p2_mat_);

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
      agile::multiplyElementwise(zeta1_mat_,zeta1_mat_,safe_);
      agile::multiplyElementwise(zeta2_mat_,zeta2_mat_,safe2_);
      agile::multiplyElementwise(zeta3_mat_,zeta3_mat_,safe3_);
      agile::scale(typename TType::value_type(2),safe3_,safe3_);//2*x3.^2
      agile::addMatrix(safe_,safe2_,safe2_);
      agile::addMatrix(safe2_,safe3_,safe_);

      agile::sqrt(safe_,safe_);

      typename TType::value_type gamma = (typename TType::value_type(2))*beta;   //gamma = 2*beta
      if(gamma == 0)
        gamma=std::numeric_limits<typename TType::value_type>::min();
      agile::scale(one_/gamma,safe_,safe_);

      agile::max(this->ones_complex_mat_nxns_,safe_,safe_); //safe_==mz

      //qx = zx./mz;
      agile::divideElementwise(zeta1_mat_,safe_,q1_mat_);
      agile::divideElementwise(zeta2_mat_,safe_,q2_mat_);
      agile::divideElementwise(zeta3_mat_,safe_,q3_mat_);


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
      IRGN<TType,TType2>::ApplyM(&Mu_mat_,Mc_mat_,&ub_mat_,cb_mat_);

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

      for(unsigned i=0; i<IRGN<TType,TType2>::num_coils_;++i)
      {
        agile::addMatrix(c[i],cb_mat_[i],safe_);

        agile::scale(alpha,safe_,safe_);
        agile::addMatrix(rhs[i+1],safe_,safe_);
        agile::addMatrix(Mc_mat_[i],safe_,eta4_mat_[i]);

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

          /*
    IRGN<TType>::timer_value = IRGN<TType>::Timer.stop();
    IRGN<TType>::timer_value = IRGN<TType>::timer_value/1000;
    std::cout << "\n TGVSolve nach iter: " << std::setprecision(5)<< IRGN<TType>::timer_value << "[s]  "<<IRGN<TType>:: timer_value/60 << "[min]";
    IRGN<TType>::Timer.start();
        */
    /*IRGN<TType>::timer_value = IRGN<TType>::Timer.stop();
    IRGN<TType>::timer_value = IRGN<TType>::timer_value/1000;
    std::cout << "\n TGVSolve: " << std::setprecision(3)<< IRGN<TType>::timer_value << "[s]  ";*/

  }
} //namespace
