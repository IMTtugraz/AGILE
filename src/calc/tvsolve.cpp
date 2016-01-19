#include "agile/calc/tvsolve.hpp"

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
  void TVSolve<TType, TType2>::Init()
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
    p1_mat_.resize(num_rows, num_columns);
    p2_mat_.resize(num_rows, num_columns);

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
  //! \brief TVSolve()
  //!  calculates TVSolve
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void TVSolve<TType, TType2>::Solve(const GPUMatrix<TType>* u, const GPUMatrix<TType>* c,
                                 const GPUMatrix<TType>* rhs, const GPUMatrix<TType>* u0,
                                 unsigned maxits, typename TType::value_type alpha, typename TType::value_type beta,
                                 GPUMatrix<TType>* du, GPUMatrix<TType>* dc)
  {

    // initialize iterates du, dc and ub, cb
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,*du);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,ub_mat_);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,p1_mat_);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,p2_mat_);
    for(unsigned i=0; i<IRGN<TType,TType2>::num_coils_;++i)
    {
      agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,dc[i]);
      agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,cb_mat_[i]);
    }

    typename agile::to_real_type<TType>::type L;
    L = IRGN<TType,TType2>::calcLipschitz();

    typename agile::to_real_type<TType>::type tau;
    tau = 1/std::sqrt(10+L);

    typename agile::to_real_type<TType>::type sigma;
    sigma = 1/std::sqrt(10+L);

    //std::cout<<"\n sigma: "<<sigma;

          /*
    IRGN<TType>::timer_value = IRGN<TType>::Timer.stop();
    IRGN<TType>::timer_value = IRGN<TType>::timer_value/1000;
    std::cout << "\nTVSolve bis iter: " << std::setprecision(5)<< IRGN<TType>::timer_value << "[s]  "<<IRGN<TType>:: timer_value/60 << "[min]";
    IRGN<TType>::Timer.start();
        */

    /*IRGN<TType>::timer_value = IRGN<TType>::Timer.stop();
    IRGN<TType>::timer_value = IRGN<TType>::timer_value/1000;
    std::cout << "\nTVSolve bis iter: " << std::setprecision(3)<< IRGN<TType>::timer_value << "[s]  ";
    IRGN<TType>::Timer.start();*/


    for (unsigned iter = 1; iter <= maxits; ++iter)
    {
  /*
      % update dual (p): (eta1,eta2) = nabla(u+du-u0)
      ukp = ub + u - u0;
      eta1 = dx(ukp);
      eta2 = dy(ukp);
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
      //eta2 = dy(ukp);
      agile::diff(2,ukp_mat_,eta2_mat_);

      //y1 = p1 + sigma * eta1;
      agile::scale(sigma,eta1_mat_,eta1_mat_);
      agile::addMatrix(p1_mat_,eta1_mat_,y1p_mat_);
      //y2 = p2 + sigma * eta2;
      agile::scale(sigma,eta2_mat_,eta2_mat_);
      agile::addMatrix(p2_mat_,eta2_mat_,y2p_mat_);

      //proj = @(y1,y2) max(1,sqrt(y1*conj(y1)+y2*conj(y2))/beta);
      agile::multiplyConjElementwise(y1p_mat_,y1p_mat_,y1p_mat_);
      agile::multiplyConjElementwise(y2p_mat_,y2p_mat_,y2p_mat_);

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

    }

/*
    IRGN<TType>::timer_value = IRGN<TType>::Timer.stop();
    IRGN<TType>::timer_value = IRGN<TType>::timer_value/1000;
    std::cout << "\nTVSolve nach iter: " << std::setprecision(5)<< IRGN<TType>::timer_value << "[s]  "<<IRGN<TType>:: timer_value/60 << "[min]";
    IRGN<TType>::Timer.start();
        */

    /*IRGN<TType>::timer_value = IRGN<TType>::Timer.stop();
    IRGN<TType>::timer_value = IRGN<TType>::timer_value/1000;
    std::cout << "\nTVSolve: " << std::setprecision(3)<< IRGN<TType>::timer_value << "[s]  ";*/

  }
}  //namespace agile
