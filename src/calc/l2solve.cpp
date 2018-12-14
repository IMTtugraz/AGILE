#include "agile/calc/l2solve.hpp"

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
  void L2Solve<TType, TType2>::Init()
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
    Mu_mat_.resize(num_rows, num_columns);
    eta3_mat_.resize(num_rows, num_columns);
    du_old_mat_.resize(num_rows, num_columns);
    safe_.resize(num_rows, num_columns);

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
  //! \brief L2solve()
  //!  calculates L2solve
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void L2Solve<TType, TType2>::Solve(const GPUMatrix<TType>* u, const GPUMatrix<TType>* c,
                                 const GPUMatrix<TType>* rhs, const GPUMatrix<TType>* u0,
                                 unsigned maxits, typename TType::value_type alpha, typename TType::value_type beta,
                                 GPUMatrix<TType>* du, GPUMatrix<TType>* dc)
  {

    // initialize iterates du, dc and ub, cb
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,*du);
    agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,ub_mat_);
    for(unsigned i=0; i<IRGN<TType,TType2>::num_coils_;++i)
    {
      agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,dc[i]);
      agile::copy(IRGN<TType,TType2>::zero_complex_mat_nxns_,cb_mat_[i]);
    }

    typename agile::to_real_type<TType>::type L;
    L = IRGN<TType,TType2>::calcLipschitz();

    typename agile::to_real_type<TType>::type tau;
    tau = 1/std::sqrt(8+L);

    /*IRGN<TType>::timer_value = IRGN<TType>::Timer.stop();
    IRGN<TType>::timer_value = IRGN<TType>::timer_value/1000;
    std::cout << "\n L2Solve bis iter: " << std::setprecision(3)<< IRGN<TType>::timer_value << "[s]  ";
    IRGN<TType>::Timer.start();*/

    for (unsigned iter = 1; iter <= maxits; ++iter)
    {
      IRGN<TType,TType2>::ApplyM(&Mu_mat_,Mc_mat_,&ub_mat_,cb_mat_);

      //eta3= Mu+rhsu+beta*(u+ub-u0)        primal variable: image
      agile::addMatrix(*u,ub_mat_,safe_);
      agile::subMatrix(safe_,*u0,safe_);
      agile::scale(beta,safe_,safe_);
      agile::addMatrix(rhs[0],safe_,safe_);
      agile::addMatrix(Mu_mat_,safe_,eta3_mat_);

      agile::copy(*du, du_old_mat_);
      agile::scale(tau,eta3_mat_,safe_);
//      test_max(*du,"du");
//      test_max(safe_,"safe_");
      agile::subMatrix(*du,safe_,*du);
//      test_max(*du,"du");

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

    /*IRGN<TType>::timer_value = IRGN<TType>::Timer.stop();
    IRGN<TType>::timer_value = IRGN<TType>::timer_value/1000;
    std::cout << "\n L2Solve: " << std::setprecision(3)<< IRGN<TType>::timer_value << "[s]  ";*/

  }
} //namespace agile
