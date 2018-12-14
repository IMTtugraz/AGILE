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
//                    L2SOLVE.HPP
//
//==================================================================================

#ifndef AGILE_IRGN_L2SOLVE_HPP
#define AGILE_IRGN_L2SOLVE_HPP


#include "irgn.hpp"

using agile::GPUMatrix;
using agile::GPUVector;


template <typename TType>
class L2Solve : public IRGN<TType>
{
  public:

    //! \brief Default constructor.
    //!
    //! The default constructor creates an empty L2Solve Class.
    L2Solve()
      : IRGN<TType>()
    {
    }


    //! \brief Constructor.
    //!
    //! The Constructor creates an L2Solve.
    L2Solve(GPUMatrix<TType>* coil, unsigned int num_coils, IRGN_Params param)
      : IRGN<TType>(coil, num_coils, param)
    {
      Init();
    }

    //! \brief Constructor.
    //!
    //! The Constructor creates an initialized L2Solve.
    L2Solve(GPUMatrix<TType>* coil, unsigned int num_coils)
        : IRGN<TType>(coil, num_coils)
    {
      Init();
    }

    //! \brief Destructor.
    virtual ~L2Solve()
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
    GPUMatrix<TType> Mu_mat_;
    GPUMatrix<TType> eta3_mat_;
    GPUMatrix<TType> du_old_mat_;
    GPUMatrix<TType> safe_;

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
void L2Solve<TType>::Init()
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
  Mu_mat_.resize(num_rows, num_columns);
  eta3_mat_.resize(num_rows, num_columns);
  du_old_mat_.resize(num_rows, num_columns);
  safe_.resize(num_rows, num_columns);

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
}

// ---------------------------------------------------------------
//! \brief L2solve()
//!  calculates L2solve
// ---------------------------------------------------------------
template <typename TType>
void L2Solve<TType>::Solve(const GPUMatrix<TType>* u, const GPUMatrix<TType>* c,
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
  for(unsigned i=0; i<IRGN<TType>::num_coils_;i++)
  {
    agile::copy(IRGN<TType>::random2_mat_[i],x2_mat_[i]);  //achtung - size nxnx - zero_complex_mat_nxnx_

    agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,dc[i]);  //achtung - size nxnx - zero_complex_mat_nxnx_
    agile::copy(IRGN<TType>::zero_complex_mat_nxnx_,cb_mat_[i]);  //achtung - size nxnx - zero_complex_mat_nxnx_
  }

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

  for (unsigned iter = 1; iter <= maxits; iter++)
  {
    ApplyM(&Mu_mat_,Mc_mat_,&ub_mat_,cb_mat_);

    //eta3= Mu+rhsu+beta*(u+ub-u0)        primal variable: image
    agile::addMatrix(*u,ub_mat_,safe_);
    agile::subMatrix(safe_,*u0,safe_);
    agile::scale(beta,safe_,safe_);
    agile::addMatrix(rhs[0],safe_,safe_);
    agile::addMatrix(Mu_mat_,safe_,eta3_mat_);

    agile::copy(*du, du_old_mat_);

    agile::scale(tau,eta3_mat_,safe_);
    agile::subMatrix(*du,safe_,*du);

    //update primal leading points - ub
    agile::scale(two_,*du,safe_);
    agile::subMatrix(safe_,du_old_mat_,ub_mat_);


    //eta4=Mc+rhsc+aspha*(c+cb)           primal variable: coils
    for(unsigned i=0; i<IRGN<TType>::num_coils_;i++)
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
}


#endif // AGILE_IRGN_L2SOLVE_HPP
