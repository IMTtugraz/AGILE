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
//                    IRGN.CPP
//
//==================================================================================


#include "irgn.hpp"

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
void IRGN<TType>::Init()
{
  init_matrixlog(myfile,"matrixlog.txt");

  cufftPlan2d(&fftplan_, num_rows_, num_columns_, CUFFT_C2C);

  typename agile::to_real_type<TType>::type zero=0;
  typename agile::to_real_type<TType>::type one=1;
  TType zero_complex;
  zero_complex.real(0);
  zero_complex.imag(0);
  TType one_complex;     //val = 1+0i
  one_complex.real(1);
  one_complex.imag(0);

  for(int i=0; i<num_rows_*num_rows_; i++)
  {
    zero_vec_cpu_nxnx_.push_back(zero);
    zero_complex_vec_cpu_nxnx_.push_back(zero_complex);
    ones_vec_cpu_nxnx_.push_back(one);
    ones_complex_vec_cpu_nxnx_.push_back(one_complex);
  }

  for(int i=0; i<num_rows_*num_columns_; i++)
  {
    zero_vec_cpu_nxns_.push_back(zero);
    zero_complex_vec_cpu_nxns_.push_back(zero_complex);
    ones_vec_cpu_nxns_.push_back(one);
    ones_complex_vec_cpu_nxns_.push_back(one_complex);
  }

  c_mat_ = new GPUMatrix<TType> [num_coils_];
  cw_mat_ = new GPUMatrix<TType> [num_coils_];
  cw_bar_mat_ = new GPUMatrix<TType> [num_coils_];


  //----- Helper Matrix: -----
  //nx nx  - matrix
  zero_mat_nxnx_.assignFromHost(num_rows_, num_rows_, &zero_vec_cpu_nxnx_[0]);
  ones_mat_nxnx_.assignFromHost(num_rows_, num_rows_, &ones_vec_cpu_nxnx_[0]);
  zero_complex_mat_nxnx_.assignFromHost(num_rows_, num_rows_, &zero_complex_vec_cpu_nxnx_[0]);
  ones_complex_mat_nxnx_.assignFromHost(num_rows_, num_rows_, &ones_complex_vec_cpu_nxnx_[0]); //val = 1+0i

  //nx ns  - matrix
  zero_complex_mat_nxns_.assignFromHost(num_rows_, num_columns_, &zero_complex_vec_cpu_nxns_[0]);
  ones_complex_mat_nxns_.assignFromHost(num_rows_, num_columns_, &ones_complex_vec_cpu_nxns_[0]);


  random2_mat_ = new GPUMatrix<TType> [IRGN<TType>::num_coils_];
  //generate complex random - vector:
  std::complex<float> val_complex;

  for (unsigned coil = 0; coil < IRGN<TType>::num_coils_; ++coil)
    for (unsigned row = 0; row < IRGN<TType>::num_rows_; ++row)
      for (unsigned column = 0; column < IRGN<TType>::num_rows_; ++column)
      {
        val_complex.imag(randomcalc(column*row));
        val_complex.real(randomcalc(column));
        x2_.push_back(val_complex);
        if(coil==0)
        {
          val_complex=(randomcalc(column*row),randomcalc(column));
          x1_.push_back(val_complex);
        }
      }



  //----- Set Init Matrics -----
  u_mat_.resize(num_rows_, num_rows_);
  u0_mat_.resize(num_rows_, num_rows_);
  u_bar_mat_.resize(num_rows_, num_rows_);
  w_mat_.resize(num_rows_, num_columns_);
  us_mat_.resize(num_rows_, num_rows_);
  image_mat_.resize(num_rows_, num_rows_);

  pattern_.resize(num_rows_, num_columns_);
  pattern_complex_.resize(num_rows_, num_columns_);

  random1_mat_.assignFromHost(num_rows_, num_rows_, &x1_[0]);
  for(unsigned i=0; i<num_coils_;i++)
  {
    random2_mat_[i].assignFromHost(num_rows_, num_rows_, &x2_[num_rows_ * num_rows_ * i]);
    c_mat_[i].resize(num_rows_, num_rows_);
    cw_bar_mat_[i].resize(num_rows_, num_rows_);
    cw_mat_[i].resize(num_rows_, num_rows_);
  }

  agile::copy(ones_complex_mat_nxnx_,u_mat_);

  for(unsigned i=0; i<num_coils_;i++)
    agile::copy(zero_complex_mat_nxnx_,c_mat_[i]);

  //u0 value calculation
  typename agile::to_real_type<TType>::type u0=0;
  agile::scale(u0,u_mat_,u0_mat_);

  //generate pattern
  agile::pattern(coil_[0],pattern_);

  agile::multiplyElementwise(pattern_, ones_complex_mat_nxns_, pattern_complex_);
}

// ---------------------------------------------------------------
//! \brief randomcalc()
//!  random value generator
// ---------------------------------------------------------------
template <typename TType>
float IRGN<TType>::randomcalc(int i)
{
  //rand(n,m,nc)
  srand ( time(NULL)*i );    //init random-generator
  int val = rand() % 100 + 1;
  float one=1;

  if (val==0)
    val=2;

  return one/val;
}

// ---------------------------------------------------------------
//! \brief HighFreqPenalty()
//!  calculates the High Freq Penalties
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::HighFreqPenalty()
{
  typename agile::to_real_type<TType>::type alpha=220;
  typename agile::to_real_type<TType>::type two=2;
  typename agile::to_real_type<TType>::type val=-8;

  GPUVector<typename agile::to_real_type<TType>::type> linspace_x_vec(num_rows_,NULL);
  GPUVector<typename agile::to_real_type<TType>::type> linspace_y_vec(num_rows_,NULL);
  GPUMatrix<typename agile::to_real_type<TType>::type> xi_mat(num_rows_,num_rows_,NULL);
  GPUMatrix<typename agile::to_real_type<TType>::type> eta_mat(num_rows_,num_rows_,NULL);

  agile::linspace(linspace_x_vec, -0.5f , 0.5f);
  agile::linspace(linspace_y_vec, -0.5f , 0.5f);
  agile::meshgrid(xi_mat,eta_mat,linspace_x_vec,linspace_y_vec);

  // w = (1+220*(xi.²+eta.²)).^(-8)
  agile::pow(two,xi_mat,xi_mat);
  agile::pow(two,eta_mat,eta_mat);
  agile::addMatrix(xi_mat,eta_mat,w_mat_);
  agile::scale(alpha,w_mat_,w_mat_);
  agile::addMatrix(ones_mat_nxnx_,w_mat_,w_mat_);
  agile::pow(val,w_mat_,w_mat_);
}

// ---------------------------------------------------------------
//! \brief Normalize()
//!  Normalize
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::Normalize()
{
  typename agile::to_real_type<TType>::type alpha=100;

  dscale_ = 0;

  GPUMatrix<typename agile::to_real_type<TType>::type> abs_coil;
  abs_coil.resize(num_rows_,num_columns_);

  for(unsigned i=0; i<num_coils_;i++)
  {
    agile::absMatrix(coil_[i],abs_coil);
    dscale_ = dscale_ + std::pow(agile::norm2(abs_coil),2);
  }
  dscale_ = std::sqrt(dscale_);
  AGILE_ASSERT(dscale_ == 0,
                StandardException::ExceptionMessage(
             "divide 0 - error - Normalize"));

  dscale_ = alpha/dscale_;

  for(unsigned i=0; i<num_coils_;i++)
  {
    agile::scale(dscale_,coil_[i],coil_[i]);
  }
}


// ---------------------------------------------------------------
//! \brief CenterdFFT(GPUMatrix<TType>* fft_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
//!  calculates the Centerd FFT for given matrix
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::CenterdFFTpattern(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
{
  typename TType::value_type val_sqrt;

  this->CopyMatrixZ(in_mat,out_mat,num_z);

  for(unsigned i=0; i<num_z;i++)
  {
    // fft
    agile::fftshift(out_mat[i]);
    cufftExecC2C(fftplan_,
                 (cufftComplex*)out_mat[i].data(),
                 (cufftComplex*)out_mat[i].data(),
                 CUFFT_FORWARD);
    agile::ifftshift(out_mat[i]);

    // fft./sqrt(num_columns*num_rows)
    AGILE_ASSERT(std::sqrt(out_mat[i].getNumRows() * out_mat[i].getNumColumns()) == 0,
                  StandardException::ExceptionMessage(
               "CenterdFFTpattern - divide 0 - error - check num_columns and num_rows"));

    val_sqrt = typename TType::value_type(1)/std::sqrt(out_mat[i].getNumRows() * out_mat[i].getNumColumns());

    agile::scale(val_sqrt, out_mat[i], out_mat[i]);


    agile::multiplyElementwise(pattern_complex_,out_mat[i],out_mat[i]);

  }
}


// ---------------------------------------------------------------
//! \brief CenterdIFFT(GPUMatrix<TType>* ifft_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
//!  calculates the Centerd IFFT for given matrix
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::CenterdIFFTpattern(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
{
  typename TType::value_type val_sqrt;

  this->CopyMatrixZ(in_mat,out_mat,num_z);

  for(unsigned i=0; i<num_z;i++)
  {
    agile::multiplyElementwise(pattern_complex_,out_mat[i],out_mat[i]);

    // ifft
    agile::fftshift(out_mat[i]);

    cufftExecC2C(fftplan_,
                 (cufftComplex*)out_mat[i].data(),
                 (cufftComplex*)out_mat[i].data(),
                 CUFFT_INVERSE);
    agile::ifftshift(out_mat[i]);

    // fft./sqrt(num_columns*num_rows)
    AGILE_ASSERT(std::sqrt(out_mat[i].getNumRows() * out_mat[i].getNumColumns()) == 0,
                  StandardException::ExceptionMessage(
               "CenterdIFFTpattern - divide 0 - error - check num_columns and num_rows"));

    val_sqrt = typename TType::value_type(1)/std::sqrt(out_mat[i].getNumRows() * out_mat[i].getNumColumns());

    agile::scale(val_sqrt, out_mat[i], out_mat[i]);

  }
}

// ---------------------------------------------------------------
//! \brief CenterdFFT(GPUMatrix<TType>* fft_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
//!  calculates the Centerd FFT for given matrix
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::CenterdFFT(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
{
  typename TType::value_type val_sqrt;

  this->CopyMatrixZ(in_mat,out_mat,num_z);

  for(unsigned i=0; i<num_z;i++)
  {

    // fft
    agile::fftshift(out_mat[i]);
    cufftExecC2C(fftplan_,
                 (cufftComplex*)out_mat[i].data(),
                 (cufftComplex*)out_mat[i].data(),
                 CUFFT_FORWARD);
    agile::ifftshift(out_mat[i]);

    // fft./sqrt(num_columns*num_rows)
    AGILE_ASSERT(std::sqrt(out_mat[i].getNumRows() * out_mat[i].getNumColumns()) == 0,
                  StandardException::ExceptionMessage(
               "CenterdFFT - divide 0 - error - check num_columns and num_rows"));

    val_sqrt = typename TType::value_type(1)/std::sqrt(out_mat[i].getNumRows() * out_mat[i].getNumColumns());

    agile::scale(val_sqrt, out_mat[i], out_mat[i]);

  }
}


// ---------------------------------------------------------------
//! \brief CenterdIFFT(GPUMatrix<TType>* ifft_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
//!  calculates the Centerd IFFT for given matrix
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::CenterdIFFT(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
{
  typename TType::value_type val_sqrt;

  this->CopyMatrixZ(in_mat,out_mat,num_z);

  for(unsigned i=0; i<num_z;i++)
  {
    // ifft
    agile::fftshift(out_mat[i]);

    cufftExecC2C(fftplan_,
                 (cufftComplex*)out_mat[i].data(),
                 (cufftComplex*)out_mat[i].data(),
                 CUFFT_INVERSE);
    agile::ifftshift(out_mat[i]);

    // fft./sqrt(num_columns*num_rows)
    AGILE_ASSERT(std::sqrt(out_mat[i].getNumRows() * out_mat[i].getNumColumns()) == 0,
                  StandardException::ExceptionMessage(
               "CenterdIFFT - divide 0 - error - check num_columns and num_rows"));

    val_sqrt = typename TType::value_type(1)/std::sqrt(out_mat[i].getNumRows() * out_mat[i].getNumColumns());

    agile::scale(val_sqrt, out_mat[i], out_mat[i]);

  }
}


// ---------------------------------------------------------------
//! \brief ApplyW(GPUMatrix<TType>* ifft_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
//!  calculates the w-scaled Centerd IFFT for given matrix
//!  W = @(x) cifft(w.*x)
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::ApplyW(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
{
  for(unsigned i=0; i<num_z;i++)
  {
    agile::multiplyElementwise(w_mat_, in_mat[i], out_mat[i]);
  }

  this->CenterdIFFT(out_mat,out_mat,num_z);
}

// ---------------------------------------------------------------
//! \brief ApplyWH(GPUMatrix<TType>* ifft_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
//!  calculates the wbar-scaled Centerd IFFT for given matrix
//!  WH = @(x) wbar.*cfft(x)
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::ApplyWH(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
{
  CenterdFFT(&in_mat[0],&out_mat[0],num_z);

  //wbar=conj(w)   - w is real -> no conj needed
  for(unsigned i=0; i<num_z;i++)
    //agile::multiplyConjElementwise(w_mat_,out_mat[i],out_mat[i]);
    agile::multiplyElementwise(w_mat_,out_mat[i],out_mat[i]);
}


// ---------------------------------------------------------------
//! \brief Iteration()
//!  Calculates Gauss-Newton Iteration
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::Iteration()
{
  float nr_k=0;
  unsigned tvits = this->tvits_;
  float alpha = this->alpha0_;
  float beta = this->beta0_;

  GPUMatrix<TType>* res_mat;
  res_mat = new GPUMatrix<TType> [num_coils_];

  GPUMatrix<TType>* rhs_mat;
  rhs_mat = new GPUMatrix<TType> [num_coils_+1];

  GPUMatrix<TType>* dc_mat;
  dc_mat = new GPUMatrix<TType> [num_coils_];

  GPUMatrix<TType> du_mat(num_rows_, num_rows_, NULL);
  GPUMatrix<TType> safe_mat(num_rows_, num_rows_, NULL);

  for(unsigned i=0; i<num_coils_;i++)
  {
    dc_mat[i].resize(num_rows_, num_rows_);
    res_mat[i].resize(num_rows_, num_rows_);
    rhs_mat[i+1].resize(num_rows_, num_rows_);
  }
  rhs_mat[0].resize(num_rows_, num_rows_);

  agile::copy(zero_complex_mat_nxnx_,rhs_mat[0]);

  //---- Start Iteration ----
  for(int k=1; k<=this->maxit_; k++)
  {
    // apply inverse weight to sensitivity iterates
    this->ApplyW(c_mat_,cw_mat_,num_coils_);

    agile::conjMatrix(u_mat_,u_bar_mat_);
    for(unsigned i=0; i<num_coils_;i++)
    {
      agile::conjMatrix(cw_mat_[i],cw_bar_mat_[i]);

      //calculate residual
      agile::multiplyElementwise(u_mat_,cw_mat_[i],safe_mat);

      CenterdFFTpattern(&safe_mat,&safe_mat,1);

      agile::subMatrix(safe_mat,coil_[i],res_mat[i]);     // UNTERSCHIEDLICHE GROESSEN!!!!
    }

    //Calculate residual norm   nr_k:
    for(unsigned i=0; i<num_coils_;i++)
    {
      nr_k = nr_k + std::pow(agile::norm2(res_mat[i]),2);
    }
    nr_k = std::sqrt(nr_k);
    nr_k_.push_back(nr_k);                                //write residual norm to vector

    std::cout<<"\nIteration "<<k<<": tvits = "<<tvits<<", alpha = "<<alpha<<", beta = "<<beta<<", residual norm = "<<nr_k;

    ApplyDFH(rhs_mat, res_mat);   //result in rhs_mat[ 0 , value[1..nc+1] ]

    //l2sove
    Solve(&u_mat_, c_mat_, rhs_mat, &u0_mat_, tvits, alpha, beta, &du_mat, dc_mat);


    agile::addMatrix(u_mat_,du_mat,u_mat_);

    for(unsigned i=0; i<num_coils_;i++)
      agile::addMatrix(c_mat_[i],dc_mat[i],c_mat_[i]);

    //reduce parameter

    alpha = (this->alpha_min_ >= (alpha_q_*alpha)) ? alpha_min_ : (alpha_q_*alpha);
    beta = (this->beta_min_ >= (beta_q_*beta)) ? beta_min_ : (beta_q_*beta);
    float tvits_help = 1.5 * float(tvits);

    tvits = ((int(tvits_help)) <= this->tvmax_) ? int(tvits_help) : tvmax_;
  }

  delete[] res_mat;
  res_mat = 0;
  delete[] rhs_mat;
  rhs_mat = 0;
  delete[] dc_mat;
  dc_mat = 0;
}


// ---------------------------------------------------------------
//! \brief ApplyDFH()
//!  calculates ApplyDFH
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::ApplyDFH(GPUMatrix<TType>* rhs_mat, const GPUMatrix<TType>* dx)
{
  GPUMatrix<TType>* fdy;
  fdy = new GPUMatrix<TType> [num_coils_];
  for(unsigned i=0; i<num_coils_;i++)
    fdy[i].resize(num_rows_, num_rows_);

  GPUMatrix<TType> safe_mat(num_rows_, num_rows_,NULL);
  GPUMatrix<TType> tmp(num_rows_, num_rows_,NULL);

  agile::copy(zero_complex_mat_nxnx_,tmp);   //set tmp to 0

  CenterdIFFTpattern(dx,fdy,num_coils_);

  for(unsigned i=0; i<num_coils_;i++)
  {
    //fdy.*cconj    (cconj=cw_bar_mat_)
    agile::multiplyElementwise(fdy[i],cw_bar_mat_[i],safe_mat);
    agile::addMatrix(safe_mat,tmp,tmp);

    //fdy.*uconj    (uconj=u_bar_mat_)
    agile::multiplyElementwise(fdy[i],u_bar_mat_,safe_mat);

    this->ApplyWH(&safe_mat,&rhs_mat[i+1],1);
  }
  agile::copy(tmp, rhs_mat[0]);

  delete[] fdy;
  fdy = 0;
}


// ---------------------------------------------------------------
//! \brief ApplyM()
//!  calculates ApplyM
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::ApplyM(GPUMatrix<TType>* gu, GPUMatrix<TType>* gc,
                               const GPUMatrix<TType>* du, const GPUMatrix<TType>* dc)
{
  GPUMatrix<TType> safe_mat(num_rows_, num_rows_,NULL);
  GPUMatrix<TType> fdy(num_rows_, num_rows_,NULL);

  agile::copy(zero_complex_mat_nxnx_, *gu);    //set gu-matrix to zero

  for(unsigned i=0; i<num_coils_;i++)
  {
    //W(dc)
    ApplyW(&dc[i],&safe_mat,1);

    //u*W(dc)
    agile::multiplyElementwise(u_mat_ , safe_mat, safe_mat);
    //c*du
    agile::multiplyElementwise(cw_mat_[i] , *du, fdy);
    //u*W(dc) + c*du
    agile::addMatrix(safe_mat, fdy, fdy);

    //fdy*cconj
    agile::multiplyElementwise(fdy, cw_bar_mat_[i], safe_mat);

    //gu + fdy*cconj
    agile::addMatrix(*gu, safe_mat, *gu);

    //fdy*uconj
    agile::multiplyElementwise(fdy, u_bar_mat_, gc[i]);
  }
  //WH(fdy*uconj)
  ApplyWH(gc,gc,num_coils_);
}


// ---------------------------------------------------------------
//! \brief Postprocess()
//!  calculates postprocess
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::Postprocess()
{
  GPUMatrix<typename agile::to_real_type<TType>::type> cscale_mat(num_rows_, num_rows_, NULL);
  GPUMatrix<typename agile::to_real_type<TType>::type> absval(num_rows_, num_rows_, NULL);

  GPUMatrix<TType> safe_mat(num_rows_, num_rows_, NULL);

  agile::copy(this->zero_mat_nxnx_,cscale_mat);   //init with 0

  for(unsigned i=0; i<num_coils_;i++)
  {
    this->ApplyW(&c_mat_[i],&safe_mat,1);
    agile::absMatrix(safe_mat,absval);
    agile::pow(float(2),absval,absval);
    agile::addMatrix(absval,cscale_mat,cscale_mat);
  }
  agile::sqrt(cscale_mat,cscale_mat);

  agile::multiplyElementwise(u_mat_, cscale_mat, us_mat_);

  agile::scale(dscale_,us_mat_,us_mat_);

  //Postprocess to Image-Data:
  agile::absMatrix(us_mat_,image_mat_);
}


// ---------------------------------------------------------------
//! \brief CopyMatrixZ()
//!  Copy number "num_z" of Matrix "in_mat" to "out_mat"
// ---------------------------------------------------------------
template <typename TType>
void IRGN<TType>::CopyMatrixZ(const GPUMatrix<TType>* in_mat,
                 GPUMatrix<TType>* out_mat,
                 unsigned int num_z)
{
  for(unsigned i=0; i<num_z;++i)
  {
    agile::copy(in_mat[i],out_mat[i]);
  }
}

template <typename TType>
bool IRGN<TType>::irgn_param_test(IRGN_Params &param)
{
  bool ok = true;

  if(param.maxit < 0)
  {
    param.maxit = 1;
    ok = false;
  }

  if((param.tvits < 0) || (param.tvits > param.tvmax))
  {
    param.tvits = 20;
    ok = false;
  }

  if((param.tvmax <= 0) || (param.tvits > param.tvmax))
  {
    param.tvmax = 1000;
    ok = false;
  }

  if(param.alpha_min < 0)
  {
    param.alpha_min = 0;
    ok = false;
  }

  if(param.beta_min < 0)
  {
    param.beta_min = 0;
    ok = false;
  }

  if(param.alpha0 <= 0)
  {
    param.alpha0 = 1;
    ok = false;
  }

  if(param.beta0 <= 0)
  {
    param.beta0 = 1;
    ok = false;
  }

  if(param.alpha_q <= 0)
  {
    param.alpha_q = 0.1;
    ok = false;
  }

  if(param.beta_q <= 0)
  {
    param.beta_q = 0.2;
    ok = false;
  }

  return ok;
}
