//==================================================================================
//
//                    IRGN.CPP
//
//==================================================================================


#include "agile/calc/irgn.hpp"

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
  void IRGN<TType, TType2>::Init()
  {
    //init_matrixlog(myfile,"toolkit_41_matrixlog_TV.txt");

    //Timer.start();

    _break_calc=false;
    fftobj_ = new agile::FFT<TType>(num_rows_, num_columns_);

    typename agile::to_real_type<TType>::type zero=0;
    typename agile::to_real_type<TType>::type one=1;
    TType zero_complex;
    zero_complex.real(0);
    zero_complex.imag(0);
    TType one_complex;     //val = 1+0i
    one_complex.real(1);
    one_complex.imag(0);

    for(int i=0; i<num_rows_*num_columns_; ++i)
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
    zero_mat_nxns_.assignFromHost(num_rows_, num_columns_, &zero_vec_cpu_nxns_[0]);
    ones_mat_nxns_.assignFromHost(num_rows_, num_columns_, &ones_vec_cpu_nxns_[0]);

    //nx ns  - matrix
    zero_complex_mat_nxns_.assignFromHost(num_rows_, num_columns_, &zero_complex_vec_cpu_nxns_[0]);
    ones_complex_mat_nxns_.assignFromHost(num_rows_, num_columns_, &ones_complex_vec_cpu_nxns_[0]);


    random2_mat_ = new GPUMatrix<TType> [num_coils_];
    //generate complex random - vector:
    std::complex<typename TType::value_type> val_complex;

    for (unsigned coil = 0; coil < num_coils_; ++coil)
      for (unsigned row = 0; row < num_rows_; ++row)
        for (unsigned column = 0; column < num_columns_; ++column)
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
    u_mat_.resize(num_rows_, num_columns_);
    u0_mat_.resize(num_rows_, num_columns_);
    u_bar_mat_.resize(num_rows_, num_columns_);
    w_mat_.resize(num_rows_, num_columns_);
    us_mat_.resize(num_rows_, num_columns_);
    image_mat_.resize(num_rows_, num_columns_);
    image_mat_final_.resize(num_rows_, num_columns_);

    random1_mat_.assignFromHost(num_rows_, num_columns_, &x1_[0]);
    for(unsigned i=0; i<num_coils_;++i)
    {
      random2_mat_[i].assignFromHost(num_rows_, num_columns_, &x2_[num_rows_ * num_columns_ * i]);
      c_mat_[i].resize(num_rows_, num_columns_);
      cw_bar_mat_[i].resize(num_rows_, num_columns_);
      cw_mat_[i].resize(num_rows_, num_columns_);
    }

    agile::copy(ones_complex_mat_nxns_,u_mat_);

    for(unsigned i=0; i<num_coils_;++i)
      agile::copy(zero_complex_mat_nxns_,c_mat_[i]);

    //u0 value calculation
    typename agile::to_real_type<TType>::type u0=0;
    agile::scale(u0,u_mat_,u0_mat_);

    fftobj_->calc_pattern(coil_[0]);


    /*
    IRGN<TType>::timer_value = IRGN<TType>::Timer.stop();
    IRGN<TType>::timer_value = timer_value/1000;
    std::cout << "\n\n Init: " << std::setprecision(5)<< IRGN<TType>::timer_value << "[s]  "<<IRGN<TType>:: timer_value/60 << "[min]";
    IRGN<TType>::Timer.start();
    */
  }


  // ---------------------------------------------------------------
  //! \brief randomcalc()
  //!  random value generator
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  typename TType::value_type IRGN<TType, TType2>::randomcalc(int i)
  {
    //rand(n,m,nc)
    srand ( time(NULL)*i );    //init random-generator
    int val = rand() % 100 + 1;
    typename TType::value_type one=1;

    if (val==0)
      val=2;

    //return one/val;
    TType2 zahl=0.5;
    return zahl;
  }

  // ---------------------------------------------------------------
  //! \brief HighFreqPenalty()
  //!  calculates the High Freq Penalties
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::HighFreqPenalty()
  {
    typename agile::to_real_type<TType>::type alpha=220;
    typename agile::to_real_type<TType>::type two=2;
    typename agile::to_real_type<TType>::type val=-8;

    GPUVector<typename agile::to_real_type<TType>::type> linspace_x_vec(num_columns_,0);
    GPUVector<typename agile::to_real_type<TType>::type> linspace_y_vec(num_rows_,0);
    GPUMatrix<typename agile::to_real_type<TType>::type> xi_mat(num_rows_,num_columns_,NULL);
    GPUMatrix<typename agile::to_real_type<TType>::type> eta_mat(num_rows_,num_columns_,NULL);

    agile::linspace(linspace_x_vec, -0.5f , 0.5f);
    agile::linspace(linspace_y_vec, -0.5f , 0.5f);
    agile::meshgrid(xi_mat,eta_mat,linspace_x_vec,linspace_y_vec);

    // w = (1+220*(xi.²+eta.²)).^(-8)
    agile::pow(two,xi_mat,xi_mat);
    agile::pow(two,eta_mat,eta_mat);
    agile::addMatrix(xi_mat,eta_mat,w_mat_);
    agile::scale(alpha,w_mat_,w_mat_);
    agile::addMatrix(ones_mat_nxns_,w_mat_,w_mat_);
    agile::pow(val,w_mat_,w_mat_);
  }

  // ---------------------------------------------------------------
  //! \brief Normalize()
  //!  Normalize
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::Normalize()
  {
    //typename agile::to_real_type<TType>::type alpha=100;
    typename agile::to_real_type<TType>::type alpha=100;
    //alpha=std::numeric_limits<typename TType::value_type>::max();

    dscale_ = 0;

    GPUMatrix<TType2> abs_coil;
    abs_coil.resize(num_rows_,num_columns_);

    for(unsigned i=0; i<num_coils_;++i)
    {
      agile::absMatrix(coil_[i],abs_coil);
      dscale_ = dscale_ + std::pow(agile::norm2(abs_coil),2);
    }
    dscale_ = std::sqrt(dscale_);
    AGILE_ASSERT(dscale_ > 0,
                  StandardException::ExceptionMessage(
               "divide 0 - error - Normalize"));

    dscale_ = alpha/dscale_;

    //std::cout<<"\ndscale_: "<<dscale_;

    for(unsigned i=0; i<num_coils_;++i)
    {
      agile::scale(dscale_,coil_[i],coil_[i]);
    }
  }


  // ---------------------------------------------------------------
  //! \brief CenterdFFT(GPUMatrix<TType>* fft_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
  //!  calculates the Centerd FFT for given matrix
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::CenterdFFTpattern(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
  {
    for(unsigned i=0; i<num_z;++i)
    {
      fftobj_->CenterdFFTpattern(in_mat[i],out_mat[i]);
    }
  }


  // ---------------------------------------------------------------
  //! \brief CenterdIFFTpattern(GPUMatrix<TType>* ifft_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
  //!  calculates the Centerd IFFT for given matrix
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::CenterdIFFTpattern(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
  {
    for(unsigned i=0; i<num_z;++i)
    {
      fftobj_->CenterdIFFTpattern(in_mat[i],out_mat[i]);
    }
  }

  // ---------------------------------------------------------------
  //! \brief CenterdFFT(GPUMatrix<TType>* fft_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
  //!  calculates the Centerd FFT for given matrix
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::CenterdFFT(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
  {
    for(unsigned i=0; i<num_z;++i)
    {
      fftobj_->CenterdFFT(in_mat[i],out_mat[i]);
    }
  }


  // ---------------------------------------------------------------
  //! \brief CenterdIFFT(GPUMatrix<TType>* ifft_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
  //!  calculates the Centerd IFFT for given matrix
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::CenterdIFFT(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
  {
    for(unsigned i=0; i<num_z;++i)
    {
      fftobj_->CenterdIFFT(in_mat[i],out_mat[i]);
    }
  }


  // ---------------------------------------------------------------
  //! \brief ApplyW(GPUMatrix<TType>* ifft_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
  //!  calculates the w-scaled Centerd IFFT for given matrix
  //!  W = @(x) cifft(w.*x)
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::ApplyW(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
  {
    for(unsigned i=0; i<num_z;++i)
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
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::ApplyWH(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z)
  {
    CenterdFFT(&in_mat[0],&out_mat[0],num_z);

    //wbar=conj(w)   - w is real -> no conj needed
    for(unsigned i=0; i<num_z;++i)
      //agile::multiplyConjElementwise(w_mat_,out_mat[i],out_mat[i]);
      agile::multiplyElementwise(w_mat_,out_mat[i],out_mat[i]);
  }


  // ---------------------------------------------------------------
  //! \brief Iteration()
  //!  Calculates Gauss-Newton Iteration
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::Iteration()
  {
    float nr_k=0;
    unsigned tvits = this->tvits_;
    TType2 alpha = this->alpha0_;
    TType2 beta= this->beta0_;

    GPUMatrix<TType>* res_mat;
    res_mat = new GPUMatrix<TType> [num_coils_];

    GPUMatrix<TType>* rhs_mat;
    rhs_mat = new GPUMatrix<TType> [num_coils_+1];

    GPUMatrix<TType>* dc_mat;
    dc_mat = new GPUMatrix<TType> [num_coils_];

    GPUMatrix<TType> du_mat(num_rows_, num_columns_, NULL);
    GPUMatrix<TType> safe_mat(num_rows_, num_columns_, NULL);

    //u_mat_old  and  c_mat_old: if res_norm==nan -> old matrics are taken.
    GPUMatrix<TType>* c_mat_old;
    c_mat_old = new GPUMatrix<TType> [num_coils_];

    GPUMatrix<TType> u_mat_old(num_rows_, num_columns_, NULL);


    for(unsigned i=0; i<num_coils_;++i)
    {
      dc_mat[i].resize(num_rows_, num_columns_);
      res_mat[i].resize(num_rows_, num_columns_);
      rhs_mat[i+1].resize(num_rows_, num_columns_);
      c_mat_old[i].resize(num_rows_, num_columns_);
    }
    rhs_mat[0].resize(num_rows_, num_columns_);

    agile::copy(zero_complex_mat_nxns_,rhs_mat[0]);

    nr_k_.push_back(100);   //start with 100

    std::cout<<std::endl;

    //---- Start Iteration ----
    for(int k=1; k<=this->maxit_; k++)
    {
      // apply inverse weight to sensitivity iterates
      this->ApplyW(c_mat_,cw_mat_,num_coils_);

      agile::conjMatrix(u_mat_,u_bar_mat_);
      for(unsigned i=0; i<num_coils_;++i)
      {
        agile::conjMatrix(cw_mat_[i],cw_bar_mat_[i]);

        //calculate residual
        agile::multiplyElementwise(u_mat_,cw_mat_[i],safe_mat);

        CenterdFFTpattern(&safe_mat,&safe_mat,1);

        agile::subMatrix(safe_mat,coil_[i],res_mat[i]);
      }

      //Calculate residual norm   nr_k:
      for(unsigned i=0; i<num_coils_;++i)
      {
        nr_k = nr_k + std::pow(agile::norm2(res_mat[i]),2);
      }
      nr_k = std::sqrt(nr_k);

      std::cout<<"Iteration "<<k<<": tvits = "<<tvits<<", alpha = "<<alpha<<", beta = "<<beta<<", residual norm = "<<nr_k<<std::endl;

      //break if nr_k is nan; and restore old u_mat and c_mat
      if ((nr_k != nr_k) || (nr_k>(nr_k_.back()+10)) || _break_calc)   //if nan or nr_k is higher than the last nr_k+10
      {
        agile::copy(u_mat_old, u_mat_);
        CopyMatrixZ(c_mat_old, c_mat_, num_coils_);
        std::cerr<<"Break at Iteration "<<k;
        nr_k_.push_back(nr_k);                                //-1 norm to vector for error
        nr_k_.push_back(float(-1));                                //-1 norm to vector for error
        break;
      }
      nr_k_.push_back(nr_k);                                //write residual norm to vector

      ApplyDFH(rhs_mat, res_mat);   //result in rhs_mat[ 0 , value[1..nc+1] ]

      //L2-/TV-/TGV-Solve
      Solve(&u_mat_, c_mat_, rhs_mat, &u0_mat_, tvits, alpha, beta, &du_mat, dc_mat);

      agile::addMatrix(u_mat_,du_mat,u_mat_);

      for(unsigned i=0; i<num_coils_;++i)
        agile::addMatrix(c_mat_[i],dc_mat[i],c_mat_[i]);

      //safe actual u_mat and c_mat to u_mat_old, c_mat_old
      agile::copy(u_mat_, u_mat_old);
      CopyMatrixZ(c_mat_, c_mat_old, num_coils_);

      //reduce parameter

      alpha = (this->alpha_min_ >= (alpha_q_*alpha)) ? alpha_min_ : (alpha_q_*alpha);
      beta = (this->beta_min_ >= (beta_q_*beta)) ? beta_min_ : (beta_q_*beta);
      float tvits_help = 2 * float(tvits);

      tvits = ((int(tvits_help)) <= this->tvmax_) ? int(tvits_help) : tvmax_;
    }

    delete[] res_mat;
    res_mat = 0;
    delete[] rhs_mat;
    rhs_mat = 0;
    delete[] dc_mat;
    dc_mat = 0;
    delete[] c_mat_old;
    c_mat_old = 0;
  }


  // ---------------------------------------------------------------
  //! \brief ApplyDFH()
  //!  calculates ApplyDFH
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::ApplyDFH(GPUMatrix<TType>* rhs_mat, const GPUMatrix<TType>* dx)
  {
    GPUMatrix<TType>* fdy;
    fdy = new GPUMatrix<TType> [num_coils_];
    for(unsigned i=0; i<num_coils_;++i)
      fdy[i].resize(num_rows_, num_columns_);

    GPUMatrix<TType> safe_mat(num_rows_, num_columns_,NULL);
    GPUMatrix<TType> tmp(num_rows_, num_columns_,NULL);

    agile::copy(zero_complex_mat_nxns_,tmp);   //set tmp to 0

    CenterdIFFTpattern(dx,fdy,num_coils_);

    for(unsigned i=0; i<num_coils_;++i)
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
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::ApplyM(GPUMatrix<TType>* gu, GPUMatrix<TType>* gc,
                                 const GPUMatrix<TType>* du, const GPUMatrix<TType>* dc)
  {
    GPUMatrix<TType> safe_mat(num_rows_, num_columns_,NULL);
    GPUMatrix<TType> fdy(num_rows_, num_columns_,NULL);

    agile::copy(zero_complex_mat_nxns_, *gu);    //set gu-matrix to zero

    for(unsigned i=0; i<num_coils_;++i)
    {
      //W(dc)
      ApplyW(&dc[i],&safe_mat,1);

      //u*W(dc)
      agile::multiplyElementwise(u_mat_ , safe_mat, safe_mat);
      //c*du
      agile::multiplyElementwise(cw_mat_[i] , *du, fdy);
      //u*W(dc) + c*du
      agile::addMatrix(safe_mat, fdy, fdy);

      CenterdFFTpattern(&fdy,&fdy,1);
      CenterdIFFTpattern(&fdy,&fdy,1);

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
  //! \brief calcLipschitz()
  //!  calculates Lipschitz constant - initial solve-calculation
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  TType2 IRGN<TType, TType2>::calcLipschitz()
  {

    std::complex<typename agile::to_real_type<TType>::type> l1;
    std::complex<typename agile::to_real_type<TType>::type> l2;
    std::complex<typename agile::to_real_type<TType>::type> l2_safe;
    l2.real(0);
    l2.imag(0);
    typename agile::to_real_type<TType>::type norm_;
    typename agile::to_real_type<TType>::type one_;
    typename agile::to_real_type<TType>::type two_;
    norm_=0;
    one_ =1;
    two_ = 2;
    typename agile::to_real_type<TType>::type L;

    //Matrix-declaration
    GPUMatrix<TType> x1_mat;
    GPUMatrix<TType> y1_mat;
    GPUMatrix<TType>* x2_mat = new GPUMatrix<TType> [IRGN<TType,TType2>::num_coils_];
    GPUMatrix<TType>* y2_mat = new GPUMatrix<TType> [IRGN<TType,TType2>::num_coils_];

    x1_mat.resize(num_rows_, num_columns_);
    y1_mat.resize(num_rows_, num_columns_);
    for(unsigned i=0; i<num_coils_;++i)
    {
      x2_mat[i].resize(num_rows_, num_columns_);
      y2_mat[i].resize(num_rows_, num_columns_);
    }
    //Matrix-definition with random matrix
    agile::copy(random1_mat_,x1_mat);
    for(unsigned i=0; i<IRGN<TType,TType2>::num_coils_;++i)
    {
      agile::copy(random2_mat_[i],x2_mat[i]);
    }

    //---start calculation:
    ApplyM(&y1_mat,y2_mat,&x1_mat,x2_mat);

    for (unsigned ii = 0; ii < 10; ++ii)
    {
      //x1=y1./norm(y1)
      norm_ = agile::norm2(y1_mat);
      if (norm_ == 0)
        agile::copy(y1_mat, x1_mat);
      else
        agile::scale(one_/norm_,y1_mat,x1_mat);

      //x2=y2./norm(y2)
      norm_ = 0;
      for(unsigned i=0; i<IRGN<TType,TType2>::num_coils_;++i)
      {
        norm_ = norm_ + std::pow(agile::norm2(y2_mat[i]),2);
      }
      norm_ = std::sqrt(norm_);

      if (norm_ == 0)
        IRGN<TType,TType2>::CopyMatrixZ(y2_mat, x2_mat, IRGN<TType,TType2>::num_coils_);
      else
      {
        for(unsigned i=0; i<IRGN<TType,TType2>::num_coils_;++i)
          agile::scale(one_/norm_,y2_mat[i],x2_mat[i]);
      }
      ApplyM(&y1_mat,y2_mat,&x1_mat,x2_mat);
    }

    agile::dotProduct(x1_mat,y1_mat,l1);

    for(unsigned i=0; i<IRGN<TType,TType2>::num_coils_;++i)
    {
      agile::dotProduct(x2_mat[i],y2_mat[i],l2_safe);
      l2 = l2 + l2_safe;
    }

    L = 2 * ((std::abs(l1) >= std::abs(l2)) ? std::abs(l1) : std::abs(l2));

    delete[] x2_mat;
    x2_mat = 0;
    delete[] y2_mat;
    y2_mat = 0;

    return L;
  }


  // ---------------------------------------------------------------
  //! \brief Postprocess()
  //!  calculates postprocess
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::Postprocess()
  {
    GPUMatrix<typename agile::to_real_type<TType>::type> cscale_mat(num_rows_, num_columns_, NULL);
    GPUMatrix<typename agile::to_real_type<TType>::type> absval(num_rows_, num_columns_, NULL);

    GPUMatrix<TType> safe_mat(num_rows_, num_columns_, NULL);

    agile::copy(this->zero_mat_nxns_,cscale_mat);   //init with 0

    for(unsigned i=0; i<num_coils_;++i)
    {
      this->ApplyW(&c_mat_[i],&safe_mat,1);
      agile::absMatrix(safe_mat,absval);
      agile::pow(typename TType::value_type(2),absval,absval);
      agile::addMatrix(absval,cscale_mat,cscale_mat);
    }
    agile::sqrt(cscale_mat,cscale_mat);

    agile::multiplyElementwise(u_mat_, cscale_mat, us_mat_);

    agile::scale(dscale_,us_mat_,us_mat_);

    //Postprocess to Image-Data:
    agile::absMatrix(us_mat_,image_mat_);

    //test
    //test_max(us_mat_,"usmat");
  }


  // ---------------------------------------------------------------
  //! \brief CopyMatrixZ()
  //!  Copy number "num_z" of Matrix "in_mat" to "out_mat"
  // ---------------------------------------------------------------
  template <typename TType, typename TType2>
  void IRGN<TType, TType2>::CopyMatrixZ(const GPUMatrix<TType>* in_mat,
                   GPUMatrix<TType>* out_mat,
                   unsigned int num_z)
  {
    for(unsigned i=0; i<num_z;++i)
    {
      agile::copy(in_mat[i],out_mat[i]);
    }
  }

  template <typename TType, typename TType2>
  bool IRGN<TType, TType2>::irgn_param_test(IRGN_Params &param)
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

} //namespace agile
