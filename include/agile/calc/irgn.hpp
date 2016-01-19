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




//==================================================================================
//
//                    IRGN.HPP
//
//==================================================================================

#ifndef AGILE_IRGN_HPP
#define AGILE_IRGN_HPP

#include <algorithm>
#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_matrix.hpp"
#include "agile/gpu_timer.hpp"

#include "agile/exception.hpp"

#include "agile/io/file.hpp"
#include "agile/calc/fft.hpp"

#include <iostream>
#include <iomanip>
#include <cufft.h>

//#include "matrixhelper.h"           //my include for simple vector, matrix console ouput

namespace agile
{

//_____________ test ___________

template <typename TType>
void test_max(const agile::GPUMatrix<TType>& in_mat, std::string text)
{

  std::vector<float> mat_host;
  agile::GPUMatrix<float> mat_abs(in_mat.getNumRows(),in_mat.getNumColumns(),NULL);
  agile::absMatrix(in_mat, mat_abs);
  mat_abs.copyToHost(mat_host);
  float mat_host_max = *std::max_element(mat_host.begin(),mat_host.end());
  std::cout<<"\n "<<text<<": "<<mat_host_max;

  /*
  std::vector<float> u_mat_host;
  std::vector<float> c_mat_host;
  std::vector<float> us_mat_host;
  agile::GPUMatrix<float> u_mat_abs(u_mat_.getNumRows(),u_mat_.getNumColumns(),NULL);
  agile::GPUMatrix<float> c_mat_abs(c_mat_[0].getNumRows(),c_mat_[0].getNumColumns(),NULL);;
  agile::GPUMatrix<float> us_mat_abs(us_mat_.getNumRows(),us_mat_.getNumColumns(),NULL);;
  agile::absMatrix(u_mat_, u_mat_abs);
  agile::absMatrix(c_mat_[0], c_mat_abs);
  agile::absMatrix(us_mat_, us_mat_abs);

  u_mat_abs.copyToHost(u_mat_host);
  c_mat_abs.copyToHost(c_mat_host);
  us_mat_abs.copyToHost(us_mat_host);
  float u_mat_host_max = *std::max_element(u_mat_host.begin(),u_mat_host.end());
  float c_mat_host_max = *std::max_element(c_mat_host.begin(),c_mat_host.end());
  float us_mat_host_max = *std::max_element(us_mat_host.begin(),us_mat_host.end());

  std::cout<<"\n u_mat_host_max: "<<u_mat_host_max;
  std::cout<<"\n c_mat_host_max: "<<c_mat_host_max;
  std::cout<<"\n us_mat_host_max: "<<us_mat_host_max;
  */
}

//_____________ testend ___________



  using agile::GPUMatrix;
  using agile::GPUVector;


  struct IRGN_Params
   {
      unsigned int maxit;         // maximum number of IRGN iterations
      unsigned char tvtype;       // regularization term: 0: L2, 1: TV, 2: TGV
      unsigned int tvits;         // initial number of gradient steps
      unsigned int tvmax;         // upper bound on number of gradient steps
      float alpha_min;            // final value of alpha
      float beta_min;             // final value of beta: 0: no T(G)V effect, >0 effect
      float alpha0;               // initial penalty alpha_0 (L2, sensitivites)
      float beta0;                // initial penalty beta_0 (image)
      float alpha_q;              // reduction factor for alpha
      float beta_q;               // reduction factor for beta
    };


  template <typename TType, typename TType2>
  class IRGN
  {
    public:

      typedef agile::to_real_type<TType> TType_real;

      //! \brief Default constructor.
      //!
      //! The default constructor creates an empty IRGN Class.
      IRGN()
        : maxit_(0), tvtype_(0), tvits_(0), tvmax_(0),alpha_min_(0),
          beta_min_(0),alpha0_(0),beta0_(0),alpha_q_(0),beta_q_(0)
      {
        num_coils_= 0;
        num_rows_= 0;
        num_columns_ = 0;
      }


      //! \brief Constructor.
      //!
      //! The Constructor creates an IRGN.
      IRGN(GPUMatrix<TType>* coil, unsigned int num_coils, IRGN_Params param)
      {
        coil_ = coil;
        set_num_coils(num_coils);
        set_num_rows(coil->getNumRows());
        set_num_columns(coil->getNumColumns());

        set_param(param);

        this->Init();
      }

      //! \brief Constructor.
      //!
      //! The Constructor creates an initialized IRGN.
      IRGN(GPUMatrix<TType>* coil, unsigned int num_coils)
        : maxit_(6), tvtype_(1), tvits_(20), tvmax_(1000),alpha_min_(0),
          beta_min_(0),alpha0_(1),beta0_(1),alpha_q_(0.1),beta_q_(0.2)
      {
        coil_ = coil;
        set_num_coils(num_coils);
        set_num_rows(coil->getNumRows());
        set_num_columns(coil->getNumColumns());

        this->Init();
      }

      //! \brief Destructor.
      virtual ~IRGN()
      {
        //close_matrixlog(myfile);

        delete fftobj_;
        fftobj_ = 0;

        delete[] c_mat_;
        c_mat_ = 0;
        delete[] cw_mat_;
        cw_mat_ = 0;
        delete[] cw_bar_mat_;
        cw_bar_mat_ = 0;
        delete[] random2_mat_;
        random2_mat_ = 0;
      }

      //! sets the IRGN parameters.
      void set_param(IRGN_Params param)
      {
        if( !irgn_param_test(param) )
          AGILE_ASSERT(true, StandardException::ExceptionMessage(
                     "Check your IRGN-Parameter"));

        maxit_ = param.maxit;
        tvtype_ = param.tvtype;
        tvits_ = param.tvits;
        tvmax_ = param.tvmax;
        alpha_min_ = param.alpha_min;
        beta_min_ = param.beta_min;
        alpha0_ = param.alpha0;
        beta0_ = param.beta0;
        alpha_q_ = param.alpha_q;
        beta_q_ = param.beta_q;
      }

      void set_break_calc(bool break_calc)
      {
        _break_calc = break_calc;
      }


      void set_num_coils(unsigned int num_coils)
      {

        AGILE_ASSERT(num_coils > 0,
                      StandardException::ExceptionMessage(
                   "num_coils == 0 - error"));

        num_coils_=num_coils;
      }
      void set_num_rows(unsigned int num_rows)
      {
        AGILE_ASSERT(num_rows > 0,
                      StandardException::ExceptionMessage(
                   "num_rows == 0 - error"));

        num_rows_=num_rows;
      }
      void set_num_columns(unsigned int num_columns)
      {
        AGILE_ASSERT(num_columns > 0,
                      StandardException::ExceptionMessage(
                   "num_columns == 0 - error"));

        num_columns_=num_columns;
      }

      void set_coil(GPUMatrix<TType>* coil, unsigned int num_coils)
      {
        AGILE_ASSERT(coil != NULL,
                      StandardException::ExceptionMessage(
                   "coil == NULL - error"));
        coil_= coil;
        set_num_coils(num_coils);
      }

      GPUMatrix<TType>* get_coil()
      {
        return coil_;
      }

      unsigned int get_numcoils()
      {
        return num_coils_;
      }

       std::vector<float> get_nr_k()
      {
        return nr_k_;
      }

      GPUMatrix<TType>* get_us_mat()
      {
        return &us_mat_;
      }

      GPUMatrix<typename agile::to_real_type<TType>::type>* get_image()
      {
        return &image_mat_;
      }

      void HighFreqPenalty();
      void Normalize();
      void Iteration();
      void Postprocess();

    protected:

      void CenterdFFT(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);
      void CenterdIFFT(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);
      void CenterdFFTpattern(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);
      void CenterdIFFTpattern(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);


      std::fstream myfile;
      agile::GPUTimer Timer;
      double timer_value;

      unsigned int num_coils_;     // number of coils
      unsigned int num_rows_;      // number of rows
      unsigned int num_columns_;   // number of columns

      void ApplyW(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);
      void ApplyWH(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);

      void ApplyDFH(GPUMatrix<TType>* rhs_mat, const GPUMatrix<TType>* dx);

      void ApplyM(GPUMatrix<TType>* gu, GPUMatrix<TType>* gc,
                  const GPUMatrix<TType>* du, const GPUMatrix<TType>* dc);

      void CopyMatrixZ(const GPUMatrix<TType>* in_mat,
                       GPUMatrix<TType>* out_mat,
                       unsigned int num_z);

      TType2 calcLipschitz();

      //Pure Virtual Solve-Methode
      virtual void Solve(const GPUMatrix<TType>* u, const GPUMatrix<TType>* c,
                                     const GPUMatrix<TType>* rhs, const GPUMatrix<TType>* u0,
                                     unsigned maxits, TType2 alpha, TType2 beta,
                                     GPUMatrix<TType>* du, GPUMatrix<TType>* dc) = 0;

      //nx ns  - matrix
      GPUMatrix<typename agile::to_real_type<TType>::type> zero_mat_nxns_;
      GPUMatrix<typename agile::to_real_type<TType>::type> ones_mat_nxns_;
      GPUMatrix<TType> zero_complex_mat_nxns_;
      GPUMatrix<TType> ones_complex_mat_nxns_;

      GPUMatrix<TType> random1_mat_;
      GPUMatrix<TType>* random2_mat_;

    private:

      void Init();
      typename TType::value_type randomcalc(int i);

      bool irgn_param_test(IRGN_Params &param);

      unsigned int maxit_;         // maximum number of IRGN iterations
      unsigned char tvtype_;       // regularization term: 0: L2, 1: TV, 2: TGV
      //unsigned int tvits_;         // initial number of gradient steps
      unsigned int tvmax_;         // upper bound on number of gradient steps
      float alpha_min_;            // final value of alpha
      float beta_min_;             // final value of beta: 0: no T(G)V effect, >0 effect
      float alpha0_;               // initial penalty alpha_0 (L2, sensitivites)
      float beta0_;                // initial penalty beta_0 (image)
      float alpha_q_;              // reduction factor for alpha
      float beta_q_;               // reduction factor for beta
      unsigned int tvits_;         // initial number of gradient steps

      typename agile::to_real_type<TType>::type dscale_;

      GPUMatrix<TType>* coil_;

      cufftHandle fftplan_;

      GPUMatrix<typename agile::to_real_type<TType>::type> w_mat_;

      GPUMatrix<TType> u_mat_;
      GPUMatrix<TType>* c_mat_;
      GPUMatrix<TType> u0_mat_;
      GPUMatrix<TType> u_bar_mat_;
      GPUMatrix<TType>* cw_mat_;
      GPUMatrix<TType>* cw_bar_mat_;
      GPUMatrix<TType> us_mat_;
      GPUMatrix<typename agile::to_real_type<TType>::type> image_mat_;
      GPUMatrix<typename agile::to_real_type<TType>::type> image_mat_final_;

      std::vector<float> nr_k_;    //residual norm;

      //vector size nx ns  (with nx = num_rows_, ns = num_columns_)
      std::vector<typename agile::to_real_type<TType>::type> zero_vec_cpu_nxns_;
      std::vector<TType> zero_complex_vec_cpu_nxns_;
      std::vector<typename agile::to_real_type<TType>::type> ones_vec_cpu_nxns_;
      std::vector<TType> ones_complex_vec_cpu_nxns_;

      //random - vector
      std::vector<TType> x1_;
      std::vector<TType> x2_;

      bool _break_calc;

      agile::FFT<TType>* fftobj_;

/*      unsigned int num_delta_yo_;       // number of added zeros rows up
      unsigned int num_delta_yu_;       // number of added zeros rows bottom
      unsigned int num_delta_xo_;       // number of added zeros rows up
      unsigned int num_delta_xu_;       // number of added zeros rows bottom
*/
  };
}
#endif // AGILE_IRGN_HPP



