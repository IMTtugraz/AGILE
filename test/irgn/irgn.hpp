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


#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_matrix.hpp"

#include "agile/exception.hpp"

#include "agile/io/file.hpp"

#include <iostream>
#include <iomanip>
#include <cufft.h>

#include "l2solve.hpp"
//#include "tvsolve.hpp"
//#include "tgvsolve.hpp"

#include "matrixhelper.h"           //my include for simple vector, matrix console ouput

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


template <typename TType>
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
      close_matrixlog(myfile);

      cufftDestroy(fftplan_);

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
        AGILE_ASSERT(false, StandardException::ExceptionMessage(
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


    void set_num_coils(unsigned int num_coils)
    {
      AGILE_ASSERT(num_coils <= 0,
                    StandardException::ExceptionMessage(
                 "num_coils == 0 - error"));

      num_coils_=num_coils;
    }
    void set_num_rows(unsigned int num_rows)
    {
      AGILE_ASSERT(num_rows <= 0,
                    StandardException::ExceptionMessage(
                 "num_rows == 0 - error"));

      num_rows_=num_rows;
    }
    void set_num_columns(unsigned int num_columns)
    {
      AGILE_ASSERT(num_columns <= 0,
                    StandardException::ExceptionMessage(
                 "num_columns == 0 - error"));

      num_columns_=num_columns;
    }

    void set_coil(GPUMatrix<TType>* coil, unsigned int num_coils)
    {
      AGILE_ASSERT(coil == NULL,
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
      return *us_mat_;
    }

    GPUMatrix<typename agile::to_real_type<TType>::type>* get_image()
    {
      return &image_mat_;
    }

    void HighFreqPenalty();
    void Normalize();
    void Iteration();
    void Postprocess();

    void ApplyM(GPUMatrix<TType>* gu, GPUMatrix<TType>* gc,
                const GPUMatrix<TType>* du, const GPUMatrix<TType>* dc);



    void CenterdFFT(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);
    void CenterdIFFT(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);
    void CenterdFFTpattern(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);
    void CenterdIFFTpattern(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);


  protected:

    std::fstream myfile;


    unsigned int num_coils_;     // number of coils
    unsigned int num_rows_;      // number of rows
    unsigned int num_columns_;   // number of columns

//    void CenterdFFT(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);
//    void CenterdIFFT(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);
//    void CenterdFFTpattern(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);
//    void CenterdIFFTpattern(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);


    void ApplyW(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);
    void ApplyWH(const GPUMatrix<TType>* in_mat, GPUMatrix<TType>* out_mat, unsigned int num_z);

    void ApplyDFH(GPUMatrix<TType>* rhs_mat, const GPUMatrix<TType>* dx);

//    void ApplyM(GPUMatrix<TType>* gu, GPUMatrix<TType>* gc,
//                const GPUMatrix<TType>* du, const GPUMatrix<TType>* dc);

    void CopyMatrixZ(const GPUMatrix<TType>* in_mat,
                     GPUMatrix<TType>* out_mat,
                     unsigned int num_z);

    //Pure Virtual Solve-Methode
    virtual void Solve(const GPUMatrix<TType>* u, const GPUMatrix<TType>* c,
                                   const GPUMatrix<TType>* rhs, const GPUMatrix<TType>* u0,
                                   unsigned maxits, float alpha, float beta,
                                   GPUMatrix<TType>* du, GPUMatrix<TType>* dc) = 0;

    //nx nx  - matrix
    GPUMatrix<typename agile::to_real_type<TType>::type> zero_mat_nxnx_;
    GPUMatrix<typename agile::to_real_type<TType>::type> ones_mat_nxnx_;
    GPUMatrix<TType> zero_complex_mat_nxnx_;
    GPUMatrix<TType> ones_complex_mat_nxnx_;

    //nx ns  - matrix
    GPUMatrix<TType> zero_complex_mat_nxns_;
    GPUMatrix<TType> ones_complex_mat_nxns_;

    GPUMatrix<TType> random1_mat_;
    GPUMatrix<TType>* random2_mat_;

  private:

    void Init();
    float randomcalc(int i);

    bool irgn_param_test(IRGN_Params &param);

    unsigned int maxit_;         // maximum number of IRGN iterations
    unsigned char tvtype_;       // regularization term: 0: L2, 1: TV, 2: TGV
    unsigned int tvits_;         // initial number of gradient steps
    unsigned int tvmax_;         // upper bound on number of gradient steps
    float alpha_min_;            // final value of alpha
    float beta_min_;             // final value of beta: 0: no T(G)V effect, >0 effect
    float alpha0_;               // initial penalty alpha_0 (L2, sensitivites)
    float beta0_;                // initial penalty beta_0 (image)
    float alpha_q_;              // reduction factor for alpha
    float beta_q_;               // reduction factor for beta

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
    GPUMatrix<typename agile::to_real_type<TType>::type> pattern_;
    GPUMatrix<TType> pattern_complex_;

    std::vector<float> nr_k_;    //residual norm;

    //vector size nx nx  (with nx = num_rows_)
    std::vector<typename agile::to_real_type<TType>::type> zero_vec_cpu_nxnx_;
    std::vector<TType> zero_complex_vec_cpu_nxnx_;
    std::vector<typename agile::to_real_type<TType>::type> ones_vec_cpu_nxnx_;
    std::vector<TType> ones_complex_vec_cpu_nxnx_;

    //vector size nx ns  (with nx = num_rows_, ns = num_columns_)
    std::vector<typename agile::to_real_type<TType>::type> zero_vec_cpu_nxns_;
    std::vector<TType> zero_complex_vec_cpu_nxns_;
    std::vector<typename agile::to_real_type<TType>::type> ones_vec_cpu_nxns_;
    std::vector<TType> ones_complex_vec_cpu_nxns_;

    //random - vector
    std::vector<TType> x1_;
    std::vector<TType> x2_;
};

#endif // AGILE_IRGN_HPP
