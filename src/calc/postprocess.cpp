#include "agile/calc/postprocess.hpp"

namespace agile
{
  // ---------------------------------------------------------------
  //! \brief init()
  //!  initialice values
  // ---------------------------------------------------------------
  template <typename TType, typename TType_real>
  void PostProcess<TType, TType_real>::init()
  {
    for(int i=0; i<num_columns_*num_rows_;  ++i)
      zero_vec_cpu_.push_back(TType_real(0));

    zero_vec_gpu_.assignFromHost(num_rows_,num_columns_,&zero_vec_cpu_[0]);
  }


  // ---------------------------------------------------------------
  //! \brief calc_abs(const GPUMatrix<FloatType> &in_mat, GPUMatrix<FloatType> &out_mat)
  //!  calculates the absvalues for given matrix
  // ---------------------------------------------------------------
  template <typename TType, typename TType_real>
  void PostProcess<TType, TType_real>::calc_abs(const GPUMatrix<TType>* in_mat, GPUMatrix<TType_real>& out_mat)
  {
    agile::GPUMatrix<TType_real> abs_val;
    abs_val.resize(num_rows_,num_columns_);
    out_mat.resize(num_rows_,num_columns_);
    agile::copy(zero_vec_gpu_, abs_val);
    agile::copy(zero_vec_gpu_, out_mat);

    for(int i=0; i<num_coils_; ++i)
    {
      agile::absMatrix(in_mat[i],abs_val);                //calculate abs-val for each Coil Matrix
      agile::pow(typename TType::value_type(2),abs_val,abs_val);                 //square all values
      agile::addMatrix(abs_val,out_mat,out_mat);            //add matrix*/
    }
    agile::sqrt(out_mat,out_mat);                   //calulate square root for final solution
  }



  // ---------------------------------------------------------------
  //! \brief calc_phase(const GPUMatrix<FloatType> &in_mat, GPUMatrix<FloatType> &out_mat)
  //!  calculates the phase-values for given matrix
  // ---------------------------------------------------------------
  template <typename TType, typename TType_real>
  void PostProcess<TType, TType_real>::calc_phase(const GPUMatrix<TType>* in_mat, GPUMatrix<TType_real> &out_mat)
  {
    agile::GPUMatrix<TType_real> phase_val;
    phase_val.resize(num_rows_,num_columns_);
    out_mat.resize(num_rows_,num_columns_);
    agile::copy(zero_vec_gpu_, phase_val);
    agile::copy(zero_vec_gpu_, out_mat);

    for(int i=0; i<num_coils_; ++i)
    {
      agile::phase(in_mat[i],phase_val);                //calculate abs-val for each Coil Matrix
      agile::addMatrix(phase_val,out_mat,out_mat);            //add matrix*/
    }
  }

  // ---------------------------------------------------------------
  //! \brief calc_real(const GPUMatrix<FloatType> &in_mat, GPUMatrix<FloatType> &out_mat)
  //!  calculates the real-values for given matrix
  // ---------------------------------------------------------------
  template <typename TType, typename TType_real>
  void PostProcess<TType, TType_real>::calc_real(const GPUMatrix<TType>* in_mat, GPUMatrix<TType_real> &out_mat)
  {
    agile::GPUMatrix<typename to_real_type<TType>::type> real_val;
    real_val.resize(num_rows_,num_columns_);
    out_mat.resize(num_rows_,num_columns_);
    agile::copy(zero_vec_gpu_, real_val);
    agile::copy(zero_vec_gpu_, out_mat);

    for(int i=0; i<num_coils_;++i)
    {
      agile::real(in_mat[i],real_val);
      agile::addMatrix(real_val,out_mat,out_mat);
    }

  }

  // ---------------------------------------------------------------
  //! \brief calc_real(const GPUMatrix<FloatType> &in_mat, GPUMatrix<FloatType> &out_mat)
  //!  calculates the imag-values for given matrix
  // ---------------------------------------------------------------
  template <typename TType, typename TType_real>
  void PostProcess<TType, TType_real>::calc_imag(const GPUMatrix<TType>* in_mat, GPUMatrix<TType_real> &out_mat)
  {
    agile::GPUMatrix<typename to_real_type<TType>::type> imag_val;
    imag_val.resize(num_rows_,num_columns_);
    out_mat.resize(num_rows_,num_columns_);
    agile::copy(zero_vec_gpu_, imag_val);
    agile::copy(zero_vec_gpu_, out_mat);

    for(int i=0; i<num_coils_;++i)
    {
      agile::imag(in_mat[i],imag_val);
      agile::addMatrix(imag_val,out_mat,out_mat);
    }
  }

}// namespace agile
