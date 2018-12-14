#include "agile/calc/fft.hpp"

namespace agile
{
  // ---------------------------------------------------------------
  //! \brief setpattern(GPUMatrix<TType>* in_mat)
  //!  creates the patternmatrix from given in_mat
  // ---------------------------------------------------------------
  template <typename TType>
  void FFT<TType>::calc_pattern(const GPUMatrix<TType>& in_mat)
  {

    unsigned num_rows = in_mat.getNumRows();
    unsigned num_cols = in_mat.getNumColumns();

    pattern_.resize(num_rows, num_cols);
    pattern_complex_.resize(num_rows, num_cols);

    //generate pattern
    agile::pattern(in_mat,pattern_);

    agile::multiplyElementwise(pattern_, ones_complex_mat_nxns_, pattern_complex_);
  }

//=========================================================
//============ Template - Helper - setfftplan(unsigned num_rows, unsigned num_columns) =============
//=========================================================
  template <typename TType>
  struct setfftplan_Helper;

  template <>
  struct setfftplan_Helper<std::complex<float> >
  {
      static void setfftpl(unsigned num_rows, unsigned num_columns, cufftHandle* fftplan)
      {
          cufftPlan2d(fftplan, num_rows, num_columns, CUFFT_C2C);
      }
  };

  template <>
  struct setfftplan_Helper<std::complex<double> >
  {
      static void setfftpl(unsigned num_rows, unsigned num_columns, cufftHandle* fftplan)
      {
          cufftPlan2d(fftplan, num_rows, num_columns, CUFFT_Z2Z);
      }
  };

  //! \brief Multiply a GPU matrix with a scalar (host function).
  template <typename TType>
  void FFT<TType>::setfftplan(unsigned num_rows, unsigned num_columns)
  {

      setfftplan_Helper<TType>::setfftpl(num_rows, num_columns, &fftplan_);
      num_rows_ = num_rows;
      num_columns_ = num_columns;

      Init();
  }

//=========================================================
//============ Template - Helper - CenterdIFFT(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat)  =============
//=========================================================
  template <typename TType>
  struct CenterdIFFT_Helper;

  template <>
  struct CenterdIFFT_Helper<std::complex<float> >
  {
    static cufftResult_t centerdifft(const std::complex<float>* in_mat, std::complex<float>* out_mat,
                                     cufftHandle* fftplan)
    {
      cufftResult_t cufftResult;

      cufftResult = cufftExecC2C(*fftplan,
                   (cufftComplex*)in_mat,
                   (cufftComplex*)out_mat,
                   CUFFT_INVERSE);
      
      AGILE_ASSERT(result == CUFFT_SUCCESS,
                        StandardException::ExceptionMessage(
                          "Error during FFT procedure"));

      return cufftResult;
    }
  };

  template <>
  struct CenterdIFFT_Helper<std::complex<double> >
  {
      static cufftResult_t centerdifft(const std::complex<double>* in_mat, std::complex<double>* out_mat,
                                       cufftHandle* fftplan)
      {
        cufftResult_t cufftResult;

        cufftResult = cufftExecZ2Z(*fftplan,
                     (cufftDoubleComplex*)in_mat,
                     (cufftDoubleComplex*)out_mat,
                     CUFFT_INVERSE);
      
        AGILE_ASSERT(result == CUFFT_SUCCESS,
                        StandardException::ExceptionMessage(
                          "Error during FFT procedure"));

        return cufftResult;
      }
  };

  // ---------------------------------------------------------------
  //! \brief CenterdIFFT(GPUMatrix<TType>* ifft_mat, GPUMatrix<TType>* out_mat)
  //!  calculates the Centerd IFFT for given matrix
  // ---------------------------------------------------------------
  template <typename TType>
  int FFT<TType>::CenterdIFFT(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat)
  {
    cufftResult_t cufftResult;
    agile::copy(in_mat,out_mat);
    agile::fftshift(out_mat);

    cufftResult = CenterdIFFT_Helper<TType>::centerdifft(out_mat.data(), out_mat.data(), &fftplan_);
    
    real_type val_sqrt = (real_type)1.0/std::sqrt(in_mat.getNumRows() * in_mat.getNumColumns());
    agile::scale(val_sqrt, out_mat, out_mat);

    agile::ifftshift(out_mat);

    return int(cufftResult);
  }
  
  // ---------------------------------------------------------------
  //! \brief CenteredInverse(GPUVector<TType>* in_vec, GPUVector<TType>* out_vec, unsigned in_offset, unsigned out_offset)
  //!  calculates the Centerd IFFT for a given vector
  //!  The vector is assumed to be a flattened matrix of the defined 
  //!  dimensions 
  // ---------------------------------------------------------------
  template <typename TType>
  int FFT<TType>::CenteredInverse(const GPUVector<TType>& in_vec, GPUVector<TType>& out_vec, unsigned in_offset, unsigned out_offset)
  {
    if (out_offset == -1)
      out_offset = in_offset;

    cufftResult_t cufftResult;
    TType* data = out_vec.data()+out_offset;

    //do not copy the whole vector
    //-> offset
    //agile::copy(in_vec,out_vec);
    agile::lowlevel::get_content(in_vec.data(), 1, num_rows_ * num_columns_, 0, 
        in_offset, data, 1, num_rows_*num_columns_);

    agile::lowlevel::fftshift(data, num_rows_, num_columns_);

    cufftResult = CenterdIFFT_Helper<TType>::centerdifft(data, data, &fftplan_);
    
    real_type val_sqrt = (real_type)1.0/std::sqrt(num_rows_ * num_columns_);
    //do not scale the whole vector
    //-> offset
    // agile::scale(val_sqrt, out_vec, out_vec);
    lowlevel::scale(val_sqrt, data, data, num_rows_ * num_columns_);

    agile::lowlevel::ifftshift(data, num_rows_, num_columns_);

    return int(cufftResult);
  }
  
  // ---------------------------------------------------------------
  //! \brief Inverse(GPUVector<TType>* in_vec, GPUVector<TType>* out_vec, unsigned in_offset, unsigned out_offset)
  //!  calculates the IFFT for a given vector
  //!  The vector is assumed to be a flattened matrix of the defined 
  //!  dimensions 
  // ---------------------------------------------------------------
  template <typename TType>
  int FFT<TType>::Inverse(const GPUVector<TType>& in_vec, GPUVector<TType>& out_vec, unsigned in_offset, unsigned out_offset)
  {
    if (out_offset == -1)
      out_offset = in_offset;

    cufftResult_t cufftResult;
    const TType* in_data = in_vec.data()+in_offset;
    TType* out_data = out_vec.data()+out_offset;

    cufftResult = CenterdIFFT_Helper<TType>::centerdifft(in_data, out_data, &fftplan_);
    
    real_type val_sqrt = (real_type)1.0/std::sqrt(num_rows_ * num_columns_);
    
    lowlevel::scale(val_sqrt, out_data, out_data, num_rows_ * num_columns_);

    return int(cufftResult);
  }

//=========================================================
//============ Template - Helper - CenterdFFT(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat)  =============
//=========================================================
  template <typename TType>
  struct CenterdFFT_Helper;

  template <>
  struct CenterdFFT_Helper<std::complex<float> >
  {
    static cufftResult_t centerdfft(const std::complex<float>* in_vec, std::complex<float>* out_vec,
                                     cufftHandle* fftplan)
    {
      cufftResult_t cufftResult;

      cufftResult = cufftExecC2C(*fftplan,
                   (cufftComplex*)in_vec,
                   (cufftComplex*)out_vec,
                   CUFFT_FORWARD);
      
      AGILE_ASSERT(result == CUFFT_SUCCESS,
                        StandardException::ExceptionMessage(
                          "Error during FFT procedure"));

      return cufftResult;
    };
  };

  template <>
  struct CenterdFFT_Helper<std::complex<double> >
  {
      static cufftResult_t centerdfft(const std::complex<double>* in_vec, std::complex<double>* out_vec,
                                       cufftHandle* fftplan)
      {
        cufftResult_t cufftResult;

        cufftResult = cufftExecZ2Z(*fftplan,
                     (cufftDoubleComplex*)in_vec,
                     (cufftDoubleComplex*)out_vec,
                     CUFFT_FORWARD);

        AGILE_ASSERT(result == CUFFT_SUCCESS,
                        StandardException::ExceptionMessage(
                          "Error during FFT procedure"));

        return cufftResult;
      };
  };

  // ---------------------------------------------------------------
  //! \brief CenterdFFT(GPUMatrix<TType>* ifft_mat, GPUMatrix<TType>* out_mat)
  //!  calculates the Centerd IFFT for given matrix
  // ---------------------------------------------------------------
  template <typename TType>
  int FFT<TType>::CenterdFFT(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat)
  {
    cufftResult_t cufftResult;
    agile::copy(in_mat,out_mat);
    agile::ifftshift(out_mat);

    cufftResult = CenterdFFT_Helper<TType>::centerdfft(out_mat.data(), out_mat.data(), &fftplan_);
    
    real_type val_sqrt = (real_type)1.0/std::sqrt(in_mat.getNumRows() * in_mat.getNumColumns());
    agile::scale(val_sqrt, out_mat, out_mat);

    agile::fftshift(out_mat);

    return int(cufftResult);
  }

  // ---------------------------------------------------------------
  //! \brief CenteredForward(GPUVector<TType>* in_vec, GPUVector<TType>* out_vec, unsigned in_offset, unsigned out_offset)
  //!  calculates the centered FFT for a given vector
  //!  The vector is assumed to be a flattened matrix of the defined 
  //!  dimensions 
  // ---------------------------------------------------------------
  template <typename TType>
  int FFT<TType>::CenteredForward(const GPUVector<TType>& in_vec, GPUVector<TType>& out_vec, unsigned in_offset, unsigned out_offset)
  {
    if (out_offset == -1)
      out_offset = in_offset;

    cufftResult_t cufftResult;
    TType* data = out_vec.data()+out_offset;

    //do not copy the whole vector
    //-> offset
    //agile::copy(in_vec,out_vec);
    agile::lowlevel::get_content(in_vec.data(), 1, num_rows_ * num_columns_, 0, 
        in_offset, data, 1, num_rows_*num_columns_);

    agile::lowlevel::ifftshift(data, num_rows_, num_columns_);

    cufftResult = CenterdFFT_Helper<TType>::centerdfft(data, data, &fftplan_);

    real_type val_sqrt = (real_type)1.0/std::sqrt(num_rows_ * num_columns_);

    //do not scale the whole vector
    //-> offset
    //agile::scale(val_sqrt, out_vec, out_vec);
    lowlevel::scale(val_sqrt, data, data, num_rows_ * num_columns_);

    agile::lowlevel::fftshift(data, num_rows_, num_columns_);

    return int(cufftResult);
  }
  
  // ---------------------------------------------------------------
  //! \brief Forward(GPUVector<TType>* in_vec, GPUVector<TType>* out_vec, unsigned in_offset, unsigned out_offset)
  //!  calculates the FFT for a given vector
  //!  The vector is assumed to be a flattened matrix of the defined 
  //!  dimensions 
  // ---------------------------------------------------------------
  template <typename TType>
  int FFT<TType>::Forward(const GPUVector<TType>& in_vec, GPUVector<TType>& out_vec, unsigned in_offset, unsigned out_offset)
  {
    if (out_offset == -1)
      out_offset = in_offset;

    cufftResult_t cufftResult;
    const TType* in_data = in_vec.data()+in_offset;
    TType* out_data = out_vec.data()+out_offset;

    cufftResult = CenterdFFT_Helper<TType>::centerdfft(in_data, out_data, &fftplan_);

    real_type val_sqrt = (real_type)1.0/std::sqrt(num_rows_ * num_columns_);

    lowlevel::scale(val_sqrt, out_data, out_data, num_rows_ * num_columns_);

    return int(cufftResult);
  }
//=========================================================
//============ Template CenterdIFFTpattern(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat)  =============
//=========================================================

// ---------------------------------------------------------------
//! \brief CenterdIFFTpattern(GPUMatrix<TType>* ifft_mat, GPUMatrix<TType>* out_mat)
//!  calculates the Centerd IFFT for given matrix with pattern
// ---------------------------------------------------------------
template <typename TType>
int FFT<TType>::CenterdIFFTpattern(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat)
{
  if((in_mat.getNumColumns()*in_mat.getNumRows()) != (pattern_complex_.getNumColumns()*pattern_complex_.getNumRows()))
  {
    std::cerr<<"\n in_mat size not equal pattern size";
    return -1;
  }

  cufftResult_t cufftResult;
  agile::copy(in_mat,out_mat);
  agile::multiplyElementwise(pattern_complex_,out_mat,out_mat);
  agile::fftshift(out_mat);

  cufftResult = CenterdIFFT_Helper<TType>::centerdifft(out_mat.data(), out_mat.data(), &fftplan_);
  
  real_type val_sqrt = (real_type)1.0/std::sqrt(in_mat.getNumRows() * in_mat.getNumColumns());
  agile::scale(val_sqrt, out_mat, out_mat);

  agile::ifftshift(out_mat);

  return int(cufftResult);
}

//=========================================================
//============ Template CenterdFFTpattern(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat)  =============
//=========================================================

// ---------------------------------------------------------------
//! \brief CenterdFFTpattern(GPUMatrix<TType>* ifft_mat, GPUMatrix<TType>* out_mat)
//!  calculates the Centerd FFT for given matrix with pattern
// ---------------------------------------------------------------
template <typename TType>
int FFT<TType>::CenterdFFTpattern(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat)
{
  if((in_mat.getNumColumns()*in_mat.getNumRows()) != (pattern_complex_.getNumColumns()*pattern_complex_.getNumRows()))
  {
    std::cerr<<"\n in_mat size not equal pattern size";
    return -1;
  }

  cufftResult_t cufftResult;
  agile::copy(in_mat,out_mat);
  agile::ifftshift(out_mat);

  cufftResult = CenterdFFT_Helper<TType>::centerdfft(out_mat.data(), out_mat.data(), &fftplan_);

  //Scale
  real_type val_sqrt = (real_type)1.0/std::sqrt(in_mat.getNumRows() * in_mat.getNumColumns());
  agile::scale(val_sqrt, out_mat, out_mat);

  agile::fftshift(out_mat);

  agile::multiplyElementwise(pattern_complex_,out_mat,out_mat);

  return int(cufftResult);
}


}// namespace agile
