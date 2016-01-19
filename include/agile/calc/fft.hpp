#ifndef AGILE_FFT_HPP
#define AGILE_FFT_HPP


#include "agile/gpu_matrix.hpp"
#include "agile/gpu_vector.hpp"
#include <cufft.h>
#include <iostream>

namespace agile
{
  template <typename TType>
  class FFT
  {
    public:
      typedef typename agile::to_real_type<TType>::type real_type;

      //! \brief Default constructor.
      //!
      //! The default constructor creates an empty FFT Class.
      FFT()
      {}


      //! \brief Constructor.
      //!
      //! The Constructor creates an FFT.
      FFT(unsigned num_rows, unsigned num_columns)
      {
        setfftplan(num_rows,num_columns);
      }

      //! \brief Destructor.
      virtual ~FFT()
      {
        cufftDestroy(fftplan_);
      }

      agile::GPUMatrix<TType>* get_pattern()
      {
        return &pattern_complex_;
      }

      void set_pattern(const agile::GPUMatrix<TType>& pattern_mat)
      {
        pattern_complex_.resize(pattern_mat.getNumRows(),pattern_mat.getNumColumns());
        agile::copy(pattern_mat,pattern_complex_);
      }

      // ---------------------------------------------------------------
      //! \brief Init()
      //!  Initialisation
      // ---------------------------------------------------------------
      void Init()
      {
        TType one_complex;     //val = 1+0i
        one_complex.real(1);
        one_complex.imag(0);

        ones_complex_vec_cpu_nxns_.clear();
        for(int i=0; i<num_rows_*num_columns_; ++i)
        {
          ones_complex_vec_cpu_nxns_.push_back(one_complex);
        }
        ones_complex_mat_nxns_.assignFromHost(num_rows_, num_columns_, &ones_complex_vec_cpu_nxns_[0]);
      }

      void setfftplan(unsigned num_rows, unsigned num_columns);
      int CenterdIFFT(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat);
      int CenterdFFT(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat);
      int CenterdIFFTpattern(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat);
      int CenterdFFTpattern(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat);

      void calc_pattern(const GPUMatrix<TType>& in_mat);

      int CenteredForward(const GPUVector<TType>& in_vec, GPUVector<TType>& out_vec, unsigned in_offset = 0, unsigned out_offset = -1);
      int CenteredInverse(const GPUVector<TType>& in_vec, GPUVector<TType>& out_vec, unsigned in_offset = 0, unsigned out_offset = -1);
      
      int Forward(const GPUVector<TType>& in_vec, GPUVector<TType>& out_vec, unsigned in_offset = 0, unsigned out_offset = -1);
      int Inverse(const GPUVector<TType>& in_vec, GPUVector<TType>& out_vec, unsigned in_offset = 0, unsigned out_offset = -1);

    private:

      cufftHandle fftplan_;
      unsigned num_rows_;
      unsigned num_columns_;

      agile::GPUMatrix<TType> pattern_complex_;
      agile::GPUMatrix<typename agile::to_real_type<TType>::type> pattern_;
      std::vector<TType> ones_complex_vec_cpu_nxns_;
      GPUMatrix<TType> ones_complex_mat_nxns_;

  };
}
#endif //AGILE_FFT_HPP
