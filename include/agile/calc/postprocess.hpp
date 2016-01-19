#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include "agile/gpu_matrix.hpp"
#include "agile/calc/fft.hpp"
#include <iostream>

namespace agile
{
  template <typename TType, typename TType_real>
  class PostProcess
  {
    public:

      //! \brief Default Constructor.
      PostProcess()
      {}

      //! \brief Constructor.
      //!
      //! The Constructor
      PostProcess(unsigned num_rows, unsigned num_columns, unsigned num_coils = 1)
      {
        num_rows_ = num_rows;
        num_columns_ = num_columns;
        num_coils_ = num_coils;

        init();
      }

      //! \brief Destructor.
      virtual ~PostProcess()
      {

      }

      void set_size(unsigned num_rows, unsigned num_columns, unsigned num_coils = 1)
      {
        if((num_rows != num_rows_) || (num_columns != num_columns_))
        {
          num_rows_ = num_rows;
          num_columns_ = num_columns;
          init();
        }
        num_coils_ = num_coils;


      }

      void calc_abs(const GPUMatrix<TType>* in_mat, GPUMatrix<TType_real>& out_mat);
      void calc_phase(const GPUMatrix<TType>* in_mat, GPUMatrix<TType_real>& out_mat);
      void calc_imag(const GPUMatrix<TType>* in_mat, GPUMatrix<TType_real>& out_mat);
      void calc_real(const GPUMatrix<TType>* in_mat, GPUMatrix<TType_real>& out_mat);
      void init();

    private:

      agile::FFT<TType>* fftobj;
      unsigned num_coils_;
      unsigned num_rows_;
      unsigned num_columns_;

      std::vector<typename to_real_type<TType>::type> zero_vec_cpu_;
      agile::GPUMatrix<typename agile::to_real_type<TType>::type> zero_vec_gpu_;

  };
}

#endif //POSTPROCESS_H
