#ifndef GENKSPACEFOV_H
#define GENKSPACEFOV_H

#include "agile/gpu_matrix.hpp"
#include "agile/calc/fft.hpp"
#include <cufft.h>
#include <iostream>

namespace agile
{
  template <typename TType>
  class KSpaceFOV
  {
    public:

      //! \brief Default constructor.
      //!
      //! The default constructor creates an empty KSpaceFOV Class.
      KSpaceFOV()
      {
        fftobj_ = new agile::FFT<TType>();
      }

      //! \brief Destructor.
      virtual ~KSpaceFOV()
      {
        delete fftobj_;
        fftobj_ = 0;
      }

      int genkspace_fov(const GPUMatrix<TType>& in_mat, GPUMatrix<TType>& out_mat);


    private:

      agile::FFT<TType>* fftobj_;

  };
}

#endif // GENKSPACEFOV_H
