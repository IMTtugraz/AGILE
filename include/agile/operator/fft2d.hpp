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

// $Id: fft2d.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_FFT2D_HPP
#define AGILE_OPERATOR_FFT2D_HPP

#include "agile/gpu_config.hpp"
#include "agile/exception.hpp"
#include "agile/gpu_type_traits.hpp"
#include "agile/gpu_vector.hpp"
#include <cufft.h>

namespace agile
{
  //! \brief Struct that changes a C++ type into the corresponding CuFFT type.
  template <typename TType>
  struct to_cufft_type
  {
  };

  //! Partial specialisation for float (--> cufftReal).
  template <>
  struct to_cufft_type<float>
  {
    typedef cufftReal type;
  };

  //! Partial specialisation for std::complex<float> (--> cufftComplex).
  template <>
  struct to_cufft_type<std::complex<float> >
  {
    typedef cufftComplex type;
  };

  //! \brief Orthonormal FFT in two dimensions.
  template <typename TSignalType, typename TSpectralType,
            bool TSignalComplex = is_complex<TSignalType>::value,
            bool TSpectralComplex = is_complex<TSpectralType>::value>
  class FFT2D;

  //! \brief Orthonormal FFT in two dimensions.
  //!
  //! Signal-type is real, spectral-type is complex.
  template <typename TSignalType, typename TSpectralType>
  class FFT2D<TSignalType, TSpectralType, false, true>
  {
    public:
      typedef typename to_cufft_type<TSignalType>::type cufft_in_type;
      typedef typename to_cufft_type<TSpectralType>::type cufft_out_type;

      FFT2D() : m_nx(-1), m_ny(-1) {}

      ~FFT2D()
      {
        cufftDestroy(m_plan);
      }

      void createPlan(unsigned nx, unsigned ny)
      {
        if (nx != m_nx || ny != m_ny)
        {
          m_nx = nx;
          m_ny = ny;
// TODO: CUFFT_R2C for backtransform, too??
          cufftResult result = cufftPlan2d(&m_plan, m_nx, m_ny, CUFFT_R2C);
          AGILE_ASSERT(result == CUFFT_SUCCESS,
                        StandardException::ExceptionMessage(
                          "Could not create FFT plan"));
        }
      }

      void forward(const GPUMatrixPitched<TSignalType>& x,
                   GPUMatrixPitched<TSpectralType>& fft_x)
      {
        // A GPU matrix is aligned to CRS_BLOCK_SIZE elements but the CuFFT
        // takes the elements in a contiguous form in row-major order. Thus,
        // we can only be sure that the CuFFT works as intended, if the aligned
        // and the contiguous storage are the same.
        AGILE_ASSERT(x.getNumColumns() == x.getPitchElements(),
                      StandardException::ExceptionMessage(
                        "FFT of non-aligned matrices not implemented, yet"));

        // create a plan
        createPlan(x.getNumRows(), x.getNumColumns());
        cufftResult result = cufftExecR2C(
          m_plan,
          (typename to_cufft_type<TSignalType>::type*)x.data(),
          (typename to_cufft_type<TSpectralType>::type*)fft_x.data());
          AGILE_ASSERT(result == CUFFT_SUCCESS,
                        StandardException::ExceptionMessage(
                          "Error during FFT procedure"));

//         scale(TSpectralType(1./std::sqrt(x.getNumRows() * x.getNumColumns())),
//               fft_x, fft_x);
      }

      void inverse(const GPUMatrixPitched<TSpectralType>& fft_x,
                   GPUMatrixPitched<TSignalType>& x)
      {
        // assert that we have no zero-padded lines
        AGILE_ASSERT(x.getNumColumns() == x.getPitchElements(),
                      StandardException::ExceptionMessage(
                        "FFT of non-aligned matrices not implemented, yet"));

        // create a plan
        createPlan(fft_x.getNumRows(), fft_x.getNumColumns());
        cufftResult result = cufftExecC2R(
          m_plan,
          (typename to_cufft_type<TSpectralType>::type*)fft_x.data(),
          (typename to_cufft_type<TSignalType>::type*)x.data());
        AGILE_ASSERT(result == CUFFT_SUCCESS,
                        StandardException::ExceptionMessage(
                          "Error during FFT procedure"));

        scale(TSignalType( 1./(fft_x.getNumRows()
                                        * fft_x.getNumColumns())),
              x, x);
      }

    private:
      cufftHandle m_plan;
      unsigned m_nx;
      unsigned m_ny;
  };

  //! \brief Orthonormal FFT in two dimensions.
  //!
  //! Signal-type is complex, spectral-type is complex.
  template <typename TSignalType, typename TSpectralType>
  class FFT2D<TSignalType, TSpectralType, true, true>
  {
    public:
      typedef typename to_cufft_type<TSignalType>::type cufft_in_type;
      typedef typename to_cufft_type<TSpectralType>::type cufft_out_type;

      FFT2D() : m_nx(-1), m_ny(-1) {}

      ~FFT2D()
      {
        cufftDestroy(m_plan);
      }

      void createPlan(unsigned nx, unsigned ny)
      {
        if (nx != m_nx || ny != m_ny)
        {
          m_nx = nx;
          m_ny = ny;
          cufftResult result = cufftPlan2d(&m_plan, m_nx, m_ny, CUFFT_C2C);
          AGILE_ASSERT(result == CUFFT_SUCCESS,
                        StandardException::ExceptionMessage(
                          "Could not create FFT plan"));
        }
      }

      void forward(const GPUMatrixPitched<TSignalType>& x,
                   GPUMatrixPitched<TSpectralType>& fft_x)
      {
        // assert that we have no zero-padded lines
        AGILE_ASSERT(x.getNumColumns() == x.getPitchElements(),
                      StandardException::ExceptionMessage(
                        "FFT of non-aligned matrices not implemented, yet"));

        // create a plan
        createPlan(x.getNumRows(), x.getNumColumns());
        cufftResult result = cufftExecC2C(
          m_plan,
          (typename to_cufft_type<TSignalType>::type*)x.data(),
          (typename to_cufft_type<TSpectralType>::type*)fft_x.data(),
          CUFFT_FORWARD);
        AGILE_ASSERT(result == CUFFT_SUCCESS,
                        StandardException::ExceptionMessage(
                          "Error during FFT procedure"));

      //scale(TSpectralType(1. /std::sqrt(x.getNumRows() * x.getNumColumns())),              fft_x, fft_x);
      }

      void inverse(const GPUMatrixPitched<TSpectralType>& fft_x,
                   GPUMatrixPitched<TSignalType>& x)
      {
        // assert that we have no zero-padded lines
        AGILE_ASSERT(x.getNumColumns() == x.getPitchElements(),
                      StandardException::ExceptionMessage(
                        "FFT of non-aligned matrices not implemented, yet"));

        // create a plan
        createPlan(x.getNumRows(), x.getNumColumns());
        cufftResult result = cufftExecC2C(
          m_plan,
          (typename to_cufft_type<TSpectralType>::type*)fft_x.data(),
          (typename to_cufft_type<TSignalType>::type*)x.data(),
          CUFFT_INVERSE);
        AGILE_ASSERT(result == CUFFT_SUCCESS,
                        StandardException::ExceptionMessage(
                          "Error during FFT procedure"));

        scale(TSignalType(1./(fft_x.getNumRows() * fft_x.getNumColumns())), x, x);
      }

    private:
      cufftHandle m_plan;
      unsigned m_nx;
      unsigned m_ny;
  };

} // namespace agile

#endif // AGILE_OPERATOR_FFT2D_HPP

// End of $Id: fft2d.hpp 476 2011-06-16 08:54:14Z freiberger $.
