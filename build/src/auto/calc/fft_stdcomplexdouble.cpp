
/* This file was generated automatically by CMake. You have to modify '/home2/GIT/AGILE/src/calc/fft.cpp.in' if you want to make changes. */

#define TType1IsComplex 1
#define TType2IsComplex 

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

#include "calc/fft.cpp"

namespace agile
{
    template void FFT<std::complex<double> >::setfftplan(
    unsigned num_rows, unsigned num_columns);

    template int FFT<std::complex<double> >::CenterdIFFT(
    const GPUMatrix<std::complex<double> >& in_mat, GPUMatrix<std::complex<double> >& out_mat);

    template int FFT<std::complex<double> >::CenterdFFT(
    const GPUMatrix<std::complex<double> >& in_mat, GPUMatrix<std::complex<double> >& out_mat);
    
    template int FFT<std::complex<double> >::CenterdIFFTpattern(
    const GPUMatrix<std::complex<double> >& in_mat, GPUMatrix<std::complex<double> >& out_mat);

    template int FFT<std::complex<double> >::CenterdFFTpattern(
    const GPUMatrix<std::complex<double> >& in_mat, GPUMatrix<std::complex<double> >& out_mat);

    template void FFT<std::complex<double> >::calc_pattern(
    const GPUMatrix<std::complex<double> >& in_mat);
    
    template int FFT<std::complex<double> >::Forward(
    const GPUVector<std::complex<double> >& in_vec, GPUVector<std::complex<double> >& out_vec, 
    unsigned in_offset, unsigned out_offset);
  
    template int FFT<std::complex<double> >::Inverse(
    const GPUVector<std::complex<double> >& in_vec, GPUVector<std::complex<double> >& out_vec, 
    unsigned in_offset, unsigned out_offset);
    
    template int FFT<std::complex<double> >::CenteredForward(
    const GPUVector<std::complex<double> >& in_vec, GPUVector<std::complex<double> >& out_vec, 
    unsigned in_offset, unsigned out_offset);
  
    template int FFT<std::complex<double> >::CenteredInverse(
    const GPUVector<std::complex<double> >& in_vec, GPUVector<std::complex<double> >& out_vec, 
    unsigned in_offset, unsigned out_offset);
} // namespace agile

