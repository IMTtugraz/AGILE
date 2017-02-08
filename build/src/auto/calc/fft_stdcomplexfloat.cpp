
/* This file was generated automatically by CMake. You have to modify '/home/dieheart/workspace/AGILE/src/calc/fft.cpp.in' if you want to make changes. */

#define TType1IsComplex 1
#define TType2IsComplex 

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

#include "calc/fft.cpp"

namespace agile
{
    template void FFT<std::complex<float> >::setfftplan(
    unsigned num_rows, unsigned num_columns);

    template int FFT<std::complex<float> >::CenterdIFFT(
    const GPUMatrix<std::complex<float> >& in_mat, GPUMatrix<std::complex<float> >& out_mat);

    template int FFT<std::complex<float> >::CenterdFFT(
    const GPUMatrix<std::complex<float> >& in_mat, GPUMatrix<std::complex<float> >& out_mat);
    
    template int FFT<std::complex<float> >::CenterdIFFTpattern(
    const GPUMatrix<std::complex<float> >& in_mat, GPUMatrix<std::complex<float> >& out_mat);

    template int FFT<std::complex<float> >::CenterdFFTpattern(
    const GPUMatrix<std::complex<float> >& in_mat, GPUMatrix<std::complex<float> >& out_mat);

    template void FFT<std::complex<float> >::calc_pattern(
    const GPUMatrix<std::complex<float> >& in_mat);
    
    template int FFT<std::complex<float> >::Forward(
    const GPUVector<std::complex<float> >& in_vec, GPUVector<std::complex<float> >& out_vec, 
    unsigned in_offset, unsigned out_offset);
  
    template int FFT<std::complex<float> >::Inverse(
    const GPUVector<std::complex<float> >& in_vec, GPUVector<std::complex<float> >& out_vec, 
    unsigned in_offset, unsigned out_offset);
    
    template int FFT<std::complex<float> >::CenteredForward(
    const GPUVector<std::complex<float> >& in_vec, GPUVector<std::complex<float> >& out_vec, 
    unsigned in_offset, unsigned out_offset);
  
    template int FFT<std::complex<float> >::CenteredInverse(
    const GPUVector<std::complex<float> >& in_vec, GPUVector<std::complex<float> >& out_vec, 
    unsigned in_offset, unsigned out_offset);
} // namespace agile

