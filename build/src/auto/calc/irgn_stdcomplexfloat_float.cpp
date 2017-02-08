
/* This file was generated automatically by CMake. You have to modify '/home/dieheart/workspace/AGILE/src/calc/irgn.cpp.in' if you want to make changes. */

#define TType1IsComplex 1
#define TType2IsComplex 0

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

#include "calc/irgn.cpp"

namespace agile
{
    template void IRGN<std::complex<float>, float >::Init();

    template void IRGN<std::complex<float> ,float >::HighFreqPenalty();
    template void IRGN<std::complex<float> ,float >::Normalize();
    template void IRGN<std::complex<float> ,float >::Iteration();
    template void IRGN<std::complex<float> ,float >::Postprocess();

    template void IRGN<std::complex<float> ,float >::CenterdIFFT(
    const GPUMatrix<std::complex<float> >* in_mat, GPUMatrix<std::complex<float> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<float> ,float >::CenterdFFT(
    const GPUMatrix<std::complex<float> >* in_mat, GPUMatrix<std::complex<float> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<float> ,float >::CenterdIFFTpattern(
    const GPUMatrix<std::complex<float> >* in_mat, GPUMatrix<std::complex<float> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<float> ,float >::CenterdFFTpattern(
    const GPUMatrix<std::complex<float> >* in_mat, GPUMatrix<std::complex<float> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<float> ,float >::ApplyW(
    const GPUMatrix<std::complex<float> >* in_mat, GPUMatrix<std::complex<float> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<float> ,float >::ApplyWH(
    const GPUMatrix<std::complex<float> >* in_mat, GPUMatrix<std::complex<float> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<float> ,float >::ApplyDFH(
    GPUMatrix<std::complex<float> >* rhs_mat, const GPUMatrix<std::complex<float> >* dx);

    template void IRGN<std::complex<float> ,float >::ApplyM(GPUMatrix<std::complex<float> >* gu, GPUMatrix<std::complex<float> >* gc,
    const GPUMatrix<std::complex<float> >* du, const GPUMatrix<std::complex<float> >* dc);

    template void IRGN<std::complex<float> ,float >::CopyMatrixZ(const GPUMatrix<std::complex<float> >* in_mat,
    GPUMatrix<std::complex<float> >* out_mat, unsigned int num_z);

    template float IRGN<std::complex<float> ,float >::calcLipschitz();

/*
    //Pure Virtual Solve-Methode
    template virtual void  IRGN<std::complex<float> ,float >::Solve(const GPUMatrix<std::complex<float> >* u, const GPUMatrix<std::complex<float> >* c,
      const GPUMatrix<std::complex<float> >* rhs, const GPUMatrix<std::complex<float> >* u0,
      unsigned maxits, float alpha, float beta,
      GPUMatrix<std::complex<float> >* du, GPUMatrix<std::complex<float> >* dc) = 0;
*/

    template float IRGN<std::complex<float> ,float >::randomcalc(int i);

    template bool IRGN<std::complex<float> ,float >::irgn_param_test(IRGN_Params &param);


} // namespace agile

