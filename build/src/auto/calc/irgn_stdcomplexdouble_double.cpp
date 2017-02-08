
/* This file was generated automatically by CMake. You have to modify '/home/dieheart/workspace/AGILE/src/calc/irgn.cpp.in' if you want to make changes. */

#define TType1IsComplex 1
#define TType2IsComplex 0

#include "agile/gpu_config.hpp"
#include "agile/gpu_type_traits.hpp"
#include <cuda.h>

#include "calc/irgn.cpp"

namespace agile
{
    template void IRGN<std::complex<double>, double >::Init();

    template void IRGN<std::complex<double> ,double >::HighFreqPenalty();
    template void IRGN<std::complex<double> ,double >::Normalize();
    template void IRGN<std::complex<double> ,double >::Iteration();
    template void IRGN<std::complex<double> ,double >::Postprocess();

    template void IRGN<std::complex<double> ,double >::CenterdIFFT(
    const GPUMatrix<std::complex<double> >* in_mat, GPUMatrix<std::complex<double> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<double> ,double >::CenterdFFT(
    const GPUMatrix<std::complex<double> >* in_mat, GPUMatrix<std::complex<double> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<double> ,double >::CenterdIFFTpattern(
    const GPUMatrix<std::complex<double> >* in_mat, GPUMatrix<std::complex<double> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<double> ,double >::CenterdFFTpattern(
    const GPUMatrix<std::complex<double> >* in_mat, GPUMatrix<std::complex<double> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<double> ,double >::ApplyW(
    const GPUMatrix<std::complex<double> >* in_mat, GPUMatrix<std::complex<double> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<double> ,double >::ApplyWH(
    const GPUMatrix<std::complex<double> >* in_mat, GPUMatrix<std::complex<double> >* out_mat, unsigned int num_z);

    template void IRGN<std::complex<double> ,double >::ApplyDFH(
    GPUMatrix<std::complex<double> >* rhs_mat, const GPUMatrix<std::complex<double> >* dx);

    template void IRGN<std::complex<double> ,double >::ApplyM(GPUMatrix<std::complex<double> >* gu, GPUMatrix<std::complex<double> >* gc,
    const GPUMatrix<std::complex<double> >* du, const GPUMatrix<std::complex<double> >* dc);

    template void IRGN<std::complex<double> ,double >::CopyMatrixZ(const GPUMatrix<std::complex<double> >* in_mat,
    GPUMatrix<std::complex<double> >* out_mat, unsigned int num_z);

    template double IRGN<std::complex<double> ,double >::calcLipschitz();

/*
    //Pure Virtual Solve-Methode
    template virtual void  IRGN<std::complex<double> ,double >::Solve(const GPUMatrix<std::complex<double> >* u, const GPUMatrix<std::complex<double> >* c,
      const GPUMatrix<std::complex<double> >* rhs, const GPUMatrix<std::complex<double> >* u0,
      unsigned maxits, double alpha, double beta,
      GPUMatrix<std::complex<double> >* du, GPUMatrix<std::complex<double> >* dc) = 0;
*/

    template double IRGN<std::complex<double> ,double >::randomcalc(int i);

    template bool IRGN<std::complex<double> ,double >::irgn_param_test(IRGN_Params &param);


} // namespace agile

