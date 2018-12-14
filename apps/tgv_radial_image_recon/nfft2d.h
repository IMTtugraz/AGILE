#include <memory.h>
#include <assert.h>
#include <math.h>
//#include <cutil.h>
#include <math_constants.h>
//#include <sm_11_atomic_functions.h>
#include <device_atomic_functions.h>
#include <cufft.h>
#include <stdio.h>

//TODO quick fix
//extracted from cutil.h
#define CUDA_SAFE_CALL(call) {                                    \
     cudaError err = call;                                                    \
     if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                 __FILE__, __LINE__, cudaGetErrorString( err) );              \
         exit(EXIT_FAILURE);                                                  \
     } }

#define CUFFT_SAFE_CALL( call) {                                           \
     cufftResult err = call;                                                  \
     if( CUFFT_SUCCESS != err) {                                              \
         fprintf(stderr, "CUFFT error in file '%s' in line %i.\n",            \
                 __FILE__, __LINE__);                                         \
         exit(EXIT_FAILURE);                                                  \
     } }

#define CUT_CHECK_ERROR(errorMessage) {                                    \
     cudaError_t err = cudaGetLastError();                                    \
     if( cudaSuccess != err) {                                                \
         fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                 errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
         exit(EXIT_FAILURE);                                                  \
     } }       
//end

typedef float complex_type[2];

typedef struct nfft_plan_t
{
  int alloc_flags;       // flags for memory allocation

  int N[2];              // source dimensions
  int n[2];              // fft dimensions
  int log2n[2];          // log2 of fft dimensions
  int M_total, coils;    // total number of nodes and coils
  size_t coil_stride;    // coil stride (in bytes)
  int N_total, n_total;  // total number of (FFT) coefficients
  
  complex_type *f_hat;  // space data
  complex_type *g;      // fft data
  complex_type *f;      // Fourier data
  float *x;              // node data

  float sigma[2], b[2];                 // oversampling factors
  float *(phi_hut_inv[2]);              // precomputed phi_hut_inv
  
  cufftHandle cufft_plan;    // CUFFT plan
  
  int stages;                // number of stages
  int2 *stageinfo;           // stage info data
  int bins;                  // number of bins
  int2 *bincoords;           // array containing bin coordinates
  int *binindex;             // array containing bin indices
} nfft_plan;



void *nfft_malloc(size_t n);
void nfft_free(void *p);
void nfft_init_2d(nfft_plan *ths, int N1, int N2, int M_total, int coils,
		  size_t coil_stride, complex_type *f_hat, complex_type *f, 
		  complex_type *g, float *x);
void nfft_finalize_2d(nfft_plan *ths);
void nfft_trafo_2d(nfft_plan *ths);
void nfft_adjoint_2d(nfft_plan *ths);
