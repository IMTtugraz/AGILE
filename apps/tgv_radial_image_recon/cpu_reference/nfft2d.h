#include <memory.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define PRE_PHI_HUT      (1U<< 0)
#define PRE_FG_PSI       (1U<< 3)
#define PRE_PSI          (1U<< 4)
#define MALLOC_X         (1U<< 6)
#define FFT_OUT_OF_PLACE (1U<< 9)
#define FFTW_INIT        (1U<< 10)


/** Macros for public members inherited by all plan structures. */
#define MACRO_MV_PLAN(float_type)                                             \
  int N_total;                          /**< Total number of Fourier          \
					     coefficients                   */\
  int M_total;                          /**< Total number of samples        */\
                                                                              \
  float_type *f_hat;                    /**< Vector of Fourier coefficients,  \
                                             size is N_total float_types    */\
  float_type *f;                        /**< Vector of samples,               \
				             size is M_total float types    */\
  void (*mv_trafo)(void*);              /**< Pointer to the own transform   */\
  void (*mv_adjoint)(void*);            /**< Pointer to the own adjoint     */\

/** Structure for a NFFT plan */
typedef struct
{
  /* api */
  
  MACRO_MV_PLAN(fftwf_complex)

  int d;                                /**< Dimension, rank                 */
  int *N;                               /**< Multi bandwidth                 */
  float *sigma;                        /**< Oversampling-factor             */
  int *n;                               /**< FFTW length, equal to sigma*N,
					     default is the power of 2 such
					     that \f$2\le\sigma<4\f$         */
  int n_total;                          /**< Total size of FFTW              */
  int m;                                /**< Cut-off parameter of the window
                                             function, default value is
					     6 (KAISER_BESSEL),
					     9 (SINC_POWER),
					     11 (B_SPLINE),
					     12 (GAUSSIAN)                   */
  float *b;                            /**< Shape parameter of the window
                                             function                        */
  int K;                                /**< Number of equispaced samples of
					     the window function for \ref
					     PRE_LIN_PSI                     */

  unsigned nfft_flags;                  /**< Flags for precomputation,
                                             (de)allocation, and FFTW usage,
					     default setting is
					     PRE_PHI_HUT| PRE_PSI|
					     MALLOC_X| MALLOC_F_HAT| MALLOC_F|
					     FFTW_INIT| FFT_OUT_OF_PLACE     */

  unsigned fftw_flags;                  /**< Flags for the FFTW, default is
					     FFTW_ESTIMATE| FFTW_DESTROY_INPUT
					                                     */

  float *x;                            /**< Nodes in time/spatial domain,
					     size is \f$dM\f$ doubles        */

  float MEASURE_TIME_t[3];             /**< Measured time for each step if
					     MEASURE_TIME is set             */

  /* internal*/
  fftwf_plan  my_fftw_plan1;             /**< Forward FFTW plan               */
  fftwf_plan  my_fftw_plan2;             /**< Backward FFTW plan              */

  float **c_phi_inv;                   /**< Precomputed data for the
                                             diagonal matrix \f$D\f$, size is
					     \f$N_0+\hdots+N_{d-1}\f$ doubles*/
  float *psi;                          /**< Precomputed data for the sparse
                                             matrix \f$B\f$, size depends on
					     precomputation scheme           */
  int *psi_index_g;                     /**< Indices in source/target vector
					     for \ref PRE_FULL_PSI           */
  int *psi_index_f;                     /**< Indices in source/target vector
					     for \ref PRE_FULL_PSI           */

  fftwf_complex *g;                    /**< Oversampled vector of samples,
					     size is \ref n_total double
					     complex                         */
  fftwf_complex *g_hat;                /**< Zero-padded vector of Fourier
                                             coefficients, size is
					     \ref n_total fftw_complex     */
  fftwf_complex *g1;                   /**< Input of fftw                   */
  fftwf_complex *g2;                   /**< Output of fftw                  */

  float *spline_coeffs;                /**< Input for de Boor algorithm if
                                             B_SPLINE or SINC_POWER is
                                             defined                         */
} nfft_plan;

void *nfft_malloc(size_t n);
void nfft_free(void *p);
void nfft_init_2d(nfft_plan *ths, int N1, int N2, int M_total);
void nfft_fftw_init(nfft_plan *ths);
void nfft_trafo_2d(nfft_plan *ths);
void nfft_adjoint_2d(nfft_plan *ths);
