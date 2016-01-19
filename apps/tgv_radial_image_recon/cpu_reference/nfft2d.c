#include <nfft2d.h>

// bessel i0 function

/*              chbevl.c
 *
 *  Evaluate Chebyshev series
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * double x, y, coef[N], chebevl();
 *
 * y = chbevl( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the series
 *
 *        N-1
 *         - '
 *  y  =   >   coef[i] T (x/2)
 *         -            i
 *        i=0
 *
 * of Chebyshev polynomials Ti at argument x/2.
 *
 * Coefficients are stored in reverse order, i.e. the zero
 * order term is last in the array.  Note N is the number of
 * coefficients, not the order.
 *
 * If coefficients are for the interval a to b, x must
 * have been transformed to x -> 2(2x - b - a)/(b-a) before
 * entering the routine.  This maps x from (a, b) to (-1, 1),
 * over which the Chebyshev polynomials are defined.
 *
 * If the coefficients are for the inverted interval, in
 * which (a, b) is mapped to (1/b, 1/a), the transformation
 * required is x -> 2(2ab/x - b - a)/(b-a).  If b is infinity,
 * this becomes x -> 4a/x - 1.
 *
 *
 *
 * SPEED:
 *
 * Taking advantage of the recurrence properties of the
 * Chebyshev polynomials, the routine requires one more
 * addition per loop than evaluating a nested polynomial of
 * the same degree.
 *
 */

/*              chbevl.c  */

/*
  Cephes Math Library Release 2.0:  April, 1987
  Copyright 1985, 1987 by Stephen L. Moshier
  Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

static float chbevl(float x, float array[], int n)
{
  float b0, b1, b2, *p;
  int i;
  
  p = array;
  b0 = *p++;
  b1 = 0.0f;
  i = n - 1;
  
  do
    {
      b2 = b1;
      b1 = b0;
      b0 = x * b1  -  b2  + *p++;
    }
  while( --i );
  
  return( 0.5f*(b0-b2) );
}

#define UNK 1
/* If you define UNK, then be sure to set BIGENDIAN properly. */
#ifdef FLOAT_WORDS_BIGENDIAN
#define BIGENDIAN 1
#else
#define BIGENDIAN 0
#endif

/*              i0.c
 *
 *  Modified Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, i0();
 *
 * y = i0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order zero of the
 * argument.
 *
 * The function is defined as i0(x) = j0( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,30         6000       8.2e-17     1.9e-17
 *    IEEE      0,30        30000       5.8e-16     1.4e-16
 *
 */

/*
  Cephes Math Library Release 2.8:  June, 2000
  Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

/* Chebyshev coefficients for exp(-x) I0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I0(x) } = 1.
 */

#ifdef UNK
static float A[] =
  {
    -4.41534164647933937950E-18,
    3.33079451882223809783E-17,
    -2.43127984654795469359E-16,
    1.71539128555513303061E-15,
    -1.16853328779934516808E-14,
    7.67618549860493561688E-14,
    -4.85644678311192946090E-13,
    2.95505266312963983461E-12,
    -1.72682629144155570723E-11,
    9.67580903537323691224E-11,
    -5.18979560163526290666E-10,
    2.65982372468238665035E-9,
    -1.30002500998624804212E-8,
    6.04699502254191894932E-8,
    -2.67079385394061173391E-7,
    1.11738753912010371815E-6,
    -4.41673835845875056359E-6,
    1.64484480707288970893E-5,
    -5.75419501008210370398E-5,
    1.88502885095841655729E-4,
    -5.76375574538582365885E-4,
    1.63947561694133579842E-3,
    -4.32430999505057594430E-3,
    1.05464603945949983183E-2,
    -2.37374148058994688156E-2,
    4.93052842396707084878E-2,
    -9.49010970480476444210E-2,
    1.71620901522208775349E-1,
    -3.04682672343198398683E-1,
    6.76795274409476084995E-1
  };
#endif

/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
 */

#ifdef UNK
static float B[] =
  {
    -7.23318048787475395456E-18,
    -4.83050448594418207126E-18,
    4.46562142029675999901E-17,
    3.46122286769746109310E-17,
    -2.82762398051658348494E-16,
    -3.42548561967721913462E-16,
    1.77256013305652638360E-15,
    3.81168066935262242075E-15,
    -9.55484669882830764870E-15,
    -4.15056934728722208663E-14,
    1.54008621752140982691E-14,
    3.85277838274214270114E-13,
    7.18012445138366623367E-13,
    -1.79417853150680611778E-12,
    -1.32158118404477131188E-11,
    -3.14991652796324136454E-11,
    1.18891471078464383424E-11,
    4.94060238822496958910E-10,
    3.39623202570838634515E-9,
    2.26666899049817806459E-8,
    2.04891858946906374183E-7,
    2.89137052083475648297E-6,
    6.88975834691682398426E-5,
    3.36911647825569408990E-3,
    8.04490411014108831608E-1
  };
#endif

/** Modified Bessel function of order zero.
 *  Cephes Math Library Release 2.8:  June, 2000
 *  Copyright 1984, 1987, 2000 by Stephen L. Moshier
 */
float nfft_i0(float x)
{
  float y;
  
  if( x < 0 )
    x = -x;
  if( x <= 8.0f )
    {
      y = (x/2.0f) - 2.0f;
      return( expf(x) * chbevl( y, A, 30 ) );
    }
  
  return(  expf(x) * chbevl( 32.0f/x - 2.0f, B, 25 ) / sqrtf(x) );
}

// define Kaiser-Bessel window functions

/** Formerly known to be an irrational number.
 */
#define PI 3.141592653589793238462643383279502884197169399375105820974944592f
#define PI2 6.283185307179586476925286766559005768394338798750211641949889185f
#define PI4 12.56637061435917295385057353311801153678867759750042328389977837f

#define THS_M 6.0f
#define THS_M_SQR 36.0f

inline float sqrf(float x) {
  return (x*x);
}

#define THS_B(d) ((float)PI*(2.0f-1.0f/ths->sigma[d]))

inline float phi_hut(float k, int d, float b, int n) {
  return ((float)nfft_i0(THS_M*sqrtf(sqrf(b)-sqrf(PI2*k/n))));
}

#define PHI_HUT(k,d) phi_hut(k,d,ths_b##d,n##d)

inline float phi(float x, int d, float b, int n) {
  float val;
  
  val = THS_M_SQR - sqrf(x*n);
  if (val > 0) return (sinhf(b*sqrtf(val))/(PI*sqrtf(val)));
  
  val = n*n - sqrf(x*THS_M);
  if (val > 0) return (sinf(b*sqrtf(val))/(PI*sqrtf(val)));
  
  return(1.0);
}

#define PHI(x,d) phi(x,d,ths_b##d,n##d)

/// modified nfft_init routines (takes from nfft source code)

/** Computes /f$n\ge N/f$ such that /f$n=2^j,\, j\in\mathhb{N}_0/f$.
 */
int nfft_next_power_of_2(int N)
{
  int n,i,logn;
  int N_is_not_power_of_2=0;

  if (N == 0)
  {
    return 1;
  }
  else
  {
    n=N;
    logn=0;
    while (n != 1)
    {
      if (n%2 == 1)
      {
          N_is_not_power_of_2=1;
      }
      n = n/2;
      logn++;
    }

    if (!N_is_not_power_of_2)
    {
      logn--;
    }

    for (i = 0; i <= logn; i++)
    {
      n = n*2;
    }

    return n;
  }
}

/** Computes integer /f$\prod_{t=0}^{d-1} v_t/f$.
 */
int nfft_prod_int(int *vec, int d)
{
  int t, prod;

  prod=1;
  for(t=0; t<d; t++)
    prod *= vec[t];

  return prod;
}

void *nfft_malloc(size_t n)
{
  void *p;

  /*if (nfft_malloc_hook)
    return;  nfft_malloc_hook(n);*/

  if (n == 0)
    n = 1;

  p = fftwf_malloc(n);

  if (!p)
  {
    fprintf(stderr,"nfft_malloc: n = %d.\n", (int)n);
    //nfft_die("nfft_malloc: out of memory\n");
  }

  return p;
}

void nfft_free(void *p)
{
  if (p)
  {
    /*if (nfft_free_hook)
    {
      nfft_free_hook(p);
      return;
      }*/
    fftwf_free(p);
  }
}


static void nfft_init_help(nfft_plan *ths)
{
  int t;                                /**< index over all dimensions       */
  int lprod;                            /**< 'bandwidth' of matrix B         */

  ths->N_total=nfft_prod_int(ths->N, ths->d);
  ths->n_total=nfft_prod_int(ths->n, ths->d);
  
  ths->sigma = (float *) nfft_malloc(ths->d*sizeof(float));
  for(t = 0;t < ths->d; t++)
    ths->sigma[t] = ((float)ths->n[t])/ths->N[t];
  
  if(ths->nfft_flags & MALLOC_X)
    ths->x = (float *)nfft_malloc(ths->d*ths->M_total*sizeof(float));
   
  ths->mv_trafo = (void (*) (void* ))nfft_trafo_2d;
  ths->mv_adjoint = (void (*) (void* ))nfft_adjoint_2d;
}

void nfft_fftw_init(nfft_plan *ths) {
  if(ths->nfft_flags & FFTW_INIT)
    { 
      ths->my_fftw_plan1 = fftwf_plan_dft(ths->d, ths->n, ths->g1, ths->g2, 
					  FFTW_FORWARD, ths->fftw_flags);
      ths->my_fftw_plan2 = fftwf_plan_dft(ths->d, ths->n, ths->g2, ths->g1,
					  FFTW_BACKWARD, ths->fftw_flags);
    }
}

void nfft_init(nfft_plan *ths, int d, int *N, int M_total)
{
  int t;                                /**< index over all dimensions       */

  ths->d = d;

  ths->N=(int*) nfft_malloc(d*sizeof(int));
  for(t = 0;t < d; t++)
    ths->N[t] = N[t];

  ths->M_total = M_total;

  ths->n = (int*) nfft_malloc(d*sizeof(int));
  for(t = 0;t < d; t++)
    ths->n[t] = 2*nfft_next_power_of_2(ths->N[t]);

  ths->nfft_flags = PRE_PHI_HUT| PRE_PSI| MALLOC_X| 
    FFTW_INIT| FFT_OUT_OF_PLACE;
  ths->fftw_flags= FFTW_ESTIMATE| FFTW_DESTROY_INPUT;

  nfft_init_help(ths);
}

void nfft_init_2d(nfft_plan *ths, int N1, int N2, int M_total)
{
  int N[2];

  N[0]=N1;
  N[1]=N2;
  nfft_init(ths,2,N,M_total);
}

/// end modified nfft


/** computes 2m+2 indices for the matrix B
 */
inline void nfft_uo(float xj, int n, int *up, int *op)
{
  int c = lrint(xj*n);

  if(xj < 0)
    {
      (*up) = c-1-THS_M;
      (*op) = c  +THS_M;
    }
  else
    {
      (*up) = c  - THS_M;
      (*op) = c+1+ THS_M;
    }
}

// wrap with n
inline void nfft_uo2(float x, int n, int *u, int *o)
{
  int c = lrint(x * n);

  if(x < 0)
    {
      *u=(c-1-(int)THS_M+n)%n;
      *o=(c+  (int)THS_M+n)%n;
    }
  else
    {
      *u=(c  -(int)THS_M+n)%n;
      *o=(c+1+(int)THS_M+n)%n;
    }
}

// evaluate convolution in xj
static void nfft_trafo_2d_compute(float _Complex *fj, const float _Complex *g,
				  const float *psij_const0, const float *psij_const1,
				  const float xj0, const float xj1,
				  const int n0, const int n1)
{
  int u0,o0,l0,u1,o1,l1;
  int m;
  const float _Complex *gj;
  const float *psij0, *psij1;

  m = THS_M;
  
  psij0=psij_const0;
  psij1=psij_const1;

  nfft_uo2(xj0,n0,&u0,&o0);
  nfft_uo2(xj1,n1,&u1,&o1);
  
  *fj=0;
  
  if(u0<o0)
    if(u1<o1)
      // no periodic wrapping
      for(l0=0; l0<=2*m+1; l0++,psij0++)
	{
	  psij1=psij_const1;
	  gj=g+(u0+l0)*n1+u1;
	  for(l1=0; l1<=2*m+1; l1++)
	    (*fj) += (*psij0) * (*psij1++) * (*gj++);
	}
    else
      // periodic wrapping in dim 1
      for(l0=0; l0<=2*m+1; l0++,psij0++)
	{
	  psij1=psij_const1;
	  gj=g+(u0+l0)*n1+u1;
	  for(l1=0; l1<2*m+1-o1; l1++)
	    (*fj) += (*psij0) * (*psij1++) * (*gj++);
	  gj=g+(u0+l0)*n1;
	  for(l1=0; l1<=o1; l1++)
	    (*fj) += (*psij0) * (*psij1++) * (*gj++);
	}
  else
    if(u1<o1)
      {
	// periodic wrapping in dim 0
	for(l0=0; l0<2*m+1-o0; l0++,psij0++)
	  {
	    psij1=psij_const1;
	    gj=g+(u0+l0)*n1+u1;
	    for(l1=0; l1<=2*m+1; l1++)
	      (*fj) += (*psij0) * (*psij1++) * (*gj++);
	  }
	for(l0=0; l0<=o0; l0++,psij0++)
	  {
	    psij1=psij_const1;
	    gj=g+l0*n1+u1;
	    for(l1=0; l1<=2*m+1; l1++)
	      (*fj) += (*psij0) * (*psij1++) * (*gj++);
	  }
      }
    else
      {
	// periodic wrapping in dim 0 and 1
	for(l0=0; l0<2*m+1-o0; l0++,psij0++)
	  {
	    psij1=psij_const1;
	    gj=g+(u0+l0)*n1+u1;
	    for(l1=0; l1<2*m+1-o1; l1++)
	      (*fj) += (*psij0) * (*psij1++) * (*gj++);
	    gj=g+(u0+l0)*n1;
	    for(l1=0; l1<=o1; l1++)
	      (*fj) += (*psij0) * (*psij1++) * (*gj++);
	  }
	for(l0=0; l0<=o0; l0++,psij0++)
	  {
	    psij1=psij_const1;
	    gj=g+l0*n1+u1;
	    for(l1=0; l1<2*m+1-o1; l1++)
	      (*fj) += (*psij0) * (*psij1++) * (*gj++);
	    gj=g+l0*n1;
	    for(l1=0; l1<=o1; l1++)
	      (*fj) += (*psij0) * (*psij1++) * (*gj++);
	  }
      }
}

static void nfft_trafo_2d_B(nfft_plan *ths)
{
  int n0,N0,n1,N1,u,o,j,M,l,m;
  float _Complex *fj, *g;
  float ths_b0, ths_b1;
  float psij_const[2*(2*(int)THS_M+2)];
  float *xj, xj0, xj1;

  // dimensions
  N0=ths->N[0];
  n0=ths->n[0];
  N1=ths->N[1];
  n1=ths->n[1];
  M=ths->M_total;
  m=THS_M;

  // "set" ths->b
  ths_b0 = THS_B(0);
  ths_b1 = THS_B(1);

  g=ths->g;

  /* no precomputed psi at all */
  // for each node
  for(j=0,fj=ths->f,xj=ths->x; j<M; j++,fj++,xj+=2)
    {
      xj0 = xj[0]; xj1 = xj[1]; 

      // dim 0
      // get margins 
      nfft_uo(xj0,n0,&u,&o);
      // compute constants
      for(l=0;l<=2*m+1;l++)
	psij_const[l]=(PHI(xj0-((float)((u+l)))/n0,0));
      
      // dim 1
      // get margins 
      nfft_uo(xj1,n1,&u,&o);
      // compute constants
      for(l=0;l<=2*m+1;l++)
	psij_const[2*m+2+l]=(PHI(xj1-((float)((u+l)))/n1,1));

      // evaluate convolution 
      nfft_trafo_2d_compute(fj, g, psij_const, psij_const+2*(int)THS_M+2, 
			    xj0, xj1, n0, n1);
    }
}

void nfft_trafo_2d(nfft_plan *ths)
{
  int k0,k1,n0,n1,N0,N1;
  float _Complex *g_hat; 
  float _Complex *f_hat;
  float ths_b0, ths_b1;
  float *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12;
  float ck01, ck02, ck11, ck12;
  float _Complex *g_hat11,*f_hat11,*g_hat21,*f_hat21,
    *g_hat12,*f_hat12,*g_hat22,*f_hat22;

  // pointers for FFTW-data
  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  // dimensions
  N0=ths->N[0];
  N1=ths->N[1];
  n0=ths->n[0];
  n1=ths->n[1];
  
  // "set" ths->b
  ths_b0 = THS_B(0);
  ths_b1 = THS_B(1);

  f_hat=ths->f_hat;
  g_hat=ths->g_hat;

  memset(ths->g_hat,0,ths->n_total*sizeof(float _Complex));
  
  // copy/multiply data with inverse phi_hut 
  // with fftshift
  for(k0=0;k0<N0/2;k0++)
    {
      ck01=1.0f/(PHI_HUT(k0-N0/2,0));
      ck02=1.0f/(PHI_HUT(k0,0));
      for(k1=0;k1<N1/2;k1++)
	{
	  ck11=1.0f/(PHI_HUT(k1-N1/2,1));
	  ck12=1.0f/(PHI_HUT(k1,1));
	  g_hat[(n0-N0/2+k0)*n1+n1-N1/2+k1] = f_hat[k0*N1+k1]             * ck01 * ck11;
	  g_hat[k0*n1+n1-N1/2+k1]           = f_hat[(N0/2+k0)*N1+k1]      * ck02 * ck11;
	  g_hat[(n0-N0/2+k0)*n1+k1]         = f_hat[k0*N1+N1/2+k1]        * ck01 * ck12;
	  g_hat[k0*n1+k1]                   = f_hat[(N0/2+k0)*N1+N1/2+k1] * ck02 * ck12;
	}
    }

  // compute fft
  fftwf_execute(ths->my_fftw_plan1);

  // convolve with phi each node
  nfft_trafo_2d_B(ths);
}

static void nfft_adjoint_2d_compute(const float _Complex *fj, float _Complex *g,
				    const float *psij_const0, const float *psij_const1,
				    const float xj0, const float xj1,
				    const int n0, const int n1)
{
  int u0,o0,l0,u1,o1,l1,m;
  float _Complex *gj;
  const float *psij0, *psij1;

  psij0=psij_const0;
  psij1=psij_const1;

  nfft_uo2(xj0,n0,&u0,&o0);
  nfft_uo2(xj1,n1,&u1,&o1);
  
  m = THS_M;

  if(u0<o0)
    if(u1<o1)
      // no periodic wrapping 
      for(l0=0; l0<=2*m+1; l0++,psij0++)
	{
	  psij1=psij_const1;
	  gj=g+(u0+l0)*n1+u1;
	  for(l1=0; l1<=2*m+1; l1++)
	    (*gj++) += (*psij0) * (*psij1++) * (*fj);
	}
    else
      // periodic wrapping in dim 1
      for(l0=0; l0<=2*m+1; l0++,psij0++)
	{
	  psij1=psij_const1;
	  gj=g+(u0+l0)*n1+u1;
	  for(l1=0; l1<2*m+1-o1; l1++)
	    (*gj++) += (*psij0) * (*psij1++) * (*fj);
	  gj=g+(u0+l0)*n1;
	  for(l1=0; l1<=o1; l1++)
	    (*gj++) += (*psij0) * (*psij1++) * (*fj);
	}
  else
    if(u1<o1)
      {
	// periodic wrapping in dim 0
	for(l0=0; l0<2*m+1-o0; l0++,psij0++)
	  {
	    psij1=psij_const1;
	    gj=g+(u0+l0)*n1+u1;
	    for(l1=0; l1<=2*m+1; l1++)
	      (*gj++) += (*psij0) * (*psij1++) * (*fj);
	  }
	for(l0=0; l0<=o0; l0++,psij0++)
	  {
	    psij1=psij_const1;
	    gj=g+l0*n1+u1;
	    for(l1=0; l1<=2*m+1; l1++)
	      (*gj++) += (*psij0) * (*psij1++) * (*fj);
	  }
      }
    else
      {
	// periodic wrapping in dim 0 and dim 1
	for(l0=0; l0<2*m+1-o0; l0++,psij0++)
	  {
	    psij1=psij_const1;
	    gj=g+(u0+l0)*n1+u1;
	    for(l1=0; l1<2*m+1-o1; l1++)
	      (*gj++) += (*psij0) * (*psij1++) * (*fj);
	    gj=g+(u0+l0)*n1;
	    for(l1=0; l1<=o1; l1++)
	      (*gj++) += (*psij0) * (*psij1++) * (*fj);
	  }
	for(l0=0; l0<=o0; l0++,psij0++)
	  {
	    psij1=psij_const1;
	    gj=g+l0*n1+u1;
	    for(l1=0; l1<2*m+1-o1; l1++)
	      (*gj++) += (*psij0) * (*psij1++) * (*fj);
	    gj=g+l0*n1;
	    for(l1=0; l1<=o1; l1++)
	      (*gj++) += (*psij0) * (*psij1++) * (*fj);
	  }
      }
}

static void nfft_adjoint_2d_B(nfft_plan *ths)
{
  int n0,N0,n1,N1,u,o,j,M,l,m;
  float _Complex *fj,*g;
  float *xj;
  float psij_const[2*(2*(int)THS_M+2)];
  float ths_b0, ths_b1;
  float xj0, xj1;

  N0=ths->N[0];
  n0=ths->n[0];
  N1=ths->N[1];
  n1=ths->n[1];
  M=ths->M_total;
  m=THS_M;

  // "set" ths->b
  ths_b0 = THS_B(0);
  ths_b1 = THS_B(1);

  g=ths->g;
  memset(g,0,ths->n_total*sizeof(float _Complex));

  /* no precomputed psi at all */
  // for each node
  for(j=0,fj=ths->f,xj=ths->x;j<M;j++,fj++,xj+=2)
    {
      xj0 = xj[0]; xj1 = xj[1];

      // compute weights dim 0
      nfft_uo(xj0,n0,&u,&o);
      for(l=0;l<=2*m+1;l++)
	psij_const[l]=(PHI(xj0-((float)((u+l)))/n0,0));
      
      // compute weights dim 1
      nfft_uo(xj1,n1,&u,&o);
      for(l=0;l<=2*m+1;l++)
	psij_const[2*m+2+l]=(PHI(xj1-((float)((u+l)))/n1,1));
      
      nfft_adjoint_2d_compute(fj, g, psij_const, psij_const+2*m+2, xj0, xj1, 
			      n0, n1);
    }
}

void nfft_adjoint_2d(nfft_plan *ths)
{
  int k0,k1,n0,n1,N0,N1;
  float _Complex *g_hat,*f_hat;
  float *c_phi_inv01, *c_phi_inv02, *c_phi_inv11, *c_phi_inv12;
  float ths_b0, ths_b1;
  float ck01, ck02, ck11, ck12;
  float _Complex *g_hat11,*f_hat11,*g_hat21,*f_hat21,
    *g_hat12,*f_hat12,*g_hat22,*f_hat22;

  // FFTW data
  ths->g_hat=ths->g1;
  ths->g=ths->g2;

  // dimensions
  N0=ths->N[0];
  N1=ths->N[1];
  n0=ths->n[0];
  n1=ths->n[1];

  // "set" ths->b
  ths_b0 = THS_B(0);
  ths_b1 = THS_B(1);

  f_hat=ths->f_hat;
  g_hat=ths->g_hat;

  // convolve adjoint
  nfft_adjoint_2d_B(ths);
  
  // adjoint (inverse) fft
  fftwf_execute(ths->my_fftw_plan2);
  
  // copy multiply with phi_hut
  // inverse fftshift
  for(k0=0;k0<N0/2;k0++)
    {
      ck01=1.0f/(PHI_HUT(k0-N0/2,0));
      ck02=1.0f/(PHI_HUT(k0,0));
      for(k1=0;k1<N1/2;k1++)
	{
	  ck11=1.0f/(PHI_HUT(k1-N1/2,1));
	  ck12=1.0f/(PHI_HUT(k1,1));
	  f_hat[k0*N1+k1]             = g_hat[(n0-N0/2+k0)*n1+n1-N1/2+k1] * ck01 * ck11;
	  f_hat[(N0/2+k0)*N1+k1]      = g_hat[k0*n1+n1-N1/2+k1]           * ck02 * ck11;
	  f_hat[k0*N1+N1/2+k1]        = g_hat[(n0-N0/2+k0)*n1+k1]         * ck01 * ck12;
	  f_hat[(N0/2+k0)*N1+N1/2+k1] = g_hat[k0*n1+k1]                   * ck02 * ck12;
	}
    }
}
