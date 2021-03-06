#include <nfft2d.h>

// constants

#define KB_M 6.0f
#define KB_M_SQR 36.0f

#define THS_M 6
#define THS_WINDOW_LEN (2*THS_M+2)
#define THS_WINDOW_SPACE 256
#define THS_BLOCK_LEN 16
#define PHI_HUT_SAFE_FRAME 32

#define NFFT_NODE_BIN_SIZE 32

#define NFFT_INIT_ALLOC_F 1
#define NFFT_INIT_ALLOC_F_HAT 2
#define NFFT_INIT_ALLOC_G 4
#define NFFT_INIT_ALLOC_BIN 8
#define NFFT_INIT_ALLOC_CUFFT 16

typedef float float_complex[2];

// bessel i0 function

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*/

/* Chebyshev coefficients for exp(-x) I0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) I0(x) } = 1.
 */
double bessel_i0_A[] =
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

double bessel_i0_B[] =
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

/*
  Cephes Math Library Release 2.0:  April, 1987
  Copyright 1985, 1987 by Stephen L. Moshier
  Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

double chbevl(double x, double array[], int n)
{
  double b0, b1, b2;
  int i, j;
  
  b0 = array[0];
  b1 = 0.0;
  i = n - 1;
  j = 1;
  
  do
    {
      b2 = b1;
      b1 = b0;
      b0 = x * b1  -  b2  + array[j++];
    }
  while( --i );
  
  return( 0.5*(b0-b2) );
}

/** Modified Bessel function of order zero.
 *  Cephes Math Library Release 2.8:  June, 2000
 *  Copyright 1984, 1987, 2000 by Stephen L. Moshier
 */

double nfft_i0(double x)
{
  float y;
  
  x = fabs(x);
  if( x <= 8.0 )
    {
      y = (x/2.0) - 2.0;
      return( exp(x) * chbevl( y, bessel_i0_A, 30 ) );
    }
  
  return(  exp(x) * chbevl( 32.0/x - 2.0, bessel_i0_B, 25 ) / sqrt(x) );
}

/* END: Cephes Math Library Release 2.8 */


// define Kaiser-Bessel window functions

/*
 * Bases on NFFT3 library source code 
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
 */

__device__ inline float sqr(float x) {
  return (x*x);
}

#define THS_B(d) (CUDART_PI*(2.0-1.0/ths->sigma[d]))

double phi_hut(double k, int d, double b, int n) {
  return (nfft_i0(KB_M*sqrt(b*b-4.0*CUDART_PI*CUDART_PI*k*k/((float)(n*n)))));
}

#define PHI_HUT(k,d) phi_hut(k,d,(CUDART_PI*(2.0-1.0/ths->sigma[d])),ths->n[d])

__device__ inline float phi(float x, int d, float b, int n) {
  float val;
  
  val = KB_M_SQR - sqr(x*n);
  if (val > 0) return (sinh(b*sqrt(val))/(CUDART_PI_F*sqrt(val)));
  
  val = n*n - sqr(x*KB_M);
  if (val > 0) return (sin(b*sqrt(val))/(CUDART_PI_F*sqrt(val)));
  
  return(1.0);
}

#define PHI(x,d) phi(x,d,ths_b[d],n[d])

/*
 * Modified nfft_init routines
 * 
 * Base on NFFT3 library source code 
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
 */

/** Computes /f$n\ge N/f$ such that /f$n=2^j,\, j\in\mathhb{N}_0/f$.
 */
int nfft_next_power_of_2(int N, int *logn)
{
  int n,i;
  int N_is_not_power_of_2=0;
  
  if (N == 0)
    {
      return 1;
    }
  else
    {
      n=N;
      *logn = 0;
      while (n != 1)
	{
	  if (n%2 == 1)
	    {
	      N_is_not_power_of_2=1;
	    }
	  n = n/2;
	  (*logn)++;
	}
      
      if (!N_is_not_power_of_2)
	{
	  (*logn)--;
	}
      
      for (i = 0; i <= *logn; i++)
	{
	  n = n*2;
	}
      
      (*logn)++;
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

  p = malloc(n);

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
    free(p);
  }
}

/** initialisation of direct transform
 */
static void nfft_precompute_phi_hut(nfft_plan *ths)
{
  float *phi_hut_inv_host;
  int ks, t, numel;
  
  // for each dimension
  for (t=0;t<2;t++) {
    numel = (ths->N[t] + PHI_HUT_SAFE_FRAME);
    
    // allocate temporary space
    phi_hut_inv_host = (float *)nfft_malloc(numel*sizeof(float));
    assert(phi_hut_inv_host != NULL);
    
    // fill with coefficients
    for(ks=0; ks<ths->N[t]; ks++) {
      phi_hut_inv_host[ks] = 1.0/(PHI_HUT(ks-ths->N[t]/2,0));
      //printf("%g ", phi_hut_inv_host[ks]);
    }
    //printf("\n");

    // allocate and copy to device
    CUDA_SAFE_CALL( cudaMalloc( (void **)&ths->phi_hut_inv[t], 
				numel*sizeof(float)) );
    CUDA_SAFE_CALL( cudaMemcpy(ths->phi_hut_inv[t], phi_hut_inv_host, 
			       numel*sizeof(float), cudaMemcpyHostToDevice) );
    nfft_free(phi_hut_inv_host);
  } 
} /* nfft_phi_hut */

/* END: NFFT3 library source code */

// node data creation routines

// finds bin with lowest level
int nfft_node_find_bin(float x[2], int rx[2], int n[2], int numbins,
		       int2 *bincoords, int *binlevel) {
  int i, k, found, u[2], o[2], u_bin[2], o_bin[2];
  
  // get limits for current node
  for (k=0; k<2; k++) {
    u[k] = (rx[k] - THS_M - (x[k] < 0)) & (n[k] - 1);
    o[k] = (rx[k] + THS_M + (x[k] >= 0)) & (n[k] - 1);
    if (o[k] < u[k]) o[k] += n[k];
  }

  found = -1;
  for (i=0; i<numbins; i++) {
    // get limits for current bin
    u_bin[0] = (bincoords[i].x - 16) & (n[0] - 1);
    u_bin[1] = (bincoords[i].y - 16) & (n[1] - 1);
    o_bin[0] = (bincoords[i].x + 15) & (n[0] - 1);
    o_bin[1] = (bincoords[i].y + 15) & (n[1] - 1);
    for (k=0; k<2; k++)
      if (o_bin[k] < u_bin[k]) o_bin[k] += n[k];

    // check whether node goes into bin
    if ( ( ((u[0] >= u_bin[0]) && (o[0] <= o_bin[0])) ||
	   ((u[0]+n[0] >= u_bin[0]) && (o[0]+n[0] <= o_bin[0])) ) &&
	 ( ((u[1] >= u_bin[1]) && (o[1] <= o_bin[1])) ||
	   ((u[1]+n[1] >= u_bin[0]) && (o[1]+n[1] <= o_bin[1])) ) ) {
      // take bin with smallest level
      if (found < 0) 
	found = i;
      else if (binlevel[i] < binlevel[found])
	found = i;
    }
  }
  
  return(found);
}

// creates bins for each node
void nfft_node_compute_bins(nfft_plan *ths, int *p_numbins, int2 **p_bincoords,
			    int (**p_binindex)[NFFT_NODE_BIN_SIZE]) {
  int i, k;
  float *ths_x;
  float x[2];
  int rx[2];
  
  int numbins, bin;
  int2 *bincoords;
  int *binlevel;
  int (*binindex)[NFFT_NODE_BIN_SIZE];

  // copy coordinates to host memory
  ths_x = (float *)nfft_malloc(2*sizeof(float)*ths->M_total);
  CUDA_SAFE_CALL( cudaMemcpy(ths_x, ths->x, 2*ths->M_total*sizeof(float), 
			     cudaMemcpyDeviceToHost) );
  
  // alloc memory according to the worst case
  numbins = 0;
  bincoords = (int2 *)nfft_malloc(sizeof(int2)*ths->M_total);
  binlevel = (int *)nfft_malloc(sizeof(int)*ths->M_total);
  binindex = (int (*)[NFFT_NODE_BIN_SIZE])nfft_malloc(sizeof(int)*ths->M_total*NFFT_NODE_BIN_SIZE);
  
  // for each node
  for(i=0;i<ths->M_total;i++) {
    x[0] = ths_x[(i << 1)];
    x[1] = ths_x[(i << 1) + 1];

    // round coordinates
    for (k=0;k<2;k++)
      rx[k] = rintf(x[k]*(float)ths->n[k]);

    // find bin for node
    bin = nfft_node_find_bin(x, rx, ths->n, numbins, bincoords, binlevel);

    // create new bin if necessary
    if ((bin < 0) || (binlevel[bin] >= NFFT_NODE_BIN_SIZE)) {
      for(k=0;k<NFFT_NODE_BIN_SIZE;k++)
	binindex[numbins][k] = -1;
      binlevel[numbins] = 0;
      bin = numbins;

      // round coordinates to multiples of 8 
      for (k=0;k<2;k++)
	rx[k] = rintf((float)rx[k]/8)*8;
      bincoords[numbins].x = rx[0];
      bincoords[numbins].y = rx[1];

      numbins++;
    } 
    
    // insert index
    binindex[bin][binlevel[bin]] = i;
    binlevel[bin]++;
  }

  // return the significant values
  nfft_free(binlevel);
  nfft_free(ths_x);
  *p_numbins = numbins;
  *p_bincoords = bincoords;
  *p_binindex = binindex;
}

// fills one stage, returns number of elements
int nfft_node_fill_stage(nfft_plan *ths, int *map, int startindex, 
			 int numbins, int2 *bincoords,
			 int (*binindex)[NFFT_NODE_BIN_SIZE]) {
  int i, k1, k2, u[2];
  int intersect;
  int tmp;
  int stage_size;

  stage_size = 0;

  // clear map
  bzero(map, sizeof(int)*ths->n[0]*ths->n[1]);
  
  for (i=startindex; i<numbins; i++) {
    // get bin extend 
    u[0] = bincoords[i].x - 16;
    u[1] = bincoords[i].y - 16;
    
    // count number of intersections with map
    intersect = 0;
    for(k1=0; k1<32; k1++)
      for(k2=0; k2<32; k2++)
	intersect += map[(((u[0]+k1) & (ths->n[0] - 1)) 
			  << ths->log2n[1]) + ((u[1]+k2) & (ths->n[1] - 1))];
    
    // if region still free, insert bin into stage
    if (intersect == 0) {
      // occupy region
      for(k1=0; k1<32; k1++)
	for(k2=0; k2<32; k2++)
	  map[(((u[0]+k1) & (ths->n[0] - 1)) << ths->log2n[1]) 
	      + ((u[1]+k2) & (ths->n[1] - 1))] = 1;
      
      // swap bin, if neccessary
      if (i != startindex+stage_size) {
	tmp = bincoords[startindex+stage_size].x;
	bincoords[startindex+stage_size].x = bincoords[i].x;
	bincoords[i].x = tmp;
	
	tmp = bincoords[startindex+stage_size].y;
	bincoords[startindex+stage_size].y = bincoords[i].y;
	bincoords[i].y = tmp;
	
	for (k1=0; k1<NFFT_NODE_BIN_SIZE; k1++) {
	  tmp = binindex[startindex+stage_size][k1];
	  binindex[startindex+stage_size][k1] = binindex[i][k1];
	  binindex[i][k1] = tmp;
	}
      }
      
      // increase stage size
      stage_size++;
    }
  }

  return(stage_size);
}

// copies contents of bincoords/binindex to device and stores pointers
void nfft_node_copy_to_device(int numstages, int numbins, int2 *bincoords, 
			      int (*binindex)[NFFT_NODE_BIN_SIZE], 
			      int2 **dp_bincoords, 
			      int **dp_binindex) {
  int size;
  char *dp;
  
  size = numbins*(sizeof(int2) + NFFT_NODE_BIN_SIZE*sizeof(int));
  CUDA_SAFE_CALL( cudaMalloc( (void **)&dp, size) );
  
  CUDA_SAFE_CALL( cudaMemcpy(dp, bincoords, numbins*sizeof(int2), 
			     cudaMemcpyHostToDevice) );
  *dp_bincoords = (int2 *)dp;
  
  CUDA_SAFE_CALL( cudaMemcpy(dp+numbins*sizeof(int2), binindex, 
			     numbins*sizeof(int)*NFFT_NODE_BIN_SIZE, 
			     cudaMemcpyHostToDevice) );
  *dp_binindex = (int *)(dp + numbins*sizeof(int2));
}

void nfft_create_node_data(nfft_plan *ths) {
  int numbins;
  int2 *bincoords;
  int (*binindex)[NFFT_NODE_BIN_SIZE];
  
  int numstages, startindex, stage_size;
  int2 *stageinfo;
  int *map;
  
  // return if plan contains no nodes
  if (ths->x == NULL) {
    ths->M_total = 0;
    return;
  }
  if (ths->M_total <= 0)
    return;
  
  nfft_node_compute_bins(ths, &numbins, &bincoords, &binindex);
  
  numstages = 0;
  stageinfo = (int2 *)nfft_malloc(sizeof(int2)*ths->M_total);
  map = (int *)nfft_malloc(sizeof(int)*ths->n[0]*ths->n[1]);
  
  startindex = 0;
  // fill stage
  while (stage_size = nfft_node_fill_stage(ths, map, startindex, 
					   numbins, bincoords, binindex)) {
    // write info about stage
    stageinfo[numstages].y = startindex;
    stageinfo[numstages].x = stage_size;
    startindex += stage_size;
    numstages++;
  }
  
  // stage info and bin output
  if (0) {
    for(int i = 0; i < numstages; i++) {
      printf("Stage %d: start %d, length %d\n", i, stageinfo[i].y, 
	     stageinfo[i].x);
    }
    for(int i = 0; i < numbins; i++) {
      printf("Bin %d at (%d,%d): Nodes ", i, bincoords[i].x, bincoords[i].y);
      for(int j = 0; j < NFFT_NODE_BIN_SIZE; j++) 
	printf("%d ", binindex[i][j]);
      printf("\n");
    }
  }
  
  // write info to plan
  nfft_free(map);
  stageinfo = (int2 *)realloc(stageinfo, sizeof(int2)*numstages);
  ths->stageinfo = stageinfo;
  ths->stages = numstages;
  // copy bin info to device
  ths->bins = numbins;
  nfft_node_copy_to_device(numstages, numbins, bincoords, binindex,
			   &ths->bincoords, &ths->binindex);
}

void nfft_init_2d(nfft_plan *ths, int N1, int N2, int M_total, int coils,
		  size_t coil_stride, float_complex *f_hat, float_complex *f, 
		  float_complex *g, float *x) {
  //cufftResult cufftRes;

  // cannot init a null-pointer plan
  if (!ths) return;
  
  // reset memory allocation flags
  ths->alloc_flags = 0;

  // set dimensions
  ths->N[0] = N1; ths->N[1] = N2;
  ths->n[0] = nfft_next_power_of_2(2*N1, &ths->log2n[0]);
  ths->n[1] = nfft_next_power_of_2(2*N2, &ths->log2n[1]);
  ths->M_total = M_total;
  ths->coils = coils;
  ths->coil_stride = coil_stride;
  
  // set products
  ths->N_total = N1*N2;
  ths->n_total = ths->n[0]*ths->n[1];
  
  // set oversampling factors
  ths->sigma[0] = ((float)ths->n[0])/ths->N[0];
  ths->sigma[1] = ((float)ths->n[1])/ths->N[1];
  ths->b[0] = THS_B(0);
  ths->b[1] = THS_B(1);
  
  // set pointers
  if (f_hat) 
    ths->f_hat = f_hat;
  else {
    CUDA_SAFE_CALL( cudaMalloc( (void **)&ths->f_hat, 
				coils*ths->N_total*sizeof(float_complex)) );
    ths->alloc_flags |= NFFT_INIT_ALLOC_F_HAT;
  }

  if (f) 
    ths->f = f;
  else {
    CUDA_SAFE_CALL( cudaMalloc( (void **)&ths->f, 
				coils*ths->M_total*sizeof(float_complex)) );
    ths->coil_stride = sizeof(float_complex)*ths->M_total;
    ths->alloc_flags |= NFFT_INIT_ALLOC_F;
  }

  if (x)
    ths->x = x;
  else {
    ths->x = 0;
    ths->M_total = 0;
  }
  
  // create bin coordinate/index pointers and stages
  nfft_create_node_data(ths);
  ths->alloc_flags |= NFFT_INIT_ALLOC_BIN;
  
  // allocate memory for g if necessary
  if (g) 
    ths->g = g;
  else {
    CUDA_SAFE_CALL( cudaMalloc( (void **)&ths->g, 
				coils*ths->n_total*sizeof(float_complex)) );
    ths->alloc_flags |= NFFT_INIT_ALLOC_G;
  }
  
  // get CUFFT plan
  CUFFT_SAFE_CALL( cufftPlan2d(&ths->cufft_plan, ths->n[0], 
			       ths->n[1], CUFFT_C2C) );
  ths->alloc_flags |= NFFT_INIT_ALLOC_CUFFT;

  // precompute phi_hut
  nfft_precompute_phi_hut(ths);
}

void nfft_finalize_2d(nfft_plan *ths) {
  if (!ths) return;

  if (ths->alloc_flags & NFFT_INIT_ALLOC_CUFFT) 
    CUFFT_SAFE_CALL( cufftDestroy(ths->cufft_plan) );

  if (ths->alloc_flags & NFFT_INIT_ALLOC_BIN) {
    CUDA_SAFE_CALL( cudaFree(ths->bincoords) );
    nfft_free(ths->stageinfo);
  }
  
  if (ths->alloc_flags & NFFT_INIT_ALLOC_G) 
    CUDA_SAFE_CALL( cudaFree(ths->g) );
  
  if (ths->alloc_flags & NFFT_INIT_ALLOC_F_HAT) 
    CUDA_SAFE_CALL( cudaFree(ths->f_hat) );
  
  if (ths->alloc_flags & NFFT_INIT_ALLOC_F) 
    CUDA_SAFE_CALL( cudaFree(ths->f) );
}

///////////////////////////
// computation for nfft

/*
 * Nfft helper routines
 * 
 * Base on NFFT3 library source code 
 * Copyright (c) 2002, 2009 Jens Keiner, Stefan Kunis, Daniel Potts
 */

/** computes 2m+2 indices for the matrix B
 */
__device__ inline void nfft_uo(float xj, int n, int *up, int *op)
{
  int c = rintf(xj*n);
  
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
__device__ inline void nfft_uo2(float x, int n, int *u, int *o)
{
  int c = rintf(x*n);
  
  if(x < 0)
    {
      *u=(c-1-THS_M+n)%n;
      *o=(c+  THS_M+n)%n;
    }
  else
    {
      *u=(c  -THS_M+n)%n;
      *o=(c+1+THS_M+n)%n;
    }
}

/* END: NFFT3 library source code */


__device__ void setScaleComplex(float_complex g, float_complex f, 
				       float fac) {
  //printf("g= %012X f=%012X fac=%g\n", g, f, fac);
  g[0] = f[0]*fac; g[1] = f[1]*fac;
}

#define F_HAT_IDX(k0,k1) f_hat[(coil*N[0] + (k0))*N[1] + (k1)]  
#define G_IDX_A(k0,k1) g[(((coil << log2n0) + ((((int)(k0)) & (n[0] - 1)))) \
			   << log2n1) + (((int)(k1)) & (n[1] - 1))]

__global__ void nfft_trafo_2d_A(int N0, int N1, int n0, int n1,
				int log2n0, int log2n1,
				float *phi_hut_inv_0,
				float *phi_hut_inv_1,
				float_complex *f_hat, 
				float_complex *g,
				int coil) 
{
  __shared__ int N[2], n[2]; 
  __shared__ float phi_hut_inv1[2][THS_BLOCK_LEN], 
    phi_hut_inv2[2][THS_BLOCK_LEN];
  __shared__ float_complex f_hat_local[2][2][THS_BLOCK_LEN][THS_BLOCK_LEN];
  int k0, k1, tid, l0, l1;
  
  // set up indices
  tid = threadIdx.y*THS_BLOCK_LEN + threadIdx.x;
  k0 = blockIdx.y*blockDim.y + threadIdx.y;
  k1 = blockIdx.x*blockDim.x + threadIdx.x;
  l0 = threadIdx.y; 
  l1 = threadIdx.x;

  // copy structure data into shared memory
  if (tid == 0) {
    N[0] = N0; N[1] = N1;
    n[0] = n0; n[1] = n1;
  }
  __syncthreads();
  
  // copy phi_hut_inv into shared memory
  if (l0 == 0) {
    phi_hut_inv1[0][l1] = phi_hut_inv_0[k0 + l1];
    phi_hut_inv2[0][l1] = phi_hut_inv_0[k0 + l1 + (N[0] >> 1)];
  }
  if (l0 == 2) {
    phi_hut_inv1[1][l1] = phi_hut_inv_1[k1];
    phi_hut_inv2[1][l1] = phi_hut_inv_1[k1 + (N[1] >> 1)];
  }
  __syncthreads();
  
  // read f_hat
  if ((k0 < (N[0] >> 1)) && (k1 < (N[1] >> 1))) {
    *((float2 *)&f_hat_local[0][0][l0][l1]) 
      = *((float2 *)&F_HAT_IDX(k0, k1));
    *((float2 *)&f_hat_local[0][1][l0][l1]) 
      = *((float2 *)&F_HAT_IDX(k0, k1+(N[1] >> 1)));
    *((float2 *)&f_hat_local[1][0][l0][l1]) 
      = *((float2 *)&F_HAT_IDX(k0+(N[0] >> 1), k1));
    *((float2 *)&f_hat_local[1][1][l0][l1]) 
      = *((float2 *)&F_HAT_IDX(k0+(N[0] >> 1), k1+(N[1] >> 1)));
  }
  
  // compute products
  f_hat_local[0][0][l0][l1][0] *= phi_hut_inv1[0][l0]*phi_hut_inv1[1][l1];
  f_hat_local[0][0][l0][l1][1] *= phi_hut_inv1[0][l0]*phi_hut_inv1[1][l1];
  f_hat_local[0][1][l0][l1][0] *= phi_hut_inv1[0][l0]*phi_hut_inv2[1][l1];
  f_hat_local[0][1][l0][l1][1] *= phi_hut_inv1[0][l0]*phi_hut_inv2[1][l1];
  f_hat_local[1][0][l0][l1][0] *= phi_hut_inv2[0][l0]*phi_hut_inv1[1][l1];
  f_hat_local[1][0][l0][l1][1] *= phi_hut_inv2[0][l0]*phi_hut_inv1[1][l1];
  f_hat_local[1][1][l0][l1][0] *= phi_hut_inv2[0][l0]*phi_hut_inv2[1][l1];
  f_hat_local[1][1][l0][l1][1] *= phi_hut_inv2[0][l0]*phi_hut_inv2[1][l1];
  
  // write in g
  if ((k0 < (N[0] >> 1)) && (k1 < (N[1] >> 1))) {
    *((float2 *)&G_IDX_A(k0-(N[0] >> 1), k1-(N[1] >> 1))) 
      = *((float2 *)&f_hat_local[0][0][l0][l1]);
    *((float2 *)&G_IDX_A(k0-(N[0] >> 1), k1)) 
      = *((float2 *)&f_hat_local[0][1][l0][l1]);
    *((float2 *)&G_IDX_A(k0, k1-(N[1] >> 1))) 
      = *((float2 *)&f_hat_local[1][0][l0][l1]);
    *((float2 *)&G_IDX_A(k0, k1)) 
      = *((float2 *)&f_hat_local[1][1][l0][l1]);
  }
}

__global__ void nfft_trafo_2d_B(int n0, int n1, int log2n0, int log2n1,
				int M_total_src, float ths_b0, float ths_b1, 
				float *x, float_complex *g, 
				float_complex *f, size_t coil_stride)
{
  __shared__ int log2n[2];
  __shared__ int n[2];
  __shared__ float ths_b[2];
  
  __shared__ float phiarg[2*THS_BLOCK_LEN];
  __shared__ float val[2*THS_BLOCK_LEN];
  __shared__ float sqrtval[2*THS_BLOCK_LEN];
  __shared__ float psij_const[2][THS_BLOCK_LEN];
  __shared__ float_complex acc[THS_WINDOW_SPACE];
  __shared__ float xj[2];
  __shared__ int u[2], o[2];
  int i, j, coil, l0, l1, tid;
  int gidx;
  
  // block index corresponds to node/coil
  j = blockIdx.y;
  coil = blockIdx.x;
  // thread ids correspond coordinates in window
  l0 = threadIdx.x; l1 = threadIdx.y;
  tid = blockDim.x*l1 + l0;
  
  // store parameters in shared memory
  if (tid == 0) {
    n[0] = n0; log2n[0] = log2n0; 
    n[1] = n1; log2n[1] = log2n1;
    ths_b[0] = ths_b0; ths_b[1] = ths_b1;
    //x = x_src; g = g_src; f = f_src;
  }
  __syncthreads();

  // load node point and lower/upper bounds
  if (tid < 2) {
    xj[tid] = x[(j << 1) + tid];
    nfft_uo(xj[tid], n[tid], &u[tid], &o[tid]);
  } 
  __syncthreads();

  // compute psij
  if (l1 < 2) {
    // argument for phi
    phiarg[tid] = xj[l1] - ((float)(u[l1]+l0))/n[l1];
    
    // case 1
    val[tid] = KB_M_SQR - sqr(phiarg[tid]*((float)n[l1]));
    sqrtval[tid] = sqrtf(fabs(val[tid]));
    if (val[tid] > 0) {
      psij_const[l1][l0] = sinhf(ths_b[l1]*sqrtval[tid])/
	(CUDART_PI_F*sqrtval[tid]); 
    } else {
      // case 2
      if (val[tid] < 0) {
	psij_const[l1][l0] = __sinf(ths_b[l1]*sqrtval[tid])/
	  (CUDART_PI_F*sqrtval[tid]); 
      } else {
	// case 3
	psij_const[l1][l0] = 1.0f;
      }
    }
  }
  __syncthreads();
  
  // get indices
  gidx = (((coil << log2n[0]) + ((u[0] + l0) & (n[0]-1))) 
	  << log2n[1]) + ((u[1] + l1) & (n[1]-1));
  // write scaled product into accumulator
  if ((l0 < THS_WINDOW_LEN) && (l1 < THS_WINDOW_LEN)) {
    acc[tid][0] = g[gidx][0]*psij_const[0][l0]*psij_const[1][l1];
    acc[tid][1] = g[gidx][1]*psij_const[0][l0]*psij_const[1][l1];
  } else {
    acc[tid][0] = acc[tid][1] = 0;
  }
  
  // sum up
  for (i = THS_WINDOW_SPACE >> 1; i > 0; i >>= 1)
    {
      __syncthreads();
      if (tid < i) {
	acc[tid][0] += acc[tid + i][0];
	acc[tid][1] += acc[tid + i][1];
      }
    }
  __syncthreads();
  
  // write out
  if (tid < 2) {
    ((float_complex *)(((char *)f) + coil*coil_stride))[j][tid] = acc[0][tid];
  }
}

#define G_IDX(l0,l1) g[(((coil << log2n0) + ((u_local[0] + l0) & (n[0] - 1))) \
			<< log2n1) + ((u_local[1] + l1) & (n[1] - 1))]  

__global__ void nfft_trafo_2d_B2(int n0, int n1, int log2n0, int log2n1,
				 float ths_b0, float ths_b1, 
				 float *x, float_complex *g, 
				 float_complex *f, size_t coil_stride,
				 int *bincoords, int *binindex)
{
  __shared__ int n[2];
  __shared__ float ths_b[2];
  
  __shared__ float phiarg[2*THS_BLOCK_LEN];
  __shared__ float val[2*THS_BLOCK_LEN];
  __shared__ float sqrtval[2*THS_BLOCK_LEN];
  __shared__ float psij_const[2][THS_BLOCK_LEN];

  __shared__ float_complex g_local[32][32];
  __shared__ int u_local[2];

  __shared__ float_complex acc[THS_WINDOW_SPACE];
  __shared__ float xj[2];
  __shared__ int u[2], o[2];

  __shared__ int bidx, binidx;
  
  int i, j, coil, l0, l1, tid;
  //int bidx, binidx;
  
  // block index corresponds to coil/node in stage
  j = blockIdx.y;
  coil = blockIdx.x;
  // thread ids correspond coordinates in window
  l0 = threadIdx.y; l1 = threadIdx.x;
  tid = blockDim.x*l0 + l1;

  // store parameters in shared memory
  if (tid == 0) {
    n[0] = n0; n[1] = n1;
    ths_b[0] = ths_b0; ths_b[1] = ths_b1;
  }
  __syncthreads();

  // copy lower index into u_local
  if (tid < 2) {
    u_local[tid] = bincoords[(j << 1) + tid] - 16;
  }
  __syncthreads();
  
  // copy g into shared memory
  *((float *)&g_local[l0][0] + l1) = *((float *)&G_IDX(l0,0) + l1);
  *((float *)&g_local[l0][8] + l1) = *((float *)&G_IDX(l0,8) + l1);
  *((float *)&g_local[l0][16] + l1) = *((float *)&G_IDX(l0,16) + l1);
  *((float *)&g_local[l0][24] + l1) = *((float *)&G_IDX(l0,24) + l1);
  *((float *)&g_local[l0+16][0] + l1) = *((float *)&G_IDX(l0+16,0) + l1);
  *((float *)&g_local[l0+16][8] + l1) = *((float *)&G_IDX(l0+16,8) + l1);
  *((float *)&g_local[l0+16][16] + l1) = *((float *)&G_IDX(l0+16,16) + l1);
  *((float *)&g_local[l0+16][24] + l1) = *((float *)&G_IDX(l0+16,24) + l1);
  __syncthreads();
  
  // initialize bidx
  if (tid == 0) {
    bidx = 0;
    binidx = binindex[(j << 5) + bidx];
  }
  __syncthreads();
  
  while ((bidx < NFFT_NODE_BIN_SIZE) && (binidx >= 0)) {
    // load node point and lower/upper bounds
    if (tid < 2) {
      xj[tid] = x[(binidx << 1) + tid];
      nfft_uo(xj[tid], n[tid], &u[tid], &o[tid]);
    }
    __syncthreads();
    
    // compute psij
    if (l0 < 2) {
      // argument for phi
      phiarg[tid] = xj[l0] - ((float)(u[l0]+l1))/n[l0];
      
      // case 1
      val[tid] = KB_M_SQR - sqr(phiarg[tid]*((float)n[l0]));
      sqrtval[tid] = sqrtf(fabs(val[tid]));
      if (val[tid] > 0) {
  	psij_const[l0][l1] = sinhf(ths_b[l0]*sqrtval[tid])/
  	  (CUDART_PI_F*sqrtval[tid]);
      } else {
  	// case 2
  	if (val[tid] < 0) {
  	  psij_const[l0][l1] = __sinf(ths_b[l0]*sqrtval[tid])/
  	    (CUDART_PI_F*sqrtval[tid]);
  	} else {
  	  // case 3
  	  psij_const[l0][l1] = 1.0f;
  	}
      }
    }
    __syncthreads();
    
    // write scaled product into accumulator
    if ((l0 < THS_WINDOW_LEN) && (l1 < THS_WINDOW_LEN)) {
      acc[tid][0] = g_local[ (l0 + u[0] - u_local[0]) & 31 ]
  	[ (l1 + u[1] - u_local[1]) & 31 ][0]
  	*psij_const[0][l0]*psij_const[1][l1];
      acc[tid][1] = g_local[ (l0 + u[0] - u_local[0]) & 31 ]
  	[ (l1 + u[1] - u_local[1]) & 31 ][1]
  	*psij_const[0][l0]*psij_const[1][l1];
    } else {
      acc[tid][0] = acc[tid][1] = 0;
    }

    // sum up
    for (i = THS_WINDOW_SPACE >> 1; i > 0; i >>= 1)
      {
  	__syncthreads();
  	if (tid < i) {
  	  acc[tid][0] += acc[tid + i][0];
  	  acc[tid][1] += acc[tid + i][1];
  	}
      }
    
    // write out
    __syncthreads();
    if (tid < 2) {
      ((float_complex *)(((char *)f) + coil*coil_stride))[binidx][tid] 
	= acc[0][tid];
    }

    //if (tid == 0)
    //  printf("%d:(%g, %g) ", binidx, acc[0][0], acc[0][1]);
    
    // increment counter
    if (tid == 0) {
      bidx++;
      binidx = binindex[(j << 5) + bidx];
    }
    __syncthreads();
  }
}

void nfft_trafo_2d(nfft_plan *ths)
{
  dim3 dimBlock, dimGrid;
  //cufftResult cufftRes;
  int i;

  // clear g
  CUDA_SAFE_CALL( cudaMemset(ths->g, 0, ths->n_total*ths->coils
			     *sizeof(float_complex)) );
  
  // multiply with phi_hut_inv and fftshift
  dimBlock.x = THS_BLOCK_LEN;
  dimBlock.y = THS_BLOCK_LEN;
  dimBlock.z = 1;
  dimGrid.x = (ths->N[1]/2 + THS_BLOCK_LEN - 1)/dimBlock.x;
  dimGrid.y = (ths->N[0]/2 + THS_BLOCK_LEN - 1)/dimBlock.y;
  dimGrid.z = 1;
  //printf("starting part A...\n");
  for (i=0; i<ths->coils; i++) {
    nfft_trafo_2d_A<<<dimGrid, dimBlock>>>(ths->N[0], ths->N[1], ths->n[0],
					   ths->n[1], ths->log2n[0], 
					   ths->log2n[1],
					   ths->phi_hut_inv[0],
					   ths->phi_hut_inv[1], ths->f_hat,
					   ths->g, i);
    CUT_CHECK_ERROR("nfft_trafo_2d_A() failed.");
  }
  
  // perform fft on g
  //printf("performing FFT...\n");
  for (i=0; i<ths->coils; i++) {
    CUFFT_SAFE_CALL( cufftExecC2C(ths->cufft_plan, 
				  (cufftComplex *)(ths->g+i*ths->n[0]*ths->n[1]), 
				  (cufftComplex *)(ths->g+i*ths->n[0]*ths->n[1]), 
				  CUFFT_FORWARD) );
  }
 
  // do convolution steps
  for(i=0; i < ths->stages; i++) {
    dimBlock.x = THS_BLOCK_LEN;
    dimBlock.y = THS_BLOCK_LEN;
    dimBlock.z = 1;
    dimGrid.x = ths->coils;
    dimGrid.y = ths->stageinfo[i].x;
    dimGrid.z = 1;
    
    nfft_trafo_2d_B2<<<dimGrid, dimBlock>>>(ths->n[0], ths->n[1], 
					    ths->log2n[0], ths->log2n[1],
					    ths->b[0], ths->b[1],
					    ths->x, ths->g, ths->f, 
					    ths->coil_stride, 
					    (int *)(ths->bincoords 
						    + ths->stageinfo[i].y),
					    ths->binindex + 
					    NFFT_NODE_BIN_SIZE
					    *ths->stageinfo[i].y);
    CUT_CHECK_ERROR("nfft_trafo_2d_B() failed.");
  }
}

// adjoint transformation

__device__ void setScaleComplexAd(float_complex g, float_complex f, 
				  float fac) {
  //printf("g= %012X f=%012X fac=%g\n", g, f, fac);
  f[0] = g[0]*fac; f[1] = g[1]*fac;
}

__global__ void nfft_adjoint_2d_A(int N0, int N1, int n0, int n1,
				  int log2n0, int log2n1,
				  float *phi_hut_inv_0,
				  float *phi_hut_inv_1,
				  float_complex *f_hat, 
				  float_complex *g,
				  int coil) 
{
  __shared__ int N[2], n[2]; 
  __shared__ float phi_hut_inv1[2][THS_BLOCK_LEN], 
    phi_hut_inv2[2][THS_BLOCK_LEN];
  __shared__ float_complex g_local[2][2][THS_BLOCK_LEN][THS_BLOCK_LEN];
  int k0, k1, tid, l0, l1;
  
  // set up indices
  tid = threadIdx.y*THS_BLOCK_LEN + threadIdx.x;
  k0 = blockIdx.y*blockDim.y + threadIdx.y;
  k1 = blockIdx.x*blockDim.x + threadIdx.x;
  l0 = threadIdx.y; 
  l1 = threadIdx.x;

  // copy structure data into shared memory
  if (tid == 0) {
    N[0] = N0; N[1] = N1;
    n[0] = n0; n[1] = n1;
  }
  __syncthreads();
  
  // copy phi_hut_inv into shared memory
  if (l0 == 0) {
    phi_hut_inv1[0][l1] = phi_hut_inv_0[k0 + l1];
    phi_hut_inv2[0][l1] = phi_hut_inv_0[k0 + l1 + (N[0] >> 1)];
  }
  if (l0 == 2) {
    phi_hut_inv1[1][l1] = phi_hut_inv_1[k1];
    phi_hut_inv2[1][l1] = phi_hut_inv_1[k1 + (N[1] >> 1)];
  }
  __syncthreads();

  // read from g
  if ((k0 < (N[0] >> 1)) && (k1 < (N[1] >> 1))) {
    *((float2 *)&(g_local[0][0][l0][l1]))
      = *((float2 *)&G_IDX_A(k0-(N[0] >> 1), k1-(N[1] >> 1)));
    *((float2 *)&g_local[0][1][l0][l1])
      = *((float2 *)&G_IDX_A(k0-(N[0] >> 1), k1));
    *((float2 *)&g_local[1][0][l0][l1])
      = *((float2 *)&G_IDX_A(k0, k1-(N[1] >> 1)));
    *((float2 *)&g_local[1][1][l0][l1])
      = *((float2 *)&G_IDX_A(k0, k1));
  }
   
  // compute products
  g_local[0][0][l0][l1][0] *= phi_hut_inv1[0][l0]*phi_hut_inv1[1][l1]; 
  g_local[0][0][l0][l1][1] *= phi_hut_inv1[0][l0]*phi_hut_inv1[1][l1]; 
  g_local[0][1][l0][l1][0] *= phi_hut_inv1[0][l0]*phi_hut_inv2[1][l1]; 
  g_local[0][1][l0][l1][1] *= phi_hut_inv1[0][l0]*phi_hut_inv2[1][l1]; 
  g_local[1][0][l0][l1][0] *= phi_hut_inv2[0][l0]*phi_hut_inv1[1][l1]; 
  g_local[1][0][l0][l1][1] *= phi_hut_inv2[0][l0]*phi_hut_inv1[1][l1]; 
  g_local[1][1][l0][l1][0] *= phi_hut_inv2[0][l0]*phi_hut_inv2[1][l1]; 
  g_local[1][1][l0][l1][1] *= phi_hut_inv2[0][l0]*phi_hut_inv2[1][l1]; 
  
  // write to f_hat
  if ((k0 < (N[0] >> 1)) && (k1 < (N[1] >> 1))) {
    *((float2 *)&F_HAT_IDX(k0, k1))
      = *((float2 *)&g_local[0][0][l0][l1]);
    *((float2 *)&F_HAT_IDX(k0, k1+(N[1] >> 1)))
      = *((float2 *)&g_local[0][1][l0][l1]);
    *((float2 *)&F_HAT_IDX(k0+(N[0] >> 1), k1))
      = *((float2 *)&g_local[1][0][l0][l1]);
    *((float2 *)&F_HAT_IDX(k0+(N[0] >> 1), k1+(N[1] >> 1))) 
      = *((float2 *)&g_local[1][1][l0][l1]);
  }
}

__device__ void atomicAddFloat(float *g, float f) {
  while (f)
    f = atomicExch(g, f + atomicExch(g, 0.0f));
}
  
__global__ void nfft_adjoint_2d_B(int n0, int n1, int log2n0, int log2n1,
				  int M_total_src, float ths_b0, float ths_b1, 
				  float *x, float_complex *g, 
				  float_complex *f, size_t coil_stride)
{
  __shared__ int n[2];
  __shared__ int log2n[2];
  __shared__ float ths_b[2];
  //__shared__ float *x;
  //__shared__ float_complex *g, *f;

  __shared__ float phiarg[2*THS_BLOCK_LEN], val[2*THS_BLOCK_LEN],
    sqrtval[2*THS_BLOCK_LEN];
  __shared__ float psij_const[2][THS_BLOCK_LEN];
  __shared__ float fj[2];
  __shared__ float xj[2];
  __shared__ int u[2], o[2];

  int j, coil, l0, l1, tid;
  int gidx;
  
  // block index corresponds to node/coil
  j = blockIdx.y;
  coil = blockIdx.x;
  // thread ids correspond coordinates in window
  l0 = threadIdx.x; l1 = threadIdx.y;
  tid = blockDim.x*l1 + l0;
  
  // store parameters in shared memory
  if (tid == 0) {
    n[0] = n0; n[1] = n1;
    log2n[0] = log2n0; log2n[1] = log2n1;
    ths_b[0] = ths_b0; ths_b[1] = ths_b1;
    //x = x_src; g = g_src; f = f_src;
  }
  __syncthreads();
  
  // load value at current node
  if ((tid >= 128) && (tid < 130)) {
    fj[tid-128] = ((float_complex *)
		   (((char *)f) + coil*coil_stride))[j][tid-128];
  }
  
  // load node point and lower/upper bounds
  if (tid < 2) {
    xj[tid] = x[(j << 1) + tid];
    nfft_uo(xj[tid], n[tid], &u[tid], &o[tid]);
  } 
  __syncthreads();

  // compute psij
  if (l1 < 2) {
    // argument for phi
    phiarg[tid] = xj[l1] - ((float)(u[l1]+l0))/n[l1];
    
    // case 1
    val[tid] = KB_M_SQR - sqr(phiarg[tid]*((float)n[l1]));
    sqrtval[tid] = sqrtf(fabs(val[tid]));
    if (val[tid] > 0) {
      psij_const[l1][l0] = sinhf(ths_b[l1]*sqrtval[tid])/
	(CUDART_PI_F*sqrtval[tid]); 
    } else {
      // case 2
      if (val[tid] < 0) {
	psij_const[l1][l0] = __sinf(ths_b[l1]*sqrtval[tid])/
	  (CUDART_PI_F*sqrtval[tid]); 
      } else {
	// case 3
	psij_const[l1][l0] = 1.0f;
      }
    }
  }
  __syncthreads();
  
  // get indices
  gidx = (((coil << log2n[0]) + ((u[0] + l0) & (n[0]-1))) 
	  << log2n[1]) + ((u[1] + l1) & (n[1]-1));
  // write scaled product into right coefficient of g
  if ((l0 < THS_WINDOW_LEN) && (l1 < THS_WINDOW_LEN)) {
    atomicAddFloat(&g[gidx][0], fj[0]*psij_const[0][l0]*psij_const[1][l1]);
    atomicAddFloat(&g[gidx][1], fj[1]*psij_const[0][l0]*psij_const[1][l1]);
  }
}


__global__ void nfft_adjoint_2d_B2(int n0, int n1, int log2n0, int log2n1,
				   float ths_b0, float ths_b1, 
				   float *x, float_complex *g, 
				   float_complex *f, size_t coil_stride,
				   int *bincoords, int *binindex)
{
  __shared__ int n[2];
  __shared__ float ths_b[2];
  
  __shared__ float phiarg[2*THS_BLOCK_LEN];
  __shared__ float val[2*THS_BLOCK_LEN];
  __shared__ float sqrtval[2*THS_BLOCK_LEN];
  __shared__ float psij_const[2][THS_BLOCK_LEN];

  __shared__ float_complex g_local[32][32];
  __shared__ int u_local[2];

  __shared__ float xj[2];
  __shared__ float fj[2];
  __shared__ int u[2], o[2];

  __shared__ int bidx, binidx;
  
  int j, coil, l0, l1, tid;
  //int bidx, binidx;
  
  // block index corresponds to coil/node in stage
  j = blockIdx.y;
  coil = blockIdx.x;
  // thread ids correspond coordinates in window
  l0 = threadIdx.y; l1 = threadIdx.x;
  tid = blockDim.x*l0 + l1;

  // store parameters in shared memory
  if (tid == 0) {
    n[0] = n0; n[1] = n1;
    ths_b[0] = ths_b0; ths_b[1] = ths_b1;
  }
  __syncthreads();

  // copy lower index into u_local
  if (tid < 2) {
    u_local[tid] = bincoords[(j << 1) + tid] - 16;
  }
  __syncthreads();

  // clear g_local
  g_local[l0][l1][0]       = g_local[l0][l1][1]       = 0.0f; 
  g_local[l0+16][l1][0]    = g_local[l0+16][l1][1]    = 0.0f;
  g_local[l0][l1+16][0]    = g_local[l0][l1+16][1]    = 0.0f;
  g_local[l0+16][l1+16][0] = g_local[l0+16][l1+16][1] = 0.0f;
  __syncthreads();
  
  // initialize bidx
  if (tid == 0) {
    bidx = 0;
    binidx = binindex[(j << 5) + bidx];
  }
  __syncthreads();
  
  while ((bidx < NFFT_NODE_BIN_SIZE) && (binidx >= 0)) {
    // load node point and lower/upper bounds
    if (tid < 2) {
      xj[tid] = x[(binidx << 1) + tid];
      nfft_uo(xj[tid], n[tid], &u[tid], &o[tid]);
    }

    // load value at current node
    if ((tid >= 128) && (tid < 130)) {
      fj[tid-128] = ((float_complex *)(((char *)f) + coil*coil_stride))
	[binidx][tid-128];
    }
    __syncthreads();
    
    // compute psij
    if (l0 < 2) {
      // argument for phi
      phiarg[tid] = xj[l0] - ((float)(u[l0]+l1))/n[l0];
      
      // case 1
      val[tid] = KB_M_SQR - sqr(phiarg[tid]*((float)n[l0]));
      sqrtval[tid] = sqrtf(fabs(val[tid]));
      if (val[tid] > 0) {
  	psij_const[l0][l1] = sinhf(ths_b[l0]*sqrtval[tid])/
  	  (CUDART_PI_F*sqrtval[tid]);
      } else {
  	// case 2
  	if (val[tid] < 0) {
  	  psij_const[l0][l1] = __sinf(ths_b[l0]*sqrtval[tid])/
  	    (CUDART_PI_F*sqrtval[tid]);
  	} else {
  	  // case 3
  	  psij_const[l0][l1] = 1.0f;
  	}
      }
    }
    __syncthreads();
    
    // write scaled fj into local g cache
    if ((l0 < THS_WINDOW_LEN) && (l1 < THS_WINDOW_LEN)) {
      g_local[ (l0 + u[0] - u_local[0]) & 31 ]
	[ (l1 + u[1] - u_local[1]) & 31 ][0] 
	+= fj[0]*psij_const[0][l0]*psij_const[1][l1];
      g_local[ (l0 + u[0] - u_local[0]) & 31 ]
	[ (l1 + u[1] - u_local[1]) & 31 ][1] 
	+= fj[1]*psij_const[0][l0]*psij_const[1][l1];
    }
  
    // increment counter
    if (tid == 0) {
      bidx++;
      binidx = binindex[(j << 5) + bidx];
    }
    __syncthreads();
  }

  // add g from shared memory
  *((float *)&G_IDX(l0,0) + l1) += *((float *)&g_local[l0][0] + l1);
  *((float *)&G_IDX(l0,8) + l1) += *((float *)&g_local[l0][8] + l1);
  *((float *)&G_IDX(l0,16) + l1) += *((float *)&g_local[l0][16] + l1);
  *((float *)&G_IDX(l0,24) + l1) += *((float *)&g_local[l0][24] + l1);
  *((float *)&G_IDX(l0+16,0) + l1) += *((float *)&g_local[l0+16][0] + l1);
  *((float *)&G_IDX(l0+16,8) + l1) += *((float *)&g_local[l0+16][8] + l1);
  *((float *)&G_IDX(l0+16,16) + l1) += *((float *)&g_local[l0+16][16] + l1);
  *((float *)&G_IDX(l0+16,24) + l1) += *((float *)&g_local[l0+16][24] + l1);
}


void nfft_adjoint_2d(nfft_plan *ths)
{
  dim3 dimBlock, dimGrid;
  //cufftResult cufftRes;
  int i;

  // clear g
  CUDA_SAFE_CALL( cudaMemset(ths->g, 0, ths->n_total*ths->coils*
			     sizeof(float_complex)) );
  
  // do adjoint convolution steps
  for(i=0; i < ths->stages; i++) {
    dimBlock.x = THS_BLOCK_LEN;
    dimBlock.y = THS_BLOCK_LEN;
    dimBlock.z = 1;
    dimGrid.x = ths->coils;
    dimGrid.y = ths->stageinfo[i].x;
    dimGrid.z = 1;
    //printf("starting part B...\n");
    
    //printf("stage: %d (%d,%d)\n", i, ths->stageinfo[i].x, ths->stageinfo[i].y);
    
    nfft_adjoint_2d_B2<<<dimGrid, dimBlock>>>(ths->n[0], ths->n[1], 
					      ths->log2n[0], ths->log2n[1],
					      ths->b[0], ths->b[1],
					      ths->x, ths->g, ths->f, 
					      ths->coil_stride, 
					      (int *)(ths->bincoords 
						      + ths->stageinfo[i].y),
					      ths->binindex + 
					      NFFT_NODE_BIN_SIZE
					      *ths->stageinfo[i].y);
    CUT_CHECK_ERROR("nfft_adjoint_2d_B() failed.");
  }

  // perform adjoint fft on g
  //printf("performing FFT...\n");
  for (i=0; i<ths->coils; i++) {
    CUFFT_SAFE_CALL( cufftExecC2C(ths->cufft_plan, 
			    (cufftComplex *)(ths->g+i*ths->n[0]*ths->n[1]), 
			    (cufftComplex *)(ths->g+i*ths->n[0]*ths->n[1]), 
					     CUFFT_INVERSE) );
    //assert(cufftRes == CUFFT_SUCCESS);
  }
  
  // multiply with phi_hut_inv and fftshift
  dimBlock.x = THS_BLOCK_LEN;
  dimBlock.y = THS_BLOCK_LEN;
  dimBlock.z = 1;
  dimGrid.x = (ths->N[1]/2 + THS_BLOCK_LEN - 1)/dimBlock.x;
  dimGrid.y = (ths->N[0]/2 + THS_BLOCK_LEN - 1)/dimBlock.y;
  dimGrid.z = 1;
  //printf("starting part A...\n");
  
  for (i=0; i<ths->coils; i++) {
    nfft_adjoint_2d_A<<<dimGrid, dimBlock>>>(ths->N[0], ths->N[1], 
					     ths->n[0], ths->n[1], 
					     ths->log2n[0], ths->log2n[1],
					     ths->phi_hut_inv[0],
					     ths->phi_hut_inv[1], ths->f_hat,
					     ths->g, i);
    CUT_CHECK_ERROR("nfft_adjoint_2d_A() failed.");
  }
}
