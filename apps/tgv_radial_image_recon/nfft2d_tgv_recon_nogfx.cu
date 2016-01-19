#include <cuda_runtime_api.h>
#include <stdio.h>
#include <time.h>
#include <poll.h>
#include <nfft2d.h>
#include <ParseParam.h>

// problem data
#define PROBLEM_N1 512
#define PROBLEM_N2 512
#define READOUT_LEN 512
#define NUM_SPOKES 96
#define COILS 32
#define NODES READOUT_LEN*NUM_SPOKES

#define PARAM_H1_ALPHA 1.0f
#define PARAM_H1_ITER 8
#define PARAM_H1_REGSPAN 1.0f

#define PARAM_SMOOTH_ITER 40

#define PARAM_ALPHA1 8e-5f
#define PARAM_ALPHA0 2.0f*PARAM_ALPHA1
#define PARAM_TGV_ITER 100
#define PARAM_REGSPAN 4.0f

#define PARAM_STRING_LEN 512

// gradient operator norm squared
#define GRAD_NORM_SQR 8.0f
// norm estimate for total operator (assuming norm(NFFT) <= 1)
#define TOTAL_NORM_EST 1.0f + sqrtf(GRAD_NORM_SQR) + GRAD_NORM_SQR
#define STEP_LENGTH 1.0f/sqrtf(TOTAL_NORM_EST)
#define PARAM_SIGMA STEP_LENGTH
#define PARAM_TAU STEP_LENGTH
#define PARAM_NORMEST 1.0f

#define PARAM_H1_SIGMA STEP_LENGTH
#define PARAM_H1_TAU STEP_LENGTH
#define PARAM_H1_NORMEST 1.0f

#define CUDA_BLOCK_LEN 16
#define CUDA_LIN_BLOCK 512
#define CUDA_LIN_GRID 512
#define CUDA_2D_BLOCK_X 32
#define CUDA_2D_BLOCK_Y 16
#define CUDA_2D_SKIP_X 16
#define CUDA_2D_SKIP_Y 15

#define IDX2(A,s1,s2) A[(k2 + (s2))*N1 + (k1 + (s1))]
#define L_IDX2(A,k1,k2) A[threadIdx.y+(k2)][threadIdx.x+(k1)]

#define IS_RANGE (k1 < N1) && (k2 < N2)
#define IS_RANGE2(s1,s2) (k1 + (s1) >= 0) && (k2 + (s2) >= 0) \
                        && (k1 + (s1) < N1) && (k2 + (s2) < N2)

#define IS_L_RANGE (threadIdx.x < CUDA_2D_SKIP_X) \
                 && (threadIdx.y < CUDA_2D_SKIP_Y) 

#define getf(_foo,_i,_j,_k) (_foo.data.f[_i+_foo.dim[0]*(_j+_foo.dim[1]*(_k))])
#define getfc(_foo,_i,_j,_k) (_foo.data.fc[_i+_foo.dim[0]*(_j+_foo.dim[1]*(_k))])
#define sizef(_foo) (_foo.dim[0]*_foo.dim[1]*_foo.dim[2]*sizeof(float))
#define sizefc(_foo) (_foo.dim[0]*_foo.dim[1]*_foo.dim[2]*sizeof(float_complex))

typedef float float_complex[2];

typedef struct array_t {
  unsigned int dim[3];
  union data_t {
    unsigned char *uc;
    float *f;
    float_complex *fc;
    double *d;
  } data;
} array;

typedef struct problem_t {
  // filenames
  char data_file[PARAM_STRING_LEN], coord_file[PARAM_STRING_LEN],
    weight_file[PARAM_STRING_LEN],  result_file[PARAM_STRING_LEN];
  // size constants
  unsigned int N1, N2, nodes, coils;  
  // tgv parameters
  int tgv_iter;
  float normest;
  float sigma, tau;      // tgv denoising norm parameters
  float alpha0, alpha1;  // tgv denoising regularization parameters
  float regspan, regreduction;
  // h1 iteration parameters
  int h1_iter;
  float h1_normest, h1_alpha, h1_sigma, h1_tau; 
  float h1_regspan, h1_regreduction;
  // smoothing iteration parameter
  int smooth_iter;
  // given data
  array img, sens, data, coords, weights;
  // data associated with transformation
  array nfft_src, nfft_dest;
  nfft_plan plan;
  // data associated with primal-dual algorithm
  // additional primal variables
  array p_x, p_y, data_grad;
  // dual variables
  array dual_v, xi_x, xi_y, eta_xx, eta_xy, eta_yy; 
  array h1_xi_x, h1_xi_y;
  // shadow variables
  array img_shadow, p_x_shadow, p_y_shadow, sens_shadow;
} problem;

void display_activate(problem *);

// allocates host memory for a 3D array
void new_array(array *ary, unsigned int N1, unsigned int N2, 
	       unsigned int N3, size_t size) {
  
  assert(ary != NULL);
  
  // save dimensions
  ary->dim[0] = N1;
  ary->dim[1] = N2;
  ary->dim[2] = N3;
  
  // allocate space
  ary->data.f = (float *)nfft_malloc(size*N1*N2*N3);
  assert(ary->data.f != NULL);
}

// allocates memory for a 3D device array
void device_new_array(array *ary, unsigned int N1, unsigned int N2, 
	       unsigned int N3, size_t size) {
  
  assert(ary != NULL);
  
  // save dimensions
  ary->dim[0] = N1;
  ary->dim[1] = N2;
  ary->dim[2] = N3;
  
  // allocate space
  CUDA_SAFE_CALL( cudaMalloc( (void **)&ary->data.f, size*N1*N2*N3) );
  printf("alloc %12d bytes: %012X\n", size*N1*N2*N3, ary->data.f);
  CUDA_SAFE_CALL( cudaMemset( ary->data.f, 0, size*N1*N2*N3) );
}

void device_dump_array_fc(array *ary) {
  float_complex *hostmem;
  hostmem = (float_complex *)malloc(sizeof(float_complex)*
				    ary->dim[0]*ary->dim[1]*ary->dim[2]);
  CUDA_SAFE_CALL( cudaMemcpy(hostmem, ary->data.fc,
			     sizeof(float_complex)*
			     ary->dim[0]*ary->dim[1]*ary->dim[2],
			     cudaMemcpyDeviceToHost) );
  for(int j = 0; j < ary->dim[2]; j++)
    for(int k = 0; k < ary->dim[1]; k++)
      for(int l = 0; l < ary->dim[0]; l++)
	{
	  printf("(%d,%d,%d):(%g,%g) ", l,k,j,
		 hostmem[(j*ary->dim[1]+k)*ary->dim[0]+l][0],
		 hostmem[(j*ary->dim[1]+k)*ary->dim[0]+l][1]);
	}
  free(hostmem);
}

// reads the data in filename and copies it into device memory
void device_reader(void *data, char *filename, unsigned long size)	
{
  FILE *file;
  void *host_data;
  unsigned long check;
  unsigned long fileLength;
  
  // open the file
  file = fopen (filename, "rb");
  assert (file != NULL);
  
  // check that the file size matches how much we want to read
  fseek (file, 0, SEEK_END);
  fileLength = ftell (file);
  printf("%s %lu %lu\n", filename, fileLength, size);
  // assert(fileLength == size);
  fseek (file, 0, SEEK_SET);

  // allocate space on host
  host_data = nfft_malloc(size);
  assert(host_data != NULL);
  
  // read the data
  check = fread(host_data, 1, size, file);
  assert (check == size);
  fclose(file);
  
  // copy host data to device
  CUDA_SAFE_CALL( cudaMemcpy(data, host_data, size,  
			     cudaMemcpyHostToDevice) ); 

  nfft_free(host_data);
}

// writes float_complex array into a file
void device_writer (char *filename, array *ary)	{
  FILE *file;
  array host_ary;
  unsigned long check;
  int i;

  assert(ary != NULL);

  // allocate host memory and copy data
  new_array(&host_ary, ary->dim[0], ary->dim[1], ary->dim[2], 
	    sizeof(float));
  CUDA_SAFE_CALL( cudaMemcpy(host_ary.data.f, ary->data.f, 
			     ary->dim[0]*ary->dim[1]*ary->dim[2]
			     *sizeof(float), cudaMemcpyDeviceToHost) );
  
  // open the file
  file = fopen (filename, "wb");
  assert (file != NULL);
  
  // write the data
  for (i=0;i<ary->dim[2];i++) 
    {
      check = fwrite(&getf(host_ary,0,0,i), 
		     ary->dim[0]*ary->dim[1], 
		     sizeof(float), file);
      //assert(check == ary->dim[0]);
    }
  fclose(file);
  
  // Free memory
  nfft_free(host_ary.data.f);
}

void read_problem_param(problem *p, char *fname) {
  int ImageDimensionX, ImageDimensionY;
  char ImageFile[PARAM_STRING_LEN];
  int DataNodes, DataCoils;
  char DataFile[PARAM_STRING_LEN], CoordinateFile[PARAM_STRING_LEN],
    WeightFile[PARAM_STRING_LEN];
  float SensAlpha, SensRegspan;
  int SensIterations, SensSmoothing;
  float TGVAlpha, TGVAlphaFac, TGVRegspan;
  int TGVIterations;

  ParseParamInt(fname, ImageDimensionX);
  ParseParamInt(fname, ImageDimensionY);
  ParseParamString(fname, ImageFile, PARAM_STRING_LEN);

  ParseParamInt(fname, DataNodes);
  ParseParamInt(fname, DataCoils);
  ParseParamString(fname, DataFile, PARAM_STRING_LEN);
  ParseParamString(fname, CoordinateFile, PARAM_STRING_LEN);
  ParseParamString(fname, WeightFile, PARAM_STRING_LEN);

  ParseParamFloat(fname, SensAlpha);
  ParseParamFloat(fname, SensRegspan);
  ParseParamInt(fname, SensIterations);
  ParseParamInt(fname, SensSmoothing);
  
  ParseParamFloat(fname, TGVAlpha);
  ParseParamFloat(fname, TGVAlphaFac);
  ParseParamFloat(fname, TGVRegspan);
  ParseParamInt(fname, TGVIterations);

  p->N1 = ImageDimensionX > 0 ? ImageDimensionX : 1;
  p->N2 = ImageDimensionY > 0 ? ImageDimensionY : 1;
  strncpy(p->result_file, ImageFile, PARAM_STRING_LEN);
  
  p->nodes = DataNodes > 0 ? DataNodes : 1;
  p->coils = DataCoils > 0 ? DataCoils : 1;
  strncpy(p->data_file, DataFile, PARAM_STRING_LEN);
  strncpy(p->coord_file, CoordinateFile, PARAM_STRING_LEN);
  strncpy(p->weight_file, WeightFile, PARAM_STRING_LEN);
  
  p->h1_alpha = SensAlpha > 0 ? SensAlpha : 1e-6f;
  p->h1_regspan = SensRegspan >= 1.0f ? SensRegspan : 1.0f;
  p->h1_iter = SensIterations > 0 ? SensIterations : 1;
  p->smooth_iter = SensSmoothing > 0 ? SensSmoothing : 1;

  p->alpha1 = TGVAlpha > 0 ? TGVAlpha : 1e-6f;
  p->alpha0 = (TGVAlphaFac > 0 ? TGVAlphaFac : 2.0f)*TGVAlpha;
  p->regspan = TGVRegspan >= 1.0f ? TGVRegspan : 1.0f;
  p->tgv_iter = TGVIterations > 0 ? TGVIterations : 1;
}

// initialize problem data
void init_problem(problem *p, char *paramfile) {
  assert(p != NULL);

  // get parameters
  read_problem_param(p, paramfile);

  // h1 iteration parameters
  p->h1_sigma = PARAM_H1_SIGMA;
  p->h1_tau = PARAM_H1_TAU;
  p->h1_normest = PARAM_H1_NORMEST;
  p->h1_regreduction = powf(p->h1_regspan, 1.0f/((float)p->h1_iter));

  // tgv iteration parameters
  p->sigma = PARAM_SIGMA;
  p->tau = PARAM_TAU;
  p->normest = PARAM_NORMEST;
  p->regreduction = powf(p->regspan, 1.0f/((float)p->tgv_iter));
  
  // allocate space for data
  device_new_array(&p->img, p->N1, p->N2, 1, sizeof(float));
  device_new_array(&p->sens, p->N1, p->N2, p->coils, sizeof(float_complex));
  device_new_array(&p->data, p->nodes, p->coils, 1, sizeof(float_complex));
  device_new_array(&p->coords, 2, p->nodes, 1, sizeof(float));
  device_new_array(&p->weights, p->nodes, 1, 1, sizeof(float));
  
  // allocate space for transformation
  device_new_array(&p->nfft_src, p->N2, p->N1, p->coils, sizeof(float_complex));
  device_new_array(&p->nfft_dest, p->nodes, p->coils, 1, sizeof(float_complex));

  // variables for primal-dual algorithm

  // additional primal variables
  device_new_array(&p->p_x, p->N1, p->N2, 1, sizeof(float));
  device_new_array(&p->p_y, p->N1, p->N2, 1, sizeof(float));
  device_new_array(&p->data_grad, p->N1, p->N2, 1, sizeof(float));

  // dual variables
  device_new_array(&p->dual_v, p->nodes, p->coils, 1, sizeof(float_complex));
  device_new_array(&p->xi_x, p->N1, p->N2, 1, sizeof(float));
  device_new_array(&p->xi_y, p->N1, p->N2, 1, sizeof(float));
  device_new_array(&p->eta_xx, p->N1, p->N2, 1, sizeof(float));
  device_new_array(&p->eta_xy, p->N1, p->N2, 1, sizeof(float));
  device_new_array(&p->eta_yy, p->N1, p->N2, 1, sizeof(float));
  device_new_array(&p->h1_xi_x, p->N1, p->N2, p->coils, sizeof(float_complex));
  device_new_array(&p->h1_xi_y, p->N1, p->N2, p->coils, sizeof(float_complex));

  // shadow variables
  device_new_array(&p->img_shadow, p->N1, p->N2, 1, sizeof(float));
  device_new_array(&p->p_x_shadow, p->N1, p->N2, 1, sizeof(float));
  device_new_array(&p->p_y_shadow, p->N1, p->N2, 1, sizeof(float));
  device_new_array(&p->sens_shadow, p->N1, p->N2, p->coils, 
		   sizeof(float_complex));
  
}

void mult_data_weights(array *, array *);

void read_problem_data(problem *p) {
  assert(p != NULL);
  
  device_reader(p->data.data.fc, p->data_file, sizefc(p->data));
  device_reader(p->coords.data.f, p->coord_file, sizef(p->coords));
  device_reader(p->weights.data.f, p->weight_file, sizef(p->weights));
  mult_data_weights(&p->data, &p->weights);
}

void init_nfft_plan(problem *p) {
  // init plan
  nfft_init_2d(&p->plan, p->N1, p->N2, p->nodes, p->coils,
               sizeof(float_complex)*p->nodes,
               p->nfft_src.data.fc, p->nfft_dest.data.fc, NULL,
               p->coords.data.f);
}

/////////////////////////////////////////////////
// general purpose routines
/////////////////////////////////////////////////

// norm square routine
static float *device_normsqr = 0;

__global__ void normsqr_A(float *ptr, int numel, float *device_normsqr) {
  __shared__ float acc[CUDA_LIN_BLOCK];
  float val;
  int bid, tid, i;

  // init ids
  tid = threadIdx.x;
  bid = blockIdx.x;
  acc[tid] = 0;

  i = tid + bid*CUDA_LIN_BLOCK;
  ptr += i;
  
  // sum every BLOCK*GRID element in acc[tid]
  for (; i < numel; i += CUDA_LIN_BLOCK*CUDA_LIN_GRID) {
    val = *ptr;
    acc[tid] += val*val;
    ptr += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
  }
  
  // reduce accumulator
  for (i = CUDA_LIN_BLOCK >> 1; i > 0; i >>= 1)
    {
      __syncthreads();
      if (tid < i) 
	acc[tid] += acc[tid + i];
    }
  __syncthreads();
  
  // store
  if (tid == 0)
    device_normsqr[blockIdx.x] = acc[0];   
}

float normsqr_f(array *ary) {
  dim3 dimGrid, dimBlock;
  float normsqr[CUDA_LIN_GRID], res;
  int i, numel;
  
  // allocate device memory for device_normsqr
  if (!device_normsqr)
    CUDA_SAFE_CALL( cudaMalloc( (void **)&device_normsqr, 
				sizeof(float)*CUDA_LIN_GRID) );
  
  // init dimensions
  numel = ary->dim[0]*ary->dim[1]*ary->dim[2]; 
  dimBlock.x = CUDA_LIN_BLOCK;
  dimBlock.y = 1;
  dimBlock.z = 1;
  dimGrid.x = CUDA_LIN_GRID;
  dimGrid.y = 1;
  dimGrid.z = 1;
  normsqr_A<<<dimGrid, dimBlock>>>(ary->data.f, numel, device_normsqr);
  CUT_CHECK_ERROR("normsqr_A() failed.");
  
  // copy results to host memory
  CUDA_SAFE_CALL( cudaMemcpy(&normsqr, device_normsqr, 
			     CUDA_LIN_GRID*sizeof(float), 
			     cudaMemcpyDeviceToHost) );
  
  // sum up the bloody rest
  res = 0.0f;
  for(i=0; i<CUDA_LIN_GRID; i++)
    res += normsqr[i];
  
  return(res);
}

float normsqr_fc(array *ary) {
  dim3 dimGrid, dimBlock;
  float normsqr[CUDA_LIN_GRID], res;
  int i, numel;
  
  // allocate device memory for device_normsqr
  if (!device_normsqr)
    CUDA_SAFE_CALL( cudaMalloc( (void **)&device_normsqr, 
				sizeof(float)*CUDA_LIN_GRID) );
  
  // init dimensions
  numel = 2*ary->dim[0]*ary->dim[1]*ary->dim[2]; 
  dimBlock.x = CUDA_LIN_BLOCK;
  dimBlock.y = 1;
  dimBlock.z = 1;
  dimGrid.x = CUDA_LIN_GRID;
  dimGrid.y = 1;
  dimGrid.z = 1;
  normsqr_A<<<dimGrid, dimBlock>>>(ary->data.f, numel, device_normsqr);
  CUT_CHECK_ERROR("normsqr_A() failed.");
  
  // copy results to host memory
  CUDA_SAFE_CALL( cudaMemcpy(&normsqr, device_normsqr, 
			     CUDA_LIN_GRID*sizeof(float), 
			     cudaMemcpyDeviceToHost) );
  
  // sum up the bloody rest
  res = 0.0f;
  for(i=0; i<CUDA_LIN_GRID; i++)
    res += normsqr[i];
  
  return(res);
}

/////////////////////////////////////////////////
// transformation-related routines
/////////////////////////////////////////////////

__global__ void mult_data_weights_A(float *data, float *weight, int nodes) 
{
  __shared__ float buf0[2*CUDA_LIN_BLOCK], buf1[CUDA_LIN_BLOCK];
  int i, j;
  
  // init id
  i = threadIdx.x + CUDA_LIN_BLOCK*blockIdx.x;
  data += 2*nodes*blockIdx.y;
  j = 2*CUDA_LIN_BLOCK*blockIdx.x + threadIdx.x;

  // read
  if (j < 2*nodes)
    buf0[threadIdx.x] = data[j];
  if (j + CUDA_LIN_BLOCK < 2*nodes)
    buf0[threadIdx.x+CUDA_LIN_BLOCK] = data[j+CUDA_LIN_BLOCK];
  if (i < nodes)
    buf1[threadIdx.x] = weight[i];
  __syncthreads();
    
  // compute 
  buf0[threadIdx.x*2] *= buf1[threadIdx.x];
  buf0[threadIdx.x*2+1] *= buf1[threadIdx.x];
  __syncthreads();
  
  // write
  if (j < 2*nodes)
    data[j] = buf0[threadIdx.x];
  if (j + CUDA_LIN_BLOCK < 2*nodes)
    data[j+CUDA_LIN_BLOCK] = buf0[threadIdx.x+CUDA_LIN_BLOCK];
}

void mult_data_weights(array *fc, array *f) {
  dim3 dimBlock, dimGrid;
  int nodes, coils;
 
  nodes = fc->dim[0];
  coils = fc->dim[1];

  // init dimensions
  dimBlock.x = CUDA_LIN_BLOCK; dimBlock.y = 1; dimBlock.z = 1;
  dimGrid.x = (nodes + CUDA_LIN_BLOCK-1)/CUDA_LIN_BLOCK; 
  dimGrid.y = coils; dimGrid.z = 1;
  mult_data_weights_A<<<dimGrid, dimBlock>>>(fc->data.f, f->data.f, nodes);
  CUT_CHECK_ERROR("mult_data_weights_A() failed.");
}

// NFFT transform img -> nfft_dest
__global__ void transform_A(int n0, int n1, 
			    float_complex *dest, float *src, 
			    float_complex *sens, int z_stride) {
  int bidxy, bidxz;
  int k0, k1, coil, l0, l1;
  __shared__ float src_l[CUDA_BLOCK_LEN][CUDA_BLOCK_LEN];
  __shared__ float_complex sens_l[CUDA_BLOCK_LEN][CUDA_BLOCK_LEN];
  __shared__ float_complex dest_l[CUDA_BLOCK_LEN][CUDA_BLOCK_LEN];

  // compute indices
  bidxz = floorf(((float)blockIdx.y)/((float)z_stride));
  bidxy = blockIdx.y-z_stride*bidxz;
  
  l0 = threadIdx.x;
  l1 = threadIdx.y;
  k0 = blockIdx.x*CUDA_BLOCK_LEN;
  k1 = bidxy*CUDA_BLOCK_LEN;
  coil = bidxz;
  
  // adjust pointers according to coil
  sens += coil*n0*n1;
  dest += coil*n0*n1;

  // read image and sensitivities
  if ((k0+l0 < n0) && (k1+l1 < n1))
    src_l[l0][l1] = src[(k1+l1)*n0 + (k0+l0)];
  if ( ((k0+l0/2) < n0) && (k1+l1 < n1) ) 
    sens_l[l0/2][l1][l0 & 1] = ((float *)&sens[(k1+l1)*n0 + k0])[l0];
  if ( ((k0+(l0+CUDA_BLOCK_LEN)/2) < n0) && (k1+l1 < n1) ) 
    sens_l[(l0+CUDA_BLOCK_LEN)/2][l1][l0 & 1] 
      = ((float *)&sens[(k1+l1)*n0 + k0])[l0+CUDA_BLOCK_LEN];
  __syncthreads();
  
  // compute product
  dest_l[l0][l1][0] = src_l[l0][l1]*sens_l[l0][l1][0]; 
  dest_l[l0][l1][1] = src_l[l0][l1]*sens_l[l0][l1][1];
  __syncthreads();

  // write result
  if ( ((k0+l0/2) < n0) && (k1+l1 < n1) ) 
    ((float *)&dest[(k1+l1)*n0 + k0])[l0] = dest_l[l0/2][l1][l0 & 1];
  if ( ((k0+(l0+CUDA_BLOCK_LEN)/2) < n0) && (k1+l1 < n1) ) 
    ((float *)&dest[(k1+l1)*n0 + k0])[l0+CUDA_BLOCK_LEN] = 
      dest_l[(l0+CUDA_BLOCK_LEN)/2][l1][l0 & 1];
}

void transform(problem *p) {
  dim3 dimBlock, dimGrid;
  int N1, N2, coils, z_stride;
  float srcnorm, destnorm;
  
  // get source norm
  srcnorm = sqrtf(normsqr_f(&p->img));
  
  // init dimensions
  N1 = p->N1; N2 = p->N2; coils = p->coils;
  dimBlock.x = CUDA_BLOCK_LEN;
  dimBlock.y = CUDA_BLOCK_LEN;
  dimBlock.z = 1;
  dimGrid.x = (N1 + CUDA_BLOCK_LEN - 1)/dimBlock.x;
  z_stride = (N2 + CUDA_BLOCK_LEN - 1)/dimBlock.y;
  dimGrid.y = coils*z_stride;
  dimGrid.z = 1;
  // multiply with sensitivities
  transform_A<<<dimGrid, dimBlock>>>(N1, N2,  
				     p->nfft_src.data.fc, p->img.data.f,
				     p->sens.data.fc, z_stride);
  CUT_CHECK_ERROR("transform_A() failed.");

  // perform transform
  nfft_trafo_2d(&p->plan);
  
  // multiply with weights
  mult_data_weights(&p->nfft_dest, &p->weights);
  
  // get destination norm
  destnorm = sqrtf(normsqr_fc(&p->nfft_dest));

  // adjust estimate if necessary
  if ((srcnorm > 0) && (destnorm > p->normest*srcnorm)) {
    printf("srcnorm = %g, destnorm = %g\n", srcnorm, destnorm);
    p->normest = destnorm/srcnorm;
    printf("WARNING: adjusting norm estimate to %g.\n", p->normest);
  }
}

__global__ void transform_adjoint_A(int n0, int n1, int coils, 
				    float *dest, float_complex *src, 
				    float_complex *sens) {
  int k0, k1, l0, l1, i;
  __shared__ float_complex src_l[CUDA_BLOCK_LEN][CUDA_BLOCK_LEN];
  __shared__ float_complex sens_l[CUDA_BLOCK_LEN][CUDA_BLOCK_LEN];
  __shared__ float dest_l[CUDA_BLOCK_LEN][CUDA_BLOCK_LEN];

  // compute indices
  l0 = threadIdx.x;
  l1 = threadIdx.y;
  k0 = blockIdx.x*CUDA_BLOCK_LEN;
  k1 = blockIdx.y*CUDA_BLOCK_LEN;

  // clear buffer
  dest_l[l0][l1] = 0.0f;

  // for each coil
  for(i=0; i<coils; i++) {
    // read nfft data
    if ( ((k0+l0/2) < n0) && (k1+l1 < n1) ) 
      src_l[l0/2][l1][l0 & 1] = ((float *)&src[(k1+l1)*n0 + k0])[l0];
    if ( ((k0+(l0+CUDA_BLOCK_LEN)/2) < n0) && (k1+l1 < n1) ) 
      src_l[(l0+CUDA_BLOCK_LEN)/2][l1][l0 & 1] = 
	((float *)&src[(k1+l1)*n0 + k0])[l0+CUDA_BLOCK_LEN];
    
    // read sensitivities
    if ( ((k0+l0/2) < n0) && (k1+l1 < n1) ) 
      sens_l[l0/2][l1][l0 & 1] = ((float *)&sens[(k1+l1)*n0 + k0])[l0];
    if ( ((k0+(l0+CUDA_BLOCK_LEN)/2) < n0) && (k1+l1 < n1) ) 
      sens_l[(l0+CUDA_BLOCK_LEN)/2][l1][l0 & 1] 
	= ((float *)&sens[(k1+l1)*n0 + k0])[l0+CUDA_BLOCK_LEN];
    __syncthreads();

    // add scalar product to local buffer
    dest_l[l0][l1] +=  src_l[l0][l1][0]*sens_l[l0][l1][0]
      + src_l[l0][l1][1]*sens_l[l0][l1][1];
    __syncthreads();
    
    src += n0*n1;
    sens += n0*n1;
  }
  
  // write image 
  if ((k0+l0 < n0) && (k1+l1 < n1))
    dest[(k1+l1)*n0 + (k0+l0)] = dest_l[l0][l1];
}

// NFFT adjoint dual_v -> data_grad
void transform_adjoint(problem *p) {
  dim3 dimBlock, dimGrid;
  int N1, N2, coils;
  float srcnorm, destnorm;
  
  // get source norm
  srcnorm = sqrtf(normsqr_fc(&p->dual_v));

  // copy into nfft_dest
  CUDA_SAFE_CALL( cudaMemcpy(p->nfft_dest.data.fc, p->dual_v.data.fc,
                             p->dual_v.dim[0]*p->dual_v.dim[1]*p->dual_v.dim[2]
                             *sizeof(float_complex),
                             cudaMemcpyDeviceToDevice) );
  
  // multiply with weights
  mult_data_weights(&p->nfft_dest, &p->weights);
  
  // perform transform adjoint
  nfft_adjoint_2d(&p->plan);
  
  // init dimensions
  N1 = p->N1; N2 = p->N2; coils = p->coils;
  dimBlock.x = CUDA_BLOCK_LEN;
  dimBlock.y = CUDA_BLOCK_LEN;
  dimBlock.z = 1;
  dimGrid.x = (N1 + CUDA_BLOCK_LEN - 1)/dimBlock.x;
  dimGrid.y = (N2 + CUDA_BLOCK_LEN - 1)/dimBlock.y;
  dimGrid.z = 1;
  // multiply with sensitivities and accumulate
  transform_adjoint_A<<<dimGrid, dimBlock>>>(N1, N2, coils,  
					     p->data_grad.data.f,
					     p->nfft_src.data.fc,
					     p->sens.data.fc);
  CUT_CHECK_ERROR("transform_adjoint_A() failed.");
  
  // get destination norm
  destnorm = sqrtf(normsqr_f(&p->data_grad));
  
  // adjust estimate if necessary
  if ((srcnorm > 0) && (destnorm > p->normest*srcnorm)) {
    printf("srcnorm = %g, destnorm = %g\n", srcnorm, destnorm);
    p->normest = destnorm/srcnorm;
    printf("WARNING: adjusting norm estimate to %g.\n", p->normest);
  }
}

void h1_transform(problem *p) {
  float srcnorm, destnorm;
  
  // get source norm
  srcnorm = sqrtf(normsqr_fc(&p->sens));
  
  // copy sens -> nfft_src
  CUDA_SAFE_CALL( cudaMemcpy(p->nfft_src.data.fc, p->sens.data.fc, 
			     p->sens.dim[0]*p->sens.dim[1]*p->sens.dim[2]
			     *sizeof(float_complex),  
			     cudaMemcpyDeviceToDevice) ); 
  
  // perform transform
  nfft_trafo_2d(&p->plan);
  
  // multiply with weights
  mult_data_weights(&p->nfft_dest, &p->weights);
  
  // get destination norm
  destnorm = sqrtf(normsqr_fc(&p->nfft_dest));

  // adjust estimate if necessary
  if ((srcnorm > 0) && (destnorm > p->h1_normest*srcnorm)) {
    printf("srcnorm = %g, destnorm = %g\n", srcnorm, destnorm);
    p->h1_normest = destnorm/srcnorm;
    printf("WARNING: adjusting norm estimate to %g.\n", p->h1_normest);
  }
}

void h1_transform_adjoint(problem *p) {
  float srcnorm, destnorm;
  
  // get source norm
  srcnorm = sqrtf(normsqr_fc(&p->dual_v));

  // copy into nfft_dest
  CUDA_SAFE_CALL( cudaMemcpy(p->nfft_dest.data.fc, p->dual_v.data.fc,
                             p->dual_v.dim[0]*p->dual_v.dim[1]*p->dual_v.dim[2]
                             *sizeof(float_complex),
                             cudaMemcpyDeviceToDevice) );
  
  // multiply with weights
  mult_data_weights(&p->nfft_dest, &p->weights);
  
  // perform transform adjoint
  nfft_adjoint_2d(&p->plan);
   
  // get destination norm
  destnorm = sqrtf(normsqr_fc(&p->nfft_src));
  
  // adjust estimate if necessary
  if ((srcnorm > 0) && (destnorm > p->h1_normest*srcnorm)) {
    printf("srcnorm = %g, destnorm = %g\n", srcnorm, destnorm);
    p->h1_normest = destnorm/srcnorm;
    printf("WARNING: adjusting norm estimate to %g.\n", p->h1_normest);
  }
}

//////////////////////////////////////////////////////////////////////
// h1 minimization algorithm related functions
//////////////////////////////////////////////////////////////////////

// swaps shadow image
void h1_swap_shadow_sens(problem *p) {
  float *f;
  
  // swap sens
  f = p->sens_shadow.data.f;
  p->sens_shadow.data.f = p->sens.data.f;
  p->sens.data.f = f;
}

// updates dual variable v according to
// v^{n+1} = (v^n + sigma_v*(K img - f))/(1 + sigma)
__global__ void h1_update_dual_v_A(float *dest, float *src0, float *src1, 
				int numel, float sigma, float sigmap1inv,
				float normestinv) {
  int i;
  
  // init id
  i = threadIdx.x + CUDA_LIN_BLOCK*blockIdx.x;
  dest += i; src0 += i; src1 += i;
  
  for (; i < numel; i += CUDA_LIN_BLOCK*CUDA_LIN_GRID) {
    // compute
    *dest = (*dest + sigma*((*src0)*normestinv - (*src1)))*sigmap1inv;
    // increment
    dest += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    src0 += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    src1 += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
  }
}

void h1_update_dual_v(problem *p) {
  dim3 dimBlock, dimGrid;
  int numel;
  float sigma, normestinv;
  
  sigma = p->h1_sigma;
  normestinv = 1.0f/p->h1_normest;

  // compute transform
  h1_transform(p);

  // compute update
  numel = 2*p->dual_v.dim[0]*p->dual_v.dim[1]*p->dual_v.dim[2];
  // init dimensions
  dimBlock.x = CUDA_LIN_BLOCK; dimBlock.y = 1; dimBlock.z = 1;
  dimGrid.x = CUDA_LIN_GRID; dimGrid.y = 1; dimGrid.z = 1;
  h1_update_dual_v_A<<<dimGrid, dimBlock>>>(p->dual_v.data.f, 
					    p->nfft_dest.data.f, 
					    p->data.data.f, numel, sigma,
					    1.0f/(sigma + 1.0f), normestinv);
  CUT_CHECK_ERROR("h1_update_dual_v_A() failed.");
}

__global__ void h1_update_xi_A(int N1, int N2, 
			       float *destx, float *desty, 
			       float *src, float sigma, 
			       float sigma_alphap1inv,
			       int z_stride) {
  int bidxy, bidxz, k1, k2;
  __shared__ float src_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X];
  float valx, valy;

  bidxz = floorf(((float)blockIdx.y)/((float)z_stride));
  bidxy = blockIdx.y-z_stride*bidxz;

  k1 = blockIdx.x*CUDA_2D_SKIP_X + threadIdx.x;
  k2 = bidxy*CUDA_2D_SKIP_Y + threadIdx.y;
  
  destx += N1*N2*bidxz; 
  desty += N1*N2*bidxz; 
  src += N1*N2*bidxz; 

  // get source buffer
  if (IS_RANGE)
    L_IDX2(src_l,0,0) = IDX2(src,0,0);
  __syncthreads();

  if (IS_L_RANGE) {
    // derivative w.r.t x
    valx = (k1 < N1-2) ? sigma*(L_IDX2(src_l,2,0) - L_IDX2(src_l,0,0)) : 0.0f; 
    
    // derivative w.r.t y
    valy = (k2 < N2-1) ? sigma*(L_IDX2(src_l,0,1) - L_IDX2(src_l,0,0)) : 0.0f; 
    
    if (IS_RANGE) {
      IDX2(destx,0,0) = (IDX2(destx,0,0) + valx)*sigma_alphap1inv;
      IDX2(desty,0,0) = (IDX2(desty,0,0) + valy)*sigma_alphap1inv;
    }
  }
}

void h1_update_xi(problem *p) {
  dim3 dimBlock, dimGrid;
  int z_stride, coils;
  float sigma, alpha;
  
  sigma = p->h1_sigma;
  alpha = p->h1_alpha;
  coils = p->coils;

  // xi = xi + sigma*\grad img
  dimBlock.x = CUDA_2D_BLOCK_X;
  dimBlock.y = CUDA_2D_BLOCK_Y;
  dimBlock.z = 1;
  dimGrid.x = (2*p->N1 + CUDA_2D_SKIP_X - 1)/CUDA_2D_SKIP_X;
  z_stride = (p->N2 + CUDA_2D_SKIP_Y - 1)/CUDA_2D_SKIP_Y;
  dimGrid.y = coils*z_stride;
  dimGrid.z = 1;
  h1_update_xi_A<<<dimGrid, dimBlock>>>(2*p->N1, p->N2, 
					p->h1_xi_x.data.f, p->h1_xi_y.data.f,
					p->sens.data.f, sigma, 
					1.0f/(1.0f + sigma/alpha), 
					z_stride);
  CUT_CHECK_ERROR("h1_update_xi_A() failed.");
}

__global__ void h1_update_sens_A(int N1, int N2, 
				 float *dest, float *srcimg, float *srcv,
				 float *srcx, float *srcy, 
				 float tau, float normestinv, int z_stride) {
  int bidxy, bidxz, k1, k2;
  __shared__ float
    srcimg_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X],
    srcv_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X],
    srcx_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X],
    srcy_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X];
  float val;

  // compute indices
  bidxz = floorf(((float)blockIdx.y)/((float)z_stride));
  bidxy = blockIdx.y-z_stride*bidxz;
  
  k1 = blockIdx.x*CUDA_2D_SKIP_X + threadIdx.x;
  k2 = bidxy*CUDA_2D_SKIP_Y + threadIdx.y;
  
  dest += N1*N2*bidxz; 
  srcimg += N1*N2*bidxz; 
  srcv += N1*N2*bidxz; 
  srcx += N1*N2*bidxz;
  srcy += N1*N2*bidxz;  

  // get source buffers
  if (IS_RANGE2(-16,-1)) {
    L_IDX2(srcv_l,0,0) = IDX2(srcv,-16,-1);
    L_IDX2(srcimg_l,0,0) = IDX2(srcimg,-16,-1);
    L_IDX2(srcx_l,0,0) = IDX2(srcx,-16,-1);
    L_IDX2(srcy_l,0,0) = IDX2(srcy,-16,-1);
  }
  __syncthreads();
  
  if (IS_L_RANGE) {
    // compute divergence
    val = (((float)(k1 < N1-2))*L_IDX2(srcx_l,16,1)
	   - ((float)(k1 > 1))*L_IDX2(srcx_l,14,1));
    val += (((float)(k2 < N2-1))*L_IDX2(srcy_l,16,1)
	    - ((float)(k2 > 0))*L_IDX2(srcy_l,16,0));
    
    // dest = img + tau*(div - v)
    val = L_IDX2(srcimg_l,16,1) 
		+ tau*(val - normestinv*L_IDX2(srcv_l,16,1));
    
    if (IS_RANGE)
      IDX2(dest,0,0) = val;
  }
}

void h1_update_sens(problem *p) {
  dim3 dimBlock, dimGrid;
  int z_stride;
  float tau, normestinv;
  
  tau = p->h1_tau;  
  normestinv = 1.0f/p->h1_normest;
  
  // compute adjoint K^*v
  h1_transform_adjoint(p);
  
  // img = shadow_img + tau*(div xi - K^*v)
  dimBlock.x = CUDA_2D_BLOCK_X;
  dimBlock.y = CUDA_2D_BLOCK_Y;
  dimBlock.z = 1;
  dimGrid.x = (2*p->N1 + CUDA_2D_SKIP_X - 1)/CUDA_2D_SKIP_X;
  z_stride = (p->N2 + CUDA_2D_SKIP_Y - 1)/CUDA_2D_SKIP_Y;
  dimGrid.y = p->coils*z_stride;
  dimGrid.z = 1;
  h1_update_sens_A<<<dimGrid, dimBlock>>>(2*p->N1, p->N2,
					  p->sens.data.f, p->sens_shadow.data.f,
					  p->nfft_src.data.f, p->h1_xi_x.data.f,
					  p->h1_xi_y.data.f, tau, normestinv,
					  z_stride);
  CUT_CHECK_ERROR("h1_update_sens_A() failed.");
}

__global__ void h1_update_shadow_sens_A(float *dest, float *src, int numel) {
  int i;
  
  // init id
  i = threadIdx.x + CUDA_LIN_BLOCK*blockIdx.x;
  dest += i; src += i;
  
  for (; i < numel; i += CUDA_LIN_BLOCK*CUDA_LIN_GRID) {
    // compute
    *dest = 2.0f*(*src) - (*dest);
    // increment
    dest += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    src += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
  }
}

void h1_update_shadow_sens(problem *p) {
  dim3 dimBlock, dimGrid;
  int numel;

  // img part
  numel = 2*p->sens.dim[0]*p->sens.dim[1]*p->sens.dim[2];
  // init dimensions
  dimBlock.x = CUDA_LIN_BLOCK; dimBlock.y = 1; dimBlock.z = 1;
  dimGrid.x = CUDA_LIN_GRID; dimGrid.y = 1; dimGrid.z = 1;
  h1_update_shadow_sens_A<<<dimGrid, dimBlock>>>(p->sens_shadow.data.f, 
						 p->sens.data.f, numel);
  CUT_CHECK_ERROR("h1_update_shadow_sens_A() failed.");
}


__global__ void h1_smooth_sens_A(int N1, int N2, 
				 float *destx, float *desty, 
				 float *src, int z_stride) {
  int bidxy, bidxz, k1, k2;
  __shared__ float src_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X];
  float valx, valy;

  bidxz = floorf(((float)blockIdx.y)/((float)z_stride));
  bidxy = blockIdx.y-z_stride*bidxz;

  k1 = blockIdx.x*CUDA_2D_SKIP_X + threadIdx.x;
  k2 = bidxy*CUDA_2D_SKIP_Y + threadIdx.y;
  
  destx += N1*N2*bidxz; 
  desty += N1*N2*bidxz; 
  src += N1*N2*bidxz; 

  // get source buffer
  if (IS_RANGE)
    L_IDX2(src_l,0,0) = IDX2(src,0,0);
  __syncthreads();

  if (IS_L_RANGE) {
    // derivative w.r.t x
    valx = (k1 < N1-2) ? 0.5f*(L_IDX2(src_l,2,0) + L_IDX2(src_l,0,0)) : 0.0f; 
    
    // derivative w.r.t y
    valy = (k2 < N2-1) ? 0.5f*(L_IDX2(src_l,0,1) + L_IDX2(src_l,0,0)) : 0.0f; 
    
    if (IS_RANGE) {
      IDX2(destx,0,0) = valx;
      IDX2(desty,0,0) = valy;
    }
  }
}

__global__ void h1_smooth_sens_B(int N1, int N2, 
				 float *dest, float *srcx, float *srcy, 
				 int z_stride) {
  int bidxy, bidxz, k1, k2;
  __shared__ float
    srcx_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X],
    srcy_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X];
  float val;
  
  // compute indices
  bidxz = floorf(((float)blockIdx.y)/((float)z_stride));
  bidxy = blockIdx.y-z_stride*bidxz;
  
  k1 = blockIdx.x*CUDA_2D_SKIP_X + threadIdx.x;
  k2 = bidxy*CUDA_2D_SKIP_Y + threadIdx.y;
  
  dest += N1*N2*bidxz; 
  srcx += N1*N2*bidxz;
  srcy += N1*N2*bidxz;  

  // get source buffers
  if (IS_RANGE2(-16,-1)) {
    L_IDX2(srcx_l,0,0) = IDX2(srcx,-16,-1);
    L_IDX2(srcy_l,0,0) = IDX2(srcy,-16,-1);
  }
  __syncthreads();
  
  if (IS_L_RANGE) {
    // compute average
    val = (((float)(k1 < N1-2))*L_IDX2(srcx_l,16,1)
	   + ((float)(k1 > 1))*L_IDX2(srcx_l,14,1));
    val += (((float)(k2 < N2-1))*L_IDX2(srcy_l,16,1)
	    + ((float)(k2 > 0))*L_IDX2(srcy_l,16,0));
    val *= 0.5f;
    
    if (IS_RANGE)
      IDX2(dest,0,0) = val;
  }
}

void h1_smooth_sens(problem *p) {
  dim3 dimBlock, dimGrid;
  int z_stride;

  // first smoothing step
  dimBlock.x = CUDA_2D_BLOCK_X;
  dimBlock.y = CUDA_2D_BLOCK_Y;
  dimBlock.z = 1;
  dimGrid.x = (2*p->N1 + CUDA_2D_SKIP_X - 1)/CUDA_2D_SKIP_X;
  z_stride = (p->N2 + CUDA_2D_SKIP_Y - 1)/CUDA_2D_SKIP_Y;
  dimGrid.y = p->coils*z_stride;
  dimGrid.z = 1;
  h1_smooth_sens_A<<<dimGrid, dimBlock>>>(2*p->N1, p->N2, 
					  p->h1_xi_x.data.f, p->h1_xi_y.data.f,
					  p->sens.data.f, z_stride);
  CUT_CHECK_ERROR("h1_smooth_sens_A() failed.");
  
  h1_smooth_sens_B<<<dimGrid, dimBlock>>>(2*p->N1, p->N2, 
					  p->sens.data.f, p->h1_xi_x.data.f, 
					  p->h1_xi_y.data.f, z_stride);
  CUT_CHECK_ERROR("h1_smooth_sens_B() failed.");
}


__global__ void h1_sum_of_squares_A(int n0, int n1, int coils, 
				    float *dest, float_complex *sens) 
{
  int k0, k1, l0, l1, i;
  __shared__ float_complex sens_l[CUDA_BLOCK_LEN][CUDA_BLOCK_LEN];
  __shared__ float dest_l[CUDA_BLOCK_LEN][CUDA_BLOCK_LEN];

  // compute indices
  l0 = threadIdx.x;
  l1 = threadIdx.y;
  k0 = blockIdx.x*CUDA_BLOCK_LEN;
  k1 = blockIdx.y*CUDA_BLOCK_LEN;

  // clear buffer
  dest_l[l0][l1] = 0.0f;

  // for each coil
  for(i=0; i<coils; i++) {
    // read sensitivities
    if ( ((k0+l0/2) < n0) && (k1+l1 < n1) ) 
      sens_l[l0/2][l1][l0 & 1] = ((float *)&sens[(k1+l1)*n0 + k0])[l0];
    if ( ((k0+(l0+CUDA_BLOCK_LEN)/2) < n0) && (k1+l1 < n1) ) 
      sens_l[(l0+CUDA_BLOCK_LEN)/2][l1][l0 & 1] 
	= ((float *)&sens[(k1+l1)*n0 + k0])[l0+CUDA_BLOCK_LEN];
    __syncthreads();

    // add square product to local buffer
    dest_l[l0][l1] +=  sens_l[l0][l1][0]*sens_l[l0][l1][0]
      + sens_l[l0][l1][1]*sens_l[l0][l1][1];
    __syncthreads();
    
    sens += n0*n1;
  }
  
  // write image 
  if ((k0+l0 < n0) && (k1+l1 < n1))
    dest[(k1+l1)*n0 + (k0+l0)] = sqrtf(dest_l[l0][l1]);
}

__global__ void h1_sum_of_squares_B(int n0, int n1, 
				    float_complex *sens, float *src, 
				    int z_stride) {
  int bidxy, bidxz;
  int k0, k1, coil, l0, l1;
  __shared__ float src_l[CUDA_BLOCK_LEN][CUDA_BLOCK_LEN];
  __shared__ float_complex sens_l[CUDA_BLOCK_LEN][CUDA_BLOCK_LEN];

  // compute indices
  bidxz = floorf(((float)blockIdx.y)/((float)z_stride));
  bidxy = blockIdx.y-z_stride*bidxz;
  
  l0 = threadIdx.x;
  l1 = threadIdx.y;
  k0 = blockIdx.x*CUDA_BLOCK_LEN;
  k1 = bidxy*CUDA_BLOCK_LEN;
  coil = bidxz;
  
  // adjust pointers according to coil
  sens += coil*n0*n1;

  // read inverse image and sensitiviity
  if ((k0+l0 < n0) && (k1+l1 < n1))
    src_l[l0][l1] = 1.0f/src[(k1+l1)*n0 + (k0+l0)];
  if ( ((k0+l0/2) < n0) && (k1+l1 < n1) ) 
    sens_l[l0/2][l1][l0 & 1] = ((float *)&sens[(k1+l1)*n0 + k0])[l0];
  if ( ((k0+(l0+CUDA_BLOCK_LEN)/2) < n0) && (k1+l1 < n1) ) 
    sens_l[(l0+CUDA_BLOCK_LEN)/2][l1][l0 & 1] 
      = ((float *)&sens[(k1+l1)*n0 + k0])[l0+CUDA_BLOCK_LEN];
  __syncthreads();
  
  // perform division
  sens_l[l0][l1][0] *= src_l[l0][l1]; 
  sens_l[l0][l1][1] *= src_l[l0][l1];
  __syncthreads();

  // write result
  if ( ((k0+l0/2) < n0) && (k1+l1 < n1) ) 
    ((float *)&sens[(k1+l1)*n0 + k0])[l0] = sens_l[l0/2][l1][l0 & 1];
  if ( ((k0+(l0+CUDA_BLOCK_LEN)/2) < n0) && (k1+l1 < n1) ) 
    ((float *)&sens[(k1+l1)*n0 + k0])[l0+CUDA_BLOCK_LEN] = 
      sens_l[(l0+CUDA_BLOCK_LEN)/2][l1][l0 & 1];
}

// NFFT adjoint dual_v -> data_grad
void h1_sum_of_squares(problem *p) {
  dim3 dimBlock, dimGrid;
  int N1, N2, coils, z_stride;
   
  // init dimensions
  N1 = p->N1; N2 = p->N2; coils = p->coils;
  dimBlock.x = CUDA_BLOCK_LEN;
  dimBlock.y = CUDA_BLOCK_LEN;
  dimBlock.z = 1;
  dimGrid.x = (N1 + CUDA_BLOCK_LEN - 1)/dimBlock.x;
  dimGrid.y = (N2 + CUDA_BLOCK_LEN - 1)/dimBlock.y;
  dimGrid.z = 1;
  // multiply with sensitivities and accumulate
  h1_sum_of_squares_A<<<dimGrid, dimBlock>>>(N1, N2, coils, p->img.data.f,
					     p->sens.data.fc);
  CUT_CHECK_ERROR("h1_sum_of_squares_A() failed.");

  // init dimensions
  dimGrid.x = (N1 + CUDA_BLOCK_LEN - 1)/dimBlock.x;
  z_stride = (N2 + CUDA_BLOCK_LEN - 1)/dimBlock.y;
  dimGrid.y = coils*z_stride;
  dimGrid.z = 1;
  // multiply with sensitivities
  h1_sum_of_squares_B<<<dimGrid, dimBlock>>>(N1, N2, p->sens.data.fc, 
					     p->img.data.f, z_stride);
  CUT_CHECK_ERROR("h1_sum_of_squares_B() failed.");
}


void h1_minimize(problem *p) {
  int i;

  //h1_init_shadow_variable(p);

  for(i=0;i<p->h1_iter;i++) {
    
    printf("-----------------------\n");
    printf("starting iteration %4d\n",i);
    printf("-----------------------\n");
    
    p->h1_alpha /= p->h1_regreduction;
    printf("alpha = %g\n", p->h1_alpha);

    // update dual variables
    h1_swap_shadow_sens(p);

    printf("update dual variable v...\n");
    h1_update_dual_v(p);
    
    printf("update dual variable xi...\n");
    h1_update_xi(p);
    
    // update primal variable
    printf("update primal variable sens...\n");
    h1_update_sens(p);
    
    // update shadow variable
    printf("update shadow variable...\n");
    h1_update_shadow_sens(p);
    
  }

  for (i=0; i<p->smooth_iter; i++) {
    // smoothing step 
    h1_smooth_sens(p);
  }

  h1_sum_of_squares(p);
}


//////////////////////////////////////////////////////////////////////
// tgv minimization algorithm related functions
//////////////////////////////////////////////////////////////////////

// copies img -> shadow_img, p_x -> shadow_p_x etc.
void init_shadow_variables(problem *p) {
  // img
  CUDA_SAFE_CALL( cudaMemcpy(p->img_shadow.data.f, p->img.data.f, 
			     p->img.dim[0]*p->img.dim[1]*p->img.dim[2]
			     *sizeof(float),  
			     cudaMemcpyDeviceToDevice) ); 

  // p_x
  CUDA_SAFE_CALL( cudaMemcpy(p->p_x_shadow.data.f, p->p_x.data.f,
			     p->p_x.dim[0]*p->p_x.dim[1]*p->p_x.dim[2]
			     *sizeof(float),  
			     cudaMemcpyDeviceToDevice) ); 

  // p_y
  CUDA_SAFE_CALL( cudaMemcpy(p->p_y_shadow.data.f, p->p_y.data.f, 
			     p->p_y.dim[0]*p->p_y.dim[1]*p->p_y.dim[2]
			     *sizeof(float),  
			     cudaMemcpyDeviceToDevice) );
}

// swaps shadow image
void swap_shadow_img(problem *p) {
  float *f;
  
  // swap image
  f = p->img_shadow.data.f;
  p->img_shadow.data.f = p->img.data.f;
  p->img.data.f = f;

  // swap p_x
  f = p->p_x_shadow.data.f;
  p->p_x_shadow.data.f = p->p_x.data.f;
  p->p_x.data.f = f;
  
   // swap p_y
  f = p->p_y_shadow.data.f;
  p->p_y_shadow.data.f = p->p_y.data.f;
  p->p_y.data.f = f;
}

// updates dual variable v according to
// v^{n+1} = (v^n + sigma_v*(K img - f))/(1 + sigma)
__global__ void update_dual_v_A(float *dest, float *src0, float *src1, 
				int numel, float sigma, float sigmap1inv,
				float normestinv) {
  int i;
  
  // init id
  i = threadIdx.x + CUDA_LIN_BLOCK*blockIdx.x;
  dest += i; src0 += i; src1 += i;
  
  for (; i < numel; i += CUDA_LIN_BLOCK*CUDA_LIN_GRID) {
    // compute
    *dest = (*dest + sigma*((*src0)*normestinv - (*src1)))*sigmap1inv;
    // increment
    dest += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    src0 += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    src1 += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
  }
}

void update_dual_v(problem *p) {
  dim3 dimBlock, dimGrid;
  int numel;
  float sigma, normestinv;
  
  sigma = p->sigma;
  normestinv = 1.0f/p->normest;

  // compute transform
  transform(p);
  
  // compute update
  numel = 2*p->dual_v.dim[0]*p->dual_v.dim[1]*p->dual_v.dim[2];
  // init dimensions
  dimBlock.x = CUDA_LIN_BLOCK; dimBlock.y = 1; dimBlock.z = 1;
  dimGrid.x = CUDA_LIN_GRID; dimGrid.y = 1; dimGrid.z = 1;
  update_dual_v_A<<<dimGrid, dimBlock>>>(p->dual_v.data.f, 
					 p->nfft_dest.data.f, 
					 p->data.data.f, numel, sigma,
					 1.0f/(sigma + 1.0f), normestinv);
  CUT_CHECK_ERROR("update_dual_v_A() failed.");
}

// updates dual variable xi according to
// xi^{n+1} = P(\xi^n + sigma_xi*grad img)
// P is the projection on ||*||_\infty \leq alpha1

__global__ void update_xi_A(int N1, int N2, 
			    float *destx, float *desty, 
			    float *src, float sigma) {
  int k1, k2;
  __shared__ float src_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X];
  float valx, valy;
  
  k1 = blockIdx.x*CUDA_2D_SKIP_X + threadIdx.x;
  k2 = blockIdx.y*CUDA_2D_SKIP_Y + threadIdx.y;

  // get source buffer
  if (IS_RANGE)
    L_IDX2(src_l,0,0) = IDX2(src,0,0);
  __syncthreads();

  if (IS_L_RANGE) {
    // derivative w.r.t x
    valx = (k1 < N1-1) ? sigma*(L_IDX2(src_l,1,0) - L_IDX2(src_l,0,0)) : 0.0f; 
    
    // derivative w.r.t y
    valy = (k2 < N2-1) ? sigma*(L_IDX2(src_l,0,1) - L_IDX2(src_l,0,0)) : 0.0f; 
    
    if (IS_RANGE) {
      IDX2(destx,0,0) += valx;
      IDX2(desty,0,0) += valy;
    }
  }
}

__global__ void update_xi_B(float *destx, float *desty, 
			    float *srcx, float *srcy, 
			    int numel, float sigma, float alpha) {
  int i, tid;
  __shared__ float valx[CUDA_LIN_BLOCK], valy[CUDA_LIN_BLOCK];
  float abs;
  
  // init id
  tid = threadIdx.x;
  i = tid + CUDA_LIN_BLOCK*blockIdx.x;
  destx += i; desty += i; 
  srcx += i; srcy += i;
  
  for (; i < numel; i += CUDA_LIN_BLOCK*CUDA_LIN_GRID) {
    // compute
    valx[tid] = (*destx) - sigma*(*srcx);
    valy[tid] = (*desty) - sigma*(*srcy);
    // project
    abs = sqrtf(valx[tid]*valx[tid] + valy[tid]*valy[tid]);
    if (abs > alpha) {
      valx[tid] *= alpha/abs;
      valy[tid] *= alpha/abs;
    }
    *destx = valx[tid];
    *desty = valy[tid];
    
    // increment
    destx += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    desty += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    srcx += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    srcy += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
  }
}

void update_xi(problem *p) {
  dim3 dimBlock, dimGrid;
  int numel;
  float sigma_xi, alpha1;
  
  sigma_xi = p->sigma;
  alpha1 = p->alpha1;

  printf("alpha1 = %g \n", alpha1);

  // xi = xi + sigma*\grad img
  dimBlock.x = CUDA_2D_BLOCK_X;
  dimBlock.y = CUDA_2D_BLOCK_Y;
  dimBlock.z = 1;
  dimGrid.x = (p->N1 + CUDA_2D_SKIP_X - 1)/CUDA_2D_SKIP_X;
  dimGrid.y = (p->N2 +  CUDA_2D_SKIP_Y - 1)/CUDA_2D_SKIP_Y;
  dimGrid.z = 1;
  update_xi_A<<<dimGrid, dimBlock>>>(p->N1, p->N2, 
				     p->xi_x.data.f, p->xi_y.data.f,
				     p->img.data.f, sigma_xi);
  CUT_CHECK_ERROR("update_xi_A() failed.");

  // xi = P_{||*||_\infty <= alpha}(xi - sigma*p)
  numel = p->xi_x.dim[0]*p->xi_x.dim[1]*p->xi_x.dim[2];
  dimBlock.x = CUDA_LIN_BLOCK; dimBlock.y = 1; dimBlock.z = 1;
  dimGrid.x = CUDA_LIN_GRID; dimGrid.y = 1; dimGrid.z = 1;
  update_xi_B<<<dimGrid, dimBlock>>>(p->xi_x.data.f, p->xi_y.data.f,
				     p->p_x.data.f, p->p_y.data.f,
				     numel, sigma_xi, alpha1);
  CUT_CHECK_ERROR("update_xi_B() failed.");
}

// update eta according to
// eta^{n+1} = P(\eta^n + sigma_eta*symgrad p)
// P is the projection on ||*||_\infty \leq alpha0
__global__ void update_eta_A(int N1, int N2, 
			     float *exx, float *exy, float *eyy, 
			     float *srcx, float *srcy,
			     float sigma) {
  int k1, k2;
  __shared__ float 
    srcx_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X],
    srcy_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X];
  float vxx, vxy, vyy;
  
  k1 = blockIdx.x*CUDA_2D_SKIP_X + threadIdx.x;
  k2 = blockIdx.y*CUDA_2D_SKIP_Y + threadIdx.y;
  
  // get source buffer
  if (IS_RANGE2(-16,-1)) {
    L_IDX2(srcx_l,0,0) = IDX2(srcx,-16,-1);
    L_IDX2(srcy_l,0,0) = IDX2(srcy,-16,-1);
  }
  __syncthreads();
  
  if (IS_L_RANGE) {
    // derivative w.r.t x
    vxx = (k1 > 0) ? (L_IDX2(srcx_l,16,1) - L_IDX2(srcx_l,15,1)) : 0.0f;
    vxy = (k1 > 0) ? 0.5f*(L_IDX2(srcy_l,16,1) - L_IDX2(srcy_l,15,1)) : 0.0f;
    // derivative w.r.t y
    vxy += (k2 > 0) ? 0.5f*(L_IDX2(srcx_l,16,1) - L_IDX2(srcx_l,16,0)) : 0.0f;
    vyy = (k2 > 0) ? (L_IDX2(srcy_l,16,1) - L_IDX2(srcy_l,16,0)) : 0.0f;
    
    // write out
    if (IS_RANGE) {
      IDX2(exx,0,0) += sigma*vxx;
      IDX2(exy,0,0) += sigma*vxy;
      IDX2(eyy,0,0) += sigma*vyy;
    }
  }
}

__global__ void update_eta_B(float *exx, float *exy, float *eyy,  
			     int numel, float alpha) {
  int i, tid;
  float vxx, vxy, vyy;
  float abs;
  
  // init id
  tid = threadIdx.x;
  i = tid + CUDA_LIN_BLOCK*blockIdx.x;
  exx += i; exy += i; eyy += i; 
  
  for (; i < numel; i += CUDA_LIN_BLOCK*CUDA_LIN_GRID) {
    // project
    vxx = (*exx); abs = vxx*vxx;
    vxy = (*exy); abs += 2.0f*vxy*vxy;
    vyy = (*eyy); abs += vyy*vyy;

    abs = sqrtf(abs);
    if (abs > alpha) {
      vxx *= alpha/abs; (*exx) = vxx;
      vxy *= alpha/abs; (*exy) = vxy;
      vyy *= alpha/abs; (*eyy) = vyy;
    }
    
    // increment
    exx += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    exy += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    eyy += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
  }
}

void update_eta(problem *p) {
  dim3 dimBlock, dimGrid;
  int numel;
  float sigma_eta, alpha0;
  
  sigma_eta = p->sigma;
  alpha0 = p->alpha0;

  // eta = eta + sigma_eta*symgrad(p)
  dimBlock.x = CUDA_2D_BLOCK_X;
  dimBlock.y = CUDA_2D_BLOCK_Y;
  dimBlock.z = 1;
  dimGrid.x = (p->N1 + CUDA_2D_SKIP_X - 1)/CUDA_2D_SKIP_X;
  dimGrid.y = (p->N2 + CUDA_2D_SKIP_Y - 1)/CUDA_2D_SKIP_Y;
  dimGrid.z = 1;
  update_eta_A<<<dimGrid, dimBlock>>>(p->N1, p->N2, 
				      p->eta_xx.data.f, p->eta_xy.data.f, 
				      p->eta_yy.data.f,
				      p->p_x.data.f, p->p_y.data.f,
				      sigma_eta);
  CUT_CHECK_ERROR("update_eta_A() failed.");

  // eta = P_{||*||_\infty <= alpha}(eta)
  numel = p->eta_xx.dim[0]*p->eta_xx.dim[1]*p->eta_xx.dim[2];
  dimBlock.x = CUDA_LIN_BLOCK; dimBlock.y = 1; dimBlock.z = 1;
  dimGrid.x = CUDA_LIN_GRID; dimGrid.y = 1; dimGrid.z = 1;
  update_eta_B<<<dimGrid, dimBlock>>>(p->eta_xx.data.f, p->eta_xy.data.f, 
				      p->eta_yy.data.f, numel, alpha0);
  CUT_CHECK_ERROR("update_eta_B() failed.");
}

// update img according to
// img^{n+1} = shadow_img^{n} + tau*div xi^{n+1} - tau*K^*v^{n+1}
// assumes that K^*v^{n+1} is in data_grad^{n+1}
__global__ void update_img_A(int N1, int N2, 
			     float *dest, float *srcimg, float *srcv,
			     float *srcx, float *srcy, 
			     float tau, float normestinv) {
  int k1, k2;
  __shared__ float
    srcimg_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X],
    srcv_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X],
    srcx_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X],
    srcy_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X];
  float val;
  
  k1 = blockIdx.x*CUDA_2D_SKIP_X + threadIdx.x;
  k2 = blockIdx.y*CUDA_2D_SKIP_Y + threadIdx.y;

  // get source buffers
  if (IS_RANGE2(-16,-1)) {
    L_IDX2(srcv_l,0,0) = IDX2(srcv,-16,-1);
    L_IDX2(srcimg_l,0,0) = IDX2(srcimg,-16,-1);
    L_IDX2(srcx_l,0,0) = IDX2(srcx,-16,-1);
    L_IDX2(srcy_l,0,0) = IDX2(srcy,-16,-1);
  }
  __syncthreads();
  
  if (IS_L_RANGE) {
    // compute divergence
    val = (((float)(k1 < N1-1))*L_IDX2(srcx_l,16,1)
	   - ((float)(k1 > 0))*L_IDX2(srcx_l,15,1));
    val += (((float)(k2 < N2-1))*L_IDX2(srcy_l,16,1)
	    - ((float)(k2 > 0))*L_IDX2(srcy_l,16,0));
    
    // dest = img + tau*(div - v)
    val = fmaxf(0.0f, L_IDX2(srcimg_l,16,1) 
		+ tau*(val - normestinv*L_IDX2(srcv_l,16,1)));
    
    if (IS_RANGE)
      IDX2(dest,0,0) = val;
  }
}

void update_img(problem *p) {
  dim3 dimBlock, dimGrid;
  float tau, normestinv;
  
  tau = p->tau;  
  normestinv = 1.0f/p->normest;
  
  // compute adjoint K^*v
  transform_adjoint(p);
  
  // img = shadow_img + tau*(div xi - K^*v)
  dimBlock.x = CUDA_2D_BLOCK_X;
  dimBlock.y = CUDA_2D_BLOCK_Y;
  dimBlock.z = 1;
  dimGrid.x = (p->N1 + CUDA_2D_SKIP_X - 1)/CUDA_2D_SKIP_X;
  dimGrid.y = (p->N2 + CUDA_2D_SKIP_Y - 1)/CUDA_2D_SKIP_Y;
  dimGrid.z = 1;
  update_img_A<<<dimGrid, dimBlock>>>(p->N1, p->N2,
				      p->img.data.f, p->img_shadow.data.f,
				      p->data_grad.data.f, p->xi_x.data.f,
				      p->xi_y.data.f, tau, normestinv);
  CUT_CHECK_ERROR("update_img_A() failed.");
}

// update p according to
// p^{n+1} = p_shadow^n + tau_xi*xi^{n+1} + tau_eta*div eta^{n+1}

// divergence into dest
__global__ void update_p_A(int N1, int N2, 
			   float *destx, float *desty,
			   float *exx, float *exy, float *eyy) {
  int k1, k2;
  __shared__ float 
    exx_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X],
    exy_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X],
    eyy_l[CUDA_2D_BLOCK_Y][CUDA_2D_BLOCK_X];
  float valx, valy;

  k1 = blockIdx.x*CUDA_2D_SKIP_X + threadIdx.x;
  k2 = blockIdx.y*CUDA_2D_SKIP_Y + threadIdx.y;
  
  // get source buffer
  if (IS_RANGE) {
    L_IDX2(exx_l,0,0) = IDX2(exx,0,0);
    L_IDX2(exy_l,0,0) = IDX2(exy,0,0);
    L_IDX2(eyy_l,0,0) = IDX2(eyy,0,0);
  }
  __syncthreads();
  
  if (IS_L_RANGE) {
    // derivative w.r.t. x
    valx = (((float)(k1 < N1-1))*L_IDX2(exx_l,1,0) 
	    -((float)(k1 > 0))*L_IDX2(exx_l,0,0));
    valy = (((float)(k1 < N1-1))*L_IDX2(exy_l,1,0)
	    -((float)(k1 > 0))*L_IDX2(exy_l,0,0));
    
    // derivative w.r.t. y
    valx += (((float)(k2 < N2-1))*L_IDX2(exy_l,0,1)
	     -((float)(k2 > 0))*L_IDX2(exy_l,0,0));
    valy += (((float)(k2 < N2-1))*L_IDX2(eyy_l,0,1)
	     -((float)(k2 > 0))*L_IDX2(eyy_l,0,0));

    // write out
    if (IS_RANGE) {
      IDX2(destx,0,0) = valx;
      IDX2(desty,0,0) = valy;
    }
  }
}

// dest = src0 + tau*(src1 + dest)
__global__ void update_p_B(float *destx, float *desty,
			   float *src0x, float *src0y,  
			   float *src1x, float *src1y,  
			   int numel, float tau) {
  int i, tid;
  float valx, valy; 
  
  // init id
  tid = threadIdx.x;
  i = tid + CUDA_LIN_BLOCK*blockIdx.x;

  destx += i; desty += i; 
  src0x += i; src0y += i; 
  src1x += i; src1y += i; 
  
  for (; i < numel; i += CUDA_LIN_BLOCK*CUDA_LIN_GRID) {
    // compute
    valx = (*src0x) + tau*((*src1x) + (*destx));
    valy = (*src0y) + tau*((*src1y) + (*desty));
   
    *destx = valx;
    *desty = valy;
    
    // increment
    destx += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    desty += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    src0x += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    src0y += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    src1x += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    src1y += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
  }
}

void update_p(problem *p) {
  dim3 dimBlock, dimGrid;
  int numel;
  float tau;

  tau = p->tau;

  // p = div eta
  dimBlock.x = CUDA_2D_BLOCK_X;
  dimBlock.y = CUDA_2D_BLOCK_Y;
  dimBlock.z = 1;
  dimGrid.x = (p->N1 + CUDA_2D_SKIP_X - 1)/CUDA_2D_SKIP_X;
  dimGrid.y = (p->N2 + CUDA_2D_SKIP_Y - 1)/CUDA_2D_SKIP_Y;
  dimGrid.z = 1;
  update_p_A<<<dimGrid, dimBlock>>>(p->N1, p->N2, 
				    p->p_x.data.f, p->p_y.data.f, 
				    p->eta_xx.data.f, p->eta_xy.data.f,  
				    p->eta_yy.data.f);
  CUT_CHECK_ERROR("update_p_A() failed.");

  // p =  p_shadow + tau*(xi + p)
  numel = p->eta_xx.dim[0]*p->eta_xx.dim[1]*p->eta_xx.dim[2];
  dimBlock.x = CUDA_LIN_BLOCK; dimBlock.y = 1; dimBlock.z = 1;
  dimGrid.x = CUDA_LIN_GRID; dimGrid.y = 1; dimGrid.z = 1;
  update_p_B<<<dimGrid, dimBlock>>>(p->p_x.data.f, p->p_y.data.f, 
				    p->p_x_shadow.data.f, p->p_y_shadow.data.f, 
				    p->xi_x.data.f, p->xi_y.data.f,
				    numel, tau); 
  CUT_CHECK_ERROR("update_p_B() failed.");
}

// update shadow_img^{n+1} = 2*img^{n+1} - shadow_img^n
__global__ void update_shadow_img_A(float *dest, float *src, int numel) {
  int i;
  
  // init id
  i = threadIdx.x + CUDA_LIN_BLOCK*blockIdx.x;
  dest += i; src += i;
  
  for (; i < numel; i += CUDA_LIN_BLOCK*CUDA_LIN_GRID) {
    // compute
    *dest = 2.0f*(*src) - (*dest);
    // increment
    dest += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
    src += CUDA_LIN_BLOCK*CUDA_LIN_GRID;
  }
}

void update_shadow_img(problem *p) {
  dim3 dimBlock, dimGrid;
  int numel;

  // img part
  numel = p->img.dim[0]*p->img.dim[1]*p->img.dim[2];
  // init dimensions
  dimBlock.x = CUDA_LIN_BLOCK; dimBlock.y = 1; dimBlock.z = 1;
  dimGrid.x = CUDA_LIN_GRID; dimGrid.y = 1; dimGrid.z = 1;
  update_shadow_img_A<<<dimGrid, dimBlock>>>(p->img_shadow.data.f, 
					     p->img.data.f, numel);
  CUT_CHECK_ERROR("update_shadow_img_A() (1) failed.");

  // p_x part
  // init dimensions
  dimBlock.x = CUDA_LIN_BLOCK; dimBlock.y = 1; dimBlock.z = 1;
  dimGrid.x = CUDA_LIN_GRID; dimGrid.y = 1; dimGrid.z = 1;
  update_shadow_img_A<<<dimGrid, dimBlock>>>(p->p_x_shadow.data.f, 
					     p->p_x.data.f, numel);
  CUT_CHECK_ERROR("update_shadow_img_A() (2) failed.");

  // p_y part
  // init dimensions
  dimBlock.x = CUDA_LIN_BLOCK; dimBlock.y = 1; dimBlock.z = 1;
  dimGrid.x = CUDA_LIN_GRID; dimGrid.y = 1; dimGrid.z = 1;
  update_shadow_img_A<<<dimGrid, dimBlock>>>(p->p_y_shadow.data.f, 
					     p->p_y.data.f, numel);
  CUT_CHECK_ERROR("update_shadow_img_A() (3) failed.");
}

void tgv_minimize(problem *p) {
  int i;

  init_shadow_variables(p);

  for(i=0;i<p->tgv_iter;i++) {
    
    printf("-----------------------\n");
    printf("starting iteration %4d\n",i);
    printf("-----------------------\n");
    
    p->alpha0 /= p->regreduction;
    p->alpha1 /= p->regreduction;
    
    // update dual variables
    swap_shadow_img(p);
    printf("update dual variable v...\n");
    update_dual_v(p);
    
    printf("update dual variable xi...\n");
    update_xi(p);
    printf("update dual variable eta...\n");
    update_eta(p);
    
    // update primal variable
    printf("update primal variable img...\n");
    update_img(p);
    printf("update primal variable p...\n");
    update_p(p);
    
    // update shadow variable
    printf("update shadow variable...\n");
    update_shadow_img(p);
  }
}

int main (int argc, char *argv[]) {
  problem p;
  int pthread_result;
  char *param_file;
  double ta, tb;
     
  // no buffering for stdout
  setvbuf(stdout, NULL, _IONBF, 0);
  
  // init device
  //TODO quick fix
  //CUT_DEVICE_INIT(argc, argv);
  int dev = 0;
  CUDA_SAFE_CALL(cudaSetDevice(dev));

  param_file = 0;
  //cutGetCmdLineArgumentstr(argc, (const char **)argv, "param", &param_file);
  //TODO quick fix
  param_file = (char *)argv[1];
  param_file = (char *)&param_file[7];
  printf("debug param file name: %s\n",param_file);

  // get parameter file name
  if (!param_file) {
    printf("please specify a parameter file via -param. exiting.\n");
    return(255);
  }

  init_problem(&p, param_file);
  read_problem_data(&p);
  init_nfft_plan(&p);
  
  // start iteration
  printf("starting computations...\n");
  ta = clock();
  h1_minimize(&p);

  CUDA_SAFE_CALL( cudaMemset( p.img.data.f, 0, sizeof(float)*
			      p.img.dim[0]*p.img.dim[1]*p.img.dim[2]) );

  tgv_minimize(&p);

  // write the result back out
  device_writer(p.result_file, &p.img);

  // wait
  tb = clock();
  printf("computations finished (time = %g seconds).\n",
	 (tb - ta)/CLOCKS_PER_SEC);

  return(0);
}
