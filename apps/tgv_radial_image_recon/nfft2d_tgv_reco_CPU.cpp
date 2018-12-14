#include <stdio.h>
#include <poll.h>
#include <nfft2d.h>
#include <ParseParam.h>

#include <cstring>
#include <cmath>
#include <cassert>

#define GFX_OUTPUT 0

#if GFX_OUTPUT
#include <pthread.h>
#include <X11/Xlib.h>
#include <gtkimageviewer-2.0/gtk-image-viewer.h>
#include <gdk/gdkx.h>
#include <gdk/gdkkeysyms.h>
#endif

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
  // display variables
  array display_img;
  array display_array;
  int display_active, redraw, curframe, display_frames, display_complex;
#if GFX_OUTPUT
  pthread_t display_thread;
  GtkWidget *gtk_window, *gtk_image_viewer;
  GdkPixbuf *gdk_pixbuf;
#endif
} problem;

void display_activate(problem *);

// allocates host memory for a 3D array
void new_array(array *ary, unsigned int N1, unsigned int N2,
               unsigned int N3, size_t size)
{
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
                      unsigned int N3, size_t size)
{
    assert(ary != NULL);

    // save dimensions
    ary->dim[0] = N1;
    ary->dim[1] = N2;
    ary->dim[2] = N3;

    // allocate space
    ary->data.f = (float*)malloc(size * N1 * N2 * N3);
    memset(ary->data.f, 0, size * N1 * N2 * N3);
}

void device_dump_array_fc(array *ary)
{
    size_t size
            = sizeof(float_complex) * ary->dim[0] * ary->dim[1] * ary->dim[2];
    float_complex *hostmem;
    hostmem = (float_complex *)malloc(size);
    memcpy(hostmem, ary->data.fc, size);
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
    memcpy(data, host_data, size);

    nfft_free(host_data);
}

// writes float_complex array into a file
void device_writer(char *filename, array *ary)
{
    FILE *file;
    array host_ary;
    unsigned long check;
    int i;

    assert(ary != NULL);

    // allocate host memory and copy data
    new_array(&host_ary, ary->dim[0], ary->dim[1], ary->dim[2],
              sizeof(float));
    memcpy(host_ary.data.f, ary->data.f,
           ary->dim[0] * ary->dim[1] * ary->dim[2] * sizeof(float));

    // open the file
    file = fopen (filename, "wb");
    assert (file != NULL);

    // write the data
    for (i=0;i<ary->dim[2];i++)
    {
        check = fwrite(&getf(host_ary,0,0,i),
                       ary->dim[0] * ary->dim[1],
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

#if GFX_OUTPUT
  // display variables
  new_array(&p->display_img, p->N1, p->N2, p->coils, sizeof(float_complex));
  p->display_array.data.fc = 0;
  p->gdk_pixbuf = 0;
  p->gtk_image_viewer = 0;
  p->redraw = FALSE;
#endif
}

void mult_data_weights(array *, array *);

void read_problem_data(problem *p) {
  assert(p != NULL);

  device_reader(p->data.data.fc, p->data_file, sizefc(p->data));
  device_reader(p->coords.data.f, p->coord_file, sizef(p->coords));
  device_reader(p->weights.data.f, p->weight_file, sizef(p->weights));
  mult_data_weights(&p->data, &p->weights);
}

void update_display_img(problem *p, array *src, int complex)
{
    // copy frame into host memory
    memcpy(p->display_img.data.f, src->data.f,
           src->dim[0] * src->dim[1] * src->dim[2] * sizeof(float)
           * (complex ? 2 : 1));
    p->display_complex = complex;
    p->display_frames = src->dim[2];
    if (p->curframe >= src->dim[2])
        p->curframe = src->dim[2] - 1;
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

float normsqr_f(array *ary)
{
    int num_elements = ary->dim[0] * ary->dim[1] * ary->dim[2];
    float* ptr = ary->data.f;
    float sum = 0;
    for (int counter = 0; counter < num_elements; ++counter)
    {
        sum += (*ptr) * (*ptr);
        ++ptr;
    }

    return sum;
}

float normsqr_fc(array *ary)
{
    int num_elements = 2 * ary->dim[0] * ary->dim[1] * ary->dim[2];
    float* ptr = ary->data.f;
    float sum = 0;
    for (int counter = 0; counter < num_elements; ++counter)
    {
        sum += (*ptr) * (*ptr);
        ++ptr;
    }

    return sum;
}

/////////////////////////////////////////////////
// transformation-related routines
/////////////////////////////////////////////////

void mult_data_weights(array *fc, array *f)
{
    int num_nodes = fc->dim[0];
    int num_coils = fc->dim[1];

    float* data = fc->data.f;
    for (int coil_counter = 0; coil_counter < num_coils; ++coil_counter)
    {
        float* w_ptr = f->data.f;
        for (int node_counter = 0; node_counter < num_nodes; ++node_counter)
        {
            // Note that data is complex, i.e. we have to multiply twice.
            (*data) *= (*w_ptr);
            ++data;
            (*data) *= (*w_ptr);
            ++data;
            ++w_ptr;
        }
    }
}

void transform(problem *p)
{
    int N1, N2, num_coils;
    float srcnorm, destnorm;

    // get source norm
    srcnorm = sqrtf(normsqr_f(&p->img));

    // init dimensions
    N1 = p->N1; N2 = p->N2; num_coils = p->coils;

    // multiply with sensitivities
    int img_size = N1 * N2;
    float* dest = p->nfft_src.data.f; // complex
    float* src = p->img.data.f;
    float* sens = p->sens.data.f; // complex
    for (int coil_counter = 0; coil_counter < num_coils; ++coil_counter)
    {
        float* src_ptr = src;
        for (int counter = 0; counter < img_size; ++counter)
        {
            // Note that dest and sens are complex-valued.
            *dest = (*sens) * (*src_ptr);
            ++dest;
            ++sens;
            *dest = (*sens) * (*src_ptr);
            ++dest;
            ++sens;
            ++src_ptr;
        }
    }

    // perform transform
    nfft_trafo_2d(&p->plan);

    // multiply with weights
    mult_data_weights(&p->nfft_dest, &p->weights);

    // get destination norm
    destnorm = sqrtf(normsqr_fc(&p->nfft_dest));

    // adjust estimate if necessary
    if ((srcnorm > 0) && (destnorm > p->normest*srcnorm))
    {
        printf("srcnorm = %g, destnorm = %g\n", srcnorm, destnorm);
        p->normest = destnorm/srcnorm;
        printf("WARNING: adjusting norm estimate to %g.\n", p->normest);
    }
}

// NFFT adjoint dual_v -> data_grad
void transform_adjoint(problem *p)
{
    int N1, N2, num_coils;
    float srcnorm, destnorm;

    // get source norm
    srcnorm = sqrtf(normsqr_fc(&p->dual_v));

    // copy into nfft_dest
    memcpy(p->nfft_dest.data.fc, p->dual_v.data.fc,
           p->dual_v.dim[0] * p->dual_v.dim[1] * p->dual_v.dim[2]
           * sizeof(float_complex));

    // multiply with weights
    mult_data_weights(&p->nfft_dest, &p->weights);

    // perform transform adjoint
    nfft_adjoint_2d(&p->plan);

    // init dimensions
    N1 = p->N1; N2 = p->N2; num_coils = p->coils;

    // multiply with sensitivities and accumulate
    int img_size = N1 * N2;
    float* dest = p->data_grad.data.f;
    float* src = p->nfft_src.data.f; // complex
    float* sens = p->sens.data.f; // complex
    for (int counter = 0; counter < img_size; ++counter)
    {
        float sum = 0;
        float* src_ptr = src;
        float* sens_ptr = sens;
        for (int coil_counter = 0; coil_counter < num_coils; ++coil_counter)
        {
            sum += (*src_ptr) * (*sens_ptr)
                   + (*(src_ptr + 1)) * (*(sens_ptr + 1));
            // Go to the next coil.
            src_ptr += 2 * img_size;
            sens_ptr += 2 * img_size;
        }
        *dest = sum;

        // Advance to the next node.
        src += 2;
        sens += 2;
        ++dest;
    }

    // get destination norm
    destnorm = sqrtf(normsqr_f(&p->data_grad));

    // adjust estimate if necessary
    if ((srcnorm > 0) && (destnorm > p->normest*srcnorm))
    {
        printf("srcnorm = %g, destnorm = %g\n", srcnorm, destnorm);
        p->normest = destnorm/srcnorm;
        printf("WARNING: adjusting norm estimate to %g.\n", p->normest);
    }
}

void h1_transform(problem *p)
{
    float srcnorm, destnorm;

    // get source norm
    srcnorm = sqrtf(normsqr_fc(&p->sens));

    // copy sens -> nfft_src
    memcpy(p->nfft_src.data.fc, p->sens.data.fc,
           p->sens.dim[0] * p->sens.dim[1] * p->sens.dim[2]
           * sizeof(float_complex));

    // perform transform
    nfft_trafo_2d(&p->plan);

    // multiply with weights
    mult_data_weights(&p->nfft_dest, &p->weights);

    // get destination norm
    destnorm = sqrtf(normsqr_fc(&p->nfft_dest));

    // adjust estimate if necessary
    if ((srcnorm > 0) && (destnorm > p->h1_normest*srcnorm))
    {
        printf("srcnorm = %g, destnorm = %g\n", srcnorm, destnorm);
        p->h1_normest = destnorm/srcnorm;
        printf("WARNING: adjusting norm estimate to %g.\n", p->h1_normest);
    }
}

void h1_transform_adjoint(problem *p)
{
    float srcnorm, destnorm;

    // get source norm
    srcnorm = sqrtf(normsqr_fc(&p->dual_v));

    // copy into nfft_dest
    memcpy(p->nfft_dest.data.fc, p->dual_v.data.fc,
           p->dual_v.dim[0] * p->dual_v.dim[1] * p->dual_v.dim[2]
           * sizeof(float_complex));

    // multiply with weights
    mult_data_weights(&p->nfft_dest, &p->weights);

    // perform transform adjoint
    nfft_adjoint_2d(&p->plan);

    // get destination norm
    destnorm = sqrtf(normsqr_fc(&p->nfft_src));

    // adjust estimate if necessary
    if ((srcnorm > 0) && (destnorm > p->h1_normest*srcnorm))
    {
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

void h1_update_dual_v(problem *p)
{
    float sigma = p->h1_sigma;
    float normestinv = 1.0f/p->h1_normest;

    // compute transform
    h1_transform(p);

    // compute update according to
    // v^{n+1} = (v^n + sigma_v*(K img - f))/(1 + sigma)
    int numel = 2 * p->dual_v.dim[0] * p->dual_v.dim[1] * p->dual_v.dim[2];
    float* dest = p->dual_v.data.f;
    float* src0 = p->nfft_dest.data.f;
    float* src1 = p->data.data.f;
    float sigmap1inv = 1.0f / (sigma + 1.0f);
    for (int counter = 0; counter < numel; ++counter)
    {
        *dest = (*dest + sigma * ((*src0) * normestinv - (*src1))) * sigmap1inv;
        ++dest;
        ++src0;
        ++src1;
    }
}

void h1_update_xi(problem *p)
{
    float sigma = p->h1_sigma;
    float alpha = p->h1_alpha;
    float sigma_alphap1inv = 1.0f/(1.0f + sigma/alpha);

    int x_size = p->N1 * 2; // Note: doubled because the image is complex
    int img_size = x_size * p->N2;
    for (int coil_counter = 0; coil_counter < p->coils; ++coil_counter)
    {
        int coil_offset = img_size * coil_counter;

        // Position input pointer at first pixel of current sensitivity.
        const float* img_ptr = p->sens.data.f + coil_offset;
        // Position output pointers at first pixel of the current sensitivity.
        float* destx = p->h1_xi_x.data.f + coil_offset;
        float* desty = p->h1_xi_y.data.f + coil_offset;

        // Loop over all rows except the last one.
        for (int y = 0; y < p->N2 - 1; ++y)
        {
            // Loop over all columns except the last one.
            for (int x = 0; x < x_size - 2; ++x)
            {
                // Note: stride of 2 because data is complex
                float x_derivative = *(img_ptr + 2) - *img_ptr;
                *destx = (*destx + sigma * x_derivative) * sigma_alphap1inv;

                float y_derivative = *(img_ptr + x_size) - *img_ptr;
                *desty = (*desty + sigma * y_derivative) * sigma_alphap1inv;

                ++destx;
                ++desty;
                ++img_ptr;
            }
            // No x-derivative in the last column.
            destx += 2;
            // y-derivative in the last column.
            float y_derivative = *(img_ptr + x_size) - *img_ptr;
            *desty = (*desty + sigma * y_derivative) * sigma_alphap1inv;
            ++desty;
            ++img_ptr;
            y_derivative = *(img_ptr + x_size) - *img_ptr;
            *desty = (*desty + sigma * y_derivative) * sigma_alphap1inv;
            ++desty;
            ++img_ptr;
        }

        // x-derivative in the last row.
        for (int x = 0; x < x_size - 2; ++x)
        {
            // Note: stride of 2 because data is complex
            float x_derivative = *(img_ptr + 2) - *img_ptr;
            *destx = (*destx + sigma * x_derivative) * sigma_alphap1inv;
            ++destx;
            ++img_ptr;
        }
        // No y-derivative in the last row.
    }
}

void h1_update_sens(problem* p)
{
    float tau = p->h1_tau;
    float normestinv = 1.0f/p->h1_normest;

    // compute adjoint K^*v
    h1_transform_adjoint(p);

    int x_size = p->N1 * 2; // Note: complex data
    int img_size = x_size * p->N2;
    for (int coil_counter = 0; coil_counter < p->coils; ++coil_counter)
    {
        int coil_offset = img_size * coil_counter;

        float* dest = p->sens.data.f + coil_offset;
        const float* srcimg = p->sens_shadow.data.f + coil_offset;
        const float* srcv = p->nfft_src.data.f + coil_offset;
        const float* srcx = p->h1_xi_x.data.f + coil_offset;
        const float* srcy = p->h1_xi_y.data.f + coil_offset;

        for (int y = 0; y < p->N2; ++y)
        {
            for (int x = 0; x < x_size; ++x)
            {
                // Note: Could be optimized by moving if-branches outside loops
                float divergence = *srcx;
                if (x >= 2)
                    divergence -= *(srcx - 2);

                divergence += *srcy;
                if (y >= 1)
                    divergence -= *(srcy - x_size);

                // dest = img + tau * (div - v)
                *dest = *srcimg + tau * (divergence - normestinv * *srcv);

                ++srcx;
                ++srcy;
                ++srcimg;
                ++srcv;
                ++dest;
            }
        }
    }
}

void h1_update_shadow_sens(problem *p)
{
    int numel = 2*p->sens.dim[0]*p->sens.dim[1]*p->sens.dim[2];
    float* dest = p->sens_shadow.data.f;
    float* src = p->sens.data.f;
    for (int counter = 0; counter < numel; ++counter)
    {
        *dest = 2.0f * (*src) - (*dest);
        ++dest;
        ++src;
    }
}

void h1_smooth_sens(problem *p)
{
    int x_size = p->N1 * 2; // Note: complex data
    int img_size = x_size * p->N2;

    // ----- Step 1 -----
    for (int coil_counter = 0; coil_counter < p->coils; ++coil_counter)
    {
        int coil_offset = img_size * coil_counter;

        float* destx = p->h1_xi_x.data.f + coil_offset;
        float* desty = p->h1_xi_y.data.f + coil_offset;
        float* src = p->sens.data.f + coil_offset;

        // Loop over all rows except the last one.
        for (int y = 0; y < p->N2 - 1; ++y)
        {
            // Loop over all columns except the last one.
            for (int x = 0; x < x_size - 2; ++x)
            {
                // Note: stride of 2 because data is complex
                float val_x = *(src + 2) + *src;
                *destx = 0.5f * val_x;

                float val_y = *(src + x_size) + *src;
                *desty = 0.5f * val_y;

                ++destx;
                ++desty;
                ++src;
            }
            // No x-sum in the last column.
            destx += 2;
            // y-sum in the last column.
            float val_y = *(src + x_size) + *src;
            *desty = 0.5f * val_y;
            ++desty;
            ++src;
            val_y = *(src + x_size) + *src;
            *desty = 0.5f * val_y;
            ++desty;
            ++src;
        }

        // x-sum in the last row.
        for (int x = 0; x < x_size - 2; ++x)
        {
            // Note: stride of 2 because data is complex
            float val_x = *(src + 2) + *src;
            *destx = 0.5f * val_x;
            ++destx;
            ++src;
        }
        // No y-sum in the last row.
    }


    // ----- Step 2 -----
    for (int coil_counter = 0; coil_counter < p->coils; ++coil_counter)
    {
        int coil_offset = img_size * coil_counter;

        float* dest = p->sens.data.f + coil_offset;
        const float* srcx = p->h1_xi_x.data.f + coil_offset;
        const float* srcy = p->h1_xi_y.data.f + coil_offset;

        for (int y = 0; y < p->N2; ++y)
        {
            for (int x = 0; x < x_size; ++x)
            {
                // Note: Could be optimized by moving if-branches outside loops
                float sum = *srcx;
                if (x >= 2)
                    sum += *(srcx - 2);

                sum += *srcy;
                if (y >= 1)
                    sum += *(srcy - x_size);

                *dest = 0.5f * sum;

                ++srcx;
                ++srcy;
                ++dest;
            }
        }
    }
}

// NFFT adjoint dual_v -> data_grad
void h1_sum_of_squares(problem *p)
{
    int img_size = p->N1 * p->N2;
    int num_coils = p->coils;

    float* dest = p->img.data.f;
    float* sens = p->sens.data.f; // complex
    for (int counter = 0; counter < img_size; ++counter)
    {
        // Compute sum of squares of sensitivities.
        float sqr_sum = 0;
        float* sens_ptr = sens;
        for (int coil_counter = 0; coil_counter < num_coils; ++coil_counter)
        {
            sqr_sum += (*sens_ptr) * (*sens_ptr)
                       + (*(sens_ptr + 1)) * (*(sens_ptr + 1));
            // Go to the next coil.
            // Note: doubled because of complex data
            sens_ptr += 2 * img_size;
        }
        *dest = sqrtf(sqr_sum);

        // Divide every sensitivity by inverse of image.
        float inv_img = 1.0f / *dest;
        sens_ptr = sens;
        for (int coil_counter = 0; coil_counter < num_coils; ++coil_counter)
        {
            *sens_ptr *= inv_img;
            *(sens_ptr + 1) *= inv_img;
            // Go to the next coil.
            sens_ptr += 2 * img_size;
        }

        // Advance to the next node.
        sens += 2;
        ++dest;
    }
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

#if GFX_OUTPUT
    update_display_img(p, &p->sens, 1);
    p->redraw = TRUE;
    display_activate(p);
#endif
  }

  for (i=0; i<p->smooth_iter; i++) {
    // smoothing step
    h1_smooth_sens(p);
  }

  h1_sum_of_squares(p);

#if GFX_OUTPUT
  update_display_img(p, &p->sens, 1);
  p->redraw = TRUE;
  display_activate(p);
#endif
}


//////////////////////////////////////////////////////////////////////
// tgv minimization algorithm related functions
//////////////////////////////////////////////////////////////////////

// copies img -> shadow_img, p_x -> shadow_p_x etc.
void init_shadow_variables(problem *p)
{
    // img
    memcpy(p->img_shadow.data.f, p->img.data.f,
           p->img.dim[0] * p->img.dim[1] * p->img.dim[2] * sizeof(float));

    // p_x
    memcpy(p->p_x_shadow.data.f, p->p_x.data.f,
           p->p_x.dim[0] * p->p_x.dim[1] * p->p_x.dim[2] * sizeof(float));

    // p_y
    memcpy(p->p_y_shadow.data.f, p->p_y.data.f,
           p->p_y.dim[0] * p->p_y.dim[1] * p->p_y.dim[2] * sizeof(float));
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

void update_dual_v(problem *p)
{
    float sigma, normestinv;

    sigma = p->sigma;
    normestinv = 1.0f/p->normest;

    // compute transform
    transform(p);

    // compute update according to
    // v^{n+1} = (v^n + sigma_v*(K img - f))/(1 + sigma)
    int numel = 2*p->dual_v.dim[0]*p->dual_v.dim[1]*p->dual_v.dim[2];
    float* dest = p->dual_v.data.f;
    float* src0 = p->nfft_dest.data.f;
    float* src1 = p->data.data.f;
    float sigmap1inv = 1.0f/(sigma + 1.0f);
    for (int counter = 0; counter < numel; ++counter)
    {
        *dest = (*dest + sigma * ((*src0) * normestinv - (*src1))) * sigmap1inv;
        ++dest;
        ++src0;
        ++src1;
    }
}

// updates dual variable xi according to
// xi^{n+1} = P(\xi^n + sigma_xi*grad img)
// P is the projection on ||*||_\infty \leq alpha1
void update_xi(problem *p)
{
    float sigma = p->sigma;

    float* destx = p->xi_x.data.f;
    float* desty = p->xi_y.data.f;
    float* img = p->img.data.f;

    // xi = xi + sigma*grad img

    // Loop over all rows except the last one.
    for (int y = 0; y < p->N2 - 1; ++y)
    {
        // Loop over all columns except the last one.
        for (int x = 0; x < p->N1 - 1; ++x)
        {
            float x_derivative = *(img + 1) - *img;
            *destx += sigma * x_derivative;

            float y_derivative = *(img + p->N1) - *img;
            *desty += sigma * y_derivative;

            ++destx;
            ++desty;
            ++img;
        }
        // No x-derivative in the last column.
        ++destx;
        // y-derivative in the last column.
        float y_derivative = *(img + p->N1) - *img;
        *desty += sigma * y_derivative;
        ++desty;
        ++img;
    }

    // x-derivative in the last row.
    for (int x = 0; x < p->N1 - 1; ++x)
    {
        float x_derivative = *(img + 1) - *img;
        *destx += sigma * x_derivative;
        ++destx;
        ++img;
    }
    // No y-derivative in the last row.

    // ----- Projection -----
    float alpha = p->alpha1;

    destx = p->xi_x.data.f;
    desty = p->xi_y.data.f;
    float* srcx = p->p_x.data.f;
    float* srcy = p->p_y.data.f;
    int numel = p->xi_x.dim[0] * p->xi_x.dim[1] * p->xi_x.dim[2];
    for (int counter = 0; counter < numel; ++counter)
    {
        float valx = (*destx) - sigma * (*srcx);
        float valy = (*desty) - sigma * (*srcy);
        // project
        float abs = sqrtf(valx * valx + valy * valy);
        if (abs > alpha)
        {
            valx *= alpha / abs;
            valy *= alpha / abs;
        }
        *destx = valx;
        *desty = valy;

        ++destx;
        ++desty;
        ++srcx;
        ++srcy;
    }
}

// update eta according to
// eta^{n+1} = P(\eta^n + sigma_eta*symgrad p)
// P is the projection on ||*||_\infty \leq alpha0
void update_eta(problem *p)
{
    float sigma = p->sigma;

    float* exx = p->eta_xx.data.f;
    float* exy = p->eta_xy.data.f;
    float* eyy = p->eta_yy.data.f;
    float* srcx = p->p_x.data.f;
    float* srcy = p->p_y.data.f;

    // ----- Symmetrized gradient -----
    // eta = eta + sigma_eta*symgrad(p)
    for (int y = 0; y < p->N2; ++y)
    {
        for (int x = 0; x < p->N1; ++x)
        {
            // Note: Could be optimized by moving if-branches outside loops

            // derivative w.r.t. x
            float vxx = 0;
            float vxy = 0;
            float vyy = 0;
            if (x > 0)
            {
                vxx =         *srcx - *(srcx - 1);
                vxy = 0.5f * (*srcy - *(srcy - 1));
            }
            // derivative w.r.t. y
            if (y > 0)
            {
                vxy += 0.5f * (*srcx - *(srcx - p->N1));
                vyy =          *srcy - *(srcy - p->N1);
            }

            *exx += sigma * vxx;
            *exy += sigma * vxy;
            *eyy += sigma * vyy;

            ++exx;
            ++exy;
            ++eyy;
            ++srcx;
            ++srcy;
        }
    }

    // ----- Projection -----
    // eta = P_{||*||_\infty <= alpha}(eta)
    float alpha = p->alpha0;

    exx = p->eta_xx.data.f;
    exy = p->eta_xy.data.f;
    eyy = p->eta_yy.data.f;

    int numel = p->eta_xx.dim[0] * p->eta_xx.dim[1] * p->eta_xx.dim[2];
    for (int counter = 0; counter < numel; ++counter)
    {
        float vxx = (*exx);
        float vxy = (*exy);
        float vyy = (*eyy);
        float abs = sqrtf(vxx * vxx + 2.0f * vxy * vxy + vyy * vyy);

        if (abs > alpha)
        {
            float t = alpha / abs;
            *exx = vxx * t;
            *exy = vxy * t;
            *eyy = vyy * t;
        }

        ++exx;
        ++exy;
        ++eyy;
    }
}

// update img according to
// img^{n+1} = shadow_img^{n} + tau*div xi^{n+1} - tau*K^*v^{n+1}
// assumes that K^*v^{n+1} is in data_grad^{n+1}
void update_img(problem *p)
{
    // compute adjoint K^*v
    transform_adjoint(p);

    float tau = p->tau;
    float normestinv = 1.0f / p->normest;

    float* srcimg = p->img_shadow.data.f;
    float* srcv = p->data_grad.data.f;
    float* srcx = p->xi_x.data.f;
    float* srcy = p->xi_y.data.f;
    float* dest = p->img.data.f;

    for (int y = 0; y < p->N2; ++y)
    {
        for (int x = 0; x < p->N1; ++x)
        {
            // Note: Could be optimized by moving if-branches outside loops
            float divergence = *srcx;
            if (x >= 1)
                divergence -= *(srcx - 1);

            divergence += *srcy;
            if (y >= 1)
                divergence -= *(srcy - p->N1);

            // dest = img + tau*(div - v)
            *dest = fmaxf(0.0f,
                          *srcimg + tau * (divergence - normestinv * (*srcv)));

            ++srcimg;
            ++srcv;
            ++srcx;
            ++srcy;
            ++dest;
        }
    }
}

// update p according to
// p^{n+1} = p_shadow^n + tau_xi*xi^{n+1} + tau_eta*div eta^{n+1}
void update_p(problem *p)
{
    // p = div eta
    float* destx = p->p_x.data.f;
    float* desty = p->p_y.data.f;
    float* exx = p->eta_xx.data.f;
    float* exy = p->eta_xy.data.f;
    float* eyy = p->eta_yy.data.f;

    for (int y = 0; y < p->N2; ++y)
    {
        for (int x = 0; x < p->N1; ++x)
        {
            // derivative w.r.t. x
            float x_derivative = 0;
            float y_derivative = 0;
            if (x < p->N1 - 1)
            {
                x_derivative += *(exx + 1);
                y_derivative += *(exy + 1);
            }
            if (x > 0)
            {
                x_derivative -= *exx;
                y_derivative -= *exy;
            }

            if (y < p->N2 - 1)
            {
                x_derivative += *(exy + p->N1);
                y_derivative += *(eyy + p->N1);
            }
            if (y > 0)
            {
                x_derivative -= *exy;
                y_derivative -= *eyy;
            }

            *destx = x_derivative;
            *desty = y_derivative;

            ++exx;
            ++exy;
            ++eyy;
            ++destx;
            ++desty;
        }
    }

    float tau = p->tau;

    // p = p_shadow + tau * (xi + p)
    int numel = p->eta_xx.dim[0] * p->eta_xx.dim[1] * p->eta_xx.dim[2];
    destx = p->p_x.data.f;
    desty = p->p_y.data.f;
    float *src_p_x = p->p_x_shadow.data.f;
    float *src_p_y = p->p_y_shadow.data.f;
    float *src_xi_x = p->xi_x.data.f;
    float *src_xi_y = p->xi_y.data.f;
    for (int counter = 0; counter < numel; ++counter)
    {
        *destx = (*src_p_x) + tau * ((*src_xi_x) + (*destx));
        *desty = (*src_p_y) + tau * ((*src_xi_y) + (*desty));

        ++src_p_x;
        ++src_p_y;
        ++src_xi_x;
        ++src_xi_y;
        ++destx;
        ++desty;
    }
}

// update shadow_img^{n+1} = 2*img^{n+1} - shadow_img^n
void update_shadow_img_helper(float* dest, float* src, int numel)
{
    for (int counter = 0; counter < numel; ++counter)
    {
        *dest = 2.0f * (*src) - (*dest);
        ++dest;
        ++src;
    }
}

void update_shadow_img(problem *p)
{
    int numel = p->img.dim[0]*p->img.dim[1]*p->img.dim[2];
    // img part
    update_shadow_img_helper(p->img_shadow.data.f, p->img.data.f, numel);
    // p_x part
    update_shadow_img_helper(p->p_x_shadow.data.f, p->p_x.data.f, numel);
    // p_y part
    update_shadow_img_helper(p->p_y_shadow.data.f, p->p_y.data.f, numel);
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

#if GFX_OUTPUT
    update_display_img(p, &p->img, 0);
    p->redraw = TRUE;
    display_activate(p);
#endif
  }
}

#if GFX_OUTPUT
// display thread

// computes, scales abs of img to range unsigned char
// and stores in ary (RGB0)
void update_display(problem *p, array *ary, int frame) {
  int numel, i;
  float curabs, maxabs;
  float *img;
  unsigned char v;
  unsigned char *display;

  if (p->display_complex) {
    // complex case
    // find maximal abs
    maxabs = 0;
    numel = p->img.dim[0]*p->img.dim[1];
    img = (float *)(p->display_img.data.fc + frame*numel);
    for(i=0;i<numel;i++) {
      curabs = sqrtf(img[0]*img[0] + img[1]*img[1]);
      if (curabs > maxabs)
        maxabs = curabs;
      img += 2;
    }

    // scale and copy
    numel = p->N1*p->N2;
    img = (float *)(p->display_img.data.fc + frame*numel);
    display = ary->data.uc;
    for(i=0;i<numel;i++) {
      v = 255.0f*sqrtf(img[0]*img[0] + img[1]*img[1])/maxabs;
      *display = v; *(display+1) = v; *(display+2) = v;
      *(display+3) = 255;
      img+=2; display+=4;
    }
  } else {
    // real case
    // find maximal abs
    maxabs = 0;
    numel = p->img.dim[0]*p->img.dim[1];
    img = p->display_img.data.f + frame*numel;
    for(i=0;i<numel;i++) {
      curabs = fabs(*img);
      if (curabs > maxabs)
        maxabs = curabs;
      img++;
    }

    // scale and copy
    numel = p->N1*p->N2;
    img = p->display_img.data.f + frame*numel;
    display = ary->data.uc;
    for(i=0;i<numel;i++) {
      v = 255.0f*fabs(*img)/maxabs;
      *display = v; *(display+1) = v; *(display+2) = v;
      *(display+3) = 255;
      img++; display+=4;
    }
  }
}

gint cb_key_press_event(GtkWidget *widget, GdkEventKey *event,
                        gpointer *func_data) {
  char titlestring[150];
  problem *p;

  p = (problem *)func_data;
  if (event->keyval == GDK_Left) {
    if (p->curframe > 0) {
      p->curframe--;
      sprintf(titlestring, "Reconstruction %d/%d",
              p->curframe+1, p->display_img.dim[2]);
      gtk_window_set_title (GTK_WINDOW (p->gtk_window), titlestring);
      p->redraw = TRUE;
    }
    return(1);
  };
  if (event->keyval == GDK_Right) {
    if (p->curframe < p->display_img.dim[2]-1) {
      p->curframe++;
      sprintf(titlestring, "Reconstrution %d/%d",
              p->curframe+1, p->display_img.dim[2]);
      gtk_window_set_title (GTK_WINDOW (p->gtk_window), titlestring);
      p->redraw = TRUE;
    }
    return(1);
  };

  return(0);
}

gboolean display_idle_func(gpointer *func_data) {
  struct timespec sleeptime={0, 50000000}, remainingtime;
  problem *p;

  p = (problem *)func_data;
  if (p->redraw) {
    update_display(p, &p->display_array, p->curframe);
    if (p->display_active)
      gtk_image_viewer_redraw(GTK_IMAGE_VIEWER(p->gtk_image_viewer), TRUE);
    p->redraw = FALSE;
  }

  nanosleep(&sleeptime, &remainingtime);
  return(1);
}

gint display_leave(GtkWidget *widget, void *dummy, gpointer *func_data) {
  problem *p;

  p = (problem *)func_data;
  gtk_widget_hide_all (p->gtk_window);
  gtk_main_quit();

  return(1);
}

void *display_thread_func(void *args) {
  problem *p;
  GdkPixbuf *pixbuf;
  GtkWidget *window, *image_viewer;
  struct timespec sleeptime={0, 50000000}, remainingtime;

  p = (problem *)args;
  printf("starting display thread...\n");
  p->display_active = 1;

  // create display data
  p->curframe = 0;
  p->redraw = TRUE;
  if (!p->display_array.data.fc) {
    new_array(&p->display_array, 4, p->N1, p->N2, sizeof(unsigned char));
    assert(p->display_array.data.fc != NULL);
  }

  // create GdkPixbuf out of display data
  if (!p->gdk_pixbuf) {
    p->gdk_pixbuf
      = gdk_pixbuf_new_from_data((guchar *)p->display_array.data.uc,
                                 GDK_COLORSPACE_RGB, TRUE,
                                 8, p->N1, p->N2, p->N1*4, NULL,
                                 NULL);
    assert(p->gdk_pixbuf != NULL);
  }
  pixbuf = p->gdk_pixbuf;

  p->gtk_window = window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  gtk_window_set_policy(GTK_WINDOW(window), TRUE, TRUE, FALSE);
  gtk_window_set_title (GTK_WINDOW (window), "Reconstruction");
  g_signal_connect (window, "key_press_event",
                    G_CALLBACK(cb_key_press_event), (gpointer)p);

  g_signal_connect (window, "delete_event", G_CALLBACK(display_leave),
                    (gpointer)p);

  p->gtk_image_viewer = gtk_image_viewer_new(pixbuf);
  gtk_container_add (GTK_CONTAINER (window), p->gtk_image_viewer);

  g_idle_add((GSourceFunc)display_idle_func, (gpointer)p);
  do {
    gtk_widget_show_all (window);
    gtk_main();

    XSetCloseDownMode (GDK_DISPLAY(), RetainPermanent);
    XCloseDisplay(GDK_DISPLAY());

    printf("display window closed.\n");
    p->display_active = 0;

    while (!p->display_active)
      nanosleep(&sleeptime, &remainingtime);

    GDK_DISPLAY() = XOpenDisplay(0);
  } while (1);

  return(NULL);
}

void display_activate(problem *p) {
  struct timespec timeout={0, 1};
  struct pollfd fds={0, POLLIN, 0};
  char s;

  if (ppoll(&fds, 1, &timeout, 0)) {
    s = fgetc(stdin);
    //printf("got %d\n", s);
    if (s == '\n')
      p->display_active = 1;
  }
}
#endif

int main (int argc, char *argv[])
{
  problem p;
  char *param_file;
#if GFX_OUTPUT
  int pthread_result;

  // init gtk
  gtk_init(&argc, &argv);
#endif

  // no buffering for stdout
  setvbuf(stdout, NULL, _IONBF, 0);

  // init device
  param_file = 0;
  for (int counter = 1; counter < argc; ++counter)
  {
      if (strcmp(argv[counter], "param") == 0 && counter < argc - 1)
      {
          param_file = argv[counter + 1];
          break;
      }
  }

  // get parameter file name
  if (!param_file) {
    printf("please specify a parameter file via -param. exiting.\n");
    return(255);
  }

  init_problem(&p, param_file);
  read_problem_data(&p);
  init_nfft_plan(&p);

#if GFX_OUTPUT
  // create display thread
  pthread_result = pthread_create(&p.display_thread, NULL,
                                  display_thread_func, (void *)&p);
  assert(pthread_result == 0);
#endif

  // start iteration
  h1_minimize(&p);

  memset(p.img.data.f, 0,
         sizeof(float) * p.img.dim[0] * p.img.dim[1] * p.img.dim[2]);

  tgv_minimize(&p);

  // write the result back out
  device_writer(p.result_file, &p.img);

  // wait
  printf("computations finished. press any key to exit.");
  int foo = getchar();

  return(0);
}
