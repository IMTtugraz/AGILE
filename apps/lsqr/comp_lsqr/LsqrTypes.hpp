#ifndef __W_TYPES_H__
#define __W_TYPES_H__
#include <stdio.h>

// integer quantities are either signed or unsigned, and may have
// 8, 16, 32, or 64 bits.  Define the ones you need here.

typedef signed char        int8;
typedef unsigned char      uint8;

typedef short              int16;
typedef unsigned short     uint16;

typedef int                int32;
typedef unsigned int       uint32;

typedef long long          long64;
typedef unsigned long long ulong64;
#ifndef __64BIT__
typedef long64             int64;
#endif
typedef ulong64            uint64;


/*------------------*/
/* Vector type definitions */
/*------------------*/

typedef struct LONG_VECTOR {
  long     length;
  long     *elements;
} lvec;
typedef struct FLOAT_VECTOR {
  long     length;
  float   *elements;
} fvec;
typedef struct LSQR_INPUTS {
  long     num_rows;
  long     num_cols;
  float    damp_val;
  float    rel_mat_err;
  float    rel_rhs_err;
  float    cond_lim;
  long     max_iter;
  FILE     *lsqr_fp_out;
  fvec     *rhs_vec;
  fvec     *sol_vec;
} lsqr_input;

typedef struct LSQR_OUTPUTS {
  long     term_flag;
  long     num_iters;
  float    frob_mat_norm;
  float    mat_cond_num;
  float    resid_norm;
  float    mat_resid_norm;
  float    sol_norm;
  fvec     *sol_vec;
  fvec     *std_err_vec;
} lsqr_output;

typedef struct LSQR_WORK {
  fvec     *bidiag_wrk_vec;
  fvec     *srch_dir_vec;
} lsqr_work;

#endif
