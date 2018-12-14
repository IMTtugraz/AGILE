/*
* test_prog.c
* Tests C++ version of LSQR.
*
* 27 June 2009: Gerald Buchgraber <gerald.buchgraber@student.tugraz.at>
*/

#include <vector>

#include "agile/io/file.hpp"

#include "Model.hpp"
#include "Lsqr.hpp"

#define SHOW_CPU_LSQR_DETAILS 0

#include "../config.hpp"

int main( int argc, char* argv[])
{
  char *A_matrix_file = NULL;
  char *AT_matrix_file = NULL;
  char *vector_file = NULL;
  char *result_file = NULL;

  // input checks
  if (argc <= 4)
  {
    std::cout << "usage: " << argv[0] << " <A matrix file> <A' matrix file> "
              << "<rhs vector file> <result file>" << std::endl;
    return -1;
  } else {
    A_matrix_file = argv[1];
    AT_matrix_file = argv[2];
    vector_file = argv[3];
    result_file = argv[4];
  }

  // read in crs matrix
  //---------------------------------------------------------------------------
  unsigned A_num_rows, A_num_columns;
  std::vector<unsigned> A_row_nnz;
  std::vector<unsigned> A_column_index;
  std::vector<float> A_data;
  bool success = agile::readCSMatrixFile(A_matrix_file,
                                           A_num_rows, A_num_columns,
                                           A_row_nnz, A_column_index,
                                           A_data);
  if (!success)
  {
    std::cerr << "error: not able to load matrix A: " << A_matrix_file
              << std::endl;
    exit(-1);
  }

  // read in crs matrix A'
  //---------------------------------------------------------------------------
  unsigned AT_num_rows, AT_num_columns;
  std::vector<unsigned> AT_row_nnz;
  std::vector<unsigned> AT_column_index;
  std::vector<float> AT_data;
  success = agile::readCSMatrixFile(AT_matrix_file,
                                      AT_num_rows, AT_num_columns,
                                      AT_row_nnz, AT_column_index,
                                      AT_data);
  if (!success)
  {
    std::cerr << "error: not able to load matrix A': " << AT_matrix_file
              << std::endl;
    exit(-1);
  }

  // read in vector
  //---------------------------------------------------------------------------
  std::vector<float> y_host;
  success = agile::readVectorFile(vector_file, y_host);
  if (!success)
  {
    std::cerr << "error: not able to load vector: " << vector_file << std::endl;
    exit(-1);
  }

  // dimension check
  if (y_host.size() != A_num_rows) {
    std::cerr << "error: incompatible dimensions" << std::endl;
  }


  //---------------------------------------------------------------------------
  // INITILIZATION FROM HERE
  //---------------------------------------------------------------------------

  // Instantiate the model
  Model *model = new Model(A_num_rows, A_num_columns);
  model->A_row_nnz = &A_row_nnz;
  model->A_column_index = &A_column_index;
  model->A_data = &A_data;
  model->AT_row_nnz = &AT_row_nnz;
  model->AT_column_index = &AT_column_index;
  model->AT_data = &AT_data;

  // and then Lsqr
  Lsqr *lsqr = new Lsqr(model);

  // Allocate memory for lsqr
  lsqr->allocLsqrMem();

  // Specify output unit
#if SHOW_CPU_LSQR_DETAILS
  lsqr->input->lsqr_fp_out = stdout;
#else
  lsqr->input->lsqr_fp_out = NULL;
#endif

/*
*     Copy the right-hand side vector into the right-hand side vector for LSQR.
*/
  fvec fvec_y_host;
  fvec_y_host.length = y_host.size();
  fvec_y_host.elements = &y_host[0];
  fvec_copy( &fvec_y_host, lsqr->input->rhs_vec );

/*
*  Set the initial guess for LSQR.
*/
  for( unsigned i = 0; i < A_num_columns; ++i)
  {
    lsqr->input->sol_vec->elements[i] = 0.0;
  }

/*
*     Set the input parameters for LSQR.
*/
  lsqr->input->num_rows = A_num_rows;
  lsqr->input->num_cols = A_num_columns;
  lsqr->input->damp_val = 0.0; // @see Lsqr.hpp
  lsqr->input->rel_mat_err = LSQR_ABS_TOLERANCE;
  lsqr->input->rel_rhs_err = LSQR_ABS_TOLERANCE;
  lsqr->input->cond_lim = 100.0/* 10.0 * act_mat_cond_num*/; //TODO check this
  lsqr->input->max_iter = LSQR_MAX_ITERATIONS;

  //---------------------------------------------------------------------------
  // COMPUTATION FROM HERE
  //---------------------------------------------------------------------------


#if WITH_TIMER
  struct timeval st, et;
  gettimeofday(&st, NULL);
#endif

  lsqr->do_lsqr( model);

#if WITH_TIMER
  gettimeofday(&et, NULL);
  float elapsed_time = ((et.tv_sec-st.tv_sec)*1000.0
    + (et.tv_usec - st.tv_usec)/1000.0);
  std::cout << "lsqr (cpu):  " << std::setprecision(5.9)
            << elapsed_time << "ms" << std::endl;
#endif

#if SHOW_LSQR_DETAILS
  std::cout << "iterations: " << lsqr->output->num_iters << std::endl;
  std::cout << "final residual: " << lsqr->output->resid_norm << std::endl;
#endif

  // transfer result from gpu to cpu
  std::vector<float> x_host(A_num_columns);
  //x_host.assign(A_num_columns, *(lsqr->output->sol_vec->elements));
  for (unsigned i=0; i<A_num_columns; ++i)
    x_host[i] = lsqr->output->sol_vec->elements[i];

  // write result to file
  agile::writeVectorFile(result_file, x_host);

  //---------------------------------------------------------------------------
  // CLEANUPS FROM HERE
  //---------------------------------------------------------------------------
/*
*     Free the memory allocated for LSQR.
*/
  lsqr->freeLsqrMem();

  delete lsqr;
  delete model;
}
