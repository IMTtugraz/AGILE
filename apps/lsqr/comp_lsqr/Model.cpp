/*
* Model.cpp
*    Code for test model
*/

#include "Model.hpp"

void Model::matVecMult( long  mode, fvec  *x, fvec  *y) {

/*
*     Compute  Y = Y + A*X
*/
  if ( mode == 0 )
  {
    unsigned row, nz_in_row, idx = 0;
    for (row=0; row<A_row_nnz->size(); ++row)
    {
      nz_in_row = (*A_row_nnz)[row];
      if (nz_in_row > 0) {
        while (nz_in_row--)
        {
          y->elements[row] += (*A_data)[idx] * x->elements[(*A_column_index)[idx]];
          idx++;
        }
      }
    }
  }

/*
*     Compute  X = X + A^T*Y
*/
  if ( mode == 1 )
  {
    unsigned row, nz_in_row, idx = 0;
    for (row=0; row<AT_row_nnz->size(); ++row)
    {
      nz_in_row = (*AT_row_nnz)[row];
      if (nz_in_row > 0) {
        while (nz_in_row--)
        {
          x->elements[row] += (*AT_data)[idx] * y->elements[(*AT_column_index)[idx]];
          idx++;
        }
      }
    }
  }

  return;
}
