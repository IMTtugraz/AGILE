/*
* Model.hpp
* Defines the Model class for use by Lsqr
*
* 16 May 2003: First C++ version by J.A. Tomlin
*/

#ifndef _MATRIX_
#define _MATRIX_

#include "LsqrTypes.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <vector>

/**
 * This class implements a (sparse) matrix for pdco
 */
class Model {

public:

  uint32 nrow;
  uint32 ncol;

  std::vector<unsigned> *A_row_nnz;
  std::vector<unsigned> *A_column_index;
  std::vector<float> *A_data;

  std::vector<unsigned> *AT_row_nnz;
  std::vector<unsigned> *AT_column_index;
  std::vector<float> *AT_data;

/*
  uint32 numElts;
  int32 *rowIndex;
  uint32 *colStarts;
  double  *values;
  double  *rowUpper;
  double  *rowLower;
  double  *colUpper;
  double  *colLower;
  double  *rhs;
  prod_data *prod;
*/
  
  Model(uint32 rowSiz, uint32 colSiz) {
    nrow = rowSiz;
    ncol = colSiz;
  } 

  ~Model(){
  }

uint32 rowSize(){
    return nrow;
  }  
  
uint32 colSize(){
    return ncol;
  }  

void matVecMult( long, fvec *, fvec *);
  
};
#endif
