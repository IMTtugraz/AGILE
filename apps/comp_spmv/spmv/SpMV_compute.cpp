/*
 * IBM Sparse Matrix-Vector Multiplication Toolkit for Graphics Processing Units
 *  (c) Copyright IBM Corp. 2008, 2009.  All Rights Reserved.
 */ 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "SpMV.h"

void computeSpMV(float *x, const float *val, const unsigned int *rowIndices, const unsigned int *indices, 
		 const float *y, const unsigned int numRows)
{
    for (unsigned int i=0; i<numRows; i++) {
	float t = 0;
	unsigned int lb = rowIndices[i]; 
	unsigned int ub = rowIndices[i+1];
	for (unsigned int j=lb; j<ub; j++) {
	    unsigned int ind = indices[j];
	    t += val[j] * y[ind];
	}
	x[i] = t;
    }
}

void computeSpMV_BCSR(float *x, const float *val, const unsigned int *rowIndices,
                 const unsigned int *indices, const float *y, const unsigned int numRows,
		 const unsigned int numCols, const unsigned int bsx, const unsigned int bsy)
{
    float *t = (float *)malloc(sizeof(float)*bsx);

    for (unsigned int i=0; i< ceild(numRows,bsx); i++) {
        memset(t,0,sizeof(float)*bsx);
        unsigned int lb = rowIndices[i];
        unsigned int ub = rowIndices[i+1];
        unsigned int r = i*bsx;
        for (unsigned int j=lb; j<ub; j++) {
            unsigned int ind = indices[j];
	    unsigned int c = ind*bsy;
	    unsigned int commonInd = j*bsx*bsy;

	    if (((c+bsy) > numCols) || ((r+bsx) > numRows)) {
                for (unsigned int bi=0; (r+bi) < min(numRows,r+bsx); bi++) {
                    float tb=0;
                    for (unsigned int bj=0; (c+bj) < min(numCols,c+bsy); bj++)
                    	tb += val[commonInd+bi*bsy+bj] * y[c+bj];
                    t[bi]+=tb;
                }
	    }
	    else {
		#if BCSR_c8 
		// Assuming bsy as fixed (8)
	    	float y0 = y[c];
	    	float y1 = y[c+1];
	    	float y2 = y[c+2];
	    	float y3 = y[c+3];
	    	float y4 = y[c+4];
	    	float y5 = y[c+5];
	    	float y6 = y[c+6];
	    	float y7 = y[c+7];
		#if BCSR_r2
            	t[0]   +=   val[commonInd] * y0 + val[commonInd+1] * y1 +\
            	           val[commonInd+2] * y2 + val[commonInd+3] * y3 +\
			   val[commonInd+4] * y4 + val[commonInd+5] * y5 +\
			   val[commonInd+6] * y6 + val[commonInd+7] * y7;
            	t[1] +=   val[commonInd+bsy] * y0 + val[commonInd+bsy+1] * y1 +\
            	           val[commonInd+bsy+2] * y2 + val[commonInd+bsy+3] * y3 +\
			   val[commonInd+bsy+4] * y4 + val[commonInd+bsy+5] * y5 +\
			   val[commonInd+bsy+6] * y6 + val[commonInd+bsy+7] * y7;
		#else
	    	for (unsigned int bi=0; (r+bi) < min(numRows,r+bsx); bi++) {
            	    t[bi] += val[commonInd+bi*bsy] * y0 + val[commonInd+bi*bsy+1] * y1 +\
            	           val[commonInd+bi*bsy+2] * y2 + val[commonInd+bi*bsy+3] * y3 +\
			   val[commonInd+bi*bsy+4] * y4 + val[commonInd+bi*bsy+5] * y5 +\
			   val[commonInd+bi*bsy+6] * y6 + val[commonInd+bi*bsy+7] * y7;
		}
		#endif
		#endif
		#if BCSR_c4 
		// Assuming bsy as fixed (4)
	    	float y0 = y[c];
	    	float y1 = y[c+1];
	    	float y2 = y[c+2];
	    	float y3 = y[c+3];
		#if BCSR_r2
            	t[0]   +=   val[commonInd] * y0 + val[commonInd+1] * y1 +\
            	           val[commonInd+2] * y2 + val[commonInd+3] * y3;
            	t[1] +=   val[commonInd+bsy] * y0 + val[commonInd+bsy+1] * y1 +\
            	           val[commonInd+bsy+2] * y2 + val[commonInd+bsy+3] * y3;
		#else
	    	for (unsigned int bi=0; (r+bi) < min(numRows,r+bsx); bi++) {
            	    t[bi] += val[commonInd+bi*bsy] * y0 + val[commonInd+bi*bsy+1] * y1 +\
            	           val[commonInd+bi*bsy+2] * y2 + val[commonInd+bi*bsy+3] * y3;
		}
		#endif
		#endif
		#if BCSR_c2
                // Assuming bsy as fixed (2)
                float y0 = y[c];
                float y1 = y[c+1];
		#if BCSR_r2
            	t[0]   +=   val[commonInd] * y0 + val[commonInd+1] * y1;
            	t[1] +=   val[commonInd+bsy] * y0 + val[commonInd+bsy+1] * y1; 
		#else
                for (unsigned int bi=0; (r+bi) < min(numRows,r+bsx); bi++) {
                    t[bi] += val[commonInd+bi*bsy] * y0 + val[commonInd+bi*bsy+1] * y1;
                }
		#endif
		#endif
	    }
        }
	for (unsigned int bi=0; (r+bi) < min(numRows,r+bsx); bi++) 
            x[r+bi] = t[bi];
    }
    free(t);
}

