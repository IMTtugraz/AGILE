// Copyright (C) 2010-2011 Institute of Medical Engineering,
// Graz University of Technology
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses>.

// $Id: interp2d.cpp 476 2011-06-16 08:54:14Z freiberger $
#include "agile/gpu_environment.hpp"
#include "agile/network/gpu_communicator.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_matrix_pitched.hpp"
#include "agile/io/file.hpp"

#include "agile/quadtree.hpp"
#include <vector>

#include <time.h>

#include <iomanip>
#include <sys/time.h>

#include "../../test/test_defines.h"

// helper defines
#define MAX(a, b) (a > b ? a : b)
#define LD(x) (log(x) / log(2))

#define NUM_ITER 1000
#define NUM_RESORT_TYPES 6

// depends on source data to be stored in tex cache
// e.g. sizeof(float) * BLOCKED_LINE_BY_LINE_MAX_BLOCK_WIDTH / 1024
//      is the need for memory for one block row in texture cache
#define BLOCKED_LINE_BY_LINE_MAX_BLOCK_WIDTH 4096


// Line-By-Line reorganization
bool compPosLineByLine( const std::complex<float>& a, const std::complex<float>& b ) {
  unsigned aCol = floor(a.real());
  unsigned aRow = floor(a.imag());
  unsigned bCol = floor(b.real());
  unsigned bRow = floor(b.imag());
  return aRow == bRow ? aCol < bCol : aRow < bRow;
}

// Blockwise Line-By-Line reorganisation
bool compPosBlockedLineByLine( const std::complex<float>& a, const std::complex<float>& b ) {
  // devide the source pic in vertical strips and go through them line-by-line

  //TODO: interpolation positions outside the source data (<0 and >width)
  // are not handled correctly.
  // in CLAMP addressing mode:
  // -> X positions lower than 0 should be within the first strip
  // -> X positions greater than source width should be within the last strip
  //
  // in WRAP addressing mode:
  // X positions lower than 0 and greater than source width should be mapped into
  // the correct strip by: x = frac(x)
  // where frac(x) = x − floor(x) and floor(x) is the largest integer not greater than x
  // @see NVIDIA CUDA Guide 2.3.1 Appendix E "Texture Fetching"

  unsigned aCol = floor(a.real());
  unsigned aRow = floor(a.imag());
  unsigned bCol = floor(b.real());
  unsigned bRow = floor(b.imag());

  // we asume interpolation points only within the source data (between 0 and width or height)
  int aStripId = aCol / BLOCKED_LINE_BY_LINE_MAX_BLOCK_WIDTH;
  int bStripId = bCol / BLOCKED_LINE_BY_LINE_MAX_BLOCK_WIDTH;
  if (aStripId == bStripId)
  {
    // within a strip order the positions line-by-line
    return aRow == bRow ? aCol < bCol : aRow < bRow;
  }
  return aStripId < bStripId;
}

// Advanced blockwise Line-By-Line reorganisation with alternating rows within block
bool compPosBlockedAltLineByLine( const std::complex<float>& a, const std::complex<float>& b ) {
  // devide the source pic in vertical strips and go through them line-by-line

  //TODO: interpolation positions outside the source data (<0 and >width)
  // are not handled correctly.
  // in CLAMP addressing mode:
  // -> X positions lower than 0 should be within the first strip
  // -> X positions greater than source width should be within the last strip
  //
  // in WRAP addressing mode:
  // X positions lower than 0 and greater than source width should be mapped into
  // the correct strip by: x = frac(x)
  // where frac(x) = x − floor(x) and floor(x) is the largest integer not greater than x
  // @see NVIDIA CUDA Guide 2.3.1 Appendix E "Texture Fetching"

  unsigned aCol = floor(a.real());
  unsigned aRow = floor(a.imag());
  unsigned bCol = floor(b.real());
  unsigned bRow = floor(b.imag());

  // we asume interpolation points only within the source data (between 0 and width or height)
  unsigned aStripId = aCol / BLOCKED_LINE_BY_LINE_MAX_BLOCK_WIDTH;
  unsigned bStripId = bCol / BLOCKED_LINE_BY_LINE_MAX_BLOCK_WIDTH;

  // within same block?
  if (aStripId == bStripId)
  {
    // within a strip order the positions line-by-line
    if (aRow == bRow)
    {
      if (aRow % 2)
        return aCol > bCol;
      return aCol < bCol;
    }
    return aRow < bRow;
  }
  return aStripId < bStripId;
}



int main(int argc, char* argv[])
{
  struct timeval st, et;
  unsigned numRows = 0;
  unsigned numCols = 0;
  unsigned numInterpPositions = 0;

  if (argc != 4)
  {
    std::cout << "USAGE: " << argv[0] << " <src cols> <src rows> <positions>" << std::endl;
    std::cout << "ARGUMENTS:" << std::endl;
    std::cout << "<src cols> .... number of columns in source data" << std::endl;
    std::cout << "<src rows> .... number of rows in source data" << std::endl;
    std::cout << "<positions> ... number of randomized interpolation positions" << std::endl;
    std::cout << std::endl;
    return -1;
  } else {
    numCols = atoi(argv[1]);
    numRows = atoi(argv[2]);
    numInterpPositions = atoi(argv[3]);
  }

  std::cout << "#-----------------------------------------------------------------------" << std::endl;
  std::cout << "# source rows: " << numRows << ", columns: " << numCols << std::endl;
  std::cout << "# interpolation positions: " << numInterpPositions << std::endl;

  // init the network
  agile::NetworkEnvironment environment(argc, argv);

  // allocate a GPU
  typedef agile::GPUCommunicator<unsigned, float, float> communicator_type;
  communicator_type com;
  com.allocateGPU();


  //Test matrix and test position vector for interpolation
  //-------------------------------------------------------
  typedef std::vector<std::complex<float> > cpu_pos_vector_type;
  typedef agile::GPUVector<std::complex<float> > gpu_pos_vector_type;
  typedef agile::GPUMatrixPitched<float> gpu_source_matrix_type;

  float *tmp_matrix = new float[numRows * numCols];

  srand(time(NULL));

  // initialize matrix with random data
  //----------------------------------------------------------------
  for (unsigned i=(numRows * numCols); i--; )
    tmp_matrix[i] = (float)rand() / (float)rand();

/*
  for (unsigned row=0; row<numRows; ++row)
    for (unsigned col=0; col<numCols; ++col)
    {
      tmp_matrix[row * numCols + col] = rand() / rand();
    }
*/


  // initialize randomized interpolation positions
  //----------------------------------------------------------------
  agile::vec2_t bounds[2] = {
    agile::vec2_t(0.0f, 0.0f),
    agile::vec2_t((float)numCols, (float)numRows)
  };

  // QUADTREE for reordering interpolation positions
  // using logarithmus dualis for max depth computation
  agile::Quadtree quadtree(bounds, floor(MAX(LD(numCols), LD(numRows))));

  cpu_pos_vector_type pos_host(numInterpPositions);
  float rx, ry;
  for (unsigned i=0; i<numInterpPositions; ++i)
  {
    rx = (rand() % numCols) + (rand() / (float)RAND_MAX);
    ry = (rand() % numRows) + (rand() / (float)RAND_MAX);
    pos_host[i] = std::complex<float>(rx, ry);
    quadtree.addItem(new agile::QuadtreeItem(pos_host[i].real(), pos_host[i].imag(), i));
  }

  std::vector<agile::QuadtreeItem *> reorderdPositionList;
  cpu_pos_vector_type pos_host_resorted(numInterpPositions);

  gpu_pos_vector_type pos(pos_host.size());

  // Init GPU source matrix
  gpu_source_matrix_type source(numRows, numCols, tmp_matrix);
  // free host matrix data
  delete[] tmp_matrix;

  // init result vector on gpu (same size as position vector)
  typedef agile::GPUVector<gpu_source_matrix_type::value_type >
    gpu_result_vector_type;
  gpu_result_vector_type res(pos_host.size());

  float times[NUM_RESORT_TYPES];

  // iterate over all interpolation types
  // (no reordering, z-order curve, hilbert curve)
  for (unsigned reorderType=0; reorderType<NUM_RESORT_TYPES; ++reorderType)
  {

    if (reorderType == 1)
    {
      std::cout << "# Z-ORDER CURVE REORDERING ............ ";
      reorderdPositionList.clear();
      quadtree.getLinearZOrderList(reorderdPositionList);
      // real reordering of positions
      for (unsigned i=0; i < reorderdPositionList.size(); ++i)
        pos_host_resorted[i] = pos_host[reorderdPositionList[i]->orgIndex];
    }
    else if (reorderType == 2)
    {
      std::cout << "# HILBERT CURVE REORDERING ............ ";
      reorderdPositionList.clear();
      quadtree.getLinearHilbertOrderList(reorderdPositionList);
      // real reordering of positions
      for (unsigned i=0; i < reorderdPositionList.size(); ++i)
        pos_host_resorted[i] = pos_host[reorderdPositionList[i]->orgIndex];
    }
    else if (reorderType == 3)
    {
      std::cout << "# LINE-BY-LINE REORDERING ............. ";
      pos_host_resorted.assign(pos_host.begin(), pos_host.end());
      std::sort(pos_host_resorted.begin(), pos_host_resorted.end(), compPosLineByLine);
    }
    else if (reorderType == 4)
    {
      std::cout << "# BLOCK LINE-BY-LINE REORDERING ....... ";
      pos_host_resorted.assign(pos_host.begin(), pos_host.end());
      std::sort(pos_host_resorted.begin(), pos_host_resorted.end(), compPosBlockedLineByLine);
    }
    else if (reorderType == 5)
    {
      std::cout << "# BLOCK L-BY-L ALT ROWs REORDERING .... ";
      pos_host_resorted.assign(pos_host.begin(), pos_host.end());
      std::sort(pos_host_resorted.begin(), pos_host_resorted.end(), compPosBlockedAltLineByLine);
    }
    else
    {
      std::cout << "# NO REORDERING ....................... ";
      // use random positions as they are
      pos_host_resorted.assign(pos_host.begin(), pos_host.end());
    }

    //PRINT_VEC("as is   ", pos_host);
    //PRINT_VEC("resorted", pos_host_resorted);

    // Init GPU position vector
    pos.assignFromHost(pos_host_resorted.begin(), pos_host_resorted.end());

    // do the interpolation
    gettimeofday(&st, NULL);
    cudaThreadSynchronize();
    for (unsigned i=0; i<NUM_ITER; ++i)
    {
      agile::interp2d(source, pos, res);
      cudaThreadSynchronize();
    }
    gettimeofday(&et, NULL);
    times[reorderType] = ((et.tv_sec-st.tv_sec)*1000.0 + (et.tv_usec - st.tv_usec)/1000.0);

    std::cout << (times[reorderType] / (float)NUM_ITER) << "[ms]";
    if (reorderType > 0)
      std::cout << " (speedup: " << (times[0] /  times[reorderType]) << ", " << (floor((times[0] /  times[reorderType] - 1.0f) * 10000.0f + .5f) / 100.0f) << "%)";
    std::cout << std::endl;
  }

  std::cout << "# iterations;numRows;numCols;numInterpPositions;time(no-reordering);time(z-order);time(hilbert);time(lbl);time(blocked lbl);time(blocked lbl alt)" << std::endl;
  std::cout << NUM_ITER << ";" << numRows << ";" << numCols << ";" << numInterpPositions << ";"
            << times[0] << ";" << times[1] << ";" << times[2] << ";" << times[3] << ";" << times[4] << ";" << times[5] << std::endl;

  // free quadtree items
  // quadtree.getLinearZOrderList(reorderdPositionList);
  for (unsigned i=0; i < reorderdPositionList.size(); ++i)
    delete reorderdPositionList[i];

  return 0;
}

// End of $Id: interp2d.cpp 476 2011-06-16 08:54:14Z freiberger $.
