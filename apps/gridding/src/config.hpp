// $Id: config.hpp 476 2011-06-16 08:54:14Z freiberger $
#ifndef _GRIDDING_CONFIG_HPP_
#define _GRIDDING_CONFIG_HPP_

//--------------------------------------------------
// C O M M O N   S E T T I N G S
//--------------------------------------------------

// data type of transformation matrix, e.g. float or std::complex<float>
//#define TRANSFORMATION_MATRIX_DATATYPE float
#define TRANSFORMATION_MATRIX_DATATYPE std::complex<float>

// data type of vector x and y, e.g. float or std::complex<float>
#define DATA_VECTOR_DATATYPE std::complex<float>

// show timer output or not
#define WITH_TIMER 1

// LSQR stopping criteria,
//  - max. number of iterations
#define LSQR_MAX_ITERATIONS 100
//  - absolute tolerance
#define LSQR_ABS_TOLERANCE 1e-3

// show lsqr details, iterations, norm of residual, ...
#define WITH_LSQR_DETAILS 1


//--------------------------------------------------
// F O R W A R D   G R I D D I N G
// (radial data -> grid data)
//--------------------------------------------------

// no explicit settings at the moment

//--------------------------------------------------
// B A C K W A R D   G R I D D I N G
// (grid data -> radial data)
//--------------------------------------------------

// the backward gridding could be done by interpolation or multiplication
#define BACKWARD_GRIDDING_WITH_INTERPOLATION 1

//--------------------------------------------------
// L O O P   G R I D D I N G
//--------------------------------------------------

// predefined number of loops
#define LOOP_GRIDDING_NUMBER_OF_LOOPS 100

// the backward gridding could be done by interpolation or multiplication
#define LOOP_BACKWARD_GRIDDING_WITH_INTERPOLATION 0

// reset x vector within loop to zero vector
#define LOOP_GRIDDING_RESET_X_TO_ZERO 1


#endif // _GRIDDING_CONFIG_HPP_

// End of $Id: config.hpp 476 2011-06-16 08:54:14Z freiberger $.
