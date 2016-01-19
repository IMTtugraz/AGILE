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

// $Id: agile.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_AGILE_HPP
#define AGILE_AGILE_HPP

// I added this file to include more or less the whole AGILE library.
// For some quick and dirty examples, it is just convenient to have a do-it-all
// header :-)

// constants, definitions, configuration
#include "agile/gpu_config.hpp"

#include "agile/exception.hpp"

// CPU vectors and matrices
#include "agile/cpu_matrix.hpp"
#include "agile/cpu_vector.hpp"

// GPU vectors and matrices
#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_vector_range.hpp"
#include "agile/gpu_host_vector.hpp"
#include "agile/gpu_matrix_pitched.hpp"
#include "agile/gpu_matrix.hpp"
#include "agile/gpu_cs_matrix.hpp"

// other basic objects
#include "agile/gpu_timer.hpp"
#include "agile/radix_exchange_sort.hpp"
#include "agile/ublas_bindings.hpp"

// networking
#include "agile/network/communicator.hpp"
#include "agile/network/gpu_communicator.hpp"
#include "agile/network/network.hpp"

// basic operators
#include "agile/operator/forward_operator.hpp"
#include "agile/operator/inverse_operator.hpp"

// solvers (inverse operators)
#include "agile/operator/diagonal_solver.hpp"
#include "agile/operator/cg.hpp"
#include "agile/operator/pcg.hpp"
#include "agile/operator/gmres.hpp"
#include "agile/operator/lsqr.hpp"
#include "agile/operator/minres.hpp"

// preconditioners (TODO: remove this)
#include "agile/operator/jacobi.hpp"

// IO
#include "agile/io.hpp"

#endif // AGILE_AGILE_HPP

// End of $Id: agile.hpp 476 2011-06-16 08:54:14Z freiberger $.
