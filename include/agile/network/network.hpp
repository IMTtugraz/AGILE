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

// $Id: network.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_NETWORK_NETWORK_HPP
#define AGILE_NETWORK_NETWORK_HPP

#include "agile/gpu_config.hpp"

// choose the correct network implementation depending on the availability of
// the MPI library
#ifdef HAVE_MPI
#include "agile/network/mpi_network.hpp"
#else
#include "agile/network/single_network.hpp"
#endif

#endif // AGILE_NETWORK_NETWORK_HPP

// End of $Id: network.hpp 476 2011-06-16 08:54:14Z freiberger $.
