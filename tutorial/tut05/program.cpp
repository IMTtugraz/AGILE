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

// $Id: program.cpp 476 2011-06-16 08:54:14Z freiberger $

#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"

// In order to use the message passing interface (MPI), we have to include the
// header for a communicator.
#include "agile/network/communicator.hpp"

#include <iostream>
#include <vector>

// Output a vector to \p std::cout
void output(const char* string, const std::vector<double>& v)
{
  // We output the vector together with the rank of the process. More on this
  // later...
  std::cout << "Process " << agile::Network::rank() << ": " << string;
  for (unsigned counter = 0; counter < v.size(); ++counter)
    std::cout << v[counter] << " ";
  std::cout << std::endl;
}

// The main routine takes also the command line arguments. These are passed on
// to the MPI later on.
int main(int argc, char* argv[])
{
  // The first step is to initialize the network communication. If the message
  // passing interface is installed, the \p NetworkEnvironment will
  // initialize it in the constructor. If MPI is not found, this will just
  // create an environment with a single process (the current one) in it.
  // You can use the network as long as the network environment object lives.
  // The destructor of \p NetworkEnvironment will finalize MPI.
  agile::NetworkEnvironment environment(argc, argv);

  // The singleton \p Network can be used to access very basic properties of
  // the parallel network. For example, you can use \p Network::size() to
  // get the number of processes that are running this program or you can
  // get the rank of this process by calling \p Network::rank().
  // We use the \p Network singleton to output some basic process information.
  std::cout << "Process " << agile::Network::rank() << " of "
            << agile::Network::size() << " is running on "
            << agile::Network::processor_name() << std::endl;

  // We introduce a barrier here. This is just a "meeting point" for all the
  // threads. The processes will only continue once all processes have reached
  // this point. We do this because the threads are not synchronized and the
  // output might be horribly disordered otherwise. If you care for speed, it
  // is possibly a good idea to remove the barriers.
  agile::Network::barrier();

  // This is a parallel example which is why we need at least two processes.
  // If there are less, we print a message to \p std::cout and quit right now.
  if (agile::Network::size() < 2)
  {
    if (agile::Network::rank() == 0)
    {
      std::cout << "Please run more than only a single process!" << std::endl;
      std::cout << "This is usually done using mpirun. ";
      std::cout << "The call should look similar to: " << std::endl;
      std::cout << "  mpirun -np 2 " << argv[0] << std::endl;
    }
    return 0;
  }

  // The next step is to create a communicator. This object is responsible
  // for transfering data over the network. Further, it eases the handling of
  // vectors that are shared by more than one process and stored in
  // distributed memory.
  // The communicator takes two template arguments. The first one is the
  // type used for vector indices (usually int or unsigned). The second one
  // is a floating point type that is compatible to the vectors used later on.
  // This type is needed to store the reciprocal (thus, floating point type)
  // of the number of nodes that share a vector entry.
  typedef agile::Communicator<unsigned, double> communicator_type;
  communicator_type com;

  // One of the easiest things is to collect a value from all processes and
  // to compute its sum. The communicator provides the \p collect() method for
  // doing so.
  // We generate one value per process and collect it:
  unsigned value = agile::Network::rank() + 1;
  std::cout << "Process " << agile::Network::rank()
            << " has got the local value " << value << std::endl;
  com.collect(value);
  // Every process has now the sum of all local value as can be seen by
  // printing it:
  std::cout << "Process " << agile::Network::rank()
            << " has got the global value " << value << std::endl;

  // Synchronize the threads again.
  agile::Network::barrier();

  // We want to demonstrate the abilities of the communicator by creating a
  // vector that lives partly on the first and partly on the second process.
  // The first process shall belong to indices \f$ 0 \le i \le 6 \f$ and the
  // second one to \f$ 3 \le i \le 9 \f$. Thus, the indices
  // \f$ 3 \le i \le 6 \f$ are shared by both processes. If a vector entry
  // does not belong to a process we do not store it. Thus, we need a way
  // to tell the library which global elements we store locally.

  // This is achieved using a local->global index map. This is simply a vector
  // whose \p i-th element holds the global index. For process 1, this has to
  // be the vector [0, 1, 2, 3, 4, 5, 6] and for process 2 the vector
  // [3, 4, 5, 6, 7, 8, 9].
  std::vector<unsigned> local_to_global_map;
  if (agile::Network::rank() == 0)
  {
    for (unsigned counter = 0; counter <= 6; ++counter)
      local_to_global_map.push_back(counter);
  }
  else if (agile::Network::rank() == 1)
  {
    for (unsigned counter = 3; counter <= 9; ++counter)
      local_to_global_map.push_back(counter);
  }
  // This local to global map we pass to the communicator.
  com.setLocal2GlobalMap(local_to_global_map.begin(),
                         local_to_global_map.end());

  // Now we construct a global vector \f$ v(i) = 2*i \f$.
  std::vector<double> global_v;
  for (unsigned counter = 0; counter <= 9; ++counter)
    global_v.push_back(2 * counter);

  // Create the local vectors from the global one and print them. Note that
  // these vectors are \p accumulated (because every process knows the full
  // value of the shared nodes).
  std::vector<double> local_v;
  for (unsigned counter = 0; counter < local_to_global_map.size(); ++counter)
    local_v.push_back(global_v[local_to_global_map[counter]]);
  if (agile::Network::rank() < 2)
    output("accumulated v_i: ", local_v);

  // Synchronize the threads again.
  agile::Network::barrier();

  // The communicator provides means to convert an \p accumulated vector to a
  // \p distributed one (method \p distribute()) and vice versa (method
  // \p accumulate()). We distribute the vector and output it:
  com.distribute(local_v);
  if (agile::Network::rank() < 2)
    output("distributed v_i: ", local_v);

  // Synchronize threads again.
  agile::Network::barrier();

  // The \p distributed vector can be accumulated again:
  com.accumulate(local_v);
  if (agile::Network::rank() < 2)
    output("accumulated v_i: ", local_v);
  agile::Network::barrier();

  // To compute the scalar product, we need a \p distributed and an
  // \p accumulated vector. \p local_v is still \p accumulated, so we copy and
  // distribute it. Note: It is generally more efficient to distribute a
  // vector than to accumulate one. The reason is that distribution can be
  // done without communicating with other processes because the shared vector
  // elements have to be divided only. Accumulation is more expensive because
  // the value from the other processes is needed which requires parallel
  // communication.
  std::vector<double> local_v_dist(local_v);
  com.distribute(local_v_dist);

  // Compute the local (partial) scalar product and print it to \p std::cout.
  double scalar_product = 0;
  for (unsigned counter = 0; counter < local_v.size(); ++counter)
    scalar_product += local_v[counter] * local_v_dist[counter];
  if (agile::Network::rank() < 2)
    std::cout << "Process " << agile::Network::rank()
              << ": local scalar product = " << scalar_product << std::endl;
  agile::Network::barrier();

  // Now we collect all the local scalar products and print the result.
  com.collect(scalar_product);
  std::cout << "Process " << agile::Network::rank()
            << ": global scalar product = " << scalar_product << std::endl;
  agile::Network::barrier();

  // The process with rank 0 shall compute the reference solution from the
  // global vector.
  if (agile::Network::rank() == 0)
  {
    double reference = 0;
    for (unsigned counter = 0; counter < global_v.size(); ++counter)
      reference += global_v[counter] * global_v[counter];
    std::cout << "Reference solution = " << reference << std::endl;
  }

  // This concludes the introduction to parallel programming. Once again we
  // want to state, that the \p barrier methods are not needed at all. The
  // program will still run as expected without failures.
  return 0;
}

// End of $Id: program.cpp 476 2011-06-16 08:54:14Z freiberger $.
