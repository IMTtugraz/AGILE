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

// We want to use MPI together with the GPU, which is why we include a
// GPUCommunicator via \p gpu_communicator.hpp.
#include "agile/network/gpu_communicator.hpp"

#include <iostream>
#include <vector>

// Output a vector to \p std::cout
void output(const char* string, const std::vector<double>& v)
{
  std::cout << "Process " << agile::Network::rank() << ": " << string;
  for (unsigned counter = 0; counter < v.size(); ++counter)
    std::cout << v[counter] << " ";
  std::cout << std::endl;
}

// The main routine takes also the command line arguments. These are passed on
// to the MPI.
int main(int argc, char* argv[])
{
  // The initialisation of the network is the same as in tutorial 5.
  agile::NetworkEnvironment environment(argc, argv);
  std::cout << "Process " << agile::Network::rank() << " of "
            << agile::Network::size() << " is running on "
            << agile::Network::processor_name() << std::endl;
  agile::Network::barrier();

  // Again, as this is a parallel example, we need some processes. If there is
  // only one process, we cannot do anything.
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

  // Instead of creating a \p Communicator as in the last tutorial, we create
  // a \p GPUCommunicator. The GPU communicator has the same abilities as
  // the ordinary communicator but it can also handle GPU vectors.
  // The \p GPUCommuncator uses <b>three</b> template arguments. The first two
  // are identical to the \p Communicator: the type used for vector indices
  // (in our case \p unsigned) and the floating point type which has to be
  // compatible to a CPU vector (we use a \p float) which is used to store
  // the reciprocal of the number of processes an element belongs to.
  // The third type is a floating point type which has to be compatible to a
  // <b>GPU</b> vector.
  // It is a good idea to make this type the same as the second argument
  // except if you know very well what you are doing!
  typedef agile::GPUCommunicator<unsigned, float, float> communicator_type;
  communicator_type com;

  // Note that we did not initialize the GPU until now. The GPUcommunicator
  // offers another method for doing so: \p allocateGPU(). This method
  // searches the best GPU available and allocates it to this process.
  // Processes with a lower rank will have the first choice. The method
  // will fail if there are not enough GPUs in the PC.
  if (!com.allocateGPU())
  {
    std::cerr << "Could not allocate a GPU for process "
              << agile::Network::rank() << "." << std::endl
              << "Are you sure you have enough GPUs?" << std::endl;
    return -1;
  }

  // Next we output the GPU that was allocated to the current process and
  // synchronize the processes afterwards.
  std::cout << "Process " << agile::Network::rank() << " uses GPU no. "
            << agile::GPUEnvironment::getDeviceID() << std::endl;
  agile::Network::barrier();

  // The rest of the tutorial is just a copy-and-paste of the previous one.
  // If you are familiar with using a communicator, you can just skip over
  // from here on.

  // Warm up by summing the ranks.
  unsigned value = agile::Network::rank() + 1;
  std::cout << "Process " << agile::Network::rank()
            << " has got the local value " << value << std::endl;
  com.collect(value);
  std::cout << "Process " << agile::Network::rank()
            << " has got the global value " << value << std::endl;
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

  return 0;
}

// End of $Id: program.cpp 476 2011-06-16 08:54:14Z freiberger $.
