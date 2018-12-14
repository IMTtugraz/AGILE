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

// $Id: mpi_network.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_NETWORK_MPI_NETWORK_HPP
#define AGILE_NETWORK_MPI_NETWORK_HPP

#include <string>
#include <vector>
#include <mpi.h>

namespace agile
{
  //! \brief Class used to initialise the network infrastructure.
  //!
  //! You can use the network infrastructure as long as an object of this
  //! class is alive. After calling the constructor of the \p NetworkEnvironment
  //! it is not allowed to use the network any longer.
  //! \author Manuel Freiberger
  class NetworkEnvironment
  {
    public:
      //! \brief Constructor.
      //!
      //! The constructor initialises the network.
      NetworkEnvironment(int& argc, char**& argv)
      {
        MPI_Init(&argc, &argv);
      }

      //! \brief Desctructor.
      //!
      //! The destructor closes the network infrastructure. No more calls to
      //! network functions are allowed after that point.
      ~NetworkEnvironment()
      {
        MPI_Finalize();
      }
  };

  /**
   * The network class provides global network communication functions using
   * the Message Passing Interface (MPI).
   * \author Manfred Liebmann
   */
  class Network
  {
    public:
      /**
       * The time procedure initializes the network infrastructure.
       * \param _t Output: Time in seconds since an arbitrary time in the past.
       */
      static void time(double &_t)
      {
        _t = MPI_Wtime();
      }

      /**
      * The rank procedure returns the rank of the current process.
      * \param _rank Output: Rank of the current process.
      */
      static void rank(unsigned &_rank)
      {
        MPI_Comm_rank(MPI_COMM_WORLD, (int*)&_rank);
      }

      //! \brief Get the rank of the current process.
      static unsigned rank()
      {
        unsigned temp;
        Network::rank(temp);
        return temp;
      }

      /**
       * The size procedure returns the total number of processes.
       * \param _size Output: Total number of processes.
       */
      static void size(unsigned &_size)
      {
        MPI_Comm_size(MPI_COMM_WORLD, (int*)&_size);
      }

      //! \brief Get the total number of processes.
      static unsigned size()
      {
        unsigned temp;
        Network::size(temp);
        return temp;
      }

      /**
       * The barrier procedure syncronizes all processes.
       */
      static void barrier()
      {
        MPI_Barrier(MPI_COMM_WORLD);
      }

      /**
       * The alltoallv procedure initiates an all-to-all communication using data blocks of varying size.
       * \param _s Input: Address of the send buffer.
       * \param _slen Input: Array of the number of bytes to send to each process.
       * \param _soff Input: Array of data displacements in the send buffer for each process.
       * \param _r Input: Address of the receive buffer .
       * \param _rlen input: Array of the number of bytes to receive from any process.
       * \param _roff input: Array of data displacements in the receive buffer for any process.
       */
      static void alltoallv(void* _s, int* _slen, int* _soff, void* _r, int* _rlen, int* _roff)
      {
        MPI_Alltoallv(_s, _slen, _soff, MPI_CHAR, _r, _rlen, _roff, MPI_CHAR, MPI_COMM_WORLD);
      }

      /**
       * The allgather procedure initiates an all-gather communication using
       * data blocks of the same size.
       * \param _s Input: Address of the send buffer.
       * \param _slen Input: Number of bytes to send to each process.
       * \param _r Input: Address of the receive buffer .
       * \param _rlen input: Number of bytes to receive from any process.
       */
      static void allgather(void* _s, int _slen, void* _r, int _rlen)
      {
        MPI_Allgather(_s, _slen, MPI_CHAR, _r, _rlen, MPI_CHAR, MPI_COMM_WORLD);
      }

      //! \brief Gather a value from all processes.
      //!
      //! This method gathers a value from all processes and stores the value
      //! of the \p i-th process in \p y[i].
      //! \param[in] x Value from the current process.
      //! \param[out] y Vector containing the value from all processes.
      template <typename TType>
      static void allgather(const TType& x, std::vector<TType>& y)
      {
        y.resize(Network::size());
        MPI_Allgather((void*)&x, sizeof(TType), MPI_CHAR, &y[0], sizeof(TType),
                      MPI_CHAR, MPI_COMM_WORLD);
      }

      //! \brief Gather a string.
      //!
      //! \param[in] x The string from the current process.
      //! \param[out] y A vector containing the strings from all processes.
      //template <>
      static void allgather(const std::string& x, std::vector<std::string>& y)
      {
        // communicate the string lengths
        std::vector<unsigned> byte_sizes;
        allgather(unsigned(x.size()), byte_sizes);

        // create a vector of offsets
        std::vector<unsigned> byte_offset(Network::size() + 1, 0);
        for (unsigned counter = 1; counter <= Network::size(); ++counter)
          byte_offset[counter] = byte_offset[counter - 1]
                                 + byte_sizes[counter - 1];

        // a vector of char big enough for all the strings
        std::vector<char> receive_buffer(byte_offset[Network::size()]);
        // gather all the strings
        MPI_Allgatherv((void*)x.c_str(), x.size(), MPI_CHAR,
                       &receive_buffer[0], (int*)&byte_sizes[0],
                       (int*)&byte_offset[0], MPI_CHAR, MPI_COMM_WORLD);

        // fill the strings into a vector of strings
        y.resize(Network::size());
        for (unsigned counter = 0; counter < Network::size(); ++counter)
          y[counter] = std::string(&receive_buffer[0] + byte_offset[counter],
                                   byte_sizes[counter]);
      }

      //! \brief Gather a vector.
      //!
      //! \param[in] x The vector from the current process.
      //! \param[out] y A vector of vectors containing the vector from all
      //! processes.
      template <typename TType>
      static void allgather(const std::vector<TType>& x,
                            std::vector<std::vector<TType> >& y)
      {
        // communicate the vector sizes
        std::vector<unsigned> byte_sizes;
        allgather<unsigned>(unsigned(x.size()) * sizeof(TType), byte_sizes);

        // create a vector of offsets
        std::vector<unsigned> byte_offset(Network::size() + 1, 0);
        for (unsigned counter = 1; counter <= Network::size(); ++counter)
          byte_offset[counter] = byte_offset[counter - 1]
                                 + byte_sizes[counter - 1];

        // a vector of char big enough for all the vectors
        std::vector<char> receive_buffer(byte_offset[Network::size()]);
        // gather all the vectors
        MPI_Allgatherv((void*)&x[0], x.size() * sizeof(TType), MPI_CHAR,
                       &receive_buffer[0], (int*)&byte_sizes[0],
                       (int*)&byte_offset[0], MPI_CHAR, MPI_COMM_WORLD);

        // fill the vectors into a vector of vectors
        y.resize(Network::size());
        for (unsigned counter = 0; counter < Network::size(); ++counter)
        {
          y[counter].resize(byte_sizes[counter] / sizeof(TType));
          memcpy(&y[counter][0], &receive_buffer[0] + byte_offset[counter],
                 byte_sizes[counter]);
        }
      }

      /**
       * The allgatherv procedure initiates an all-gather communication using
       * data blocks of varying size.
       * \param _s Input: Address of the send buffer.
       * \param _slen Input: Number of bytes to send to each process.
       * \param _r Input: Address of the receive buffer .
       * \param _rlen input: Array of the number of bytes to receive from any
       * process.
       * \param _roff input: Array of data displacements in the receive buffer
       * for any process.
       */
      static void allgatherv(void* _s, int _slen, void* _r, int* _rlen,
                             int* _roff)
      {
        MPI_Allgatherv(_s, _slen, MPI_CHAR, _r, _rlen, _roff, MPI_CHAR,
                       MPI_COMM_WORLD);
      }

      //! \brief Get the name of the processor.
      //!
      //! This method returns a string containing the name of the processor
      //! this process is running on.
      //! \return The name of the processor.
      static std::string processor_name()
      {
        char name[MPI_MAX_PROCESSOR_NAME];
        int length;
        MPI_Get_processor_name(name, &length);
        return std::string(name, length);
      }
  };

} // namespace agile

#endif // AGILE_NETWORK_MPI_NETWORK_HPP

// End of $Id: mpi_network.hpp 476 2011-06-16 08:54:14Z freiberger $.
