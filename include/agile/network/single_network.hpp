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

// $Id: single_network.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_NETWORK_SINGLE_NETWORK_HPP
#define AGILE_NETWORK_SINGLE_NETWORK_HPP

// This code is based on work done by Manfred Liebmann.

#include <ctime>
#include <cstring>  // memcpy
#include <string>
#include <vector>

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
      }

      //! \brief Desctructor.
      //!
      //! The destructor closes the network infrastructure. No more calls to
      //! network functions are allowed after that point.
      ~NetworkEnvironment()
      {
      }
  };

  /**
   * The network class provides global network communication functions for a
   * single process not using the Message Passing Interface (MPI).
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
        _t = 1.0/double(CLOCKS_PER_SEC)*clock();
      }

      /**
       * The rank procedure returns the rank of the current process.
       * \param _rank Output: Rank of the current process.
       */
      static void rank(unsigned &_rank)
      {
        _rank = 0;
      }

      //! \brief Get the rank of the current process.
      static unsigned rank()
      {
        return 0;
      }

      /**
       * The size procedure returns the total number of processes.
       * \param _size Output: Total number of processes.
       */
      static void size(unsigned &_size)
      {
        _size = 1;
      }

      //! \brief Get the total number of processes.
      static unsigned size()
      {
        return 1;
      }

      /**
       * The barrier procedure syncronizes all processes.
       */
      static void barrier()
      {
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
        memcpy(_r, _s, _slen[0]);
      }

      /**
       * The allgather procedure initiates an all-gather communication using data blocks of the same size.
       * \param _s Input: Address of the send buffer.
       * \param _slen Input: Number of bytes to send to each process.
       * \param _r Input: Address of the receive buffer .
       * \param _rlen input: Number of bytes to receive from any process.
       */
      static void allgather(void* _s, int _slen, void* _r, int _rlen)
      {
        memcpy(_r, _s, _slen);
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
        y.assign(1, x);
      }

      /**
       * The allgatherv procedure initiates an all-gather communication using data blocks of varying size.
       * \param _s Input: Address of the send buffer.
       * \param _slen Input: Number of bytes to send to each process.
       * \param _r Input: Address of the receive buffer .
       * \param _rlen input: Array of the number of bytes to receive from any process.
       * \param _roff input: Array of data displacements in the receive buffer for any process.
       */
      static void allgatherv(void* _s, int _slen, void* _r, int* _rlen, int* _roff)
      {
        memcpy(_r, _s, _slen);
      }

      //! \brief Get the name of the processor.
      static std::string processor_name()
      {
        // TODO: find out the correct name
        return "TODO";
      }
    };

} // namespace agile

#endif // AGILE_NETWORK_SINGLE_NETWORK_HPP

// End of $Id: single_network.hpp 476 2011-06-16 08:54:14Z freiberger $.
