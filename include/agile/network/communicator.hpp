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

// $Id: communicator.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_NETWORK_COMMUNICATOR_HPP
#define AGILE_NETWORK_COMMUNICATOR_HPP

// This code is based on previous work done by Manfred Liebmann.
// Clean-up, some re-structuring and GPU code by Manuel Freiberger.

#include <vector>
#include <algorithm>
#include <numeric>

#include "agile/gpu_config.hpp"
#include "agile/network/network.hpp"
#include "agile/radix_exchange_sort.hpp"

#include <vector>

namespace agile
{
  //! \brief Class for parallel communication.
  //!
  //! This class takes the global node numbers associated with the current
  //! process and establishes all the necessary communication with other
  //! processes.
  //! \p TIndexType is the type used for node indices. Usually this is
  //! either \p int or \p unsigned.
  //! \p TRealType is the type used to store the reciprocal of the amount
  //! of processes a shared node belongs to. Usually this is a \p double.
  //! \author Manfred Liebmann, Manuel Freiberger
  template <typename TIndexType, typename TRealType>
  class Communicator
  {
    public:
      //! \brief Default constructor.
      //!
      //! The communicator is not usable before \p setLocal2GlobalMap was
      //! called.
      Communicator()
      {
      }

      //! \brief Constructor.
      //!
      //! This constructor takes the global node numbers associated with the
      //! current process in vector form. This vector also establishes a map
      //! between the local node numbering and the global numbering (local
      //! node 0 is mapped to global[0], local node 1 to global[1] and so on).
      Communicator(const std::vector<TIndexType>& global)
      {
        setLocal2GlobalMap(global);
      }

      //! \brief Supply a vector of global nodes.
      //!
      //! Use this method to tell the communicator which global nodes are to
      //! be computed on the current process. The order of the nodes is
      //! important because the communicator assumes that the local node
      //! 0 corresponds to the global node _gnod[0], local node 1 to
      //! the global node _gnod[1] and so on.
      //! \param[in] _gnod The local->global mapping of node numbers.
      void setLocal2GlobalMap(const std::vector<TIndexType>& global)
      {
        setLocal2GlobalMap(global.begin(), global.end());
      }

      //! \brief Supply a vector of global nodes.
      //!
      //! Use this method to supply the global nodes as an iterator range.
      //! \param[in] start The start iterator for the global nodes.
      //! \param[in] end The end iterator for the global nodes.
      template <typename TIterator>
      void setLocal2GlobalMap(const TIterator& start, const TIterator& end)
      {
        // if this is the only process, we do not have to do anything
        if (Network::size() == 1)
          return;

        // copy the global nodes into a temporary vector that can be sorted
        std::vector<TIndexType> global_nodes;
        global_nodes.assign(start, end);
        unsigned num_nodes = global_nodes.size();

        // create a vector of local nodes; then sort the global nodes and
        // copy the corresponding local indices such that the mapping is
        // preserved
        std::vector<TIndexType> local_nodes(num_nodes);
        for (TIndexType counter = 0; counter < num_nodes; ++counter)
          local_nodes[counter] = counter;
        RadixExchangeSortCopy<typename std::vector<TIndexType>::iterator>::sort(
          global_nodes.begin(), global_nodes.end(),
          local_nodes.begin(), local_nodes.end());

        // communicate the amount of nodes
        Network::allgather(num_nodes, m_shared_node_count);
        // sum the amount of nodes on all processors
        unsigned total_num_nodes
          = std::accumulate(m_shared_node_count.begin(),
                            m_shared_node_count.end(), 0);

        // create the displacement vector
        m_shared_node_offset.resize(Network::size());
        m_shared_node_offset[0] = 0;
        for (unsigned counter = 1; counter < m_shared_node_offset.size();
             ++counter)
          m_shared_node_offset[counter] = m_shared_node_offset[counter - 1]
                                          + m_shared_node_count[counter - 1];

        // communicate all global node indices
        m_shared_node.resize(total_num_nodes);
        {
          std::vector<unsigned> block_size(m_shared_node_count.size());
          std::vector<unsigned> block_offset(m_shared_node_count.size());
          for(unsigned counter = 0; counter < m_shared_node_count.size();
              ++counter)
          {
            block_size[counter] = m_shared_node_count[counter]
                                  * sizeof(TIndexType);
            block_offset[counter] = m_shared_node_offset[counter]
                                    * sizeof(TIndexType);
          }
          Network::allgatherv(&global_nodes[0], num_nodes * sizeof(TIndexType),
                              &m_shared_node[0], (int*)&block_size[0],
                              (int*)&block_offset[0]);
        }

        // store the amount of global nodes
        m_num_global_nodes = 1 + *std::max_element(m_shared_node.begin(),
                                                   m_shared_node.end());

        // at this point we have
        // - the sorted global nodes of the current process (global_nodes)
        // - the corresponding local nodes of the current process (local_nodes)
        // - the global nodes of all processes (m_shared_node)
        // - the amount of nodes on each process (m_shared_node_count)
        // - the offset in m_shared_node for each process (m_shared_node_offset)
        // we intersect the global nodes of our process with the global nodes
        // from all other processes and store the corresponding local node
        // of this process
        // the intersection is done in-place overwriting m_shared_node which is
        // legal as m_shared_node is only traversed once
        unsigned result_index = 0;
        for (unsigned other_process = 0; other_process < Network::size();
             ++other_process)
        {
          // point to the beginning of our global nodes
          unsigned our_index = 0;
          // point to the beginning and past-the-end of the global nodes of the
          // other processor
          unsigned other_index = m_shared_node_offset[other_process];
          unsigned other_end = other_index + m_shared_node_count[other_process];
          // set the new displacement for the intersection result
          m_shared_node_offset[other_process] = result_index;

          if (other_process != Network::rank())
          {
            while ((our_index < global_nodes.size())
                   && (other_index < other_end))
            {
              if (global_nodes[our_index] < m_shared_node[other_index])
                ++our_index;
              else if (m_shared_node[other_index] < global_nodes[our_index])
                ++other_index;
              else
              {
                // this is an intersection because the values match
                // copy the local node number to the final vector
                m_shared_node[result_index] = local_nodes[our_index];
                ++result_index;
                ++our_index;
                ++other_index;
              }
            }
          }
          // set the amount of intersecting nodes
          m_shared_node_count[other_process]
            = result_index - m_shared_node_offset[other_process];
        }
        // reduce the vector such that only the intersection result is stored
        m_shared_node.resize(result_index);

        // clear all the vectors that we do not need any longer
        global_nodes.clear();
        local_nodes.clear();

        // create a vector holding only the intersecting local nodes and
        // a second one to store how often each node occurs
        m_unique_shared_node = m_shared_node;
        m_unique_shared_node_process_count.assign(m_shared_node.size(), 0);
        RadixExchangeSort<typename std::vector<TIndexType>::iterator>::sort(
          m_unique_shared_node.begin(), m_unique_shared_node.end());
        // make the intersection result unique and count the occurences
        result_index = 0;
        for (typename std::vector<TIndexType>::const_iterator iter
               = m_unique_shared_node.begin();
             iter != m_unique_shared_node.end(); ++iter)
        {
          if (*iter != m_unique_shared_node[result_index])
          {
            ++result_index;
            m_unique_shared_node[result_index] = *iter;
          }
          ++m_unique_shared_node_process_count[result_index];
        }
        // clip the non-unique values
        if (m_shared_node.size())
        {
          m_unique_shared_node.resize(result_index + 1);
          m_unique_shared_node_process_count.resize(result_index + 1);
        }

        // store the inverse of the occurences as this is needed later on all
        // the time; we also have to correct the number by one to account for
        // the fact that the local node also belongs to the current process
        // which is not stored in the vector m_shared_node
        for (typename std::vector<TRealType>::iterator iter
               = m_unique_shared_node_process_count.begin();
             iter != m_unique_shared_node_process_count.end(); ++iter)
          *iter = 1.0 / (1.0 + *iter);

        // frem: store the local->global index map as it is needed to
        // re-assemble vectors from all processes
        m_local_to_global_map.assign(start, end);
      }

      //! \brief Calculate the global sum of the input value.
      //!
      //! This method communicates the input value to all other processes and
      //! sums them afterwards.
      //! \param[in,out] value Input: The value of this process to sum;
      //! Output: The sum of all input values.
      template <typename TType>
      void collect(TType& value)
      {
        // if this is the only process, we do not have to do anything
        if (Network::size() == 1)
          return;

        unsigned network_size = Network::size();
        m_receive_buffer.resize(sizeof(TType) * network_size);
        // gather all the values from all processes
        Network::allgather(&value, sizeof(TType), &m_receive_buffer[0],
                           sizeof(TType));
        // sum all the values
        value = TType(0);
        TType* ptr = (TType*)(&m_receive_buffer[0]);
        for (unsigned counter = 0; counter < network_size; ++counter)
          value += ptr[counter];
      }

      //! \brief Calculate the global sum of the input vectors.
      //!
      //! This method communicates the input vector to all other processes and
      //! sums them afterwards.
      //! \param[in,out] value Input: The vector of this process to sum;
      //! Output: The sum of all input vectors.
      template <typename TVectorType>
      void collect(TVectorType& vector,
                   typename TVectorType::const_iterator())
      {
        // if this is the only process, we do not have to do anything
        if (Network::size() == 1)
          return;

        typedef typename TVectorType::value_type value_type;
        unsigned vector_size = vector.size();
        unsigned network_size = Network::size();
        // resize the receive buffer so that it can hold the values of all
        // vectors
        m_receive_buffer.resize(
          sizeof(value_type) * vector_size * network_size);
        Network::allgather(&vector[0], vector_size * sizeof(value_type),
                           &m_receive_buffer[0],
                           vector_size * sizeof(value_type));

        // for every element of the vector sum the values accross the processes
        value_type* ptr = (value_type*)(&m_receive_buffer[0]);
        for (unsigned element_counter = 0; element_counter < vector_size;
             ++element_counter)
        {
          value_type sum = 0;
          for (unsigned counter = 0; counter < network_size; ++counter)
            sum += ptr[counter * vector_size + element_counter];
          vector[element_counter] = sum;
        }
      }

      /**
       * The accumulate procedure converts a distributed vector to an
       * accumulated vector.
       * \param[in,out] x Input: Distributed vector, Output: Accumulated vector.
       */
      template <typename TVectorType>
      void accumulate(TVectorType& x)
      {
        // if this is the only process, we do not have to do anything
        if (Network::size() == 1)
          return;

        typedef typename TVectorType::value_type value_type;

        // set up the buffers for the transaction
        m_send_buffer.resize(m_shared_node.size() * sizeof(value_type));
        m_receive_buffer.resize(m_shared_node.size() * sizeof(value_type));
        value_type* ptr = (value_type*)(&m_send_buffer[0]);
        for (unsigned counter = 0; counter < m_shared_node.size(); ++counter)
          ptr[counter] = x[m_shared_node[counter]];

        // calculate the amount of bytes to transfer and the offset
        std::vector<unsigned> block_size(m_shared_node_count.size());
        std::vector<unsigned> block_offset(m_shared_node_count.size());
        for(unsigned counter = 0; counter < m_shared_node_count.size();
            ++counter)
        {
          block_size[counter] = m_shared_node_count[counter]
                                * sizeof(value_type);
          block_offset[counter] = m_shared_node_offset[counter]
                                  * sizeof(value_type);
        }
        Network::alltoallv(&m_send_buffer[0], (int*)&block_size[0],
                           (int*)&block_offset[0], &m_receive_buffer[0],
                           (int*)&block_size[0], (int*)&block_offset[0]);

        // sum the individual contributions
        ptr = (value_type*)(&m_receive_buffer[0]);
        for (unsigned counter = 0; counter < m_shared_node.size(); ++counter)
          x[m_shared_node[counter]] += ptr[counter];
      }

      /**
       * The distribute procedure converts an accumulated vector to a
       * distributed vector.
       * \param _x Input: Accumulated vector, Output: Distributed vector.
       */
      template <typename TVectorType>
      void distribute(TVectorType& x)
      {
        // if this is the only process, we do not have to do anything
        if (Network::size() == 1)
          return;

        // distributing is easy as we only have to loop over our local vector
        // and divide the shared nodes by the number of processes that node
        // belongs to; note that the division is implemented via a
        // multiplication with the reciprocal
        for (unsigned counter = 0; counter < m_unique_shared_node.size();
             ++counter)
          x[m_unique_shared_node[counter]]
            *= m_unique_shared_node_process_count[counter];
      }

      //! \brief Gather an accumulated vector from all processes.
      //!
      //! This method gathers an accumulated vector from all processes and
      //! stores the result in one huge vector.
      //! \param[in,out] x Input: the local accumulated vector; output: the
      //! global gathered vector.
      template <typename TVectorType>
      void allgatherAccumulatedVector(TVectorType& x)
      {
        // make a distributed vector and sum all distributed vectors
        distribute(x);
        allgatherDistributedVector(x);
      }

      //! \brief Gather a distributed vector from all processes.
      //!
      //! This method gathers a distributed vector from all processes and
      //! stores the result in one huge vector.
      //! \param[in,out] x Input: the local distributed vector; output: the
      //! global gathered vector.
      template <typename TVectorType>
      void allgatherDistributedVector(TVectorType& x)
      {
        // if this is the only process, we do not have to do anything
        if (Network::size() == 1)
          return;

        // fill our local values into a huge global vector
        TVectorType global_vector(m_num_global_nodes, 0);
        for (unsigned counter = 0; counter < x.size(); ++counter)
          global_vector[m_local_to_global_map[counter]] = x[counter];
        // sum all global vectors
        collect(global_vector);
        // copy the result
        x = global_vector;
      }

    private:
      //! \brief The total amount of global nodes.
      unsigned m_num_global_nodes;

      //! \brief Vector containing the local -> global node index map.
      std::vector<TIndexType> m_local_to_global_map;

      //! \brief Vector containing the amount of shared nodes.
      //!
      //! The i-th value in this vector is the amount of nodes that the current
      //! process shares with the i-th process.
      std::vector<unsigned> m_shared_node_count;

      //! \brief Vector containing the offset for shared nodes.
      //!
      //! The i-th value in this vector is the offset of the first node that
      //! is shared with the i-th process in the vector \p m_shared_node.
      std::vector<unsigned> m_shared_node_offset;

      //! \brief Vector containing the shared local node indices.
      //!
      //! This vector contains the local indices of nodes that this process
      //! shares with other processes. First the indices of nodes shared with
      //! process 0 are stored, followed by those shared with process 1...
      //! The first shared node of every process can be determined using
      //! the vector \p m_shared_node_offset and the amount of nodes shared
      //! with each process is stored in \p m_shared_node_count.
      std::vector<TIndexType> m_shared_node;

      //! \brief Vector holding the sorted shared local node indices.
      //!
      //! This vector is nothing but the sorted and "unique-ified" version of
      //! \p m_shared_node. This vector is needed for the \p distribute method
      //! as we need to traverse the shared nodes only once in there.
      std::vector<TIndexType> m_unique_shared_node;

      //! \brief Amount of processes a shared node belongs to.
      //!
      //! This vector stores the amount of processes to which a shared node in
      //! \p m_unique_shared_node belongs to. The value is stored in reciprocal
      //! form (i.e. each entry contains 1/(number of processes) as this is
      //! more convenient for multiplication purposes.
      std::vector<TRealType> m_unique_shared_node_process_count;

      //! \brief Buffer for receiving over the network.
      std::vector<unsigned char> m_receive_buffer;

      //! \brief Buffer for sending over the network.
      std::vector<unsigned char> m_send_buffer;
  };

} // namespace agile

#endif // AGILE_NETWORK_COMMUNICATOR_HPP

// End of $Id: communicator.hpp 476 2011-06-16 08:54:14Z freiberger $.
