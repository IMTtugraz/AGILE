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

// $Id: gpu_communicator.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_NETWORK_GPU_COMMUNICATOR_HPP
#define AGILE_NETWORK_GPU_COMMUNICATOR_HPP

// This code is based on previous work done by Manfred Liebmann.
// Clean-up, some re-structuring and GPU code by Manuel Freiberger.

#include <vector>
#include <algorithm>
#include <numeric>

#include "agile/gpu_config.hpp"
#include "agile/network/network.hpp"
#include "agile/gpu_environment.hpp"
#include "agile/gpu_vector.hpp"
#include "agile/gpu_host_vector.hpp"
#include "agile/radix_exchange_sort.hpp"
#include "agile/exception.hpp"

#include <cuda.h>

namespace agile
{
  // forward declarations
  template <typename TType1, typename TType2>
  void extractSharedNodes(const GPUVector<TType1>& map_shared_to_vector,
                          const GPUVector<TType2>& v, void* shared);
  template <typename TType1, typename TType2>
  void insertSharedNodes(const GPUVector<TType1>& map_vector_to_shared,
                         const void* shared, unsigned num_shared_indices,
                         GPUVector<TType2>& v);
  template <typename TType1, typename TType2, typename TType3>
  void multiplySharedNodes(
    const GPUVector<TType1>& shared_node_indices,
    const GPUVector<TType1>& vector_to_shared_node_indices,
    const GPUVector<TType2>& shared_node_multiplicand,
    void* gpu_buffer, GPUVector<TType3>& v);

  //! \brief Class for parallel communication.
  //!
  //! This class takes the global node numbers associated with the current
  //! process and establishes all the necessary communication with other
  //! processes.
  //! \p TIndexType is the type used for node indices. Usually this is
  //! either \p int or \p unsigned.
  //! \p TRealType is the type used to store the reciprocal of the amount
  //! of processes a shared node belongs to. Usually this is a \p double.
  //! \p TGPURealType is the type used to store the reciprocal of the amount
  //! of processes a shared node belongs to on the GPU.
  //! \author Manfred Liebmann, Manuel Freiberger
  template <typename TIndexType, typename TRealType, typename TGPURealType>
  class GPUCommunicator
  {
    public:
      //! \brief Default constructor.
      //!
      //! The communicator is not usable before \p setLocal2GlobalMap was
      //! called.
      GPUCommunicator()
        : m_gpu_part_initialized(false)
      {
      }

      //! \brief Constructor.
      //!
      //! This constructor takes the global node numbers associated with the
      //! current process in vector form. This vector also establishes a map
      //! between the local node numbering and the global numbering (local
      //! node 0 is mapped to global[0], local node 1 to global[1] and so on).
      GPUCommunicator(const std::vector<TIndexType>& global)
        : m_gpu_part_initialized(false)
      {
        setLocal2GlobalMap(global);
      }

      //! \brief Allocate the best GPU that is available.
      //!
      //! This method loops over all GPUs in the current PC and communicates
      //! with other processes on this PC. Every process gets its own GPU.
      //! Processes with lower ranks can choose first. Currently the number
      //! of multiprocessors is used as quality criterion.
      //! \return True, if successful.
      bool allocateGPU()
      {
        // from the device name we determine the number of PCs involved
        std::string processor_name = Network::processor_name();
        // communicate the names
        std::vector<std::string> all_processor_names;
        Network::allgather(processor_name, all_processor_names);
        // have a look how oftern our name occurs in all the names
        // we store the amount of processes on the same PC that have a higher
        // rank than we
        unsigned higher_rank_processors_with_same_name = 0;
        // and we also store the ranks of all the processes on the same PC
        std::vector<unsigned> share_name_with_process;
        for (unsigned counter = 0; counter < Network::size(); ++counter)
          if (all_processor_names[counter] == processor_name)
          {
            if (counter != Network::rank())
              share_name_with_process.push_back(counter);
            if (counter < Network::rank())
              ++higher_rank_processors_with_same_name;
          }
        // free memory for names
        all_processor_names.clear();

        // determine the amount of GPUs in this PC
        unsigned num_gpus = GPUEnvironment::getNumGPUs();
        // there have to be more (or at least the same number of) GPUs than
        // processes
        if (num_gpus < share_name_with_process.size() + 1)
          return false;

        // a vector holding the GPU number assigned to each process
        std::vector<unsigned> assigned_gpus(Network::size(), unsigned(-1));
        // loop until all processors have got a GPU
        unsigned our_gpu = unsigned(-1);
        unsigned iteration = 0;
        while (std::find(assigned_gpus.begin(), assigned_gpus.end(), -1)
               != assigned_gpus.end())
        {
          // if all processes with a higher rank have got their GPU already,
          // it is our turn now
          if (iteration == higher_rank_processors_with_same_name)
          {
            // loop through all the GPUs on this PC
            unsigned max_num_multiprocessors = 0;
            unsigned gpu_with_max_num_multiprocessors = unsigned(-1);
            for (unsigned gpu_counter = 0; gpu_counter < num_gpus;
                 ++gpu_counter)
            {
              // if another process on this PC has already taken this GPU,
              // we neglect it
              bool gpu_already_taken = false;
              for (unsigned counter = 0;
                   counter < share_name_with_process.size(); ++counter)
                if (assigned_gpus[share_name_with_process[counter]]
                    == gpu_counter)
                {
                  gpu_already_taken = true;
                  break;
                }
              if (gpu_already_taken)
                continue;

              // this GPU is still available -> check the properties
              cudaDeviceProp properties;
              cudaGetDeviceProperties(&properties, gpu_counter);
              // look for the one GPU with the maximum amount of multiprocessors
              // at this point we could also implement different quality
              // criteria like the maximum amount of memory, for example
              if (properties.multiProcessorCount > int(max_num_multiprocessors))
              {
                max_num_multiprocessors = properties.multiProcessorCount;
                gpu_with_max_num_multiprocessors = gpu_counter;
              }
            }
            // right now we should have found a good GPU
            if (gpu_with_max_num_multiprocessors == unsigned(-1))
              return false;
            // this GPU is ours!!
            our_gpu = gpu_with_max_num_multiprocessors;
          }

          // communicate the assigned GPUs
          Network::allgather(our_gpu, assigned_gpus);
          ++iteration;
        }

        // set the GPU
        GPUEnvironment::allocateGPU(our_gpu);
        return true;
      }

      //! \brief Set the local -> global mapping.
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

      //! \brief Set the local -> global mapping.
      //!
      //! Use this method to supply the global nodes which have to be handled
      //! by the current process as an iterator range.
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

        // ----- set up the GPU vectors -----

        // if no GPU was initialized, we do not set up the GPU vectors
        if (!GPUEnvironment::initialized())
          return;

        m_shared_node_gpu.assignFromHost(
          m_shared_node.begin(), m_shared_node.end());
        m_unique_shared_node_gpu.assignFromHost(
          m_unique_shared_node.begin(), m_unique_shared_node.end());
        m_unique_shared_node_process_count_gpu.assignFromHost(
          m_unique_shared_node_process_count.begin(),
          m_unique_shared_node_process_count.end());

        // we also need the inverse mapping of m_shared_node_gpu; determine the
        // last index that is shared, make a vector of that size and construct
        // the inverse mapping
        TIndexType max_shared_index = *std::max_element(m_shared_node.begin(),
                                                        m_shared_node.end());
        std::vector<TIndexType> vector_to_shared_map(
          max_shared_index + 1, max_shared_index + 1);
        for (unsigned counter = 0; counter < m_shared_node.size(); ++counter)
          vector_to_shared_map[m_shared_node[counter]] = counter;
        m_map_vector_to_shared_node_gpu.assignFromHost(
          vector_to_shared_map.begin(), vector_to_shared_map.end());

        // create the same inverse mapping from m_unique_shared_node_gpu
        max_shared_index = *std::max_element(m_unique_shared_node.begin(),
                                             m_unique_shared_node.end());
        vector_to_shared_map.assign(max_shared_index + 1, max_shared_index + 1);
        for (unsigned counter = 0; counter < m_unique_shared_node.size();
             ++counter)
          vector_to_shared_map[m_unique_shared_node[counter]] = counter;
        m_map_vector_to_unique_shared_node_gpu.assignFromHost(
          vector_to_shared_map.begin(), vector_to_shared_map.end());

        // remember that the GPU part was initialized
        m_gpu_part_initialized = true;
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
      //! \param[in,out] vector Input: The vector of this process to sum.
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

      //! \brief Calculate the global sum of the input vectors.
      //!
      //! This method communicates the input vectors to all processes and
      //! builds their sum.
      //! \param[in,out] vector Input: The GPU vector from this process.
      //! Output: The sum of all vectors.
      template <typename TType>
      void collect(GPUVector<TType>& vector)
      {
        // if this is the only process, we do not have to do anything
        if (Network::size() == 1)
          return;

        // the GPU part has to be initialized
        AGILE_ASSERT(m_gpu_part_initialized,
                      StandardException::ExceptionMessage(
                        "The GPU part of this communicator was not "
                        "initialized."));

        typedef typename GPUVector<TType>::value_type value_type;
        unsigned vector_size = vector.size();
        unsigned network_size = Network::size();
        // first copy the GPUVector to the host (the vector will be resized
        // automatically)
        vector.copyToHost(m_gpu_host_buffer);
        // resize the receive buffer so that it can hold the values of all
        // vectors
        m_receive_buffer.resize(
          sizeof(value_type) * vector_size * network_size);
        // communicate the vectors
        Network::allgather(
          &m_gpu_host_buffer[0], vector_size * sizeof(value_type),
          &m_receive_buffer[0], vector_size * sizeof(value_type));

        // for every element of the vector sum the values accross the processes
        value_type* ptr = (value_type*)(&m_receive_buffer[0]);
        for (unsigned element_counter = 0; element_counter < vector_size;
             ++element_counter)
        {
          value_type sum = 0;
          for (unsigned counter = 0; counter < network_size; ++counter)
            sum += ptr[counter * vector_size + element_counter];
          m_gpu_host_buffer[element_counter] = sum;
        }

        // upload the vector to the GPU again
        vector.assignFromHost(m_gpu_host_buffer.begin(),
                              m_gpu_host_buffer.end());
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

      //! \brief Accumulation of a GPU vector.
      //!
      //! \param[in,out] x in: Distributed vector, out: Accumulated vector.
      template <typename TType>
      void accumulate(GPUVector<TType>& x)
      {
        // if this is the only process, we do not have to do anything
        if (Network::size() == 1)
          return;

        AGILE_ASSERT(m_gpu_part_initialized,
                      StandardException::ExceptionMessage(
                        "The GPU part of this communicator was not "
                        "initialized."));

        // extract the shared elements from the GPU vector
        unsigned byte_size = m_shared_node_gpu.size() * sizeof(TType);
        m_gpu_buffer.resize(byte_size);
        m_gpu_host_buffer.resize(byte_size);
        m_receive_buffer.resize(byte_size);
        extractSharedNodes(m_shared_node_gpu, x, m_gpu_buffer.data());
        // copy them to the host
        CUDA_SAFE_CALL(cudaMemcpy(&m_gpu_host_buffer[0], m_gpu_buffer.data(),
                                  byte_size, cudaMemcpyDeviceToHost));
        // communicate the shared nodes
        std::vector<unsigned> block_size(m_shared_node_count.size());
        std::vector<unsigned> block_offset(m_shared_node_count.size());
        for(unsigned counter = 0; counter < m_shared_node_count.size();
            ++counter)
        {
          block_size[counter] = m_shared_node_count[counter] * sizeof(TType);
          block_offset[counter] = m_shared_node_offset[counter] * sizeof(TType);
        }
        Network::alltoallv(&m_gpu_host_buffer[0], (int*)&block_size[0],
                           (int*)&block_offset[0], &m_receive_buffer[0],
                           (int*)&block_size[0], (int*)&block_offset[0]);
        // sum the individual contributions
        TType* src_ptr = (TType*)(&m_receive_buffer[0]);
        TType* dest_ptr = (TType*)(&m_gpu_host_buffer[0]);
        for (unsigned counter = 0; counter < m_shared_node_gpu.size();
             ++counter)
          dest_ptr[counter] += src_ptr[counter];
        // transfer back to the GPU
        CUDA_SAFE_CALL(cudaMemcpy(m_gpu_buffer.data(), &m_gpu_host_buffer[0],
                                  byte_size, cudaMemcpyDeviceToHost));
        // and insert into the original vector
        insertSharedNodes(m_map_vector_to_shared_node_gpu,
                          m_gpu_buffer.data(), m_shared_node_gpu.size(), x);
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

      //! \brief Distribution of a GPU vector.
      //!
      //! \param[in,out] x in: Accumulated vector, out: Distributed vector.
      template <typename TType>
      void distribute(GPUVector<TType>& x)
      {
        // if this is the only process, we do not have to do anything
        if (Network::size() == 1)
          return;

        AGILE_ASSERT(m_gpu_part_initialized,
                      StandardException::ExceptionMessage(
                        "The GPU part of this communicator was not "
                        "initialized."));

        unsigned byte_size = m_unique_shared_node_gpu.size() * sizeof(TType);
        m_gpu_buffer.resize(byte_size);
        multiplySharedNodes(m_unique_shared_node_gpu,
                            m_map_vector_to_unique_shared_node_gpu,
                            m_unique_shared_node_process_count_gpu,
                            m_gpu_buffer.data(), x);
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

      //! \brief A buffer on GPU memory to extract/insert vector elements.
      GPUVector<unsigned char> m_gpu_buffer;

      //! \brief A host vector for fast GPU <-> host communication.
      GPUHostVector<unsigned char> m_gpu_host_buffer;

      //! \brief A GPU vector holding the shared local node indices.
      //!
      //! This is a GPU vector holding the content of m_shared_node in GPU
      //! memory. It defines a mapping from the shared nodes of this process to
      //! a local vector.
      GPUVector<TIndexType> m_shared_node_gpu;

      //! \brief A GPU vector holding m_unique_shared_node.
      GPUVector<TIndexType> m_unique_shared_node_gpu;

      //! \brief A GPU vector holding m_unique_shared_node_process_count.
      GPUVector<TGPURealType> m_unique_shared_node_process_count_gpu;

      //! \brief A GPU vector holding a map from a vector to shared nodes.
      //!
      //! This GPU vector defines to which shared node a local vector index
      //! belongs to. If a vector element is not to be shared, the destination
      //! index in this mapping has to be larger than the size of this
      //! vector. This is more or less the inverse mapping of m_shared_node_gpu.
      GPUVector<TIndexType> m_map_vector_to_shared_node_gpu;

      //! \brief The inverse mapping of m_unique_shared_node_gpu.
      GPUVector<TIndexType> m_map_vector_to_unique_shared_node_gpu;

      //! \brief A flag to store if the GPU part was initialized.
      bool m_gpu_part_initialized;
  };

} // namespace agile

#endif // AGILE_NETWORK_GPU_COMMUNICATOR_HPP

// End of $Id: gpu_communicator.hpp 476 2011-06-16 08:54:14Z freiberger $.
