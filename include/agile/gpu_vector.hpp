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

// $Id: gpu_vector.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_VECTOR_HPP
#define AGILE_GPU_VECTOR_HPP

#include "agile/gpu_config.hpp"
#include "agile/gpu_memory.hpp"
#include "agile/gpu_vector_base.hpp"

#include <complex>
#include <vector>
#include <cuda_runtime.h>

// in debug mode we need the iostream library for output
#ifdef AGILE_DEBUG
#include <iostream>
#endif

namespace agile
{
  //! \brief A vector on the GPU.
  template <typename TType>
  class GPUVector : public GPUVectorBase<GPUVector<TType> >
  {
    public:
      //! The type of the elements in this vector.
      typedef TType value_type;
      //! The iterator type.
      typedef TType* iterator;
      //! The const-iterator type.
      typedef const TType* const_iterator;

      //! \brief Default constructor.
      //!
      //! The default constructor creates an empty GPU vector which does not
      //! allocate any memory.
      GPUVector()
        : m_size(0), m_capacity(0), m_data(0)
      {
      }

      //! \brief Copy constructor.
      //!
      //! This copy constructor implements deep-copying of the GPU vector \p v.
      //! \param[in] v The GPU vector to be copied.
      GPUVector(const GPUVector& v)
        : m_size(v.size()), m_capacity(v.size()), m_data(0)
      {
        if (m_size)
        {
          allocateMemory((void**)&m_data, m_size * sizeof(TType));
          CUDA_SAFE_CALL(cudaMemcpy(m_data, v.data(), m_size * sizeof(TType),
                                    cudaMemcpyDeviceToDevice));
        }
      }

      //! \brief Construct a GPU vector with a given size and initial value.
      //!
      //! Create a GPU vector with \p size elements and set all elements to
      //! \p value.
      //! \param[in] size The amount of elements the new vector shall contain.
      //! \param[in] value The value used for the vector elements.
      GPUVector(unsigned size, const TType& value = TType())
        : m_size(size), m_capacity(size), m_data(0)
      {
        if (m_size)
        {
          allocateMemory((void**)&m_data, m_size * sizeof(TType));
          lowlevel::setVectorConstant(value, m_data, m_size);
        }
      }

      //! \brief Construct from an iterator range.
      //!
      //! This constructor creates a new GPU vector by copying the range
      //! specified by [begin, end). The iterators have to point to GPU
      //! memory which means that it is illegal to copy from the host to
      //! the device using this constructor.
      //! \param[in] begin Iterator to the start of the vector to be copied.
      //! \param[in] end Past-the-end iterator of the vector to be copied.
      GPUVector(const const_iterator& begin, const const_iterator& end)
        : m_size(end - begin), m_capacity(end - begin), m_data(0)
      {
        if (m_size)
        {
          allocateMemory((void**)&m_data, m_size * sizeof(TType));
          // copy from device to the device
          CUDA_SAFE_CALL(cudaMemcpy(m_data, &(*begin), m_size * sizeof(TType),
                                    cudaMemcpyDeviceToDevice));
        }
      }

      //! \brief Destructor.
      ~GPUVector()
      {
        freeMemory();
      }

      //! \brief Assignment operator.
      //!
      //! Performs a deep-copy of the GPU vector \p rhs.
      //! \param[in] rhs The GPU vector to be copied.
      //! \return A reference to the current GPU vector.
      GPUVector& operator= (const GPUVector& rhs)
      {
        if (rhs.size() == 0)
        {
          m_size = 0;
          return *this;
        }

        // check if the new vector fits in the old memory block; if this is not
        // the case, free the old memory and allocate a new block
        if (rhs.size() > m_capacity)
        {
          freeMemory();
          allocateMemory((void**)&m_data, rhs.size() * sizeof(TType));
          m_capacity = rhs.size();
        }

        m_size = rhs.size();
        // copy the data
        CUDA_SAFE_CALL(cudaMemcpy(m_data, rhs.data(), m_size * sizeof(TType),
                                  cudaMemcpyDeviceToDevice));

        return *this;
      }

      //! \brief Assign \p size copies of \p value to the vector.
      //!
      //! The vector will be resized to hold \p size elements and all elements
      //! are initialized to \p value.
      //! \param[in] size The new vector size.
      //! \param[in] value The value to which all vector elements are set.
      void assign(unsigned size, const TType& value)
      {
        if (size == 0)
        {
          m_size = 0;
          return;
        }

        // check if the new vector fits in the old memory block; if this is not
        // the case, free the old memory and allocate a new block
        if (size > m_capacity)
        {
          freeMemory();
          allocateMemory((void**)&m_data, size * sizeof(TType));
          m_capacity = size;
        }

        m_size = size;
        lowlevel::setVectorConstant(value, m_data, m_size);
      }

      //! \brief Assign a vector on the host to a GPU vector.
      //!
      //! The iterators have to point to a consecutive memory area such that
      //! memcpy can be used to copy the content. This means that is illegal
      //! to use list iterators, for example.
      //! \param[in] begin Iterator to the start of the vector to be copied.
      //! \param[in] end Past-the-end iterator of the vector to be copied.
      template <typename TIterator>
      void assignFromHost(const TIterator& begin, const TIterator& end)
      {
        unsigned size = end - begin;
        // a vector of size 0 is quite easy to handle
        if (size == 0)
        {
          m_size = 0;
          return;
        }

        // check if the new vector fits in the old memory block; if this is not
        // the case, free the old memory and allocate a new block
        if (size > m_capacity)
        {
          freeMemory();
          allocateMemory((void**)&m_data, size * sizeof(TType));
          m_capacity = size;
        }

        // at this point we have a data block that is large enough to hold
        // the new vector so perform a memcopy from the host to the device
        CUDA_SAFE_CALL(cudaMemcpy(m_data, &(*begin), size * sizeof(TType),
                                  cudaMemcpyHostToDevice));
        m_size = size;
      }

      //! \brief Get an iterator to the first vector element.
      iterator begin()
      {
        return m_data;
      }

      //! \brief Get a const-iterator to the first vector element.
      const_iterator begin() const
      {
        return m_data;
      }

      //! \brief Get the maximum numer of elements the vector can hold.
      //!
      //! This method returns the maximum amount of elements that the vector
      //! can hold before it needs to allocate new GPU memory.
      //! \return The maxmimum size of the GPU vector that fits in the currently
      //! allocated GPU memory.
      unsigned capacity() const
      {
        return m_capacity;
      }

      //! \brief Remove all elements from the vector.
      void clear()
      {
        m_size = 0;
      }

      //! \brief Copy the vector from the GPU memory to the host.
      //!
      //! This method copies the GPU memory holding the current GPU vector
      //! back to the host. The host vector is resized such that its byte-size
      //! is large enough to hold the GPU memory allocated by the GPU vector.
      //! This means that if you copy a GPUVector<float> to a
      //! std::vector<double>, for example, the host vector will be resized
      //! to half of the size of the GPU vector because sizeof(double) ==
      //! 2*sizeof(float). Then two floats from the GPU are copied into
      //! one double of the host meaning that you have to explicitly cast
      //! the host vector in order to get the content of the GPU vector.
      //! Therefore it is recommended to perform copies between the same
      //! GPU vector and host vector types only.
      //! \param[out] host_vector The host vector into which the GPU vector
      //! will be copied.
      template <typename THostVector>
      void copyToHost(THostVector& host_vector) const
      {
        typedef typename THostVector::value_type host_type;
        unsigned byte_size = m_size * sizeof(TType);
        host_vector.resize((byte_size + sizeof(host_type) - 1)
                           / sizeof(host_type));
        CUDA_SAFE_CALL(cudaMemcpy(&host_vector[0],
                                  m_data, m_size * sizeof(TType),
                                  cudaMemcpyDeviceToHost));
      }

      //! \brief Get a pointer to the data.
      //!
      //! This function can be seen as replacement for operator[].
      //! \return A pointer to the beginning of the GPU vector's memory.
      TType* data()
      {
        return m_data;
      }

      //! \brief Get a constant pointer to the data.
      //!
      //! This function can be seen as replacement for operator[].
      //! \return A pointer to the beginning of the GPU vector's memory.
      const TType* data() const
      {
        return m_data;
      }

      //! \brief Return true if the vector has no elements.
      bool empty() const
      {
        return m_size == 0;
      }

      //! \brief Get an iterator pointing past the last vector element.
      iterator end()
      {
        return m_data + m_size;
      }

      //! \brief Get a const-iterator pointing past the last vector element.
      const_iterator end() const
      {
        return m_data + m_size;
      }

      //! \brief Resize the vector.
      //!
      //! This method resizes the GPU vector to \p size elements. If the vector
      //! needs to allocate new memory (which happens if \p size >
      //! \p capacity()), these new elements are initialized to \p value.
      //! No initialization is done for elements with index
      //! [size(), capacity()).
      //! \param[in] size The new vector size.
      //! \param[in] value The value for newly created elements.
      //! \sa assign
      void resize(unsigned size, const TType& value = TType())
      {
        // if the vector is still large enough, we can simply change the size
        // and everything is done
        if (size <= m_capacity)
        {
          m_size = size;
          return;
        }

        // hm... the new vector does not fit into the allocated memory
        if (m_data)
        {
          // allocate more memory and copy the old data
          TType* new_data;
          allocateMemory((void**)&new_data, size * sizeof(TType));
          CUDA_SAFE_CALL(cudaMemcpy(new_data, m_data, m_size * sizeof(TType),
                                    cudaMemcpyDeviceToDevice));
          // the old space can be freed now
          freeMemory();
          m_data = new_data;
        }
        else
        {
          // we only need a new data block
          allocateMemory((void**)&m_data, size * sizeof(TType));
        }

        // the new elements have to be initialised correctly
        lowlevel::setVectorConstant(value, m_data + m_size, size - m_size);
        m_size = size;
        m_capacity = size;
      }

      //! \brief Get the number of vector elements.
      unsigned size() const
      {
        return m_size;
      }

    private:
      //! \brief Number of vector elements.
      unsigned m_size;

      //! \brief Amount of elements that can be stored in the GPU memory.
      unsigned m_capacity;

      //! \brief Pointer to GPU memory containing the elements.
      TType* m_data;

      //! \brief Allocate GPU memory
      inline void allocateMemory(void **devPtr, size_t count)
      {
        /*
#ifdef AGILE_DEBUG
        std::cout << "GPUVector mem alloc: " << (count / sizeof(TType))
          << " elements" << std::endl << "                     " << count
          << " [Byte] (" << (count / 1024.0 / 1024.0) << " [MB])" << std::endl;
#endif
*/
        CUDA_SAFE_CALL(cudaMalloc(devPtr, count));
      }

      //! \brief Free GPU memory if allocated.
      void freeMemory()
      {
        if (m_data)
        {
          /*
#ifdef AGILE_DEBUG
          std::cout << "GPUVector mem free:  " << m_capacity << " elements"
            << std::endl << "                     "
            << (m_capacity * sizeof(TType)) << " [Byte] ("
            << (m_capacity * sizeof(TType) / 1024.0 / 1024.0)
            << " [MB])" << std::endl;
#endif
*/
          CUDA_SAFE_CALL(cudaFree(m_data));
          m_data = 0;
        }
      }
  };

  //CUBLAS functions
  //
  
  template <typename TType>
  void copy(const GPUVector<TType>& x, GPUVector<TType>& y);

  template <typename TType>
  void maxElement(const GPUVector<TType>& x, int* maxVal);

} // namespace agile

#endif // AGILE_GPU_VECTOR_HPP

// End of $Id: gpu_vector.hpp 476 2011-06-16 08:54:14Z freiberger $.
