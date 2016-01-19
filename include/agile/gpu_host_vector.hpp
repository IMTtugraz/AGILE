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

// $Id: gpu_host_vector.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_HOST_VECTOR_HPP
#define AGILE_GPU_HOST_VECTOR_HPP

#include <string.h>

#include <vector>
#include <cuda_runtime.h>

#include "agile/gpu_config.hpp"

namespace agile
{
  //! \brief A GPU vector using pinned memory in the host's memory area.
  template <typename TType>
  class GPUHostVector
  {
    public:
      //! The type of the elements in this vector.
      typedef TType value_type;
      //! The iterator type.
      typedef TType* iterator;
      //! The const-iterator type.
      typedef const TType* const_iterator;

      //! \brief Constructor.
      //!
      //! The default constructor creates an empty GPU host vector which does
      //! not allocate any memory.
      GPUHostVector()
        : m_size(0), m_capacity(0), m_data(0)
      {
      }

      //! \brief Copy constructor.
      //!
      //! This copy constructor implements deep-copying of the GPU vector \p v.
      //! \param[in] v The GPU vector to be copied.
      GPUHostVector(const GPUHostVector& v)
        : m_size(v.size()), m_capacity(v.size()), m_data(0)
      {
        if (m_size)
        {
          CUDA_SAFE_CALL(cudaMallocHost((void**)&m_data,
                                        m_size * sizeof(TType)));
          memcpy(m_data, &v[0], m_size * sizeof(TType));
        }
      }

      //! \brief Construct with a given size and initial value.
      //!
      //! Create a GPU host vector with \p size elements and set all elements
      //! to \p value.
      //! \param[in] size The amount of elements the new vector shall contain.
      //! \param[in] value The value used for the vector elements.
      GPUHostVector(unsigned size, const TType& value = TType())
        : m_size(size), m_capacity(size), m_data(0)
      {
        if (m_size)
        {
          CUDA_SAFE_CALL(cudaMallocHost((void**)&m_data,
                                        m_size * sizeof(TType)));
          initialize(m_data, m_size, value);
        }
      }

      //! \brief Construct from an iterator range.
      //!
      //! This constructor creates a new GPU host vector by copying the range
      //! specified by [begin, end). The iterators have to point to host
      //! memory which means that it is illegal to copy from the device to
      //! the host using this constructor.
      //! \param[in] begin Iterator to the start of the vector to be copied.
      //! \param[in] end Past-the-end iterator of the vector to be copied.
      GPUHostVector(const const_iterator& begin, const const_iterator& end)
        : m_size(end - begin), m_capacity(end - begin), m_data(0)
      {
        if (m_size)
        {
          CUDA_SAFE_CALL(cudaMallocHost((void**)&m_data,
                                        m_size * sizeof(TType)));
          memcpy(m_data, &(*begin), m_size * sizeof(TType));
        }
      }

      //! \brief Destructor.
      ~GPUHostVector()
      {
        freeMemory();
      }

      //! \brief Assignment operator.
      //!
      //! Performs a deep-copy of the GPU host vector \p rhs.
      //! \param[in] rhs The GPU host vector to be copied.
      //! \return A reference to the current GPU host vector.
      GPUHostVector& operator= (const GPUHostVector& rhs)
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
          CUDA_SAFE_CALL(cudaMallocHost((void**)&m_data,
                                        rhs.size() * sizeof(TType)));
          m_capacity = rhs.size();
        }

        m_size = rhs.size();
        // copy the data
        memcpy(m_data, &rhs[0], m_size * sizeof(TType));

        return *this;
      }

      //! \brief Access an element.
      TType& operator[](unsigned index)
      {
        return *(m_data + index);
      }

      //! \brief Access a constant element.
      const TType& operator[](unsigned index) const
      {
        return *(m_data + index);
      }

      //! \brief Assign \p size copies of \p value to the vector.
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
          CUDA_SAFE_CALL(cudaMallocHost((void**)&m_data, size * sizeof(TType)));
          m_capacity = size;
        }

        m_size = size;
        initialize(m_data, m_size, value);
      }

      //! \brief Assign the elements between \p start and \p end to the vector.
      template <typename TIterator>
      void assign(const TIterator& begin, const TIterator& end,
                  typename TIterator::value_type = 0)
      {
        resizeMemory(end - begin);
        TType* ptr = m_data;
        for (TIterator iter = begin; iter != end; ++iter)
          *ptr++ = *iter;
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
      //! can hold before it needs to allocate new memory.
      //! \return The maxmimum size of the GPU host vector that fits in the
      //! currently allocated memory.
      unsigned capacity() const
      {
        return m_capacity;
      }

      //! \brief Remove all elements from the vector.
      void clear()
      {
        m_size = 0;
      }

      //! \brief Return true if the vector has no elements.
      bool empty() const
      {
        return m_size == 0;
      }

      //! \brief Get an iterator pointing pas the last vector element.
      iterator end()
      {
        return m_data + m_size;
      }

      //! \brief Get a const-iterator pointing pas the last vector element.
      const_iterator end() const
      {
        return m_data + m_size;
      }

      //! \brief Resize the vector.
      //!
      //! This method resizes the GPU host vector to \p size elements. If the
      //! vector needs to allocate new memory (which happens if \p size >
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
          CUDA_SAFE_CALL(cudaMallocHost((void**)&new_data,
                                        size * sizeof(TType)));
          memcpy(new_data, m_data, m_size * sizeof(TType));
          // the old space can be freed now
          freeMemory();
          m_data = new_data;
        }
        else
        {
          // we only need a new data block
          CUDA_SAFE_CALL(cudaMallocHost((void**)&m_data, size * sizeof(TType)));
        }

        // the new elements have to be initialised correctly
        initialize(m_data + m_size, size - m_size, value);
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

      //! \brief Amount of elements that can be stored in the allocated memory.
      unsigned m_capacity;

      //! \brief Pointer to the host memory containing the elements.
      TType* m_data;

      //! \brief Free memory if it is allocated.
      void freeMemory()
      {
        if (m_data)
        {
          cudaFreeHost(m_data);
          m_data = 0;
        }
      }

      //! \brief Initialize memory.
      //!
      //! This method initializes the host memory located at \p data with
      //! size \p size to \p value.
      void initialize(TType* data, unsigned size, const TType& value)
      {
        for (unsigned counter = 0; counter < size; ++counter)
          data[counter] = value;
      }
  };

}  // namespace agile

#endif // AGILE_GPU_HOST_VECTOR_HPP

// End of $Id: gpu_host_vector.hpp 476 2011-06-16 08:54:14Z freiberger $.
