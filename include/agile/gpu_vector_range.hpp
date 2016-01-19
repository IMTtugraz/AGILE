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

// $Id: gpu_vector_range.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_VECTOR_RANGE_HPP
#define AGILE_GPU_VECTOR_RANGE_HPP

#include "agile/gpu_config.hpp"
#include "agile/gpu_vector_base.hpp"

#include <complex>
#include <vector>
#include <cuda_runtime.h>

namespace agile
{
  //! \brief A range of a GPU vector.
  //!
  //! The GPU vector range is a contiguous block of elements of a \p GPUVector
  //! and is used to address sub-vectors, for example. The vector range is
  //! valid as long as the memory location of the corresponding GPU vector
  //! does not change. This means the vector range loses validity, when the
  //! vector it is attached to is either deleted or resized. There is no
  //! mechanism implemented to check these cases.
  template <typename TType>
  class GPUVectorRange : public GPUVectorBase<GPUVectorRange<TType> >
  {
    private:
      //! Hidden default constructor.
      GPUVectorRange();

      //! Hidden copy constructor.
      GPUVectorRange(const GPUVectorRange&);

      //! Hidden assignment operator.
      GPUVectorRange& operator= (const GPUVectorRange&);

    public:
      //! The type of the elements in this vector.
      typedef TType value_type;
      //! The iterator type.
      typedef TType* iterator;
      //! The const-iterator type.
      typedef const TType* const_iterator;

      //! \brief Default constructor.
      //!
      //! This constructor creates a vector range taking \p size elements from
      //! \p vector starting with the \p start-th element.
      GPUVectorRange(GPUVector<TType>& vector, unsigned start, unsigned size)
        : m_data(vector.m_data + start), m_size(size)
      {
      }

      //! \brief Construction from a pointer and a size.
      //!
      //! This method constructs a \p GPUVectorRange from a pointer in GPU
      //! memory and the size of the memory area. This enables the user to
      //! treat any memory area as \p GPUVector with the additional flexibility
      //! of its own memory management.
      GPUVectorRange(TType* data, unsigned size)
        : m_data(data), m_size(size)
      {
      }

#if 0
TODO: is this method needed for the vector range?
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
#endif

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

      //! \brief Copy the vector range from the GPU memory to the host.
      //!
      //! This method copies the GPU memory holding the current GPU vector range
      //! back to the host. The host vector is resized such that its byte-size
      //! is large enough to hold the GPU memory of the GPU vector range.
      //! \see GPUVector::copyToHost
      //! \param[out] host_vector The host vector into which the GPU vector
      //! range will be copied.
      template <typename THostVector>
      void copyToHost(THostVector& host_vector) const
      {
        typedef typename THostVector::value_type host_type;
        size_t byte_size = m_size * sizeof(TType);
        host_vector.resize((byte_size + sizeof(host_type) - 1)
                           / sizeof(host_type));
        CUDA_SAFE_CALL(cudaMemcpy(&host_vector[0],
                                  m_data, m_size * sizeof(TType),
                                  cudaMemcpyDeviceToHost));
      }

      //! \brief Get a pointer to the data.
      //!
      //! This function can be seen as replacement for operator[].
      //! \return A pointer to the beginning of the memory area the GPU
      //! vector range points to.
      TType* data()
      {
        return m_data;
      }

      //! \brief Get a constant pointer to the data.
      //!
      //! This function can be seen as replacement for operator[].
      //! \return A pointer to the beginning of the memory area the GPU
      //! vector range points to.
      const TType* data() const
      {
        return m_data;
      }

      //! \brief Return true if the vector range has no elements.
      bool empty() const
      {
        return m_size == 0;
      }

      //! \brief Get an iterator pointing past the last range element.
      iterator end()
      {
        return m_data + m_size;
      }

      //! \brief Get a const-iterator pointing past the last range element.
      const_iterator end() const
      {
        return m_data + m_size;
      }

      //! \brief Get the number of vector elements in the range.
      unsigned size() const
      {
        return m_size;
      }

    private:
      //! \brief Pointer to GPU memory containing the elements.
      TType* m_data;

      //! \brief Number of vector elements.
      unsigned m_size;
  };

} // namespace agile

#endif // AGILE_GPU_VECTOR_RANGE_HPP

// End of $Id: gpu_vector_range.hpp 476 2011-06-16 08:54:14Z freiberger $.
