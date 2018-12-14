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

// $Id: cpu_vector.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_CPU_VECTOR_HPP
#define AGILE_CPU_VECTOR_HPP

#include "agile/cpu_vector_base.hpp"

#include <cstring>  // memcpy
#include <vector>

namespace agile
{
  //! \brief A vector on the CPU.
  template <typename TType>
  class CPUVector : public CPUVectorBase<CPUVector<TType> >
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
      //! The default constructor creates an empty CPU vector which does not
      //! allocate any memory.
      CPUVector()
      {
      }

      //! \brief Copy constructor.
      //!
      //! This copy constructor implements deep-copying of the CPU vector \p v.
      //! \param[in] v The CPU vector to be copied.
      CPUVector(const CPUVector& v)
        : m_data(v.m_data)
      {
      }

      //! \brief Construct a CPU vector with a given size and initial value.
      //!
      //! Create a CPU vector with \p size elements and set all elements to
      //! \p value.
      //! \param[in] size The amount of elements the new vector shall contain.
      //! \param[in] value The value used for the vector elements.
      CPUVector(unsigned size, const TType& value = TType())
        : m_data(size, value)
      {
      }

      //! \brief Construct from an iterator range.
      //!
      //! This constructor creates a new CPU vector by copying the range
      //! specified by [begin, end). The iterators have to point to CPU
      //! memory which means that it is illegal to copy from the host to
      //! the device using this constructor.
      //! \param[in] begin Iterator to the start of the vector to be copied.
      //! \param[in] end Past-the-end iterator of the vector to be copied.
      CPUVector(const const_iterator& begin, const const_iterator& end)
        : m_data(begin, end)
      {
      }

      //! \brief Destructor.
      ~CPUVector()
      {
      }

      //! \brief Assign \p size copies of \p value to the vector.
      //!
      //! The vector will be resized to hold \p size elements and all elements
      //! are initialized to \p value.
      //! \param[in] size The new vector size.
      //! \param[in] value The value to which all vector elements are set.
      void assign(unsigned size, const TType& value)
      {
        m_data.assign(size, value);
      }

      //! \brief Assign a vector on the host to a CPU vector.
      //!
      //! The iterators have to point to a consecutive memory area such that
      //! memcpy can be used to copy the content. This means that is illegal
      //! to use list iterators, for example.
      //! \param[in] begin Iterator to the start of the vector to be copied.
      //! \param[in] end Past-the-end iterator of the vector to be copied.
      template <typename TIterator>
      void assignFromHost(const TIterator& begin, const TIterator& end)
      {
        m_data.assign(begin, end);
      }

      //! \brief Get an iterator to the first vector element.
      iterator begin()
      {
        return (m_data.size() ? &m_data[0] : 0);
      }

      //! \brief Get a const-iterator to the first vector element.
      const_iterator begin() const
      {
        return (m_data.size() ? &m_data[0] : 0);
      }

      //! \brief Get the maximum numer of elements the vector can hold.
      //!
      //! This method returns the maximum amount of elements that the vector
      //! can hold before it needs to allocate new CPU memory.
      //! \return The maxmimum size of the CPU vector that fits in the currently
      //! allocated CPU memory.
      unsigned capacity() const
      {
        return m_data.capacity();
      }

      //! \brief Remove all elements from the vector.
      void clear()
      {
        m_data.clear();
      }

      //! \brief Copy the vector from the CPU memory to the host.
      //!
      //! This method copies the CPU memory holding the current CPU vector
      //! back to the host. The host vector is resized such that its byte-size
      //! is large enough to hold the CPU memory allocated by the CPU vector.
      //! This means that if you copy a CPUVector<float> to a
      //! std::vector<double>, for example, the host vector will be resized
      //! to half of the size of the CPU vector because sizeof(double) ==
      //! 2*sizeof(float). Then two floats from the CPU are copied into
      //! one double of the host meaning that you have to explicitly cast
      //! the host vector in order to get the content of the CPU vector.
      //! Therefore it is recommended to perform copies between the same
      //! CPU vector and host vector types only.
      //! \param[out] host_vector The host vector into which the CPU vector
      //! will be copied.
      template <typename THostVector>
      void copyToHost(THostVector& host_vector) const
      {
        typedef typename THostVector::value_type host_type;
        size_t byte_size = m_data.size() * sizeof(TType);
        host_vector.resize((byte_size + sizeof(host_type) - 1)
                           / sizeof(host_type));
        memcpy(&host_vector[0], &m_data[0], byte_size);
      }

      //! \brief Get a pointer to the data.
      //!
      //! This function can be seen as replacement for operator[].
      //! \return A pointer to the beginning of the CPU vector's memory.
      TType* data()
      {
        return (m_data.size() ? &m_data[0] : 0);
      }

      //! \brief Get a constant pointer to the data.
      //!
      //! This function can be seen as replacement for operator[].
      //! \return A pointer to the beginning of the CPU vector's memory.
      const TType* data() const
      {
        return (m_data.size() ? &m_data[0] : 0);
      }

      //! \brief Return true if the vector has no elements.
      bool empty() const
      {
        return m_data.empty();
      }

      //! \brief Get an iterator pointing past the last vector element.
      iterator end()
      {
        return (m_data.size() ? &m_data[0] + m_data.size() : 0);
      }

      //! \brief Get a const-iterator pointing past the last vector element.
      const_iterator end() const
      {
        return (m_data.size() ? &m_data[0] + m_data.size() : 0);
      }

      //! Append an element to the vector.
      void push_back(const TType& value)
      {
        m_data.push_back(value);
      }

      //! \brief Resize the vector.
      //!
      //! This method resizes the CPU vector to \p size elements. If the vector
      //! needs to allocate new memory (which happens if \p size >
      //! \p capacity()), these new elements are initialized to \p value.
      //! No initialization is done for elements with index
      //! [size(), capacity()).
      //! \param[in] size The new vector size.
      //! \param[in] value The value for newly created elements.
      //! \sa assign
      void resize(unsigned size, const TType& value = TType())
      {
        m_data.resize(size, value);
      }

      //! \brief Get the number of vector elements.
      unsigned size() const
      {
        return m_data.size();
      }

      //! \brief Assignment operator.
      //!
      //! Performs a deep-copy of the CPU vector \p rhs.
      //! \param[in] rhs The CPU vector to be copied.
      //! \return A reference to the current CPU vector.
      CPUVector& operator= (const CPUVector& rhs)
      {
        m_data = rhs.m_data;
        return *this;
      }

      //! Get a certain element from the vector.
      //!
      //! \param[in] index The index of the element to fetch.
      TType& operator[] (unsigned index)
      {
        return m_data[index];
      }

      //! Get a certain element from the vector.
      //!
      //! \param[in] index The index of the element to fetch.
      const TType& operator[] (unsigned index) const
      {
        return m_data[index];
      }

    private:
      //! The vector holding our elements.
      std::vector<TType> m_data;
  };

} // namespace agile

#endif // AGILE_CPU_VECTOR_HPP

// End of $Id: cpu_vector.hpp 476 2011-06-16 08:54:14Z freiberger $.
