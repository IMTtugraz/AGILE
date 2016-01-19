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

// $Id: radix_exchange_sort.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_NETWORK_RADIX_EXCHANGE_SORT_HPP
#define AGILE_NETWORK_RADIX_EXCHANGE_SORT_HPP

#include "agile/exception.hpp"

namespace agile
{
  using namespace StandardException;

  // forward declarations
  template <typename TIterator1, typename TIterator2,
            typename TValueType1, typename TValueType2>
  class RadixExchangeSortCopy;

  template <typename TIterator1, typename TIterator2, typename TIterator3,
            typename TValueType1, typename TValueType2, typename TValueType3>
  class RadixExchangeSortCopyCopy;

  //! \brief Radix exchange sort.
  template <typename TIterator,
            typename TValueType = typename TIterator::value_type>
  class RadixExchangeSort
  {
    public:
      //! \brief Sort a vector given by [begin, end) using radix exchange sort.
      static void sort(const TIterator& begin, const TIterator& end)
      {
        // vectors with a size lower than 2 do not have to be sorted
        if (end - begin < 2)
          return;

        doSort(begin, end, getMaxBit(begin, end));
      }

    private:
      //! \brief The actual sorting procedure.
      //!
      //! \param[in] begin Iterator to the first element to sort.
      //! \param[in] end Iterator pointing past-the-end of the element range
      //!            to sort.
      //! \param[in] bit Bit mask used for sorting. This shall be a positive
      //!            value with only a single bit set.
      static void doSort(const TIterator& begin, const TIterator& end,
                         TValueType bit)
      {
        TIterator left = begin;
        TIterator right = end;

        // the right iterator points past-the-end, so decrease it
        --right;

        while (left != right)
        {
          // advance the left iterator until it points to an element with a set
          // bit
          while ((left != right) && ((*left & bit) == 0))
            ++left;
          // decrease the right iterator until it points to an element with
          // a zero bit
          while ((left != right) && (*right & bit))
            --right;
          // swap the elements
          TValueType temp = *left;
          *left = *right;
          *right = temp;
        }

        // if all the elements had a zero-bit, we have to increase the right
        // iterator such that it points past-the-end again because otherwise
        // we would miss the last element in the next recursion
        if ((*right & bit) == 0)
          ++right;

        // advance to the next bit
        bit >>= 1;
        if (bit)
        {
          if (right - begin > 1)
            doSort(begin, right, bit);
          if (end - right > 1)
            doSort(right, end, bit);
        }
      }

      //! \brief Get the key length for radix sort.
      //!
      //! This function calculates 2^floor(ld(x)) for all input values and
      //! returns the largest result. Note that this function only works for
      //! positive values.
      static TValueType getMaxBit(const TIterator& begin, const TIterator& end)
      {
        // first or all bit patterns
        TValueType ored_bit_pattern = 0;
        for (TIterator iter = begin; iter != end; ++iter)
          ored_bit_pattern |= *iter;

        // shift the bit pattern to the right until it is zero and at the same
        // time shift the result to the left
        TValueType result = 1;
        while (ored_bit_pattern >>= 1)
          result <<= 1;

        return result;
      }

      // declare friends
      template <typename TIterator1, typename TIterator2,
                typename TValueType1, typename TValueType2>
      friend class RadixExchangeSortCopy;

      template <typename TIterator1, typename TIterator2, typename TIterator3,
                typename TValueType1, typename TValueType2,
                typename TValueType3>
      friend class RadixExchangeSortCopyCopy;
  };

  //! \brief Radix exchange sort with copying of a second vector.
  template <typename TIterator1, typename TIterator2 = TIterator1,
            typename TValueType1 = typename TIterator1::value_type,
            typename TValueType2 = typename TIterator2::value_type>
  class RadixExchangeSortCopy
  {
    public:
      //! \brief Sort a vector while copying a second one.
      //!
      //! This method sorts the elements in the vector given by
      //! [sort_begin, sort_end) and copies the corresponding elements in the
      //! vector [copy_begin, copy_end).
      //! The vectors have to have the same length!
      static void sort(const TIterator1& sort_begin, const TIterator1& sort_end,
                       const TIterator2& copy_begin, const TIterator2& copy_end)
      {
        AGILE_ASSERT(sort_end - sort_begin == copy_end - copy_begin,
                      ExceptionSizeMismatch("sort vector", "copy vector",
                                            sort_end - sort_begin,
                                            copy_end - copy_begin));

        // vectors with a size lower than 2 do not have to be sorted
        if (sort_end - sort_begin < 2)
          return;

        doSort(sort_begin, sort_end, copy_begin, copy_end,
               RadixExchangeSort<TIterator1, TValueType1>::getMaxBit(
                 sort_begin, sort_end));
      }

    private:
      //! \brief The actual sorting procedure.
      static void doSort(const TIterator1& sort_begin,
                         const TIterator1& sort_end,
                         const TIterator2& copy_begin,
                         const TIterator2& copy_end,
                         TValueType1 bit)
      {
        TIterator1 sort_left = sort_begin;
        TIterator1 sort_right = sort_end;
        TIterator2 copy_left = copy_begin;
        TIterator2 copy_right = copy_end;

        // the right iterator point past-the-end, so decrease them
        --sort_right;
        --copy_right;

        while (sort_left != sort_right)
        {
          // advance the left iterator until it points to an element with a set
          // bit
          while ((sort_left != sort_right) && ((*sort_left & bit) == 0))
          {
            ++sort_left;
            ++copy_left;
          }
          // decrease the right iterator until it points to an element with
          // a zero bit
          while ((sort_left != sort_right) && (*sort_right & bit))
          {
            --sort_right;
            --copy_right;
          }
          // swap the elements
          TValueType1 sort_temp = *sort_left;
          *sort_left = *sort_right;
          *sort_right = sort_temp;
          TValueType2 copy_temp = *copy_left;
          *copy_left = *copy_right;
          *copy_right = copy_temp;
        }

        // if all the elements had a zero-bit, we have to increase the right
        // iterator such that it points past-the-end again because otherwise
        // we would miss the last element in the next recursion
        if ((*sort_right & bit) == 0)
        {
          ++sort_right;
          ++copy_right;
        }

        // advance to the next bit
        bit >>= 1;
        if (bit)
        {
          if (sort_right - sort_begin > 1)
            doSort(sort_begin, sort_right, copy_begin, copy_right, bit);
          if (sort_end - sort_right > 1)
            doSort(sort_right, sort_end, copy_right, copy_end, bit);
        }
      }
  };

  //! \brief Radix exchange sort with copying of a two other vectors.
  template <typename TIterator1, typename TIterator2 = TIterator1,
            typename TIterator3 = TIterator2,
            typename TValueType1 = typename TIterator1::value_type,
            typename TValueType2 = typename TIterator2::value_type,
            typename TValueType3 = typename TIterator3::value_type>
  class RadixExchangeSortCopyCopy
  {
    public:
      //! \brief Sort a vector while copying two other ones.
      //!
      //! This method sorts the elements in the vector given by
      //! [sort_begin, sort_end) and copies the corresponding elements in the
      //! vectors [copy1_begin, copy1_end) and [copy2_begin, copy2_end).
      //! The vectors have to have the same length!
      static void sort(
        const TIterator1& sort_begin, const TIterator1& sort_end,
        const TIterator2& copy1_begin, const TIterator2& copy1_end,
        const TIterator3& copy2_begin, const TIterator3& copy2_end)
      {
        AGILE_ASSERT(sort_end - sort_begin == copy1_end - copy1_begin,
                      ExceptionSizeMismatch("sort vector", "1st copy vector",
                                            sort_end - sort_begin,
                                            copy1_end - copy1_begin));
        AGILE_ASSERT(sort_end - sort_begin == copy2_end - copy2_begin,
                      ExceptionSizeMismatch("sort vector", "2nd copy vector",
                                            sort_end - sort_begin,
                                            copy2_end - copy2_begin));

        // vectors with a size lower than 2 do not have to be sorted
        if (sort_end - sort_begin < 2)
          return;

        doSort(sort_begin, sort_end, copy1_begin, copy1_end,
               copy2_begin, copy2_end,
               RadixExchangeSort<TIterator1, TValueType1>::getMaxBit(
                 sort_begin, sort_end));
      }

    private:
      //! \brief The actual sorting procedure.
      static void doSort(const TIterator1& sort_begin,
                         const TIterator1& sort_end,
                         const TIterator2& copy1_begin,
                         const TIterator2& copy1_end,
                         const TIterator3& copy2_begin,
                         const TIterator3& copy2_end,
                         TValueType1 bit)
      {
        TIterator1 sort_left = sort_begin;
        TIterator1 sort_right = sort_end;
        TIterator2 copy1_left = copy1_begin;
        TIterator2 copy1_right = copy1_end;
        TIterator3 copy2_left = copy2_begin;
        TIterator3 copy2_right = copy2_end;

        // the right iterator point past-the-end, so decrease them
        --sort_right;
        --copy1_right;
        --copy2_right;

        while (sort_left != sort_right)
        {
          // advance the left iterator until it points to an element with a set
          // bit
          while ((sort_left != sort_right) && ((*sort_left & bit) == 0))
          {
            ++sort_left;
            ++copy1_left;
            ++copy2_left;
          }
          // decrease the right iterator until it points to an element with
          // a zero bit
          while ((sort_left != sort_right) && (*sort_right & bit))
          {
            --sort_right;
            --copy1_right;
            --copy2_right;
          }
          // swap the elements
          TValueType1 sort_temp = *sort_left;
          *sort_left = *sort_right;
          *sort_right = sort_temp;
          TValueType2 copy1_temp = *copy1_left;
          *copy1_left = *copy1_right;
          *copy1_right = copy1_temp;
          TValueType3 copy2_temp = *copy2_left;
          *copy2_left = *copy2_right;
          *copy2_right = copy2_temp;
        }

        // if all the elements had a zero-bit, we have to increase the right
        // iterator such that it points past-the-end again because otherwise
        // we would miss the last element in the next recursion
        if ((*sort_right & bit) == 0)
        {
          ++sort_right;
          ++copy1_right;
          ++copy2_right;
        }

        // advance to the next bit
        bit >>= 1;
        if (bit)
        {
          if (sort_right - sort_begin > 1)
            doSort(sort_begin, sort_right, copy1_begin, copy1_right,
                   copy2_begin, copy2_right, bit);
          if (sort_end - sort_right > 1)
            doSort(sort_right, sort_end, copy1_right, copy1_end,
                   copy2_right, copy2_end, bit);
        }
      }
  };

} // namespace agile

#endif // AGILE_NETWORK_RADIX_EXCHANGE_SORT_HPP

// End of $Id: radix_exchange_sort.hpp 476 2011-06-16 08:54:14Z freiberger $.
