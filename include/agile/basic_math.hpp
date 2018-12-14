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

// $Id: basic_math.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_BASIC_MATH_HPP
#define AGILE_BASIC_MATH_HPP

#include <cmath>

namespace agile
{
  //! Complex signum function.
  //!
  //! sgn(x) = 1       if x = 0
  //!          x/|x|   otherwise
  template <typename TType>
  inline
  std::complex<TType> sgn(const std::complex<TType>& x)
  {
    TType abs_x = std::abs(x);
    if (abs_x != TType(0))
      return x / abs_x;
    else
      return TType(1);
  }

} // namespace agile

#endif // AGILE_BASIC_MATH_HPP

// End of $Id: basic_math.hpp 476 2011-06-16 08:54:14Z freiberger $.
