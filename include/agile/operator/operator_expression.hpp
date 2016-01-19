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

// $Id: operator_expression.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_OPERATOR_EXPRESSION_HPP
#define AGILE_OPERATOR_OPERATOR_EXPRESSION_HPP

#include "agile/gpu_config.hpp"

namespace agile
{
  //! \brief Base class for operator expressions.
  //!
  //! This is the base class for all operator expressions. It uses the
  //! curiously recurring template pattern (CRTP) to know its derived class.
  template <typename TExpression>
  class OperatorExpression
  {
    public:
      typedef TExpression expression_type;

      //! \brief Downcast to derived type.
      inline
      const expression_type& operator() () const
      {
        return *static_cast<const expression_type*> (this);
      }

      //! \brief Downcast to derived type.
      inline
      expression_type& operator() ()
      {
        return *static_cast<expression_type*> (this);
      }

    protected:
      //! \brief Protected constructor.
      OperatorExpression() {}

      //! \brief Protected destructor.
      ~OperatorExpression() {}

    private:
      //! \brief Hidden assignment operator.
      const OperatorExpression& operator= (const OperatorExpression&);
  };

  // ============================= free functions =============================

  //! \brief Creator function for adjoints.
  template <typename TExpression>
  typename TExpression::adjoint_type adjoint(const TExpression& expression)
  {
    return expression.getAdjoint();
  }

} // namespace agile

#endif // AGILE_OPERATOR_OPERATOR_EXPRESSION_HPP

// End of $Id: operator_expression.hpp 476 2011-06-16 08:54:14Z freiberger $.
