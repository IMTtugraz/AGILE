//$Id: binary_operator_expression.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_OPERATOR_BINARY_OPERATOR_EXPRESSION_HPP
#define AGILE_OPERATOR_BINARY_OPERATOR_EXPRESSION_HPP

#include "agile/gpu_config.hpp"
#include "agile/operator/forward_operator.hpp"
#include "agile/operator/inverse_operator.hpp"

namespace agile
{
  //! \brief Base class for a binary operator applied on forward operator
  //! expressions.
  //!
  //! With this class it is possible to define operators acting on operator
  //! expressions. The class takes the types of the first and the second
  //! expression and the operator class implementing the mathematical
  //! operations as template arguments.
  template <typename TExpr1, typename TExpr2, class TOperator>
  class BinaryForwardOperatorExpression
    : public ForwardOperatorExpression<
        BinaryForwardOperatorExpression<TExpr1, TExpr2, TOperator> >
  {
      //! The first expression.
      TExpr1 m_expression1;
      //! The second expression.
      TExpr2 m_expression2;

    public:
      // this base class cannot know the type of the adjoint as it depends
      // on TOperator, so take it from there
      typedef typename TOperator::forward_adjoint_type adjoint_type;

      //! Constructor for a binary operator expression acting on operator
      //! expressions.
      BinaryForwardOperatorExpression(const TExpr1& expr1, const TExpr2& expr2)
        : m_expression1(expr1), m_expression2(expr2)
      {
      }

      //! Applies the two expressions and the binary operator to the vector
      //! \p x. The result will be stored in \p y. It is the user's
      //! responsibility to resize the vectors as needed.
      template <typename TVector>
      void operator() (const TVector& x, TVector& y)
      {
        TOperator::apply(m_expression1, m_expression2, x, y);
      }

      adjoint_type getAdjoint() const
      {
        return TOperator::getForwardAdjoint(m_expression1, m_expression2);
      }
  };

  //! \brief Base class for a binary operator applied on inverse operator
  //! expressions.
  //!
  //! \see BinaryForwardOperatorExpression
  template <typename TExpr1, typename TExpr2, class TOperator>
  class BinaryInverseOperatorExpression
    : public InverseOperatorExpression<
        BinaryInverseOperatorExpression<TExpr1, TExpr2, TOperator> >
  {
      //! The first expression.
      TExpr1 m_expression1;
      //! The second expression.
      TExpr2 m_expression2;

    public:
      // this base class cannot know the type of the adjoint as it depends
      // on TOperator, so take it from there
      typedef typename TOperator::inverse_adjoint_type adjoint_type;

      //! Constructor for a binary operator expression acting on operator
      //! expressions.
      BinaryInverseOperatorExpression(const TExpr1& expr1, const TExpr2& expr2)
        : m_expression1(expr1), m_expression2(expr2)
      {
      }

      //! Applies the two expressions and the binary operator to the vector
      //! \p x. The result will be stored in \p y. It is the user's
      //! responsibility to resize the vectors as needed.
      template <typename TVector>
      void operator() (const TVector& x, TVector& y)
      {
        TOperator::apply(m_expression1, m_expression2, x, y);
      }

      adjoint_type getAdjoint() const
      {
        return TOperator::getInverseAdjoint(m_expression1, m_expression2);
      }
  };

  //! \brief Binary operator to add operator expressions.
  //!
  //! This is the implementation of the operator addition. Let \f$ F \f$ and
  //! \f$ G \f$ denote two operator expressions. The addition is defined by
  //! \f$ (F + G)(x) = F(x) + G(x) \f$. The adjoint type is defined as
  //! \f$ (F + G)^*(x) = F^*(x) + G^*(x) \f$.
  template <typename TExpr1, typename TExpr2>
  struct OperatorExpressionAddition
  {
    typedef BinaryForwardOperatorExpression<typename TExpr1::adjoint_type,
                                            typename TExpr2::adjoint_type,
                                            OperatorExpressionAddition<
                                              typename TExpr1::adjoint_type,
                                              typename TExpr2::adjoint_type> >
      forward_adjoint_type;

    typedef BinaryInverseOperatorExpression<typename TExpr1::adjoint_type,
                                            typename TExpr2::adjoint_type,
                                            OperatorExpressionAddition<
                                              typename TExpr1::adjoint_type,
                                              typename TExpr2::adjoint_type> >
      inverse_adjoint_type;

    //! Apply the addition operator.
    //!
    //! (F + G)(x) = F(x) + G(x)
    template <typename TVector>
    static void apply(TExpr1& F, TExpr2& G, const TVector& x, TVector& y)
    {
      // (F + G)(x) = F(x) + G(x)
      TVector y_temp(y.size());
      F(x, y_temp);
      G(x, y);
      addVector(y, y_temp, y);
    }

    static forward_adjoint_type getForwardAdjoint(TExpr1& F, TExpr2& G)
    {
      // (F + G)^*(x) = F^*(x) + G^*(x)
      return forward_adjoint_type(F.getAdjoint(), G.getAdjoint());
    }

    static inverse_adjoint_type getInverseAdjoint(TExpr1& F, TExpr2& G)
    {
      // (F + G)^*(x) = F^*(x) + G^*(x)
      return inverse_adjoint_type(F.getAdjoint(), G.getAdjoint());
    }

  };

  //! \brief Binary operator to concatenate operator expressions.
  //!
  //! This class allows the user to chain operators. For \f$ F \f$ and
  //! \f$ G \f$ being operator expressions, the concatenation is implemented
  //! as \f$ (F o G)(x) = F(G(x)) \f$. For matrices this is equivalent to
  //! multiplying \f$ x \f$ with the product \f$ F \cdot G \f$.
  template <typename TExpr1, typename TExpr2>
  struct OperatorExpressionConcatenation
  {
    typedef BinaryForwardOperatorExpression<typename TExpr2::adjoint_type,
                                            typename TExpr1::adjoint_type,
                                              OperatorExpressionConcatenation<
                                                typename TExpr2::adjoint_type,
                                                typename TExpr1::adjoint_type> >
      forward_adjoint_type;

    typedef BinaryInverseOperatorExpression<typename TExpr2::adjoint_type,
                                            typename TExpr1::adjoint_type,
                                            OperatorExpressionConcatenation<
                                              typename TExpr2::adjoint_type,
                                              typename TExpr1::adjoint_type> >
      inverse_adjoint_type;

    //! \brief Apply the concatenation operator.
    //!
    //! (F o G)(x) = F(G(x))
    template <typename TVector>
    static void apply(TExpr1& F, TExpr2& G, const TVector& x, TVector& y)
    {
      // (F o G)(x) = F(G(x))
      TVector y_temp(y.size());
      G(x, y_temp);
      F(y_temp, y);
    }

    //! \brief Get the adjoint of two concatenated operators.
    //!
    //! (F o G)^*(x) = (G^* o F^*)(x)
    static forward_adjoint_type getForwardAdjoint(TExpr1& F, TExpr2& G)
    {
      // (F o G)^*(x) = (G^* o F^*)(x)
      return forward_adjoint_type(G.getAdjoint(), F.getAdjoint());
    }

    //! \brief Get the adjoint of two concatenated operators.
    //!
    //! (F o G)^*(x) = (G^* o F^*)(x)
    static inverse_adjoint_type getInverseAdjoint(TExpr1& F, TExpr2& G)
    {
      // (F o G)^*(x) = (G^* o F^*)(x)
      return inverse_adjoint_type(G.getAdjoint(), F.getAdjoint());
    }
  };

  //! \brief Binary operator to multiply an operator expression with a scalar.
  template <typename TScalar, typename TExpr>
  struct OperatorExpressionScalarMultiplication
  {
    typedef BinaryForwardOperatorExpression<
      TScalar, typename TExpr::adjoint_type,
      OperatorExpressionScalarMultiplication<
        TScalar, typename TExpr::adjoint_type> >
    forward_adjoint_type;

    typedef BinaryInverseOperatorExpression<
      TScalar, typename TExpr::adjoint_type,
      OperatorExpressionScalarMultiplication<
        TScalar, typename TExpr::adjoint_type> >
    inverse_adjoint_type;

    //! \brief Apply the scalar multiplication operator.
    //!
    //! (aF)(x) = aF(x)
    template <typename TVector>
    static void apply(TScalar& a, TExpr& F, const TVector& x, TVector& y)
    {
      // (aF)(x) = aF(x)
      F(x, y);
      scale(a, y, y);
    }

    //! \brief Get the adjoint of an operator multiplied with a scalar.
    //!
    //! (aF)^*(x) = a^*F^*(x)
    static forward_adjoint_type getForwardAdjoint(TScalar& a, TExpr& F)
    {
      // (aF)^*(x) = a^*F^*(x)
      return forward_adjoint_type(conj(a), F.getAdjoint());
    }

    //! \brief Get the adjoint of an operator multiplied with a scalar.
    //!
    //! (aF)^*(x) = a^*F^*(x)
    static inverse_adjoint_type getInverseAdjoint(TScalar& a, TExpr& F)
    {
      // (aF)^*(x) = a^*F^*(x)
      return inverse_adjoint_type(conj(a), F.getAdjoint());
    }
  };

  // ============================= free functions =============================

  //! Creator function for adding two forward operator expressions.
  template <typename TExpression1, typename TExpression2>
  BinaryForwardOperatorExpression<
    TExpression1, TExpression2,
    OperatorExpressionAddition<TExpression1,
                               TExpression2> >
  operator+ (const ForwardOperatorExpression<TExpression1>& expression1,
             const ForwardOperatorExpression<TExpression2>& expression2)
  {
    return BinaryForwardOperatorExpression<
             TExpression1, TExpression2,
             OperatorExpressionAddition<TExpression1, TExpression2> >(
               expression1(), expression2());
  }

  //! Creator function for adding two inverse operator expressions.
  template <typename TExpression1, typename TExpression2>
  BinaryForwardOperatorExpression<
    TExpression1, TExpression2,
    OperatorExpressionAddition<TExpression1,
                               TExpression2> >
  operator+ (const InverseOperatorExpression<TExpression1>& expression1,
             const InverseOperatorExpression<TExpression2>& expression2)
  {
    return BinaryInverseOperatorExpression<
             TExpression1, TExpression2,
             OperatorExpressionAddition<TExpression1, TExpression2> >(
               expression1(), expression2());
  }

  //! Creator function for forward operator concatenation.
  template <typename TExpression1, typename TExpression2>
  BinaryForwardOperatorExpression<
    TExpression1, TExpression2,
    OperatorExpressionConcatenation<TExpression1,
                                    TExpression2> >
  operator% (const ForwardOperatorExpression<TExpression1>& expression1,
             const ForwardOperatorExpression<TExpression2>& expression2)
  {
    return BinaryForwardOperatorExpression<
             TExpression1, TExpression2,
             OperatorExpressionConcatenation<TExpression1, TExpression2> >(
               expression1(), expression2());
  }

  //! Creator function for inverse operator concatenation.
  template <typename TExpression1, typename TExpression2>
  BinaryInverseOperatorExpression<
    TExpression1, TExpression2,
    OperatorExpressionConcatenation<TExpression1,
                                    TExpression2> >
  operator% (const InverseOperatorExpression<TExpression1>& expression1,
             const InverseOperatorExpression<TExpression2>& expression2)
  {
    return BinaryInverseOperatorExpression<
             TExpression1, TExpression2,
             OperatorExpressionConcatenation<TExpression1, TExpression2> >(
               expression1(), expression2());
  }

  //! Creator function for operator multiplication with a scalar.
  template <typename TScalar, typename TExpression>
  BinaryForwardOperatorExpression<
    TScalar, TExpression,
    OperatorExpressionScalarMultiplication<TScalar, TExpression> >
  operator* (const TScalar& scalar,
             const ForwardOperatorExpression<TExpression>& expression)
  {
    return BinaryForwardOperatorExpression<
             TScalar, TExpression,
             OperatorExpressionScalarMultiplication<TScalar, TExpression> >(
               scalar, expression());
  }

} // namespace agile

#endif // AGILE_OPERATOR_BINARY_OPERATOR_EXPRESSION_HPP

// End of $Id: binary_operator_expression.hpp 476 2011-06-16 08:54:14Z freiberger $.
