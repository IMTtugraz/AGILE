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

// $Id: exception.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_EXCEPTION_HPP
#define AGILE_EXCEPTION_HPP

// Exception modelled after the corresponding class in Deal.II.

#include "agile/gpu_config.hpp"

#include <iostream>
#include <cstdlib> // std::abort

namespace agile
{
  //! \brief A class defining an AGILE exception.
  class AGILEException : public std::exception
  {
    public:
      //! \brief Default constructor.
      AGILEException() throw() {}

      //! \brief Constructor.
      //!
      //! \param[in] file_name The name of the file in which the exception
      //! was thrown.
      //! \param[in] line The line number in which the exception was thrown.
      //! \param[in] function_name The name of the function in which the
      //! exception was thrown.
      //! \param[in] exception_call The code for calling the exception. This
      //! is typically a constructor of a class derived from
      //! \p AGILEException.
      AGILEException(
        const char* file_name, unsigned int line,
        const char* function_name, const char* exception_call) throw()
        : m_file_name(file_name), m_line(line), m_function_name(function_name),
          m_exception_call(exception_call)
      {
      }

      //! \brief Copy constructor.
      AGILEException(const AGILEException& exception) throw()
        : m_file_name(exception.m_file_name), m_line(exception.m_line),
          m_function_name(exception.m_function_name),
          m_exception_call(exception.m_exception_call)
      {
      }

      //! \brief Destructor.
      virtual
      ~AGILEException() throw() {}

      //! \brief Method to set the parameters.
      //!
      //! \param[in] file_name The name of the file in which the exception
      //! was thrown.
      //! \param[in] line The line number in which the exception was thrown.
      //! \param[in] function_name The name of the function in which the
      //! exception was thrown.
      //! \param[in] exception_call The code for calling the exception. This
      //! is typically a constructor of a class derived from
      //! \p AGILEException.
      void set(const char* file_name, unsigned int line,
               const char* function_name, const char* exception_call)
      {
        m_file_name = file_name;
        m_line = line;
        m_function_name = function_name;
        m_exception_call = exception_call;
      }

      //! \brief Method to print information on where the exception was raised.
      //!
      //! \param[out] out_stream A stream to output the information.
      void printOccurence(std::ostream& out_stream) const
      {
        out_stream << "In file '" << m_file_name << "'" << std::endl
                   << "in line " << m_line << std::endl
                   << "in function '" << m_function_name << "'" << std::endl
                   << "a '" << m_exception_call << "'" << std::endl
                   << "was raised." << std::endl;
      }

      //! \brief Method to print an additional message.
      //!
      //! Derived classes can override this method to output some additional
      //! message which should help the user to understand why the exception
      //! was thrown.
      virtual
      void printMessage(std::ostream& out_stream) const {}

      //! \brief Get the error message.
      //!
      //! \return The code for calling the exception.
      virtual
      const char* what() const throw()
      {
        return m_exception_call;
      }

    private:
      //! \brief The file in which this exception was raised.
      const char* m_file_name;

      //! \brief The line in which this exception was raised.
      unsigned int m_line;

      //! \brief The function in which this exception was thrown.
      const char* m_function_name;

      //! \brief The exception thrown as string.
      const char* m_exception_call;
  };

  //! \brief Print an exception and exit.
  //!
  //! This function is called by the \p AGILE_ASSERT macro.
  //! It sets the members of \p imt_cuda_exception to the values
  //! given, prints the exception's messages and aborts the program.
  template <typename TException>
  void handleException(const char* file_name, unsigned int line,
                       const char* function_name, const char* exception_call,
                       TException imt_cuda_exception)
  {
    imt_cuda_exception.set(file_name, line, function_name, exception_call);

    std::cerr << "----------------------------------------"
                 "----------------------------------------" << std::endl;
    imt_cuda_exception.printOccurence(std::cerr);
    std::cerr << "                    --------------------"
                 "--------------------" << std::endl;
    imt_cuda_exception.printMessage(std::cerr);
    std::cerr << "----------------------------------------"
                 "----------------------------------------" << std::endl;

    // halt the program
    std::abort();
  }
/*
  //! \brief A helper function to wrap the exception.
  //!
  //! The function \p handleException has to take a non-const reference to the
  //! base class in order to set the
  template <typename TException>
  void handleExceptionHelper(const char* file_name, unsigned int line,
                       const char* function_name, const char* exception_call,
                       AGILEException& imt_cuda_exception)
*/
}  // namespace agile

//! \brief Declare an exception that does not take any parameter.
#define DeclareAGILEException0(AGILEException0) \
class AGILEException0 :  public AGILEException {}

//! \brief Declare an exception that takes a message as parameter.
#define DeclareAGILEException1(AGILEException1, message) \
class AGILEException1 : public AGILEException            \
{                                                            \
  public:                                                    \
    virtual                                                  \
    ~AGILEException1() throw() {}                          \
                                                             \
    virtual                                                  \
    void printMessage(std::ostream &out_stream) const        \
    {                                                        \
      out_stream << message << std::endl;                    \
    }                                                        \
};

//! \brief Declare an exception that takes a message and another parameter.
#define DeclareAGILEException2(AGILEException2, message, type1) \
class AGILEException2 : public AGILEException                   \
{                                                                   \
  private:                                                          \
    type1 arg1;                                                     \
                                                                    \
  public:                                                           \
    AGILEException2(type1 argument1)                              \
      : arg1(argument1)                                             \
    {}                                                              \
                                                                    \
    virtual                                                         \
    ~AGILEException2() throw() {}                                 \
                                                                    \
                                                                    \
    virtual                                                         \
    void printMessage(std::ostream &out_stream) const               \
    {                                                               \
      out_stream << message << std::endl;                           \
    }                                                               \
};

//! \brief Declare an exception that takes a message and 2 other parameters.
#define DeclareAGILEException3(AGILEException3, message, type1, type2) \
class AGILEException3 : public AGILEException                          \
{                                                                          \
  private:                                                                 \
    type1 arg1;                                                            \
    type2 arg2;                                                            \
                                                                           \
  public:                                                                  \
    AGILEException3(type1 argument1, type2 argument2)                    \
      : arg1(argument1), arg2(argument2)                                   \
    {}                                                                     \
                                                                           \
    virtual                                                                \
    ~AGILEException3() throw() {}                                        \
                                                                           \
                                                                           \
    virtual                                                                \
    void printMessage(std::ostream &out_stream) const                      \
    {                                                                      \
      out_stream << message << std::endl;                                  \
    }                                                                      \
};

//! \brief Declare an exception that takes a message and 3 other parameters.
#define DeclareAGILEException4(AGILEException4, message, type1, type2, \
                                 type3)                                    \
class AGILEException4 : public AGILEException                          \
{                                                                          \
  private:                                                                 \
    type1 arg1;                                                            \
    type2 arg2;                                                            \
    type3 arg3;                                                            \
                                                                           \
  public:                                                                  \
    AGILEException4(type1 argument1, type2 argument2, type3 argument3)   \
      : arg1(argument1), arg2(argument2), arg3(argument3)                  \
    {}                                                                     \
                                                                           \
                                                                           \
    virtual                                                                \
    ~AGILEException4() throw() {}                                        \
                                                                           \
    virtual                                                                \
    void printMessage(std::ostream &out_stream) const                      \
    {                                                                      \
      out_stream << message << std::endl;                                  \
    }                                                                      \
};

//! \brief Declare an exception that takes a message and 4 other parameters.
#define DeclareAGILEException5(AGILEException5, message, type1, type2, \
                                 type3, type4)                             \
class AGILEException5 : public AGILEException                          \
{                                                                          \
  private:                                                                 \
    type1 arg1;                                                            \
    type2 arg2;                                                            \
    type3 arg3;                                                            \
    type4 arg4;                                                            \
                                                                           \
  public:                                                                  \
    AGILEException5(type1 argument1, type2 argument2,                    \
                      type3 argument3, type4 argument4)                    \
      : arg1(argument1), arg2(argument2), arg3(argument3),                 \
        arg4(argument4)                                                    \
    {}                                                                     \
                                                                           \
                                                                           \
    virtual                                                                \
    ~AGILEException5() throw() {}                                        \
                                                                           \
    virtual                                                                \
    void printMessage(std::ostream &out_stream) const                      \
    {                                                                      \
      out_stream << message << std::endl;                                  \
    }                                                                      \
};

namespace agile
{
  namespace StandardException
  {
    //! \brief A general internal exception.
    DeclareAGILEException1(ExceptionInternal,
                             "An internal exception occured.");

    //! \brief A simple exception with a variable message.
    //!
    //! Call this exception with <tt>ExceptionMessage("Some text")</tt>.
    DeclareAGILEException2(ExceptionMessage, arg1, const char*);

    //! \brief Exception for index out of bounds.
    //!
    //! This exception can be thrown if a variable/a value is outside the
    //! range <tt>[lower, upper)</tt>. Use it like
    //! <tt>ExceptionIndexOutOfBounds(variable_name, lower, upper)</tt>.
    DeclareAGILEException4(ExceptionIndexOutOfBounds,
                             arg1 << " is not in [" << arg2 << ", " << arg3
                             << ")", const char*, int, int);

    //! \brief Exception for invalid program flow.
    //!
    //! Throw this exception whenever the program comes to a point that
    //! should not be reached under normal conditions. An example might
    //! be a code like the following:
    //! \code
    //! switch(expression)
    //! {
    //!   case 1: return ...;
    //!   case 2: return ...;
    //!   case 3: return ...;
    //! }
    //! AGILE_ASSERT(false, ExceptionInvalidProgramFlow());
    //! \endcode
    //! In the example, the switch-case-statement should handle all possible
    //! conditions and the program should never reach the point past the
    //! switch statement.
    DeclareAGILEException1(ExceptionInvalidProgramFlow,
                             "The program flow reached an invalid point.");

    //! \brief Exception if a value is below a lower bound.
    //!
    //! Throw this exception whenever a value is too small, i.e. if
    //! value < lower_bound. Use this exception like
    //! <tt>ExceptionValueBelowLowerBound(variable_name, lower_bound)</tt>.
    DeclareAGILEException3(ExceptionValueTooSmall,
                             arg1 << " is smaller than " << arg2,
                             const char*, int);

    //! \brief Exception if a value exceeded an upper bound.
    //!
    //! This exception can be used to indicate that a value is larger than
    //! or equal to an upper limit.
    //! Use this exception like
    //! <tt>ExceptionValueAboveUpperBound(variable_name, upper_bound)</tt>.
    DeclareAGILEException3(ExceptionValueTooLarge,
                             arg1 << " is equal to or larger than " << arg2,
                             char*, int);

    //! \brief Exception for a size mismatch.
    //!
    //! Throw this exception if the sizes of two objects do not match. The
    //! call should be
    //! <tt>ExceptionSizeMismatch(variable1, variable2, size1, size2)</tt>.
    DeclareAGILEException5(ExceptionSizeMismatch,
                             arg1 << " has size " << arg3 << " but "
                             << arg2 << " has size " << arg4,
                             const char*, const char*, int, int);

    //! \brief Exception which is thrown in case an in/output operation fails.
    DeclareAGILEException1(ExceptionIO, "IO operation failed.");

  }  // namespace StandardException

}  // namespace agile


#ifdef AGILE_DEBUG

// in debug mode, the assertion has to be enabled
// here the condition is tests for truth rather than testing its negation;
// the reason is that this does not disturb C++ shortcuts during evaluation
#define AGILE_ASSERT(condition, exception)                         \
{                                                                   \
  if (condition)                                                    \
  {                                                                 \
  }                                                                 \
  else                                                              \
  {                                                                 \
    handleException (__FILE__, __LINE__, __PRETTY_FUNCTION__,  \
                          #exception, (exception));                 \
  }                                                                 \
}

#else

// if we are not in debug mode, disable the assertions
#define AGILE_ASSERT(condition, exception) /* do nothing */
#endif


#endif // AGILE_EXCEPTION_HPP

// End of $Id: exception.hpp 476 2011-06-16 08:54:14Z freiberger $.
