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

// $Id: gpu_complex.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_COMPLEX_HPP
#define AGILE_GPU_COMPLEX_HPP

#include "gpu_type_traits.hpp"
#include <limits>

#ifdef __CUDACC__
#define DEVICEHOST __inline__ __device__ __host__
#else
#define DEVICEHOST inline
#endif

namespace agile
{
  //! A complex class on the GPU.
  template <typename TType>
  class GPUComplex
  {
    public:
      //! Default constructor. The value of the resultant variable is undefined.
      DEVICEHOST GPUComplex() { }

      //! Constructor taking real and imaginary part.
      DEVICEHOST GPUComplex(const TType& real, const TType& imag = 0.0)
        : x(real), y(imag)
      {
      }

      //! Constructor taking a GPU tuple type (e.g. float2, double2) to
      //! initialize the complex value.
      DEVICEHOST GPUComplex(const typename make_gpu_tuple<TType>::type& a)
        : x(a.x), y(a.y)
      {
      }

      //! Get the conjugate the complex number.
      //!
      //! \return The conjugated number (x, -y).
      DEVICEHOST GPUComplex conj() const
      {
        return GPUComplex(x, -y);
      }

      //! Get the imaginary part.
      DEVICEHOST const TType& imag() const
      {
        return y;
      }

      //! Get the square of the absolute value of a complex number.
      //!
      //! \param[in] a a complex number (x, y).
      //! \return x*x + y*y.
      DEVICEHOST TType norm() const
      {
        return x * x + y * y;
      }

      //! Get the real part.
      DEVICEHOST const TType& real() const
      {
        return x;
      }

      //! Get the sqrt part.
      DEVICEHOST GPUComplex complex_sqrt() const
      {
        TType vz=1;
        if(y<0)
          vz=-1;

        TType real = TType(sqrtf( (sqrtf(x*x+y*y)+x)/2 ));
        TType imag = TType(vz* sqrtf( (sqrtf(x*x+y*y)-x)/2 ));

        if (real != real) //if very low -> nan
//          real = 0;
          real = std::numeric_limits<TType>::min();
        if (imag != imag) //if very low -> nan
//          imag = 0;
          imag = std::numeric_limits<TType>::min();


        return GPUComplex(real,imag);

      }


      //! Assignment from a tuple type.
      DEVICEHOST GPUComplex& operator= (const float2& a)
      {
        x = a.x;
        y = a.y;
        return *this;
      }

      //! Assignment from a real number. The imaginary part is set to 0.
      DEVICEHOST GPUComplex& operator= (const TType& real)
      {
        x = real;
        y = TType(0);
        return *this;
      }

      //! Add a real value to the current complex number.
      DEVICEHOST GPUComplex<TType>& operator+= (const TType& a)
      {
        x += a;
        return *this;
      }

      //! Subtract a real value from the current complex number.
      DEVICEHOST GPUComplex<TType>& operator-= (const TType& a)
      {
        x -= a;
        return *this;
      }

      //! Multiply the current complex number with a real value.
      DEVICEHOST GPUComplex<TType>& operator*= (const TType& a)
      {
        x *= a;
        y *= a;
        return *this;
      }

      //! Divide the current complex number by a real value.
      DEVICEHOST GPUComplex<TType>& operator/= (const TType& a)
      {
        x /= a;
        y /= a;
        return *this;
      }

      //! Assignment from a complex number.
      DEVICEHOST GPUComplex& operator= (const GPUComplex& a)
      {
        x = a.x;
        y = a.y;
        return *this;
      }

      //! Add a complex value to the current one.
      DEVICEHOST GPUComplex<TType>& operator+= (const GPUComplex<TType>& a)
      {
        x += a.x;
        y += a.y;
        return *this;
      }

      //! Subtract a complex value from the current one.
      DEVICEHOST GPUComplex<TType>& operator-= (const GPUComplex<TType>& a)
      {
        x -= a.x;
        y -= a.y;
        return *this;
      }

      //! Multiply the current complex number with a complex value.
      DEVICEHOST GPUComplex<TType>& operator*= (const GPUComplex<TType>& a)
      {
        TType new_x = x * a.x - y * a.y;
        y =           y * a.x + x * a.y;
        x = new_x;
        return *this;
      }

      // =====================================================================
      // friend functions
      // =====================================================================
      //! Get the absolute value (i.e. the magnitude) of the complex number.
      //!
      //! \param[in] a a complex number (x, y).
      //! \return The magnitude sqrt(x*x + y*y).
   /*   friend DEVICEHOST TType abs(GPUComplex<TType>& a)
      {
        return sqrt(a.x * a.x + a.y * a.y);
      }
*/
//   TODO: friend double arg(complex a);

// TODO:   friend complex polar(double r, double t);

      //! Add two complex values.
      friend DEVICEHOST GPUComplex<TType> operator+ (const GPUComplex<TType>& a,
                                                     const GPUComplex<TType>& b)
      {
        return GPUComplex<TType>(a.x + b.x, a.y + b.y);
      }

      //! Add a real and a complex value.
      friend DEVICEHOST GPUComplex<TType> operator+ (const TType& a,
                                                     const GPUComplex<TType>& b)
      {
        return GPUComplex<TType>(a + b.x, b.y);
      }

      //! Add a complex and a real value.
      friend DEVICEHOST GPUComplex<TType> operator+ (const GPUComplex<TType>& a,
                                                     const TType& b)
      {
        return GPUComplex<TType>(a.x + b, a.y);
      }

      //! Subtract two complex values.
      friend DEVICEHOST GPUComplex<TType> operator- (const GPUComplex<TType>& a,
                                                     const GPUComplex<TType>& b)
      {
        return GPUComplex<TType>(a.x - b.x, a.y - b.y);
      }

      //! Subtract a complex from a real value.
      friend DEVICEHOST GPUComplex<TType> operator- (const TType& a,
                                                     const GPUComplex<TType>& b)
      {
        return GPUComplex<TType>(a - b.x, -b.y);
      }

      //! Subtract a real from a complex value.
      friend DEVICEHOST GPUComplex<TType> operator- (const GPUComplex<TType>& a,
                                                     const TType& b)
      {
        return GPUComplex<TType>(a.x - b, a.y);
      }

      //! Multiply two complex values.
      friend DEVICEHOST GPUComplex<TType> operator* (const GPUComplex<TType>& a,
                                                     const GPUComplex<TType>& b)
      {
        return GPUComplex<TType>(a.x * b.x - a.y * b.y,
                                 a.y * b.x + a.x * b.y);
      }

      //! Multiply a complex value with a real value.
      friend DEVICEHOST GPUComplex<TType> operator* (const GPUComplex<TType>& a,
                                                     const TType& b)
      {
        return GPUComplex<TType>(a.x * b, a.y * b);
      }

      //! Multiply a real value with a complex one.
      friend DEVICEHOST GPUComplex<TType> operator* (const TType& a,
                                                     const GPUComplex<TType>& b)
      {
        return GPUComplex<TType>(a * b.x, a * b.y);
      }

      //! Divide two complex values.
      friend DEVICEHOST GPUComplex<TType> operator/ (const GPUComplex<TType>& a,
                                                     const GPUComplex<TType>& b)
      {
        TType denom = b.x * b.x + b.y * b.y;
        return GPUComplex<TType>((a.x * b.x + a.y * b.y) / denom,
                                 (a.y * b.x - a.x * b.y) / denom);
      }

      //! Divide a complex by a complex value.
      friend DEVICEHOST GPUComplex<TType> operator/ (const GPUComplex<TType>& a,
                                                     const TType& b)
      {
        return GPUComplex<TType>(a.x / b, a.y / b);
      }

      //! Divide a real by a complex value.
      friend DEVICEHOST GPUComplex<TType> operator/ (const TType& a,
                                                     const GPUComplex<TType>& b)
      {
        TType factor = a / (b.x * b.x + b.y * b.y);
        return GPUComplex<TType>(b.x * factor, -b.y * factor);
      }

    private:
      //! The real part.
      TType x;

      //! The imaginary part.
      TType y;
  };
  
  template <typename TType>
  struct exp_Helper;

  template <>
  struct exp_Helper<float>
  {
      static DEVICEHOST float exp(float a)
      {
        return expf(a);
      }
  };
  template <>
  struct exp_Helper<double>
  {
      static DEVICEHOST double exp(double a)
      {
        return exp(a);
      }
  };

  template <typename TType>
  struct sincos_Helper;

  template <>
  struct sincos_Helper<float>
  {
      static DEVICEHOST void sincos(float a, float* s, float* c)
      {
        sincosf(a, s, c);
      }
  };
  template <>
  struct sincos_Helper<double>
  {
      static DEVICEHOST void sincos(double a, double* s, double* c)
      {
        sincos(a, s, c);
      }
  };

  namespace detail
  {
    //! The next struct is needed for dispatching between real and complex
    //! types. This is the definition for a complex type.
    template <typename TType, bool complex = is_complex<TType>::value>
    struct real_complex_dispatcher
    {
      typedef typename to_real_type<TType>::type real_type;

      //! Get the conjugated number.
      //!
      //! \param[in] a a complex number.
      //! \return \p a conjugated.
      static DEVICEHOST TType conj(const TType& a)
      {
        return a.conj();
      }

      //! Get the imaginary part of a complex number.
      //!
      //! \param[in] a a complex number (x, y).
      //! \return y
      static DEVICEHOST real_type imag(const TType& a)
      {
        return a.imag();
      }

      //! Get the norm of a complex number.
      //!
      //! \param[in] a a complex number (x, y).
      //! \return x^2 + y^2
      static DEVICEHOST real_type norm(const TType& a)
      {
        return a.norm();
      }
      
      //! Get the exp of a complex number.
      //!
      //! \param[in] a a complex number (x, y).
      //! \return Return exp(z)  =exp(x)*(cos(y)+i*sin(y)).
      static DEVICEHOST TType myExp(const TType& a)
      {
        real_type s, c;
        real_type e = exp_Helper<real_type>::exp(a.real());
        sincos_Helper<real_type>::sincos(a.imag(),&s,&c);
        return TType(e * c, e * s);
      }

      //! Get the real part of a complex number.
      //!
      //! \param[in] a a complex number (x, y).
      //! \return x
      static DEVICEHOST real_type real(const TType& a)
      {
        return a.real();
      }

      //! Get the sqrt of a complex number.
      //!
      //! \param[in] a a complex number (x, y).
      //! \return sqrt of a
      static DEVICEHOST TType sqroot(const TType& a)
      {
        return a.complex_sqrt();
      }
    };

    //! The partial specialisation for a real type.
    template <typename TType>
    struct real_complex_dispatcher<TType, false>
    {
      //! Get the conjugated real number.
      //!
      //! \param[in] a a real number.
      //! \return \p a itself.
      static DEVICEHOST TType conj(const TType& a)
      {
        return a;
      }

      //! Get the imaginary part of a real number.
      //!
      //! \param[in] a a real number.
      //! \return 0.
      static DEVICEHOST TType imag(const TType& /*a*/)
      {
        return 0;
      }

      //! Get the norm of a real number.
      //!
      //! \param[in] a a real number.
      //! \return Analoguously to the norm of a complex number a^2 is returned.
      static DEVICEHOST TType norm(const TType& a)
      {
        return a*a;
      }
      
      //! Get the exp of a real number.
      //!
      //! \param[in] a a real number.
      //! \return Return exp(a).
      static DEVICEHOST TType myExp(const TType& a)
      {
        return TType(exp(a));
      }

      //! Get the real part of a real number.
      //!
      //! \param[in] a a real number.
      //! \return The real number \p a.
      static DEVICEHOST TType real(const TType& a)
      {
        return a;
      }

      //! Get the sqrt of a real number.
      //!
      //! \param[in] a a real number.
      //! \return The sqrt number of \p a.
      static DEVICEHOST TType sqroot(const TType& a)
      {
        return std::sqrt(a);
      }
    };

  } // namespace detail

  //! \brief Conjugate a number.
  //!
  //! \param[in] x a real or complex number.
  //! \return The conjugated number.
  template <typename TType>
  DEVICEHOST TType conj(const TType& x)
  {
    return detail::real_complex_dispatcher<TType>::conj(x);
  }

  //! \brief Get the imag part of a number.
  //!
  //! \param[in] x a real or complex number.
  //! \return The real part of the number.
  template <typename TType>
  DEVICEHOST typename to_real_type<TType>::type imag(const TType& x)
  {
    return detail::real_complex_dispatcher<TType>::imag(x);
  }

  //! \brief Get the norm of a number.
  //!
  //! Get the norm of a real or complex number. For a real number this
  //! method will return the squared value, for a complex number it will
  //! return the sum of squares of the real and the imaginary part (and the
  //! square of the absolute value, thus).
  template <typename TType>
  DEVICEHOST typename to_real_type<TType>::type norm(const TType& x)
  {
    return detail::real_complex_dispatcher<TType>::norm(x);
  }
  
  //! \brief Get the exp of a number.
  //!
  //! Get the exp of a real or complex number. 
  template <typename TType>
  DEVICEHOST TType exp(const TType& x)
  {
    return detail::real_complex_dispatcher<TType>::myExp(x);
  }

  //! \brief Get the real part of a number.
  //!
  //! \param[in] x a real or complex number.
  //! \return The imaginary part of the number. If the number is a real type,
  //! zero is returned.
  template <typename TType>
  DEVICEHOST typename to_real_type<TType>::type real(const TType& x)
  {
    return detail::real_complex_dispatcher<TType>::real(x);
  }

  //! \brief Get the sqrt number.
  //!
  //! \param[in] x a real or complex number.
  //! \return The sqrt of the number.
  template <typename TType>
  DEVICEHOST TType cusqrt(const TType& x)
  {
    return detail::real_complex_dispatcher<TType>::sqroot(x);
  }


} // namespace agile

#endif // AGILE_GPU_COMPLEX_HPP

// End of $Id: gpu_complex.hpp 476 2011-06-16 08:54:14Z freiberger $.
