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

// $Id: gpu_type_traits.hpp 476 2011-06-16 08:54:14Z freiberger $

#ifndef AGILE_GPU_TYPE_TRAITS_HPP
#define AGILE_GPU_TYPE_TRAITS_HPP

#include "agile/gpu_config.hpp"
#include <complex>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cuComplex.h>

#if __CUDA_ARCH__ >= 130
#include <sm_13_double_functions.h>
#endif

namespace agile
{
  // forward declaration
  template <typename TType>
  class GPUComplex;


  // ----------------------------------------------------------------------
  // to_real_type
  // ----------------------------------------------------------------------
  //! \brief Template that returns the base type of a complex type.
  template <typename TType>
  struct to_real_type
  {
    typedef TType type;
  };

  // Partial specialisation for complex types.
  template <typename TType>
  struct to_real_type<std::complex<TType> >
  {
    typedef TType type;
  };

  // Partial specialisation for GPUComplex types.
  template <typename TType>
  struct to_real_type<GPUComplex<TType> >
  {
    typedef TType type;
  };


  // ----------------------------------------------------------------------
  // to_tuple_type
  // ----------------------------------------------------------------------
  //! \brief Template that returns the CUDA tuple type for a complex type.
  //!
  //! If the input is a real type, this method will return the input type.
  //! If the input is either a <tt>std::complex<T></tt> or a
  //! <tt>GPUComplex<T></tt> the struct will return the appropriate
  //! CUDA tuple type.
  template <typename TType>
  struct to_tuple_type
  {
    typedef TType type;
    typedef TType texture_type;

#ifdef __CUDACC__
    __inline__ __device__ static
    type texture2type(const texture_type& value)
    {
      return value;
    }
#endif
  };

  // Specialisation for double returning int2 for textures.
  template <>
  struct to_tuple_type<double>
  {
    typedef double type;
    typedef int2 texture_type;

#if __CUDA_ARCH__ >= 130
#ifdef __CUDACC__
    __inline__ __device__ static
    type texture2type(const texture_type& value)
    {
      return __hiloint2double(value.y, value.x);
    }
#endif
#endif /* __CUDA_ARCH__ > 130 */
  };

  // Partial specialisation for complex types.
  template <typename TType>
  struct to_tuple_type<std::complex<TType> >
  {
    typedef typename to_tuple_type<GPUComplex<TType> >::type type;
    typedef typename to_tuple_type<GPUComplex<TType> >::texture_type texture_type;

#ifdef __CUDACC__
    __inline__ __device__ static
    type texture2type(const texture_type& value)
    {
      return to_tuple_type<GPUComplex<TType> >::texture2type(value);
    }
#endif
  };

  // Specialisation for a GPUComplex<float> returning float2.
  template <>
  struct to_tuple_type<GPUComplex<float> >
  {
    typedef float2 type;
    typedef float2 texture_type;

#ifdef __CUDACC__
    __inline__ __device__ static
    type texture2type(const texture_type& value)
    {
      return value;
    }
#endif
  };

  // Specialisation for a GPUComplex<double> returning double2.
  template <>
  struct to_tuple_type<GPUComplex<double> >
  {
    typedef double2 type;
    typedef int4 texture_type;

#if __CUDA_ARCH__ >= 130
#ifdef __CUDACC__
    __inline__ __device__ static
    type texture2type(const texture_type& value)
    {
      return make_double2(__hiloint2double(value.y, value.x),
                          __hiloint2double(value.w, value.z));
    }
#endif
#endif /* __CUDA_ARCH__ > 130 */
  };


  // ----------------------------------------------------------------------
  // make_gpu_tuple
  // ----------------------------------------------------------------------
  //! \brief Template that returns the CUDA tuple type for a given input type.
  template <typename TType>
  struct make_gpu_tuple { };

  // Specialisation for char.
  template <>
  struct make_gpu_tuple<char> { typedef char2 type; };
  // Specialisation for unsigned char.
  template <>
  struct make_gpu_tuple<unsigned char> { typedef uchar2 type; };
  // Specialisation for int.
  template <>
  struct make_gpu_tuple<int> { typedef int2 type; };
  // Specialisation for unsigned.
  template <>
  struct make_gpu_tuple<unsigned> { typedef uint2 type; };
  // Specialisation for float.
  template <>
  struct make_gpu_tuple<float> { typedef float2 type; };
  // Specialisation for double.
  template <>
  struct make_gpu_tuple<double> { typedef double2 type; };


  // ----------------------------------------------------------------------
  // is_complex
  // ----------------------------------------------------------------------
  //! \brief Template that returns true, if \p TType is a complex (GPU) type.
  template <typename TType>
  struct is_complex { enum { value = false }; };

  // Specialisation for std::complex<T>.
  template <typename TType>
  struct is_complex<std::complex<TType> > { enum { value = true }; };

  // Specialisation for GPUComplex<T>.
  template <typename TType>
  struct is_complex<GPUComplex<TType> > { enum { value = true }; };


  // Specialisation for CublasComplexFloat
  template <typename TType>
  struct CUBLASComplexHelper
  {
      typedef cuComplex type;
  };

  // Specialisation for CublasComplexDouble
  template <>
  struct CUBLASComplexHelper<double>
  {
      typedef cuDoubleComplex type;
  };


  // ----------------------------------------------------------------------
  // substitute_gpu_complex
  // ----------------------------------------------------------------------
  //! \brief Substitute complex types by its GPU equivalents.
  //!
  //! If the input type is a \p std::complex<TType>, the output type will be
  //! a \p GPUComplex<TType>. For all other types the input type is returned.
  template <typename TType>
  struct substitute_gpu_complex
  {
    typedef TType type;
    typedef TType cublas_type;
  };

  // Partial specialisation for a complex type.
  template <typename TType>
  struct substitute_gpu_complex<std::complex<TType> >
  {
    typedef GPUComplex<TType> type;
    typedef typename CUBLASComplexHelper<TType>::type cublas_type;
  };


  // ----------------------------------------------------------------------
  // promote
  // ----------------------------------------------------------------------
  //! \brief Template that returns a complex type if one of the input types is
  //! complex.
  //!
  //! The two base types have to be equal which means that both \p TType1 and
  //! \p TType2 can only be TType or \p std::complex<TType>. In other words:
  //! Do not instantiate this template with mixed types as \p unsigned and
  //! \p double, for example. However, this promotion should be added in the
  //! future.
  template <typename TType1, typename TType2>
  struct promote
  {
  };

  // Partial specialisation for equal types.
  template <typename TType>
  struct promote<TType, TType>
  {
    typedef TType type;
  };

  // Partial specialisation if the first type is the complex version of the 2nd
  // type.
  template <typename TType>
  struct promote<std::complex<TType>, TType>
  {
    typedef std::complex<TType> type;
  };

  // Partial specialisation if the second type is the complex version of the 1st
  // type.
  template <typename TType>
  struct promote<TType, std::complex<TType> >
  {
    typedef std::complex<TType> type;
  };

  // Partial specialisation if the first type is a GPU complex type.
  template <typename TType>
  struct promote<GPUComplex<TType>, TType>
  {
    typedef GPUComplex<TType> type;
  };

  // Partial specialisation if the second type is a GPU complex type.
  template <typename TType>
  struct promote<TType, GPUComplex<TType> >
  {
    typedef GPUComplex<TType> type;
  };

  // ----------------------------------------------------------------------
  // to_float_type
  // ----------------------------------------------------------------------
  //! \brief Template that returns a floating point type.
  //!
  //! This method returns always a floating point type. If \p TType is
  //! already floating point, \p TType is returned. Otherwise, the default
  //! floating point type (which is \p float) is returned.
  template <typename TType>
  struct to_float_type { typedef float type; };

  // Specialisation for double.
  template <>
  struct to_float_type<double> { typedef double type; };

} // namespace agile

#endif // AGILE_GPU_TYPE_TRAITS_HPP

// End of $Id: gpu_type_traits.hpp 476 2011-06-16 08:54:14Z freiberger $.
