#ifndef __HLIB_BLAS_TRAITS_HH
#define __HLIB_BLAS_TRAITS_HH
//
// Project     : HLib
// File        : traits.hh
// Description : provide type traits for BLAS functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/base/types.hh"

namespace HLIB
{

//!
//! signals integer types
//!
template <typename T> struct is_integer                   { static constexpr bool value = false; };
template <>           struct is_integer< char >           { static constexpr bool value = true; };
template <>           struct is_integer< unsigned char >  { static constexpr bool value = true; };
template <>           struct is_integer< short >          { static constexpr bool value = true; };
template <>           struct is_integer< unsigned short > { static constexpr bool value = true; };
template <>           struct is_integer< int >            { static constexpr bool value = true; };
template <>           struct is_integer< unsigned int >   { static constexpr bool value = true; };
template <>           struct is_integer< long >           { static constexpr bool value = true; };
template <>           struct is_integer< unsigned long >  { static constexpr bool value = true; };

#if __cplusplus >= 201703L
template <typename T> inline constexpr bool is_integer_v = is_integer< T >::value;
#endif

//!
//! signals floating point types
//!
template <typename T> struct is_float           { static constexpr bool value = false; };
template <>           struct is_float< float >  { static constexpr bool value = true; };
template <>           struct is_float< double > { static constexpr bool value = true; };

#if __cplusplus >= 201703L
template <typename T> inline constexpr bool is_float_v = is_float< T >::value;
#endif

//!
//! signals complex valued types
//!
template <typename T> struct is_complex_type                            { static constexpr bool value = false; };
template <>           struct is_complex_type< std::complex< float > >   { static constexpr bool value = true; };
template <>           struct is_complex_type< std::complex< double > >  { static constexpr bool value = true; };

#if __cplusplus >= 201703L
template <typename T> inline constexpr bool is_complex_type_v = is_complex_type< T >::value;
#endif

//!
//! signals single/double precision
//!
template <typename T> struct is_single_prec                            { static constexpr bool value = false; };
template <>           struct is_single_prec< float >                   { static constexpr bool value = true; };
template <>           struct is_single_prec< std::complex< float > >   { static constexpr bool value = true; };

template <typename T> struct is_double_prec                            { static constexpr bool value = false; };
template <>           struct is_double_prec< double >                  { static constexpr bool value = true; };
template <>           struct is_double_prec< std::complex< double > >  { static constexpr bool value = true; };

#if __cplusplus >= 201703L
template <typename T> inline constexpr bool is_single_prec_v = is_single_prec< T >::value;
template <typename T> inline constexpr bool is_double_prec_v = is_double_prec< T >::value;
#endif

//!
//! provide real valued type forming base of T
//!
template <typename T> struct real_type                       { using  type_t = T; };
template <typename T> struct real_type< std::complex< T > >  { using  type_t = T; };

template <typename T> using real_type_t = typename real_type< T >::type_t;

//!
//! convert from C++ type to value_type_t
//!
template <typename T> struct value_type                       { static constexpr value_type_t value = real_valued; };
template <typename T> struct value_type< std::complex< T > >  { static constexpr value_type_t value = complex_valued; };

#if __cplusplus >= 201703L
template <typename T> inline constexpr value_type_t value_type_v = value_type< T >::value;
#endif

//
// test if T1 and T2 are of the same type
//
template < typename T1, typename T2 > struct is_same_type           { static constexpr bool  value = false; };
template < typename T1 >              struct is_same_type< T1, T1 > { static constexpr bool  value = true;  };

#if __cplusplus >= 201703L
template < typename T1, typename T2 > inline constexpr bool is_same_type_v = is_same_type< T1, T2 >::value;
#endif

}// namespace HLIB

#endif  // __HLIB_BLAS_TRAITS_HH
