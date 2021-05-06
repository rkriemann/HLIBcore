#ifndef __HLIB_UNORDERED_MAP_HH
#define __HLIB_UNORDERED_MAP_HH
//
// Project     : HLIBpro
// File        : unordered_map.hh
// Description : provide map class (unordered_map preferred)
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <tuple>
#include <unordered_map>

#include "hpro/cluster/TIndexSet.hh"

/////////////////////////////////////////////////////////////////////////////
//
// index set hash functions
//
/////////////////////////////////////////////////////////////////////////////

namespace std
{

template <>
struct hash< std::pair< int, int > >
{
    size_t operator () ( const std::pair< int, int > &  arg ) const
    {
        const auto  a = arg.first;
        const auto  b = arg.second;
        
        return  ( a >= b
                  ? a * a + a + b
                  : a + b * b );
    }
};

template <>
struct hash< HLIB::TIndexSet >
{
    size_t operator () ( const HLIB::TIndexSet &  is ) const
    {
        return ( std::hash< HLIB::idx_t >()( is.first() ) +
                 std::hash< HLIB::idx_t >()( is.last()  ) );
    }
};

template <>
struct hash< HLIB::TBlockIndexSet >
{
    size_t operator () ( const HLIB::TBlockIndexSet &  is ) const
    {
        return ( std::hash< HLIB::TIndexSet >()( is.row_is() ) +
                 std::hash< HLIB::TIndexSet >()( is.col_is() ) +
                 std::hash< HLIB::idx_t >()( is.row_is().first() ) ); // ← to make it unique w.r.t. transposition
    }
};

}// namespace std

#endif  // __HLIB_UNORDERED_MAP_HH
