#ifndef __HLIB_BLAS_RANGE_HH
#define __HLIB_BLAS_RANGE_HH
//
// Project     : HLib
// File        : Range.hh
// Description : implements class for index set with stride
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <ostream>

#include "hpro/base/error.hh"
#include "hpro/cluster/TIndexSet.hh"

namespace HLIB
{

namespace BLAS
{

//!
//! \ingroup  BLAS_Module
//! \class    Range
//! \brief    indexset with modified ctors
//!
class Range : public TIndexSet
{
public:
    //
    // constructors
    // - last < first => indexset = ∅
    //

    //! create index set { \a afirst ... \a last } with stride \a astride
    Range ( const idx_t   afirst,
            const idx_t   alast ) noexcept
            : TIndexSet( afirst, alast )
    {
    }

    //! create index set { \a pos }
    Range ( const idx_t  pos ) noexcept
            : TIndexSet( pos, pos )
    {}

    //! copy constructor for TIndexSet objects
    Range ( const TIndexSet &  is ) noexcept
            : TIndexSet( is )
    {}
            
    //
    // access data
    //

    //! return stride of index set
    size_t  stride () const noexcept { return 1; }

    //
    // output
    //

    //! stream output
    friend std::ostream & operator << ( std::ostream & os, const Range & r )
    {
        os << "[ " << r.first() << " , " << r.last() << " ]";
        return os;
    }

public:
    static Range all;
};

}// namespace BLAS

}// namespace HLIB

#endif  // __HLIB_BLAS_RANGE_HH
