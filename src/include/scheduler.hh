#ifndef __HLIB_SCHEDULER_HH
#define __HLIB_SCHEDULER_HH
//
// Project     : HLib
// File        : TListSched.hh
// Description : provides List scheduling algorithm
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/base/types.hh"

namespace HLIB
{

//
// MultiFit scheduling based upon bin-packing
//
class TMFitSched
{
protected:
    // depth of binary search
    const uint  _search_depth;
    
public:
    ///////////////////////////////////////////////
    //
    // constructor
    //

    TMFitSched ( const uint d = 10 ) : _search_depth(d) {}
    
    virtual ~TMFitSched () {}
    
    ///////////////////////////////////////////////
    //
    // scheduling algorithms
    //

    virtual bool schedule ( const uint                     p,
                            std::vector< int > &           sched,
                            const std::vector< double > &  costs ) const;

protected:
    // bin-packing "first fit decreasing" (FFD) algorithm
    uint ffd ( double                         cmax,
               uint                           p,
               const std::vector< double > &  costs,
               std::vector< int > &           sched ) const;

    DISABLE_COPY_OP( TMFitSched );
};

}// namespace HLIB

#endif  // __HLIB_SCHEDULER_HH
