#ifndef __HLIB_TBSHMBUILDER_HH
#define __HLIB_TBSHMBUILDER_HH
//
// Project     : HLib
// File        : TBSHMBuilder.hh
// Description : class for building h-matrices out of bytestreams
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/base/TByteStream.hh"

#include "hpro/matrix/TMatrix.hh"

namespace HLIB
{

//
// takes a bytestream and reads in matrix
//
class TBSHMBuilder
{
public:
    //
    // constructor and destructor
    //

    TBSHMBuilder () {}

    virtual ~TBSHMBuilder () {}

    //
    // construct matrix out of given bytestream
    //

    virtual std::unique_ptr< TMatrix >  build ( TByteStream & bs ) const;

protected:
    //
    // special functions
    //

    // return matrix corresponding to given type
    std::unique_ptr< TMatrix >  build_matrix ( uint t ) const;
};

}// namespace

#endif // __HLIB_TBSHMBUILDER_HH
