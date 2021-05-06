//
// Project     : HLib
// File        : TBSHMBuilder.cc
// Description : class for building h-matrices out of bytestreams
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <iostream>

#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/matrix/TBlockMatrix.hh"
#include "hpro/matrix/THMatrix.hh"

#include "hpro/matrix/TBSHMBuilder.hh"

namespace HLIB
{

using std::make_unique;

/////////////////////////////////////////////////////////////////////////////
//
// local defines
//
/////////////////////////////////////////////////////////////////////////////

// define to do print detailed informations about computations
#define LOG( msg )   // std::cout << msg << std::endl << std::flush;

/////////////////////////////////////////////////////////////////////////////
//
// TBSHMBuilder
//
/////////////////////////////////////////////////////////////////////////////

//
// construct matrix out of given bytestream
//

std::unique_ptr< TMatrix >
TBSHMBuilder::build ( TByteStream & bs ) const
{
    //
    // read type and build corresponding sub-matrix
    //

    uint  t;

    bs.get( & t, sizeof(uint) );

    std::unique_ptr< TMatrix >  A( build_matrix( t ) );

    // check if type of matrix is known
    if ( A.get() == nullptr )
        HERROR( ERR_NULL, "(TBSHMBuilder) build", "" );

    // use member function to build matrix
    A->build( bs );

    return A;
}

//
// return matrix corresponding to given type
//
std::unique_ptr< TMatrix >
TBSHMBuilder::build_matrix ( uint t ) const
{
    std::unique_ptr< TMatrix >  M;

    if      ( t == TYPE_ID( TDenseMatrix ) ) M = make_unique< TDenseMatrix >();
    else if ( t == TYPE_ID( TRkMatrix )    ) M = make_unique< TRkMatrix >();
    else if ( t == TYPE_ID( TBlockMatrix ) ) M = make_unique< TBlockMatrix >();
    else if ( t == TYPE_ID( THMatrix )     ) M = make_unique< THMatrix >();
    else
        HERROR( ERR_MAT_TYPE, "(TBSHMBuilder) build_matrix", RTTI::id_to_type( t ) );

    LOG( "created " << M->typestr() );

    return M;
}

}// namespace
