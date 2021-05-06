//
// Project     : HLib
// File        : TUpdateAccumulator.cc
// Description : class for handling updates to matrix blocks
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TBlockMatrix.hh"

#include "TRNG.hh"
#include "list.hh"

#include "hpro/matrix/structure.hh"

namespace HLIB
{

namespace B = BLAS;

using std::deque;
using std::pair;
using std::list;
using std::make_pair;
using std::unique_ptr;

//
// initialise matrix for accumulated updates
//
void
TUpdateAccumulator::init ( const TMatrix *  M )
{
    TScopedLock  mlock( *this );

    if ( _accumulated.get() == nullptr )
    {
        if ( is_dense( M ) )
        {
            _accumulated = std::make_unique< TDenseMatrix >( M->row_is(), M->col_is(), M->value_type() );
        }// if
        else
        {
            _accumulated = std::make_unique< TRkMatrix >( M->row_is(), M->col_is(), M->value_type() );
        }// else
    }// if
}

//
// return true if accumulator holds any updates
//
bool
TUpdateAccumulator::has_updates () const
{
    if ( pending_direct().empty() && pending_recursive().empty() && ( accumulated_updates() == nullptr ))
        return false;

    return true;
}
    
//
// compute sum of local direct updates
//
void
TUpdateAccumulator::apply_direct ( const TTruncAcc &  /* acc */,
                                   TMatrix *          /* M */,
                                   const bool         /* update_M */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}
    

//
// remove list of pending updates
//
void
TUpdateAccumulator::clear_pending ()
{
    TScopedLock  mlock( *this );
    
    for ( auto  upd : _pending_direct )
        delete upd;
    
    _pending_direct.clear();

    for ( auto  upd : _pending_recursive )
        delete upd;
    
    _pending_recursive.clear();
}

//
// add update matrix
//
void
TUpdateAccumulator::add_update ( const TMatrix *    /* M */,
                                 const TTruncAcc &  /* acc */,
                                 const TMatrix *    /* dest */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// add update from parent matrix
//
void
TUpdateAccumulator::add_parent_update ( const TMatrix *    /* M */,
                                        const TTruncAcc &  /* acc */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

}// namespace HLIB
