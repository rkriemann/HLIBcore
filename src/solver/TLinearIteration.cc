//
// Project     : HLib
// File        : TLinearIteration.cc
// Description : base-class for all iterative solvers
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <cmath>
#include <sstream>

#include "hpro/matrix/TMatrix.hh"
#include "hpro/solver/TLinearIteration.hh"

namespace HLIB
{

////////////////////////////////////////////////
//
// constructor and destructor
//

TLinearIteration::TLinearIteration ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
        , _damping( 1.0 )
{
}

TLinearIteration::TLinearIteration ( const real              damping,
                                     const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
        , _damping( damping )
{
}

TLinearIteration::~TLinearIteration ()
{}

////////////////////////////////////////////////
//
// solving the system
//

//
// the real solving
//
void
TLinearIteration::solve ( const TLinearOperator *  A,
                          TVector *                x,
                          const TVector *          b,
                          const TLinearOperator *  W,
                          TSolverInfo *            info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TLinearIteration) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TLinearIteration) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TLinearIteration) solve", "A = nullptr" );
    
    real  norm, norm0;
    auto  r = b->copy();

    if ( info != nullptr )
        info->set_solver( "Lin. Iter." );
    
    set_start_value( x, b, W );

    // r = Ax - b
    apply( A, x, r.get() );
    r->axpy( real(-1), b );
    
    norm = norm0 = r->norm2();
    
    if ( info != nullptr )
        info->append( 0, norm );

    //
    // standard linear iteration
    //

    uint  it = 0;
    
    while ( true )
    {
        // iteration step : x_new = x_old - W r
        if ( W != nullptr )
        {
            apply_add( - _damping, W, r.get(), x );
        }// if
        else
        {
            axpy<real>( -_damping, r.get(), x );
        }// else

        it++;

        // update residual
        apply( A, x, r.get() );
        r->axpy( real(-1), b );

        norm = r->norm2();

        if ( info != nullptr )
            info->append( it, norm );

        if ( stopped( it, norm, norm0, info ) )
            break;
    }// while
}

//
// the real solving
//
void
TLinearIteration::solve ( const TLinearOperator *  /* A */,
                          TMatrix *                /* X */,
                          const TMatrix *          /* B */,
                          const TLinearOperator *  /* W */,
                          TSolverInfo *            /* info */ ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

}// namespace HLIB
