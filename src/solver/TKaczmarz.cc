//
// Project     : HLib
// File        : TKaczmarz.cc
// Description : iterative solver based on Kaczmarz algorithm
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TKaczmarz.hh"

namespace HLIB
{

////////////////////////////////////////////////
//
// constructor and destructor
//

TKaczmarz::TKaczmarz ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TKaczmarz::~TKaczmarz ()
{}

////////////////////////////////////////////////
//
// solving the system
//

//
// the actual Kaczmarz iteration
// (see Hackbusch: "Large iterative Systems", Sec 8.2)
//
void
TKaczmarz::solve ( const TLinearOperator *  A,
                   TVector *                x,
                   const TVector *          b,
                   const TLinearOperator *  W,
                   TSolverInfo *            info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TKaczmarz) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TKaczmarz) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TKaczmarz) solve", "A = nullptr" );
    
    real                        norm, norm0;
    std::unique_ptr< TVector >  r( b->copy() );
    std::unique_ptr< TVector >  a( b->copy() );
    std::unique_ptr< TVector >  e( b->copy() );
    std::unique_ptr< TVector >  t( b->copy() );

    if ( info != nullptr )
        info->set_solver( "Kaczmarz" );

    set_start_value( x, b, W );    

    // r = A x - b
    apply( A, x, t.get() );
    t->axpy( real(-1), b );
    if ( W != nullptr ) apply( W, t.get(), r.get() );
    else                r->assign( real(1), t.get() );

    norm = norm0 = r->norm2();

    if ( info != nullptr )
        info->append( 0, norm );

    e->scale( real(0) );

    //
    // Kaczmarz iteration
    //

    uint        it = 0;
    const uint  n  = uint(e->size());
    
    while ( ! stopped( it, norm, norm0, info ) )
    {
        // determine <it>'th column of A
        e->set_entry( it % n, real(1) );

#if 1
        //
        // left transformation
        //
        
        // a = A e
        apply( A, e.get(), t.get() );
        if ( W != nullptr ) apply( W, t.get(), a.get() );
        else                a->assign( real(1), t.get() );
        
        // <Ax-b, a>
        const complex r_dot_a = dot( r.get(), a.get() );

        // <a, a>
        const complex a_dot_a = dot( a.get(), a.get() );

        // if <a,a> == 0, then A is not regular
        if ( abs( a_dot_a ) < 1e-32 )
            HERROR( ERR_MAT_SINGULAR, "(TKaczmarz) solve",
                    to_string( "e_%d in kern(A)", it % n ) );
        
        // x = x - e <r,a>/<a,a>
        x->caxpy( - r_dot_a / a_dot_a, e.get() );
#else
        //
        // right transformation
        //
        
        // a = A^H e
        apply( A, e.get(), a.get(), MATOP_ADJ );
        
        // <Ax-b, e>
        const complex r_dot_e = dot( r.get(), e.get() );

        // <a, a>
        const complex a_dot_a = dot( a.get(), a.get() );

        // if <a,a> == 0, then A is not regular
        if ( abs( a_dot_a ) < 1e-32 )
            HERROR( ERR_MAT_SINGULAR, "(TKaczmarz) solve",
                    to_string( "e_%d in kern(A)", it+1 ) );
        
        // x = x - a <r,e>/<a,a>
        x->caxpy( - r_dot_e / a_dot_a, a.get() );
#endif
        
        e->set_entry( it % n, real(0) );

        // update residual
        apply( A, x, t.get() );
        t->axpy( real(-1), b );
        if ( W != nullptr ) apply( W, t.get(), r.get() );
        else             r->assign( real(1), t.get() );

        norm = r->norm2();
        
        // Sys::Out::printfln( "|r_%d| = %.8e", it+1, norm );

        it++;
        if ( info != nullptr )
            info->append( it+1, norm );
    }// while
}

}// namespace
