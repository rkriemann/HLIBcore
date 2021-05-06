//
// Project     : HLib
// File        : TCGS.cc
// Description : class implementing CGS
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TCGS.hh"

namespace HLIB
{

namespace
{

//
// actual CGS iteration
//
template <typename T>
void
cgs ( const TLinearOperator *   A,
      TVector *                 x,
      const TVector *           b,
      const TLinearOperator *   W,
      TSolverInfo *             info,
      const TSolver *           solver )
{
    using  real_t = typename  real_type< T >::type_t;

    uint    it;
    real_t  norm = 0, norm_0 = 0;

    if ( info != nullptr )
        info->set_solver( "CGS" );
    
    solver->set_start_value( x, b, W );

    auto  r  = b->copy();
    auto  r0 = b->copy();
    auto  rt = b->copy(); // for exact residual

    // r0 = W( b - A x )
    apply_add( real(-1), A, x, r0.get() );

    if ( solver->use_exact_residual() )
        norm_0 = r0->norm2();

    if ( W != nullptr )
    {
        apply( W, r0.get(), r.get() );
        r0->assign( real(1), r.get() );
    }// if
    else
        r->assign( real(1), r0.get() );

    if ( ! solver->use_exact_residual() )
        norm_0 = r0->norm2();

    auto  p  = r->copy();
    auto  u  = r->copy();
    auto  q  = r->copy();
    auto  t  = r->copy();
    
    norm = norm_0;

    if ( info != nullptr )
        info->append( 0, norm );
    
    it = 0;
    while ( true )
    {
        // t = Ap
        if ( W != nullptr )
        {
            apply( A, p.get(), q.get() );
            apply( W, q.get(), t.get() );
        }// if
        else
        {
            apply( A, p.get(), t.get() );
        }// else

        // α = <rᵢ, r₀> / <Apᵢ, r₀>
        auto  sigma = tdot<T>( t.get(), r0.get() );
        auto  rho   = tdot<T>( r.get(), r0.get() );
        auto  alpha = rho / sigma;

        // q = u - αAp
        q->assign( 1.0, u.get() );
        axpy<T>( -alpha, t.get(), q.get() );

        // u = u + q
        axpy<T>( 1.0, q.get(), u.get() );
        
        // x = x + α(u + q)
        axpy<T>( alpha, u.get(), x );

        // r = r - α A (u + q)
        if ( W != nullptr )
        {
            apply( A, u.get(), t.get() );
            apply_add( -alpha, W, t.get(), r.get() );
        }// if
        else
        {
            apply_add( -alpha, A, u.get(), r.get() );
        }// else

        // one step per MVM, so two steps here
        it += 2;

        if ( solver->use_exact_residual() && ( W != nullptr ))
        {
            apply( A, x, rt.get() );
            rt->axpy( real(-1), b );
            norm = rt->norm2();
        }// if
        else
        {
            norm = r->norm2();
        }// else
        
        if ( info != nullptr )
            info->append( it, norm );
        
        if ( solver->stopped( it, norm, norm_0, info ) )
            break;

        
        auto  beta = tdot<T>( r.get(), r0.get() ) / rho;

        // u = r + βq
        u->assign( 1.0, r.get() );
        axpy<T>( beta, q.get(), u.get() );
        
        // p = u + βq + β²p
        scale<T>( beta*beta, p.get() );
        axpy<T>( beta, q.get(), p.get() );
        axpy<T>( T(1), u.get(), p.get() );
    }// while
}

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TCGS::TCGS ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TCGS::~TCGS ()
{}

////////////////////////////////////////////////
//
// solving the system
//

//
// the real solving
//
void
TCGS::solve ( const TLinearOperator *  A,
             TVector *                 x,
             const TVector *           b,
             const TLinearOperator *   W,
             TSolverInfo *             info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TCGS) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TCGS) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TCGS) solve", "A = nullptr" );
    
    if ( A->is_complex() )
        cgs<complex>( A, x, b, W, info, this );
    else
        cgs<real>( A, x, b, W, info, this );
}

}// namespace HLIB
