//
// Project     : HLib
// File        : TTFQMR.cc
// Description : class implementing TFQMR
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TTFQMR.hh"

namespace HLIB
{

namespace
{

//
// actual TFQMR iteration
//
template <typename T>
void
tfqmr ( const TLinearOperator *  A,
        TVector *                x,
        const TVector *          b,
        const TLinearOperator *  W,
        TSolverInfo *            info,
        const TSolver *          solver )
{
    using  real_t = typename  real_type< T >::type_t;

    uint    it;
    real_t  norm   = 0;
    real_t  norm_0 = 0;

    if ( info != nullptr )
        info->set_solver( "TFQMR" );
    
    solver->set_start_value( x, b, W );

    auto  t  = b->copy();
    auto  t2 = b->copy();
    auto  rr = b->copy(); // for exact residual
    
    // r0 = W( b - A x )
    auto  r0 = b->copy();

    if ( W != nullptr )
    {
        apply_add( real(-1), A, x, r0.get() );

        if ( solver->use_exact_residual() )
            norm = norm_0 = r0->norm2();
        
        apply( W, r0.get(), t.get() );
        r0->assign( real(1), t.get() );

        if ( ! solver->use_exact_residual() )
            norm = norm_0 = r0->norm2();
    }// if
    else
    {
        apply_add( real(-1), A, x, r0.get() );
        norm = norm_0 = r0->norm2();
    }// else

    auto  w  = r0->copy();
    auto  y  = r0->copy();
    auto  d  = x->copy();

    // v = Ay
    auto  v  = x->copy();

    if ( W != nullptr )
    {
        apply( A, y.get(), d.get() );
        apply( W, d.get(), v.get() );
    }// if
    else
    {
        apply( A, y.get(), v.get() );
    }// else

    t->assign( real(1), v.get() );
    
    d->scale( 0 );
    
    real  tau   = r0->norm2();
    real  theta = 0;
    T     eta   = 0;
    T     rho   = tdot<T>( r0.get(), r0.get() );

    if ( info != nullptr )
        info->append( 0, norm );

    bool  converged = solver->stopped( 0, norm, norm_0, info );
    
    it = 0;
    while ( ! converged )
    {
        auto  sigma = tdot<T>( r0.get(), v.get() );

        if ( sigma == T(0) )
        {
            if ( info != nullptr )
                info->set_status( failed );
            else
                HERROR( ERR_DIV_ZERO, "(TTFQMR) solve", "breakdown" );
            
            return;
        }// if
        
        auto  alpha = rho / sigma;

        if ( alpha == T(0) )
        {
            if ( info != nullptr )
                info->set_status( failed );
            else
                HERROR( ERR_DIV_ZERO, "(TTFQMR) solve", "breakdown" );
            
            return;
        }// if
            
        for ( uint  m = 1; m <= 2; ++m )
        {
            // w = w - α A y
            if ( m == 1 )
            {
                // reuse computed product from last step
            }// if
            else
            {
                // compute Ay (also reused below)
                if ( W != nullptr )
                {
                    apply( A, y.get(), t2.get() );
                    apply( W, t2.get(), t.get() );
                }// if
                else
                {
                    apply( A, y.get(), t.get() );
                }// else
            }// else

            axpy<T>( -alpha, t.get(), w.get() );
            
            auto  theta_old = theta;
            auto  eta_old   = eta;

            if ( tau == real(0) )
            {
                if ( info != nullptr )
                    info->set_status( failed );
                else
                    HERROR( ERR_DIV_ZERO, "(TTFQMR) solve", "breakdown" );
                
                return;
            }// if
            
            theta   = w->norm2() / tau;
            auto  c = real(1) / std::sqrt( 1 + theta * theta );
            tau     = tau * theta * c;
            eta     = c*c * alpha;

            // d = y + (θ²η/α)d
            scale<T>( Math::square( theta_old ) * eta_old / alpha, d.get() );
            d->axpy( real(1), y.get() );

            // x = x + η d
            axpy<T>( eta, d.get(), x );

            // norm (approximation)
            if ( solver->use_exact_residual() )
            {
                apply( A, x, rr.get() );
                rr->axpy( real(-1), b );
                norm = rr->norm2();
            }// if
            else
            {
                norm = std::sqrt( 2*it + m ) * tau;
            }// else

            if ( info != nullptr )
                info->append( 2*it + m, norm );
            
            if ( solver->stopped( 2*it + m, norm, norm_0, info ) )
            {
                converged = true;
                break;
            }// if

            // update y after first step
            if ( m == 1 )
            {
                // y = y - αv
                axpy<T>( -alpha, v.get(), y.get() );
            }// if
        }// for

        if ( rho == T(0) )
        {
            if ( info != nullptr )
                info->set_status( failed );
            else
                HERROR( ERR_DIV_ZERO, "(TTFQMR) solve", "breakdown" );
            
            return;
        }// if
        
        auto  rho_new = tdot<T>( r0.get(), w.get() );
        auto  beta    = rho_new / rho;

        rho = rho_new;
        
        // y' = w + βy
        scale<T>( beta, y.get() );
        axpy<T>( T(1), w.get(), y.get() );
        
        // v = Ay' + β( Ay + β v ) = Ay' + βAy + β²v
        scale<T>( beta*beta, v.get() );
        axpy<T>( beta, t.get(), v.get() );

        if ( W != nullptr )
        {
            apply( A, y.get(), t2.get() );
            apply( W, t2.get(), t.get() );
        }// if
        else
        {
            apply( A, y.get(), t.get() );
        }// else

        axpy<T>( real(1), t.get(), v.get() );
        
        ++it;
    }// while
}

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TTFQMR::TTFQMR ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TTFQMR::~TTFQMR ()
{}

////////////////////////////////////////////////
//
// solving the system
//

//
// the real solving
//
void
TTFQMR::solve ( const TLinearOperator *  A,
                TVector *                x,
                const TVector *          b,
                const TLinearOperator *  W,
                TSolverInfo *            info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TTFQMR) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TTFQMR) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TTFQMR) solve", "A = nullptr" );
    
    if ( A->is_complex() )
        tfqmr<complex>( A, x, b, W, info, this );
    else
        tfqmr<real>( A, x, b, W, info, this );
}

}// namespace HLIB
