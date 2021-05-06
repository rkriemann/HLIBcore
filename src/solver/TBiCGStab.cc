//
// Project     : HLib
// File        : TBiCGStab.cc
// Description : class implementing BiCG-Stab
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/base/System.hh"

#include "hpro/solver/TBiCGStab.hh"

namespace HLIB
{

namespace
{

//
// actual BiCG-Stab algorithm
//
template <typename T>
void
bicgstab ( const TLinearOperator *   A,
           TVector *                 x,
           const TVector *           b,
           const TLinearOperator *   W,
           TSolverInfo *             info,
           const TSolver *           solver );

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TBiCGStab::TBiCGStab ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TBiCGStab::~TBiCGStab ()
{}

////////////////////////////////////////////////
//
// solving the system
//

//
// the real solving
//
void
TBiCGStab::solve ( const TLinearOperator *  A,
                   TVector *                x,
                   const TVector *          b,
                   const TLinearOperator *  W,
                   TSolverInfo *            info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TBiCGStab) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TBiCGStab) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TBiCGStab) solve", "A = nullptr" );

    if ( A->is_complex() )
        bicgstab<complex>( A, x, b, W, info, this );
    else
        bicgstab<real>( A, x, b, W, info, this );
}

namespace
{

//
// actual BiCG-Stab algorithm
//
template <typename T>
void
bicgstab ( const TLinearOperator *   A,
           TVector *                 x,
           const TVector *           b,
           const TLinearOperator *   W,
           TSolverInfo *             info,
           const TSolver *           solver )
{
    using  real_t = typename real_type< T >::type_t;
    
    const real_t  zero_limit = Math::square( Limits::epsilon< real_t >() );
    
    if ( info != nullptr )
        info->set_solver( "BiCGStab" );
    
    solver->set_start_value( x, b, W );

    // r = b - Ax
    auto  r  = b->copy();
    auto  t  = b->copy(); // for exact residual

    real  norm0 = real(0);

    if ( W != nullptr )
    {
        apply_add( real(-1), A, x, t.get() );

        if ( solver->use_exact_residual() )
            norm0 = t->norm2();
        
        apply( W, t.get(), r.get() );

        if ( ! solver->use_exact_residual() )
            norm0 = r->norm2();
    }// if
    else
    {
        apply_add( real(-1), A, x, r.get() );
        norm0 = r->norm2();
    }// else

    auto  norm  = norm0;

    // choose r0, with <r0,r> != 0
    auto  r0   = r->copy();

    // create/set up remaining aux. vectors
    auto  p  = r0->copy();
    auto  Ap = p->copy();
    auto  Ar = r->copy();
    
    //
    // iterate until stop-criterion reached
    //

    if ( info != nullptr )
        info->append( 0, norm0 );

    if ( solver->stopped( 0, norm, norm0, info ) )
        return;
    
    uint  it = 1;
    
    while ( true )
    {
        // α = <r_j, r_0> / <Ap_j, r_0>
        if ( W != nullptr )
        {
            apply( A, p.get(), t.get() );
            apply( W, t.get(), Ap.get() );
        }// if
        else
        {
            apply( A, p.get(), Ap.get() );
        }// else

        auto  r_r0  = tdot<T>( r.get(),  r0.get() );
        auto  Ap_r0 = tdot<T>( Ap.get(), r0.get() );
        
        if (( Math::abs( Ap_r0 ) < zero_limit ) || ( Math::abs( r_r0 ) < zero_limit ))
        {
            if ( info != nullptr ) info->set_status( failed );
            else                   HERROR( ERR_DIV_ZERO, "(TBiCGStab) solve", "breakdown" );
            
            return;
        }// if
        
        auto  alpha = r_r0 / Ap_r0;

        // x = x + αp
        axpy<T>( alpha, p.get(), x );

        // r = r - alpha A p_j
        axpy<T>( -alpha, Ap.get(), r.get() );

        if ( solver->use_exact_residual() && ( W != nullptr ))
        {
            apply( A, x, t.get() );
            t->axpy( real(-1), b );
            norm = t->norm2();
        }// if
        else
        {
            norm = r->norm2();
        }// else
        
        if ( info != nullptr )
            info->append( it, norm );
        
        if ( solver->stopped( it, norm, norm0, info ) )
            break;

        it++;

        // ω = <Ar_j, s_j> / <Ar_j, Ar_j>
        if ( W != nullptr )
        {
            apply( A, r.get(), t.get() );
            apply( W, t.get(), Ar.get() );
        }// if
        else
        {
            apply( A, r.get(), Ar.get() );
        }// else

        auto  Ar_r  = tdot<T>( Ar.get(), r.get() );
        auto  Ar_Ar = tdot<T>( Ar.get(), Ar.get() );
        auto  omega = Ar_r / Ar_Ar;

        if ( Math::abs( omega ) < zero_limit )
        {
            if ( info != nullptr ) info->set_status( failed );
            else                   HERROR( ERR_DIV_ZERO, "(TBiCGStab) solve", "breakdown" );
            
            return;
        }// if
        
        // x = x + ωr
        axpy<T>( omega, r.get(), x );

        // r = r - ωAr
        axpy<T>( - omega, Ar.get(), r.get() );

        if ( solver->use_exact_residual() && ( W != nullptr ))
        {
            apply( A, x, t.get() );
            t->axpy( real(-1), b );
            norm = t->norm2();
        }// if
        else
        {
            norm = r->norm2();
        }// else
        
        if ( info != nullptr )

            info->append( it, norm );
        
        if ( solver->stopped( it, norm, norm0, info ) )
            break;

        it++;

        // β = <r_j+1,r0> / <r_j,r0> · alpha_j / omega_j
        auto  beta = ( tdot<T>( r.get(), r0.get() ) / r_r0 ) * ( alpha / omega );

        // p = r_j+1 + β( p_j - ωAp_j )
        scale<T>( beta, p.get() );
        axpy<T>( -omega * beta, Ap.get(), p.get() );
        axpy<T>( T(1), r.get(), p.get() );
    }// while
}

}// namespace anonymous

}// namespace HLIB
