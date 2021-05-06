//
// Project     : HLib
// File        : TMINRES.cc
// Description : class implementing MINRES
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TMINRES.hh"

namespace HLIB
{

namespace
{

//
// epsilon for check of (almost) linear dependency in Krylov subspace
//
const real MINRES_EPS = Limits::epsilon<real>() * real(5e1);

//
// apply plane rotation (  cs   sn ) to (a)
//                      ( -sn'  cs )    (b)
//
template <typename T>
void
apply_rot ( const T  cs,
            const T  sn,
            T &      a,
            T &      b )
{
    const T  ta = a;
    const T  tb = b;

    a =   cs             * ta + sn * tb;
    b = - Math::conj(sn) * ta + cs * tb;
}

//
// MINRES algorithm
//
template <typename T>
void
minres ( const TLinearOperator *   A,
         TVector *                 x,
         const TVector *           b,
         const TLinearOperator *   W,
         TSolverInfo *             info,
         const TSolver *           solver )
{
    using  real_t = typename real_type< T >::type_t;
    
    solver->set_start_value( x, b, W );

    auto  p    = b->copy();
    auto  p1   = b->copy();
    auto  p2   = b->copy();
    auto  v    = b->copy();
    auto  v1   = b->copy();
    auto  vbar = b->copy();
    auto  rr   = b->copy(); // for exact residual

    //
    // v_0 = r / | r | = W ( Ax - b ) / |...|
    //
        
    real_t   norm  = 0;
    real_t   norm0 = 0;
    T        e;
    
    if ( W != nullptr )
    {
        apply( A, x, p.get() );
        p->axpy( real(-1), b );

        if ( solver->use_exact_residual() )
            norm0 = p->norm2();
        
        apply( W, p.get(), v.get() );

        if ( ! solver->use_exact_residual() )
            norm0 = v->norm2();
    }// if
    else
    {
        apply( A, x, v.get() );
        v->axpy( real(-1), b );
        norm0 = v->norm2();
    }// else

    norm = norm0;

    if ( info != nullptr )
        info->append( 0, norm );
        
    // test stop condition
    if ( solver->stopped( 0, norm, norm0, info ) )
        return;

    auto  norm_v = norm;

    // adjust norm_v if neccessary
    if ( solver->use_exact_residual() && ( W != nullptr ))
        norm_v = v->norm2();
    
    v->scale( real(1) / norm_v );

    e = norm_v;
        
    //
    // MINRES iteration
    //
    
    uint  it  = 0;
    T     cs1 = real(0), sn1 = real(0);
    T     cs2 = real(0), sn2 = real(0);

    while ( ! solver->stopped( it, norm, norm0, info ) )
    {
        //
        // v_j+1 = W A v_j
        //
            
        if ( W != nullptr )
        {
            apply( A, v.get(), p.get()  );
            apply( W, p.get(), vbar.get() );
        }// if
        else
            apply( A, v.get(), vbar.get() );

        norm = vbar->norm2();
            
        //
        // compute scalar products with previous directions
        //

        T  H;
        T  H1 = real(0);
        T  H2 = real(0);

        if ( it > 0 )
            H1 = tdot<T>( v1.get(), vbar.get() );
        H = tdot<T>( v.get(), vbar.get() );

        //
        // orthogonalize wrt. previously computed basis-vectors
        //
        
        if ( it > 0 ) axpy( -H1, v1.get(), vbar.get() );
        axpy( -H, v.get(), vbar.get() );
        
        T  Hbar = vbar->norm2();

        //
        // check if new direction is in span of previous directions
        // and stop inner iteration if so
        //
        
        if ( Math::abs( Hbar ) < norm * MINRES_EPS )
            break;
        
        scale( T(1) / Hbar, vbar.get() );
            
        //
        // apply old Givens rotations to new column
        //
        
        if ( it > 1 ) apply_rot( cs2, sn2, H2, H1 );
        if ( it > 0 ) apply_rot( cs1, sn1, H1, H  );
            
        //
        // compute and apply Givens rotation for new element
        //

        T  cs, sn;
        
        Math::givens( H, Hbar, cs, sn );
        apply_rot( cs, sn, H, Hbar );

        //
        // apply Givens rotation to RHS (new pos is 0 in e1)
        //

        T  ebar;

        ebar = - sn * e;
        e    =   cs * e;

        //
        // adjust x
        //
        
        p->assign( real(1), v.get() );
        if ( it > 0 ) axpy( - H1, p1.get(), p.get() );
        if ( it > 1 ) axpy( - H2, p2.get(), p.get() );
        scale( T(1) / H, p.get() );

        axpy( -e, p.get(), x );

        if ( solver->use_exact_residual() && ( W != nullptr ))
        {
            apply( A, x, rr.get() );
            rr->axpy( real(-1), b );
            norm = rr->norm2();
        }// if
        else
        {
            // last coefficient is actual residual norm
            norm = Math::abs( ebar );
        }// else
                
        //
        // copy values for next step
        //
        
        v1->assign( real(1), v.get()  );
        v->assign(  real(1), vbar.get() );
        p2->assign( real(1), p1.get() );
        p1->assign( real(1), p.get()  );

        cs2 = cs1;
        sn2 = sn1;
        cs1 = cs;
        sn1 = sn;

        e  = ebar;

        it++;
        
        if ( info != nullptr )
            info->append( it, norm );
    }// for
}

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TMINRES::TMINRES ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TMINRES::~TMINRES ()
{}

////////////////////////////////////////////////
//
// solving the system
//

//
// MINRES iteration
// - symmetric and optimized version if GMRES
//
void
TMINRES::solve ( const TLinearOperator *  A,
                 TVector *                x,
                 const TVector *          b,
                 const TLinearOperator *  W,
                 TSolverInfo *            info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TMINRES) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TMINRES) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TMINRES) solve", "A = nullptr" );
    
    if ( info != nullptr )
        info->set_solver( "MINRES" );

    if ( ! A->is_self_adjoint() )
        HWARNING( "(TMINRES) solve : operator not self adjoint" );

    if (( W != nullptr ) && ! W->is_self_adjoint() )
        HWARNING( "(TMINRES) solve : precondition operator not self adjoint" );

    if ( A->is_complex() )
        minres<complex>( A, x, b, W, info, this );
    else
        minres<real>( A, x, b, W, info, this );
}

namespace
{

}// namespace anonymous

}// namespace HLIB
