//
// Project     : HLib
// File        : TCG.cc
// Description : class implementing CG
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TCG.hh"

namespace HLIB
{

namespace
{

//
// actual CG iteration
//
template <typename T>
void
cg ( const TLinearOperator *  A,
     TVector *                x,
     const TVector *          b,
     const TLinearOperator *  W,
     TSolverInfo *            info,
     const TSolver *          solver )
{
    using  real_t = typename  real_type< T >::type_t;

    uint    it;
    T       lambda, rho, rho_new;
    real_t  norm, norm_0;

    if ( info != nullptr )
        info->set_solver( "CG" );
    
    solver->set_start_value( x, b, W );

    auto  r = b->copy();
    auto  p = x->copy();
    auto  a = x->copy();
    auto  q = x->copy();

    // r = b - A x
    apply_add( real(-1), A, x, r.get() );

    // p = W r
    if ( W != nullptr ) apply( W, r.get(), p.get() );
    else                p->assign( real(1), r.get() );

    // rho = < r, p >
    rho = tdot<T>( r.get(), p.get() );
    norm = norm_0 = r->norm2(); 

    if ( info != nullptr )
        info->append( 0, norm );
    
    it = 0;
    while ( ! solver->stopped( it, norm, norm_0, info ) )
    {
        // a = A p
        apply( A, p.get(), a.get() );
        
        // lambda = rho / < a, p >
        lambda = rho / tdot<T>( a.get(), p.get() );

        // x = x + lambda p; r = r - lambda a
        axpy(  lambda, p.get(), x );
        axpy( -lambda, a.get(), r.get() );

        it++;

        // q = W r
        if ( W != nullptr ) apply( W, r.get(), q.get() );
        else                q->assign( real(1), r.get() );
        
        // rho' = < q, r >
        rho_new = tdot<T>( q.get(), r.get() );

        // check stop criteria
        if ( solver->use_exact_residual() && ( W != nullptr ))
        {
            norm = r->norm2();
        }// if
        else
        {
            norm = Math::abs( Math::sqrt( rho_new ) );
        }// else

        // p = r + (rho' / rho) p
        scale( rho_new / rho, p.get() );
        axpy( T(1), q.get(), p.get() );

        rho = rho_new;

        if ( info != nullptr )
            info->append( it, norm );
    }// while
}

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TCG::TCG ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TCG::~TCG ()
{}

////////////////////////////////////////////////
//
// solving the system
//

//
// the real solving
//
void
TCG::solve ( const TLinearOperator *  A,
             TVector *                x,
             const TVector *          b,
             const TLinearOperator *  W,
             TSolverInfo *            info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TCG) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TCG) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TCG) solve", "A = nullptr" );
    
    if ( ! A->is_self_adjoint() )
        HWARNING( "(TCG) solve : operator not self adjoint" );

    if (( W != nullptr ) && ! W->is_self_adjoint() )
        HWARNING( "(TCG) solve : precondition operator not self adjoint" );

    if ( A->is_complex() )
        cg<complex>( A, x, b, W, info, this );
    else
        cg<real>( A, x, b, W, info, this );
}

}// namespace HLIB
