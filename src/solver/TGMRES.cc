//
// Project     : HLib
// File        : TGMRES.cc
// Description : class implementing GMRES
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/base/System.hh"
#include "hpro/blas/Algebra.hh"

#include "hpro/solver/TGMRES.hh"

namespace HLIB
{

using std::vector;
using std::unique_ptr;

// namespace abbr.
namespace B = BLAS;

//
// abbreviation for applying preconditioner or just copy vector
//
void
apply ( const TLinearOperator *  M,
        const TVector &          v_src,
        TVector &                v_dest )
{
    if ( M != nullptr ) apply( M, & v_src, & v_dest );
    else                v_dest.assign( real(1), & v_src );
}

namespace
{

//
// epsilon for check of (almost) linear dependency in Krylov subspace
//
const real GMRES_EPS = Limits::epsilon<real>() * real(5e1);

////////////////////////////////////////////////
//
// apply plane rotation (  cs  sn ) to (a)
//                      ( -sn  cs )    (b)
//

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
// the real solving
//
template <typename T>
void
gmres ( const uint               restart,
        const TLinearOperator *  A,
        TVector *                x,
        const TVector *          b,
        const TLinearOperator *  W,
        TSolverInfo *            info,
        const TSolver *          solver )
{
    using  real_t = typename real_type< T >::type_t;

    uint                             j, it;
    real_t                           norm = 0, norm_0 = 0, norm_old = 0;
    auto                             d = b->copy();
    B::Vector< T >                   y( restart ), e1( restart+1 );
    B::Vector< T >                   givc( restart+1 ), givs( restart+1 );
    B::Matrix< T >                   H( restart+1, restart );
    vector< unique_ptr< TVector > >  V( restart+1 );

    if ( info != nullptr )
        info->set_solver( "GMRES(" + to_string( "%d", restart ) + ")" );
    
    solver->set_start_value( x, b, W );

    //
    // GMRES iteration
    //
    
    int  ndiv = 0;     // counter for divergent inner iterations
        
    it = 0;
    
    while ( true )
    {
        //
        // v_0 = r / | r | = W ( b - Ax ) / |...|
        //
        
        if ( V[0] == nullptr )
            V[0] = b->copy();
        
        d->assign( real(1), b );
        apply_add( real(-1), A, x, d.get() );

        if ( solver->use_exact_residual() )
        {
            norm = d->norm2();
        }// if
        
        apply( W, *d, *V[0] );

        if ( ! solver->use_exact_residual() )
        {
            norm = V[0]->norm2();
        }// if

        if ( it == 0 )
        {
            norm_0 = norm_old = norm;

            if ( info != nullptr )
                info->append( 0, norm );
        }// if
        
        // test stop condition
        if ( solver->stopped( it, norm, norm_0, info ) )
            break;

        // check for divergent inner iterations
        if ( norm > norm_old ) ndiv++;
        else                   ndiv = 0;

        // stop after 10 divergent inner iterations
        if ( ndiv > 10 )
        {
            HWARNING( "in (TGMRES) solve : residual continuously growing (10 inner it.)" );
            break;
        }// if

        auto  norm_V0 = norm;
        
        // adjust norm_V0 if neccessary
        if ( solver->use_exact_residual() && ( W != nullptr ))
            norm_V0 = V[0]->norm2();
        
        if ( norm_V0 != real(0) ) scale( T(1) / norm_V0, V[0].get() );
        else                      break;

        //
        // e1 = | r | (1,0,0,...)^T
        //
        
        B::fill( T(0), e1 );
        e1(0) = norm_V0;
        
        //
        // iterate in Krylov-subspace
        //

        uint  dim = restart;
        
        for ( j = 0; j < restart; j++ )
        {
            //
            // v_j+1 = W A v_j
            //
            
            if ( V[j+1] == nullptr )
                V[j+1] = b->copy();
        
            apply( A, *V[j], *d );
            apply( W, *d,    *V[j+1] );

            norm = V[j+1]->norm2();
            
            //
            // orthogonalize wrt. previously computed basis-vectors
            //
            
            for ( uint i = 0; i <= j; i++ )
            {
                H(i,j) = tdot<T>( V[i].get(), V[j+1].get() );
                axpy( -H(i,j), V[i].get(), V[j+1].get() );
            }// for
            
            H(j+1,j) = V[j+1]->norm2();

            //
            // check if new direction is in span of previous directions
            // and stop inner iteration if it is the case
            //
            
            if ( Math::abs( H(j+1,j) ) < norm * GMRES_EPS )
            {
                dim = j+1;
                it++;
                break;
            }// if
            
            scale( T(1) / H(j+1,j), V[j+1].get() );
            
            //
            // apply old givens rotations to new column
            //
            
            for ( uint i = 0; i < j; i++ )
                apply_rot( givc(i), givs(i), H(i,j), H(i+1,j) );

            //
            // compute and apply givens rotation for new element
            //

            Math::givens( H(j,j), H(j+1,j), givc(j), givs(j) );
            apply_rot( givc(j), givs(j), H(j,j), H(j+1,j) );
            H(j+1,j) = real(0);

            //
            // apply givens rotation to RHS (new pos is 0 in e1)
            //
            
            e1(j+1) = - Math::conj( givs(j) ) * e1(j);
            e1(j)   =   givc(j)               * e1(j);


            if ( solver->use_exact_residual() )
            {
                //
                // solve min_y | beta e^1 - H y |
                // (H is upper triangular)
                //

                const uint      ldim = j+1;
                B::Matrix< T >  H_dim(  H,  B::Range( 0, ldim-1 ), B::Range( 0, ldim-1 ), copy_value );
                B::Vector< T >  e1_dim( e1, B::Range( 0, ldim-1 ), copy_value );
                
                B::solve_tri( B::upper_triangular, B::general_diag, H_dim, e1_dim );

                auto  x_copy = x->copy();
                
                for ( uint i = 0; i < ldim; i++ )
                    axpy( e1_dim(i), V[i].get(), x_copy.get() );

                apply( A, x_copy.get(), d.get() );
                d->axpy( real(-1), b );
                
                norm = d->norm2();
            }// if
            else
            {
                // last coefficient is actual residual norm
                norm = Math::abs( e1(j+1) );
            }// else

            it++;
            
            if ( info != nullptr )
                info->append( it, norm );
            
            //
            // test stop condition
            //

            if ( solver->stopped( it, norm, norm_0, info ) )
            {
                dim = j+1;
                break;
            }// if

            norm_old = norm;
        }// for

        //
        // solve min_y | beta e^1 - H y |
        // (H is upper triangular)
        //

        B::Matrix< T >  H_dim(  H,  B::Range( 0, dim-1 ), B::Range( 0, dim-1 ) );
        B::Vector< T >  e1_dim( e1, B::Range( 0, dim-1 ) );
        
        B::solve_tri( B::upper_triangular, B::general_diag, H_dim, e1_dim );
        
        //
        // adjust x
        //

        for ( uint i = 0; i < dim; i++ )
            axpy( e1(i), V[i].get(), x );

        // test stop condition
        if ( solver->stopped( it, norm, norm_0, info ) )
            break;
    }// while
}

}// namespace anonymous

////////////////////////////////////////////////
//
// constructor and destructor
//

TGMRES::TGMRES ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
        , _restart( CFG::Solver::gmres_restart )
{}

TGMRES::TGMRES ( const uint              restart,
                 const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
        , _restart( restart )
{}

TGMRES::~TGMRES ()
{}

////////////////////////////////////////////////
//
// solving the system
//

//
// the real solving
//
void
TGMRES::solve ( const TLinearOperator *  A,
                TVector *                x,
                const TVector *          b,
                const TLinearOperator *  W,
                TSolverInfo *            info ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TGMRES) solve", "x = nullptr" );
    if ( b == nullptr ) HERROR( ERR_ARG, "(TGMRES) solve", "b = nullptr" );
    if ( A == nullptr ) HERROR( ERR_ARG, "(TGMRES) solve", "A = nullptr" );

    if ( A->is_complex() )
        gmres<complex>( _restart, A, x, b, W, info, this );
    else
        gmres<real>( _restart, A, x, b, W, info, this );
}

}// namespace
