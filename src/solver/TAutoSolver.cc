//
// Project     : HLib
// File        : TAutoSolver.cc
// Description : class implementing a solver which automatically decides best strategy
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hlib-solver.hh"

#include "hpro/solver/TAutoSolver.hh"

namespace HLIB
{

////////////////////////////////////////////////
//
// constructor and destructor
//

TAutoSolver::TAutoSolver ( const TStopCriterion &  stop_crit )
        : TSolver( stop_crit )
{}

TAutoSolver::~TAutoSolver ()
{}

////////////////////////////////////////////////
//
// solving the system
//

void
TAutoSolver::solve ( const TLinearOperator *  A,
                     TVector *                x,
                     const TVector *          b,
                     const TLinearOperator *  W,
                     TSolverInfo *            info ) const
{
    //
    // if symmetric or hermition, use MINRES, else GMRES
    //
    
    if ( W == nullptr )
    {
        if ( A->is_self_adjoint() )
        {
            TMINRES  solver( _stop_criterion );

            solver.initialise_start_value( _initialise_start_value );
            
            solver.solve( A, x, b, W, info );
        }// if
        else
        {
            TGMRES  solver( _stop_criterion );
            
            solver.initialise_start_value( _initialise_start_value );
            
            solver.solve( A, x, b, W, info );
        }// else
    }// if
    else
    {
        if ( A->is_self_adjoint()  &&  W->is_self_adjoint() )
        {
            TMINRES  solver( _stop_criterion );
            
            solver.initialise_start_value( _initialise_start_value );
            
            solver.solve( A, x, b, W, info );
        }// if
        else
        {
            #if 1

            TGMRES  solver( _stop_criterion );
            
            solver.initialise_start_value( _initialise_start_value );
            
            solver.solve( A, x, b, W, info );

            #else
            
            // always store history to be used as reference
            const bool  store_hist = true; // ( info != nullptr ? info->store_hist() : false );
            const bool  print_hist = ( info != nullptr ? info->print_hist() : false );
            uint        iter       = 0;

            //
            // we assume, that preconditioner is based on H-LU and start
            // with some steps of the preconditioned Richardson, i.e., H-LU
            // iteration
            //

            TStopCriterion  rich_stop( _stop_criterion );

            rich_stop.max_iter = 10;
            
            TRichardson     rich_solver( 1.0, rich_stop );
            TInfo           rich_info( store_hist, print_hist );

            rich_solver.initialise_start_value( _initialise_start_value );

            auto  rich_tic = Time::Wall::now();
            
            rich_solver.solve( A, x, b, W, & rich_info );

            auto  rich_toc = Time::Wall::since( rich_tic );

            std::cout << "t(Rich) = " << rich_toc << std::endl;
            
            // copy iteration data
            if ( info != nullptr )
            {
                info->set_print_hist( false );
                
                for ( auto  i : rich_info.history() )
                    info->append( iter++, i.res_norm );

                info->set_print_hist( print_hist );
                info->set_status( rich_info.status() );
            }// if
                      
            if ( rich_info.has_converged() )
                return;

            //
            // adjust stopping criterion based on current residual
            //

            TStopCriterion  tfqmr_stop( _stop_criterion );

            {
                auto  res_norm_start = rich_info.history().begin()->res_norm;
                auto  res_norm_end   = rich_info.res_norm();

                if ( res_norm_end != real(0) )
                    tfqmr_stop.rel_res_reduct = rich_stop.rel_res_reduct * ( res_norm_start / res_norm_end );

                tfqmr_stop.max_iter = rich_stop.max_iter - rich_info.n_iter();
            }

            //
            // Richardson did not converge, so try TFQMR
            //

            TTFQMR          tfqmr_solver( tfqmr_stop );
            TInfo           tfqmr_info( store_hist, print_hist );
            bool            has_breakdown = false;

            tfqmr_solver.initialise_start_value( false );

            auto  tfqmr_tic = Time::Wall::now();
            
            tfqmr_solver.solve( A, x, b, W, & tfqmr_info );
            
            auto  tfqmr_toc = Time::Wall::since( tfqmr_tic );
            
            std::cout << "t(TFQMR) = " << tfqmr_toc << std::endl;

            // copy iteration history
            if ( info != nullptr )
            {
                info->set_print_hist( false );

                for ( auto  i : tfqmr_info.history() )
                    info->append( iter++, i.res_norm );

                info->set_print_hist( print_hist );
                info->set_status( tfqmr_info.status() );
            }// if
                      
            if ( tfqmr_info.has_converged() )
                return;

            //
            // in case of breakdown in TFQMR, try GMRes
            //

            if ( tfqmr_info.has_failed() )
            {
                TStopCriterion  gmres_stop( _stop_criterion );

                auto  res_norm_start = tfqmr_info.history().begin()->res_norm;
                auto  res_norm_end   = tfqmr_info.res_norm();

                if ( res_norm_end != real(0) )
                    gmres_stop.rel_res_reduct = tfqmr_stop.rel_res_reduct * ( res_norm_start / res_norm_end );
                
                gmres_stop.max_iter = tfqmr_stop.max_iter - tfqmr_info.n_iter();

                
                TGMRES          gmres_solver( gmres_stop );
                TInfo           gmres_info( store_hist, print_hist );

                gmres_solver.initialise_start_value( false );

                auto  gmres_tic = Time::Wall::now();

                gmres_solver.solve( A, x, b, W, & gmres_info );

                auto  gmres_toc = Time::Wall::since( gmres_tic );

                std::cout << "t(GMRES) = " << gmres_toc << std::endl;

                // copy iteration history
                if ( info != nullptr )
                {
                    info->set_print_hist( false );

                    for ( auto  i : gmres_info.history() )
                        info->append( iter++, i.res_norm );

                    info->set_print_hist( print_hist );
                    info->set_status( gmres_info.status() );
                }// if
            }// if

            #endif
        }// else
    }// else
}

}// namespace HLIB
