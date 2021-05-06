#ifndef __HLIB_TSOLVER_HH
#define __HLIB_TSOLVER_HH
//
// Project     : HLib
// File        : TSolver.hh
// Description : base-class for all iterative solvers
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <ostream>
#include <list>

#include "hpro/matrix/TLinearOperator.hh"
#include "hpro/vector/TVector.hh"

namespace HLIB
{

enum solver_status_t
{
    iterating = 0,  // solver is still iterating
    converged,      // solver has converged
    diverged,       // solver has diverged
    failed          // solver has failed
};

//!
//! \ingroup  Solver_Module
//! \struct   TSolverInfo
//! \brief    datatype to share iteration history
//!
struct TSolverInfo
{
    //! @cond
        
public:
    //!
    //! \struct  hist_t
    //! \brief   used to store history information
    //!
    struct hist_t {
        real   res_norm;
    };
        
protected:
    //! name of solver
    std::string          _solver;

    //! status of solver
    solver_status_t      _status;
        
    //! number of iterations
    size_t               _n_iter;
        
    //! average convergance rate
    real                 _conv_rate;

    //! residual norm after iteration
    real                 _res_norm;

    //! if true, iteration history information is stored
    bool                 _store_hist;

    //! iteration history data
    std::list< hist_t >  _history;

    //! if true, print convergence history
    bool                 _print;

    //! internal flags
    const bool           _flag0, _flag1;
    
    //! @endcond
        
public:
    //
    // ctor
    //
        
    //! constructor and store history if \a store_hist is true
    //! and print history if \a print is true respectively
    TSolverInfo ( const bool store_hist = false,
                  const bool print      = false );

    //
    // access data
    //
        
    //! return true if convergence information is present
    bool has_data () const { return (n_iter() != 0) || (history().size() > 0); }

    //! return solver name
    const std::string &  solver () const { return _solver; }

    //! return true, if iteration has converged
    bool           has_converged () const { return _status == converged; }

    //! return true, if iteration has diverged
    bool           has_diverged  () const { return _status == diverged; }

    //! return true, if iteration has failed
    bool           has_failed    () const { return _status == failed; }

    //! return status of solver
    solver_status_t status        () const { return _status; }

    //! return number of iteration steps
    size_t         n_iter    () const { return _n_iter; }

    //! return convergence rate
    real           conv_rate () const { return _conv_rate; }

    //! return current residual norm
    real           res_norm  () const { return _res_norm; }

    //! return iteration history data
    const std::list< hist_t > &  history () const { return _history; }

    //! return true, if history is stored
    bool           store_hist () const { return _store_hist; }
        
    //! return true, if convergence history is printed is stored
    bool           print_hist () const { return _print; }
        
    //
    // write internal data and options
    //
        
    //! append data for single iteration step to history (if set so)
    //! and update internal iteration data
    void append         ( const uint it, const real res_norm );

    //! set solver name
    void set_solver     ( const std::string & name ) { _solver = name; }

    //! set convergence status
    void set_status     ( const solver_status_t  s  ) { _status = s; }

    //! set number of iteration steps
    void set_n_iter     ( const uint     n    ) { _n_iter    = n; }

    //! set convergence rate
    void set_conv_rate  ( const real     conv ) { _conv_rate = conv; }

    //! set current residual norm
    void set_res_norm   ( const real     norm ) { _res_norm  = norm; }

    //! turn on/off printing of history
    void set_print_hist ( const bool     b    ) { _print = b; }

    //! reset all data
    void reset ();

    //
    // misc.
    //
        
    //! convert data in information object to string
    std::string  to_string () const;

    //! stream output
    friend std::ostream & operator << ( std::ostream & os, const TSolverInfo & info )
    {
        return os << info.to_string();
    }

    //! write history in Gnuplot format to stream
    void  print_gnuplot ( std::ostream &  os );
};

//!
//! \ingroup  Solver_Module
//! \struct   TStopCriterion
//! \brief     stopping criterion for iterative solvers
//!
class TStopCriterion
{
public:
    // maximal number of iterations
    uint  max_iter;
        
    // absolute reduction of residual norm
    real  abs_res_reduct;
        
    // relative reduction of residual norm compared to start residual
    real  rel_res_reduct;
        
    // maximal relative growth of residual norm compared to start residual
    real  rel_res_growth;

    // ctors
    TStopCriterion ( const uint  max_iter       = CFG::Solver::max_iter,
                     const real  abs_res_red    = CFG::Solver::abs_res_red,
                     const real  rel_res_red    = CFG::Solver::rel_res_red,
                     const real  rel_res_growth = CFG::Solver::rel_res_growth );

    TStopCriterion ( const TStopCriterion &  stop_crit );
    TStopCriterion ( TStopCriterion &&       stop_crit );

    TStopCriterion & operator =  ( const TStopCriterion &  stop_crit );
    TStopCriterion & operator =  ( TStopCriterion &&       stop_crit );
        
    //! return true if stop condition is met
    virtual bool stopped ( const uint     it, 
                           const real     norm, 
                           const real     norm0, 
                           TSolverInfo *  info ) const;

    //! convert to string
    std::string  to_string () const;
};

//!
//! \ingroup  Solver_Module
//! \brief    joins stop criteria while overloading/preferring undefined values
//!           (interprete as logical and)
//!
TStopCriterion
operator + ( const TStopCriterion &  crit1,
             const TStopCriterion &  crit2 );

//!
//! \ingroup  Solver_Module
//! \brief    sets maximal number of iteration steps
//!
TStopCriterion 
max_steps           ( const uint  steps );

//!
//! \ingroup  Solver_Module
//! \brief    sets relative reduction of residual
//!
TStopCriterion 
relative_reduction  ( const real  red );

//!
//! \ingroup  Solver_Module
//! \brief    sets absolute reduction of residual
//!
TStopCriterion 
absolute_reduction  ( const real  red );

//!
//! \ingroup  Solver_Module
//! \class    TSolver
//! \brief    Solver base class defining interface and implementing
//!           simple solver (Richardson iteration)
//!
class TSolver
{
protected:

    //! @cond

    // stop criterion
    TStopCriterion  _stop_criterion;

    // if true, start value of iteration will be initialised by iteration method
    bool            _initialise_start_value;

    // if true, exact residual is computed (default: depends on precond./solver)
    bool            _use_exact_residual;
    
    //! @admcond
    
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct solver object with corresponding stop criteria
    TSolver ( const TStopCriterion &  stop_crit );

    //! dtor
    virtual ~TSolver ();

    //! turn initialisation of start value on/off
    void  initialise_start_value ( const bool  b );

    //! choose type of residual
    void  use_exact_residual     ( const bool  b );

    //! return choice of residual
    bool  use_exact_residual     () const { return _use_exact_residual; }
    
    ////////////////////////////////////////////////
    //
    // solving the system
    //

    //! solve A·x = b with optional preconditioner \a W
    virtual void solve         ( const TLinearOperator *  A,
                                 TVector *                x,
                                 const TVector *          b,
                                 const TLinearOperator *  W    = nullptr,
                                 TSolverInfo *            info = nullptr ) const = 0;

    //! set stop criterion
    virtual void set_stop_crit ( const TStopCriterion &  stop_crit );

    //! return true if stop condition is met
    virtual bool stopped       ( const uint     it, 
                                 const real     norm, 
                                 const real     norm0, 
                                 TSolverInfo *  info ) const;
    
    //! initialises start value of iteration
    virtual void set_start_value ( TVector *                x,
                                   const TVector *          b,
                                   const TLinearOperator *  W = nullptr ) const;

};

}// namespace HLIB

#endif  // __HLIB_TSOLVER_HH
