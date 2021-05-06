#ifndef __HLIB_TAUTOSOLVER_HH
#define __HLIB_TAUTOSOLVER_HH
//
// Project     : HLib
// File        : TAutoSolver.hh
// Description : class implementing a solver which automatically decides best strategy
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TSolver.hh"

namespace HLIB
{

//!
//! \ingroup  Solver_Module
//! \class    TAutoSolver
//! \brief    Implements an iterative solver automatically choosing
//!           appropriate algorithm based on matrix criteria.
//!
class TAutoSolver : public TSolver
{
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct auto solver object with corresponding stop criteria
    TAutoSolver ( const TStopCriterion &  stop_crit = TStopCriterion() );
    
    virtual ~TAutoSolver ();

    ////////////////////////////////////////////////
    //
    // solving the system
    //

    //! solve A·x = b with optional preconditioner \a W
    virtual void solve ( const TLinearOperator *  A,
                         TVector *                x,
                         const TVector *          b,
                         const TLinearOperator *  W    = nullptr,
                         TSolverInfo *            info = nullptr ) const;
};

//!
//! \ingroup  Solver_Module
//! \brief    Solve A·x = b with optional preconditioner \a W (functional version).
//!
inline
void
solve ( const TLinearOperator *  A,
        TVector *                x,
        const TVector *          b,
        const TLinearOperator *  W,
        TSolverInfo *            info      = nullptr,
        const TStopCriterion &   stop_crit = TStopCriterion() )
{
    TAutoSolver  solver( stop_crit );

    solver.solve( A, x, b, W, info );
}

//!
//! \ingroup  Solver_Module
//! \brief    Solve A·x = b (functional version).
//!
inline
void
solve ( const TLinearOperator *  A,
        TVector *                x,
        const TVector *          b,
        TSolverInfo *            info      = nullptr,
        const TStopCriterion &   stop_crit = TStopCriterion() )
{
    solve( A, x, b, nullptr, info, stop_crit );
}

}// namespace HLIB

#endif  // __HLIB_TAUTOSOLVER_HH
