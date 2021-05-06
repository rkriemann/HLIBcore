#ifndef __HLIB_TCG_HH
#define __HLIB_TCG_HH
//
// Project     : HLib
// File        : TCG.hh
// Description : class implementing CG
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TSolver.hh"

namespace HLIB
{

//!
//! \ingroup  Solver_Module
//! \class    TCG
//! \brief    Implements conjugate gradient iteration.
//!
class TCG : public TSolver
{
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct CG solver object
    TCG ( const TStopCriterion &  stop_crit = TStopCriterion() );
    
    virtual ~TCG ();

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
//! \fn       cg
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
inline
void
cg ( const TLinearOperator *  A,
     TVector *                x,
     const TVector *          b,
     const TLinearOperator *  W         = nullptr,
     TSolverInfo *            info      = nullptr,
     const TStopCriterion &   stop_crit = TStopCriterion() )
{
    TCG  solver( stop_crit );

    solver.solve( A, x, b, W, info );
}

}// namespace HLIB

#endif  // __HLIB_TCG_HH
