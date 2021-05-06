#ifndef __HLIB_TCGS_HH
#define __HLIB_TCGS_HH
//
// Project     : HLib
// File        : TCGS.hh
// Description : class implementing CGS
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TSolver.hh"

namespace HLIB
{

//!
//! \ingroup  Solver_Module
//! \class    TCGS
//! \brief    Implements squared conjugate gradient iteration.
//!
class TCGS : public TSolver
{
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct CG solver object with corresponding stop criteria
    TCGS ( const TStopCriterion &  stop_crit = TStopCriterion() );
    
    virtual ~TCGS ();

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
//! \fn       cgs
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
inline
void
cgs ( const TLinearOperator *  A,
      TVector *                x,
      const TVector *          b,
      const TLinearOperator *  W         = nullptr,
      TSolverInfo *            info      = nullptr,
      const TStopCriterion &   stop_crit = TStopCriterion() )
{
    TCGS  solver( stop_crit );

    solver.solve( A, x, b, W, info );
}

}// namespace HLIB

#endif  // __HLIB_TCGS_HH
