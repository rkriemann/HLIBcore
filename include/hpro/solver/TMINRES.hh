#ifndef __HLIB_TMINRES_HH
#define __HLIB_TMINRES_HH
//
// Project     : HLib
// File        : TMINRES.hh
// Description : class implementing MINRES
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TSolver.hh"

namespace HLIB
{

//!
//! \ingroup  Solver_Module
//! \class    TMINRES
//! \brief    Implements the MINRES iteration
//!
class TMINRES : public TSolver
{
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct MINRES solver object with corresponding stop criteria
    TMINRES ( const TStopCriterion &  stop_crit = TStopCriterion() );
    
    virtual ~TMINRES ();

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
//! \fn       minres
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
inline
void
minres ( const TLinearOperator *  A,
         TVector *                x,
         const TVector *          b,
         const TLinearOperator *  W         = nullptr,
         TSolverInfo *            info      = nullptr,
         const TStopCriterion &   stop_crit = TStopCriterion() )
{
    TMINRES  solver( stop_crit );

    solver.solve( A, x, b, W, info );
}

}// namespace HLIB

#endif  // __HLIB_TMINRES_HH
