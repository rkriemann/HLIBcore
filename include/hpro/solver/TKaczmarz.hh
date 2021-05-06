#ifndef __HLIB_TKACZMARZ_HH
#define __HLIB_TKACZMARZ_HH
//
// Project     : HLib
// File        : TKaczmarz.hh
// Description : iterative solver based on Kaczmarz algorithm
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TSolver.hh"

namespace HLIB
{

//! @cond

//!
//! \ingroup  Solver_Module
//! \class    TKaczmarz
//! \brief    Implements Kaczmarz iteration.
//!
class TKaczmarz : public TSolver
{
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct Kaczmarz solver object with corresponding stop criteria
    TKaczmarz ( const TStopCriterion &  stop_crit = TStopCriterion() );

    virtual ~TKaczmarz ();

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
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
inline
void
kaczmarz ( const TLinearOperator *  A,
           TVector *                x,
           const TVector *          b,
           const TLinearOperator *  W         = nullptr,
           TSolverInfo *            info      = nullptr,
           const TStopCriterion &   stop_crit = TStopCriterion() )
{
    TKaczmarz  solver( stop_crit );

    solver.solve( A, x, b, W, info );
}

//! @endcond

}// namespace HLIB

#endif  // __HLIB_TKACZMARZ_HH
