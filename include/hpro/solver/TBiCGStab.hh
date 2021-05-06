#ifndef __HLIB_TBICGSTAB_HH
#define __HLIB_TBICGSTAB_HH
//
// Project     : HLib
// File        : TBiCGStab.hh
// Description : class representing BiCG-Stab
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TSolver.hh"

namespace HLIB
{

//!
//! \ingroup  Solver_Module
//! \class    TBiCGStab
//! \brief    Implements BiCG-Stab iteration.
//!
class TBiCGStab : public TSolver
{
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct BiCG-Stab solver object with corresponding stop criteria
    TBiCGStab ( const TStopCriterion &  stop_crit = TStopCriterion() );

    virtual ~TBiCGStab ();

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
//! \fn       bicgstab
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
inline
void
bicgstab ( const TLinearOperator *  A,
           TVector *                x,
           const TVector *          b,
           const TLinearOperator *  W         = nullptr,
           TSolverInfo *            info      = nullptr,
           const TStopCriterion &   stop_crit = TStopCriterion() )
{
    TBiCGStab  solver( stop_crit );

    solver.solve( A, x, b, W, info );
}

}// namespace HLIB

#endif  // __HLIB_TBICGSTAB_HH
