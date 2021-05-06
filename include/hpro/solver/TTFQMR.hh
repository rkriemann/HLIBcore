#ifndef __HLIB_TTFQMR_HH
#define __HLIB_TTFQMR_HH
//
// Project     : HLib
// File        : TTFQMR.hh
// Description : class implementing TFQMR
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TSolver.hh"

namespace HLIB
{

//!
//! \ingroup  Solver_Module
//! \class    TTFQMR
//! \brief    Implements squared conjugate gradient iteration.
//!
class TTFQMR : public TSolver
{
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct TFQMR solver object with corresponding stop criteria
    TTFQMR ( const TStopCriterion &  stop_crit = TStopCriterion() );
    
    virtual ~TTFQMR ();

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
//! \fn       tfqmr
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
inline
void
tfqmr ( const TLinearOperator *  A,
        TVector *                x,
        const TVector *          b,
        const TLinearOperator *  W         = nullptr,
        TSolverInfo *            info      = nullptr,
        const TStopCriterion &   stop_crit = TStopCriterion() )
{
    TTFQMR  solver( stop_crit );

    solver.solve( A, x, b, W, info );
}

}// namespace HLIB

#endif  // __HLIB_TTFQMR_HH
