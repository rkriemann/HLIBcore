#ifndef __HLIB_TGMRES_HH
#define __HLIB_TGMRES_HH
//
// Project     : HLib
// File        : TGMRES.hh
// Description : class implementing GMRES
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/solver/TSolver.hh"

namespace HLIB
{

//!
//! \ingroup  Solver_Module
//! \class    TGMRES
//! \brief    Implements GMRES iteration with restart.
//!
class TGMRES : public TSolver
{
private:
    //! maximal dimension of Krylov space before restart
    const uint  _restart;
    
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct GMRES solver object with default restart
    TGMRES ( const TStopCriterion &  stop_crit = TStopCriterion() );
    
    //! construct GMRES solver object with given restart
    TGMRES ( const uint              restart,
             const TStopCriterion &  stop_crit = TStopCriterion() );
    
    virtual ~TGMRES ();

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
//! \fn       gmres
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
inline
void
gmres ( const uint               restart,
        const TLinearOperator *  A,
        TVector *                x,
        const TVector *          b,
        const TLinearOperator *  W         = nullptr,
        TSolverInfo *            info      = nullptr,
        const TStopCriterion &   stop_crit = TStopCriterion() )
{
    TGMRES  solver( restart, stop_crit );

    solver.solve( A, x, b, W, info );
}

}// namespace HLIB

#endif  // __HLIB_TGMRES_HH
