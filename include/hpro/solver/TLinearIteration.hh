#ifndef __HLIB_TLINEARITERATION_HH
#define __HLIB_TLINEARITERATION_HH
//
// Project     : HLib
// File        : TLinearIteration.hh
// Description : Linear Iteration solver
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <ostream>
#include <list>

#include "hpro/solver/TSolver.hh"

namespace HLIB
{

//!
//! \ingroup  Solver_Module
//! \class    TLinearIteration
//! \brief    Implements linear iteration \f$x_{i+1} = x_k + N (A x_i - b)\f$
//!
//! \detail   The class TLinearIteration implements a linear iteration solver
//!           based on the second normal form \f$x_{i+1} = x_k + N (A x_i - b)\f$
//!           with the iteration matrix \f$N\f$. Please note, that the iteration
//!           matrix is defined as the preconditioner \f$W\f$, e.g., for other
//!           linear iterations like Jacobi, GS or SOR, please use the correct
//!           preconditioner.
//!
class TLinearIteration : public TSolver
{
protected:
    //! @cond

    // damping factor
    real  _damping;

    //! @endcond
    
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct linear iteration solver object 
    TLinearIteration ( const TStopCriterion &  stop_crit = TStopCriterion() );
    
    //! construct linear iteration solver object with damping
    TLinearIteration ( const real              damping,
                       const TStopCriterion &  stop_crit = TStopCriterion() );
    
    //! dtor
    virtual ~TLinearIteration ();

    ////////////////////////////////////////////////
    //
    // solving the system
    //

    //! solve A·x = b with optional preconditioner \a W
    virtual void solve  ( const TLinearOperator *  A,
                          TVector *                x,
                          const TVector *          b,
                          const TLinearOperator *  W    = nullptr,
                          TSolverInfo *            info = nullptr ) const;

    //! solve A·X = B with optional preconditioner \a W
    virtual void solve  ( const TLinearOperator *  A,
                          TMatrix *                X,
                          const TMatrix *          B,
                          const TLinearOperator *  W    = nullptr,
                          TSolverInfo *            info = nullptr ) const;
};

//!
//! \ingroup  Solver_Module
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
inline
void linear_iteration ( const TLinearOperator *  A,
                        TVector *                x,
                        const TVector *          b,
                        const TLinearOperator *  W         = nullptr,
                        TSolverInfo *            info      = nullptr,
                        const TStopCriterion &   stop_crit = TStopCriterion() )
{
    TLinearIteration  solver( 1.0, stop_crit );

    solver.solve( A, x, b, W, info );
}

/////////////////////////////////////////////////////
//
// for backward compatibility
//

using TRichardson = TLinearIteration;

//!
//! \ingroup  Solver_Module
//! \brief    Solve A·x = b with optional preconditioner \a W (functional approach)
//!
inline
void richardson ( const TLinearOperator *  A,
                  TVector *                x,
                  const TVector *          b,
                  const TLinearOperator *  W         = nullptr,
                  TSolverInfo *            info      = nullptr,
                  const TStopCriterion &   stop_crit = TStopCriterion() )
{
    TLinearIteration  solver( 1.0, stop_crit );

    solver.solve( A, x, b, W, info );
}

}// namespace HLIB

#endif  // __HLIB_TRICHARDSON_HH
