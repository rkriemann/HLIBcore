#ifndef __HLIB_TLINEAROPERATOR_HH
#define __HLIB_TLINEAROPERATOR_HH
//
// Project     : HLib
// File        : TLinearOperator.hh
// Description : baseclass for all linear operators
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/base/error.hh"
#include "hpro/base/TTypeInfo.hh"
#include "hpro/blas/types.hh"
#include "hpro/blas/Vector.hh"
#include "hpro/blas/Matrix.hh"
#include "hpro/vector/TVector.hh"

namespace HLIB
{

// local RTTI types
DECLARE_TYPE( TLinearOperator );

class TMatrix;

//!
//! \ingroup Matrix_Module
//! \class   TLinearOperator
//! \brief   Base class for all linear operators mapping vectors to vectors.
//!
//!          Many standard arithmetic operations only depend upon a linear operator
//!          providing the mapping between vectors, e.g. iterativ solvers. Instead of
//!          requiring a full matrix and hence the need for an implementation of the
//!          full matrix algebra, an object of type TLinearOperator is fully sufficient
//!          in such cases.
//!
class TLinearOperator : public TTypeInfo
{
public:
    ///////////////////////////////////////////////////////////
    //
    // dtor
    //

    virtual ~TLinearOperator () {}

    ///////////////////////////////////////////////////////////
    //
    // linear operator properties
    //

    //! return true, if field type is complex
    virtual bool  is_complex      () const = 0;
    
    //! return true, of operator is self adjoint
    virtual bool  is_self_adjoint () const = 0;
    
    ///////////////////////////////////////////////////////////
    //
    // linear operator mapping
    //
    
    //!
    //! mapping function of linear operator \f$A\f$, e.g. \f$ y := A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply       ( const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const = 0;

    //!
    //! mapping function with update: \f$ y := y + \alpha A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply_add   ( const real       alpha,
                                const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const = 0;
    virtual void  capply_add  ( const complex    alpha,
                                const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const = 0;

    virtual void  apply_add   ( const real       alpha,
                                const TMatrix *  X,
                                TMatrix *        Y,
                                const matop_t    op = apply_normal ) const = 0;
    
    //! same as above but only the dimension of the vector spaces is tested,
    //! not the corresponding index sets
    virtual void  apply_add   ( const real                       alpha,
                                const BLAS::Vector< real > &     x,
                                BLAS::Vector< real > &           y,
                                const matop_t                    op = apply_normal ) const = 0;
    virtual void  apply_add   ( const complex                    alpha,
                                const BLAS::Vector< complex > &  x,
                                BLAS::Vector< complex > &        y,
                                const matop_t                    op = apply_normal ) const = 0;

    virtual void  apply_add   ( const real                       alpha,
                                const BLAS::Matrix< real > &     x,
                                BLAS::Matrix< real > &           y,
                                const matop_t                    op = apply_normal ) const = 0;
    virtual void  apply_add   ( const complex                    alpha,
                                const BLAS::Matrix< complex > &  x,
                                BLAS::Matrix< complex > &        y,
                                const matop_t                    op = apply_normal ) const = 0;

    ///////////////////////////////////////////////////////////
    //
    // access to vector space data
    //

    //! return dimension of domain
    virtual size_t  domain_dim () const = 0;
    
    //! return dimension of range
    virtual size_t  range_dim  () const = 0;
    
    //! return vector in domain space
    virtual auto    domain_vector  () const -> std::unique_ptr< TVector > = 0;

    //! return vector in range space
    virtual auto    range_vector   () const -> std::unique_ptr< TVector > = 0;

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //

    HLIB_RTTI_BASE( TLinearOperator );
};

//////////////////////////////////////////////////////////
//
// functional form of TLinearOperator methods
//

//!
//! \ingroup Matrix_Module
//! \fn      apply
//! \brief   Functional form of TLinearOperator::apply.
//!
inline
void
apply ( const TLinearOperator *  A,
        const TVector *          x,
        TVector *                y,
        const matop_t            op = apply_normal )
{
    if (( A == NULL ) || ( x == NULL ) || ( y == NULL ))
        HERROR( ERR_ARG, "apply", "argument is NULL" );

    A->apply( x, y, op );
}

//!
//! \ingroup Matrix_Module
//! \fn      apply
//! \brief   Functional form of TLinearOperator::apply.
//!
template < typename T_field >
void
apply_add ( const T_field            alpha,
            const TLinearOperator *  A,
            const TVector *          x,
            TVector *                y,
            const matop_t            op = apply_normal );

template <>
inline
void
apply_add< real > ( const real               alpha,
                    const TLinearOperator *  A,
                    const TVector *          x,
                    TVector *                y,
                    const matop_t            op )
{
    if (( A == NULL ) || ( x == NULL ) || ( y == NULL ))
        HERROR( ERR_ARG, "apply_add", "argument is NULL" );

    A->apply_add( alpha, x, y, op );
}

template <>
inline
void
apply_add< complex > ( const complex            alpha,
                       const TLinearOperator *  A,
                       const TVector *          x,
                       TVector *                y,
                       const matop_t            op )
{
    if (( A == NULL ) || ( x == NULL ) || ( y == NULL ))
        HERROR( ERR_ARG, "apply_add", "argument is NULL" );

    A->capply_add( alpha, x, y, op );
}

inline
void
apply_add ( const real               alpha,
            const TLinearOperator *  A,
            const TMatrix *          X,
            TMatrix *                Y,
            const matop_t            op = apply_normal )
{
    if (( A == NULL ) || ( X == NULL ) || ( Y == NULL ))
        HERROR( ERR_ARG, "apply_add", "argument is NULL" );

    A->apply_add( alpha, X, Y, op );
}

//! return dimension of domain
inline size_t  domain_dim ( const TLinearOperator *  op ) { return op->domain_dim(); }
inline size_t  domain_dim ( const TLinearOperator &  op ) { return op.domain_dim(); }
    
//! return dimension of range
inline size_t  range_dim  ( const TLinearOperator *  op ) { return op->range_dim(); }
inline size_t  range_dim  ( const TLinearOperator &  op ) { return op.range_dim(); }
    
//! return vector in domain space
inline std::unique_ptr< TVector >  domain_vector  ( const TLinearOperator *  op ) { return op->domain_vector(); }
inline std::unique_ptr< TVector >  domain_vector  ( const TLinearOperator &  op ) { return op.domain_vector(); }

//! return vector in range space
inline std::unique_ptr< TVector >  range_vector   ( const TLinearOperator *  op ) { return op->range_vector(); }
inline std::unique_ptr< TVector >  range_vector   ( const TLinearOperator &  op ) { return op.range_vector(); }

//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//!
//! write linear operator to file
//!
void
write ( const TLinearOperator *  M,
        const std::string &      filename,
        const std::string &      matname );

}// namespace DBG

}// namespace HLIB

#endif  // __HLIB_TLINEAROPERATOR_HH
