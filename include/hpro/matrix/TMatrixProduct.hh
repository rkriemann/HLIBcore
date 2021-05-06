#ifndef __HLIB_TMATRIXPRODUCT_HH
#define __HLIB_TMATRIXPRODUCT_HH
//
// Project     : HLib
// File        : TMatrixProduct.hh
// Description : Represents product of two matrices (linear ops)
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <deque>

#include "hpro/matrix/TLinearOperator.hh"

namespace HLIB
{

// local RTTI types
DECLARE_TYPE( TMatrixProduct );

//!
//! \ingroup Matrix_Module
//! \class   TMatrixProduct
//! \brief   Represents product α₁A₁ · α₂A₂ · α₃A₃... of matrices (linear ops)
//!
template < typename T_value >
class TMatrixProduct : public TLinearOperator
{
    //! @cond
public:

    using  value_t = T_value;
    
private:

    struct factor_t {
        const TLinearOperator *  linop;
        const matop_t            op;
        const value_t            scale;
    };
    
    // scaling factors (Π_i α_i)
    std::deque< factor_t >  _factors;

    // if true, local object has ownership
    const bool              _is_owner;
    
    //! @endcond
public:
    ///////////////////////////////////////////////////////////
    //
    // ctor/dtor
    //

    //! single matrix
    TMatrixProduct ( const value_t            alpha1,
                     const TLinearOperator *  A1,
                     const bool               is_owner = false );

    TMatrixProduct ( const value_t            alpha1,
                     const matop_t            op1,
                     const TLinearOperator *  A1,
                     const bool               is_owner = false );

    //! two matrices
    TMatrixProduct ( const value_t            alpha1,
                     const TLinearOperator *  A1,
                     const value_t            alpha2,
                     const TLinearOperator *  A2,
                     const bool               is_owner = false );

    TMatrixProduct ( const value_t            alpha1,
                     const matop_t            op1,
                     const TLinearOperator *  A1,
                     const value_t            alpha2,
                     const matop_t            op2,
                     const TLinearOperator *  A2,
                     const bool               is_owner = false );

    //! three matrices
    TMatrixProduct ( const value_t            alpha1,
                     const TLinearOperator *  A1,
                     const value_t            alpha2,
                     const TLinearOperator *  A2,
                     const value_t            alpha3,
                     const TLinearOperator *  A3,
                     const bool               is_owner = false );

    TMatrixProduct ( const value_t            alpha1,
                     const matop_t            op1,
                     const TLinearOperator *  A1,
                     const value_t            alpha2,
                     const matop_t            op2,
                     const TLinearOperator *  A2,
                     const value_t            alpha3,
                     const matop_t            op3,
                     const TLinearOperator *  A3,
                     const bool               is_owner = false );

    virtual ~TMatrixProduct ();

    ///////////////////////////////////////////////////////////
    //
    // linear operator properties
    //

    //! return true, if field type is complex
    virtual bool  is_complex      () const;
    
    //! return true, of operator is self adjoint
    virtual bool  is_self_adjoint () const;
    
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
                                const matop_t    op = apply_normal ) const;

    //!
    //! mapping function with update: \f$ y := y + \alpha A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply_add   ( const real       alpha,
                                const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;
    virtual void  capply_add  ( const complex    alpha,
                                const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;

    virtual void  apply_add   ( const real       alpha,
                                const TMatrix *  X,
                                TMatrix *        Y,
                                const matop_t    op = apply_normal ) const;
    
    //! same as above but only the dimension of the vector spaces is tested,
    //! not the corresponding index sets
    virtual void  apply_add   ( const real                       alpha,
                                const BLAS::Vector< real > &     x,
                                BLAS::Vector< real > &           y,
                                const matop_t                    op = apply_normal ) const;
    virtual void  apply_add   ( const complex                    alpha,
                                const BLAS::Vector< complex > &  x,
                                BLAS::Vector< complex > &        y,
                                const matop_t                    op = apply_normal ) const;

    virtual void  apply_add   ( const real                       alpha,
                                const BLAS::Matrix< real > &     X,
                                BLAS::Matrix< real > &           Y,
                                const matop_t                    op = apply_normal ) const;
    virtual void  apply_add   ( const complex                    alpha,
                                const BLAS::Matrix< complex > &  X,
                                BLAS::Matrix< complex > &        Y,
                                const matop_t                    op = apply_normal ) const;

    ///////////////////////////////////////////////////////////
    //
    // access to vector space elements
    //

    //! return dimension of domain
    virtual size_t  domain_dim () const;
    
    //! return dimension of range
    virtual size_t  range_dim  () const;
    
    //! return vector in domain space
    virtual auto    domain_vector  () const -> std::unique_ptr< TVector >;

    //! return vector in range space
    virtual auto    range_vector   () const -> std::unique_ptr< TVector >;

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( TMatrixProduct, TLinearOperator );
};

//
// functions to return matrix product objects
//

inline
std::unique_ptr< TMatrixProduct< real > >
matrix_product ( const TLinearOperator *  A0,
                 const bool               is_owner = false )
{
    return std::make_unique< TMatrixProduct< real > >( real(1), A0, is_owner );
}
                 
inline                 
std::unique_ptr< TMatrixProduct< real > >
matrix_product ( const real               alpha0,
                 const TLinearOperator *  A0,
                 const bool               is_owner = false )
{
    return std::make_unique< TMatrixProduct< real > >( alpha0, A0, is_owner );
}

inline                 
std::unique_ptr< TMatrixProduct< complex > >
matrix_product ( const complex            alpha0,
                 const TLinearOperator *  A0,
                 const bool               is_owner = false )
{
    return std::make_unique< TMatrixProduct< complex > >( alpha0, A0, is_owner );
}
                 
inline
std::unique_ptr< TMatrixProduct< real > >
matrix_product ( const TLinearOperator *  A0,
                 const TLinearOperator *  A1,
                 const bool               is_owner = false )
{
    return std::make_unique< TMatrixProduct< real > >( real(1), A0, real(1), A1, is_owner );
}
                 
template < typename value_t >
std::unique_ptr< TMatrixProduct< value_t > >
matrix_product ( const value_t            alpha0,
                 const TLinearOperator *  A0,
                 const value_t            alpha1,
                 const TLinearOperator *  A1,
                 const bool               is_owner = false )
{
    return std::make_unique< TMatrixProduct< value_t > >( alpha0, A0, alpha1, A1, is_owner );
}
                 
template < typename value_t >
std::unique_ptr< TMatrixProduct< value_t > >
matrix_product ( const value_t            alpha0,
                 const matop_t            op0,
                 const TLinearOperator *  A0,
                 const value_t            alpha1,
                 const matop_t            op1,
                 const TLinearOperator *  A1,
                 const bool               is_owner = false )
{
    return std::make_unique< TMatrixProduct< value_t > >( alpha0, op0, A0, alpha1, op1, A1, is_owner );
}
                 
inline
std::unique_ptr< TMatrixProduct< real > >
matrix_product ( const TLinearOperator *  A0,
                 const TLinearOperator *  A1,
                 const TLinearOperator *  A2,
                 const bool               is_owner = false )
{
    return std::make_unique< TMatrixProduct< real > >( real(1), A0, real(1), A1, real(1), A2, is_owner );
}
                 
template < typename value_t >
std::unique_ptr< TMatrixProduct< value_t > >
matrix_product ( const value_t            alpha0,
                 const TLinearOperator *  A0,
                 const value_t            alpha1,
                 const TLinearOperator *  A1,
                 const value_t            alpha2,
                 const TLinearOperator *  A2,
                 const bool               is_owner = false )
{
    return std::make_unique< TMatrixProduct< value_t > >( alpha0, A0, alpha1, A1, alpha2, A2, is_owner );
}
                 
template < typename value_t >
std::unique_ptr< TMatrixProduct< value_t > >
matrix_product ( const value_t            alpha0,
                 const matop_t            op0,
                 const TLinearOperator *  A0,
                 const value_t            alpha1,
                 const matop_t            op1,
                 const TLinearOperator *  A1,
                 const value_t            alpha2,
                 const matop_t            op2,
                 const TLinearOperator *  A2,
                 const bool               is_owner = false )
{
    return std::make_unique< TMatrixProduct< value_t > >( alpha0, op0, A0, alpha1, op1, A1, alpha2, op2, A2, is_owner );
}
                 
}// namespace HLIB

#endif  // __HLIB_TMATRIXPRODUCT_HH
