#ifndef __HLIB_TBEMRHS_HH
#define __HLIB_TBEMRHS_HH
//
// Project     : HLib
// File        : TBEMRHS.hh
// Description : classes for building RHS-vectors in BEM-applications
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/cluster/TIndexSet.hh"
#include "hpro/cluster/TPermutation.hh"

#include "hpro/bem/TFnSpace.hh"

#include "hpro/vector/TVector.hh"

namespace HLIB
{

////////////////////////////////////////////////////////
//
// baseclass defining interface
//
template < typename  T_fnspace,
           typename  T_val >
class TBEMRHS
{
public:
    //! template arguments as internal types
    using  fnspace_t = T_fnspace;
    using  value_t   = T_val;
    
public:
    ///////////////////////////////////////////////////////////////
    //
    // constructur and destructor
    //

    TBEMRHS () {}

    virtual ~TBEMRHS () {}

    ///////////////////////////////////////////////////////////////
    //!
    //! build RHS-vector for \a fnspace and RHS
    //! function \a rhs; by \a perm one can set a permutation for
    //! to final vector, e.g., from internal to external numbering
    //!
    virtual std::unique_ptr< TVector >  build ( const fnspace_t *                fnspace,
                                                const TBEMFunction< value_t > *  rhs,
                                                const TPermutation *             perm = NULL ) const = 0;
};

////////////////////////////////////////////////////////
//
// building RHS using quadrature
//
template < typename  T_fnspace,
           typename  T_val >
class TQuadBEMRHS : public TBEMRHS< T_fnspace, T_val >
{
public:
    //! template arguments as internal types
    using  fnspace_t = T_fnspace;
    using  value_t   = T_val;
    
private:
    // quadrature order
    const uint  _quad_order;
    
public:
    ///////////////////////////////////////////////////////////////
    //
    // constructur and destructor
    //

    TQuadBEMRHS ( const uint  quad_order );

    virtual ~TQuadBEMRHS () {}

    ///////////////////////////////////////////////////////////////
    //!
    //! build RHS-vector for \a fnspace and RHS
    //! function \a rhs; by \a perm one can set a permutation for
    //! to final vector, e.g., from internal to external numbering
    //!
    virtual std::unique_ptr< TVector >  build ( const fnspace_t *                fnspace,
                                                const TBEMFunction< value_t > *  rhs,
                                                const TPermutation *             perm = NULL ) const;

    DISABLE_COPY_OP( TQuadBEMRHS );
};

}// namespace

#endif  // __HLIB_TBEMRHS_HH
