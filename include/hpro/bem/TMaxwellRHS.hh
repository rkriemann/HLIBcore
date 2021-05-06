#ifndef __HLIB_TMAXWELLRHS_HH
#define __HLIB_TMAXWELLRHS_HH
//
// Project     : HLib
// File        : TMaxwellRHS.hh
// Description : classes for building RHS-vectors in Maxwell applications
// Author      : Ronald Kriemann, Jonas Ballani
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/bem/TBEMRHS.hh"

namespace HLIB
{
//
// value type of BEM function with traits
//
struct complex3_t
{
    complex  val[3];
};

template <>
struct real_type< complex3_t >
{
    using  type_t = real;
};

template <>
struct is_complex_type< complex3_t >
{
    static const bool value = true;
};

////////////////////////////////////////////////////////
//
// building RHS for Maxwell EFIE using quadrature
//

template < typename T_fnspace >
class TMaxwellEFIERHS : public TBEMRHS< T_fnspace, complex3_t >
{
public:
    //! template argument as internal type
    using  value_t   = complex3_t;
    using  fnspace_t = T_fnspace;

private:
    // quadrature order
    const uint  _quad_order;
    
public:
    ///////////////////////////////////////////////////////////////
    //
    // constructur and destructor
    //

    TMaxwellEFIERHS ( const uint  aquad_order );

    virtual ~TMaxwellEFIERHS () {}

    ///////////////////////////////////////////////////////////////
    //!
    //! build RHS-vector for \a fnspace and RHS
    //! function \a rhs; by \a perm one can set a permutation for
    //! to final vector, e.g., from internal to external numbering
    //!
    virtual std::unique_ptr< TVector >  build ( const fnspace_t *                fnspace,
                                                const TBEMFunction< value_t > *  rhs,
                                                const TPermutation *             perm = NULL ) const;

    DISABLE_COPY_OP( TMaxwellEFIERHS );
};

////////////////////////////////////////////////////////
//
// building RHS for Maxwell MFIE using quadrature
//

template < typename T_fnspace >
class TMaxwellMFIERHS : public TBEMRHS< T_fnspace, complex3_t >
{
public:
    //! template argument as internal type
    using  value_t   = complex3_t;
    using  fnspace_t = T_fnspace;
          
private:
    // quadrature order
    const uint  _quad_order;
    
public:
    ///////////////////////////////////////////////////////////////
    //
    // constructur and destructor
    //

    TMaxwellMFIERHS ( const uint  aquad_order );

    virtual ~TMaxwellMFIERHS () {}

    ///////////////////////////////////////////////////////////////
    //!
    //! build RHS-vector for \a fnspace and RHS
    //! function \a rhs; by \a perm one can set a permutation for
    //! to final vector, e.g., from internal to external numbering
    //!
    virtual std::unique_ptr< TVector >  build ( const fnspace_t *                fnspace,
                                                const TBEMFunction< value_t > *  rhs,
                                                const TPermutation *             perm = NULL ) const;

    DISABLE_COPY_OP( TMaxwellMFIERHS );
};

}// namespace

#endif  // __HLIB_TMAXWELLRHS_HH
