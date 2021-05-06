#ifndef __HLIB_BEM_TBFCOEFFFN_H
#define __HLIB_BEM_TBFCOEFFFN_H
//
// Project     : HLib
// File        : TBFCoeffFn.hh
// Description : matrix coefficient functions for bilinear forms
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/matrix/TCoeffFn.hh"
#include "hpro/bem/TBEMBF.hh"

namespace HLIB
{

////////////////////////////////////////////////////////////////////
//!
//! \class  TBFCoeffFn
//! \brief  Provide matrix coefficients defined by bilinear forms
//!
template < typename  T_bf > 
class TBFCoeffFn : public TCoeffFn< typename T_bf::value_t >
{
public:
    //
    // template types as internal types
    //
    using  bf_t       = T_bf;
    using  ansatzsp_t = typename bf_t::ansatzsp_t;
    using  testsp_t   = typename bf_t::testsp_t;
    using  value_t    = typename bf_t::value_t;

protected:
    //! bilinear form to evaluate
    const bf_t *  _bf;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct coefficient function for bilinear form \a bf
    TBFCoeffFn ( const bf_t *  bf )
            : _bf( bf )
    {
        if ( _bf == NULL )
            HERROR( ERR_ARG, "(TBFCoeffFn)", "invalid bilinear form supplied (NULL)" );
    }

    //! destructor
    virtual ~TBFCoeffFn () {}

    //! return true if function is complex valued
    bool       is_complex     () const { return _bf->is_complex(); }
        
    //! return format of matrix, e.g. symmetric or hermitian
    matform_t  matrix_format  () const { return _bf->format(); }
    
    //////////////////////////////////////
    //
    // evaluation of coefficient
    //

    //! evaluate matrix coefficients in \a rowis × \a colis and store values
    //! in \a matrix (real valued)
    virtual void eval ( const std::vector< idx_t > &  rowidxs,
                        const std::vector< idx_t > &  colidxs,
                        value_t *                     matrix ) const
    {
        //
        // call bilinear form
        //
        
        const size_t  n = rowidxs.size();
        const size_t  m = colidxs.size();

        BLAS::Matrix< value_t >  M( n, m );
        
        _bf->eval( rowidxs, colidxs, M );
            
        for ( idx_t  i = 0; i < idx_t(n*m); ++i )
            matrix[i] = M.data()[i];
    }

    using TCoeffFn< value_t >::eval;

    DISABLE_COPY_OP( TBFCoeffFn );
};

}// namespace HLIB

#endif  // __HLIB_BEM_TBFCOEFFFN_H
