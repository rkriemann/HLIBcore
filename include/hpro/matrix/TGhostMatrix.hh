#ifndef __HLIB_TGHOSTMATRIX_HH
#define __HLIB_TGHOSTMATRIX_HH
//
// \file         TGhostMatrix.hh
//
// Project     : HLib
// Description : class for representing remote matrix blocks
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TGhostMatrix );

//!
//! \ingroup Matrix_Module
//! \class   TGhostMatrix
//! \brief   The class acts as a place holder for non-local matrix blocks to
//!          access logical information, e.g. size, processor number, but
//!          can not perform any computations
//!
class TGhostMatrix : public TMatrix
{
private:
    //! number of rows in matrix
    size_t  _n_rows;
    
    //! number of columns in matrix
    size_t  _n_cols;

public:
    ///////////////////////////////////////////////////////
    //
    // ctors and dtor
    //

    //! construct zero-sized matrix
    TGhostMatrix ( const value_type_t  avalue_type = real_valued )
            : TMatrix( avalue_type )
            , _n_rows( 0 )
            , _n_cols(0)
    {}

    //! construct matrix defined over given block index set and
    //! on given processor set
    TGhostMatrix ( const TBlockIndexSet &  is,
                   const TProcSet &        ps,
                   const value_type_t      avalue_type = real_valued )
            : TMatrix( avalue_type )
            , _n_rows( is.row_is().size() )
            , _n_cols( is.col_is().size() )
    {
        set_ofs( is.row_is().first(), is.col_is().first() );
        set_procs( ps );
    }

    ////////////////////////////////////////////////////////
    //
    // structure of matrix
    //

    //! return number of rows
    virtual size_t  rows      () const { return _n_rows; }

    //! return number of columns
    virtual size_t  cols      () const { return _n_cols; }

    //! set dimension of matrix
    virtual void    set_size  ( const size_t n, const size_t m )
    {
        _n_rows = n;
        _n_cols = m;
    }

    ///////////////////////////////////////////
    //
    // management of field type
    //

    //! convert matrix data to real valued format
    virtual void to_real    () {}
    
    //! convert matrix data to complex valued format
    virtual void to_complex () {}
    
    /////////////////////////////////////////////////
    //
    // misc.
    //

    // transpose matrix
    virtual void   transpose  ()
    {
        TMatrix::transpose();

        std::swap( _n_rows, _n_cols );
    }  
    
    // conjugate matrix coefficients
    virtual void   conjugate  ()
    {}
    
    //! truncate matrix to given accuracy
    virtual void   truncate   ( const TTruncAcc & )
    {}

    //! return matrix of same class (but no content)
    virtual auto   create     () const -> std::unique_ptr< TMatrix >
    {
        return std::make_unique< TGhostMatrix >();
    }

    //! return size in bytes used by this object
    virtual size_t byte_size  () const
    {
        return TMatrix::byte_size() + sizeof(size_t) * 2;
    }

    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( TGhostMatrix, TMatrix )
};

}// namespace HLIB

#endif  // __HLIB_TGHOSTMATRIX_HH
