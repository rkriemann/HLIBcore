#ifndef __HLIB_TZEROMATRIX_HH
#define __HLIB_TZEROMATRIX_HH
//
// Project     : HLib
// File        : TZeroMatrix.hh
// Description : class for a zero matrix, i.e. with only zero coefficients
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/cluster/TPermutation.hh"

#include "hpro/matrix/TMatrix.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TZeroMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    TZeroMatrix
//! \brief    Class for a null matrix with only zero coefficients
//
class TZeroMatrix : public TMatrix
{
private:
    //! number of rows in matrix
    size_t                  _rows;

    //! number of columns in matrix
    size_t                  _cols;

public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct null matrix of size \a anrows × \a ancols
    TZeroMatrix ( const size_t  anrows,
                  const size_t  ancols )
            : _rows(0)
            , _cols(0)
    {
        set_size( anrows, ancols );
    }
    
    //! construct null matrix with size defined by block cluster \a bct
    TZeroMatrix ( const TBlockCluster *  bct         = nullptr,
                  const value_type_t     avalue_type = real_valued )
            : TMatrix( bct, avalue_type )
            , _rows(0)
            , _cols(0)
    {
        if ( bct != nullptr )
            set_cluster( bct );
    }

    virtual ~TZeroMatrix () {}

    /////////////////////////////////////////////////
    //
    // access data
    //

    //! set block cluster of matrix
    virtual void  set_cluster  ( const TBlockCluster * bct );

    //! directly set dimension of matrix
    virtual void  set_size     ( const size_t  nrows,
                                 const size_t  ncols );
    
    //! return number of rows in matrix
    size_t        rows         () const { return _rows; }

    //! return number of columns in matrix
    size_t        cols         () const { return _cols; }

    //! convert coefficients to real valued representation (if possible)
    virtual void  to_real      () {}

    //! convert coefficients to complex valued representation
    virtual void  to_complex   () {}

    //! return true, if matrix is zero
    virtual bool  is_zero      () const { return true; }
    
    /////////////////////////////////////////////////
    //
    // manage stored entries
    //

    //! return matrix coefficient a_ij (real valued)
    virtual real           entry      ( const idx_t,
                                        const idx_t ) const { return real(0); }

    //! return matrix coefficient a_ij (complex valued)
    virtual const complex  centry     ( const idx_t,
                                        const idx_t ) const { return complex(0); }

    /////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    //! compute this ≔ α·this
    virtual void scale ( const real ) {}
    
    //! compute y ≔ β·y + α·op(M)·x, with M = this
    virtual void mul_vec ( const real      alpha,
                           const TVector * x,
                           const real      beta,
                           TVector       * y,
                           const matop_t   op = MATOP_NORM ) const;

    //! compute this ≔ this + α · matrix
    virtual void add ( const real alpha, const TMatrix * matrix );
        
    //! transpose matrix
    virtual void transpose ();
    
    //! conjugate matrix coefficients
    virtual void conjugate ();
    
    //! truncate matrix to given accuracy (NOT YET IMPLEMENTED)
    virtual void truncate ( const TTruncAcc & ) {}
    
    /////////////////////////////////////////////////
    //
    // BLAS-routines (complex valued)
    //

    //! compute this ≔ α·this
    virtual void cscale ( const complex ) {}
    
    //! compute y ≔ β·y + α·op(M)·x, with M = this
    virtual void cmul_vec ( const complex   alpha,
                            const TVector * x,
                            const complex   beta,
                            TVector       * y,
                            const matop_t   op = MATOP_NORM ) const;

    //! compute this ≔ this + α · matrix
    virtual void cadd ( const complex a, const TMatrix * matrix );
        
    /////////////////////////////////////////////////
    //
    // serialisation
    //

    //! read data from stream \a s and copy to matrix
    virtual void read  ( TByteStream & s );
    
    //! use data from stream \a s to build matrix
    virtual void build ( TByteStream & s );

    //! write data to stream \a s
    virtual void write ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;

    /////////////////////////////////////////////////
    //
    // misc.
    //

    //
    // virtual ctors
    //
    
    //! return matrix of same class (but no content)
    virtual auto  create       () const -> std::unique_ptr< TMatrix > { return std::make_unique< TZeroMatrix >(); }
    
    //! return copy of matrix
    virtual auto  copy         () const -> std::unique_ptr< TMatrix >;
    using TMatrix::copy;

    //! return structural copy of matrix
    virtual auto  copy_struct  () const -> std::unique_ptr< TMatrix >;

    //! copy matrix into matrix \a A
    virtual void copy_to       ( TMatrix * A ) const;
    using TMatrix::copy_to;

    //
    // type checking
    //

    HLIB_RTTI_DERIVED( TZeroMatrix, TMatrix );

    //! return size in bytes used by this object
    virtual size_t byte_size () const;
};

}// namespace HLIB

#endif  // __HLIB_TZEROMATRIX_HH
