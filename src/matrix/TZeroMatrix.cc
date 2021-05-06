//
// Project     : HLib
// File        : TZeroMatrix.cc
// Description : class for a zero matrix, i.e. with only zero coefficients
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/matrix/TZeroMatrix.hh"

namespace HLIB
{

/////////////////////////////////////////////////
//
// access data
//

//! set block cluster of matrix
void
TZeroMatrix::set_cluster  ( const TBlockCluster * bct )
{
    TMatrix::set_cluster( bct );

    if ( bct != nullptr )
        set_size( bct->rowcl()->size(),
                  bct->colcl()->size() );
}

//! directly set dimension of matrix
void
TZeroMatrix::set_size ( const size_t  nrows,
                        const size_t  ncols )
{
    _rows = nrows;
    _cols = ncols;
}
    
/////////////////////////////////////////////////
//
// BLAS-routines (real valued)
//

//! compute y ≔ β·y + α·op(M)·x, with M = this
void
TZeroMatrix::mul_vec ( const real,
                       const TVector *,
                       const real      beta,
                       TVector       * y,
                       const matop_t ) const
{
    y->scale( beta );
}

//! compute this ≔ this + α · matrix
void
TZeroMatrix::add ( const real, const TMatrix * M )
{
    if ( ! IS_TYPE( M, TZeroMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TZeroMatrix) cadd", M->typestr() );
}

//
// transpose matrix
//
void
TZeroMatrix::transpose ()
{
    TMatrix::transpose();
    
    std::swap( _rows, _cols );
}

//
// conjugate matrix coefficients
//
void
TZeroMatrix::conjugate ()
{
}

/////////////////////////////////////////////////
//
// BLAS-routines (complex valued)
//

//! compute y ≔ β·y + α·op(M)·x, with M = this
void
TZeroMatrix::cmul_vec ( const complex,
                        const TVector *,
                        const complex   beta,
                        TVector       * y,
                        const matop_t ) const
{
    y->cscale( beta );
}

//! compute this ≔ this + α · matrix
void
TZeroMatrix::cadd ( const complex,
                    const TMatrix *  M )
{
    if ( ! IS_TYPE( M, TZeroMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TZeroMatrix) cadd", M->typestr() );
}
        
//
// serialisation
//

//! read data from stream \a s and copy to matrix
void
TZeroMatrix::read  ( TByteStream & s )
{
    TMatrix::read( s );

    s.get( _rows );
    s.get( _cols );
}
    
//! use data from stream \a s to build matrix
void
TZeroMatrix::build ( TByteStream & s )
{
    TMatrix::build( s );

    s.get( _rows );
    s.get( _cols );
}

//! write data to stream \a s
void
TZeroMatrix::write ( TByteStream & s ) const
{
    TMatrix::write( s );

    s.put( _rows );
    s.put( _cols );
}

//! returns size of object in bytestream
size_t
TZeroMatrix::bs_size () const
{
    return TMatrix::bs_size() + + sizeof(_rows) + sizeof(_cols);
}

/////////////////////////////////////////////////
//
// misc.
//

//
// return copy of matrix
//
std::unique_ptr< TMatrix >
TZeroMatrix::copy  () const
{
    auto           M = TMatrix::copy();
    TZeroMatrix *  Z = ptrcast( M.get(), TZeroMatrix );

    copy_to( Z );

    return M;
}

//
// return structural copy of matrix
//
std::unique_ptr< TMatrix >
TZeroMatrix::copy_struct  () const
{
    auto           M = TMatrix::copy_struct();
    TZeroMatrix *  Z = ptrcast( M.get(), TZeroMatrix );

    copy_to( Z );

    return M;
}

//
// copy matrix into matrix \a A
//
void
TZeroMatrix::copy_to ( TMatrix * A ) const
{
    if ( ! IS_TYPE( A, TZeroMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TRkMatrix) copy_to", A->typestr() );

    TZeroMatrix *  M = ptrcast( A, TZeroMatrix );

    M->set_complex( is_complex() );
    M->set_size( rows(), cols() );
    M->set_ofs( row_ofs(), col_ofs() );
}

//
// return size in bytes used by this object
//
size_t
TZeroMatrix::byte_size () const
{
    return TMatrix::byte_size() + sizeof(_rows) + sizeof(_cols);
}

}// namespace HLIB
