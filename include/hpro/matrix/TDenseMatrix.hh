#ifndef __HLIB_TDENSEMATRIX_HH
#define __HLIB_TDENSEMATRIX_HH
//
// Project     : HLib
// File        : TDenseMatrix.hh
// Description : class for dense matrices of arbitrary (small) size
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/base/types.hh"
#include "hpro/blas/Matrix.hh"

#include "hpro/cluster/TPermutation.hh"

#include "hpro/matrix/TMatrix.hh"

#include "hpro/vector/TScalarVector.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TDenseMatrix );

//!
//! \ingroup Matrix_Module
//! \class   TDenseMatrix
//! \brief   Represent a dense matrix
//!
class TDenseMatrix : public TMatrix
{
private:
    //! @cond
    
    //! number of rows 
    size_t                   _rows;

    //! number of columns rows 
    size_t                   _cols;

    //! real valued matrix
    BLAS::Matrix< real >     _rmat;

    //! complex valued matrix
    BLAS::Matrix< complex >  _cmat;

    //! @endcond
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct zero sized matrix
    TDenseMatrix ()
            : _rows(0), _cols(0)
    {}

    //! construct matrix of size \a n × \a m
    TDenseMatrix ( const size_t        n,
                   const size_t        m,
                   const value_type_t  avalue_type = real_valued )
            : TMatrix( avalue_type )
            , _rows( 0 )
            , _cols( 0 )
    {
        set_size(n,m);
    }
    
    //! construct matrix of size \a n × \a m
    TDenseMatrix ( const size_t        n,
                   const size_t        m,
                   const bool          acomplex )
            : TMatrix( acomplex ? complex_valued : real_valued )
            , _rows( 0 )
            , _cols( 0 )
    {
        set_size(n,m);
    }
    
    //! construct matrix with size defined by \a arow_is × \a acol_is
    TDenseMatrix ( const TIndexSet &   arow_is,
                   const TIndexSet &   acol_is,
                   const value_type_t  avalue_type = real_valued )
            : TMatrix( avalue_type )
            , _rows( 0 )
            , _cols( 0 )
    {
        set_block_is( TBlockIndexSet( arow_is, acol_is ) );
    }

    //! construct matrix with size defined by \a arow_is × \a acol_is
    TDenseMatrix ( const TIndexSet &   arow_is,
                   const TIndexSet &   acol_is,
                   const bool          acomplex )
            : TMatrix( acomplex ? complex_valued : real_valued )
            , _rows( 0 )
            , _cols( 0 )
    {
        set_block_is( TBlockIndexSet( arow_is, acol_is ) );
    }

    //! construct matrix with size defined by \a arow_is × \a acol_is
    //! and data given by \a M (real valued)
    TDenseMatrix ( const TIndexSet &             arow_is,
                   const TIndexSet &             acol_is,
                   const BLAS::Matrix< real > &  M )
            : TMatrix( real_valued )
            , _rows( arow_is.size() )
            , _cols( acol_is.size() )
            , _rmat( M )
    {
        // do not call set_block_is to avoid size initialisation
        set_ofs( arow_is.first(), acol_is.first() );
    }

    //! construct matrix with size defined by \a arow_is × \a acol_is
    //! and data given by \a M (real valued)
    TDenseMatrix ( const TIndexSet &                arow_is,
                   const TIndexSet &                acol_is,
                   const BLAS::Matrix< complex > &  M )
            : TMatrix( complex_valued )
            , _rows( arow_is.size() )
            , _cols( acol_is.size() )
            , _cmat( M )
    {
        // do not call set_block_is to avoid size initialisation
        set_ofs( arow_is.first(), acol_is.first() );
    }

    //! construct matrix with size defined by \a arow_is × \a acol_is
    //! and move data from \a M (real valued)
    TDenseMatrix ( const TIndexSet &           arow_is,
                   const TIndexSet &           acol_is,
                   BLAS::Matrix< real > &&     M )
            : TMatrix( real_valued )
            , _rows( arow_is.size() )
            , _cols( acol_is.size() )
            , _rmat( std::move( M ) )
    {
        // do not call set_block_is to avoid size initialisation
        set_ofs( arow_is.first(), acol_is.first() );
    }
    
    //! construct matrix with size defined by \a arow_is × \a acol_is
    //! and move data from \a M (complex valued)
    TDenseMatrix ( const TIndexSet &           arow_is,
                   const TIndexSet &           acol_is,
                   BLAS::Matrix< complex > &&  M )
            : TMatrix( complex_valued )
            , _rows( arow_is.size() )
            , _cols( acol_is.size() )
            , _cmat( std::move( M ) )
    {
        // do not call set_block_is to avoid size initialisation
        set_ofs( arow_is.first(), acol_is.first() );
    }

    //! copy constructor
    TDenseMatrix ( const TDenseMatrix & mat )
            : TMatrix()
            , _rows(0)
            , _cols(0)
    {
        mat.copy_to( this );
    }

    //! construct matrix with size defined by block cluster
    TDenseMatrix ( const TBlockCluster *  bct,
                   const value_type_t     avalue_type = real_valued )
            : TMatrix( bct, avalue_type )
            , _rows(0)
            , _cols(0)
    {
        set_cluster( bct );
    }

    //! destructor
    virtual ~TDenseMatrix ()
    {}

    ////////////////////////////////////////////////
    //
    // manage internal data
    //

    //! return number of rows
    virtual size_t  rows () const { return _rows; }

    //! return number of columns
    virtual size_t  cols () const { return _cols; }

    //! set size of matrix to \a n × \a m
    void            set_size ( const size_t  n,
                               const size_t  m );
    
    //! set size as defined by block cluster \a c
    virtual void    set_cluster ( const TBlockCluster * c );

    //! return true, if matrix is dense
    virtual bool    is_dense () const { return true; }
    
    //
    // return internal BLAS matrices
    //

    //! return real valued matrix
    BLAS::Matrix< real > &          blas_rmat  ()       { return _rmat; }

    //! return constant real valued matrix
    const BLAS::Matrix< real > &    blas_rmat  () const { return _rmat; }

    //! return complex valued matrix
    BLAS::Matrix< complex > &       blas_cmat ()       { return _cmat; }

    //! return constant complex valued matrix
    const BLAS::Matrix< complex > & blas_cmat () const { return _cmat; }

    //
    // return rows and columns
    //

    //! return row \a i as vector object
    const TScalarVector            row      ( const idx_t  i ) const
    {
        if ( is_complex() )
            return TScalarVector( col_is(),
                                  BLAS::Vector< complex >( _cmat, i, BLAS::Range( 0, idx_t(_cols)-1 ) ) );
        else
            return TScalarVector( col_is(),
                                  BLAS::Vector< real >( _rmat, i, BLAS::Range( 0, idx_t(_cols)-1 ) ) );
    }

    //! return column \a i as vector object
    const TScalarVector            column   ( const idx_t  i ) const
    {
        if ( is_complex() )
            return TScalarVector( row_is(),
                                  BLAS::Vector< complex >( _cmat, BLAS::Range( 0, idx_t(_rows)-1 ), i ) );
        else
            return TScalarVector( row_is(),
                                  BLAS::Vector< real >( _rmat, BLAS::Range( 0, idx_t(_rows)-1 ), i ) );
    }

    //! return row \a i as BLAS vector (real value)
    const BLAS::Vector< real >      blas_rrow      ( const idx_t  i ) const
    { return BLAS::Vector< real >( _rmat, i, BLAS::Range( 0, idx_t(_cols)-1 ) ); }

    //! return column \a i as BLAS vector (real value)
    const BLAS::Vector< real >      blas_rcol      ( const idx_t  i ) const 
    { return BLAS::Vector< real >( _rmat, BLAS::Range( 0, idx_t(_rows)-1 ), i ); }

    //! return row \a i as BLAS vector (complex value)
    const BLAS::Vector< complex >   blas_crow      ( const idx_t  i ) const
    { return BLAS::Vector< complex >( _cmat, i, BLAS::Range( 0, idx_t(_cols)-1 ) ); }

    //! return column \a i as BLAS vector (complex value)
    const BLAS::Vector< complex >   blas_ccol      ( const idx_t  i ) const
    { return BLAS::Vector< complex >( _cmat, BLAS::Range( 0, idx_t(_rows)-1 ), i ); }

    //
    // and to pointers to matrix entries
    //

    //! return pointer to data starting at coefficient (\a i, \a j) (real valued)
    real * entry_ptr ( const idx_t i, const idx_t j )
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) entry_ptr", "" );
        return & _rmat( i, j );
    }

    //! return pointer to data starting at coefficient (\a i, \a j) (real valued)
    complex * centry_ptr ( const idx_t i, const idx_t j )
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) entry_ptr", "" );
        return & _cmat( i, j );
    }

    //
    // handle real/complex representation
    //
    
    //! convert to real valued representation if possible
    virtual void  to_real    ();

    //! convert to complex valued representation if possible
    virtual void  to_complex ();
    
    ///////////////////////////////////////////////
    //
    // access coeff.
    //

    //! return coefficient (\a i, \a j) (real valued)
    real entry ( const idx_t i, const idx_t j ) const
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) entry", "" );
        return _rmat( i, j );
    }
    
    //! return coefficient (\a i, \a j) (complex valued)
    const complex centry ( const idx_t i, const idx_t j ) const
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) centry", "" );
        return _cmat( i, j );
    }

    //! set coefficient (\a i, \a j) to \a f (real valued)
    void set_entry  ( const idx_t i, const idx_t j, const real f )
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) set_entry", "" );
        _rmat( i, j ) = f;
    }
    
    //! set coefficient (\a i, \a j) to \a f (complex valued)
    void set_centry ( const idx_t i, const idx_t j, const complex f )
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) set_centry", "" );
        _cmat( i, j ) = f;
    }

    //! add \a f to coefficient (\a i, \a j) (real valued)
    void add_entry  ( const idx_t i, const idx_t j, const real f )
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) add_entry", "" );
        _rmat( i, j ) += f;
    }
    
    //! add \a f to coefficient (\a i, \a j) (complex valued)
    void add_centry ( const idx_t i, const idx_t j, const complex f )
    {
        HASSERT( (i < idx_t(_rows)) && ( j < idx_t(_cols)), ERR_ARR_BOUND, "(TDenseMatrix) add_centry", "" );
        _cmat( i, j ) += f;
    }

    /////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    // scale matrix by constant factor \a f
    virtual void       scale      ( const real            f );
    
    //! compute y ≔ α·op(this)·x + β·y
    virtual void       mul_vec    ( const real            alpha,
                                    const TVector *       x,
                                    const real            beta,
                                    TVector *             y,
                                    const matop_t         op = apply_normal ) const;
    
    //! compute this = this + α·A without truncation
    virtual void       add        ( const real            alpha,
                                    const TMatrix *       A );

    //! compute matrix product α·op_A(this)·op_B(B) 
    virtual TMatrix *  mul_right  ( const real            alpha,
                                    const TMatrix *       B,
                                    const matop_t         op_A,
                                    const matop_t         op_B ) const;

    //! compute matrix product α·op_A(A)·op_B(this) 
    virtual TMatrix *  mul_left   ( const real            alpha,
                                    const TMatrix *       A,
                                    const matop_t         op_A,
                                    const matop_t         op_B ) const;
    
    //! compute α·this + β·op(M), where op(M) is subblock of this or this a subblock of op(M))
    void               add_block  ( const real            alpha,
                                    const real            beta,
                                    const TDenseMatrix *  M,
                                    const matop_t         op = apply_normal );
    
    /////////////////////////////////////////////////
    //
    // BLAS-routines (complex valued)
    //

    // scale matrix by constant factor \a f
    virtual void      cscale      ( const complex         f );
    
    //! compute y ≔ α·op(this)·x + β·y
    virtual void      cmul_vec    ( const complex         alpha,
                                    const TVector *       x,
                                    const complex         beta,
                                    TVector *             y,
                                    const matop_t         op_A = apply_normal ) const;

    //! compute this = this + α·A without truncation
    virtual void      cadd        ( const complex         alpha,
                                    const TMatrix *       A );

    //! compute matrix product α·op_A(this)·op_B(B) 
    virtual TMatrix * cmul_right  ( const complex         alpha,
                                    const TMatrix *       B,
                                    const matop_t         op_A,
                                    const matop_t         op_B ) const;

    //! compute matrix product α·op_A(A)·op_B(this) 
    virtual TMatrix * cmul_left   ( const complex         alpha,
                                    const TMatrix *       A,
                                    const matop_t         op_A,
                                    const matop_t         op_B ) const;
    
    //! compute α·this + β·op(M), where op(M) is subblock of this or this a subblock of op(M))
    void              add_block  ( const complex          alpha,
                                   const complex          beta,
                                   const TDenseMatrix *   M,
                                   const matop_t          op = apply_normal );
    
    ///////////////////////////////////////////////////////////
    //
    // linear operator mapping
    //

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

    using TMatrix::apply_add;

    /////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! transpose matrix
    virtual void transpose ();
    
    //! conjugate matrix coefficients
    virtual void conjugate ();
    
    //! truncate matrix to given accuracy: nothing to be done
    virtual void truncate ( const TTruncAcc & ) {}

    //! copy operator
    TDenseMatrix & operator = ( const TDenseMatrix & mat );

    //
    // virtual constructors
    //

    //! return matrix object of same class as this
    virtual auto  create       () const -> std::unique_ptr< TMatrix >
    {
        return std::make_unique< TDenseMatrix >();
    }

    //! return copy of matrix
    virtual auto  copy         () const -> std::unique_ptr< TMatrix >;
    using TMatrix::copy;

    //! return structural copy of matrix
    virtual auto  copy_struct  () const -> std::unique_ptr< TMatrix >;

    //! copy matrix into A
    virtual void  copy_to      ( TMatrix * A ) const;
    using TMatrix::copy_to;
    
    //! permute rows/columns in matrix according to \a row_pern and \a col_perm
    //! - nullptr is treated as identity
    void permute ( const TPermutation *  row_perm,
                   const TPermutation *  col_perm );

    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( TDenseMatrix, TMatrix )

    //
    // serialisation
    //

    //! read matrix from byte stream
    virtual void    read     ( TByteStream & s );

    //! construct matrix based on data in byte stream
    virtual void    build    ( TByteStream & s );

    //! write matrix to byte stream
    virtual void    write    ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t  bs_size  () const;

    //
    // tests
    //

    //! test data for invalid values, e.g. INF and NAN
    virtual void check_data () const;
};

//
// wrapper for entry_ptr/centry_ptr (relative to global indexset)
//
template <typename T>
T *
entry_ptr                ( const TDenseMatrix *  A,
                           const idx_t           i,
                           const idx_t           j );

template<>
inline real *
entry_ptr< real >        ( const TDenseMatrix *  A,
                           const idx_t           i,
                           const idx_t           j )
{ return const_cast< TDenseMatrix * >( A )->entry_ptr( i - A->row_ofs(), j - A->col_ofs() ); }

template<>
inline complex *
entry_ptr< complex >     ( const TDenseMatrix *  A,
                           const idx_t           i,
                           const idx_t           j )
{ return const_cast< TDenseMatrix * >( A )->centry_ptr( i - A->row_ofs(), j - A->col_ofs() ); }

//
// wrapper for entry_ptr/centry_ptr (relative to local indexset)
//
template <typename T>
T *
entry_ptr_loc            ( const TDenseMatrix *  A,
                           const idx_t           i,
                           const idx_t           j );

template<>
inline real *
entry_ptr_loc< real >    ( const TDenseMatrix *  A,
                           const idx_t           i,
                           const idx_t           j )
{ return const_cast< TDenseMatrix * >( A )->entry_ptr( i, j ); }

template<>
inline complex *
entry_ptr_loc< complex > ( const TDenseMatrix *  A,
                           const idx_t           i,
                           const idx_t           j )
{ return const_cast< TDenseMatrix * >( A )->centry_ptr( i, j ); }

//
// wrapper for set_entry/set_centry
//
template <typename T>
void
set_entry            ( TDenseMatrix *  A,
                       const idx_t     i,
                       const idx_t     j,
                       const T         val );

template<>
inline void
set_entry< real >    ( TDenseMatrix *  A,
                       const idx_t     i,
                       const idx_t     j,
                       const real      val )
{
    return A->set_entry( i, j, val );
}

template<>
inline void
set_entry< complex > ( TDenseMatrix *  A,
                       const idx_t     i,
                       const idx_t     j,
                       const complex   val )
{
    return A->set_centry( i, j, val );
}

//
// wrapper for BLAS matrix
//
template <typename T> const BLAS::Matrix< T > &    blas_mat          ( const TDenseMatrix *  A );
template <> inline const BLAS::Matrix< real > &    blas_mat<real>    ( const TDenseMatrix *  A ) { return A->blas_rmat(); }
template <> inline const BLAS::Matrix< complex > & blas_mat<complex> ( const TDenseMatrix *  A ) { return A->blas_cmat(); }

template <typename T> BLAS::Matrix< T > &          blas_mat          ( TDenseMatrix *        A ); 
template <> inline BLAS::Matrix< real > &          blas_mat<real>    ( TDenseMatrix *        A ) { return A->blas_rmat(); }
template <> inline BLAS::Matrix< complex > &       blas_mat<complex> ( TDenseMatrix *        A ) { return A->blas_cmat(); }

template <typename T> const BLAS::Matrix< T > &    blas_mat          ( const TDenseMatrix &  A );
template <> inline const BLAS::Matrix< real > &    blas_mat<real>    ( const TDenseMatrix &  A ) { return A.blas_rmat(); }
template <> inline const BLAS::Matrix< complex > & blas_mat<complex> ( const TDenseMatrix &  A ) { return A.blas_cmat(); }

template <typename T> BLAS::Matrix< T > &          blas_mat          ( TDenseMatrix &        A ); 
template <> inline BLAS::Matrix< real > &          blas_mat<real>    ( TDenseMatrix &        A ) { return A.blas_rmat(); }
template <> inline BLAS::Matrix< complex > &       blas_mat<complex> ( TDenseMatrix &        A ) { return A.blas_cmat(); }

template <typename T> BLAS::Matrix< T > &          blas_mat          ( std::unique_ptr< TDenseMatrix > &  M ) { return blas_mat<T>( M.get() ); }

//
// wrapper for row/column (relative to local indexset)
//
template <typename T> const BLAS::Vector< T >    blas_row           ( const TDenseMatrix *  A, const idx_t  i );
template <> inline const BLAS::Vector< real >    blas_row<real>     ( const TDenseMatrix *  A, const idx_t  i ) { return A->blas_rrow( i ); }
template <> inline const BLAS::Vector< complex > blas_row<complex>  ( const TDenseMatrix *  A, const idx_t  i ) { return A->blas_crow( i ); }

template <typename T> const BLAS::Vector< T >    blas_col           ( const TDenseMatrix *  A, const idx_t  i );
template <> inline const BLAS::Vector< real >    blas_col<real>     ( const TDenseMatrix *  A, const idx_t  i ) { return A->blas_rcol( i ); }
template <> inline const BLAS::Vector< complex > blas_col<complex>  ( const TDenseMatrix *  A, const idx_t  i ) { return A->blas_ccol( i ); }

template <typename T> const BLAS::Vector< T >    blas_row           ( const TDenseMatrix &  A, const idx_t  i );
template <> inline const BLAS::Vector< real >    blas_row<real>     ( const TDenseMatrix &  A, const idx_t  i ) { return A.blas_rrow( i ); }
template <> inline const BLAS::Vector< complex > blas_row<complex>  ( const TDenseMatrix &  A, const idx_t  i ) { return A.blas_crow( i ); }

template <typename T> const BLAS::Vector< T >    blas_col           ( const TDenseMatrix &  A, const idx_t  i );
template <> inline const BLAS::Vector< real >    blas_col<real>     ( const TDenseMatrix &  A, const idx_t  i ) { return A.blas_rcol( i ); }
template <> inline const BLAS::Vector< complex > blas_col<complex>  ( const TDenseMatrix &  A, const idx_t  i ) { return A.blas_ccol( i ); }

}// namespace HLIB

#endif  // __HLIB_TDENSEMATRIX_HH
