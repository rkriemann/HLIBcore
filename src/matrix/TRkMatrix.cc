//
// Project     : HLib
// File        : TRkMatrix.cc
// Description : class for rank-k-matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/blas/Algebra.hh"
#include "hpro/vector/TBlockVector.hh"
#include "hpro/vector/vec_conv.hh"
#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/structure.hh"

#include "hpro/matrix/TRkMatrix.hh"

namespace HLIB
{

namespace B = BLAS;

using std::unique_ptr;
using std::make_unique;

namespace
{

template < typename value_t >
void
task_mulvec ( const value_t                 alpha,
              const matop_t                 op,
              const B::Matrix< value_t > &  A,
              const B::Matrix< value_t > &  B,
              const B::Vector< value_t > &  x,
              B::Vector< value_t > &        y )
{
    const bool  test_dense = false;

    B::Vector< value_t >  res_exact;
    
    if ( test_dense )
    {
        const auto  D = B::prod( value_t(1), A, B::adjoint( B ) );

        res_exact = B::mulvec( alpha, B::mat_view( op, D ), x );
        B::add( value_t(1), y, res_exact );
    }// if

    // B::Matrix< value_t >  CA( A, copy_value );
    // B::Matrix< value_t >  CB( B, copy_value );
    // B::Vector< value_t >  cx( x, copy_value );
    // B::Vector< value_t >  cy( y, copy_value );
    
    const auto            rank = A.ncols();
    B::Vector< value_t >  v( rank );
                
    if ( op == apply_normal )
    {
        B::mulvec( alpha, B::adjoint( B ), x, value_t(0), v );
        B::mulvec( value_t(1), A, v, value_t(1), y );
    }// if
    else if ( is_complex_type< value_t >::value && ( op == apply_transposed ))
    {
        //
        // (AB^H)^T x = conj(B) A^T x
        //   1. v = A^T x
        //   2. t = B conj(v)
        //   3. y = y + conj(t)
        //
        
        B::mulvec( alpha, B::transposed( A ), x, value_t(0), v );
        B::conj( v );

        {
            B::Vector< value_t >  t( y.length() );
            
            B::mulvec( value_t(1), B, v, value_t(0), t );
            B::conj( t );
            B::add( value_t(1), t, y );
        }// else
    }// else
    else // if ( op == apply_adjoint ) or value_t == real and op == transposed/adjoint
    {
        B::mulvec( alpha, B::adjoint( A ), x, value_t(0), v );
        B::mulvec( value_t(1), B, v, value_t(1), y );
    }// else

    if ( test_dense )
    {
        B::add( value_t(-1), y, res_exact );

        const auto  f = B::norm2( res_exact );
        
        if ( f > 1e-10 )
        {
            LOG::printf( "%.8e", f );
        }// if
    }// if
}

}// namespace anonymous

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// TRkMatrix
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
//
// constructor and destructor
//

TRkMatrix::TRkMatrix ()
        : _rows( 0 )
        , _cols( 0 )
        , _rank( 0 )
{}

TRkMatrix::TRkMatrix ( const size_t  r,
                       const size_t  c )
        : _rows( r )
        , _cols( c )
        , _rank( 0 )
{}

TRkMatrix::TRkMatrix ( const TIndexSet &   arow_is,
                       const TIndexSet &   acol_is,
                       const value_type_t  avalue_type )
        : TMatrix( avalue_type )
        , _rows( 0 )
        , _cols( 0 )
        , _rank( 0 )
{
    set_block_is( TBlockIndexSet( arow_is, acol_is ) );
}

TRkMatrix::TRkMatrix ( const TIndexSet &   arow_is,
                       const TIndexSet &   acol_is,
                       const bool          acomplex )
        : TMatrix( acomplex ? complex_valued : real_valued )
        , _rows( 0 )
        , _cols( 0 )
        , _rank( 0 )
{
    set_block_is( TBlockIndexSet( arow_is, acol_is ) );
}

TRkMatrix::TRkMatrix ( const TIndexSet &          arow_is,
                       const TIndexSet &          acol_is,
                       const B::Matrix< real > &  A,
                       const B::Matrix< real > &  B )
        : TMatrix( real_valued )
        , _rmat_A( A )
        , _rmat_B( B )
        , _rows( arow_is.size() )
        , _cols( acol_is.size() )
        , _rank( _rmat_A.ncols() )
{
    if ( A.ncols() != B.ncols() )
        HERROR( ERR_ARG, "(TRkMatrix) ctor", "rank of A and B differs" );

    // do not call set_block_is to avoid size initialisation
    set_ofs( arow_is.first(), acol_is.first() );
}

TRkMatrix::TRkMatrix ( const TIndexSet &             arow_is,
                       const TIndexSet &             acol_is,
                       const B::Matrix< complex > &  A,
                       const B::Matrix< complex > &  B )
        : TMatrix( complex_valued ),
          _cmat_A( A ), _cmat_B( B ),
          _rows( arow_is.size() ),
          _cols( acol_is.size() ),
          _rank( _cmat_A.ncols() )
{
    if ( A.ncols() != B.ncols() )
        HERROR( ERR_ARG, "(TRkMatrix) ctor", "rank of A and B differs" );

    // do not call set_block_is to avoid size initialisation
    set_ofs( arow_is.first(), acol_is.first() );
}

TRkMatrix::TRkMatrix ( const TIndexSet &     arow_is,
                       const TIndexSet &     acol_is,
                       B::Matrix< real > &&  A,
                       B::Matrix< real > &&  B )
        : TMatrix( real_valued )
        , _rmat_A( std::move( A ) )
        , _rmat_B( std::move( B ) )
        , _rows( arow_is.size() )
        , _cols( acol_is.size() )
        , _rank( _rmat_A.ncols() )
{
    if ( A.ncols() != B.ncols() )
        HERROR( ERR_ARG, "(TRkMatrix) ctor", "rank of A and B differs" );

    // do not call set_block_is to avoid size initialisation
    set_ofs( arow_is.first(), acol_is.first() );
}

TRkMatrix::TRkMatrix ( const TIndexSet &        arow_is,
                       const TIndexSet &        acol_is,
                       B::Matrix< complex > &&  A,
                       B::Matrix< complex > &&  B )
        : TMatrix( complex_valued )
        , _cmat_A( std::move( A ) )
        , _cmat_B( std::move( B ) )
        , _rows( arow_is.size() )
        , _cols( acol_is.size() )
        , _rank( _cmat_A.ncols() )
{
    if ( A.ncols() != B.ncols() )
        HERROR( ERR_ARG, "(TRkMatrix) ctor", "rank of A and B differs" );

    // do not call set_block_is to avoid size initialisation
    set_ofs( arow_is.first(), acol_is.first() );
}

TRkMatrix::TRkMatrix ( const TBlockIndexSet &  ablock_is,
                       const value_type_t      avalue_type )
        : TMatrix( avalue_type )
        , _rows( 0 )
        , _cols( 0 )
        , _rank( 0 )
{
    set_block_is( ablock_is );
}

TRkMatrix::TRkMatrix ( const TBlockCluster * bcl,
                       const value_type_t    avalue_type )
        : TMatrix( bcl, avalue_type )
        , _rows( 0 )
        , _cols( 0 )
        , _rank( 0 )
{
    if ( bcl != nullptr )
    {
        set_size( bcl->rowcl()->size(), bcl->colcl()->size() );
    }// if
}

TRkMatrix::TRkMatrix ( const TRkMatrix &  A )
        : TMatrix()
        , _rows( 0 )
        , _cols( 0 )
        , _rank( 0 )
{
    set_cluster( A.cluster() );
    set_complex( A.is_complex() );
    set_size( A.rows(), A.cols(), A.rank() );
    set_ofs( A.row_ofs(), A.col_ofs() );
    
    if ( is_complex() )
    {
        B::copy( A._cmat_A, _cmat_A );
        B::copy( A._cmat_B, _cmat_B );
    }// if
    else
    {
        B::copy( A._rmat_A, _rmat_A );
        B::copy( A._rmat_B, _rmat_B );
    }// if
}

//
// set rank of matrix
//

void
TRkMatrix::set_rank ( const size_t  k )
{
    if ( k == _rank )
        return;

    // allocate memory for new arrays
    set_size( _rows, _cols, k );
}

//
// usual matrix access (read-only)
//
real
TRkMatrix::entry ( const idx_t i, const idx_t j ) const
{
    if ( is_complex() )
        HERROR( ERR_REAL_CMPLX, "(TRkMatrix) entry", "matrix is complex valued" );

    real  f = real(0);
        
    for ( uint k = 0; k < _rank; k++ )
        f += _rmat_A( i, k ) * _rmat_B( j, k );
    
    return f;
}

const complex
TRkMatrix::centry ( const idx_t i, const idx_t j ) const
{
    if ( is_complex() )
    {
        complex  f = complex(0);

        for ( uint k = 0; k < _rank; k++ )
            f += _cmat_A( i, k ) * conj( _cmat_B( j, k ) );

        return f;
    }// if
    else
        return entry( i, j );
}

//
// update size of matrices if cluster changed
//
void
TRkMatrix::set_cluster ( const TBlockCluster * c )
{
    TMatrix::set_cluster( c );

    if ( c != nullptr )
        set_size( c->rowcl()->size(), c->colcl()->size(), _rank );
}

//
// switch between complex and real format
//
void
TRkMatrix::to_real ()
{
    if ( ! is_complex() )
        return;

    if ( _rank == 0 )
        return;
    
    TScopedLock  mlock( *this );

    // check if matrix has imaginary part
    for ( uint k = 0; k < _rank; ++k )
        for ( uint i = 0; i < _rows; ++i )
        {
            if ( std::imag( _cmat_A( i, k ) ) != real(0) )
            {
                TMatrix::set_complex( true );
                HERROR( ERR_COMPLEX, "(TRkMatrix) to_real", "matrix has imaginary part" );
            }// if
        }// for
        
    for ( uint k = 0; k < _rank; ++k )
        for ( uint i = 0; i < _cols; ++i )
        {
            if ( std::imag( _cmat_B( i, k ) ) != real(0) )
            {
                TMatrix::set_complex( true );
                HERROR( ERR_COMPLEX, "(TRkMatrix) to_real", "matrix has imaginary part" );
            }// if
        }// for

    B::Matrix< real > TA( _rows, _rank );
    B::Matrix< real > TB( _cols, _rank );

    for ( uint k = 0; k < _rank; ++k )
        for ( uint i = 0; i < _rows; ++i )
            TA( i, k ) = std::real( _cmat_A( i, k ) );

    for ( uint k = 0; k < _rank; ++k )
        for ( uint i = 0; i < _cols; ++i )
            TB( i, k ) = std::real( _cmat_B( i, k ) );

    _cmat_A = B::Matrix< complex >();
    _cmat_B = B::Matrix< complex >();
    _rmat_A = std::move( TA );
    _rmat_B = std::move( TB );
}

void
TRkMatrix::to_complex ()
{
    if ( is_complex() )
        return;

    if ( _rank == 0 )
        return;
    
    TScopedLock           mlock( *this );
    B::Matrix< complex >  TA( _rows, _rank );
    B::Matrix< complex >  TB( _cols, _rank );
        
    for ( uint k = 0; k < _rank; ++k )
        for ( uint i = 0; i < _rows; ++i )
            TA( i, k ) = _rmat_A( i, k );

    for ( uint k = 0; k < _rank; ++k )
        for ( uint i = 0; i < _cols; ++i )
            TB( i, k ) = _rmat_B( i, k );

    _rmat_A = B::Matrix< real >();
    _rmat_B = B::Matrix< real >();
    _cmat_A = std::move( TA );
    _cmat_B = std::move( TB );
}

/////////////////////////////////////////////////
//
// management of update accumulator
//

//
// apply stored updates U to local matrix M, e.g., M = M + U,
// with accuracy \a acc
//
void
TRkMatrix::apply_updates ( const TTruncAcc &       /* acc */,
                           const recursion_type_t )
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

/////////////////////////////////////////////////
//
// truncate the rank via bestapproximation in frobenius-norm
//

void
TRkMatrix::truncate ( const TTruncAcc & acc )
{
    if ( _rank == 0 )
        return;
        
    TScopedLock  mlock( *this );

    if ( is_complex() )
    {
        _rank = B::truncate( _cmat_A, _cmat_B, acc( this ) );

        if (( _cmat_A.ncols() != _rank ) || ( _cmat_B.ncols() != _rank ))
            HERROR( ERR_MAT_SIZE, "(TRkMatrix) truncate", "" );
    }// if
    else
    {
        _rank = B::truncate( _rmat_A, _rmat_B, acc( this ) );

        if (( _rmat_A.ncols() != _rank ) || ( _rmat_B.ncols() != _rank ))
            HERROR( ERR_MAT_SIZE, "(TRkMatrix) truncate", "" );
    }// else

}

//
// copy given dense matrix into local rank-matrix
//
void
TRkMatrix::copy_dense ( const TDenseMatrix * A, const TTruncAcc & acc )
{
    if ((A->rows() != _rows) || (A->cols() != _cols))
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) copy_dense", "argument has wrong dimensions" );

    set_complex( A->is_complex() );
    
    TScopedLock   mlock( *this );
    const size_t  n = A->rows();
    const size_t  m = A->cols();

    if ( A->is_complex() )
    {
        if ( acc.is_exact() )
        {
            //
            // copy dense matrix directly
            //

            if ( A->rows() > A->cols() )
            {
                const size_t          rk = A->cols();
                B::Matrix< complex >  TA( rows(), rk );
                B::Matrix< complex >  TB( cols(), rk );
                
                // copy A to R_A
                B::copy( A->blas_cmat(), TA );
                
                // and set R_B to e_1 ... e_rk
                B::fill( complex(0), TB );
                for ( uint i = 0; i < rk; i++ )
                    TB( i, i ) = complex(1);

                _cmat_A = std::move( TA );
                _cmat_B = std::move( TB );
                _rank   = rk;
            }// if
            else
            {
                const size_t          rk = A->rows();
                B::Matrix< complex >  TA( rows(), rk );
                B::Matrix< complex >  TB( cols(), rk );
                
                // copy A^T to R_B; conjugate since B^H == A
                B::copy( B::adjoint( A->blas_cmat() ), TB );
                
                // and set R_A to e_1 ... e_rk
                B::fill( complex(0), TA );
                for ( uint i = 0; i < rk; i++ )
                    TA( i, i ) = complex(1);

                _cmat_A = std::move( TA );
                _cmat_B = std::move( TB );
                _rank   = rk;
            }// else
        }// if
        else
        {
            //
            // compute low-rank approximation of dense matrix
            //
            
            B::Matrix< complex > D( n, m );
        
            B::copy( A->blas_cmat(), D );
            _rank = B::approx( D, acc( this ), _cmat_A, _cmat_B );
        }// else

        if (( _cmat_A.ncols() != _rank ) || ( _cmat_B.ncols() != _rank ))
            HERROR( ERR_MAT_SIZE, "(TRkMatrix) copy_dense", "" );
    }// if
    else
    {
        if ( acc.is_exact() )
        {
            //
            // copy dense matrix directly
            //

            if ( A->rows() > A->cols() )
            {
                const size_t       rk = A->cols();
                B::Matrix< real >  TA( rows(), rk );
                B::Matrix< real >  TB( cols(), rk );
                
                // copy A to R_A
                B::copy( A->blas_rmat(), TA );
                
                // and set R_B to e_1 ... e_rk
                B::fill( real(0), TB );
                for ( uint i = 0; i < rk; i++ )
                    TB( i, i ) = real(1);

                _rmat_A = std::move( TA );
                _rmat_B = std::move( TB );
                _rank   = rk;
            }// if
            else
            {
                const size_t       rk = A->rows();
                B::Matrix< real >  TA( rows(), rk );
                B::Matrix< real >  TB( cols(), rk );
                
                // copy A^T to R_B; conjugate since B^H == A
                B::copy( B::adjoint( A->blas_rmat() ), TB );
                
                // and set R_A to e_1 ... e_rk
                B::fill( real(0), TA );
                for ( uint i = 0; i < rk; i++ )
                    TA( i, i ) = real(1);

                _rmat_A = std::move( TA );
                _rmat_B = std::move( TB );
                _rank   = rk;
            }// else
        }// if
        else
        {
            //
            // compute low-rank approximation of dense matrix
            //
            
            B::Matrix< real > D( n, m );
        
            B::copy( A->blas_rmat(), D );
            _rank = B::approx( D, acc( this ), _rmat_A, _rmat_B );
        }// else

        if (( _rmat_A.ncols() != _rank ) || ( _rmat_B.ncols() != _rank ))
            HERROR( ERR_MAT_SIZE, "(TRkMatrix) copy_dense", "" );
    }// else
}

//
// set this ≔ A·B^H
//
void
TRkMatrix::set_lrmat ( const BLAS::Matrix< real > &     A,
                       const BLAS::Matrix< real > &     B )
{
    if ( is_complex() )
        HERROR( ERR_ARG, "(TRkMatrix) set_lrmat", "expecting real valued data" );
        
    if (( A.nrows() != _rmat_A.nrows() ) ||
        ( B.nrows() != _rmat_B.nrows() ) ||
        ( A.ncols() != B.ncols()))
        HERROR( ERR_ARG, "(TRkMatrix) set_lrmat", "input matrices have invalid dimension" );

    if ( _rank == A.ncols() )
    {
        B::copy( A, _rmat_A );
        B::copy( B, _rmat_B );
    }// if
    else
    {
        _rmat_A = std::move( B::Matrix< real >( A, copy_value ) );
        _rmat_B = std::move( B::Matrix< real >( B, copy_value ) );
        _rank   = A.ncols();
    }// else
}

void
TRkMatrix::set_lrmat ( const BLAS::Matrix< complex > &  A,
                       const BLAS::Matrix< complex > &  B )
{
    if ( ! is_complex() )
        HERROR( ERR_ARG, "(TRkMatrix) set_lrmat", "expecting complex valued data" );
        
    if (( A.nrows() != _cmat_A.nrows() ) ||
        ( B.nrows() != _cmat_B.nrows() ) ||
        ( A.ncols() != B.ncols()))
        HERROR( ERR_ARG, "(TRkMatrix) set_lrmat", "input matrices have invalid dimension" );

    if ( _rank == A.ncols() )
    {
        B::copy( A, _cmat_A );
        B::copy( B, _cmat_B );
    }// if
    else
    {
        _cmat_A = std::move( B::Matrix< complex >( A, copy_value ) );
        _cmat_B = std::move( B::Matrix< complex >( B, copy_value ) );
        _rank   = A.ncols();
    }// else
}

//
// add a rank-k-matrix and truncate result
//
void
TRkMatrix::add_rank ( const real                 alpha,
                      const B::Matrix< real > &  A,
                      const B::Matrix< real > &  B,
                      const TTruncAcc &          acc )
{
    if ( A.ncols() != B.ncols() )
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) add_rank", "given matrices have different ncols" );

    const size_t  k = A.ncols();
    
    // check trivial cases
    if (( k == 0 ) || ( alpha == real(0) ))
        return;

    if ( is_complex() )
        HERROR( ERR_REAL_CMPLX, "(TRkMatrix) add_rank", "" );

    //
    // build new matrix holding the sum
    //

    TScopedLock        mlock( *this );
    B::Matrix< real >  TA( _rows, _rank + k );
    B::Matrix< real >  TB( _cols, _rank + k );

    // copy this
    if ( _rank > 0 )
    {
        const B::Range     old_vecs( 0, idx_t(_rank)-1 );
        B::Matrix< real >  partA( TA, B::Range( 0, idx_t(_rows)-1 ), old_vecs );
        B::Matrix< real >  partB( TB, B::Range( 0, idx_t(_cols)-1 ), old_vecs );

        B::copy( _rmat_A, partA );
        B::copy( _rmat_B, partB );
    }// if

    // copy given matrix
    if ( k > 0 )
    {
        const B::Range     new_vecs( idx_t(_rank), idx_t(_rank + k) - 1 );
        B::Matrix< real >  partA( TA, B::Range( 0, idx_t(_rows)-1 ), new_vecs );
        B::Matrix< real >  partB( TB, B::Range( 0, idx_t(_cols)-1 ), new_vecs );

        B::copy( A, partA );
        B::copy( B, partB );

        if ( alpha != real(1) )
        {
            if ( A.nrows() < B.nrows() ) B::scale( alpha, partA );
            else                         B::scale( alpha, partB );
        }// if
    }// if

    //
    // truncate result
    //

    if ( CFG::Arith::zero_sum_trunc || ( _rank > 0 ))
        _rank = B::truncate( TA, TB, acc( this ) );
    else
        _rank = k;

    _rmat_A = std::move( TA );
    _rmat_B = std::move( TB );

    if (( _rmat_A.ncols() != _rank ) || ( _rmat_B.ncols() != _rank ))
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) add_rank", "" );
}

void
TRkMatrix::add_rank ( const complex                 alpha,
                      const B::Matrix< complex > &  A,
                      const B::Matrix< complex > &  B,
                      const TTruncAcc &             acc )
{
    if ( A.ncols() != B.ncols() )
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) add_rank", "given matrices have different ncols" );

    const size_t  k = A.ncols();
    
    // check trivial case
    if (( k == 0 ) || ( alpha == complex(0) ))
        return;

    set_complex( true );
    
    //
    // build new matrix holding the sum
    //

    TScopedLock           mlock( *this );
    B::Matrix< complex >  TA( _rows, _rank + k );
    B::Matrix< complex >  TB( _cols, _rank + k );

    // copy this
    if ( _rank > 0 )
    {
        const B::Range        old_vecs( 0, idx_t(_rank)-1 );
        B::Matrix< complex >  partA( TA, B::Range( 0, idx_t(_rows)-1 ), old_vecs );
        B::Matrix< complex >  partB( TB, B::Range( 0, idx_t(_cols)-1 ), old_vecs );

        B::copy( _cmat_A, partA );
        B::copy( _cmat_B, partB );
    }// if

    // copy given matrix
    if ( k > 0 )
    {
        const B::Range        new_vecs( idx_t(_rank), idx_t(_rank + k) - 1 );
        B::Matrix< complex >  partA( TA, B::Range( 0, idx_t(_rows)-1 ), new_vecs );
        B::Matrix< complex >  partB( TB, B::Range( 0, idx_t(_cols)-1 ), new_vecs );

        B::copy( A, partA );
        B::copy( B, partB );

        if ( alpha != complex(1) )
        {
            if ( A.nrows() < B.nrows() ) B::scale( alpha,       partA );
            else                         B::scale( conj(alpha), partB );
        }// if
    }// if

    //
    // truncate result
    //

    if ( CFG::Arith::zero_sum_trunc || ( _rank > 0 ))
        _rank = B::truncate( TA, TB, acc( this ) );
    else
        _rank = k;

    _cmat_A = std::move( TA );
    _cmat_B = std::move( TB );
}

//
// add a dense matrix and truncate
//
void
TRkMatrix::add_dense ( const real                 alpha,
                       const B::Matrix< real > &  D,
                       const TTruncAcc &          acc )
{
    if ( alpha == real(0) )
        return;
    
    TScopedLock  mlock( *this );

    if ( is_complex() )
    {
        HERROR( ERR_REAL_CMPLX, "(TRkMatrix) add_dense", "" );
    }// if
    else
    {
        B::Matrix< real >  TD( _rows, _cols );

        B::copy( D, TD );

        if ( alpha != real(1) )
            B::scale( alpha, TD );
        
        // add this to D
        if ( _rank > 0 )
            B::prod( real(1), _rmat_A, adjoint(_rmat_B), real(1), TD );

        // truncate and store result in this
        _rank = B::approx( TD, acc( this ), _rmat_A, _rmat_B );

        if (( _rmat_A.ncols() != _rank ) || ( _rmat_B.ncols() != _rank ))
            HERROR( ERR_MAT_SIZE, "(TRkMatrix) add_dense", "" );
    }// if
}
    
void
TRkMatrix::add_dense ( const complex                 alpha,
                       const B::Matrix< complex > &  D,
                       const TTruncAcc &             acc )
{
    if ( alpha == complex(0) )
        return;
    
    set_complex( true );
    
    TScopedLock           mlock( *this );
    B::Matrix< complex >  TD( _rows, _cols );
    
    B::copy( D, TD );

    if ( alpha != complex(1) )
        B::scale( alpha, TD );
    
    // add this to D
    if ( _rank > 0 )
        B::prod( complex(1), _cmat_A, adjoint(_cmat_B), complex(1), TD );
    
    // truncate and store result in this
    _rank = B::approx( TD, acc( this ), _cmat_A, _cmat_B );
}

/////////////////////////////////////////////////
//
// BLAS-routines (real valued)
//

//
// scale matrix by constant factor
//
void
TRkMatrix::scale ( const real f )
{
    if (( f == real(1) ) || ( rank() == 0 ))
        return;

    TScopedLock  mlock( *this );

    if ( is_complex() )
    {
        if ( f == real(0) )
        {
            B::fill( complex(0), blas_cmat_A() );
            B::fill( complex(0), blas_cmat_B() );
        }// if
        else
        {
            //
            // scale matrix with lower norm
            //

            if ( B::normF( blas_cmat_A() ) < B::normF( blas_cmat_B() ) )
                B::scale( complex(f), blas_cmat_A() );
            else
                B::scale( complex(f), blas_cmat_B() );
        }// else
    }// if
    else
    {
        if ( f == real(0) )
        {
            B::fill( real(0), blas_rmat_A() );
            B::fill( real(0), blas_rmat_B() );
        }// if
        else
        {
            //
            // scale matrix with lower norm
            //

            if ( B::normF( blas_rmat_A() ) < B::normF( blas_rmat_B() ) )
                B::scale( f, blas_rmat_A() );
            else
                B::scale( f, blas_rmat_B() );
        }// else
    }// else
}
    
//
// matrix-vector-multiplication
//
void
TRkMatrix::mul_vec ( const real      alpha,
                     const TVector * x,
                     const real      beta,
                     TVector       * y,
                     const matop_t   op ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TRkMatrix) mul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TRkMatrix) mul_vec", "y = nullptr" );

    if ( op == apply_normal )
    {
        if (( row_is() != y->is() ) || ( col_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TRkMatrix) mul_vec", "incompatible vector index set" );
    }// if
    else
    {
        if (( col_is() != y->is() ) || ( row_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TRkMatrix) mul_vec", "incompatible vector index set" );
    }// if
    
    //
    // check if scalar vector
    //
    
    if ( IS_TYPE( x, TScalarVector ) && IS_TYPE( y, TScalarVector ) )
    {
        TScalarVector       * sy = ptrcast( y, TScalarVector );
        const TScalarVector * sx = cptrcast( x, TScalarVector );

        //
        // we assume, that the index-sets are coherent
        // and we can access by local indices
        //

        // scale left-side
        if      ( beta == real(0) ) sy->fill( real(0) );
        else if ( beta != real(1) ) sy->scale( beta );

        if (( alpha == real(0) ) || ( _rank == 0 ))
            return;

        if ( is_complex() || x->is_complex() )
            y->set_complex( true );

        if ( is_complex() )
        {
            if ( ! x->is_complex() )
                HERROR( ERR_NOT_IMPL, "(TRkMatrix) mul_vec", "complex * real" );

            task_mulvec< complex >( alpha, op, blas_cmat_A(), blas_cmat_B(), sx->blas_cvec(),
                                    sy->blas_cvec() );
        }// if
        else
        {
            if ( x->is_complex() )
                HERROR( ERR_NOT_IMPL, "(TRkMatrix) mul_vec", "real * complex" );

            task_mulvec< real >( alpha, op, blas_rmat_A(), blas_rmat_B(), sx->blas_rvec(),
                                 sy->blas_rvec() );
        }// else
    }// if
    else if ( IS_TYPE( x, TBlockVector ) && IS_TYPE( y, TBlockVector ) )
    {
        unique_ptr< TVector >  sx( to_scalar( x ) );
        TScalarVector          sy( y->is(), y->is_complex() );

        mul_vec( alpha, sx.get(), beta, & sy, op );

        y->axpy( real(1), & sy );
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TRkMatrix) mul_vec",
                y->typestr() + " += TRkMatrix * " + x->typestr() );
}

//
// compute this ≔ this + α · M
//
void
TRkMatrix::add ( const real alpha, const TMatrix * M )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "(TRkMatrix) add", "nullptr argument" );
    
    if ( ! IS_TYPE( M, TRkMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TRkMatrix) add", M->typestr() );

    const TRkMatrix * R = cptrcast( M, TRkMatrix );

    //
    // some trivial checks
    //

    if (( alpha == real(0) ) || ( R->rank() == 0 ))
        return;

    if ((rows() != R->rows()) || (cols() != R->cols()))
        HERROR( ERR_ARG, "(TRkMatrix) add", "argument has wrong dimensions" );

    //
    // In general, when adding two rank-k-matrices, we get a rank-2k-matrix.
    // we do not truncate the rank here (just in case, you want something
    // different), so you have to do this in some other place.
    //

    const B::Range  new_vecs( idx_t(_rank), idx_t(_rank + R->rank()) - 1 );
        
    if ( R->is_complex() )
        set_complex( true );
    
    set_rank( _rank + R->rank() );

    if ( R->is_complex() )
    {
        B::Matrix< complex >  A_new( _cmat_A, B::Range( 0, idx_t(_rows)-1 ), new_vecs );
        B::Matrix< complex >  B_new( _cmat_B, B::Range( 0, idx_t(_cols)-1 ), new_vecs );
        
        // copy new vectors
        B::copy( R->blas_cmat_A(), A_new );
        B::copy( R->blas_cmat_B(), B_new );
        
        // scale vectors (choose smaller ones)
        if ( alpha != real(1) )
        {
            if (_rows < _cols) B::scale( complex(alpha), A_new );
            else               B::scale( complex(alpha), B_new );
        }// if
    }// if
    else
    {
        if ( is_complex() )
            HERROR( ERR_REAL_CMPLX, "(TRkMatrix) add", "" );
        else
        {
            B::Matrix< real >  A_new( _rmat_A, B::Range( 0, idx_t(_rows)-1 ), new_vecs );
            B::Matrix< real >  B_new( _rmat_B, B::Range( 0, idx_t(_cols)-1 ), new_vecs );
            
            // copy new vectors
            B::copy( R->blas_rmat_A(), A_new );
            B::copy( R->blas_rmat_B(), B_new );
        
            // scale vectors (choose smaller ones)
            if ( alpha != real(1) )
            {
                if (_rows < _cols) B::scale( alpha, A_new );
                else               B::scale( alpha, B_new );
            }// if
        }// else
    }// else
}

//
// matrix-matrix-multiplication
//
TRkMatrix *
TRkMatrix::mul_right ( const real       alpha,
                       const TMatrix *  B,
                       const matop_t    op_A,
                       const matop_t    op_B ) const
{
    const bool  trans_A = (op_A != apply_normal);
    const bool  trans_B = (op_B != apply_normal);
    
    if ( B == nullptr )
        HERROR( ERR_ARG, "(TRkMatrix) mul_right", "argument is nullptr" );

    if ( (trans_A ? rows() : cols()) != (trans_B ? B->cols() : B->rows()) )
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) mul_right", "incompatible dimensions" );

    if ( is_complex() != B->is_complex() )
        HERROR( ERR_NOT_IMPL, "(TRkMatrix) mul_right", "complex * real" );

    ////////////////////////////////////////////////////////////////
    //
    // build matrix holding the result
    //

    const TRkMatrix *  A = this;
    const size_t       n = ( trans_A ? A->cols() : A->rows() );
    const size_t       m = ( trans_B ? B->rows() : B->cols() );
    const size_t       k = A->rank();
    auto               C = make_unique< TRkMatrix >();

    C->set_complex( A->is_complex() || B->is_complex() );
    C->set_size( n, m, k );
    C->set_ofs( (trans_A ? A->col_ofs() : A->row_ofs()),
                (trans_B ? B->row_ofs() : B->col_ofs()) );

    if ( k == 0 )
        return C.release();

    ////////////////////////////////////////////////////////////////
    //
    // multiply
    //

    //
    // depending on op_A, copy either A- or B-vectors to A(C)
    //
    
    if ( is_complex() )
    {
        if ( op_A == apply_normal )
            B::copy( A->blas_cmat_A(), C->blas_cmat_A() );
        else if ( op_A == apply_transposed )
        {
            B::copy( A->blas_cmat_B(), C->blas_cmat_A() );
            B::conj( C->blas_cmat_A() );
        }// if
        else
            B::copy( A->blas_cmat_B(), C->blas_cmat_A() );
    }// if
    else
    {
        if ( trans_A ) B::copy( A->blas_rmat_B(), C->blas_rmat_A() );
        else           B::copy( A->blas_rmat_A(), C->blas_rmat_A() );
    }// else

    //
    // compute B-vectors of C
    //
    
    const size_t  ma   = ( trans_A ? A->rows()    : A->cols()    );
    const idx_t   xofs = ( trans_B ? B->col_ofs() : B->row_ofs() );

    if ( is_complex() )
    {
        if ( op_A == apply_normal )
        {
            if ( op_B == apply_normal )
            {
                // B(C)^H = B(A)^H * B => B(C) = B^H * B(A)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_B( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_adjoint );
                }// for
            }// if
            else if ( op_B == apply_transposed )
            {
                TScalarVector  sx( ma, xofs, true );
            
                // B(C)^H = B(A)^H * B^T => B(C) = conj(B) * B(A) = conj( B * conj(B(A)) )
                for ( uint i = 0; i < k; i++ )
                {
                    B::Vector< complex >  B_i( A->blas_cmat_B().column( i ) );
                    TScalarVector         vy = C->vec_B( i );
                    
                    B::copy( B_i, sx.blas_cvec() );
                    B::conj( sx.blas_cvec() );
                    B->mul_vec( real(1), & sx, real(0), & vy, apply_normal );
                    B::conj( vy.blas_cvec() );
                }// for
            }// if
            else if ( op_B == apply_adjoint )
            {
                // B(C)^H = B(A)^H * B^H => B(C) = B * B(A)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_B( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_normal );
                }// for
            }// if
        }// if
        else if ( op_A == apply_transposed )
        {
            if ( op_B == apply_normal )
            {
                // B(C)^H = A(A)^T * B => B(C) = B^H * conj(A(A)) = conj(B^T * A)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_transposed );
                }// for

                B::conj( C->blas_cmat_B() );
            }// if
            else if ( op_B == apply_transposed )
            {
                // B(C)^H = A(A)^T * B^T => B(C) = conj(B) * conj(A(A)) = conj(B * A(A))
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_normal );
                }// for

                B::conj( C->blas_cmat_B() );
            }// if
            else if ( op_B == apply_adjoint )
            {
                TScalarVector  sx( ma, xofs, true );
                
                // B(C)^H = A(A)^T * B^H => B(C) = B * conj(A(A))
                for ( uint i = 0; i < k; i++ )
                {
                    B::Vector< complex >  A_i( A->blas_cmat_A().column( i ) );
                    TScalarVector         vy = C->vec_B( i );
                    
                    B::copy( A_i, sx.blas_cvec() );
                    B::conj( sx.blas_cvec() );
                    B->mul_vec( real(1), & sx, real(0), & vy, apply_normal );
                }// for
            }// if
        }// if
        else if ( op_A == apply_adjoint )
        {
            if ( op_B == apply_normal )
            {
                // B(C)^H = A(A)^H * B => B(C) = B^H * A(A)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_adjoint );
                }// for
            }// if
            else if ( op_B == apply_transposed )
            {
                TScalarVector  sx( ma, xofs, true );
                
                // B(C)^H = A(A)^H * B^T => B(C) = conj(B) * A(A) = conj(B * conj(A(A)))
                for ( uint i = 0; i < k; i++ )
                {
                    B::Vector< complex >  A_i( A->blas_cmat_A().column( i ) );
                    TScalarVector         vy = C->vec_B( i );
                    
                    B::copy( A_i, sx.blas_cvec() );
                    B::conj( sx.blas_cvec() );
                    B->mul_vec( real(1), & sx, real(0), & vy, apply_normal );
                    B::conj( vy.blas_cvec() );
                }// for
            }// if
            else if ( op_B == apply_adjoint )
            {
                // B(C)^H = A(A)^H * B^H => B(C) = B * A(A)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_normal );
                }// for
            }// if
        }// if
    }// if
    else
    {
        if ( op_A == apply_normal )
        {
            if ( op_B == apply_normal )
            {
                // compute B(A)^T * B as B^T B(A) for each vector b_i
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_B( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_transposed );
                }// for
            }// if
            else
            {
                // compute B(A)^T * B^T as B B(A) for each vector b_i
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_B( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_normal );
                }// for
            }// else
        }// if
        else
        {
            if ( op_B == apply_normal )
            {
                // compute A(B)^T * B as B^T A(B)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_transposed );
                }// for
            }// if
            else
            {
                // compute A(B)^T * B^T as B A(B)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_normal );
                }// for
            }// else
        }// else
    }// else

    //
    // finally, scale C by alpha
    //
    
    if ( alpha != real(1) )
    {
        if ( C->is_complex() )
        {
            if ( n > m ) B::scale( complex(alpha), C->blas_cmat_B() );
            else         B::scale( complex(alpha), C->blas_cmat_A() );
        }// if
        else
        {
            if ( n > m ) B::scale( alpha, C->blas_rmat_B() );
            else         B::scale( alpha, C->blas_rmat_A() );
        }// else
    }// if

    return C.release();
}

TRkMatrix *
TRkMatrix::mul_left ( const real       alpha,
                      const TMatrix *  A,
                      const matop_t    op_A,
                      const matop_t    op_B ) const
{
    const bool  trans_A = (op_A != apply_normal);
    const bool  trans_B = (op_B != apply_normal);
    
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TRkMatrix) mul_left", "argument is nullptr" );

    if ( (trans_A ? A->rows() : A->cols()) != (trans_B ? cols() : rows()) )
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) mul_left", "incompatible dimensions" );

    if ( is_complex() != A->is_complex() )
        HERROR( ERR_NOT_IMPL, "(TRkMatrix) mul_left", "complex * real" );

    ////////////////////////////////////////////////////////////////
    //
    // build temporary matrix holding the result
    //

    const TRkMatrix *  B = this;
    const size_t       n = ( trans_A ? A->cols() : A->rows() );
    const size_t       m = ( trans_B ? B->rows() : B->cols() );
    const size_t       k = B->rank();
    auto               C = make_unique< TRkMatrix >();

    C->set_complex( A->is_complex() || B->is_complex() );
    C->set_size( n, m, k );
    C->set_ofs( (trans_A ? A->col_ofs() : A->row_ofs()),
                (trans_B ? B->row_ofs() : B->col_ofs()) );
    
    if ( k == 0 )
        return C.release();
    
    ////////////////////////////////////////////////////////////////
    //
    // multiply
    //

    //
    // copy B-/A-vectors
    //
    
    if ( is_complex() )
    {
        if ( op_B == apply_normal )
            B::copy( B->blas_cmat_B(), C->blas_cmat_B() );
        else if ( op_B == apply_transposed )
        {
            B::copy( B->blas_cmat_A(), C->blas_cmat_B() );
            B::conj( C->blas_cmat_B() );
        }// if
        else if ( op_B == apply_adjoint )
            B::copy( B->blas_cmat_A(), C->blas_cmat_B() );
    }// if
    else
    {
        if ( op_B == apply_normal )
            B::copy( B->blas_rmat_B(), C->blas_rmat_B() );
        else
            B::copy( B->blas_rmat_A(), C->blas_rmat_B() );
    }// else

    //
    // compute A- or B-vectors of C
    //
    
    const size_t    nb   = ( trans_B ? B->cols()    : B->rows()    );
    const idx_t     xofs = ( trans_A ? A->row_ofs() : A->col_ofs() );

    if ( is_complex() )
    {
        if ( op_B == apply_normal )
        {
            // A(C) = op(A) * A(B)
            for ( uint i = 0; i < k; i++ )
            {
                TScalarVector  vx = B->vec_A( i );
                TScalarVector  vy = C->vec_A( i );
                
                A->mul_vec( real(1), & vx, real(0), & vy, op_A );
            }// for
        }// if
        else if ( op_B == apply_transposed )
        {
            TScalarVector  sx( nb, xofs, true );
            
            // A(C) = op(A) * conj(B(B))
            for ( uint i = 0; i < k; i++ )
            {
                B::Vector< complex >  B_i( B->blas_cmat_B().column( i ) );
                TScalarVector         vy = C->vec_A( i );
                
                B::copy( B_i, sx.blas_cvec() );
                B::conj( sx.blas_cvec() );
                A->mul_vec( real(1), & sx, real(0), & vy, op_A );
            }// for
        }// if
        else if ( op_B == apply_adjoint )
        {
            // A(C) = op(A) * B(B)
            for ( uint i = 0; i < k; i++ )
            {
                TScalarVector  vx = B->vec_B( i );
                TScalarVector  vy = C->vec_A( i );
                
                A->mul_vec( real(1), & vx, real(0), & vy, op_A );
            }// for
        }// if
    }// if
    else
    {
        if ( op_B == apply_normal )
        {
            // multiply A * A(B)
            for ( uint i = 0; i < k; i++ )
            {
                TScalarVector  vx = B->vec_A( i );
                TScalarVector  vy = C->vec_A( i );
                
                A->mul_vec( real(1), & vx, real(0), & vy, op_A );
            }// for
        }// if
        else
        {
            // multiply A * B(B)^T
            for ( uint i = 0; i < k; i++ )
            {
                TScalarVector  vx = B->vec_B( i );
                TScalarVector  vy = C->vec_A( i );
                
                A->mul_vec( real(1), & vx, real(0), & vy, op_A );
            }// for
        }// else
    }// else

    //
    // finally, scale C with alpha
    //
    
    if ( alpha != real(1) )
    {
        if ( C->is_complex() )
        {
            if ( n > m ) B::scale( complex(alpha), C->blas_cmat_B() );
            else         B::scale( complex(alpha), C->blas_cmat_A() );
        }// if
        else
        {
            if ( n > m ) B::scale( alpha, C->blas_rmat_B() );
            else         B::scale( alpha, C->blas_rmat_A() );
        }// else
    }// if

    return C.release();
}
    
/////////////////////////////////////////////////
//
// BLAS-routines (complex valued)
//

//
// scale matrix by constant factor
//
void
TRkMatrix::cscale ( const complex f )
{
    if (( f == complex(1) ) || ( rank() == 0 ))
        return;

    if ( std::imag( f ) != real(0) )
        set_complex( true );

    TScopedLock  mlock( *this );

    if ( is_complex() )
    {
        if ( f == complex(0) )
        {
            B::fill( complex(0), blas_cmat_A() );
            B::fill( complex(0), blas_cmat_B() );
        }// if
        else
        {
            //
            // scale matrix with lower norm
            //

            if ( B::normF( blas_cmat_A() ) < B::normF( blas_cmat_B() ) )
                B::scale( f, blas_cmat_A() );
            else
                B::scale( f, blas_cmat_B() );
        }// else
    }// if
    else
    {
        if ( f == complex(0) )
        {
            B::fill( real(0), blas_rmat_A() );
            B::fill( real(0), blas_rmat_B() );
        }// if
        else
        {
            //
            // scale matrix with lower norm
            //

            if ( B::normF( blas_rmat_A() ) < B::normF( blas_rmat_B() ) )
                B::scale( std::real( f ), blas_rmat_A() );
            else
                B::scale( std::real( f ), blas_rmat_B() );
        }// else
    }// else
}
    
//
// matrix-vector-multiplication
//
void
TRkMatrix::cmul_vec ( const complex   alpha,
                      const TVector * x,
                      const complex   beta,
                      TVector       * y,
                      const matop_t   op ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TRkMatrix) cmul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TRkMatrix) cmul_vec", "y = nullptr" );
    
    if ( op == apply_normal )
    {
        if (( row_is() != y->is() ) || ( col_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TRkMatrix) cmul_vec", "incompatible vector index set" );
    }// if
    else
    {
        if (( col_is() != y->is() ) || ( row_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TRkMatrix) cmul_vec", "incompatible vector index set" );
    }// if
    
    //
    // check if scalar vector
    //
    
    if ( IS_TYPE( x, TScalarVector ) && IS_TYPE( y, TScalarVector ) )
    {
        TScalarVector       * sy = ptrcast( y, TScalarVector );
        const TScalarVector * sx = cptrcast( x, TScalarVector );

        //
        // we assume, that the index-sets are coherent
        // and we can access by local indices
        //

        // scale left-side
        if      ( beta == complex(0) ) sy->fill( real(0) );
        else if ( beta != complex(1) ) sy->cscale( beta );

        if (( alpha == complex(0) ) || ( _rank == 0 ))
            return;

        if ( is_complex() || x->is_complex() || ( std::imag(alpha) != real(0) ))
            y->set_complex( true );

        if ( is_complex() )
        {
            if ( ! x->is_complex() )
                HERROR( ERR_REAL_CMPLX, "(TRkMatrix) cmul_vec", "" );
                
            task_mulvec< complex >( alpha, op, blas_cmat_A(), blas_cmat_B(), sx->blas_cvec(),
                                    sy->blas_cvec() );
        }// if
        else
        {
            if ( x->is_complex() || ( std::imag( alpha ) != real(0) ))
                HERROR( ERR_REAL_CMPLX, "(TRkMatrix) cmul_vec", "" );

            task_mulvec< real >( std::real( alpha ), op, blas_rmat_A(), blas_rmat_B(), sx->blas_rvec(),
                                 sy->blas_rvec() );
        }// else
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TRkMatrix) cmul_vec",
                y->typestr() + " += TRkMatrix * " + x->typestr() );
}

//
// compute this = this + a * matrix
//
void
TRkMatrix::cadd ( const complex alpha, const TMatrix * M )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "(TRkMatrix) add", "nullptr argument" );
    
    if ( ! IS_TYPE( M, TRkMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TRkMatrix) add", M->typestr() );

    const TRkMatrix * R = cptrcast( M, TRkMatrix );

    //
    // some trivial checks
    //

    if (( alpha == complex(0) ) || ( R->rank() == 0 ))
        return;

    if ((rows() != R->rows()) || (cols() != R->cols()))
        HERROR( ERR_ARG, "(TRkMatrix) cadd", "argument has wrong dimensions" );

    //
    // In general, when adding two rank-k-matrices, we get a rank-2k-matrix.
    // we do not truncate the rank here (just in case, you want something
    // different), so you have to do this in some other place.
    //

    const B::Range  new_vecs( idx_t(_rank), idx_t(_rank + R->rank()) - 1 );

    if ( R->is_complex() || ( std::imag( alpha ) != real(0) ))
        set_complex( true );
    
    set_rank( _rank + R->rank() );

    if ( R->is_complex() )
    {
        B::Matrix< complex >  A_new( _cmat_A, B::Range( 0, idx_t(_rows)-1 ), new_vecs );
        B::Matrix< complex >  B_new( _cmat_B, B::Range( 0, idx_t(_cols)-1 ), new_vecs );
        
        // copy new vectors
        B::copy( R->blas_cmat_A(), A_new );
        B::copy( R->blas_cmat_B(), B_new );
        
        // scale vectors (choose smaller ones)
        if ( alpha != complex(1) )
        {
            if (_rows < _cols) B::scale( alpha, A_new );
            else               B::scale( alpha, B_new );
        }// if
    }// if
    else
    {
        if ( is_complex() || ( std::imag( alpha ) != real(0) ))
            HERROR( ERR_REAL_CMPLX, "(TRkMatrix) cadd", "" );
        else
        {
            B::Matrix< real >  A_new( _rmat_A, B::Range( 0, idx_t(_rows)-1 ), new_vecs );
            B::Matrix< real >  B_new( _rmat_B, B::Range( 0, idx_t(_cols)-1 ), new_vecs );
            
            // copy new vectors
            B::copy( R->blas_rmat_A(), A_new );
            B::copy( R->blas_rmat_B(), B_new );
        
            // scale vectors (choose smaller ones)
            if ( alpha != complex(1) )
            {
                if (_rows < _cols) B::scale( std::real(alpha), A_new );
                else               B::scale( std::real(alpha), B_new );
            }// if
        }// else
    }// else
}

//
// matrix-matrix-multiplication
//
TRkMatrix *
TRkMatrix::cmul_right ( const complex alpha, const TMatrix * B,
                        const matop_t op_A, const matop_t op_B ) const
{
    const bool  trans_A = (op_A != apply_normal);
    const bool  trans_B = (op_B != apply_normal);
    
    if ( B == nullptr )
        HERROR( ERR_ARG, "(TRkMatrix) cmul_right", "argument is nullptr" );

    if ( (trans_A ? rows() : cols()) != (trans_B ? B->cols() : B->rows()) )
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) cmul_right", "incompatible dimensions" );
    
    if ( is_complex() != B->is_complex() )
        HERROR( ERR_NOT_IMPL, "(TRkMatrix) cmul_right", "complex * real" );

    ////////////////////////////////////////////////////////////////
    //
    // build temporary matrix holding the result
    //

    const TRkMatrix *  A = this;
    const size_t       n = ( trans_A ? A->cols() : A->rows() );
    const size_t       m = ( trans_B ? B->rows() : B->cols() );
    const size_t       k = A->rank();
    auto               C = make_unique< TRkMatrix >();

    C->set_complex( A->is_complex() || B->is_complex());
    C->set_size( n, m, k );
    C->set_ofs( (trans_A ? A->col_ofs() : A->row_ofs()),
                (trans_B ? B->row_ofs() : B->col_ofs()) );

    if ( k == 0 )
        return C.release();
    
    ////////////////////////////////////////////////////////////////
    //
    // multiply
    //

    //
    // depending on op_A, copy either A- or B-vectors to A(C)
    //
    
    if ( is_complex() )
    {
        if ( op_A == apply_normal )
            B::copy( A->blas_cmat_A(), C->blas_cmat_A() );
        else if ( op_A == apply_transposed )
        {
            B::copy( A->blas_cmat_B(), C->blas_cmat_A() );
            B::conj( C->blas_cmat_A() );
        }// if
        else
            B::copy( A->blas_cmat_B(), C->blas_cmat_A() );
    }// if
    else
    {
        if ( trans_A ) B::copy( A->blas_rmat_B(), C->blas_rmat_A() );
        else           B::copy( A->blas_rmat_A(), C->blas_rmat_A() );
    }// else

    //
    // compute B-vectors of C
    //
    
    const size_t    ma   = ( trans_A ? A->rows()    : A->cols()    );
    const idx_t     xofs = ( trans_B ? B->col_ofs() : B->row_ofs() );
    
    if ( is_complex() )
    {
        if ( op_A == apply_normal )
        {
            if ( op_B == apply_normal )
            {
                // B(C)^H = B(A)^H * B => B(C) = B^H * B(A)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_B( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_adjoint );
                }// for
            }// if
            else if ( op_B == apply_transposed )
            {
                TScalarVector  sx( ma, xofs, true );
            
                // B(C)^H = B(A)^H * B^T => B(C) = conj(B) * B(A) = conj( B * conj(B(A)) )
                for ( uint i = 0; i < k; i++ )
                {
                    B::Vector< complex >  B_i( A->blas_cmat_B().column( i ) );
                    TScalarVector         vy = C->vec_B( i );
                    
                    B::copy( B_i, sx.blas_cvec() );
                    B::conj( sx.blas_cvec() );
                    B->mul_vec( real(1), & sx, real(0), & vy, apply_normal );
                    B::conj( vy.blas_cvec() );
                }// for
            }// if
            else if ( op_B == apply_adjoint )
            {
                // B(C)^H = B(A)^H * B^H => B(C) = B * B(A)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_B( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_normal );
                }// for
            }// if
        }// if
        else if ( op_A == apply_transposed )
        {
            if ( op_B == apply_normal )
            {
                // B(C)^H = A(A)^T * B => B(C) = B^H * conj(A(A)) = conj(B^T * A)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_transposed );
                }// for

                B::conj( C->blas_cmat_B() );
            }// if
            else if ( op_B == apply_transposed )
            {
                // B(C)^H = A(A)^T * B^T => B(C) = conj(B) * conj(A(A)) = conj(B * A(A))
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_normal );
                }// for

                B::conj( C->blas_cmat_B() );
            }// if
            else if ( op_B == apply_adjoint )
            {
                TScalarVector  sx( ma, xofs, true );
                
                // B(C)^H = A(A)^T * B^H => B(C) = B * conj(A(A))
                for ( uint i = 0; i < k; i++ )
                {
                    B::Vector< complex >  A_i( A->blas_cmat_A().column( i ) );
                    TScalarVector         vy = C->vec_B( i );
                    
                    B::copy( A_i, sx.blas_cvec() );
                    B::conj( sx.blas_cvec() );
                    B->mul_vec( real(1), & sx, real(0), & vy, apply_normal );
                }// for
            }// if
        }// if
        else if ( op_A == apply_adjoint )
        {
            if ( op_B == apply_normal )
            {
                // B(C)^H = A(A)^H * B => B(C) = B^H * A(A)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_adjoint );
                }// for
            }// if
            else if ( op_B == apply_transposed )
            {
                TScalarVector  sx( ma, xofs, true );
                
                // B(C)^H = A(A)^H * B^T => B(C) = conj(B) * A(A) = conj(B * conj(A(A)))
                for ( uint i = 0; i < k; i++ )
                {
                    B::Vector< complex >  A_i( A->blas_cmat_A().column( i ) );
                    TScalarVector         vy = C->vec_B( i );
                    
                    B::copy( A_i, sx.blas_cvec() );
                    B::conj( sx.blas_cvec() );
                    B->mul_vec( real(1), & sx, real(0), & vy, apply_normal );
                    B::conj( vy.blas_cvec() );
                }// for
            }// if
            else if ( op_B == apply_adjoint )
            {
                // B(C)^H = A(A)^H * B^H => B(C) = B * A(A)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_normal );
                }// for
            }// if
        }// if
    }// if
    else
    {
        if ( op_A == apply_normal )
        {
            if ( op_B == apply_normal )
            {
                // compute B(A)^T * B as B^T B(A) for each vector b_i
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_B( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_transposed );
                }// for
            }// if
            else
            {
                // compute B(A)^T * B^T as B B(A) for each vector b_i
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_B( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_normal );
                }// for
            }// else
        }// if
        else
        {
            if ( op_B == apply_normal )
            {
                // compute A(B)^T * B as B^T A(B)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_transposed );
                }// for
            }// if
            else
            {
                // compute A(B)^T * B^T as B A(B)
                for ( uint i = 0; i < k; i++ )
                {
                    TScalarVector  vx = A->vec_A( i );
                    TScalarVector  vy = C->vec_B( i );
                    
                    B->mul_vec( real(1), & vx, real(0), & vy, apply_normal );
                }// for
            }// else
        }// else
    }// else

    //
    // finally, scale C by alpha
    //
    
    if ( alpha != complex(1) )
    {
        // adjust field type of C in case of complex alpha
        // TODO: create C with correct field type and multiply accordingly
        if ( std::imag( alpha ) != real(0) )
            C->set_complex( true );
        
        if ( C->is_complex() )
        {
            if ( n > m ) B::scale( conj(alpha), C->blas_cmat_B() );
            else         B::scale( alpha,       C->blas_cmat_A() );
        }// if
        else
        {
            if ( n > m ) B::scale( std::real( alpha ), C->blas_rmat_B() );
            else         B::scale( std::real( alpha ), C->blas_rmat_A() );
        }// else
    }// if

    return C.release();
}

TRkMatrix *
TRkMatrix::cmul_left ( const complex alpha, const TMatrix * A,
                       const matop_t op_A, const matop_t op_B ) const
{
    const bool  trans_A = (op_A != apply_normal);
    const bool  trans_B = (op_B != apply_normal);
    
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TRkMatrix) cmul_left", "argument is nullptr" );

    if ( (trans_A ? A->rows() : A->cols()) != (trans_B ? cols() : rows()) )
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) cmul_left", "incompatible dimensions" );

    if ( is_complex() != A->is_complex() )
        HERROR( ERR_NOT_IMPL, "(TRkMatrix) cmul_left", "complex * real" );

    ////////////////////////////////////////////////////////////////
    //
    // build temporary matrix holding the result
    //

    const TRkMatrix *  B = this;
    const size_t       n = ( trans_A ? A->cols() : A->rows() );
    const size_t       m = ( trans_B ? B->rows() : B->cols() );
    const size_t       k = B->rank();
    auto               C = make_unique< TRkMatrix >();

    C->set_complex( A->is_complex() || B->is_complex() );
    C->set_size( n, m, k );
    C->set_ofs( (trans_A ? A->col_ofs() : A->row_ofs()),
                (trans_B ? B->row_ofs() : B->col_ofs()) );
    
    if ( k == 0 )
        return C.release();
    
    ////////////////////////////////////////////////////////////////
    //
    // multiply
    //

    //
    // copy B-/A-vectors
    //
    
    if ( is_complex() )
    {
        if ( op_B == apply_normal )
            B::copy( B->blas_cmat_B(), C->blas_cmat_B() );
        else if ( op_B == apply_transposed )
        {
            B::copy( B->blas_cmat_A(), C->blas_cmat_B() );
            B::conj( C->blas_cmat_B() );
        }// if
        else if ( op_B == apply_adjoint )
            B::copy( B->blas_cmat_A(), C->blas_cmat_B() );
    }// if
    else
    {
        if ( op_B == apply_normal )
            B::copy( B->blas_rmat_B(), C->blas_rmat_B() );
        else
            B::copy( B->blas_rmat_A(), C->blas_rmat_B() );
    }// else

    //
    // compute A- or B-vectors of C
    //
    
    const size_t  nb   = ( trans_B ? B->cols()    : B->rows()    );
    const idx_t   xofs = ( trans_A ? A->row_ofs() : A->col_ofs() );

    if ( is_complex() )
    {
        if ( op_B == apply_normal )
        {
            // multiply A * A(B)
            for ( uint i = 0; i < k; i++ )
            {
                TScalarVector  vx = B->vec_A( i );
                TScalarVector  vy = C->vec_A( i );
                
                A->mul_vec( real(1), & vx, real(0), & vy, op_A );
            }// for
        }// if
        else if ( op_B == apply_transposed )
        {
            TScalarVector  sx( nb, xofs, true );
            
            // multiply A * conj(B(B))
            for ( uint i = 0; i < k; i++ )
            {
                B::Vector< complex >  B_i( B->blas_cmat_B().column( i ) );
                TScalarVector         vy = C->vec_A( i );
                
                B::copy( B_i, sx.blas_cvec() );
                B::conj( sx.blas_cvec() );
                A->mul_vec( real(1), & sx, real(0), & vy, op_A );
            }// for
        }// if
        else if ( op_B == apply_adjoint )
        {
            // multiply A * B(B)
            for ( uint i = 0; i < k; i++ )
            {
                TScalarVector  vx = B->vec_B( i );
                TScalarVector  vy = C->vec_A( i );
                
                A->mul_vec( real(1), & vx, real(0), & vy, op_A );
            }// for
        }// if
    }// if
    else
    {
        if ( op_B == apply_normal )
        {
            // multiply A * A(B)
            for ( uint i = 0; i < k; i++ )
            {
                TScalarVector  vx = B->vec_A( i );
                TScalarVector  vy = C->vec_A( i );
                
                A->mul_vec( real(1), & vx, real(0), & vy, op_A );
            }// for
        }// if
        else
        {
            // multiply A * B(B)^T
            for ( uint i = 0; i < k; i++ )
            {
                TScalarVector  vx = B->vec_B( i );
                TScalarVector  vy = C->vec_A( i );
                
                A->mul_vec( real(1), & vx, real(0), & vy, op_A );
            }// for
        }// else
    }// else

    //
    // finally, scale C with alpha
    //
    
    if ( alpha != complex(1) )
    {
        // adjust field type of C in case of complex alpha
        // TODO: create C with correct field type and multiply accordingly
        if ( std::imag( alpha ) != real(0) )
            C->set_complex( true );
        
        if ( C->is_complex() )
        {
            if ( n > m ) B::scale( conj(alpha), C->blas_cmat_B() );
            else         B::scale( alpha,       C->blas_cmat_A() );
        }// if
        else
        {
            if ( n > m ) B::scale( std::real( alpha ), C->blas_rmat_B() );
            else         B::scale( std::real( alpha ), C->blas_rmat_A() );
        }// else
    }// if

    return C.release();
}
    
///////////////////////////////////////////////////////////
//
// linear operator mapping
//

void
TRkMatrix::apply_add   ( const real                       alpha,
                         const BLAS::Vector< real > &     x,
                         BLAS::Vector< real > &           y,
                         const matop_t                    op ) const
{
    HASSERT( ( x.length() == ncols( op ) ) && ( y.length() == nrows( op ) ),
             ERR_ARG, "(TRkMatrix) apply_add", "incompatible vector dimensions" );
    HASSERT( is_real(), 
             ERR_REAL_CMPLX, "(TRkMatrix) apply_add", "real vectors, complex matrix" );

    switch ( op )
    {
        case apply_normal :
        {
            const auto  t = BLAS::mulvec( alpha, BLAS::adjoint( blas_mat_B< real >( this ) ), x );
                
            BLAS::mulvec( real(1), blas_mat_A< real >( this ), t, real(1), y );
        }
        break;

        case apply_transposed :
        {
            const auto  t = BLAS::mulvec( alpha, BLAS::transposed( blas_mat_A< real >( this ) ), x );
                
            BLAS::mulvec( real(1), blas_mat_B< real >( this ), t, real(1), y );
        }
        break;

        case apply_adjoint :
        {
            const auto  t = BLAS::mulvec( alpha, BLAS::adjoint( blas_mat_A< real >( this ) ), x );
                
            BLAS::mulvec( real(1), blas_mat_B< real >( this ), t, real(1), y );
        }
        break;

        default:  // apply_conjugate
        {
            const auto  t = BLAS::mulvec( alpha, BLAS::adjoint( blas_mat_B< real >( this ) ), x );
                
            BLAS::mulvec( real(1), blas_mat_A< real >( this ), t, real(1), y );
        }
    }// switch
}

void
TRkMatrix::apply_add   ( const complex                    alpha,
                         const BLAS::Vector< complex > &  x,
                         BLAS::Vector< complex > &        y,
                         const matop_t                    op ) const
{
    HASSERT( ( x.length() == ncols( op ) ) && ( y.length() == nrows( op ) ),
             ERR_ARG, "(TRkMatrix) apply_add", "incompatible vector dimensions" );
    HASSERT( is_complex(), 
             ERR_REAL_CMPLX, "(TRkMatrix) apply_add", "complex vectors, real matrix" );

    switch ( op )
    {
        case apply_normal :
        {
            const auto  t = BLAS::mulvec( alpha, BLAS::adjoint( blas_mat_B< complex >( this ) ), x );
                
            BLAS::mulvec( complex(1), blas_mat_A< complex >( this ), t, complex(1), y );
        }
        break;

        case apply_transposed :
        {
            const auto  t = BLAS::mulvec( complex(1), BLAS::transposed( blas_mat_A< complex >( this ) ), x );
            auto        s = BLAS::mulvec( complex(1), blas_mat_B< complex >( this ), t );
                
            BLAS::conj( s );
            BLAS::add( alpha, s, y );
        }
        break;

        case apply_adjoint :
        {
            const auto  t = BLAS::mulvec( alpha, BLAS::adjoint( blas_mat_A< complex >( this ) ), x );
                
            BLAS::mulvec( complex(1), blas_mat_B< complex >( this ), t, complex(1), y );
        }
        break;

        default:  // apply_conjugate
        {
            const auto  t = BLAS::mulvec( complex(1), BLAS::transposed( blas_mat_B< complex >( this ) ), x );
            auto        s = BLAS::mulvec( complex(1), blas_mat_A< complex >( this ), t );
                
            BLAS::conj( s );
            BLAS::add( alpha, s, y );
        }
    }// switch
}

/////////////////////////////////////////////////
//
// misc methods
//

//
// transpose matrix
//
void
TRkMatrix::transpose ()
{
    //
    // (A·B^H)^T = B^H^T A^T = conj(B) conj(A)^H
    //
    
    conjugate();

    {
        TScopedLock  mlock( *this );

        if ( is_complex() )
        {
            std::swap( _cmat_A, _cmat_B );
        }// if
        else
        {
            std::swap( _rmat_A, _rmat_B );
        }// else
    
        std::swap( _rows, _cols );
    }

    TMatrix::transpose();
}
    
//
// conjugate matrix coefficients
//
void
TRkMatrix::conjugate ()
{
    if ( is_complex() && ! is_hermitian() )
    {
        TScopedLock  mlock( *this );
        
        B::conj( blas_cmat_A() );
        B::conj( blas_cmat_B() );
    }// if
}
    
//
// allocate memory (and update data)
//
void
TRkMatrix::set_size ( const size_t  n,
                      const size_t  m,
                      const size_t  new_rank )
{
    //
    // first check if we really have to do something
    //
    if (( _rows == n ) && ( _cols == m ) && ( new_rank == _rank ))
        return;

    TScopedLock  mlock( *this );

    //
    // a rank 0 means -> remove all data
    //

    if ( new_rank == 0 )
    {
        // delete old data
        if ( is_complex() )
        {
            _cmat_A = B::Matrix< complex >( n, 0 );
            _cmat_B = B::Matrix< complex >( m, 0 );
        }// if
        else
        {
            _rmat_A = B::Matrix< real >( n, 0 );
            _rmat_B = B::Matrix< real >( m, 0 );
        }// else

        _rows = n;
        _cols = m;
        _rank = 0;
    }// if
    else
    {
        const  B::Range  loc_row_is( 0, idx_t(n)-1 );
        const  B::Range  loc_col_is( 0, idx_t(m)-1 );
        const  B::Range  old_vec_is( 0, idx_t( std::min( _rank, new_rank ) ) - 1 );
            
        //
        // handle matrix A
        //

        if ( is_complex() )
        {
            if ( _rows == n )
            {
                // only copy if different rank
                if ( _rank != new_rank )
                {
                    B::Matrix< complex >  new_A( n, new_rank );
                    B::Matrix< complex >  sub_old_A( _cmat_A, loc_row_is, old_vec_is );
                    B::Matrix< complex >  sub_new_A( new_A, loc_row_is, old_vec_is );

                    B::copy( sub_old_A, sub_new_A );

                    _cmat_A = std::move( new_A );
                }// if
            }// if
            else
            {
                // simply create (nullified) array
                _rows   = n;
                _cmat_A = B::Matrix< complex >( n, new_rank );
            }// else
        }// if
        else
        {
            if ( _rows == n )
            {
                // only copy if different rank
                if ( _rank != new_rank )
                {
                    B::Matrix< real >  new_A( n, new_rank );
                    B::Matrix< real >  sub_old_A( _rmat_A, loc_row_is, old_vec_is );
                    B::Matrix< real >  sub_new_A( new_A, loc_row_is, old_vec_is );

                    B::copy( sub_old_A, sub_new_A );
                
                    _rmat_A = std::move( new_A );
                }// if
            }// if
            else
            {
                // simply create (nullified) array
                _rows   = n;
                _rmat_A = B::Matrix< real >( n, new_rank );
            }// else
        }// else

        //
        // handle matrix B
        //

        if ( is_complex() )
        {
            if ( _cols == m )
            {
                // only copy if different rank
                if ( _rank != new_rank )
                {
                    B::Matrix< complex >  new_B( m, new_rank );
                    B::Matrix< complex >  sub_old_B( _cmat_B, loc_col_is, old_vec_is );
                    B::Matrix< complex >  sub_new_B( new_B, loc_col_is, old_vec_is );

                    B::copy( sub_old_B, sub_new_B );

                    _cmat_B = std::move( new_B );
                }// if
            }// if
            else
            {
                // simply create (nullified) arrays
                _cols   = m;
                _cmat_B = B::Matrix< complex >( m, new_rank );
            }// else
        }// if
        else
        {
            if ( _cols == m )
            {
                // only copy if different rank
                if ( _rank != new_rank )
                {
                    B::Matrix< real >  new_B( m, new_rank );
                    B::Matrix< real >  sub_old_B( _rmat_B, loc_col_is, old_vec_is );
                    B::Matrix< real >  sub_new_B( new_B, loc_col_is, old_vec_is );

                    B::copy( sub_old_B, sub_new_B );

                    _rmat_B = std::move( new_B );
                }// if
            }// if
            else
            {
                // simply create (nullified) arrays
                _cols   = m;
                _rmat_B = B::Matrix< real >( m, new_rank );
            }// else
        }// else

        // adjust rank
        _rank = new_rank;
    }// else

    if ( ! is_complex() )
    {
        if ( rows() != _rmat_A.nrows() )
            HERROR( ERR_CONSISTENCY, "", "" );

        if ( cols() != _rmat_B.nrows() )
            HERROR( ERR_CONSISTENCY, "", "" );
    }// if
}

//
// virtual constructor
//
std::unique_ptr< TMatrix >
TRkMatrix::copy () const
{
    auto  M = TMatrix::copy();
    auto  R = ptrcast( M.get(), TRkMatrix );

    if ( cluster() != nullptr )
        R->set_cluster( cluster() );

    R->set_rank( _rank );

    if ( _rank > 0 )
    {
        if ( is_complex() )
        {
            B::copy( blas_cmat_A(), R->blas_cmat_A() );
            B::copy( blas_cmat_B(), R->blas_cmat_B() );
        }// if
        else
        {
            B::copy( blas_rmat_A(), R->blas_rmat_A() );
            B::copy( blas_rmat_B(), R->blas_rmat_B() );
        }// else
    }// if

    return M;
}

//
// copy matrix wrt. given accuracy
//
std::unique_ptr< TMatrix >
TRkMatrix::copy ( const TTruncAcc & acc, const bool ) const
{
    auto  M = copy();
    auto  R = ptrcast( M.get(), TRkMatrix );

    R->truncate( acc( this ) );
    
    return M;
}

//
// return structural copy
//
std::unique_ptr< TMatrix >
TRkMatrix::copy_struct () const
{
    auto  M = TMatrix::copy_struct();
    auto  R = ptrcast( M.get(), TRkMatrix );

    if ( cluster() != nullptr )
        R->set_cluster( cluster() );

    R->set_rank( 0 );

    return M;
}

//
// copy matrix into A
//
void
TRkMatrix::copy_to ( TMatrix * A ) const
{
    TMatrix::copy_to( A );
    
    if ( ! IS_TYPE( A, TRkMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TRkMatrix) copy_to", A->typestr() );
    
    TRkMatrix * R = ptrcast( A, TRkMatrix );

    if (( R->rows() != rows() ) || ( R->cols() != cols() ))
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) copy_to", "matrix has wrong dimension" );

    R->set_complex( is_complex() );
    R->set_size( rows(), cols(), _rank );
    R->set_ofs( row_ofs(), col_ofs() );

    if ( _rank > 0 )
    {
        if ( is_complex() )
        {
            B::copy( blas_cmat_A(), R->blas_cmat_A() );
            B::copy( blas_cmat_B(), R->blas_cmat_B() );
        }// if
        else
        {
            B::copy( blas_rmat_A(), R->blas_rmat_A() );
            B::copy( blas_rmat_B(), R->blas_rmat_B() );
        }// else
    }// if
}

void
TRkMatrix::copy_to ( TMatrix * A, const TTruncAcc & acc, const bool ) const
{
    TMatrix::copy_to( A );
    
    if ( ! IS_TYPE( A, TRkMatrix ) )
    {
        HERROR( ERR_NOT_IMPL, "", "" );
        // convert( this, A, acc );
        // return;
    }// if
    
    TRkMatrix * R = ptrcast( A, TRkMatrix );

    if (( R->rows() != rows() ) || ( R->cols() != cols() ))
        HERROR( ERR_MAT_SIZE, "(TRkMatrix) copy_to", "matrix has wrong dimension" );

    R->set_complex( is_complex() );
    R->set_size( rows(), cols(), _rank );
    R->set_ofs( row_ofs(), col_ofs() );

    if ( _rank > 0 )
    {
        if ( is_complex() )
        {
            B::copy( blas_cmat_A(), R->blas_cmat_A() );
            B::copy( blas_cmat_B(), R->blas_cmat_B() );
        }// if
        else
        {
            B::copy( blas_rmat_A(), R->blas_rmat_A() );
            B::copy( blas_rmat_B(), R->blas_rmat_B() );
        }// else
        
        R->truncate( acc( this ) );
    }// if
}

//
// return size in bytes used by this object
//
size_t
TRkMatrix::byte_size () const
{
    size_t  s = TMatrix::byte_size() + 3 * sizeof(size_t);

    s += sizeof(B::Matrix<real>) + sizeof(B::Matrix<complex>);
    
    if ( is_complex() )
        s += (_rank * (_rows + _cols) * sizeof(complex));
    else
        s += (_rank * (_rows + _cols) * sizeof(real));

    return s;
}

//
// serialisation
//

void
TRkMatrix::read  ( TByteStream & s )
{
    TMatrix::read( s );

    size_t  n, m, k;
    
    s.get( n );
    s.get( m );
    s.get( k );

    set_size( n, m, k );

    if ( is_complex() )
    {
        s.get( _cmat_A.data(), sizeof( complex ) * _rows * _rank );
        s.get( _cmat_B.data(), sizeof( complex ) * _cols * _rank );
    }// if
    else
    {
        s.get( _rmat_A.data(), sizeof( real ) * _rows * _rank );
        s.get( _rmat_B.data(), sizeof( real ) * _cols * _rank );
    }// else
}

void
TRkMatrix::build  ( TByteStream & s )
{
    TMatrix::build( s );

    size_t  n, m, k;
    
    s.get( n );
    s.get( m );
    s.get( k );

    set_size( n, m, k );

    if ( is_complex() )
    {
        s.get( _cmat_A.data(), sizeof( complex ) * _rows * _rank );
        s.get( _cmat_B.data(), sizeof( complex ) * _cols * _rank );
    }// if
    else
    {
        s.get( _rmat_A.data(), sizeof( real ) * _rows * _rank );
        s.get( _rmat_B.data(), sizeof( real ) * _cols * _rank );
    }// else
}

void
TRkMatrix::write ( TByteStream & s ) const
{
    TMatrix::write( s );

    s.put( _rows );
    s.put( _cols );
    s.put( _rank );

    if ( is_complex() )
    {
        s.put( _cmat_A.data(), sizeof( complex ) * _rows * _rank );
        s.put( _cmat_B.data(), sizeof( complex ) * _cols * _rank );
    }// if
    else
    {
        s.put( _rmat_A.data(), sizeof( real ) * _rows * _rank );
        s.put( _rmat_B.data(), sizeof( real ) * _cols * _rank );
    }// else
}

//
// returns size of object in bytestream
//
size_t
TRkMatrix::bs_size () const
{
    return (TMatrix::bs_size() + sizeof(_rows) + sizeof(_cols) + sizeof(_rank) +
            (_rank * (_rows + _cols) * (is_complex() ? sizeof(complex) : sizeof(real))));
}

//
// test data for invalid values, e.g. INF and NAN
//
void
TRkMatrix::check_data () const
{
    if ( is_complex() )
    {
        _cmat_A.check_data();
        _cmat_B.check_data();
    }// if
    else
    {
        _rmat_A.check_data();
        _rmat_B.check_data();
    }// else
}

}// namespace
