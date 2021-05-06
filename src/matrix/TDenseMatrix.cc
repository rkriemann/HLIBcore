//
// Project     : HLib
// File        : TDenseMatrix.cc
// Description : class for dense matrices of arbitrary (small) size
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/base/error.hh"
#include "hpro/vector/TScalarVector.hh"
#include "hpro/matrix/TRkMatrix.hh"

#include "hpro/matrix/TDenseMatrix.hh"

namespace HLIB
{

// namespace abbr.
namespace B = BLAS;

using std::unique_ptr;
using std::make_unique;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// local defines, functions
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

namespace
{

// defines maximal ratio of rows vs. columns before considered
// as low rank matrix
const size_t  MAX_DIM_RATIO = 30;

//
// dense MVM with optional parallelisation
//
template < typename value_t >
void
task_mulvec ( const value_t                 alpha,
              const matop_t                 op,
              const B::Matrix< value_t > &  M,
              const B::Vector< value_t > &  x,
              const value_t                 beta,
              B::Vector< value_t > &        y )
{
    B::mulvec( alpha, B::mat_view( op, M ), x, beta, y );
}

}// namespace anonymous

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// TDenseMatrix
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//
// set size of matrix
//
void
TDenseMatrix::set_size ( const size_t n, const size_t m )
{
    TScopedLock  mlock( *this );
    
    if (( n != _rows ) || (m != _cols))
    {
        if ( is_complex() )
        {
            if ( n * m > 0 )
            {
                _cmat = B::Matrix< complex >( n, m );
                _rows = n;
                _cols = m;
            }// if
            else
            {            
                _cmat = B::Matrix< complex >();
                _rows = 0;
                _cols = 0;
            }// else
        }// if
        else
        {
            if ( n * m > 0 )
            {
                _rmat = B::Matrix< real >( n, m );
                _rows = n;
                _cols = m;
            }// if
            else
            {            
                _rmat = B::Matrix< real >();
                _rows = 0;
                _cols = 0;
            }// else
        }// else
    }// if
}

//
// switch between complex and real format
//
void
TDenseMatrix::to_real ()
{
    if ( ! is_complex() )
        return;

    if ( _rows * _cols == 0 )
        return;
    
    TScopedLock  mlock( *this );

    for ( uint j = 0; j < _cols; ++j )
        for ( uint i = 0; i < _rows; ++i )
            if ( std::imag( _cmat( i, j ) ) != real(0) )
            {
                TMatrix::set_complex( true );
                HERROR( ERR_NREAL, "(TDenseMatrix) to_real",
                        "matrix has non-zero imaginary part" );
            }// if
    
    _rmat = B::Matrix< real >( _rows, _cols );
    
    for ( uint j = 0; j < _cols; ++j )
        for ( uint i = 0; i < _rows; ++i )
            _rmat( i, j ) = std::real( _cmat( i, j ) );

    _cmat = B::Matrix< complex >();
}

void
TDenseMatrix::to_complex ()
{
    if ( is_complex() )
        return;

    if ( _rows*_cols == 0 )
        return;
    
    TScopedLock  mlock( *this );

    _cmat = B::Matrix< complex >( _rows, _cols );
    
    for ( uint j = 0; j < _cols; ++j )
        for ( uint i = 0; i < _rows; ++i )
            _cmat( i, j ) = _rmat( i, j );

    _rmat = B::Matrix< real >();
}
    
///////////////////////////////////////////////
//
// block-handling
//

//
// do α·this + β·op(M)
// (either M is sub block of this or this a sub block of M)
//
void
TDenseMatrix::add_block ( const real           alpha,
                          const real           beta,
                          const TDenseMatrix * M,
                          const matop_t        op )
{
    if ( M->is_complex() )
        set_complex( true );
        
    TScopedLock  mlock( *this );

    if ( op == apply_normal )
    {
        if ( block_is().is_subset( M->block_is() ) )
        {
            //
            // M is sub matrix of this
            //
            
            if ( is_complex() )
            {
                B::Matrix< complex >  Rthis( _cmat,
                                             intersect( row_is(), M->row_is() ) - row_ofs(),
                                             intersect( col_is(), M->col_is() ) - col_ofs() );

                // this ≔ α·this
                if ( alpha != real(1) )
                {
                    if ( alpha == real(0) ) B::fill( complex(0), Rthis );
                    else                B::scale( complex(alpha), Rthis );
                }// if

                // this ≔ this + β·M
                if ( beta != real(0) )
                {
                    B::add( complex(beta), M->blas_cmat(), Rthis );
                }// else
            }// if
            else
            {
                B::Matrix< real >  Rthis( _rmat,
                                          intersect( row_is(), M->row_is() ) - row_ofs(),
                                          intersect( col_is(), M->col_is() ) - col_ofs() );
                
                // this ≔ α·this
                if ( alpha != real(1) )
                {
                    if ( alpha == real(0) ) B::fill( real(0), Rthis );
                    else                B::scale( alpha, Rthis );
                }// if

                // this ≔ this + β·M
                if ( beta != real(0) )
                {
                    B::add( beta, M->blas_rmat(), Rthis );
                }// else
            }// if
        }// if
        else if ( M->block_is().is_subset( block_is() ) )
        {
            //
            // this is sub matrix of M
            //
            
            if ( is_complex() )
            {
                // this ≔ α·this
                if ( alpha != real(1) )
                {
                    if ( alpha == real(0) ) B::fill( complex(0), _cmat );
                    else                    B::scale( complex(alpha), _cmat );
                }// if

                // this ≔ this + β·M
                if ( beta != real(0) )
                {
                    if ( ! M->is_complex() )
                        HERROR( ERR_REAL_CMPLX, "(TDenseMatrix) add_block", "" );
                        
                    B::Matrix< complex >  RM( M->blas_cmat(),
                                              intersect( row_is(), M->row_is() ) - M->row_ofs(),
                                              intersect( col_is(), M->col_is() ) - M->col_ofs() );

                    B::add( complex(beta), RM, _cmat );
                }// if
            }// if
            else
            {
                // this ≔ α·this
                if ( alpha != real(1) )
                {
                    if ( alpha == real(0) ) B::fill( real(0), _rmat );
                    else                B::scale( alpha, _rmat );
                }// else

                // this ≔ this + β·M
                if ( beta != real(0) )
                {
                    B::Matrix< real >  RM( M->blas_rmat(),
                                           intersect( row_is(), M->row_is() ) - M->row_ofs(),
                                           intersect( col_is(), M->col_is() ) - M->col_ofs() );

                    B::add( beta, RM, _rmat );
                }// if
            }// if
        }// if
        else
            HERROR( ERR_INDEXSET, "(TDenseMatrix) add_block",
                   "given matrix is neither a subblock nor a superblock" );
    }// if
    else if ( op == apply_transposed )
    {
        TBlockIndexSet  MT_bs( HLIB::transpose( M->block_is() ) );

        if ( block_is().is_subset( MT_bs ) )
        {
            //
            // M^T is sub matrix of this
            //
            
            if ( is_complex() )
            {
                B::Matrix< complex >  Rthis( _cmat,
                                             intersect( row_is(), MT_bs.row_is() ) - row_ofs(),
                                             intersect( col_is(), MT_bs.col_is() ) - col_ofs() );

                // this ≔ α·this
                if ( alpha != real(1) )
                {
                    if ( alpha == real(0) ) B::fill( complex(0), Rthis );
                    else                B::scale( complex(alpha), Rthis );
                }// if


                // this ≔ this + β·M
                if ( beta != real(0) )
                {
                    B::add( complex(beta), B::transposed( M->blas_cmat() ), Rthis );
                }// else
            }// if
            else
            {
                B::Matrix< real >  Rthis( _rmat,
                                          intersect( row_is(), MT_bs.row_is() ) - row_ofs(),
                                          intersect( col_is(), MT_bs.col_is() ) - col_ofs() );
                
                // this ≔ α·this
                if ( alpha != real(1) )
                {
                    if ( alpha == real(0) ) B::fill( real(0), Rthis );
                    else                B::scale( alpha, Rthis );
                }// if

                // this ≔ this + β·M
                if ( beta != real(0) )
                {
                    B::add( beta, B::transposed( M->blas_rmat() ), Rthis );
                }// else
            }// else
        }// if
        else
            HERROR( ERR_INDEXSET, "(TDenseMatrix) add_block",
                   "given matrix is neither a subblock nor a superblock" );
    }// if
    else // if ( op == apply_adjoint )
    {
        TBlockIndexSet  MT_bs( HLIB::transpose( M->block_is() ) );

        if ( block_is().is_subset( MT_bs ) )
        {
            //
            // M^T is sub matrix of this
            //
        
            if ( is_complex() )
            {
                B::Matrix< complex >  Rthis( _cmat,
                                             intersect( row_is(), MT_bs.row_is() ) - row_ofs(),
                                             intersect( col_is(), MT_bs.col_is() ) - col_ofs() );

                // this ≔ α·this
                if ( alpha != real(1) )
                {
                    if ( alpha == real(0) ) B::fill( complex(0), Rthis );
                    else                B::scale( complex(alpha), Rthis );
                }// if

                // this ≔ this + β·M
                if ( beta != real(0) )
                {
                    B::add( complex(beta), B::adjoint( M->blas_cmat() ), Rthis );
                }// else
            }// if
            else
            {
                B::Matrix< real >  Rthis( _rmat,
                                          intersect( row_is(), MT_bs.row_is() ) - row_ofs(),
                                          intersect( col_is(), MT_bs.col_is() ) - col_ofs() );
                
                // this ≔ α·this
                if ( alpha != real(1) )
                {
                    if ( alpha == real(0) ) B::fill( real(0), Rthis );
                    else                B::scale( alpha, Rthis );
                }// if

                // this ≔ this + β·M
                if ( beta != real(0) )
                {
                    B::add( beta, B::adjoint( M->blas_rmat() ), Rthis );
                }// else
            }// else
        }// if
        else
            HERROR( ERR_INDEXSET, "(TDenseMatrix) add_block",
                   "given matrix is neither a subblock nor a superblock" );
    }// else
}

///////////////////////////////////////////////
//
// handle size (via cluster)
//

void
TDenseMatrix::set_cluster ( const TBlockCluster * c )
{
    TMatrix::set_cluster( c );

    // set size to the size of the cluster
    if ( c != nullptr )
        set_size( c->rowcl()->size(), c->colcl()->size() );
}

/////////////////////////////////////////////////
//
// BLAS-routines (real valued)
//

//
// scale matrix by constant factor
//
void
TDenseMatrix::scale ( const real f )
{
    if ( f == real(1) )
        return;

    TScopedLock  mlock( *this );

    if ( is_complex() )
    {
        if ( f == real(0) ) B::fill<complex>(  real(0), _cmat );
        else            B::scale<complex>(   f, _cmat );
    }// if
    else
    {
        if ( f == real(0) ) B::fill<real>(  real(0), _rmat );
        else            B::scale<real>(   f, _rmat );
    }// else
}

//
// matrix-vector-multiplication
//
void
TDenseMatrix::mul_vec ( const real      alpha,
                        const TVector * x,
                        const real      beta,
                        TVector       * y,
                        const matop_t   op ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TDenseMatrix) mul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TDenseMatrix) mul_vec", "y = nullptr" );
    
    if ( op == apply_normal )
    {
        if (( row_is() != y->is() ) || ( col_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TDenseMatrix) mul_vec", "incompatible vector index set" );
    }// if
    else
    {
        if (( col_is() != y->is() ) || ( row_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TDenseMatrix) mul_vec", "incompatible vector index set" );
    }// if
    
    if ( ! IS_TYPE( x, TScalarVector ) || ! IS_TYPE( y, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TDenseMatrix) mul_vec",
                y->typestr() + " += TDenseMatrix * " + x->typestr() );
        
    //
    // check if scalar vector
    //

    auto  s_y = ptrcast( y, TScalarVector );
    auto  s_x = cptrcast( x, TScalarVector );

    if ( ! (( is_complex() == x->is_complex() ) && ( is_complex() == y->is_complex() )) )
        HERROR( ERR_REAL_CMPLX, "(TDenseMatrix) mul_vec", "" );
        
    if ( is_complex() )
    {
        task_mulvec< complex >( alpha, op, blas_cmat(), s_x->blas_cvec(), beta, s_y->blas_cvec() );
    }// if
    else
    {
        task_mulvec< real >( alpha, op, blas_rmat(), s_x->blas_rvec(), beta, s_y->blas_rvec() );
    }// else
}

//
// compute this ≔ this + a·M
//
void
TDenseMatrix::add ( const real a, const TMatrix * M )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "(TDenseMatrix) add", "argument is nullptr" );
    
    if ( ! IS_TYPE( M, TDenseMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TDenseMatrix) add", M->typestr() );

    const TDenseMatrix *  D = cptrcast( M, TDenseMatrix );

    //
    // some checks
    //
    
    if ( a == real(0) )
        return;

    if ( block_is() != D->block_is() )
        HERROR( ERR_MAT_SIZE, "(TDenseMatrix) add", "" );

    //
    // add matrix
    //

    if ( D->is_complex() )
        set_complex( true );

    TScopedLock  mlock( *this );
    
    if ( D->is_complex() )
        B::add<complex>( a, D->blas_cmat(), blas_cmat() );
    else
    {
        if ( is_complex() )
            HERROR( ERR_COMPLEX, "(TDenseMatrix) add", "can not mix real and complex" );
        else
            B::add( a, D->blas_rmat(), blas_rmat() );
    }// else
}

//
// matrix-matrix-multiplication: alpha op_A(this) * op_B(B)
//
TMatrix *
TDenseMatrix::mul_right ( const real       /* alpha */,
                          const TMatrix *  /* B */,
                          const matop_t    /* op_A */,
                          const matop_t    /* op_B */ ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// matrix-matrix-multiplication: alpha op_A(A)·op_B(this)
//
TMatrix *
TDenseMatrix::mul_left ( const real       /* alpha */,
                         const TMatrix *  /* A */,
                         const matop_t    /* op_A */,
                         const matop_t    /* op_B */ ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

/////////////////////////////////////////////////
//
// BLAS-routines (complex valued)
//

//
// scale matrix by constant factor
//
void
TDenseMatrix::cscale ( const complex f )
{
    if ( f == real(1) )
        return;

    TScopedLock  mlock( *this );

    if ( is_complex() )
    {
        if ( f == real(0) ) B::fill<complex>(  real(0), _cmat );
        else            B::scale<complex>(   f, _cmat );
    }// if
    else
    {
        HERROR( ERR_REAL_CMPLX, "(TDenseMatrix) cscale", "real matrix, complex factor" );
        
        // if ( f == real(0) )
        //     B::fill<real>(  real(0), _rmat );
        // else if ( std::imag( f ) != real(0) )
        // {
        //     set_complex( true );
        //     B::scale<complex>( f, _cmat );
        // }
        // else
        //     B::scale<real>( std::real(f), _rmat );
    }// else
}

//
// matrix-vector-multiplication
//
void
TDenseMatrix::cmul_vec ( const complex   alpha,
                         const TVector * x,
                         const complex   beta,
                         TVector       * y,
                         const matop_t   op ) const
{
    if ( x == nullptr ) HERROR( ERR_ARG, "(TDenseMatrix) cmul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TDenseMatrix) cmul_vec", "y = nullptr" );
    
    if ( op == apply_normal )
    {
        if (( row_is() != y->is() ) || ( col_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TDenseMatrix) cmul_vec", "incompatible vector index set" );
    }// if
    else
    {
        if (( col_is() != y->is() ) || ( row_is() != x->is() ))
            HERROR( ERR_INDEXSET, "(TDenseMatrix) cmul_vec", "incompatible vector index set" );
    }// if
    
    if ( ! IS_TYPE( x, TScalarVector ) || ! IS_TYPE( y, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TDenseMatrix) cmul_vec",
                y->typestr() + " += TDenseMatrix * " + x->typestr() );
        
    //
    // multiply
    //

    auto  s_y = ptrcast( y, TScalarVector );
    auto  s_x = cptrcast( x, TScalarVector );
                
    if ( ! (( is_complex() == x->is_complex() ) && ( is_complex() == y->is_complex() )) )
        HERROR( ERR_COMPLEX, "(TDenseMatrix) mul_vec", "can not mix real and complex" );
        
    if ( is_complex() )
    {
        task_mulvec< complex >( alpha, op, blas_cmat(), s_x->blas_cvec(), beta, s_y->blas_cvec() );
    }// if
    else
    {
        if ( ( std::imag( alpha ) != real(0) ) || ( std::imag( beta ) != real(0) ) )
            HERROR( ERR_COMPLEX, "(TDenseMatrix) mul_vec", "can not mix real and complex" );
        
        task_mulvec< real >( std::real( alpha ), op, blas_rmat(), s_x->blas_rvec(), std::real( beta ), s_y->blas_rvec() );
    }// else
}

//
// compute this = this + a * matrix
//
void
TDenseMatrix::cadd ( const complex a, const TMatrix * M )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "(TDenseMatrix) cadd", "argument is nullptr" );
    
    if ( ! IS_TYPE( M, TDenseMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TDenseMatrix) cadd", M->typestr() );

    const TDenseMatrix *  D = cptrcast( M, TDenseMatrix );

    //
    // some checks
    //
    
    if ( a == real(0) )
        return;

    if ( block_is() != D->block_is() )
        HERROR( ERR_MAT_SIZE, "(TDenseMatrix) cadd", "" );

    //
    // add matrix
    //

    if ( M->is_complex() || ( std::imag( a ) != real(0) ))
        set_complex( true );

    TScopedLock  mlock( *this );
    
    if ( is_complex() )
    {
        if ( D->is_complex() )
            B::add( a, D->blas_cmat(), blas_cmat() );
        else
            HERROR( ERR_COMPLEX, "(TDenseMatrix) cadd", "can not mix real and complex" );
    }// if
    else
        B::add( std::real(a), D->blas_rmat(), blas_rmat() );
}

//
// matrix-matrix-multiplication: alpha op_A(this) * op_B(B)
//
TMatrix *
TDenseMatrix::cmul_right ( const complex    /* alpha */,
                           const TMatrix *  /* B */,
                           const matop_t    /* op_A */,
                           const matop_t    /* op_B */ ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// matrix-matrix-multiplication: alpha op_A(A) * op_B(this)
//
TMatrix *
TDenseMatrix::cmul_left ( const complex    /* alpha */,
                          const TMatrix *  /* A */,
                          const matop_t    /* op_A */,
                          const matop_t    /* op_B */ ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}
    
//
// do α * this + β * op(M)
// (either M is subblock of this or this a subblock of M)
//
void
TDenseMatrix::add_block ( const complex        alpha,
                          const complex        beta,
                          const TDenseMatrix * M,
                          const matop_t        op )
{
    if ( M->is_complex() || ( std::imag( alpha ) != real(0) ) || ( std::imag( beta ) != real(0) ) )
        set_complex( true );
        
    TScopedLock  mlock( *this );
    
    if ( op == apply_normal )
    {
        //
        // check if bigger or smaller
        //

        if ( block_is().is_subset( M->block_is() ) )
        {
            //
            // M is sub matrix of this
            //
            
            if ( is_complex() )
            {
                B::Matrix< complex >  Rthis( _cmat,
                                             intersect( row_is(), M->row_is() ) - row_ofs(),
                                             intersect( col_is(), M->col_is() ) - col_ofs() );

                // this ≔ α·this
                if ( alpha != complex(1) )
                {
                    if ( alpha == real(0) ) B::fill( complex(0), Rthis );
                    else                B::scale( alpha, Rthis );
                }// if

                // this ≔ this + β·M
                if ( beta != complex(0) )
                {
                    B::add( beta, M->blas_cmat(), Rthis );
                }// else
            }// if
            else
            {
                B::Matrix< real >  Rthis( _rmat,
                                          intersect( row_is(), M->row_is() ) - row_ofs(),
                                          intersect( col_is(), M->col_is() ) - col_ofs() );
                
                // this ≔ α·this
                if ( alpha != complex(1) )
                {
                    if ( alpha == complex(0) ) B::fill( real(0), _rmat );
                    else                       B::scale( std::real(alpha), _rmat );
                }// if

                // this ≔ this + β·M
                if ( beta != complex(0) )
                {
                    B::add( std::real(beta), M->blas_rmat(), Rthis );
                }// else
            }// if
        }// if
        else if ( M->block_is().is_subset( block_is() ) )
        {
            //
            // this is sub matrix of M
            //
            
            if ( is_complex() )
            {
                // this ≔ α·this
                if ( alpha != complex(1) )
                {
                    if ( alpha == real(0) ) B::fill( complex(0), _cmat );
                    else                B::scale( alpha, _cmat );
                }// if

                // this ≔ this + β·M
                if ( beta != complex(0) )
                {
                    if ( ! M->is_complex() )
                        HERROR( ERR_REAL_CMPLX, "(TDenseMatrix) add_block", "" );
                        
                    B::Matrix< complex >  RM( M->blas_cmat(),
                                              intersect( row_is(), M->row_is() ) - M->row_ofs(),
                                              intersect( col_is(), M->col_is() ) - M->col_ofs() );

                    B::add( beta, RM, _cmat );
                }// if
            }// if
            else
            {
                // this ≔ α·this
                if ( alpha != complex(1) )
                {
                    if ( alpha == real(0) ) B::fill( real(0), _rmat );
                    else                    B::scale( std::real(alpha), _rmat );
                }// if

                // this ≔ this + β·M
                if ( beta != complex(0) )
                {
                    B::Matrix< real >  RM( M->blas_rmat(),
                                           intersect( row_is(), M->row_is() ) - M->row_ofs(),
                                           intersect( col_is(), M->col_is() ) - M->col_ofs() );

                    B::add( std::real(beta), RM, _rmat );
                }// if
            }// if
        }// if
        else
            HERROR( ERR_INDEXSET, "(TDenseMatrix) add_block",
                   "given matrix is neither a subblock nor a superblock" );
    }// if
    else if ( op == apply_transposed )
    {
        TBlockIndexSet  MT_bs( HLIB::transpose( M->block_is() ) );

        if ( block_is().is_subset( MT_bs ) )
        {
            //
            // M^T is sub matrix of this
            //
            
            if ( is_complex() )
            {
                B::Matrix< complex >  Rthis( _cmat,
                                             intersect( row_is(), MT_bs.row_is() ) - row_ofs(),
                                             intersect( col_is(), MT_bs.col_is() ) - col_ofs() );

                // this ≔ α·this
                if ( alpha != complex(1) )
                {
                    if ( alpha == real(0) ) B::fill( complex(0), Rthis );
                    else                    B::scale( alpha, Rthis );
                }// if

                // this ≔ this + β·M
                if ( beta != complex(0) )
                {
                    B::add( beta, B::transposed( M->blas_cmat() ), Rthis );
                }// else
            }// if
            else
            {
                B::Matrix< real >  Rthis( _rmat,
                                          intersect( row_is(), MT_bs.row_is() ) - row_ofs(),
                                          intersect( col_is(), MT_bs.col_is() ) - col_ofs() );
                
                // this ≔ α·this
                if ( alpha != complex(1) )
                {
                    if ( alpha == complex(0) ) B::fill( real(0), Rthis );
                    else                       B::scale( std::real(alpha), Rthis );
                }// if

                // this ≔ this + β·M
                if ( beta != complex(0) )
                {
                    B::add( std::real(beta), B::transposed( M->blas_rmat() ), Rthis );
                }// else
            }// else
        }// if
        else
            HERROR( ERR_INDEXSET, "(TDenseMatrix) add_block",
                   "given matrix is neither a subblock nor a superblock" );
    }// if
    else // if ( op == apply_adjoint )
    {
        TBlockIndexSet  MT_bs( HLIB::transpose( M->block_is() ) );

        if ( block_is().is_subset( MT_bs ) )
        {
            //
            // M^T is sub matrix of this
            //
            
            if ( is_complex() )
            {
                B::Matrix< complex >  Rthis( _cmat,
                                             intersect( row_is(), MT_bs.row_is() ) - row_ofs(),
                                             intersect( col_is(), MT_bs.col_is() ) - col_ofs() );

                // this ≔ α·this
                if ( alpha != complex(1) )
                {
                    if ( alpha == complex(0) ) B::fill( complex(0), Rthis );
                    else                       B::scale( alpha, Rthis );
                }// if

                // this ≔ this + β·M
                if ( beta != complex(0) )
                {
                    B::add( beta, B::adjoint( M->blas_cmat() ), Rthis );
                }// else
            }// if
            else
            {
                B::Matrix< real >  Rthis( _rmat,
                                          intersect( row_is(), MT_bs.row_is() ) - row_ofs(),
                                          intersect( col_is(), MT_bs.col_is() ) - col_ofs() );
                
                // this ≔ α·this
                if ( alpha != complex(1) )
                {
                    if ( alpha == complex(0) ) B::fill( real(0), Rthis );
                    else                       B::scale( std::real(alpha), Rthis );
                }// if

                // this ≔ this + β·M
                if ( beta != complex(0) )
                {
                    B::add( std::real(beta), B::adjoint( M->blas_rmat() ), Rthis );
                }// else
            }// else
        }// if
        else
            HERROR( ERR_INDEXSET, "(TDenseMatrix) add_block",
                   "given matrix is neither a subblock nor a superblock" );
    }// else
}

///////////////////////////////////////////////////////////
//
// linear operator mapping
//

//! same as above but only the dimension of the vector spaces is tested,
//! not the corresponding index sets
void
TDenseMatrix::apply_add   ( const real                       alpha,
                            const BLAS::Vector< real > &     x,
                            BLAS::Vector< real > &           y,
                            const matop_t                    op ) const
{
    HASSERT( ( x.length() == ncols( op ) ) && ( y.length() == nrows( op ) ),
             ERR_ARG, "(TDenseMatrix) apply_add", "incompatible vector dimensions" );
    HASSERT( is_real(), 
             ERR_REAL_CMPLX, "(TDenseMatrix) apply_add", "real vectors, complex matrix" );

    BLAS::mulvec( alpha, mat_view( op, blas_mat< real >( this ) ), x, real(1), y );
}
    
void
TDenseMatrix::apply_add   ( const complex                    alpha,
                            const BLAS::Vector< complex > &  x,
                            BLAS::Vector< complex > &        y,
                            const matop_t                    op ) const
{
    HASSERT( ( x.length() == ncols( op ) ) && ( y.length() == nrows( op ) ),
             ERR_ARG, "(TDenseMatrix) apply_add", "incompatible vector dimensions" );
    HASSERT( is_complex(), 
             ERR_REAL_CMPLX, "(TDenseMatrix) apply_add", "complex vectors, real matrix" );

    BLAS::mulvec( alpha, mat_view( op, blas_mat< complex >( this ) ), x, complex(1), y );
}

/////////////////////////////////////////////////
//
// misc methods
//

//
// transpose matrix
//
void
TDenseMatrix::transpose ()
{
    {
        TScopedLock  mlock( *this );
    
        if ( is_complex() )
        {
            B::Matrix< complex >  T( cols(), rows() );
            
            B::copy( transposed( blas_cmat() ), T );
            blas_cmat() = std::move( T );
        }// if
        else
        {
            B::Matrix< real >  T( cols(), rows() );
            
            B::copy( transposed( blas_rmat() ), T );
            blas_rmat() = std::move( T );
        }// else
        
        std::swap( _rows, _cols );
    }

    TMatrix::transpose();
}
    
//
// conjugate matrix coefficients
//
void
TDenseMatrix::conjugate ()
{
    if ( is_complex() && ! is_hermitian() )
    {
        TScopedLock  mlock( *this );
        
        B::conj( blas_cmat() );
    }// if
}

//
// virtual constructor
//
std::unique_ptr< TMatrix >
TDenseMatrix::copy () const
{
    auto  M = TMatrix::copy();
    auto  D = ptrcast( M.get(), TDenseMatrix );

    if ( is_complex() )
        B::copy( blas_cmat(), D->blas_cmat() );
    else
        B::copy( blas_rmat(), D->blas_rmat() );

    return M;
}

//
// copy matrix into A
//
void
TDenseMatrix::copy_to ( TMatrix * A ) const
{
    TMatrix::copy_to( A );
    
    if ( ! IS_TYPE( A, TDenseMatrix ) )
    {
        HERROR( ERR_NOT_IMPL, "", "" );
        
        // convert( this, A, acc_exact );
        // return;
    }// if
    
    TDenseMatrix * D = ptrcast( A, TDenseMatrix );

//     if ((D->rows() != rows()) || (D->cols() != cols()))
//         WARNING( ERR_MAT_SIZE, "(TDenseMatrix) copy_to", "" );

    D->set_complex( is_complex() );
    D->set_size( rows(), cols() );

    if ( is_complex() )
        B::copy( blas_cmat(), D->blas_cmat() );
    else
        B::copy( blas_rmat(), D->blas_rmat() );
}

//
// return structural copy
//
std::unique_ptr< TMatrix >
TDenseMatrix::copy_struct () const
{
    return TMatrix::copy_struct();
}

//
// return size in bytes used by this object
//
size_t
TDenseMatrix::byte_size () const
{
    size_t  size = TMatrix::byte_size() + 2*sizeof(size_t) + sizeof(B::Matrix<real>) + sizeof(B::Matrix<complex>);

    if ( is_complex() ) size += _rows * _cols * sizeof(complex);
    else                size += _rows * _cols * sizeof(real);

    return size;
}

//
// serialisation
//

void
TDenseMatrix::read  ( TByteStream & s )
{
    TMatrix::read( s );

    size_t  n, m;
    
    s.get( n );
    s.get( m );

    set_size( n, m );

    if ( is_complex() )
        s.get( blas_cmat().data(), sizeof( complex ) * _rows * _cols );
    else
        s.get( blas_rmat().data(), sizeof( real ) * _rows * _cols );
}

void
TDenseMatrix::build ( TByteStream & s )
{
    TMatrix::build( s );

    size_t  n, m;
    
    s.get( n );
    s.get( m );

    set_size( n, m );

    if ( is_complex() )
        s.get( blas_cmat().data(), sizeof( complex ) * _rows * _cols );
    else
        s.get( blas_rmat().data(), sizeof( real ) * _rows * _cols );
}

void
TDenseMatrix::write ( TByteStream & s ) const
{
    TMatrix::write( s );

    s.put( _rows );
    s.put( _cols );

    if ( is_complex() )
        s.put( blas_cmat().data(), sizeof( complex ) * _rows * _cols );
    else
        s.put( blas_rmat().data(), sizeof( real ) * _rows * _cols );
}

//
// returns size of object in bytestream
//
size_t
TDenseMatrix::bs_size () const
{
    return (TMatrix::bs_size() + sizeof(_rows) + sizeof(_cols) +
            (_rows * _cols * ( is_complex() ? sizeof(complex) : sizeof(real) ) ) );
}

///////////////////////////////////////////////////////////////////////////
//
// matrix permutation functions
//
///////////////////////////////////////////////////////////////////////////

namespace
{

//
// copy row/column <j> to row/column <i> in matrix M
//
template <typename T>
void
copy ( B::Matrix< T > & M,
       const bool       sort_rows,
       const idx_t      i,
       const idx_t      j )
{
    if ( sort_rows )
    {
        B::Vector< T >  row_i( M.row( i ) );
        B::Vector< T >  row_j( M.row( j ) );
            
        B::copy( row_j, row_i );
    }// if
    else
    {
        B::Vector< T >  col_i( M.column( i ) );
        B::Vector< T >  col_j( M.column( j ) );
            
        B::copy( col_j, col_i );
    }// else
}

//
// copy row/column <i> to vector <tmp>
//
template <typename T>
void
copy ( B::Matrix< T > &  M,
       const bool        sort_rows,
       const idx_t       i,
       B::Vector< T > &  tmp )
{
    if ( sort_rows )
    {
        B::Vector< T >  row_i( M.row( i ) );
            
        B::copy( row_i, tmp );
    }// if
    else
    {
        B::Vector< T >  col_i( M.column( i ) );
            
        B::copy( col_i, tmp );
    }// else
}

//
// copy vector <tmp> to row/column <i>
//
template <typename T>
void
copy ( B::Vector< T > &  tmp,
       const bool        sort_rows,
       const idx_t       i,
       B::Matrix< T > &  M )
{
    if ( sort_rows )
    {
        B::Vector< T >  row_i( M.row( i ) );
            
        B::copy( tmp, row_i );
    }// if
    else
    {
        B::Vector< T >  col_i( M.column( i ) );
            
        B::copy( tmp, col_i );
    }// else
}

//
// swap the rows/columns <i> and <j> in matrix M
//
template <typename T>
void
swap ( B::Matrix< T > &  M,
       const bool        sort_rows,
       const idx_t       i,
       const idx_t       j,
       B::Vector< T > &  tmp )
{
    if ( sort_rows )
    {
        B::Vector< T >  row_i( M.row( i ) );
        B::Vector< T >  row_j( M.row( j ) );

        B::copy( row_i, tmp );
        B::copy( row_j, row_i );
        B::copy( tmp,   row_j );
    }// if
    else
    {
        B::Vector< T >  col_i( M.column( i ) );
        B::Vector< T >  col_j( M.column( j ) );

        B::copy( col_i, tmp );
        B::copy( col_j, col_i );
        B::copy( tmp,   col_j );
    }// else
}

//
// sort rows/columns of matrix according to given permutation
//
template <typename T>
void
matrix_sort( B::Matrix< T > &  M,
             const bool        sort_rows,
             TPermutation &    perm,
             const idx_t       lb,
             const idx_t       ub,
             B::Vector< T > &  tmp )
{
    if ( lb >= ub ) return;

    if ( (ub - lb) < 20 )
    {
        //
        // apply insertion sort for small ranges
        //

        for ( idx_t  i = lb+1; i <= ub; i++ )
        {
            const idx_t  v = perm[i];
            idx_t        j = i-1;

            copy( M, sort_rows, i, tmp );
            
            while (( j >= 0 ) && ( perm[j] > v ))
            {
                copy( M, sort_rows, j+1, j );
                perm[j+1] = perm[j];
                j--;
            }// if

            copy( tmp, sort_rows, j+1, M );
            perm[j+1] = v;
        }// for
    }// if
    else
    {
        //
        // apply quick sort for standard ranges
        //

        idx_t        i         = lb;
        idx_t        j         = ub;
        const idx_t  mid       = (lb + ub) / 2;
        idx_t        choice[3] = { perm[lb], perm[mid], perm[ub] };
        idx_t        pivot;

        // choose pivot (median-of-three)
        if ( choice[0] > choice[1] ) std::swap( choice[0], choice[1] );
        if ( choice[0] > choice[2] ) std::swap( choice[0], choice[2] );
        if ( choice[1] > choice[2] ) std::swap( choice[1], choice[2] );
        pivot = choice[1];

        // partition
        while ( i < j )
        {
            while ( perm[i] < pivot   ) i++;
            while ( pivot   < perm[j] ) j--;

            if ( i < j )
            {
                swap( M, sort_rows, i, j, tmp );
                std::swap( perm[i], perm[j] );
            }// if
        }// while

        // recursion
        matrix_sort( M, sort_rows, perm, lb, i-1, tmp );
        matrix_sort( M, sort_rows, perm, i+1, ub, tmp );
    }// else
}

}// namespace anonymous

//
// permute rows/columns in matrix
//
void
TDenseMatrix::permute ( const TPermutation *  row_perm,
                        const TPermutation *  col_perm )
{
#if 1

    //
    // apply permutations by sorting the rows/columns
    // simultaneously to sorting permutations, thereby
    // reverting order
    //
    
    TPermutation   perm;

    if ( is_complex() )
    {
        B::Vector< complex >  tmp;

        if ( row_perm != nullptr )
        {
            tmp  = B::Vector< complex >( cols() );
            perm = * row_perm;
            matrix_sort( blas_cmat(), true, perm, 0, idx_t(rows())-1, tmp );
        }// if
    
        if ( col_perm != nullptr )
        {
            tmp  = B::Vector< complex >( rows() );
            perm = * col_perm;
            matrix_sort( blas_cmat(), false, perm, 0, idx_t(cols())-1, tmp );
        }// if
    }// if
    else
    {
        B::Vector< real >  tmp( std::max( rows(), cols() ) );

        if ( row_perm != nullptr )
        {
            tmp  = B::Vector< real >( cols() );
            perm = * row_perm;
            matrix_sort( blas_rmat(), true, perm, 0, idx_t(rows())-1, tmp );
        }// if
    
        if ( col_perm != nullptr )
        {
            tmp  = B::Vector< real >( rows() );
            perm = * col_perm;
            matrix_sort( blas_rmat(), false, perm, 0, idx_t(cols())-1, tmp );
        }// if
    }// else
        
#else
    
    unique_ptr< TDenseMatrix >  T( ptrcast( copy().release(), TDenseMatrix ) );
    const uint                  n = rows();
    const uint                  m = cols();
    
    //
    // apply row permutation first by copying the old rows
    // to the new position defined by the permutation
    //

    if ( row_perm != nullptr )
    {
        if ( is_complex() )
        {
            for ( uint i = 0; i < n; i++ )
            {
                const uint pi = row_perm->permute( i );

                B::Vector< complex >  this_i( blas_cmat().row( i ) );
                B::Vector< complex >  T_pi( T->blas_cmat().row( pi ) );
            
                B::copy( this_i, T_pi );
            }// for

            B::copy( T->blas_cmat(), blas_cmat() );
        }// if
        else
        {
            for ( uint i = 0; i < n; i++ )
            {
                const uint pi = row_perm->permute( i );

                B::Vector< real >  this_i( blas_rmat().row( i ) );
                B::Vector< real >  T_pi( T->blas_rmat().row( pi ) );
            
                B::copy( this_i, T_pi );
            }// for

            B::copy( T->blas_rmat(), blas_rmat() );
        }// else
    }// if

    //
    // now apply column permutation in the same way
    //

    if ( col_perm != nullptr )
    {
        if ( is_complex() )
        {
            for ( uint j = 0; j < m; j++ )
            {
                const uint pj = col_perm->permute( j );

                B::Vector< complex >  this_j( blas_cmat().column( j ) );
                B::Vector< complex >  T_pj( T->blas_cmat().column( pj ) );
            
                B::copy( this_j, T_pj );
            }// for

            B::copy( T->blas_cmat(), blas_cmat() );
        }// if
        else
        {
            for ( uint j = 0; j < m; j++ )
            {
                const uint pj = col_perm->permute( j );

                B::Vector< real >  this_j( blas_rmat().column( j ) );
                B::Vector< real >  T_pj( T->blas_rmat().column( pj ) );
            
                B::copy( this_j, T_pj );
            }// for

            B::copy( T->blas_rmat(), blas_rmat() );
        }// else
    }// if
    
#endif
}

//
// test data for invalid values, e.g. INF and NAN
//
void
TDenseMatrix::check_data () const
{
    if ( is_complex() )
    {
        _cmat.check_data();
    }// if
    else
    {
        _rmat.check_data();
    }// else
}

}// namespace
