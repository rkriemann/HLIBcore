//
// Project     : HLib
// File        : TBlockMatrix.cc
// Description : class for a matrix consisting of submatrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <list>
#include <iostream>

#include "hpro/vector/TBlockVector.hh"

#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/matrix/TBSHMBuilder.hh"
#include "hpro/matrix/structure.hh"

#include "hpro/matrix/TBlockMatrix.hh"

namespace HLIB
{

// namespace abbr.
namespace B = BLAS;

using namespace std;

///////////////////////////////////////////
//
// constructor and destructor
//

TBlockMatrix::~TBlockMatrix ()
{
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
            if ( block(i,j) != nullptr )
                delete block(i,j);

    delete[] _blocks;
}

//
// set cluster of matrix
//
void
TBlockMatrix::set_cluster ( const TBlockCluster * c )
{
    TMatrix::set_cluster( c );

    if ( c != nullptr )
    {
        set_block_struct( c->nrows(), c->ncols() );
        set_size( c->rowcl()->size(), c->colcl()->size() );
    }// if
}

//
// set block-structure
//
void
TBlockMatrix::set_block_struct ( const uint bn, const uint bm )
{
    // adjust values
    uint n = bn;
    uint m = bm;
    
    if (( n == 0 ) && ( m != 0 )) n = 1;
    if (( n != 0 ) && ( m == 0 )) m = 1;
    
    TMatrix  ** tmp = new TMatrix*[ n * m ];
    const uint  mrc = min(_block_rows*_block_cols, n*m);

    for ( uint i = 0; i < mrc; i++ )
        tmp[i] = _blocks[i];

    for ( uint i = mrc; i < _block_rows*_block_cols; i++ )
        delete _blocks[i];
    
    for ( uint i = mrc; i < n*m; i++ )
        tmp[i] = nullptr;
    
    delete[] _blocks;

    _block_rows = n;
    _block_cols = m;
    _blocks = tmp;
}

//
// set matrix format
//
void
TBlockMatrix::set_form ( const matform_t  f )
{
    TMatrix::set_form( f );

    // also change for diagonal sub blocks
    for ( uint i = 0; i < min( block_rows(), block_cols() ); ++i )
    {
        if ( block( i, i ) != nullptr )
            block( i, i )->set_form( f );
    }// for
}
    
void
TBlockMatrix::set_form ( const matform_t         f,
                         const recursion_type_t  t )
{
    if ( t == recursive )
    {
        set_form( f );
    }// if
    else
    {
        // only change locally
        TMatrix::set_form( f );
    }// else
}
    
//
// replace matrix block <A> by matrix <B> (don't delete A !)
//
void
TBlockMatrix::replace_block ( TMatrix * A, TMatrix * B )
{
    if ( A == nullptr )
        HWARNING( "in (TBlockMatrix) replace_block : replacing nullptr block" );

    TScopedLock  slock( * this );
    
    for ( uint i = 0; i < block_rows(); i++ )
    {
        for ( uint j = 0; j < block_cols(); j++ )
        {
            if ( block( i, j ) == A )
            {
                set_block( i, j, B );
                return;
            }// if
        }// for
    }// for
    
    HERROR( ERR_MAT_STRUCT, "(TBlockMatrix) replace_block", "given matrix not a sub block" );
}

//
// clear pointers to all subblocks
//
void
TBlockMatrix::clear_blocks ()
{
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
            set_block( i, j, nullptr );
}
    
//
// access single matrix coefficient
//
real
TBlockMatrix::entry ( const idx_t row, const idx_t col ) const
{
    // translate index pair to lower half in case of symmetry
    idx_t  srow = row;
    idx_t  scol = col;

    if ( ! is_nonsym() && ( srow < scol ))
        swap( srow, scol );

    const idx_t rofs = row_ofs();
    const idx_t cofs = col_ofs();
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            const TMatrix * A_ij = block(i,j);

            if ( A_ij == nullptr )
                continue;
            
            if ( A_ij->block_is().is_in( rofs + srow, cofs + scol ) )
                return A_ij->entry( srow - A_ij->row_ofs() + rofs, scol - A_ij->col_ofs() + cofs );
        }// for

    return real(0);
}

const complex
TBlockMatrix::centry ( const idx_t row, const idx_t col ) const
{
    // translate index pair to lower half in case of symmetry
    idx_t  srow     = row;
    idx_t  scol     = col;
    bool   swap_idx = false;

    if ( ! is_nonsym() && ( srow < scol ))
    {
        swap( srow, scol );
        swap_idx = true;
    }// if

    const idx_t  rofs = row_ofs();
    const idx_t  cofs = col_ofs();
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            const TMatrix * A_ij = block(i,j);

            if ( A_ij == nullptr )
                continue;
            
            if ( A_ij->block_is().is_in( rofs + srow, cofs + scol ) )
            {
                const complex val = A_ij->centry( srow - A_ij->row_ofs() + rofs,
                                                  scol - A_ij->col_ofs() + cofs );
                
                if ( swap_idx && is_hermitian() )
                    return conj(val);
                else
                    return val;
            }// if
        }// for

    return real(0);
}

//
// set processor set of matrix
//
void
TBlockMatrix::set_procs ( const TProcSet &        ps,
                          const recursion_type_t  t )
{
    TMatrix::set_procs( ps, nonrecursive );

    if ( t == recursive )
    {
        for ( uint i = 0; i < block_rows(); i++ )
            for ( uint j = 0; j < block_cols(); j++ )
            {
                TMatrix * A_ij = block(i,j);

                if ( A_ij != nullptr )
                    A_ij->set_procs( ps, t );
            }// for
    }// if
}

/////////////////////////////////////////////////
//
// update accumulator
//
void
TBlockMatrix::apply_updates ( const TTruncAcc &       acc,
                              const recursion_type_t  recursion )
{
    // DBG::printf( "apply %d", id() );
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// return true, if matrix has updates not yet applied
//
bool
TBlockMatrix::has_updates ( const recursion_type_t  recursion ) const
{
    if ( TMatrix::has_updates( recursion ) )
        return true;

    if ( recursion )
    {
        for ( uint  i = 0; i < block_rows(); ++i )
            for ( uint  j = 0; j < block_cols(); ++j )
            {
                if (( block( i, j ) != nullptr ) && block( i, j )->has_updates( recursion ) )
                    return true;
            }// for
    }// if

    return false;
}


/////////////////////////////////////////////////
//
// BLAS-routines (real valued)
//

//
// scale matrix by constant factor
//
void
TBlockMatrix::scale ( const real f )
{
    if ( f == real(1) )
        return;
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            TMatrix  * A_ij = block(i,j);
            
            if ( A_ij != nullptr )
                A_ij->scale( f );
        }// for
}
    
//
// matrix-vector-mult.
//

void
TBlockMatrix::mul_vec ( const real      alpha,
                        const TVector * x,
                        const real      beta,
                        TVector       * y,
                        const matop_t   op ) const
{
    // if ( ! is_small( this ) || ( procs().size() > 1 ))
    // {
    //     HLIB::mul_vec( procs(), alpha, this, x, beta, y, op );
    //     return;
    // }// if
    
    if ( x == nullptr ) HERROR( ERR_ARG, "(TBlockMatrix) mul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TBlockMatrix) mul_vec", "y = nullptr" );

    if (( row_is( op ) != y->is() ) || ( col_is( op ) != x->is() ))
        HERROR( ERR_INDEXSET, "(TBlockMatrix) mul_vec", "incompatible vector index set" );
    
    //
    // decide upon type
    //

    if ( IS_TYPE( x, TScalarVector ) && IS_TYPE( y, TScalarVector ) )
    {
        const TScalarVector * sx = cptrcast( x, TScalarVector );
        TScalarVector       * sy = ptrcast( y, TScalarVector );

        //
        // get pointer to actual data by comparing local indices and indices
        // in blockmatrix
        //
        
        if ( beta != real(1) )
            y->scale( beta );
        
        if ( alpha == real(0) )
            return;
        
        if ( op == apply_normal )
        {
            const matop_t  sym_op = is_symmetric() ? apply_transposed : apply_adjoint;
            
            for ( uint i = 0; i < block_rows(); i++ )
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const TMatrix  * A_ij = block(i,j);
                    
                    if ( A_ij == nullptr )
                        continue;

                    {
                        // multiply with sub block
                        const TScalarVector  tx = sub_vector( sx, A_ij->col_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->row_is() );
                
                        A_ij->mul_vec( alpha, & tx, real(1), & ty, op );
                    }

                    // in case of symmetry, the same again but transposed
                    if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        const TScalarVector  tx = sub_vector( sx, A_ij->row_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->col_is() );

                        A_ij->mul_vec( alpha, & tx, real(1), & ty, sym_op );
                    }// if
                }// for
        }// if
        else if ( op == apply_transposed )
        {
            for ( uint i = 0; i < block_rows(); i++ )
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const TMatrix  * A_ij = block(i,j);
                    
                    if ( A_ij == nullptr )
                        continue;

                    {
                        // multiply with sub block
                        const TScalarVector  tx = sub_vector( sx, A_ij->row_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->col_is() );
                    
                        A_ij->mul_vec( alpha, & tx, real(1), & ty, op );
                    }

                    // in case of symmetry, the same again but transposed
                    if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        const TScalarVector  tx = sub_vector( sx, A_ij->col_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->row_is() );
                        
                        if ( A_ij->is_complex() && is_hermitian() )
                        {
                            // multiply with conjugated but not transposed matrix B
                            // via: y = conj(B) x  ⇔  conj(y) = B conj(x)
                            TScalarVector  bx( A_ij->col_is() );
                            TScalarVector  by( A_ij->row_is(), real(0) );

                            bx.set_complex( tx.is_complex() );
                            by.set_complex( ty.is_complex() );

                            if ( tx.is_complex() )
                            {
                                B::copy( tx.blas_cvec(), bx.blas_cvec() );
                                B::conj( bx.blas_cvec() );
                            }// if
                            else
                                B::copy( tx.blas_rvec(), bx.blas_rvec() );
                                            
                            A_ij->mul_vec( alpha, & bx, real(1), & by, apply_normal );

                            if ( by.is_complex() )
                            {
                                B::conj( by.blas_cvec() );
                                B::add( complex(1), by.blas_cvec(), ty.blas_cvec() );
                            }// if
                            else 
                                B::add( real(1), by.blas_rvec(), ty.blas_rvec() );
                        }// if
                        else
                        {
                            A_ij->mul_vec( alpha, & tx, real(1), & ty, apply_normal );
                        }// else
                    }// if
                }// for
        }// if
        else if ( op == apply_adjoint )
        {
            for ( uint i = 0; i < block_rows(); i++ )
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const TMatrix  * A_ij = block(i,j);
                    
                    if ( A_ij == nullptr )
                        continue;

                    {
                        // multiply with sub block
                        const TScalarVector  tx = sub_vector( sx, A_ij->row_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->col_is() );
                    
                        A_ij->mul_vec( alpha, & tx, real(1), & ty, op );
                    }

                    // in case of symmetry, the same again but transposed
                    if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        const TScalarVector  tx = sub_vector( sx, A_ij->col_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->row_is() );
                        
                        if ( is_symmetric() && A_ij->is_complex() )
                        {
                            // multiply with conjugated but not transposed matrix B
                            // via: y = conj(B) x <=> conj(y) = B conj(x)
                            TScalarVector  bx( A_ij->col_is() );
                            TScalarVector  by( A_ij->row_is(), real(0) );

                            bx.set_complex( tx.is_complex() );
                            by.set_complex( ty.is_complex() );

                            if ( tx.is_complex() )
                            {
                                B::copy( tx.blas_cvec(), bx.blas_cvec() );
                                B::conj( bx.blas_cvec() );
                            }// if
                            else
                                B::copy( tx.blas_rvec(), bx.blas_rvec() );
                                            
                            A_ij->mul_vec( alpha, & bx, real(1), & by, apply_normal );

                            if ( by.is_complex() )
                            {
                                B::conj( by.blas_cvec() );
                                B::add( complex(1), by.blas_cvec(), ty.blas_cvec() );
                            }// if
                            else 
                                B::add( real(1), by.blas_rvec(), ty.blas_rvec() );
                        }// if
                        else
                        {
                            A_ij->mul_vec( alpha, & tx, real(1), & ty, apply_normal );
                        }// else
                    }// if
                }// for
        }// if
    }// if
    else if ( IS_TYPE( x, TBlockVector ) && IS_TYPE( y, TBlockVector ) )
    {
        const TBlockVector * bx = cptrcast( x, TBlockVector );
        TBlockVector       * by = ptrcast( y, TBlockVector );

        //
        // get pointer to actual data by comparing local indices and indices
        // in blockmatrix
        //
        
        if ( beta != real(1) )
            y->scale( beta );
        
        if ( alpha == real(0) )
            return;
        
        if ( op == apply_normal )
        {
            const matop_t  sym_op = is_symmetric() ? apply_transposed : apply_adjoint;
            
            for ( uint i = 0; i < block_rows(); i++ )
            {
                TVector *  y_i = by->block( i );
                
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const TMatrix *  A_ij = block(i,j);
                    const TVector *  x_j  = bx->block( j );
                    
                    if ( A_ij == nullptr )
                        continue;
                    
                    A_ij->mul_vec( alpha, x_j, real(1), y_i, op );

                    // in case of symmetry, the same again but transposed
                    if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        A_ij->mul_vec( alpha, bx->block( i ), real(1), by->block( j ), sym_op );
                    }// if
                }// for
            }// for
        }// if
        else if ( op == apply_transposed )
        {
            for ( uint i = 0; i < block_rows(); i++ )
            {
                const TVector *  x_i = bx->block( i );
                
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const TMatrix *  A_ij = block(i,j);
                    TVector *        y_j  = by->block( j );
                    
                    if ( A_ij == nullptr )
                        continue;
                    
                    A_ij->mul_vec( alpha, x_i, real(1), y_j, op );

                    // in case of symmetry, the same again but transposed
                    if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        if ( A_ij->is_complex() && is_hermitian() )
                        {
                            // multiply with conjugated but not transposed matrix B
                            // via: y = conj(B) x <=> conj(y) = B conj(x)
                            HERROR( ERR_NOT_IMPL, "", "" );
                        }// if
                        else
                        {
                            A_ij->mul_vec( alpha, bx->block(j), real(1), by->block(i), apply_normal );
                        }// else
                    }// if
                }// for
            }// for
        }// if
        else if ( op == apply_adjoint )
        {
            for ( uint i = 0; i < block_rows(); i++ )
            {
                const TVector *  x_i = bx->block( i );
                
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const TMatrix *  A_ij = block(i,j);
                    TVector *        y_j  = by->block( j );
                    
                    if ( A_ij == nullptr )
                        continue;
                    
                    A_ij->mul_vec( alpha, x_i, real(1), y_j, op );

                    // in case of symmetry, the same again but transposed
                    if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        if ( is_symmetric() && A_ij->is_complex() )
                        {
                            // multiply with conjugated but not transposed matrix B
                            // via: y = conj(B) x <=> conj(y) = B conj(x)
                            HERROR( ERR_NOT_IMPL, "", "" );
                        }// if
                        else
                        {
                            A_ij->mul_vec( alpha, bx->block(j), real(1), by->block(i), apply_normal );
                        }// else
                    }// if
                }// for
            }// for
        }// if
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TBlockMatrix) mul_vec",
                y->typestr() + " = TBlockMatrix * " + x->typestr() );
}

//
// compute this = this + alpha * matrix
//
void
TBlockMatrix::add ( const real       /* alpha */,
                    const TMatrix *  /* B */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
    // return HLIB::add( alpha, B, real(1), this, TTruncAcc( Limits::epsilon<real>() ) );
}

/////////////////////////////////////////////////
//
// BLAS-routines (complex valued)
//

//
// scale matrix by constant factor
//
void
TBlockMatrix::cscale ( const complex f )
{
    if ( f == real(1) )
        return;
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            TMatrix  * A_ij = block(i,j);
            
            if ( A_ij != nullptr )
                A_ij->cscale( f );
        }// for
}

//
// compute this = this + a * matrix
// (matrix must be of compatible type !)
//
void
TBlockMatrix::cadd ( const complex    /* a */,
                     const TMatrix *  /* B */ )
{
    HERROR( ERR_NOT_IMPL, "", "" );
    // HLIB::add< complex >( a, B, real(1), this, TTruncAcc( Limits::epsilon<real>() ) );
}

//
// matrix-vector-multiplication : y = alpha op(A) * x + beta * y
//
void
TBlockMatrix::cmul_vec ( const complex   alpha,
                         const TVector * x,
                         const complex   beta,
                         TVector       * y,
                         const matop_t   op ) const
{
    // if ( ! is_small( this ) )
    // {
    //     HLIB::cmul_vec( procs(), alpha, this, x, beta, y, op );
    //     return;
    // }// if
    
    if ( x == nullptr ) HERROR( ERR_ARG, "(TBlockMatrix) cmul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TBlockMatrix) cmul_vec", "y = nullptr" );

    if (( row_is( op ) != y->is() ) || ( col_is( op ) != x->is() ))
        HERROR( ERR_INDEXSET, "(TBlockMatrix) cmul_vec", "incompatible vector index set" );
    
    //
    // decide upon type
    //

    if ( IS_TYPE( x, TScalarVector ) && IS_TYPE( y, TScalarVector ) )
    {
        const TScalarVector * sx = cptrcast( x, TScalarVector );
        TScalarVector       * sy = ptrcast( y, TScalarVector );

        //
        // get pointer to actual data by comparing local indices and indices
        // in blockmatrix
        //
        
        if ( beta != real(1) )
            y->cscale( beta );
        
        if ( alpha == real(0) )
            return;
        
        if ( op == apply_normal )
        {
            for ( uint i = 0; i < block_rows(); i++ )
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const TMatrix  * A_ij = block(i,j);
                    
                    if ( A_ij == nullptr )
                        continue;

                    {
                        // multiply with sub block
                        const TScalarVector  tx = sub_vector( sx, A_ij->col_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->row_is() );
                
                        A_ij->cmul_vec( alpha, & tx, real(1), & ty, op );
                    }

                    // in case of symmetry, the same again but transposed
                    if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        const TScalarVector  tx = sub_vector( sx, A_ij->row_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->col_is() );

                        if ( is_symmetric() )
                            A_ij->cmul_vec( alpha, & tx, real(1), & ty, apply_transposed );
                        else
                            A_ij->cmul_vec( alpha, & tx, real(1), & ty, apply_adjoint );
                    }// if
                }// for
        }// if
        else if ( op == apply_transposed )
        {
            for ( uint i = 0; i < block_rows(); i++ )
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const TMatrix  * A_ij = block(i,j);
                    
                    if ( A_ij == nullptr )
                        continue;

                    {
                        // multiply with sub block
                        const TScalarVector  tx = sub_vector( sx, A_ij->row_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->col_is() );
                    
                        A_ij->cmul_vec( alpha, & tx, real(1), & ty, op );
                    }

                    // in case of symmetry, the same again but transposed
                    if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        const TScalarVector  tx = sub_vector( sx, A_ij->col_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->row_is() );
                        
                        if ( A_ij->is_complex() && is_hermitian() )
                        {
                            // multiply with conjugated but not transposed matrix B
                            // via: y = conj(B) x  ⇔  conj(y) = B conj(x)
                            TScalarVector  bx( tx.is(), tx.is_complex() );
                            TScalarVector  by( ty.is(), ty.is_complex() );

                            if ( tx.is_complex() )
                            {
                                B::copy( tx.blas_cvec(), bx.blas_cvec() );
                                B::conj( bx.blas_cvec() );
                            }// if
                            else
                                B::copy( tx.blas_rvec(), bx.blas_rvec() );
                                            
                            A_ij->cmul_vec( conj( alpha ), & bx, real(1), & by, apply_normal );

                            if ( by.is_complex() )
                            {
                                B::conj( by.blas_cvec() );
                                B::add( complex(1), by.blas_cvec(), ty.blas_cvec() );
                            }// if
                            else 
                                B::add( real(1), by.blas_rvec(), ty.blas_rvec() );
                        }// if
                        else
                        {
                            A_ij->cmul_vec( alpha, & tx, real(1), & ty, apply_normal );
                        }// else
                    }// if
                }// for
        }// if
        else if ( op == apply_adjoint )
        {
            for ( uint i = 0; i < block_rows(); i++ )
                for ( uint j = 0; j < block_cols(); j++ )
                {
                    const TMatrix  * A_ij = block(i,j);
                    
                    if ( A_ij == nullptr )
                        continue;

                    {
                        // multiply with sub block
                        const TScalarVector  tx = sub_vector( sx, A_ij->row_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->col_is() );
                    
                        A_ij->cmul_vec( alpha, & tx, real(1), & ty, op );
                    }

                    // in case of symmetry, the same again but transposed
                    if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                    {
                        const TScalarVector  tx = sub_vector( sx, A_ij->col_is() );
                        TScalarVector        ty = sub_vector( sy, A_ij->row_is() );
                        
                        if ( is_symmetric() && A_ij->is_complex() )
                        {
                            // multiply with conjugated but not transposed matrix B
                            // via: y = conj(B) x <=> conj(y) = B conj(x)
                            TScalarVector  bx( tx.is(), tx.is_complex() ); 
                            TScalarVector  by( ty.is(), ty.is_complex() );

                            if ( tx.is_complex() )
                            {
                                B::copy( tx.blas_cvec(), bx.blas_cvec() );
                                B::conj( bx.blas_cvec() );
                            }// if
                            else
                                B::copy( tx.blas_rvec(), bx.blas_rvec() );
                                            
                            A_ij->cmul_vec( conj( alpha ), & bx, real(1), & by, apply_normal );

                            if ( by.is_complex() )
                            {
                                B::conj( by.blas_cvec() );
                                B::add( complex(1), by.blas_cvec(), ty.blas_cvec() );
                            }// if
                            else 
                                B::add( real(1), by.blas_rvec(), ty.blas_rvec() );
                        }// if
                        else
                        {
                            A_ij->cmul_vec( alpha, & tx, real(1), & ty, apply_normal );
                        }// else
                    }// if
                }// for
        }// if
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TBlockMatrix) cmul_vec",
                y->typestr() + " = TBlockMatrix * " + x->typestr() );
}

///////////////////////////////////////////////////////////
//
// linear operator mapping
//

//
// same as above but only the dimension of the vector spaces is tested,
// not the corresponding index sets
//
void
TBlockMatrix::apply_add   ( const real                       alpha,
                            const BLAS::Vector< real > &     x,
                            BLAS::Vector< real > &           y,
                            const matop_t                    op ) const
{
    if ( alpha == real(0) )
        return;
        
    if ( op == apply_normal )
    {
        const matop_t  sym_op = is_symmetric() ? apply_transposed : apply_adjoint;
            
        for ( uint i = 0; i < block_rows(); i++ )
        {
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    const B::Vector< real >  tx( x, A_ij->col_is() - col_ofs() );
                    B::Vector< real >        ty( y, A_ij->row_is() - row_ofs() );
                
                    A_ij->apply_add( alpha, tx, ty, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                {
                    const B::Vector< real >  tx( x, A_ij->row_is() - row_ofs() );
                    B::Vector< real >        ty( y, A_ij->col_is() - col_ofs() );

                    A_ij->apply_add( alpha, tx, ty, sym_op );
                }// if
            }// for
        }// for
    }// if
    else if ( op == apply_transposed )
    {
        for ( uint i = 0; i < block_rows(); i++ )
        {
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    const B::Vector< real >  tx( x, A_ij->row_is() - row_ofs() );
                    B::Vector< real >        ty( y, A_ij->col_is() - col_ofs() );
                    
                    A_ij->apply_add( alpha, tx, ty, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                {
                    const B::Vector< real >  tx( x, A_ij->col_is() - col_ofs() );
                    B::Vector< real >        ty( y, A_ij->row_is() - row_ofs() );
                        
                    if ( A_ij->is_complex() && is_hermitian() )
                    {
                        // multiply with conjugated but not transposed matrix B
                        // via: y = conj(B) x  ⇔  conj(y) = B conj(x)
                        HERROR( ERR_NOT_IMPL, "", "" );
                    }// if
                    else
                    {
                        A_ij->apply_add( alpha, tx, ty, apply_normal );
                    }// else
                }// if
            }// for
        }// for
    }// if
    else if ( op == apply_adjoint )
    {
        for ( uint i = 0; i < block_rows(); i++ )
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    const B::Vector< real >  tx( x, A_ij->row_is() - row_ofs() );
                    B::Vector< real >        ty( y, A_ij->col_is() - col_ofs() );
                    
                    A_ij->apply_add( alpha, tx, ty, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr ))
                {
                    const B::Vector< real >  tx( x, A_ij->col_is() - col_ofs() );
                    B::Vector< real >        ty( y, A_ij->row_is() - row_ofs() );
                        
                    if ( is_symmetric() && A_ij->is_complex() )
                    {
                        // multiply with conjugated but not transposed matrix B
                        // via: y = conj(B) x <=> conj(y) = B conj(x)
                        HERROR( ERR_NOT_IMPL, "", "" );
                    }// if
                    else
                    {
                        A_ij->apply_add( alpha, tx, ty, apply_normal );
                    }// else
                }// if
            }// for
    }// if
}

void
TBlockMatrix::apply_add   ( const complex                    alpha,
                            const BLAS::Vector< complex > &  x,
                            BLAS::Vector< complex > &        y,
                            const matop_t                    op ) const
{
    if ( alpha == complex(0) )
        return;
        
    if ( op == apply_normal )
    {
        const matop_t  sym_op = is_symmetric() ? apply_transposed : apply_adjoint;
            
        for ( uint i = 0; i < block_rows(); i++ )
        {
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    const B::Vector< complex >  tx( x, A_ij->col_is() - col_ofs() );
                    B::Vector< complex >        ty( y, A_ij->row_is() - row_ofs() );
                
                    A_ij->apply_add( alpha, tx, ty, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                {
                    const B::Vector< complex >  tx( x, A_ij->row_is() - row_ofs() );
                    B::Vector< complex >        ty( y, A_ij->col_is() - col_ofs() );

                    A_ij->apply_add( alpha, tx, ty, sym_op );
                }// if
            }// for
        }// for
    }// if
    else if ( op == apply_transposed )
    {
        for ( uint i = 0; i < block_rows(); i++ )
        {
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    const B::Vector< complex >  tx( x, A_ij->row_is() - row_ofs() );
                    B::Vector< complex >        ty( y, A_ij->col_is() - col_ofs() );
                    
                    A_ij->apply_add( alpha, tx, ty, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr  ))
                {
                    const B::Vector< complex >  tx( x, A_ij->col_is() - col_ofs() );
                    B::Vector< complex >        ty( y, A_ij->row_is() - row_ofs() );
                        
                    if ( A_ij->is_complex() && is_hermitian() )
                    {
                        // multiply with conjugated but not transposed matrix B
                        // via: y = conj(B) x  ⇔  conj(y) = B conj(x)
                        HERROR( ERR_NOT_IMPL, "", "" );
                    }// if
                    else
                    {
                        A_ij->apply_add( alpha, tx, ty, apply_normal );
                    }// else
                }// if
            }// for
        }// for
    }// if
    else if ( op == apply_adjoint )
    {
        for ( uint i = 0; i < block_rows(); i++ )
            for ( uint j = 0; j < block_cols(); j++ )
            {
                const auto  A_ij = block(i,j);
                    
                if ( A_ij == nullptr )
                    continue;

                {
                    // multiply with sub block
                    const B::Vector< complex >  tx( x, A_ij->row_is() - row_ofs() );
                    B::Vector< complex >        ty( y, A_ij->col_is() - col_ofs() );
                    
                    A_ij->apply_add( alpha, tx, ty, op );
                }

                // in case of symmetry, the same again but transposed
                if ( ! is_nonsym() && ( i > j ) && ( block( j, i ) == nullptr ))
                {
                    const B::Vector< complex >  tx( x, A_ij->col_is() - col_ofs() );
                    B::Vector< complex >        ty( y, A_ij->row_is() - row_ofs() );
                        
                    if ( is_symmetric() && A_ij->is_complex() )
                    {
                        // multiply with conjugated but not transposed matrix B
                        // via: y = conj(B) x <=> conj(y) = B conj(x)
                        HERROR( ERR_NOT_IMPL, "", "" );
                    }// if
                    else
                    {
                        A_ij->apply_add( alpha, tx, ty, apply_normal );
                    }// else
                }// if
            }// for
    }// if
}

///////////////////////////////////////////
//
// block-handling
//

//
// return subblock of matrix corresponding to t
//
TMatrix *
TBlockMatrix::bc_block ( const TBlockCluster * t ) const
{
    if ( t == nullptr )
        HERROR( ERR_ARG, "(TBlockMatrix) bc_block", "given block cluster is nullptr" );
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            const TMatrix  * A_ij = block(i,j);
            
            if (( A_ij != nullptr ) && ( A_ij->cluster() == t ))
                return const_cast< TMatrix * >( A_ij );
        }// for
    
    return nullptr;
}

TMatrix *
TBlockMatrix::bc_block ( const TCluster * tau, const TCluster * sigma ) const
{
    if (( tau == nullptr ) || ( sigma == nullptr ))
        HERROR( ERR_ARG, "(TBlockMatrix) bc_block", "given clusters are nullptr" );
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            const TMatrix  * A_ij = block(i,j);
            
            if ( ( A_ij != nullptr ) && ( A_ij->cluster() != nullptr ) &&
                 ( A_ij->cluster()->rowcl() == tau) &&
                 ( A_ij->cluster()->colcl() == sigma ) )
                return const_cast< TMatrix * >( A_ij );
        }// for
    
    return nullptr;
}

//
// collect matrix-blocks corresponding to leaves
// on processor proc (-1 = all leaves)
//
// void
// TBlockMatrix::collect_leaves ( list< TMatrix * > &  leaf_list ) const
// {
//     for ( uint i = 0; i < block_rows(); i++ )
//         for ( uint j = 0; j < block_cols(); j++ )
//         {
//             const TMatrix  * A_ij = block(i,j);
            
//             if ( A_ij == nullptr )
//                 continue;
            
//             if ( ! IS_TYPE( A_ij, TBlockMatrix ) )
//                 leaf_list.push_back( const_cast< TMatrix * >( A_ij ) );
//             else
//                 cptrcast( A_ij, TBlockMatrix )->collect_leaves( leaf_list );
//         }// for
// }

//
// truncate all rank-blocks in matrix to given accuracy
//
void
TBlockMatrix::truncate ( const TTruncAcc & acc )
{
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
        {
            TMatrix  * A_ij = block(i,j);
            
            if ( A_ij == nullptr )
                continue;

            A_ij->truncate( acc );
        }// for
}
    
//
// switch between complex and real format
//
void
TBlockMatrix::to_real ()
{
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
            if ( block(i,j) != nullptr )
                block(i,j)->set_complex( false, true );
}

void
TBlockMatrix::to_complex ()
{
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
            if ( block(i,j) != nullptr )
                block(i,j)->set_complex( true, true );
}

//
// make value type of this and all sub blocks consistent
//
void
TBlockMatrix::adjust_value_type ()
{
    bool  all_real    = true;
    bool  all_complex = true;

    //
    // test, if all sub blocks are real or complex valued
    //
    
    for ( uint i = 0; i < block_rows(); i++ )
    {
        for ( uint j = 0; j < block_cols(); j++ )
        {
            if ( block(i,j) != nullptr )
            {
                if ( block(i,j)->is_complex() )
                    all_real    = false;
                else if ( block(i,j)->is_real() )
                    all_complex = false;
            }// if
        }// for
    }// for

    // check if no sub blocks are present
    if ( all_real && all_complex )
        return;
    
    // adjust local value type
    if ( all_real && is_complex() )
        set_complex( false );
    
    if ( all_complex && is_real() )
        set_complex( true );

    // and stop if all are equal
    if ( all_real || all_complex )
        return;
    
    //
    // otherwise make all complex (at least one is complex)
    //
    
    for ( uint i = 0; i < block_rows(); i++ )
    {
        for ( uint j = 0; j < block_cols(); j++ )
        {
            if ( block(i,j) != nullptr )
                block(i,j)->set_complex( true );
        }// for
    }// for

    set_complex( true );
}
    
//
// transpose matrix
//
void
TBlockMatrix::transpose ()
{
    if ( is_symmetric() )
    {
        // nothing to do, because transposed matrix is the same
        return;
    }// if
    else if ( is_hermitian() )
    {
        // A = A^H ⇒ A^T = conj(A)
        conjugate();
    }// if
    else
    {
        //
        // transpose sub matrices and transpose internal
        // pointer matrix
        //

        vector< TMatrix * >  T( block_cols() * block_rows() );

        for ( uint  j = 0; j < block_cols(); ++j )
            for ( uint  i = 0; i < block_rows(); ++i )
            {
                if ( block( i, j ) != nullptr )
                    block( i, j )->transpose();

                T[ i * block_cols() + j ] = block( i, j );
            }// for
        
        for ( uint  i = 0; i < block_cols() * block_rows(); ++i )
            _blocks[ i ] = T[ i ];

        swap( _block_rows, _block_cols );
        swap( _rows, _cols );
    }// else

    TMatrix::transpose();
}

//
// conjugate matrix coefficients
//
void
TBlockMatrix::conjugate ()
{
    if ( is_hermitian() )
    {
        return;
    }// if
    else
    {
        for ( uint  j = 0; j < block_cols(); ++j )
            for ( uint  i = 0; i < block_rows(); ++i )
            {
                if ( block( i, j ) != nullptr )
                    block( i, j )->conjugate();
            }// for
    }// else
}

//
// output of matrix
//
void
TBlockMatrix::print ( const uint ofs ) const
{
    for ( uint i = 0; i < ofs; i++ )
        cout << ' ';

    cout << typestr()
         << " ( " << rows() << " x " << cols()
         << ", +" << row_ofs() << "+" << col_ofs()
         << ", " << block_rows() << " x " << block_cols()
         << " )"
         << endl;
    
    for ( uint i = 0; i < block_rows(); i++ )
        for ( uint j = 0; j < block_cols(); j++ )
            if ( block(i,j) != nullptr )
            {
                block(i,j)->print( ofs + 4 );
                if (( i < block_rows()-1 ) || ( j < block_cols()-1 ))
                    cout << endl;
            }// if
}

//
// virtual constructors
//
std::unique_ptr< TMatrix >
TBlockMatrix::copy () const
{
    auto  M = TMatrix::copy();
    auto  B = ptrcast( M.get(), TBlockMatrix );

    B->set_block_struct( block_rows(), block_cols() );

    //
    // copy each of the blocks and add to new matrix
    //

    for ( uint  i = 0; i < nblock_rows(); ++i )
    {
        for ( uint  j = 0; j < nblock_cols(); ++j )
        {
            const TMatrix *  B_ij = block(i,j);

            if ( B_ij != nullptr )
            {
                auto  T = B_ij->copy();
                
                T->set_parent( B );
                B->set_block( i, j, T.release() );
            }// if
        }// for
    }// for
        
    return M;
}

//
// copy matrix wrt. given accuracy and coarsening
//
std::unique_ptr< TMatrix >
TBlockMatrix::copy ( const TTruncAcc & acc, const bool coarsen ) const
{
    auto  M = TMatrix::copy();
    auto  B = ptrcast( M.get(), TBlockMatrix );

    B->set_block_struct( block_rows(), block_cols() );

    //
    // copy each of the blocks and add to new matrix
    //

    for ( uint  i = 0; i < nblock_rows(); ++i )
    {
        for ( uint  j = 0; j < nblock_cols(); ++j )
        {
            const TMatrix *  M_ij = block(i,j);

            if ( M_ij == nullptr )
            {
                B->set_block( i, j, nullptr );
            }// if
            else
            {
                unique_ptr< TMatrix >  T_ij( M_ij->copy( acc( M_ij ), coarsen ) );
                
                // if ( coarsen )
                // {
                //     TCoarsen  mcoarsen;
                
                //     T_ij = unique_ptr< TMatrix >( mcoarsen.coarsen( T_ij.release(), acc( M_ij ) ) );
                // }// if
                
                B->set_block( i, j, T_ij.release() );
            }// else
        }// for
    }// for
    
    return M;
}

//
// return structural copy
//
std::unique_ptr< TMatrix >
TBlockMatrix::copy_struct () const
{
    auto  M = TMatrix::copy_struct();
    auto  B = ptrcast( M.get(), TBlockMatrix );

    B->set_block_struct( block_rows(), block_cols() );

    //
    // copy each of the blocks and add to new matrix
    //

    for ( uint  i = 0; i < nblock_rows(); ++i )
    {
        for ( uint  j = 0; j < nblock_cols(); ++j )
        {
            const TMatrix *  B_ij = block(i,j);

            if ( B_ij != nullptr )
            {
                auto  T = B_ij->copy_struct();
                
                T->set_parent( B );
                B->set_block( i, j, T.release() );
            }// if
        }// for
    }// for
        
    return M;
}

//
// copy matrix into A
//
void
TBlockMatrix::copy_to ( TMatrix * A ) const
{
    TMatrix::copy_to( A );
    
    if ( ! IS_TYPE( A, TBlockMatrix ) )
    {
        HERROR( ERR_NOT_IMPL, "", "" );
        // convert( this, A, acc_exact );
        // return;
    }// if

    TBlockMatrix  * b_A = ptrcast( A, TBlockMatrix );
    
    //
    // go over all blocks, get assoc. blocks in A and
    // copy them one by one
    //
    
    for ( uint  i = 0; i < nblock_rows(); ++i )
    {
        for ( uint  j = 0; j < nblock_cols(); ++j )
        {
            if ( block(i,j) == nullptr )
                continue;
            
            const TMatrix * B = block(i,j);
            TMatrix       * C = b_A->block(i,j);
            
            if ( C != nullptr )
                B->copy_to( C );
            else
                HERROR( ERR_MAT_STRUCT, "(TBlockMatrix) copy_to",
                        "could not find corresponding block in A" );
        }// for
    }// for
}

void
TBlockMatrix::copy_to ( TMatrix         * A,
                        const TTruncAcc & acc,
                        const bool        coarsen ) const
{
    TMatrix::copy_to( A );
    
    if ( ! IS_TYPE( A, TBlockMatrix ) )
    {
        HERROR( ERR_NOT_IMPL, "", "" );

        // convert( this, A, acc );
        // return;
    }// if

    TBlockMatrix  * b_A = ptrcast( A, TBlockMatrix );
    
    //
    // go over all blocks, get assoc. blocks in A and
    // copy them one by one
    //

    for ( uint  i = 0; i < nblock_rows(); ++i )
    {
        for ( uint  j = 0; j < nblock_cols(); ++j )
        {
            if ( block(i,j) == nullptr )
                continue;
            
            const TMatrix * B_ij = block(i,j);
            TMatrix       * C_ij = b_A->block(i,j);
            
            if ( C_ij == nullptr )
                HERROR( ERR_MAT_STRUCT, "(TBlockMatrix) copy_to",
                        "could not find corresponding block in A" );
            
            if ( IS_TYPE( B_ij, TBlockMatrix ) )
                cptrcast( B_ij, TBlockMatrix )->copy_to( C_ij, acc( B_ij ), coarsen );
            else if ( IS_TYPE( B_ij, TRkMatrix ) )
                cptrcast( B_ij, TRkMatrix )->copy_to( C_ij, acc( B_ij ) );
            else
                B_ij->copy_to( C_ij );
            
            // if ( coarsen && ( acc.rank() == 0 ))
            // {
            //     TCoarsen  mcoarsen;
            
            //     C_ij = mcoarsen.coarsen( C_ij, acc( B_ij ) );
            // }// if
            
            b_A->set_block( i, j, C_ij );
        }// for
    }// for
}

//
// copy complete structural information from given matrix
//
void
TBlockMatrix::copy_struct_from ( const TMatrix * M )
{
    TMatrix::copy_struct_from( M );

    if ( IS_TYPE( M, TBlockMatrix ) )
    {
        const TBlockMatrix * B = cptrcast( M, TBlockMatrix );

        set_block_struct( B->block_rows(), B->block_cols() );
    }// if
}
    
//
// return size in bytes used by this object
//
size_t
TBlockMatrix::byte_size () const
{
    size_t  size = 0;
    
    size += TMatrix::byte_size();
    size += 2*sizeof(uint) + _block_rows*_block_cols*sizeof(TMatrix*);

    //
    // compute memory consumption of sub blocks
    //

    for ( uint i = 0; i < block_rows(); i++ )
    {
        for ( uint j = 0; j < block_cols(); j++ )
        {
            const TMatrix  * A_ij = block(i,j);
            
            if ( A_ij != nullptr )
            {
                size += ulong(A_ij->byte_size());
            }// if
        }// for
    }// for

    return size;
}

//
// serialisation
//

void
TBlockMatrix::read ( TByteStream & s )
{
    TMatrix::read( s );

    size_t  n, m;

    s.get( n );
    s.get( m );

    set_size( n, m );
    
    uint  bn, bm;

    s.get( bn );
    s.get( bm );

    // compare with current setting
    if (( _block_rows != bn ) || ( _block_cols != bm ))
        HERROR( ERR_MAT_STRUCT, "(TBlockMatrix) read", "wrong block-structure in stream" );

    vector< char >  present( _block_rows * _block_cols );
    
    s.get( present.data(), _block_rows * _block_cols * sizeof(char) );
    
    // iterate over submatrices and read stream
    for ( uint i = 0; i < _block_rows*_block_cols; i++ )
    {
        if ( present[i] != 0 )
        {
            if ( _blocks[i] != nullptr )
                _blocks[i]->read( s );
            else
                HERROR( ERR_CONSISTENCY, "(TBlockMatrix) read", "nullptr submatrix but matrix in stream" );
        }// if
        else
        {
            if ( _blocks[i] != nullptr )
                HERROR( ERR_CONSISTENCY, "(TBlockMatrix) read", "non-nullptr submatrix but no matrix in stream" );
        }// else
    }// for
}

void
TBlockMatrix::build ( TByteStream & s )
{
    TMatrix::build( s );

    size_t  n, m;

    s.get( n );
    s.get( m );

    set_size( n, m );
    
    uint  bn, bm;
    
    s.get( bn );
    s.get( bm );

    set_block_struct( bn, bm );

    //
    // build submatrices
    //

    vector< char >  present( _block_rows * _block_cols );
    TBSHMBuilder    builder;
    
    s.get( present.data(), _block_rows * _block_cols * sizeof(char) );
    
    for ( uint i = 0; i < _block_rows*_block_cols; i++ )
    {
        if ( present[i] != 0 )
        {
            auto  B_ij = builder.build( s );
            
            if ( B_ij.get() == nullptr )
                HWARNING( "in (TBlockMatrix) build : submatrix is nullptr" );

            _blocks[i] = B_ij.release();
        }// if
        else
            _blocks[i] = nullptr;
    }// for
}

void
TBlockMatrix::write ( TByteStream & s ) const
{
    TMatrix::write( s );

    s.put( _rows );
    s.put( _cols );

    s.put( _block_rows );
    s.put( _block_cols );

    // write boolean array indicating non-nullptr submatrices
    vector< char >  present( _block_rows * _block_cols );
    
    for ( uint i = 0; i < _block_rows*_block_cols; i++ )
    {
        if ( _blocks[i] != nullptr )
            present[i] = 1;
        else
            present[i] = 0;
    }// for

    s.put( present.data(), _block_rows*_block_cols*sizeof(char) );
    
    // iterate over submatrices and read stream
    for ( uint i = 0; i < _block_rows*_block_cols; i++ )
    {
        if ( _blocks[i] != nullptr )
            _blocks[i]->write( s );
    }// for
}

//
// returns size of object in bytestream
//
size_t
TBlockMatrix::bs_size () const
{
    size_t  size = ( TMatrix::bs_size() +
                     sizeof(_rows) + sizeof(_cols) +
                     sizeof(_block_rows) + sizeof(_block_cols) +
                     _block_rows * _block_cols * sizeof(char) );
        
    for ( uint i = 0; i < _block_rows*_block_cols; i++ )
    {
        if ( _blocks[i] != nullptr )
            size += _blocks[i]->bs_size();
    }// for

    return size;
}

//
// test data for invalid values, e.g. INF and NAN
//
void
TBlockMatrix::check_data () const
{
    for ( uint  i = 0; i < block_rows(); ++i )
    {
        for ( uint  j = 0; j < block_cols(); ++j )
        {
            if ( block(i,j) != nullptr )
                block(i,j)->check_data();
        }// for
    }// for
}

}// namespace
