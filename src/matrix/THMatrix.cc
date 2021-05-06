//
// Project     : HLib
// File        : THMatrix.cc
// Description : class for H-matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/matrix/structure.hh"
#include "hpro/matrix/TMatrixHierarchy.hh"

#include "hpro/matrix/THMatrix.hh"

namespace HLIB
{

using std::unique_ptr;

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//
// model independent H-Matrix implementation
//
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////
//
// constructor
//

THMatrix::THMatrix ( const TBlockCluster * bct )
        : TBlockMatrix( bct )
        , _max_id( 0 )
{}

THMatrix::THMatrix ( const TBlockIndexSet &  bis )
        : TBlockMatrix( bis )
        , _max_id( 0 )
{
}

/////////////////////////////////////////////////
//
// access internal data
//

//
// access single matrix coefficient
//
real
THMatrix::entry ( const idx_t row, const idx_t col ) const
{
    idx_t  srow = row;
    idx_t  scol = col;
    
    // translate index pair to lower half in case of symmetry
    if ( ! is_nonsym() && ( srow < scol ))
        std::swap( srow, scol );

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

    return 0.0;
}

const complex
THMatrix::centry ( const idx_t row, const idx_t col ) const
{
    idx_t  srow = row;
    idx_t  scol = col;
    bool   swap_idx = false;
    
    // translate index pair to lower half in case of symmetry
    if ( ! is_nonsym() && ( srow < scol ))
    {
        std::swap( srow, scol );
        swap_idx = true;
    }// if
    
    const idx_t rofs = row_ofs();
    const idx_t cofs = col_ofs();
    
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


    return 0.0;
}

/////////////////////////////////////////////////
//
// BLAS-routines
//

//
// matrix-vector-mult.
//
void
THMatrix::mul_vec ( const real      alpha,
                    const TVector * x,
                    const real      beta,
                    TVector       * y,
                    const matop_t   op ) const
{
    // TBlockMatrix::mul_vec( alpha, x, beta, y, op );
    // return;
    
    if ( x == nullptr ) HERROR( ERR_ARG, "(THMatrix) mul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(THMatrix) mul_vec", "y = nullptr" );

    //
    // multiply
    //

    TBlockMatrix::mul_vec( alpha, x, beta, y, op );
}

/////////////////////////////////////////////////
//
// BLAS-routines (complex valued)
//

//
// matrix-vector-multiplication : y = alpha op(A) * x + beta * y
//
void
THMatrix::cmul_vec ( const complex   alpha,
                     const TVector * x,
                     const complex   beta,
                     TVector       * y,
                     const matop_t   op ) const
{
    // TBlockMatrix::cmul_vec( alpha, x, beta, y, op );
    // return;
    
    if ( x == nullptr ) HERROR( ERR_ARG, "(THMatrix) cmul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(THMatrix) cmul_vec", "y = nullptr" );

    //
    // multiply
    //
    
    TBlockMatrix::cmul_vec( alpha, x, beta, y, op );
}
    
/////////////////////////////////////////////////
//
// misc.
//

//
// return appropriate vector-types for matrix
//
auto
THMatrix::row_vector () const -> unique_ptr< TVector >
{
    return std::make_unique< TScalarVector >( rows(), row_ofs(), is_complex() );
}

auto
THMatrix::col_vector () const -> unique_ptr< TVector >
{
    return std::make_unique< TScalarVector >( cols(), col_ofs(), is_complex() );
}
    
//
// transpose matrix
//
void
THMatrix::transpose ()
{
    TBlockMatrix::transpose();

    swap( _row_perm_e2i, _col_perm_e2i );
    swap( _row_perm_i2e, _col_perm_i2e );

    comp_min_max_idx();
}
    
//
// return size in bytes used by this object
//
size_t
THMatrix::byte_size () const
{
    size_t  size = TBlockMatrix::byte_size();

    size += _row_perm_e2i.byte_size();
    size += _col_perm_e2i.byte_size();
    size += _row_perm_i2e.byte_size();
    size += _col_perm_i2e.byte_size();
    
#if 0
    size += _leaves.byte_size();
#endif

    return size;
}

//
// serialisation
//

void
THMatrix::read  ( TByteStream & s )
{
    TBlockMatrix::read( s );
}

void
THMatrix::build ( TByteStream & s )
{
    TBlockMatrix::build( s );
}

void
THMatrix::write ( TByteStream & s ) const
{
    TBlockMatrix::write( s );
}

//
// returns size of object in bytestream
//
size_t
THMatrix::bs_size () const
{
    return (TBlockMatrix::bs_size() +
            _row_perm_e2i.byte_size() +
            _col_perm_e2i.byte_size() +
            _row_perm_i2e.byte_size() +
            _col_perm_i2e.byte_size() );
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//
// model dependent H-Matrix implementation (BSP)
//
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//
// compute min/max tau/sigma indices
//
void
THMatrix::comp_min_max_idx ()
{
    //
    // set up matrix hierarchy
    //

    set_hierarchy_data();
    _max_id = HLIB::max_id( this );
}

/////////////////////////////////////////////////
//
// misc.
//

//
// virtual constructor
//
std::unique_ptr< TMatrix >
THMatrix::copy () const
{
    auto        M = TBlockMatrix::copy();
    THMatrix *  H = ptrcast( M.get(), THMatrix );

    H->comp_min_max_idx();
    H->set_row_perm( _row_perm_e2i, _row_perm_i2e );
    H->set_col_perm( _col_perm_e2i, _col_perm_i2e );
    
    return M;
}

//
// copy matrix wrt. given accuracy and coarsening
//
std::unique_ptr< TMatrix >
THMatrix::copy ( const TTruncAcc & acc, const bool coarsen ) const
{
    auto        M = TBlockMatrix::copy( acc, coarsen );
    THMatrix *  H = ptrcast( M.get(), THMatrix );

    H->comp_min_max_idx();
    H->set_row_perm( _row_perm_e2i, _row_perm_i2e );
    H->set_col_perm( _col_perm_e2i, _col_perm_i2e );
    
    return M;
}

//
// return structural copy
//
std::unique_ptr< TMatrix >
THMatrix::copy_struct () const
{
    auto        M = TBlockMatrix::copy_struct();
    THMatrix *  H = ptrcast( M.get(), THMatrix );

    H->comp_min_max_idx();
    H->set_row_perm( _row_perm_e2i, _row_perm_i2e );
    H->set_col_perm( _col_perm_e2i, _col_perm_i2e );
    
    return M;
}

//
// copy matrix into A
//
void
THMatrix::copy_to ( TMatrix * A ) const
{
    TBlockMatrix::copy_to( A );

    if ( ! IS_TYPE( A, THMatrix ) )
        HERROR( ERR_MAT_TYPE, "(THMatrix) copy_to", A->typestr() );

    THMatrix * H = ptrcast( A, THMatrix );
    
    H->comp_min_max_idx();
    H->set_row_perm( _row_perm_e2i, _row_perm_i2e );
    H->set_col_perm( _col_perm_e2i, _col_perm_i2e );
}

void
THMatrix::copy_to ( TMatrix * A, const TTruncAcc & acc, const bool coarsen ) const
{
    TBlockMatrix::copy_to( A, acc, coarsen );

    if ( ! IS_TYPE( A, THMatrix ) )
        HERROR( ERR_MAT_TYPE, "(THMatrix) copy_to", A->typestr() );

    THMatrix * H = ptrcast( A, THMatrix );
    
    H->comp_min_max_idx();
    H->set_row_perm( _row_perm_e2i, _row_perm_i2e );
    H->set_col_perm( _col_perm_e2i, _col_perm_i2e );
}

//
// copy complete structural information from given matrix
//
void
THMatrix::copy_struct_from ( const TMatrix * M )
{
    TBlockMatrix::copy_struct_from( M );

    if ( IS_TYPE( M, THMatrix ) )
    {
        const THMatrix * H = cptrcast( M, THMatrix );

        comp_min_max_idx();
        set_row_perm( H->row_perm_e2i(), H->row_perm_i2e() );
        set_col_perm( H->col_perm_e2i(), H->col_perm_i2e() );
    }// if
}
    
//
// collect leaves
//
bool
size_sort ( TMatrix * M1, TMatrix * M2 )
{
    return std::max( M1->rows(), M1->cols() ) > std::max( M2->rows(), M2->cols() );
}

// void
// THMatrix::build_leaf_list ()
// {
//     std::list< TMatrix * >  tmp;
    
//     TBlockMatrix::collect_leaves( tmp );

//     _leaves.resize( tmp.size() );

//     size_t  i = 0;
//     for ( auto M : tmp ) _leaves[i++] = M;

//     sort( _leaves.begin(), _leaves.end(), size_sort );
// }

}// namespace
