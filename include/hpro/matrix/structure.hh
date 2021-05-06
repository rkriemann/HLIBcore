#ifndef __HLIB_STRUCTURE_HH
#define __HLIB_STRUCTURE_HH
//
// Project     : HLib
// File        : structure.hh
// Description : test functions for various matrix structures
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/matrix/TSparseMatrix.hh"
#include "hpro/matrix/THMatrix.hh"
#include "hpro/matrix/TGhostMatrix.hh"

namespace HLIB
{

//!
//! return true, if given op is a matrix
//!
inline bool is_matrix       ( const TLinearOperator *  O ) noexcept { return IS_TYPE( O, TMatrix ); }

//!
//! return true, if matrix is an H-matrix
//!
inline bool is_hmat         ( const TMatrix *  A ) noexcept { return IS_TYPE( A, THMatrix ); }
inline bool is_hmat         ( const TMatrix &  A ) noexcept { return IS_TYPE( &A, THMatrix ); }

//!
//! return true, if matrix is a block matrix
//!
inline bool is_blocked      ( const TMatrix *  A ) noexcept { return (A != nullptr) && A->is_blocked(); }
inline bool is_blocked      ( const TMatrix &  A ) noexcept { return A.is_blocked(); }

//!
//! return true, if all matrices are blocked
//!
inline bool is_blocked_all  ( const TMatrix *  A ) noexcept { return is_blocked( A ); }
inline bool is_blocked_all  ( const TMatrix &  A ) noexcept { return is_blocked( A ); }

template < typename... T >
bool is_blocked_all  ( const TMatrix *  A, T...  mats ) { return is_blocked( A ) && is_blocked_all( mats... ); }
template < typename... T >
bool is_blocked_all  ( const TMatrix &  A, T...  mats ) { return is_blocked( A ) && is_blocked_all( mats... ); }


//!
//! return true, if matrix is a dense matrix
//!
inline bool is_dense        ( const TMatrix *  A ) noexcept { return IS_TYPE( A, TDenseMatrix ); }

//!
//! return true, if matrix is a lowrank matrix
//!
inline bool is_lowrank      ( const TMatrix *  A ) noexcept { return IS_TYPE( A, TRkMatrix ); }

//!
//! return true, if matrix is a sparse matrix
//!
inline bool is_sparse       ( const TMatrix *  A ) noexcept { return IS_TYPE( A, TSparseMatrix ); }

//!
//! return true, if matrix is a ghost matrix
//!
inline bool is_ghost        ( const TMatrix *  A ) noexcept { return IS_TYPE( A, TGhostMatrix ); }

//!
//! return true if \a A is block diagonal
//!
bool is_diag         ( const TMatrix *  A );

//!
//! return true if \a A has domain-decomposition format, e.g.
//! only block diagonal and non-zero last block row/column
//!
bool is_dd           ( const TMatrix *  A );

//!
//! return true if \a A is block matrix with flat hierarchy
//!
bool is_flat         ( const TMatrix *  A );

//!
//! return true if \a A is a diagonal block
//!
inline bool is_on_diag      ( const TMatrix *  A ) noexcept { return A->row_is() == A->col_is(); }

//!
//! return true if \a A is in lower left part of matrix
//!
inline bool is_in_lower_left ( const TMatrix *  A ) noexcept { return A->row_is().is_right_or_equal_to( A->col_is() ); }

//!
//! return true if \a A is in lower left part of matrix
//!
inline bool is_strictly_in_lower_left ( const TMatrix *  A ) noexcept { return A->row_is().is_strictly_right_of( A->col_is() ); }

//!
//! return true if \a A is in upper right part of matrix
//!
inline bool is_in_upper_right ( const TMatrix *  A ) noexcept { return A->col_is().is_right_or_equal_to( A->row_is() ); }

//!
//! return true if \a A is in upper right part of matrix
//!
inline bool is_strictly_in_upper_right ( const TMatrix *  A ) noexcept { return A->col_is().is_strictly_right_of( A->row_is() ); }

//!
//! return true if \a A is lower left block triangular matrix
//!
bool is_lower_left   ( const TMatrix *  A );

//!
//! return true if \a A is upper right block triangular matrix
//!
bool is_upper_right  ( const TMatrix *  A );

//
// return true if \a A is corresponding to leaf block
//
inline bool is_leaf         ( const TMatrix *  A ) noexcept { return ! A->is_blocked(); }

//!
//! return true, if matrix is considered too small for parallel treatement
//!
inline bool is_small        ( const TMatrix *  A ) noexcept { return ( std::min( A->rows(), A->cols() ) <= CFG::Arith::max_seq_size ); }

//!
//! return true, if any/all matrix is considered too small for parallel treatement
//!
inline bool is_small_any    ( const TMatrix *  A ) noexcept { return is_small( A ); }

template < typename... T >
bool is_small_any    ( const TMatrix *  A, T...  mats ) { return is_small( A ) || is_small_any( mats... ); }

//!
//! return true if matrix is zero, e.g. all entries zero
//!
bool
is_zero         ( const TMatrix *  A );

//!
//! return maximal ID in all submatrices
//!
int
max_id          ( const TMatrix *  A );

//
// return first level on diagonal on which leaf blocks are found
//
uint
get_first_diag_leaf_level ( const TMatrix *  A );

//!
//! return matrix with rowis( A_row ) × colis( A_col )
//!
TMatrix *
product_block   ( const TMatrix *  A_row,
                  const TMatrix *  A_col );

//!
//! return matrix with rowis( op( A_row ) ) × colis( op( A_col ) )
//!
TMatrix *
product_block   ( const matop_t    op_row,
                  const TMatrix *  A_row,
                  const matop_t    op_col,
                  const TMatrix *  A_col );

//!
//! return sub block of M with given id
//!
const TMatrix *
get_block       ( const TMatrix *  M,
                  const int        id );

//!
//! return sub block of M with given block index set
//!
const TMatrix *
get_block       ( const TMatrix *         M,
                  const TBlockIndexSet &  bis );

//!
//! return number of sub blocks in matrix (including inner blocks)
//!
size_t
get_nblocks        ( const TMatrix *   M );

//!
//! return number of leaf blocks in matrix
//!
size_t
get_nblocks_leaf   ( const TMatrix *   M );

//!
//! return number of low-rank blocks in matrix
//!
size_t
get_nblocks_lr     ( const TMatrix *   M );

//!
//! return number of dense blocks in matrix
//!
size_t
get_nblocks_dense  ( const TMatrix *   M );

}// namespace HLIB

#endif  // __HLIB_STRUCTURE_HH

