//
// Project     : HLib
// File        : TSparseMatrix.cc
// Description : class for a sparse matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <iostream>

#include "hpro/matrix/TSparseMatrix.hh"

namespace HLIB
{

using std::vector;
using std::unique_ptr;
using std::make_unique;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// local functions
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

namespace
{

//
// return number of non-zeroes in lower left part plus diagonal
//
size_t
lower_left_nonzeroes ( const TSparseMatrix *  S )
{
    if ( S == nullptr )
        return 0;

    size_t  llnnz = 0;

    for ( idx_t  row = 0; row < idx_t(S->rows()); row++ )
    {
        const idx_t  lb = S->rowptr(row);
        const idx_t  ub = S->rowptr(row+1);

        for ( idx_t  j = lb; j < ub; ++j )
            if ( S->colind(j) <= row )
                llnnz++;
    }// for

    return llnnz;
}

}// namespace anonymous

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TSparseMatrix
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

///////////////////////////////////////////
//
// constructor and destructor
//

TSparseMatrix::~TSparseMatrix ()
{
    _rowptr.resize( 0 );
    _colind.resize( 0 );
    _rcoeff.resize( 0 );
    _ccoeff.resize( 0 );
}

/////////////////////////////////////////////////
//
// access data
//

//
// set cluster of matrix
//
void
TSparseMatrix::set_cluster ( const TBlockCluster * c )
{
    if ( _colind.size() > 0 )
        HERROR( ERR_CONSISTENCY, "(TSparseMatrix) set_cluster",
                "changing sparse matrix in CRS format" );
    
    if ( c != nullptr )
    {
        if (( c->rowcl() != nullptr ) && ( c->colcl() != nullptr ))
            set_size( c->rowcl()->size(), c->colcl()->size() );
    }// if
}

//
// directly set dimension of matrix
//
void
TSparseMatrix::set_size ( const size_t  n,
                          const size_t  m )
{
    if (( _rows != n ) || ( _cols != m ))
    {
        if ( _colind.size() > 0 )
            HERROR( ERR_CONSISTENCY, "(TSparseMatrix) set_size",
                    "can not change non-zero sparse matrix" );
        
        _rows = n;
        _cols = m;

        // initialise with no non zeroes
        init( 0 );
    }// if
}

//
// switch between complex and real format
//
void
TSparseMatrix::to_real ()
{
    if ( ! is_complex() )
        return;

    for ( size_t  i = 0; i < _ccoeff.size(); ++i )
        if ( std::imag( _ccoeff[i] ) != real(0) )
        {
            TMatrix::set_complex( true );
            HERROR( ERR_COMPLEX, "(TSparseMatrix) to_real", "matrix has imaginary part" );
        }// if
    
    _rcoeff.resize( _ccoeff.size() );
    
    for ( size_t  i = 0; i < _ccoeff.size(); ++i )
        _rcoeff[i] = std::real( _ccoeff[i] );
    
    _ccoeff.resize( 0 );
}

void
TSparseMatrix::to_complex ()
{
    if ( is_complex() )
        return;
    
    _ccoeff.resize( _rcoeff.size() );
        
    for ( size_t  i = 0; i < _rcoeff.size(); ++i )
        _ccoeff[i] = _rcoeff[i];
    
    _rcoeff.resize( 0 );
}

//
// return a_ij
//
real
TSparseMatrix::entry ( const idx_t  ai,
                       const idx_t  aj ) const
{
    if ( ! ( row_is().is_in( ai ) && col_is().is_in( aj ) ) )
        HERROR( ERR_MAT_SIZE, "(TSparseMatrix) entry",
                to_string( "(%d,%d) not in indexset", ai, aj ) );
        
    if ( is_complex() )
        HERROR( ERR_COMPLEX, "(TSparseMatrix) entry", "matrix is complex valued" );

    const idx_t  lb = _rowptr[ai];
    const idx_t  ub = _rowptr[ai+1];
    
    for ( idx_t  pos = lb; pos < ub; pos++ )
    {
        if ( _colind[pos] == aj )
            return _rcoeff[pos];
    }// for

    return real(0);
}

const complex
TSparseMatrix::centry ( const idx_t  ai,
                        const idx_t  aj ) const
{
    if ( ! ( row_is().is_in( ai ) && col_is().is_in( aj ) ) )
        HERROR( ERR_MAT_SIZE, "(TSparseMatrix) centry",
                to_string( "(%d,%d) not in indexset", ai, aj ) );
        
    const idx_t  lb = _rowptr[ai];
    const idx_t  ub = _rowptr[ai+1];
    
    for ( idx_t  pos = lb; pos < ub; pos++ )
    {
        if ( _colind[pos] == aj )
        {
            if ( is_complex() ) return _ccoeff[pos];
            else                return _rcoeff[pos];
        }// if
    }// for

    return real(0);
}

//
// set entry a_ij (add if not present)
//
void
TSparseMatrix::set_entry ( const idx_t  ai,
                           const idx_t  aj,
                           const real   c )
{
    if ( ! ( row_is().is_in( ai ) && col_is().is_in( aj ) ) )
        HERROR( ERR_MAT_SIZE, "(TSparseMatrix) set_entry",
                to_string( "(%d,%d) not in indexset", ai, aj ) );
    
    const idx_t  lb = _rowptr[ai];
    const idx_t  ub = _rowptr[ai+1];

    for ( idx_t  pos = lb; pos < ub; pos++ )
    {
        if ( _colind[pos] == aj )
        {
            if ( is_complex() ) _ccoeff[pos] = c;
            else                _rcoeff[pos]  = c;
            
            return;
        }// if
    }// for

    HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) set_entry",
            to_string( "(%d,%d) not in sparsity pattern", ai, aj ) );
}

void
TSparseMatrix::set_entry ( const idx_t    ai,
                           const idx_t    aj,
                           const complex  c )
{
    if ( ! ( row_is().is_in( ai ) && col_is().is_in( aj ) ) )
        HERROR( ERR_MAT_SIZE, "(TSparseMatrix) set_entry",
                to_string( "(%d,%d) not in indexset", ai, aj ) );
    
    if ( ! is_complex() && ( std::imag( c ) != real(0) ))
        HERROR( ERR_REAL, "(TSparseMatrix) set_entry", "matrix is real valued" );
    
    const idx_t  lb = _rowptr[ai];
    const idx_t  ub = _rowptr[ai+1];

    for ( idx_t  pos = lb; pos < ub; pos++ )
    {
        if ( _colind[pos] == aj )
        {
            if ( is_complex() ) _ccoeff[pos] = c;
            else                _rcoeff[pos]  = std::real( c );
            
            return;
        }// if
    }// for

    HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) set_entry",
            to_string( "(%d,%d) not in sparsity pattern", ai, aj ) );
}

//
// add to entry a_ij if present
//
void
TSparseMatrix::add_entry ( const idx_t  ai,
                           const idx_t  aj,
                           const real   c )
{
    if ( ! ( row_is().is_in( ai ) && col_is().is_in( aj ) ) )
        HERROR( ERR_MAT_SIZE, "(TSparseMatrix) add_entry",
                to_string( "(%d,%d) not in indexset", ai, aj ) );
    
    const idx_t  lb = _rowptr[ai];
    const idx_t  ub = _rowptr[ai+1];

    for ( idx_t  pos = lb; pos < ub; pos++ )
    {
        if ( _colind[pos] == aj )
        {
            if ( is_complex() ) _ccoeff[pos] += c;
            else                _rcoeff[pos]  += c;
            
            return;
        }// if
    }// for

    HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) add_entry",
            to_string( "(%d,%d) not in sparsity pattern", ai, aj ) );
}

void
TSparseMatrix::add_entry ( const idx_t    ai,
                           const idx_t    aj,
                           const complex  c )
{
    if ( ! ( row_is().is_in( ai ) && col_is().is_in( aj ) ) )
        HERROR( ERR_MAT_SIZE, "(TSparseMatrix) add_entry",
                to_string( "(%d,%d) not in indexset", ai, aj ) );
    
    if ( ! is_complex() && ( std::imag( c ) != real(0) ))
        HERROR( ERR_REAL, "(TSparseMatrix) add_entry", "matrix is real valued" );
    
    const idx_t  lb = _rowptr[ai];
    const idx_t  ub = _rowptr[ai+1];

    for ( idx_t  pos = lb; pos < ub; pos++ )
    {
        if ( _colind[pos] == aj )
        {
            if ( is_complex() ) _ccoeff[pos] += c;
            else                _rcoeff[pos]  += std::real( c );
            
            return;
        }// if
    }// for

    HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) add_entry",
            to_string( "(%d,%d) not in sparsity pattern", ai, aj ) );
}

//
// add to entry a_ij (insert if not present)
//
bool
TSparseMatrix::has_entry ( const idx_t  ai,
                           const idx_t  aj ) const
{
    if ( ! ( row_is().is_in( ai ) && col_is().is_in( aj ) ) )
        HERROR( ERR_MAT_SIZE, "(TSparseMatrix) has_entry",
                to_string( "(%d,%d) not in indexset", ai, aj ) );
    
    const idx_t  lb = _rowptr[ai];
    const idx_t  ub = _rowptr[ai+1];

    for ( idx_t  pos = lb; pos < ub; pos++ )
    {
        if ( _colind[pos] == aj )
            return true;
    }// for
    
    return false;
}

//
// sort list per row wrt column
//
void
TSparseMatrix::sort_entries ()
{
    //
    // sort CRS data structure
    //
    
    for ( size_t i = 0; i < rows(); i++ )
    {
        const idx_t   lb = _rowptr[i];
        const idx_t   ub = _rowptr[i+1];
        
        for ( idx_t  j = lb; j < ub; j++ )
        {
            bool changed = false;
            
            for ( idx_t  l = lb; l < (ub-j+lb-1); l++ )
                if ( _colind[l] > _colind[l+1] )
                {
                    std::swap( _colind[l], _colind[l+1] );
                    
                    if ( is_complex() )
                        std::swap( _ccoeff[l], _ccoeff[l+1] );
                    else
                        std::swap( _rcoeff[l], _rcoeff[l+1] );
                    
                    changed = true;
                }// if
            
            if ( ! changed )
                break;
        }// for
    }// for
}

//
// initialise CRS data for n entries
//
void
TSparseMatrix::init ( const size_t nnz )
{
    if ( is_complex() ) _ccoeff.resize( nnz );
    else                _rcoeff.resize( nnz );
    
    _colind.resize( nnz );
    _rowptr.resize( rows() + 1 );

    // nullify CRS data
    for ( size_t  i = 0; i <= rows(); i++ )
        _rowptr[i] = 0;
}

//
// import CRS data
//
template <typename T_idx, typename T_val>
void
TSparseMatrix::import_crs ( const size_t   nrows,
                            const size_t   ncols,
                            const size_t   nnz,
                            const T_idx *  arowptr,
                            const T_idx *  acolind,
                            const T_val *  acoeffs )
{
    const bool  is_complex_val = is_complex_type< T_val >::value;

    init( 0 );

    set_complex( is_complex_val );

    _rows = nrows;
    _cols = ncols;

    _rowptr.resize( nrows+1 );
    _colind.resize( nnz );

    for ( size_t i = 0; i <= nrows; ++i ) _rowptr[i] = arowptr[i];
    for ( size_t i = 0; i <  nnz;   ++i ) _colind[i] = acolind[i];

    if ( is_complex_val )
    {
        _ccoeff.resize( nnz );
        _rcoeff.resize( 0 );

        for ( size_t i = 0; i < nnz; ++i )
            _ccoeff[i]  = acoeffs[i];
    }// if
    else
    {
        _rcoeff.resize( nnz );
        _ccoeff.resize( 0 );

        for ( size_t i = 0; i < nnz; ++i )
            _rcoeff[i]  = std::real( acoeffs[i] ); // use std::real(·) to prevent compiler error
                                           // in complex valued case
    }// else
}

//
// import CCS data
//
template <typename T_idx, typename T_val>
void
TSparseMatrix::import_ccs ( const size_t   nrows,
                            const size_t   ncols,
                            const size_t   nnz,
                            const T_idx *  acolptr,
                            const T_idx *  arowind,
                            const T_val *  acoeffs )
{
    const bool  is_complex_val = is_complex_type< T_val >::value;

    init( 0 );

    set_complex( is_complex_val );

    _rows = nrows;
    _cols = ncols;

    _rowptr.resize( nrows+1 );
    _colind.resize( nnz );

    if ( is_complex_val )
    {
        _rcoeff.resize( 0 );
        _ccoeff.resize( nnz );
    }// if
    else
    {
        _rcoeff.resize( nnz );
        _ccoeff.resize( 0 );
    }// else

    // count number of entries per column
    for ( idx_t  col = 0; col < idx_t(ncols); ++col )
    {
        const idx_t  lb = acolptr[ col ];
        const idx_t  ub = acolptr[ col+1 ];

        for ( idx_t  j = lb; j < ub; j++ )
            _rowptr[ arowind[j] ]++;
    }// for

    // compute real rowptr array
    idx_t  pos = 0;

    for ( size_t  row = 0; row < nrows; ++row )
    {
        const idx_t  n = _rowptr[row];

        _rowptr[row] = pos;
        pos         += n;
    }// for
    _rowptr[nrows] = pos;

    if ( pos != idx_t(nnz) )
        HERROR( ERR_CONSISTENCY, "(TSparseMatrix) import_ccs",
                "given number of non-zeroes differs from counted number" );

    // set up rowind and coeff arrays
    vector< idx_t >   rowpos( nrows );

    for ( idx_t  col = 0; col < idx_t(ncols); ++col )
    {
        const idx_t  lb = acolptr[ col ];
        const idx_t  ub = acolptr[ col+1 ];

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  row = arowind[j];
            const idx_t  idx = _rowptr[row] + rowpos[row];

            _colind[ idx ] = col;

            if ( is_complex_val )
                _ccoeff[ idx ] = acoeffs[ j ];
            else
                _rcoeff[ idx ] = std::real( acoeffs[ j ] ); // for std::real(·) see above

            rowpos[row]++;
        }// for
    }// for
}

//
// export internal data to CCS format (real valued)
//
template <typename T_val>
T_val
sparse_rcoeff ( const TSparseMatrix *  S,
               const idx_t            i );

template <>
inline real
sparse_rcoeff< real > ( const TSparseMatrix *  S,
                       const idx_t            i )
{
    return S->rcoeff( i );
}

template <>
inline complex
sparse_rcoeff< complex > ( const TSparseMatrix *  S,
                          const idx_t            i )
{
    return S->ccoeff( i );
}

template <typename T_idx, typename T_val>
void
TSparseMatrix::export_ccs  ( vector< T_idx > &  colptr,
                             vector< T_idx > &  rowind,
                             vector< T_val > &  acoeffs,
                             const bool         use_sym ) const
{
    const bool  is_complex_val = is_complex_type< T_val >::value;

    if ( is_complex() && ! is_complex_val )
        HERROR( ERR_REAL_CMPLX, "(TSparseMatrix) export_ccs", "" );

    const size_t  nrows  = rows();
    const size_t  ncols  = cols();
    const bool    is_sym = use_sym && ( is_symmetric() || is_hermitian() );
    const size_t  nnz    = ( is_sym ? lower_left_nonzeroes( this ) : n_non_zero() );

    colptr.resize( ncols+1 );

    for ( idx_t  col = 0; col <= idx_t(ncols); ++col )
        colptr[ col ] = 0;

    // count number of entries per column
    for ( idx_t  row = 0; row < idx_t(nrows); ++row )
    {
        const idx_t  lb = rowptr(idx_t(row));
        const idx_t  ub = rowptr(idx_t(row)+1);

        if ( is_sym )
        {
            for ( idx_t  j = lb; j < ub; j++ )
                if ( colind(j) <= row ) // count only lower left part
                    colptr[ colind(j) ]++;
        }// if
        else
        {
            for ( idx_t  j = lb; j < ub; j++ )
                colptr[ colind(j) ]++;
        }// else
    }// for

    // compute real colptr array
    idx_t  pos = 0;

    for ( size_t  j = 0; j < ncols; j++ )
    {
        const idx_t  n = colptr[j];

        colptr[j] = T_idx(pos);
        pos      += n;
    }// for
    colptr[ncols] = T_idx(pos);

    // set up rowind and coeff arrays
    vector< idx_t >   colpos( ncols );

    rowind.resize( nnz );
    acoeffs.resize( nnz );

    for ( idx_t  row = 0; row < idx_t(nrows); row++ )
    {
        const idx_t  lb = rowptr(row);
        const idx_t  ub = rowptr(row+1);

        if ( is_sym )
        {
            for ( idx_t  j = lb; j < ub; j++ )
            {
                const idx_t  col = colind(j);
                const idx_t  idx = colptr[col] + colpos[col];

                if ( col > row ) // only lower left part
                    continue;

                rowind[ idx ]  = T_idx( row );
                acoeffs[ idx ] = sparse_rcoeff< T_val >( this, j );

                colpos[col]++;
            }// for
        }// if
        else
        {
            for ( idx_t  j = lb; j < ub; j++ )
            {
                const idx_t  col = colind(j);
                const idx_t  idx = colptr[col] + colpos[col];

                rowind[ idx ]  = T_idx( row );
                acoeffs[ idx ] = sparse_rcoeff< T_val >( this, j );

                colpos[col]++;
            }// for
        }// else
    }// for
}
        

//
// permute entries in sparse matrix
//
void
TSparseMatrix::permute ( const TPermutation &  rowperm,
                         const TPermutation &  colperm )
{
    TPermutation  inv_rowperm( rowperm );

    inv_rowperm.invert();

    vector< idx_t >    new_rowptr;
    vector< idx_t >    new_colind;
    vector< real >     new_rcoeff;
    vector< complex >  new_ccoeff;

    if ( is_complex() )
    {
        new_ccoeff.resize( _ccoeff.size() );
    }// if
    else
    {
        new_rcoeff.resize( _rcoeff.size() );
    }// else
    
    new_colind.resize( _colind.size() );
    new_rowptr.resize( _rowptr.size() );

    //
    // copy data from old to new position
    //
    
    idx_t  cpos = 0;
    idx_t  rpos = 0;
    
    for ( size_t i = 0; i < rows(); i++ )
    {
        const idx_t  old_row = inv_rowperm.permute( idx_t(i) );
        const idx_t  lb      = _rowptr[old_row];
        const idx_t  ub      = _rowptr[old_row+1];

        new_rowptr[ i ] = rpos;
        
        for ( idx_t  j = lb; j < ub; j++, rpos++, cpos++ )
        {
            const idx_t   new_col = colperm.permute( _colind[j] );

            if ( is_complex() ) new_ccoeff[cpos] = _ccoeff[j];
            else                new_rcoeff[cpos]  = _rcoeff[j];
            
            new_colind[cpos] = new_col;
        }// for
    }// for

    new_rowptr[ rows() ] = _rowptr[ rows() ];

    //
    // copy to local data
    //
    
    for ( size_t i = 0; i < _rcoeff.size();  i++ ) _rcoeff[i]  = new_rcoeff[i];
    for ( size_t i = 0; i < _ccoeff.size(); i++ ) _ccoeff[i] = new_ccoeff[i];
    for ( size_t i = 0; i < _colind.size(); i++ ) _colind[i] = new_colind[i];
    for ( size_t i = 0; i < _rowptr.size(); i++ ) _rowptr[i] = new_rowptr[i];
}

//
// check symmetry of matrix, return true if symmetric
//
bool
TSparseMatrix::test_symmetry () 
{
    //
    // go through sparsity pattern and compare entries
    // with transposed entries
    //

    const idx_t  n       = idx_t(rows());
    const idx_t  m       = idx_t(cols());
    bool         is_sym  = true;
    bool         is_herm = true;

    set_nonsym();
    
    for ( idx_t  i = 0; i < n; i++ )
    {
        // stop as soon as not symmetric and not hermitian
        if ( ! ( is_sym || is_herm ) )
            break;
        
        const idx_t  lb = _rowptr[i];
        const idx_t  ub = _rowptr[i+1];
        
        for ( idx_t j = lb; j < ub; j++ )
        {
            const idx_t  col = _colind[j];
            
            if (( col >= n ) || ( i >= m ))
                is_herm = is_sym = false;
            
            if ( is_complex() )
            {
                const complex  f = centry( col, i );

                if      ( _ccoeff[j] != f         ) is_sym  = false;
                else if ( _ccoeff[j] != conj( f ) ) is_herm = false;
            }// if
            else
            {
                if ( _rcoeff[j] != entry( col, i ) )
                    is_herm = is_sym = false;
            }// else
        }// for
    }// for

    if      ( is_sym  ) set_symmetric();
    else if ( is_herm ) set_hermitian();
    
    return is_sym || is_herm;
}

//
// compute average/maximal number of entries per row
//
size_t
TSparseMatrix::avg_entries_per_row () const
{
    return _colind.size() / rows();
}
    
size_t
TSparseMatrix::max_entries_per_row () const
{
    size_t max_entries = 0;

    for ( size_t i = 0; i < rows(); i++ )
        max_entries = std::max( size_t( _rowptr[i+1] - _rowptr[i] ), max_entries );

    return max_entries;
}

//
// return true if matrix contains zero or elements < ε on diagonal
//
bool
TSparseMatrix::has_diag_zero ( const real  eps )
{
    const size_t  nrows = rows();
    const size_t  ncols = cols();

    if ( is_complex() )
    {
        for ( idx_t  i = 0; i < std::min( idx_t( nrows ), idx_t( ncols ) ); ++i )
        {
            if ( ! has_entry( i, i ) || ( Math::abs( centry( i, i ) ) < eps ))
                return true;
        }// for
    }// if
    else
    {
        for ( idx_t  i = 0; i < std::min( idx_t( nrows ), idx_t( ncols ) ); ++i )
        {
            if ( ! has_entry( i, i ) || ( Math::abs( entry( i, i ) ) < eps ))
                return true;
        }// for
    }// else

    return false;
}

//
// return true if S is (weakly) diagonally dominant
//

namespace
{

//
// return ∑_j |S_ij|, j≠i
//
real
comp_off_diag ( const idx_t            i,
                const TSparseMatrix *  S )
{
    real         retval = real(0);
    const idx_t  lb     = S->rowptr( i );
    const idx_t  ub     = S->rowptr( i+1 );

    if ( S->is_complex() )
    {
        for ( idx_t  j = lb; j < ub; j++ )
        {
            if ( S->colind( j ) != i )
                retval += Math::abs( S->ccoeff( j ) );
        }// for
    }// if
    else
    {
        for ( idx_t  j = lb; j < ub; j++ )
        {
            if ( S->colind( j ) != i )
                retval += Math::abs( S->rcoeff( j ) );
        }// for
    }// else

    return retval;
}

}// namespace anonymous

bool
TSparseMatrix::is_diag_dom ( const bool  weak )
{
    if ( weak )
    {
        for ( auto  i : row_is() )
        {
            const real  off_diag = comp_off_diag( i, this );
            const real  diag     = Math::abs( is_complex() ? centry( i, i ) : entry( i, i ) );
            
            if ( diag < off_diag )
                return false;
        }// for
    }// if
    else
    {
        for ( auto  i : row_is() )
        {
            const real  off_diag = comp_off_diag( i, this );
            const real  diag     = Math::abs( is_complex() ? centry( i, i ) : entry( i, i ) );
            
            if ( diag <= off_diag )
                return false;
        }// for
    }// else

    return true;
}
    
//
// return factor α, such that S + α·I is diagonally dominant (this = S)
//
real
TSparseMatrix::diag_dom_factor ()
{
    real  alpha = real(0);
    
    for ( auto  i : row_is() )
    {
        const real  off_diag = comp_off_diag( i, this );
        const real  diag     = Math::abs( is_complex() ? centry( i, i ) : entry( i, i ) );

        if ( diag < off_diag )
            alpha = std::max( alpha, off_diag - diag );
    }// for

    return alpha;
}

//
// some tests on the matrix
//
void
TSparseMatrix::check_matrix () const
{
    //
    // print number of entries per row
    //

    size_t  min_entries = 0;
    size_t  max_entries = 0;
    size_t  min_pos     = 0;
    size_t  max_pos     = 0;
    
    for ( size_t i = 0; i < rows(); i++ )
    {
        size_t  nentries = _rowptr[i+1] - _rowptr[i];
        
        if (( min_entries == 0 ) || ( nentries < min_entries ))
        {
            min_pos     = i;
            min_entries = nentries;
        }// if
        
        if ( nentries > max_entries )
        {
            max_pos     = i;
            max_entries = nentries;
        }// if
    }// for

    std::cout << "(TSparseMatrix) check_matrix : "
              << "minimal/maximal/average degree per row = "
              << min_entries << "(at " << min_pos << ") / "
              << max_entries << "(at " << max_pos << ") / "
              << avg_entries_per_row()
              << std::endl;

    //
    // look for zero rows
    //

    size_t  n_zero_row = 0;
    
    for ( size_t row = 0; row < rows(); row++ )
    {
        const idx_t lb      = _colind[row];
        const idx_t ub      = _colind[row+1];
        real       row_sum = real(0);

        for ( idx_t j = lb; j < ub; j++ )
        {
            row_sum += _rcoeff[j];
        }// for

        if ( Math::abs( row_sum ) < 1e-8 )
        {
            n_zero_row++;
        }// if
    }// for

    std::cout << "(TSparseMatrix) check_matrix : "
              << n_zero_row << " rows with (nearly) zero row-sum"
              << std::endl;
}

/////////////////////////////////////////////////
//
// BLAS-routines (real valued)
//

//
// scale matrix by constant factor
//
void
TSparseMatrix::scale ( const real f )
{
    if ( f == real(1) )
        return;

    if ( is_complex() )
    {
        for ( size_t i = 0; i < _rcoeff.size(); i++ )
            _ccoeff[i] *= f;
    }// if
    else
    {
        for ( size_t i = 0; i < _rcoeff.size(); i++ )
            _rcoeff[i] *= f;
    }// else
}
    
//
// matrix-vector-multiplication
//
void
TSparseMatrix::mul_vec ( const real      alpha,
                         const TVector * x,
                         const real      beta,
                         TVector *       y,
                         const matop_t   op ) const
{
    //
    // check for nullptr, index set and type
    //
    
    if ( x == nullptr ) HERROR( ERR_ARG, "(TSparseMatrix) mul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TSparseMatrix) mul_vec", "y = nullptr" );

    if (( is_complex() != x->is_complex() ) || ( is_complex() != y->is_complex() ))
        HERROR( ERR_REAL_CMPLX, "(TSparseMatrix) mul_vec", "" );
    
    if ( op == MATOP_NORM )
    {
        if (( row_is() != y->is() ) || ( col_is() != x->is() ))
            HERROR( ERR_VEC_STRUCT, "(TSparseMatrix) mul_vec", "incompatible vector index set" );
    }// if
    else
    {
        if (( col_is() != y->is() ) || ( row_is() != x->is() ))
            HERROR( ERR_VEC_STRUCT, "(TSparseMatrix) mul_vec", "incompatible vector index set" );
    }// if
    
    if ( ! IS_TYPE( x, TScalarVector ) || ! IS_TYPE( y, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TSparseMatrix) mul_vec",
                y->typestr() + " += TSparseMatrix * " + x->typestr() );

    //
    // apply β
    //

    if      ( beta == real(0) ) y->fill( real(0) );
    else if ( beta != real(1) ) y->scale( beta );

    //
    // stop, if nothing to do
    //
    
    if (( alpha == real(0) ) || ( n_non_zero() == 0 ))
        return;

    //
    // entry-wise multiplication
    //

    const TScalarVector *  sx = cptrcast( x, TScalarVector );
    TScalarVector *        sy = ptrcast( y, TScalarVector );
    const size_t           n  = rows();
            
    if ( is_complex() )
    {
        if ( op == MATOP_NORM )
        {
            //
            // complex matrix
            //
                
            for ( idx_t i = 0; i < idx_t(n); i++ )
            {
                const idx_t  lb = _rowptr[i];
                const idx_t  ub = _rowptr[i+1];
                complex      s  = real(0);
                    
                for ( idx_t j = lb; j < ub; j++ )
                    s += _ccoeff[j] * sx->blas_cvec()( _colind[j] );

                sy->blas_cvec()( i ) += alpha * s;
            }// for
        }// if
        else if ( op == MATOP_TRANS )
        {
            //
            // complex, transposed matrix
            //
                
            for ( idx_t  i = 0; i < idx_t(n); i++ )
            {
                const idx_t    lb = _rowptr[i];
                const idx_t    ub = _rowptr[i+1];
                const complex  a  = alpha * sx->blas_cvec()( i );
                
                for ( idx_t j = lb; j < ub; j++ )
                    sy->blas_cvec()( _colind[j] ) += _ccoeff[j] * a;
            }// for
        }// if
        else if ( op == MATOP_ADJ )
        {
            //
            // complex, adjoint matrix
            //
                
            for ( idx_t  i = 0; i < idx_t(n); i++ )
            {
                const idx_t    lb = _rowptr[i];
                const idx_t    ub = _rowptr[i+1];
                const complex  a  = alpha * sx->blas_cvec()( i );

                for ( idx_t j = lb; j < ub; j++ )
                    sy->blas_cvec()( _colind[j] ) += conj( _ccoeff[j] ) * a;
            }// for
        }// if
    }// if
    else
    {
        if ( op == MATOP_NORM )
        {
            //
            // real matrix
            //

            for ( idx_t  i = 0; i < idx_t(n); i++ )
            {
                const idx_t  lb = _rowptr[i];
                const idx_t  ub = _rowptr[i+1];
                real         s  = 0;
                    
                for ( idx_t j = lb; j < ub; j++ )
                    s += _rcoeff[j] * sx->blas_rvec()( _colind[j] );
                    
                sy->blas_rvec()( i ) += alpha * s;
            }// for
        }// if
        else
        {
            //
            // real, transposed matrix
            //
                
            for ( idx_t  i = 0; i < idx_t(n); i++ )
            {
                const idx_t  lb = _rowptr[i];
                const idx_t  ub = _rowptr[i+1];
                const real   a  = alpha * sx->blas_rvec()( i );

                for ( idx_t  j = lb; j < ub; j++ )
                    sy->blas_rvec()( _colind[j] ) += _rcoeff[j] * a;
            }// for
        }// else
    }// else
}

//
// compute this = this + a * M
//
void
TSparseMatrix::add ( const real a, const TMatrix * M )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "(TSparseMatrix) add", "matrix is nullptr" );

    if ( ! IS_TYPE( M, TSparseMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TSparseMatrix) add", M->typestr() );

    if ( a == real(0) )
        return;
    
    const TSparseMatrix * S = cptrcast( M, TSparseMatrix );

    if (( _rcoeff.size()  != S->_rcoeff.size() ) ||
        ( _colind.size() != S->_colind.size() ))
        HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) add", "argument has different CRS format" );
    
    //
    // check for same structure
    //
    
    for ( size_t i = 0; i < _rowptr.size(); i++ )
        if ( _rowptr[i] != S->_rowptr[i] )
            HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) add", "argument has different CRS format" );
    
    for ( size_t i = 0; i < _colind.size(); i++ )
        if ( _colind[i] != S->_colind[i] )
            HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) add", "argument has different CRS format" );
    
    if ( S->is_complex() != is_complex() )
        HERROR( ERR_NOT_IMPL, "(TSparseMatrix) add", "different field type (real <=> complex)" );
    
    // finally: same format, just add entries
    // (assuming same ordering !!!)
    if ( is_complex() )
    {
        for ( size_t  i = 0; i < _ccoeff.size(); i++ )
            _ccoeff[i] += a * S->_ccoeff[i];
    }// if
    else
    {
        for ( size_t  i = 0; i < _rcoeff.size(); i++ )
            _rcoeff[i] += a * S->_rcoeff[i];
    }// else
}

/////////////////////////////////////////////////
//
// BLAS-routines (complex valued)
//

//
// scale matrix by constant factor
//
void
TSparseMatrix::cscale ( const complex f )
{
    if ( f == real(1) )
        return;

    if ( is_complex() || ( std::imag( f ) != real(0) ))
    {
        set_complex( true );
        
        for ( size_t i = 0; i < _ccoeff.size(); i++ )
            _ccoeff[i] *= f;
    }// if
    else
    {
        const real  g = std::real( f );
        
        for ( size_t i = 0; i < _rcoeff.size(); i++ )
            _rcoeff[i] *= g;
    }// else
}
    
//
// matrix-vector-multiplication
//
void
TSparseMatrix::cmul_vec ( const complex   alpha,
                          const TVector * x,
                          const complex   beta,
                          TVector       * y,
                          const matop_t   op ) const
{
    //
    // check for nullptr, index set and type
    //
    
    if ( x == nullptr ) HERROR( ERR_ARG, "(TSparseMatrix) cmul_vec", "x = nullptr" );
    if ( y == nullptr ) HERROR( ERR_ARG, "(TSparseMatrix) cmul_vec", "y = nullptr" );
    
    if ( ! is_complex() || ( is_complex() != x->is_complex() ) || ( is_complex() != y->is_complex() ))
        HERROR( ERR_REAL_CMPLX, "(TSparseMatrix) cmul_vec", "" );

    if ( op == MATOP_NORM )
    {
        if (( row_is() != y->is() ) || ( col_is() != x->is() ))
            HERROR( ERR_VEC_STRUCT, "(TSparseMatrix) cmul_vec", "incompatible vector index set" );
    }// if
    else
    {
        if (( col_is() != y->is() ) || ( row_is() != x->is() ))
            HERROR( ERR_VEC_STRUCT, "(TSparseMatrix) cmul_vec", "incompatible vector index set" );
    }// if
    
    if ( ! IS_TYPE( x, TScalarVector ) || ! IS_TYPE( y, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TSparseMatrix) cmul_vec",
                y->typestr() + " += TSparseMatrix * " + x->typestr() );
    
    //
    // apply β
    //

    if      ( beta == real(0) ) y->fill( real(0) );
    else if ( beta != real(1) ) y->cscale( beta );

    //
    // stop, if nothing to do
    //
    
    if (( alpha == complex(0) ) || ( n_non_zero() == 0 ))
        return;

    //
    // entry-wise multiplication
    //
    
    const TScalarVector  * sx = cptrcast( x, TScalarVector );
    TScalarVector        * sy = ptrcast( y, TScalarVector );
    const size_t           n  = rows();

    if ( op == MATOP_NORM )
    {
        //
        // complex matrix, complex argument
        //
                
        for ( idx_t  i = 0; i < idx_t(n); i++ )
        {
            const idx_t  lb = _rowptr[i];
            const idx_t  ub = _rowptr[i+1];
            complex      s  = real(0);

            for ( idx_t j = lb; j < ub; j++ )
                s += _ccoeff[j] * sx->blas_cvec()( _colind[j] );

            sy->blas_cvec()( i ) += alpha * s;
        }// for
    }// if
    else if ( op == MATOP_TRANS )
    {
        //
        // complex, transposed matrix, complex argument
        //
                
        for ( idx_t  i = 0; i < idx_t(n); i++ )
        {
            const idx_t    lb = _rowptr[i];
            const idx_t    ub = _rowptr[i+1];
            const complex  a  = alpha * sx->blas_cvec()( i );

            for ( idx_t j = lb; j < ub; j++ )
                sy->blas_cvec()( _colind[j] ) += _ccoeff[j] * a;
        }// for
    }// if
    else if ( op == MATOP_ADJ )
    {
        //
        // complex, adjoint matrix, complex argument
        //
                
                
        for ( idx_t  i = 0; i < idx_t(n); i++ )
        {
            const idx_t    lb = _rowptr[i];
            const idx_t    ub = _rowptr[i+1];
            const complex  a  = alpha * sx->blas_cvec()( i );

            for ( idx_t j = lb; j < ub; j++ )
                sy->blas_cvec()( _colind[j] ) += conj( _ccoeff[j] ) * a;
        }// for
    }// if
}

//
// compute this = this + a * M
//
void
TSparseMatrix::cadd ( const complex a, const TMatrix * M )
{
    if ( M == nullptr )
        HERROR( ERR_ARG, "(TSparseMatrix) cadd", "argument is nullptr" );

    if ( ! IS_TYPE( M, TSparseMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TSparseMatrix) cadd", M->typestr() );

    if ( a == real(0) )
        return;
    
    if ( ! is_complex() || ! M->is_complex() )
        HERROR( ERR_REAL_CMPLX, "(TSparseMatrix) cadd", "" );
    
    const TSparseMatrix * S = cptrcast( M, TSparseMatrix );

    if (( _ccoeff.size() != S->_ccoeff.size() ) ||
        ( _colind.size() != S->_colind.size() ))
        HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) cadd", "argument has different CRS format" );
        
    //
    // check for same structure
    //

    for ( size_t i = 0; i < _rowptr.size(); i++ )
        if ( _rowptr[i] != S->_rowptr[i] )
            HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) cadd", "argument has different CRS format" );
        
    for ( size_t i = 0; i < _colind.size(); i++ )
        if ( _colind[i] != S->_colind[i] )
            HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) cadd", "argument has different CRS format" );

    // finally: same format, just add entries
    for ( size_t i = 0; i < _ccoeff.size(); i++ )
        _ccoeff[i] += a * S->_ccoeff[i];
}

/////////////////////////////////////////////////
//
// input/output related
//

//
// output of matrix
//
void
TSparseMatrix::print ( const uint ofs ) const
{
    TMatrix::print( ofs );
}

//
// serialisation
//

void
TSparseMatrix::read  ( TByteStream & s )
{
    return build( s );
}

void
TSparseMatrix::build ( TByteStream & s )
{
    TMatrix::read( s );

    s.get( & _rows, sizeof(size_t) );
    s.get( & _cols, sizeof(size_t) );

    //
    // read CRS data
    //
    
    size_t  nnz;
    
    s.get( & nnz, sizeof(size_t) );
    
    init( nnz );
    
    s.get( _rowptr.data(), _rows+1 * sizeof(idx_t) );
    s.get( _colind.data(), nnz * sizeof(idx_t) );
    
    if ( is_complex() )
    {
        s.get( _ccoeff.data(), nnz * sizeof(complex) );
    }// if
    else
    {
        s.get( _rcoeff.data(), nnz * sizeof(real) );
    }// else
}

void
TSparseMatrix::write ( TByteStream & s ) const
{
    TMatrix::write( s );

    s.put( & _rows, sizeof(size_t) );
    s.put( & _cols, sizeof(size_t) );

    const size_t  nnz  = _colind.size();
        
    s.put( & nnz, sizeof(nnz) );
    s.put( const_cast< vector< idx_t > & >(_rowptr).data(), sizeof(size_t) * _rowptr.size() );
    s.put( const_cast< vector< idx_t > & >(_colind).data(), sizeof(size_t) * _colind.size() );
    
    if ( is_complex() )
    {
        s.put( const_cast< vector< complex > & >(_ccoeff).data(), sizeof(complex) * _ccoeff.size() );
    }// if
    else
    {
        s.put( const_cast< vector< real > & >(_rcoeff).data(),  sizeof(real) * _rcoeff.size() );
    }// else
}

//
// returns size of object in bytestream
//
size_t
TSparseMatrix::bs_size () const
{
    size_t  size = TMatrix::bs_size() + 2*sizeof(size_t);

    size += sizeof(size_t);
    size += _rowptr.size() * sizeof(size_t);
    size += _colind.size() * sizeof(size_t);

    if ( is_complex() ) size += _ccoeff.size() * sizeof(complex);
    else                size += _rcoeff.size()  * sizeof(real);

    return size;
}

/////////////////////////////////////////////////
//
// misc.
//

//
// transpose matrix
//
void
TSparseMatrix::transpose ()
{
    // to be done
    HERROR( ERR_NOT_IMPL, "(TSparseMatrix) transpose", "" );

    TMatrix::transpose();
}
    
//
// conjugate matrix coefficients
//
void
TSparseMatrix::conjugate ()
{
    if ( is_complex() && ! is_hermitian() )
    {
        for ( size_t  i = 0; i < _ccoeff.size(); ++i )
            _ccoeff[i] = conj( _ccoeff[i] );
    }// if
}

//
// restrict sparse matrix to given block index set
//
unique_ptr< TSparseMatrix >
TSparseMatrix::restrict ( const TIndexSet &  rowis,
                          const TIndexSet &  colis ) const
{
    if ( ! block_is().is_sub( TBlockIndexSet( rowis, colis ) ) )
        HERROR( ERR_MAT_STRUCT, "TSparseMatrix) restrict", "given block index set is not a local sub set" );
                
    auto  T = make_unique< TSparseMatrix >( rowis.size(), colis.size() );

    T->set_block_is( TBlockIndexSet( rowis, colis ) );
    T->set_complex( is_complex() );
    
    //
    // determine total number of nonzeroes and number per row in subblock
    //

    size_t           sub_nnzero = 0;
    vector< idx_t >  sub_rowptr( rowis.size()+1 );

    for ( auto  row : rowis )
    {
        const idx_t  lb = _rowptr[ row - row_ofs() ];
        const idx_t  ub = _rowptr[ row - row_ofs() + 1 ];

        sub_rowptr[ row - rowis.first() ] = idx_t(sub_nnzero);
        
        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col = _colind[j] + col_ofs();
            
            if ( colis.is_in( col ) )
                sub_nnzero++;
        }// for
    }// for
    sub_rowptr[ rowis.size() ] = idx_t(sub_nnzero);

    T->init( sub_nnzero );

    //
    // set up CRS data in sub matrix
    //

    for ( idx_t  i = 0; i < idx_t(sub_rowptr.size()); ++i )
        T->rowptr( i ) = sub_rowptr[i];

    idx_t  pos = 0;
    
    for ( auto  row : rowis )
    {
        const idx_t  lb = _rowptr[ row - row_ofs() ];
        const idx_t  ub = _rowptr[ row - row_ofs() + 1 ];

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col = _colind[j] + col_ofs();
            
            if ( colis.is_in( col ) )
            {
                T->colind( pos ) = col - colis.first();

                if ( is_complex() )
                    T->ccoeff( pos ) = _ccoeff[j];
                else
                    T->rcoeff( pos ) = _rcoeff[j];
                pos++;
            }// if
        }// for
    }// for

    return T;
}

//
// restrict sparse matrix to given block index set
// which is given in a different ordering
//
unique_ptr< TSparseMatrix >
TSparseMatrix::restrict ( const TIndexSet &     rowis,
                          const TPermutation *  rowperm,
                          const TIndexSet &     colis,
                          const TPermutation *  colperm ) const
{
    if ( ! block_is().is_sub( TBlockIndexSet( rowis, colis ) ) )
        HERROR( ERR_MAT_STRUCT, "TSparseMatrix) restrict", "given block index set is not a local sub set" );
                
    auto  T = make_unique< TSparseMatrix >( rowis.size(), colis.size() );

    T->set_block_is( TBlockIndexSet( rowis, colis ) );
    T->set_complex( is_complex() );
    
    //
    // determine total number of nonzeroes and number per row in subblock
    //

    size_t           sub_nnzero = 0;
    vector< idx_t >  sub_rowptr( rowis.size()+1 );

    for ( auto  row : rowis )
    {
        const idx_t  prow = ( rowperm != nullptr ? rowperm->permute( row ) : row );
        const idx_t  lb   = _rowptr[ prow - row_ofs() ];
        const idx_t  ub   = _rowptr[ prow - row_ofs() + 1 ];

        sub_rowptr[ row - rowis.first() ] = idx_t(sub_nnzero);
        
        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col  = _colind[j] + col_ofs();
            const idx_t  pcol = ( colperm != nullptr ? colperm->permute( col ) : col );
            
            if ( colis.is_in( pcol ) )
                sub_nnzero++;
        }// for
    }// for
    sub_rowptr[ rowis.size() ] = idx_t(sub_nnzero);

    T->init( sub_nnzero );

    //
    // set up CRS data in sub matrix
    //

    for ( idx_t  i = 0; i < idx_t(sub_rowptr.size()); ++i )
        T->rowptr( i ) = sub_rowptr[i];

    idx_t  pos = 0;
    
    for ( auto  row : rowis )
    {
        const idx_t  prow = ( rowperm != nullptr ? rowperm->permute( row ) : row );
        const idx_t  lb   = _rowptr[ prow - row_ofs() ];
        const idx_t  ub   = _rowptr[ prow - row_ofs() + 1 ];

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col  = _colind[j] + col_ofs();
            const idx_t  pcol = ( colperm != nullptr ? colperm->permute( col ) : col );
            
            if ( colis.is_in( pcol ) )
            {
                T->colind( pos ) = pcol - colis.first();

                if ( is_complex() ) T->ccoeff( pos ) = _ccoeff[j];
                else                T->rcoeff(  pos ) = _rcoeff[j];
                pos++;
            }// if
        }// for
    }// for

    return T;
}

//
// return number of coefficients in sub block index set \a rowis × \a colis
//
size_t
TSparseMatrix::restrict_nonzeroes ( const TIndexSet &  rowis,
                                    const TIndexSet &  colis ) const
{
    if ( ! block_is().is_sub( TBlockIndexSet( rowis, colis ) ) )
        HERROR( ERR_MAT_STRUCT, "(TSparseMatrix) restrict_nonzeroes", "given block index set is not a local sub set" );
                
    //
    // determine total number of nonzeroes and number per row in subblock
    //

    size_t  sub_nnzero = 0;

    for ( auto  row : rowis )
    {
        const idx_t  lb = _rowptr[ row - row_ofs() ];
        const idx_t  ub = _rowptr[ row - row_ofs() + 1 ];

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col = _colind[j] + col_ofs();
            
            if ( colis.is_in( col ) )
                sub_nnzero++;
        }// for
    }// for

    return sub_nnzero;
}

size_t
TSparseMatrix::restrict_nonzeroes ( const TIndexSet &     rowis,
                                    const TPermutation *  rowperm,
                                    const TIndexSet &     colis,
                                    const TPermutation *  colperm ) const
{
    if ( ! block_is().is_sub( TBlockIndexSet( rowis, colis ) ) )
        HERROR( ERR_MAT_STRUCT, "TSparseMatrix) restrict", "given block index set is not a local sub set" );
                
    //
    // determine total number of nonzeroes and number per row in subblock
    //

    size_t  sub_nnzero = 0;

    for ( auto  row : rowis )
    {
        const idx_t  prow = ( rowperm != nullptr ? rowperm->permute( row ) : row );
        const idx_t  lb   = _rowptr[ prow - row_ofs() ];
        const idx_t  ub   = _rowptr[ prow - row_ofs() + 1 ];

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col  = _colind[j] + col_ofs();
            const idx_t  pcol = ( colperm != nullptr ? colperm->permute( col ) : col );
            
            if ( colis.is_in( pcol ) )
                sub_nnzero++;
        }// for
    }// for

    return sub_nnzero;
}

//
// virtual constructors
//
std::unique_ptr< TMatrix >
TSparseMatrix::copy () const
{
    auto             M = TMatrix::copy();
    TSparseMatrix *  S = ptrcast( M.get(), TSparseMatrix );

    if ( cluster() != nullptr )
        S->set_cluster( cluster() );
    else
        S->set_block_is( block_is() );

    const size_t  nnz = _colind.size();

    S->init( nnz );

    for ( size_t i = 0; i <= rows(); i++ )
        S->_rowptr[i] = _rowptr[i];

    if ( is_complex() )
    {
        for ( size_t i = 0; i < nnz; i++ )
        {
            S->_ccoeff[i]  = _ccoeff[i];
            S->_colind[i] = _colind[i];
        }// for
    }// if
    else
    {
        for ( size_t i = 0; i < nnz; i++ )
        {
            S->_rcoeff[i]   = _rcoeff[i];
            S->_colind[i] = _colind[i];
        }// for
    }// else

    return M;
}

//
// copy matrix into A
//
void
TSparseMatrix::copy_to ( TMatrix * A ) const
{
    TMatrix::copy_to( A );
    
    if ( ! IS_TYPE( A, TSparseMatrix ) )
        HERROR( ERR_MAT_TYPE, "(TSparseMatrix) copy_to", A->typestr() );
    
    TSparseMatrix * M = ptrcast( A, TSparseMatrix );

    M->set_complex( is_complex() );
    if ( cluster() != nullptr )
        M->set_cluster( cluster() );
    else
        M->set_block_is( block_is() );

    const size_t  nnz = _colind.size();
        
    M->init( nnz );

    for ( size_t i = 0; i <= rows(); i++ )
        M->_rowptr[i] = _rowptr[i];
    
    if ( is_complex() )
    {
        for ( size_t i = 0; i < nnz; i++ )
        {
            M->_colind[i] = _colind[i];
            M->_ccoeff[i] = _ccoeff[i];
        }// for
    }// if
    else
    {
        for ( size_t i = 0; i < nnz; i++ )
        {
            M->_colind[i] = _colind[i];
            M->_rcoeff[i] = _rcoeff[i];
        }// for
    }// else
}

//
// return structural copy
//
std::unique_ptr< TMatrix >
TSparseMatrix::copy_struct () const
{
    auto             M = TMatrix::copy_struct();
    TSparseMatrix *  S = ptrcast( M.get(), TSparseMatrix );

    if ( cluster() != nullptr )
        S->set_cluster( cluster() );
    else
        S->set_block_is( block_is() );

    S->init( 0 );

    return M;
}

//
// return size in bytes used by this object
//
size_t
TSparseMatrix::byte_size () const
{
    size_t  size = TMatrix::byte_size() + 2*sizeof(size_t);

    size += _rcoeff.size()  * sizeof(real)    + sizeof(_rcoeff);
    size += _ccoeff.size() * sizeof(complex) + sizeof(_ccoeff);
    size += _colind.size() * sizeof(idx_t)   + sizeof(_colind);
    size += _rowptr.size() * sizeof(idx_t)   + sizeof(_rowptr);
    
    return size;
}

//
// test data for invalid values, e.g. INF and NAN
//
void
TSparseMatrix::check_data () const
{
    if ( is_complex() )
    {
        for ( vector< complex >::const_iterator  iter = _ccoeff.begin();
              iter != _ccoeff.end();
              ++iter )
        {
            const complex  val = *iter;
            
            if ( Math::is_inf( val ) )
                HERROR( ERR_INF, "(TSparseMatrix) check_data", "" );
            
            if ( Math::is_nan( val ) )
                HERROR( ERR_NAN, "(TSparseMatrix) check_data", "" );
        }// for
    }// if
    else
    {
        for ( vector< real >::const_iterator  iter = _rcoeff.begin();
              iter != _rcoeff.end();
              ++iter )
        {
            const real  val = *iter;
            
            if ( Math::is_inf( val ) )
                HERROR( ERR_INF, "(TSparseMatrix) check_data", "" );
            
            if ( Math::is_nan( val ) )
                HERROR( ERR_NAN, "(TSparseMatrix) check_data", "" );
        }// for
    }// else
}

//
// print histogram for entries per row in GNUPLOT format
//
void
TSparseMatrix::print_pattern_hist ( std::ostream & os ) const
{
    //
    // compute histogram
    //
    
    const size_t      max_deg = max_entries_per_row();
    vector< size_t >  hist( max_deg + 1 );

    for ( idx_t  i = 0; i < idx_t( rows() ); ++i )
    {
        hist[ rowptr(i+1) - rowptr(i) ]++;
    }// for

    //
    // print to stream
    //
    
    os << "# See GNUPLOT manual for a description of the format" << std::endl
       << std::endl
       << "# <number of coeffs per row>  <number of rows in matrix>" << std::endl;

    for ( size_t  i = 0; i <= max_deg; ++i )
        os << i << " " << hist[i] << std::endl;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

#define INST_IMPORT_CRS( T_idx, T_val )                                  \
    template void                                                        \
    TSparseMatrix::import_crs< T_idx, T_val > ( const size_t   nrows,    \
                                                const size_t   ncols,    \
                                                const size_t   nnonzero, \
                                                const T_idx *  rowptr,   \
                                                const T_idx *  colind,   \
                                                const T_val *  coeffs )

INST_IMPORT_CRS( int,   real );
INST_IMPORT_CRS( idx_t, real );
INST_IMPORT_CRS( int,   complex );
INST_IMPORT_CRS( idx_t, complex );

#define INST_IMPORT_CCS( T_idx, T_val )                                  \
    template void                                                        \
    TSparseMatrix::import_ccs< T_idx, T_val > ( const size_t   nrows,    \
                                                const size_t   ncols,    \
                                                const size_t   nnonzero, \
                                                const T_idx *  colptr,   \
                                                const T_idx *  rowind,   \
                                                const T_val *  coeffs )

INST_IMPORT_CCS( int,   real );
INST_IMPORT_CCS( idx_t, real );
INST_IMPORT_CCS( int,   complex );
INST_IMPORT_CCS( idx_t, complex );

#define INST_EXPORT_CCS( T_idx, T_val )                                  \
    template void                                                        \
    TSparseMatrix::export_ccs< T_idx, T_val > ( std::vector< T_idx > &  colptr, \
                                                std::vector< T_idx > &  rowind, \
                                                std::vector< T_val > &  coeffs, \
                                                const bool              use_sym ) const

INST_EXPORT_CCS( idx_t, real );
INST_EXPORT_CCS( idx_t, complex );

}// namespace HLIB
