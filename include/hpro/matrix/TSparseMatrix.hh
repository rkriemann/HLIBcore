#ifndef __HLIB_TSPARSEMATRIX_HH
#define __HLIB_TSPARSEMATRIX_HH
//
// Project     : HLib
// File        : TSparseMatrix.hh
// Description : class for a sparse matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/cluster/TPermutation.hh"

#include "hpro/matrix/TMatrix.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TSparseMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    TSparseMatrix
//! \brief    Class for a sparse matrix stored in compressed row storage format.
//
class TSparseMatrix : public TMatrix
{
private:
    //! number of rows in matrix
    size_t                  _rows;

    //! number of columns in matrix
    size_t                  _cols;

    //
    // data for compressed row storage
    //

    //! start index in colind and coefficient arrays for each row: _rowptr[i] for row i
    std::vector< idx_t >    _rowptr;

    //! column indices for all coefficients; indices are local, i.e. ∈ [0,cols()-1]
    std::vector< idx_t >    _colind;

    //! matrix coefficients (real valued)
    std::vector< real >     _rcoeff;

    //! matrix coefficients (complex valued)
    std::vector< complex >  _ccoeff;

public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct sparse matrix of size \a anrows × \a ancols
    TSparseMatrix ( const size_t  anrows,
                    const size_t  ancols )
            : _rows(0)
            , _cols(0)
    {
        set_size( anrows, ancols );
    }
    
    //! construct sparse matrix of size \a anrows × \a ancols
    TSparseMatrix ( const TIndexSet  arow_is,
                    const TIndexSet  acol_is )
            : _rows( arow_is.first() )
            , _cols( acol_is.first() )
    {
        set_size( arow_is.size(), acol_is.size() );
    }
    
    //! construct sparse matrix with size defined by block cluster \a bct
    TSparseMatrix ( const TBlockCluster * bct = nullptr )
            : TMatrix( bct )
            , _rows(0)
            , _cols(0)
    {
        if ( bct != nullptr )
            set_cluster( bct );
    }

    virtual ~TSparseMatrix ();

    /////////////////////////////////////////////////
    //
    // access data
    //

    //! set block cluster of matrix
    virtual void  set_cluster  ( const TBlockCluster * bct );

    //! directly set dimension of matrix
    virtual void  set_size     ( const size_t  nrows,
                                 const size_t  ncols );
    
    //! initialise CRS data for \a nnz non-zero entries
    void          init         ( const size_t  nnz );
    
    //! return number of rows in matrix
    size_t        rows         () const { return _rows; }

    //! return number of columns in matrix
    size_t        cols         () const { return _cols; }

    //! return number of non-zero elements
    size_t        n_non_zero   () const { return _colind.size(); }

    //! convert coefficients to real valued representation (if possible)
    virtual void  to_real      ();

    //! convert coefficients to complex valued representation
    virtual void  to_complex   ();

    /////////////////////////////////////////////////
    //
    // manage stored entries
    //

    //! return matrix coefficient a_ij (real valued)
    virtual real           entry      ( const idx_t    i,
                                        const idx_t    j ) const;

    //! return matrix coefficient a_ij (complex valued)
    virtual const complex  centry     ( const idx_t    i,
                                        const idx_t    j ) const;

    //! set matrix coefficient a_ij if existent (real valued)
    virtual void           set_entry  ( const idx_t    i,
                                        const idx_t    j,
                                        const real     c );

    //! set matrix coefficient a_ij if existent (complex valued)
    virtual void           set_entry  ( const idx_t    i,
                                        const idx_t    j,
                                        const complex  c );
    
    //! add \a c to matrix coefficient a_ij if existent (real valued)
    virtual void           add_entry  ( const idx_t    i,
                                        const idx_t    j,
                                        const real     c );

    //! add \a c to matrix coefficient a_ij if existent (complex valued)
    virtual void           add_entry  ( const idx_t    i,
                                        const idx_t    j,
                                        const complex  c );

    //! return true if entry (\a i, \a j) exists
    virtual bool           has_entry  ( const idx_t    i,
                                        const idx_t    j ) const;


    //! sort matrix coefficients per row with respect to column
    void sort_entries ();

    
    //! return \a i'th row pointer (constant)
    idx_t          rowptr ( const idx_t  i ) const { return _rowptr[i]; }
    //! return \a i'th row pointer 
    idx_t &        rowptr ( const idx_t  i ) { return _rowptr[i]; }
    
    //! return \a i'th column index (constant)
    idx_t          colind ( const idx_t  i ) const { return _colind[i]; }
    //! return \a i'th column index
    idx_t &        colind ( const idx_t  i ) { return _colind[i]; }
    
    //! return \a i'th coefficient (real valued, constant)
    real           rcoeff  ( const idx_t  i ) const { return _rcoeff[i]; }
    //! return \a i'th coefficient (real valued)
    real &         rcoeff  ( const idx_t  i ) { return _rcoeff[i]; }
    
    //! return \a i'th coefficient (complex valued, constant)
    const complex  ccoeff ( const idx_t  i ) const { return _ccoeff[i]; }
    //! return \a i'th coefficient (complex valued)
    complex &      ccoeff ( const idx_t  i ) { return _ccoeff[i]; }

    
    //! permute entries in sparse matrix according to \a rowperm (for rows)
    //! and \a colperm (for columns)
    void permute ( const TPermutation &  rowperm,
                   const TPermutation &  colperm );

    //! return true if matrix is symmetric (really check data)
    bool test_symmetry ();


    //! compute average number of entries per row
    size_t avg_entries_per_row () const;
    
    //! compute maximal number of entries per row
    size_t max_entries_per_row () const;

    //! return true if matrix contains zero or elements < ε on diagonal
    bool   has_diag_zero       ( const real  eps = 0.0 );

    //! return true if S is (weakly) diagonally dominant
    bool   is_diag_dom         ( const bool  weak = false );
    
    //! return factor α, such that S + α·I is diagonally dominant (this = S)
    real   diag_dom_factor     ();
        
    //! perform some tests on the matrix (for debugging)
    void   check_matrix        () const;
    
    /////////////////////////////////////////////////
    //
    // data conversion
    //

    //! import CRS data into local matrix
    template <typename T_idx, typename T_val>
    void import_crs  ( const size_t   nrows,
                       const size_t   ncols,
                       const size_t   nnonzero,
                       const T_idx *  rowptr,
                       const T_idx *  colind,
                       const T_val *  coeffs );

    //! import CCS data into local matrix
    template <typename T_idx, typename T_val>
    void import_ccs  ( const size_t   nrows,
                       const size_t   ncols,
                       const size_t   nnonzero,
                       const T_idx *  colptr,
                       const T_idx *  rowind,
                       const T_val *  coeffs );

    //! export internal data to CCS format (real valued);
    //! if \a use_sym is true, only the lower triangular part is exported
    template <typename T_idx, typename T_val>
    void export_ccs  ( std::vector< T_idx > &  colptr,
                       std::vector< T_idx > &  rowind,
                       std::vector< T_val > &  coeffs,
                       const bool              use_sym ) const;

    /////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    //! compute this ≔ α·this
    virtual void scale ( const real alpha );
    
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
    
    /////////////////////////////////////////////////
    //
    // BLAS-routines (complex valued)
    //

    //! compute this ≔ α·this
    virtual void cscale ( const complex  alpha );
    
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
    // truncation and restriction
    //

    //! truncate matrix to given accuracy (NOT YET IMPLEMENTED)
    virtual void truncate ( const TTruncAcc & )
    {
        // TODO: apply sparse matrix truncation
    }

    //! restrict sparse matrix to block index set \a rowis × \a colis
    std::unique_ptr< TSparseMatrix >  restrict ( const TIndexSet &  rowis,
                                                 const TIndexSet &  colis ) const;
    
    //! restrict sparse matrix to block index set \a rowis × \a colis
    //! which is given in a different ordering
    //! - \a rowperm permutes the row indices to the local ordering of the sparse matrix,
    //!   while \a colperm permutes the local column indices to the ordering of colis
    //! - the returned matrix has the same ordering as \a rowis and \a colis
    //! 
    std::unique_ptr< TSparseMatrix >  restrict ( const TIndexSet &     rowis,
                                                 const TPermutation *  rowperm,
                                                 const TIndexSet &     colis,
                                                 const TPermutation *  colperm ) const;
        
    //! return number of coefficients in sub block index set \a rowis × \a colis
    size_t           restrict_nonzeroes ( const TIndexSet &  rowis,
                                          const TIndexSet &  colis ) const;
    
    //! return number of coefficients in sub block index set \a rowis × \a colis
    size_t           restrict_nonzeroes ( const TIndexSet &     rowis,
                                          const TPermutation *  rowperm,
                                          const TIndexSet &     colis,
                                          const TPermutation *  colperm ) const;
    
    /////////////////////////////////////////////////
    //
    // input/output related
    //

    //! print matrix to stdout
    virtual void print ( const uint ofs = 0 ) const;
    
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
    virtual auto  create       () const -> std::unique_ptr< TMatrix >
    {
        return std::make_unique< TSparseMatrix >();
    }
    
    //! return copy of matrix
    virtual auto  copy         () const -> std::unique_ptr< TMatrix >;
    using TMatrix::copy;

    //! return structural copy of matrix
    virtual auto  copy_struct  () const -> std::unique_ptr< TMatrix >;

    //! copy matrix into matrix \a A
    virtual void  copy_to      ( TMatrix * A ) const;
    using TMatrix::copy_to;

    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( TSparseMatrix, TMatrix )

    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //
    // tests
    //

    //! test data for invalid values, e.g. INF and NAN
    virtual void check_data () const;

    //! print histogram for entries per row in GNUPLOT format
    virtual void print_pattern_hist ( std::ostream & os ) const;
};

//////////////////////////////////////////////////////////
//
// template wrappers
//

//
// wrapper for rcoeff/ccoeff
//
template <typename T>
T
coeff ( const TSparseMatrix *  S,
        const idx_t            i );

template <>
inline real
coeff<real>  ( const TSparseMatrix *  S,
               const idx_t            i )
{
    return S->rcoeff( i );
}

template <>
inline complex
coeff<complex>  ( const TSparseMatrix *  S,
                  const idx_t            i )
{
    return S->ccoeff( i );
}

}// namespace HLIB

#endif  // __HLIB_TSPARSEMATRIX_HH
