#ifndef __HLIB_TRKMATRIX_HH
#define __HLIB_TRKMATRIX_HH
//
// Project     : HLib
// File        : TRkMatrix.hh
// Description : class for rank-k-matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/base/TTruncAcc.hh"
#include "hpro/blas/Matrix.hh"
#include "hpro/matrix/TMatrix.hh"
#include "hpro/vector/TScalarVector.hh"

namespace HLIB
{

// forward decl.
class TDenseMatrix;

// local matrix type
DECLARE_TYPE( TRkMatrix );

//!
//! \ingroup Matrix_Module
//! \class   TRkMatrix
//! \brief   Represents low rank matrices in factored form: \f$ M = A B^H \f$.
//!
class TRkMatrix : public TMatrix
{
protected:
    //! @cond
    
    // factors of low-rank representation
    BLAS::Matrix< real >     _rmat_A, _rmat_B;
    BLAS::Matrix< complex >  _cmat_A, _cmat_B;

    // size of the vectors in A and B
    size_t                   _rows, _cols;
    
    // current rank of the matrix
    size_t                   _rank;
    
    //! @endcond
    
public:
    /////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct zero sized low-rank matrix
    TRkMatrix ();

    //! construct low-rank matrix of size \a rows × \a cols
    TRkMatrix ( const size_t                     rows,
                const size_t                     cols );

    //! construct low-rank matrix of size defined by block index set
    TRkMatrix ( const TBlockIndexSet &           block_is,
                const value_type_t               avalue_type = real_valued );
    
    //! construct low-rank matrix of size defined by block index set
    TRkMatrix ( const TIndexSet &                arow_is,
                const TIndexSet &                acol_is,
                const value_type_t               avalue_type = real_valued );
    
    //! construct low-rank matrix of size defined by block index set
    TRkMatrix ( const TIndexSet &                arow_is,
                const TIndexSet &                acol_is,
                const bool                       acomplex );
    
    //! construct low-rank matrix of size defined by block index set
    //! and real factors \a A and \a B
    TRkMatrix ( const TIndexSet &                arow_is,
                const TIndexSet &                acol_is,
                const BLAS::Matrix< real > &     A,
                const BLAS::Matrix< real > &     B );

    //! construct low-rank matrix of size defined by block index set
    //! and complex factors \a A and \a B
    TRkMatrix ( const TIndexSet &                arow_is,
                const TIndexSet &                acol_is,
                const BLAS::Matrix< complex > &  A,
                const BLAS::Matrix< complex > &  B );

    //! construct low-rank matrix of size defined by block index set
    //! and real factors \a A and \a B (move version)
    TRkMatrix ( const TIndexSet &                arow_is,
                const TIndexSet &                acol_is,
                BLAS::Matrix< real > &&          A,
                BLAS::Matrix< real > &&          B );

    //! construct low-rank matrix of size defined by block index set
    //! and complex factors \a A and \a B (move version)
    TRkMatrix ( const TIndexSet &                arow_is,
                const TIndexSet &                acol_is,
                BLAS::Matrix< complex > &&       A,
                BLAS::Matrix< complex > &&       B );

    //! construct low-rank matrix of size defined by block cluster
    TRkMatrix ( const TBlockCluster *            bc,
                const value_type_t               avalue_type = real_valued );

    //! copy constructor
    TRkMatrix ( const TRkMatrix &                A );

    //! destructor
    ~TRkMatrix ()
    {}

    /////////////////////////////////////////////////
    //
    // access local variables
    //

    //! return true, if matrix is zero
    virtual bool    is_zero   () const { return (( _rank == 0 ) && ! accumulator().has_updates()); }
    
    //
    // access actual data of matrix factors as BLAS::Matrix
    //
    
    //! return pointer to internal matrix data of matrix A
    BLAS::Matrix< real > &           blas_rmat_A  ()       { return _rmat_A; }
    
    //! return const pointer to internal matrix data of matrix A
    const BLAS::Matrix< real > &     blas_rmat_A  () const { return _rmat_A; }
    
    //! return pointer to internal matrix data of matrix A
    BLAS::Matrix< complex > &        blas_cmat_A  ()       { return _cmat_A; }
    
    //! return const pointer to internal matrix data of matrix A
    const BLAS::Matrix< complex > &  blas_cmat_A  () const { return _cmat_A; }

    //! return pointer to internal matrix data of matrix B
    BLAS::Matrix< real > &           blas_rmat_B  ()       { return _rmat_B; }
    
    //! return const pointer to internal matrix data of matrix B
    const BLAS::Matrix< real > &     blas_rmat_B  () const { return _rmat_B; }
    
    //! return pointer to internal matrix data of matrix B
    BLAS::Matrix< complex > &        blas_cmat_B  ()       { return _cmat_B; }
    
    //! return const pointer to internal matrix data of matrix B
    const BLAS::Matrix< complex > &  blas_cmat_B  () const { return _cmat_B; }

    //
    // access individual vectors in A and B as BLAS vectors
    //
    
    //! return vector A_i
    const BLAS::Vector< real >       blas_rvec_A  ( const idx_t  i ) const { return _rmat_A.column( i ); }

    //! return vector B_i
    const BLAS::Vector< real >       blas_rvec_B  ( const idx_t  i ) const { return _rmat_B.column( i ); }
    
    //! return vector A_i
    const BLAS::Vector< complex >    blas_cvec_A  ( const idx_t  i ) const { return _cmat_A.column( i ); }

    //! return vector B_i
    const BLAS::Vector< complex >    blas_cvec_B  ( const idx_t  i ) const { return _cmat_B.column( i ); }
    
    //
    // access individual vectors in A and B as H vectors
    //
    
    //! return vector A_i
    const TScalarVector             vec_A   ( const idx_t  i ) const
    {
        if ( is_complex() ) return TScalarVector( row_is(), _cmat_A.column( i ) );
        else                return TScalarVector( row_is(), _rmat_A.column( i ) );
    }

    //! return vector B_i
    const TScalarVector             vec_B   ( const idx_t  i ) const
    {
        if ( is_complex() ) return TScalarVector( col_is(), _cmat_B.column( i ) );
        else                return TScalarVector( col_is(), _rmat_B.column( i ) );
    }
    
    //
    // manage rank of matrix
    //
    
    //! return rank of matrix
    size_t          rank         () const { return _rank; }

    //! set rank of matrix without truncating data
    void            set_rank     ( const size_t  k );

    //! compute actual rank of matrix (remove zero vectors)
    void            comp_rank    ();

    //
    // access size information
    //
    
    //! return number of rows 
    virtual size_t  rows         () const { return _rows; }

    //! return number of columns 
    virtual size_t  cols         () const { return _cols; }

    //
    // change size
    //
    
    //! set new cluster over which matrix is defined and change size accordingly
    virtual void    set_cluster  ( const TBlockCluster * c );

    //! set size and rank of matrix (if zero == true fill new memory with zeros)
    void            set_size     ( const size_t  n,
                                   const size_t  m,
                                   const size_t  k );
    
    //! set new size but keep current rank of matrix
    void            set_size     ( const size_t  n,
                                   const size_t  m )
    {
        return set_size( n, m, _rank );
    }

    //
    // access size information
    //
    
    //! return real valued coefficent (\a i, \a j) of matrix
    real            entry        ( const idx_t  i, const idx_t j ) const;

    //! return complex valued coefficent (\a i, \a j) of matrix
    const complex   centry       ( const idx_t  i, const idx_t j ) const;

    //
    // change field type 
    //
    
    //! switch to real valued storage if possible, e.g. all imaginary components zero
    virtual void    to_real      ();

    //! switch to complex valued storage
    virtual void    to_complex   ();
    
    ///////////////////////////////////////////
    //
    // management of update accumulator
    //

    //! apply stored updates U to local matrix M, e.g., M = M + U,
    //! with accuracy \a acc
    virtual void    apply_updates ( const TTruncAcc &       acc,
                                    const recursion_type_t  recursion );
    
    /////////////////////////////////////////////////
    //
    // manage Rk-matrices
    //

    //! transpose matrix
    virtual void    transpose    ();
    
    //! conjugate matrix coefficients
    virtual void    conjugate    ();
    
    //! truncate matrix  w.r.t. accuracy \a acc
    virtual void    truncate     ( const TTruncAcc & acc );

    //! set this ≔ A·B^H
    void            set_lrmat    ( const BLAS::Matrix< real > &     A,
                                   const BLAS::Matrix< real > &     B );
    void            set_lrmat    ( const BLAS::Matrix< complex > &  A,
                                   const BLAS::Matrix< complex > &  B );
    
    void            set_lrmat    ( BLAS::Matrix< real > &&     A,
                                   BLAS::Matrix< real > &&     B )
    {
        if ( is_complex() )
            HERROR( ERR_ARG, "(TRkMatrix) set_lrmat", "expecting real valued data" );
        
        if (( A.nrows() != _rmat_A.nrows() ) ||
            ( B.nrows() != _rmat_B.nrows() ) ||
            ( A.ncols() != B.ncols()))
            HERROR( ERR_ARG, "(TRkMatrix) set_lrmat", "input matrices have invalid dimension" );

        if ( _rmat_A.ncols() == A.ncols() )
        {
            BLAS::copy( A, _rmat_A );
            BLAS::copy( B, _rmat_B );
        }// if
        else
        {
            _rank   = A.ncols();
            _rmat_A = std::move( A );
            _rmat_B = std::move( B );
        }// else
    }
    
    void            set_lrmat    ( BLAS::Matrix< complex > &&  A,
                                   BLAS::Matrix< complex > &&  B )
    {
        if ( ! is_complex() )
            HERROR( ERR_ARG, "(TRkMatrix) set_lrmat", "expecting complex valued data" );
        
        if (( A.nrows() != _cmat_A.nrows() ) ||
            ( B.nrows() != _cmat_B.nrows() ) ||
            ( A.ncols() != B.ncols()))
            HERROR( ERR_ARG, "(TRkMatrix) set_lrmat", "input matrices have invalid dimension" );
        
        if ( _cmat_A.ncols() == A.ncols() )
        {
            BLAS::copy( A, _cmat_A );
            BLAS::copy( B, _cmat_B );
        }// if
        else
        {
            _rank   = A.ncols();
            _cmat_A = std::move( A );
            _cmat_B = std::move( B );
        }// else
    }
    
    //! compute this ≔ this + α·A·B^H and truncate result w.r.t. \a acc (real valued variant)
    void            add_rank     ( const real                       alpha,
                                   const BLAS::Matrix< real > &     A,
                                   const BLAS::Matrix< real > &     B,
                                   const TTruncAcc &                acc );

    //! compute this ≔ this + α·A·B^H and truncate result w.r.t. \a acc (complex valued variant)
    void            add_rank     ( const complex                    alpha,
                                   const BLAS::Matrix< complex > &  A,
                                   const BLAS::Matrix< complex > &  B,
                                   const TTruncAcc &                acc );
    
    //! add a dense matrix and truncate w.r.t. \a acc (real valued variant)
    void            add_dense    ( const real                       alpha,
                                   const BLAS::Matrix< real > &     D,
                                   const TTruncAcc &                acc );

    //! add a dense matrix and truncate w.r.t. \a acc (complex valued variant)
    void            add_dense    ( const complex                    alpha,
                                   const BLAS::Matrix< complex > &  D,
                                   const TTruncAcc &                acc );
    
    //! copy a densematrix (nxm) as low-rank matrix (rank = min{n,m})
    void            copy_dense   ( const TDenseMatrix * A,
                                   const TTruncAcc &    acc );
    
    /////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    //! scale matrix by constant factor
    virtual void         scale      ( const real       f );
    
    //! compute y ≔ α·M·x + β·y, where M is either this, this^T or this^H depending on \a op
    virtual void         mul_vec    ( const real       alpha,
                                      const TVector *  x,
                                      const real       beta,
                                      TVector *        y,
                                      const matop_t    op = apply_normal ) const;

    //! compute this = this + α·A without truncation
    virtual void         add        ( const real       alpha,
                                      const TMatrix *  A );
        
    //! compute matrix product α·op_A(this)·op_B(B) 
    virtual TRkMatrix *  mul_right  ( const real       alpha,
                                      const TMatrix *  B,
                                      const matop_t    op_A,
                                      const matop_t    op_B ) const;

    //! compute matrix product α·op_A(A)·op_B(this) 
    virtual TRkMatrix *  mul_left   ( const real       alpha,
                                      const TMatrix *  A,
                                      const matop_t    op_A,
                                      const matop_t    op_B ) const;

    /////////////////////////////////////////////////
    //
    // BLAS-routines (complex valued)
    //

    //! scale matrix by constant factor
    virtual void        cscale      ( const complex    f );
    
    //! compute y ≔ α·M·x + β·y, where M is either this, this^T or this^H depending on \a op
    virtual void        cmul_vec    ( const complex    alpha,
                                      const TVector *  x,
                                      const complex    beta,
                                      TVector *        y,
                                      const matop_t    op_A = apply_normal ) const;

    //! compute this = this + α·A without truncation
    virtual void        cadd        ( const complex    a,
                                      const TMatrix *  matrix );
        
    //! compute matrix product α·op_A(this)·op_B(B) 
    virtual TRkMatrix * cmul_right  ( const complex    alpha,
                                      const TMatrix *  B,
                                      const matop_t    op_A,
                                      const matop_t    op_B ) const;
    
    //! compute matrix product α·op_A(A)·op_B(this) 
    virtual TRkMatrix * cmul_left   ( const complex    alpha,
                                      const TMatrix *  A,
                                      const matop_t    op_A,
                                      const matop_t    op_B ) const;

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
    // misc methods
    //

    //
    // virtual constructors
    //

    //! return object of same type
    virtual auto   create       () const -> std::unique_ptr< TMatrix >
    {
        return std::make_unique< TRkMatrix >();
    }

    //! return copy of matrix
    virtual auto   copy         () const -> std::unique_ptr< TMatrix >;

    //! return copy matrix wrt. given accuracy; if \a do_coarsen is set, perform coarsening
    virtual auto   copy         ( const TTruncAcc &  acc,
                                  const bool         do_coarsen = false ) const -> std::unique_ptr< TMatrix >;

    //! return structural copy of matrix
    virtual auto   copy_struct  () const -> std::unique_ptr< TMatrix >;

    // copy matrix data to \a A
    virtual void   copy_to      ( TMatrix *          A ) const;

    // copy matrix data to \a A and truncate w.r.t. \acc with optional coarsening
    virtual void   copy_to      ( TMatrix *          A,
                                  const TTruncAcc &  acc,
                                  const bool         do_coarsen = false ) const;
    
    //! return size in bytes used by this object
    virtual size_t byte_size    () const;

    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( TRkMatrix, TMatrix )

    //
    // serialisation
    //

    //! read matrix from byte stream
    virtual void read  ( TByteStream & s );

    //! construct matrix from byte stream
    virtual void build ( TByteStream & s );
    
    //! write matrix into byte stream
    virtual void write ( TByteStream & s ) const;

    //! return size of object in a bytestream
    virtual size_t bs_size () const;

    //
    // tests
    //

    //! test data for invalid values, e.g. INF and NAN
    virtual void check_data () const;
};

//
// wrapper for matrix_A/B and cmatrix_A/B to avoid
// if clauses ( if ( is_complex() ) )
//

template <typename T>  BLAS::Matrix< T > &        blas_mat_A  ( TRkMatrix *        A );
template <typename T>  const BLAS::Matrix< T > &  blas_mat_A  ( const TRkMatrix *  A );
template <typename T>  BLAS::Matrix< T > &        blas_mat_B  ( TRkMatrix *        A );
template <typename T>  const BLAS::Matrix< T > &  blas_mat_B  ( const TRkMatrix *  A );

template <> inline
BLAS::Matrix< real > &
blas_mat_A< real >    ( TRkMatrix *  A )
{ return A->blas_rmat_A();  }

template <> inline
BLAS::Matrix< complex > &
blas_mat_A< complex > ( TRkMatrix *  A )
{ return A->blas_cmat_A(); }

template <> inline
BLAS::Matrix< real > &
blas_mat_B< real >    ( TRkMatrix *  A )
{ return A->blas_rmat_B();  }

template <> inline
BLAS::Matrix< complex > &
blas_mat_B< complex > ( TRkMatrix *  A )
{ return A->blas_cmat_B(); }

template <> inline
const BLAS::Matrix< real > &
blas_mat_A< real >    ( const TRkMatrix *  A )
{ return A->blas_rmat_A();  }

template <> inline
const BLAS::Matrix< complex > &
blas_mat_A< complex > ( const TRkMatrix *  A )
{ return A->blas_cmat_A(); }

template <> inline
const BLAS::Matrix< real > &
blas_mat_B< real >    ( const TRkMatrix *  A )
{ return A->blas_rmat_B();  }

template <> inline
const BLAS::Matrix< complex > &
blas_mat_B< complex > ( const TRkMatrix *  A )
{ return A->blas_cmat_B(); }

template <typename T>  BLAS::Matrix< T > &        blas_mat_A  ( TRkMatrix &        A );
template <typename T>  const BLAS::Matrix< T > &  blas_mat_A  ( const TRkMatrix &  A );
template <typename T>  BLAS::Matrix< T > &        blas_mat_B  ( TRkMatrix &        A );
template <typename T>  const BLAS::Matrix< T > &  blas_mat_B  ( const TRkMatrix &  A );

template <> inline
BLAS::Matrix< real > &
blas_mat_A< real >    ( TRkMatrix &  A )
{ return A.blas_rmat_A();  }

template <> inline
BLAS::Matrix< complex > &
blas_mat_A< complex > ( TRkMatrix &  A )
{ return A.blas_cmat_A(); }

template <> inline
BLAS::Matrix< real > &
blas_mat_B< real >    ( TRkMatrix &  A )
{ return A.blas_rmat_B();  }

template <> inline
BLAS::Matrix< complex > &
blas_mat_B< complex > ( TRkMatrix &  A )
{ return A.blas_cmat_B(); }

template <> inline
const BLAS::Matrix< real > &
blas_mat_A< real >    ( const TRkMatrix &  A )
{ return A.blas_rmat_A();  }

template <> inline
const BLAS::Matrix< complex > &
blas_mat_A< complex > ( const TRkMatrix &  A )
{ return A.blas_cmat_A(); }

template <> inline
const BLAS::Matrix< real > &
blas_mat_B< real >    ( const TRkMatrix &  A )
{ return A.blas_rmat_B();  }

template <> inline
const BLAS::Matrix< complex > &
blas_mat_B< complex > ( const TRkMatrix &  A )
{ return A.blas_cmat_B(); }

template <typename T>  BLAS::Matrix< T > &  blas_mat_A  ( std::unique_ptr< TRkMatrix > & M ) { return blas_mat_A<T>( M.get() ); }
template <typename T>  BLAS::Matrix< T > &  blas_mat_B  ( std::unique_ptr< TRkMatrix > & M ) { return blas_mat_B<T>( M.get() ); }

}// namespace HLIB

#endif  // __HLIB_TRKMATRIX_HH
