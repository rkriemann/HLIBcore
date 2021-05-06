#ifndef __HLIB_THMATRIX_HH
#define __HLIB_THMATRIX_HH
//
// Project     : HLib
// File        : THMatrix.hh
// Description : class for H-matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <list>

#include "hpro/cluster/TPermutation.hh"

#include "hpro/matrix/TBlockMatrix.hh"

#include "hpro/parallel/NET.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( THMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    THMatrix
//! \brief    Class for an H-matrix, which extends block matrices
//!           with additional functionality, e.g. permutations.
//!
class THMatrix : public TBlockMatrix
{
private:

    //! @cond

    // maximal ID of all subblocks
    int            _max_id;
    
    // mappings of indices from external to
    // internal numbering (for mul_vec)
    TPermutation   _row_perm_e2i, _col_perm_e2i;
    TPermutation   _row_perm_i2e, _col_perm_i2e;

    //! @endcond
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    THMatrix ( const TBlockCluster *   cluster = nullptr );

    THMatrix ( const TBlockIndexSet &  bis );

    virtual ~THMatrix () {}

    /////////////////////////////////////////////////
    //
    // access internal data
    //

    // access single matrix coefficient
    virtual real          entry  ( const idx_t i, const idx_t j ) const;
    virtual const complex centry ( const idx_t i, const idx_t j ) const;

    // compute min/max tau/sigma indices
    void comp_min_max_idx ();

    // access permutations
    void set_row_perm ( const TPermutation & perm_e2i,
                        const TPermutation & perm_i2e )
    {
        _row_perm_e2i = perm_e2i;
        _row_perm_i2e = perm_i2e;
    }

    void set_col_perm ( const TPermutation & perm_e2i,
                        const TPermutation & perm_i2e )
    {
        _col_perm_e2i = perm_e2i;
        _col_perm_i2e = perm_i2e;
    }

    const TPermutation & row_perm_e2i () const { return _row_perm_e2i; }
    const TPermutation & row_perm_i2e () const { return _row_perm_i2e; }
    const TPermutation & col_perm_e2i () const { return _col_perm_e2i; }
    const TPermutation & col_perm_i2e () const { return _col_perm_i2e; }

    // return true if _all_ mappings are present and of correct size
    bool has_perm () const
    {
        return (( _row_perm_e2i.size() == rows() ) &&
                ( _row_perm_i2e.size() == rows() ) &&
                ( _col_perm_e2i.size() == cols() ) &&
                ( _col_perm_i2e.size() == cols() ));
    }

    // return maximal ID of subblocks
    int  max_id () const { return _max_id; }
    
    /////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    //! compute y ≔ β·y + α·op(M)·x, with M = this
    virtual void mul_vec ( const real      alpha,
                           const TVector * x,
                           const real      beta,
                           TVector       * y,
                           const matop_t   op = MATOP_NORM ) const;

    /////////////////////////////////////////////////
    //
    // BLAS-routines (complex valued)
    //

    //! compute y ≔ β·y + α·op(M)·x, with M = this
    virtual void cmul_vec ( const complex   alpha,
                            const TVector * x,
                            const complex   beta,
                            TVector       * y,
                            const matop_t   op = MATOP_NORM ) const;
    
    /////////////////////////////////////////////////
    //
    // matrix operations
    //

#if 0
    // truncate all rank-k-blocks in matrix to given accuracy
    virtual void truncate ( const TTruncAcc & acc );
#endif
    
    //! transpose matrix
    virtual void transpose ();
    
    /////////////////////////////////////////////////
    //
    // misc.
    //

    //
    // virtual constructor
    //

    // return matrix of same class (but no content)
    virtual auto  create       () const -> std::unique_ptr< TMatrix >
    {
        return std::make_unique< THMatrix >();
    }
    
    // return copy of matrix
    virtual auto  copy         () const -> std::unique_ptr< TMatrix >;
    
    // copy matrix wrt. given accuracy and coarsening
    virtual auto  copy         ( const TTruncAcc &     acc,
                                 const bool            coarsen = false ) const -> std::unique_ptr< TMatrix >;

    // return structural copy of matrix
    virtual auto  copy_struct  () const -> std::unique_ptr< TMatrix >;
    
    // copy matrix into A
    virtual void  copy_to      ( TMatrix *          A ) const;
    virtual void  copy_to      ( TMatrix *          A,
                                 const TTruncAcc &  acc,
                                 const bool         coarsen = false ) const;

    // copy complete structural information from given matrix
    virtual void  copy_struct_from ( const TMatrix * M );
    
    // return appropriate vector-types for matrix
    virtual auto  row_vector   () const -> std::unique_ptr< TVector >;
    virtual auto  col_vector   () const -> std::unique_ptr< TVector >;
    
    //
    // size of object
    //
    
    // return size in bytes used by this object
    virtual size_t byte_size () const;

    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( THMatrix, TBlockMatrix )
    
    //
    // serialisation
    //

    virtual void read  ( TByteStream & s );
    virtual void build ( TByteStream & s );
    virtual void write ( TByteStream & s ) const;

    // returns size of object in bytestream
    virtual size_t bs_size () const;
};

}// namespace

#endif  // __HLIB_THMATRIX_HH
