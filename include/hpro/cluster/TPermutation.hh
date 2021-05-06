#ifndef __HLIB_TPERMUTATION_HH
#define __HLIB_TPERMUTATION_HH
//
// Project     : HLib
// File        : TPermutation.hh
// Description : class for a permutation
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/base/types.hh"
#include "hpro/base/TByteStream.hh"

#include "hpro/cluster/TIndexSet.hh"

#include "hpro/vector/TVector.hh"

namespace HLIB
{

//!
//! \class  TPermutation
//! \brief  Describes permutation of index sets.
//!
class TPermutation : public std::vector< idx_t >
{
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct permutation of size \a asize
    TPermutation ( const size_t                  asize = 0 );

    //! construct permutation with data defined by array \a perm
    TPermutation ( const std::vector< idx_t > &  perm );

    //! copy ctor
    TPermutation ( const TPermutation &          perm );

    //! dtor
    virtual ~TPermutation ();

    ///////////////////////////////////////////
    //
    // misc. methods for permutation handling
    //

    //! permute given index
    idx_t         permute      ( const idx_t      idx ) const { return (*this)[ idx ]; }

    //! reorder given vector \a x and write result to \a y
    void          permute      ( const TVector *  x,
                                 TVector *        y ) const;
    
    //! use inverse permutation to reorder given vector \a x and write result to \a y
    void          permute_inv  ( const TVector *  x,
                                 TVector *        y ) const;
    
    //! permute given vector \x inplace
    void          permute      ( TVector *        x ) const;
    
    //! apply inverse permutation to given vector \x
    void          permute_inv  ( TVector *        x ) const;
    
    //! invert permutation
    void          invert   ();
    
    //! return inverse permutation
    TPermutation  inverse ();
    
    /////////////////////////////////////////////////
    //
    // operators
    //

    idx_t  operator () ( const idx_t  idx ) const { return permute( idx ); }
    
    /////////////////////////////////////////////////
    //
    // serialisation
    //

    //! read object data from bytestream \a s
    virtual void read  ( TByteStream & s );

    //! write object data to bytestream \a s
    virtual void write ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;

    /////////////////////////////////////////////////
    //
    // misc.
    //

    //! copy operator
    TPermutation &  operator =  ( const TPermutation &  perm )
    {
        resize( perm.size() );
        
        for ( size_t  i = 0; i < size(); i++ )
            (*this)[i] = perm[i];

        return *this;
    }
    
    //! copy operator for vector
    TPermutation &  operator =  ( const std::vector< idx_t > &  perm )
    {
        resize( perm.size() );
        
        for ( size_t  i = 0; i < size(); i++ )
            (*this)[i] = perm[i];

        return *this;
    }
    
    //! return size in bytes used by this object
    virtual size_t byte_size () const;
};

}// namespace HLIB

#endif  // __HLIB_TPERMUTATION_HH
