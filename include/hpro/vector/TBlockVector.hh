#ifndef __HLIB_TBLOCKVECTOR_HH
#define __HLIB_TBLOCKVECTOR_HH
//
// Project     : HLib
// File        : TBlockVector.hh
// Description : class for a vector build out of other vectors
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/base/System.hh"
#include "hpro/vector/TVector.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TBlockVector );

//!
//! \ingroup Vector_Module
//! \class   TBlockVector
//! \brief   Class for a blocked, scalar vector
//!
class TBlockVector : public TVector
{
private:
    //! individual vector blocks
    std::vector< TVector * >  _blocks;

public:
    
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct block vector with \a nb sub blocks
    TBlockVector ( const uint  nb = 0 )
    {
        _blocks.resize( nb, NULL );
    }

    //! construct block vector over index set \a ais and
    //! with subblocks \a ablocks
    TBlockVector ( const TIndexSet &                 ais,
                   const std::vector< TVector * > &  ablocks )
            : TVector( ais.first() )
            , _blocks( ablocks )
    {}

    //! dtor
    virtual ~TBlockVector ();

    ///////////////////////////////////////////////
    //
    // access internal data
    //

    //! return size of vector
    virtual size_t   size              () const;

    //! return number of blocks
    virtual uint     n_blocks          () const { return uint(_blocks.size()); }

    //! access single vector block
    TVector *        block             ( const uint i )       { return _blocks[i]; }

    //! access single vector block
    const TVector *  block             ( const uint i ) const { return _blocks[i]; }

    //! set single vector block
    void             set_block         ( const uint i, TVector * v );
    
    //! setup block structure of vector
    void             set_block_struct  ( const uint i );

protected:

    //! convert data to real valued representation
    virtual void     to_real    ();

    //! convert data to complex valued representation
    virtual void     to_complex ();

public:
    ////////////////////////////////////////////////
    //
    // access entries
    //

    //! return i'th entry
    virtual real          entry  ( const idx_t  i ) const;

    //! return i'th entry
    virtual const complex centry ( const idx_t  i ) const;

    //! set i'th entry
    virtual void set_entry  ( const idx_t  i, const real      f );

    //! set i'th entry
    virtual void set_centry ( const idx_t  i, const complex   f );

    //////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    //! fill with constant
    virtual void fill ( const real f );
    
    //! fill with random numbers
    virtual void fill_rand ( const uint seed );

    //! set this ≔ this + α · x
    virtual void axpy ( const real alpha, const TVector * x );

    //! set this ≔ α · this
    virtual void scale ( const real alpha );

    //! set this ≔ α · x
    virtual void assign ( const real alpha, const TVector * x );

    //! compute ‖·‖₂
    virtual real norm2 () const { return Math::sqrt( dot( this ).real() ); }

    //! compute ‖·‖∞
    virtual real norm_inf () const;

    //////////////////////////////////////////////////
    //
    // BLAS-routines (complex valued)
    //

    //! conjugate entries
    virtual void conjugate ();
        
    //! fill with constant
    virtual void cfill     ( const complex & f );
    
    //! set this ≔ α · this
    virtual void cscale    ( const complex & alpha );

    //! set this ≔ α · x
    virtual void cassign   ( const complex & alpha, const TVector * x );

    //! set this ≔ this + α · x
    virtual void caxpy     ( const complex & f, const TVector * x );

    // return dot product <this,x> = this^H x
    virtual complex dot    ( const TVector * x ) const;

    // return dot product <this,x> = this^T x
    virtual complex dotu   ( const TVector * x ) const;

    ///////////////////////////////////////////////
    //
    // misc.
    //

    //
    // memory consumption
    //
    
    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //
    // virtual constructor
    //

    //! return copy of vector
    virtual auto  copy  () const -> std::unique_ptr< TVector >;
    
    //! return object of same class
    virtual auto  create () const -> std::unique_ptr< TVector > { return std::make_unique< TBlockVector >(); }

    //
    // output
    //

    //! print vector to stdout
    virtual void print ( const uint ofs = 0 ) const;

    
    HLIB_RTTI_DERIVED( TBlockVector, TVector )
};

}// namespace

#endif  // __HLIB_TBLOCKVECTOR_HH
