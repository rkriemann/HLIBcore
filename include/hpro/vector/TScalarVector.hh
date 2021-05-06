#ifndef __HLIB_TSCALARVECTOR_HH
#define __HLIB_TSCALARVECTOR_HH
//
// Project     : HLib
// File        : TScalarVector.hh
// Description : class for a vector of scalar type
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <mutex>

#include "hpro/base/types.hh"
#include "hpro/blas/Vector.hh"
#include "hpro/blas/Algebra.hh"

#include "hpro/cluster/TCluster.hh"
#include "hpro/cluster/TPermutation.hh"

#include "hpro/vector/TVector.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TScalarVector );

// chunk size for updates to scalar vectors
const size_t  SCALAR_CHUNK_SIZE = 64;
    
//!
//! \ingroup Vector_Module
//! \class   TScalarVector 
//! \brief   Class for a scalar vector.
//!
class TScalarVector : public TVector
{
private:
    // local mutex and lock types
    using  mutex_t = std::mutex;
    using  lock_t  = std::lock_guard< mutex_t >;
    
protected:
    //! @cond
    
    //! real valued vector data
    BLAS::Vector< real >     _rvec;
    
    //! complex valued vector data
    BLAS::Vector< complex >  _cvec;

    //! size of vector
    size_t                   _size;

    // //! mutices for each chunk of vector
    // std::vector< mutex_t >   _mutices;

    //! @endcond
    
public:
    //////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //
    // constructors without copy/move
    //
    
    //! construct zero sized vector
    TScalarVector ( const value_type_t               avalue_type = real_valued )
            : TVector( 0, avalue_type )
            , _size( 0 ) 
    {}
    TScalarVector ( const bool                       acomplex )
            : TVector( 0, acomplex ? complex_valued : real_valued )
            , _size( 0 ) 
    {}
    
    //! construct vector of size \a n with offset \a offset
    TScalarVector ( const size_t                     n,
                    const idx_t                      offset      = 0,
                    const value_type_t               avalue_type = real_valued )
            : TVector( offset, avalue_type )
            , _size( 0 )
    {
        set_size( n );
    }
    TScalarVector ( const size_t                     n,
                    const idx_t                      offset,
                    const bool                       acomplex )
            : TVector( offset, acomplex ? complex_valued : real_valued )
            , _size( 0 )
    {
        set_size( n );
    }
    
    //! construct vector with size defined by indexset \a ais
    TScalarVector ( const TIndexSet &                ais,
                    const value_type_t               avalue_type = real_valued )
            : TVector( ais.first(), avalue_type )
            , _size( 0 )
    {
        set_size( ais.size() );
    }
    TScalarVector ( const TIndexSet &                ais,
                    const bool                       acomplex )
            : TVector( ais.first(), acomplex ? complex_valued : real_valued )
            , _size( 0 )
    {
        set_size( ais.size() );
    }
    
    //
    // copy constructors
    //
    
    //! construct vector with size defined by indexset \a ais
    //! and data defined by \a bvec
    TScalarVector ( const TIndexSet &                ais,
                    const BLAS::Vector< real > &     bvec )
            : TVector( ais.first(), real_valued )
            , _rvec( bvec )
            , _size( bvec.length() )
    {
        init_chunk_mutices();
    }
    
    //! construct vector with size defined by indexset \a ais
    //! and data defined by \a bvec
    TScalarVector ( const TIndexSet &                ais,
                    const BLAS::Vector< complex > &  bvec )
            : TVector( ais.first(), complex_valued )
            , _cvec( bvec )
            , _size( bvec.length() )
    {
        init_chunk_mutices();
    }

    //! standard copy constructor
    TScalarVector ( const TScalarVector &            v )
            : TVector()
            , _size( 0 )
    {
        assign( 1.0, & v );
    }

    //! standard move constructor
    TScalarVector ( TScalarVector &&                 v )
            : TVector( v.ofs(), v.value_type() )
            , _size( 0 )
    {
        if ( v.is_complex() )
        {
            _size = v._cvec.length();
            _cvec = std::move( v._cvec );
        }// if
        else
        {
            _size = v._rvec.length();
            _rvec = std::move( v._rvec );
        }// else
    }

    //
    // move constructors
    //

    //! construct vector with size defined by indexset \a ais
    //! and data defined by \a bvec
    TScalarVector ( const TIndexSet &           ais,
                    BLAS::Vector< real > &&     bvec )
            : TVector( ais.first(), real_valued )
            , _rvec( std::move( bvec ) )
            , _size( _rvec.length() )  // use _rvec, because bvec is now "empty"
    {
        init_chunk_mutices();
    }
    
    //! construct vector with size defined by indexset \a ais
    //! and data defined by \a bvec
    TScalarVector ( const TIndexSet &           ais,
                    BLAS::Vector< complex > &&  bvec )
            : TVector( ais.first(), complex_valued )
            , _cvec( std::move( bvec ) )
            , _size( _cvec.length() )
    {
        init_chunk_mutices();
    }

    //! destructor
    virtual ~TScalarVector ()
    {}

    //////////////////////////////////////////////////
    //
    // access vector coefficients
    //

    //! return size of vector
    virtual size_t          size        () const
    {
        return _size;
    }

    //! access coefficent \a i (real valued)
    virtual real            entry       ( const idx_t  i ) const
    {
        HASSERT( i < idx_t(_size), ERR_ARR_BOUND, "(TScalarVector) entry", "" );
        return _rvec(i);
    }
    
    //! access coefficent \a i (complex valued)
    virtual const complex   centry      ( const idx_t  i ) const
    {
        HASSERT( i < idx_t(_size), ERR_ARR_BOUND, "(TScalarVector) centry", "" );
        return _cvec(i);
    }
    
    //! set coefficient \a i to \a f (real valued)
    virtual void            set_entry   ( const idx_t  i, const real  f )
    {
        HASSERT( i < idx_t(_size), ERR_ARR_BOUND, "(TScalarVector) set_entry", "" );
        _rvec(i) = f;
    }
    
    //! set coefficient \a i to \a f (complex valued)
    virtual void            set_centry  ( const idx_t  i, const complex  f )
    {
        HASSERT( i < idx_t(_size), ERR_ARR_BOUND, "(TScalarVector) set_centry", "" );
        _cvec(i) = f;
    }
    
    //! add \a f to \a i'th entry
    virtual void            add_entry   ( const idx_t  i, const real  f )
    {
        HASSERT( i < idx_t(_size), ERR_ARR_BOUND, "(TScalarVector) add_entry", "" );
        _rvec(i) += f;
    }
    
    //! add \a f to \a i'th entry (complex valued)
    virtual void            add_centry  ( const idx_t  i, const complex  f )
    {
        HASSERT( i < idx_t(_size), ERR_ARR_BOUND, "(TScalarVector) add_centry", "" );
        _cvec(i) += f;
    }

    //
    // access internal vectors as BLAS vectors
    //

    //! return real valued data
    BLAS::Vector< real > &           blas_rvec  ()       { return _rvec; }

    //! return constant real valued data
    const BLAS::Vector< real > &     blas_rvec  () const { return _rvec; }
    
    //! return complex valued data
    BLAS::Vector< complex > &        blas_cvec  ()       { return _cvec; }

    //! return constant complex valued data
    const BLAS::Vector< complex > &  blas_cvec  () const { return _cvec; }

    //
    // handle index set
    //
    
    //! set size of vector
    virtual void set_size    ( const size_t  n );
    
    //! define vector by cluster
    virtual void set_cluster ( const TCluster *  c )
    {
        if ( c != NULL )
        {
            set_size( c->size() );
            set_ofs( c->first() );
        }// if
    }

    //! define vector by indexset
    virtual void set_is ( const TIndexSet &  ais )
    {
        set_size( ais.size() );
        set_ofs( ais.first() );
    }

    //
    // copy/assign methods
    //
    
    //! set internal data directly (real valued)
    virtual void set_vector ( const BLAS::Vector< real > &  vec,
                              const idx_t                   offset )
    {
        set_complex( false );
        set_size( vec.length() );
        set_ofs( offset );
        BLAS::copy( vec, _rvec );
    }

    //! set internal data directly (complex valued)
    virtual void set_vector ( const BLAS::Vector< complex > &  vec,
                              const idx_t                      offset )
    {
        set_complex( true );
        set_size( vec.length() );
        set_ofs( offset );
        BLAS::copy( vec, _cvec );
    }

    //! copy from vector \a v
    virtual void copy_from ( const TScalarVector * v )
    {
        if ( v == NULL )
            HERROR( ERR_ARG, "(TScalarVector) copy_from", "v is NULL" );

        if ( v->is_complex() )
            set_vector( v->blas_cvec(), v->ofs() );
        else
            set_vector( v->blas_rvec(), v->ofs() );
    }

    //! copy to vector \a v
    virtual void copy_to ( TScalarVector * v ) const
    {
        if ( v == NULL )
            HERROR( ERR_ARG, "(TScalarVector) copy_to", "v is NULL" );

        v->copy_from( this );
    }

    //
    // return reference to sub vector
    //

    TScalarVector  sub_vector ( const TIndexSet &  ais )
    {
        if ( ! is().is_sub( ais ) )
            HERROR( ERR_INDEXSET, "(TScalarVector) sub_vector",
                    "given index set is NOT a sub set of local index set" );

        if ( is_complex() )
        {
            return TScalarVector( ais, BLAS::Vector< complex >( _cvec, ais - ofs() ) );
        }// if
        else
        {
            return TScalarVector( ais, BLAS::Vector< real >( _rvec, ais - ofs() ) );
        }// else
    }
    
    const TScalarVector  sub_vector ( const TIndexSet &  ais ) const
    {
        if ( ! is().is_sub( ais ) )
            HERROR( ERR_INDEXSET, "(TScalarVector) sub_vector",
                    "given index set is NOT a sub set of local index set" );

        if ( is_complex() )
        {
            return TScalarVector( ais, BLAS::Vector< complex >( _cvec, ais - ofs() ) );
        }// if
        else
        {
            return TScalarVector( ais, BLAS::Vector< real >( _rvec, ais - ofs() ) );
        }// else
    }
    
    
    //! copy from C array \a v
    virtual void copy_from ( const real * v );

    //! copy to C array \a v
    virtual void copy_to   ( real *       v );
    using TVector::copy_to;

    //! standard copy operator
    TScalarVector &  operator = ( const TScalarVector &  v )
    {
        assign( 1.0, & v );
        return * this;
    }

    //! permute entries according to \a perm
    void permute ( const TPermutation &  perm );

protected:

    //
    // switch between real and complex
    //
    
    //! switch to real valued representation if possible
    virtual void to_real    ();

    //! switch to complex valued representation
    virtual void to_complex ();

    //
    // handle chunk mutex array
    //

    //! initialise mutices for vector chunks
    void init_chunk_mutices ()
    {
        // _mutices = std::move( std::vector< mutex_t >( size() / SCALAR_CHUNK_SIZE + 1 ) );
        // // _mutices.resize( size() / SCALAR_CHUNK_SIZE + 1 ); // +1 for non-multiples
    }

public:
    //////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    //! fill vector with constant \a α
    virtual void   fill        ( const real alpha );

    //! fill vector with random numbers
    virtual void   fill_rand   ( const uint seed );

    //! set this ≔ α · this
    virtual void   scale       ( const real alpha );

    //! set this ≔ α · x
    virtual void   assign      ( const real alpha, const TVector * x );

    //! compute ‖·‖₂
    virtual real   norm2       () const;

    //! compute ‖·‖∞
    virtual real   norm_inf    () const;
    
    //! set this ≔ this + α · x
    virtual void   axpy        ( const real alpha, const TVector * x );

    //! set this ≔ this + x (thread safe, is(x) ⊆ is(this))
    virtual void   add_sub_mt  ( const TScalarVector &  x );

    //////////////////////////////////////////////////
    //
    // BLAS-routines (complex valued)
    //

    //! conjugate coefficients
    virtual void     conjugate  ()
    {
        if ( is_complex() )
            BLAS::conj( blas_cvec() );
    }
    
    //! fill vector with constant \a α
    virtual void     cfill      ( const complex &  alpha );

    //! set this ≔ α · this
    virtual void     cscale     ( const complex &  alpha );

    //! set this ≔ α · x
    virtual void     cassign    ( const complex &  alpha, const TVector * x );

    //! return inner product <this, x> = this^H · x
    virtual complex  dot        ( const TVector * x ) const;

    //! return inner product <this, x> = this^T · x
    virtual complex  dotu       ( const TVector * x ) const;

    //! set this ≔ this + α · x
    virtual void     caxpy      ( const complex & f, const TVector * x );

    //////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //
    // virtual constructor
    //

    //! return copy of vector
    virtual auto  copy   () const -> std::unique_ptr< TVector > { return std::make_unique< TScalarVector >( *this ); }
    
    //! return object of same class
    virtual auto  create () const -> std::unique_ptr< TVector > { return std::make_unique< TScalarVector >(); }
    
    //
    // restriction
    //

    //! return vector restricted to real part of coefficients
    virtual auto restrict_re  () const -> std::unique_ptr< TVector >;

    //! return vector restricted to imaginary part of coefficients
    virtual auto restrict_im  () const -> std::unique_ptr< TVector >;
    
    //! return vector restricted to absolute value of coefficients
    virtual auto restrict_abs () const -> std::unique_ptr< TVector >;
    
    //
    // stream output
    //

    //! print vector information
    virtual void print ( const uint ofs = 0 ) const;

    //
    // serialisation
    //

    //! read vector from stream
    virtual void read  ( TByteStream & s );

    //! write vector to stream
    virtual void write ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;

    //
    // parallel methods
    //

    //! pointwise summation between all vectors in \a ps
    virtual void sum ( const TProcSet & ps );
    using TVector::sum;

    HLIB_RTTI_DERIVED( TScalarVector, TVector )
};

//////////////////////////////////////////////////////////
//
// template wrappers
//

//
// blas_rvec/_cvec const
//
template <typename T> BLAS::Vector< T > &  blas_vec  ( TScalarVector *  v ); 
template <typename T> BLAS::Vector< T > &  blas_vec  ( TScalarVector &  v );

template <> inline BLAS::Vector< real > &     blas_vec<real>     ( TScalarVector *  v ) { return v->blas_rvec(); }
template <> inline BLAS::Vector< complex > &  blas_vec<complex>  ( TScalarVector *  v ) { return v->blas_cvec(); }

template <> inline BLAS::Vector< real > &     blas_vec<real>     ( TScalarVector &  v ) { return v.blas_rvec(); }
template <> inline BLAS::Vector< complex > &  blas_vec<complex>  ( TScalarVector &  v ) { return v.blas_cvec(); }

//
// blas_rvec/_cvec const
//
template <typename T> const BLAS::Vector< T > &  blas_vec  ( const TScalarVector *  v );
template <typename T> const BLAS::Vector< T > &  blas_vec  ( const TScalarVector &  v );

template <> inline const BLAS::Vector< real > &     blas_vec<real>     ( const TScalarVector *  v ) { return v->blas_rvec(); }
template <> inline const BLAS::Vector< complex > &  blas_vec<complex>  ( const TScalarVector *  v ) { return v->blas_cvec(); }

template <> inline const BLAS::Vector< real > &     blas_vec<real>     ( const TScalarVector &  v ) { return v.blas_rvec(); }
template <> inline const BLAS::Vector< complex > &  blas_vec<complex>  ( const TScalarVector &  v ) { return v.blas_cvec(); }

//!
//! functional version of TScalarVector::sub_vector
//!
inline
TScalarVector
sub_vector ( TScalarVector *    v,
             const TIndexSet &  is )
{
    return v->sub_vector( is );
}
             
inline
const TScalarVector
sub_vector ( const TScalarVector *  v,
             const TIndexSet &      is )
{
    return v->sub_vector( is );
}

//
// type checks
//

inline bool is_scalar ( const TVector &  v ) noexcept { return IS_TYPE( & v, TScalarVector ); }
inline bool is_scalar ( const TVector *  v ) noexcept { return ( v != nullptr ) && IS_TYPE( v, TScalarVector ); }

}// namespace HLIB

#endif  // __HLIB_TSCALARVECTOR_HH
