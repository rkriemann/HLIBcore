#ifndef __HLIB_TVECTOR_HH
#define __HLIB_TVECTOR_HH
//
// Project     : HLib
// File        : TVector.hh
// Description : baseclass for all vector-classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/base/System.hh"
#include "hpro/base/TStreamable.hh"
#include "hpro/base/TTypeInfo.hh"
#include "hpro/cluster/TIndexSet.hh"

namespace HLIB
{

//////////////////////////////////////////////////////////
//
// base vector class
//

// local matrix type
DECLARE_TYPE( TVector );

//!
//! \ingroup Vector_Module
//! \class   TVector
//! \brief   Base class for all vectors defining basic interface.
//!
class TVector : public TStreamable, public TTypeInfo
{
protected:
    //@cond
    
    //! first index vector represents
    idx_t  _ofs;

    //! flag indicating real of complex values
    bool   _complex;
    
    //@endcond
    
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct real or complex valued vector with first index \a offset
    TVector ( const idx_t         offset      = 0,
              const value_type_t  avalue_type = real_valued )
            : _ofs(offset)
            , _complex( avalue_type == complex_valued )
    {}

    //! copy constructor
    TVector ( const TVector &  v )
            : _ofs(0)
            , _complex(false)
    {
        assign( 1.0, & v );
    }

    //! dtor
    virtual ~TVector () {}

    ////////////////////////////////////////////////
    //! \{
    //!
    //! \name index set functionality
    //!

    //! return first index (offset)
    idx_t           ofs      () const { return _ofs; }

    //! set first index (offset)
    virtual void    set_ofs  ( const idx_t  n ) { _ofs = n; }
    
    //! return size of vector
    virtual size_t  size     () const = 0;

    //! return index set
    TIndexSet       is       () const { return TIndexSet( ofs(), idx_t(ofs() + size()) - 1 ); }

    //!\}

    ///////////////////////////////////////////
    //! \{
    //!
    //! \name management of field type
    //!

    //! return value type of vector
    value_type_t  value_type  () const { return ( is_complex() ? complex_valued : real_valued ); }
    
    //! return true if vector is complex valued
    bool          is_complex  () const { return _complex; }

    //! change between real and complex valued representation
    void set_complex ( const bool b )
    {
        if ( b == _complex )
            return;

        if ( b ) to_complex();
        else     to_real();
        
        _complex = b;
    }

    //! \}
    
protected:

    //! @cond INCLUDE_PROT_PRIV
    
    //! convert data to real representation if possible
    virtual void to_real    () = 0;

    //! convert data to complex representation if possible
    virtual void to_complex () = 0;

    //! @endcond
    
public:
    ////////////////////////////////////////////////
    //! \{
    //!
    //! \name access entries
    //!

    //! return i'th entry
    virtual real          entry  ( const idx_t  i ) const;

    //! return i'th entry
    virtual const complex centry ( const idx_t  i ) const;

    //! set \a i'th entry
    virtual void set_entry  ( const idx_t  i, const real     f );

    //! set \a i'th entry
    virtual void set_centry ( const idx_t  i, const complex  f );

    //! \}

    ////////////////////////////////////////////////
    //! \{
    //!
    //! \name BLAS functionality (real valued)
    //!

    //! fill with constant
    virtual void fill ( const real )
    { HERROR( ERR_NOT_IMPL, "(TVector) fill", "" ); }

    //! fill with random numbers
    virtual void fill_rand ( const uint )
    { HERROR( ERR_NOT_IMPL, "(TVector) fill_rand", "" ); }

    //! scale vector by constant factor
    virtual void scale ( const real )
    { HERROR( ERR_NOT_IMPL, "(TVector) scale", "" ); }

    //! this ≔ f · vector
    virtual void assign ( const real, const TVector * )
    { HERROR( ERR_NOT_IMPL, "(TVector) assign", "" ); }

    //! copy operator for all vectors
    TVector &  operator = ( const TVector & v )
    {
        assign( real(1), & v );
        return *this;
    }
    
    //! return euclidean norm
    virtual real norm2 () const { return Math::sqrt( dot( this ).real() ); }

    //! return infimum norm
    virtual real norm_inf () const { HERROR( ERR_NOT_IMPL, "(TVector) norm_inf", "" ); }
    
    //! this ≔ this + α·x
    virtual void axpy ( const real, const TVector * )
    { HERROR( ERR_NOT_IMPL, "(TVector) axpy", "" ); }
    
    //! \}
    
    ////////////////////////////////////////////////
    //! \{
    //!
    //! \name BLAS functionality (complex valued)
    //!

    //! conjugate entries
    virtual void conjugate ()
    { HERROR( ERR_NOT_IMPL, "(TVector) conjugate", "" ); }
        
    //! fill with constant
    virtual void cfill ( const complex & )
    { HERROR( ERR_NOT_IMPL, "(TVector) cfill", "" ); }

    //! scale vector by constant factor
    virtual void cscale ( const complex & )
    { HERROR( ERR_NOT_IMPL, "(TVector) cscale", "" ); }

    //! this ≔ f · vector
    virtual void cassign ( const complex &, const TVector * )
    { HERROR( ERR_NOT_IMPL, "(TVector) cassign", "" ); }

    //! return dot-product, \f$<x,y> = x^H · y\f$, where \f$x\f$ = this
    virtual complex dot  ( const TVector * ) const
    { HERROR( ERR_NOT_IMPL, "(TVector) dot", "" ); }

    //! return dot-product, \f$<x,y> = x^T · y\f$, where \f$x\f$ = this
    virtual complex dotu ( const TVector * ) const
    { HERROR( ERR_NOT_IMPL, "(TVector) dotu", "" ); }

    //! this ≔ this + α·x
    virtual void caxpy ( const complex &, const TVector * )
    { HERROR( ERR_NOT_IMPL, "(TVector) caxpy", "" ); }
    
    //! \}
    
    ////////////////////////////////////////////////
    //! \{
    //!
    //! \name misc.
    //!

    //
    // memory consumption
    //
    
    //! return size in bytes used by this object
    virtual size_t  byte_size  () const { return sizeof(_ofs) + sizeof(_complex); }

    //! return size in bytes used by this distributed object,
    //! i.e. of all distributed sub vectors
    virtual size_t  global_byte_size () const;

    //
    // virtual constructor
    //

    //! return copy of vector
    virtual auto  copy    () const -> std::unique_ptr< TVector > = 0;
    
    //! assign local values to vector \a x
    virtual void  copy_to ( TVector * x ) const
    {
        if ( x == nullptr )
            HERROR( ERR_ARG, "(TVector) copy_to", "vector is NULL" );
        x->assign( 1.0, this );
    }
    
    //! return object of same class
    virtual auto  create  () const -> std::unique_ptr< TVector > = 0;

    //
    // restriction
    //

    //! create vector restricted to real part of coefficients
    virtual auto restrict_re  () const -> std::unique_ptr< TVector >;

    //! create vector restricted to imaginary part of coefficients
    virtual auto restrict_im  () const -> std::unique_ptr< TVector >;
    
    //! create vector restricted to absolute value of coefficients
    virtual auto restrict_abs () const -> std::unique_ptr< TVector >;
    
    //
    // serialisation
    //

    //! read vector data from byte stream
    virtual void read  ( TByteStream & s );

    //! write vector data to byte stream
    virtual void write ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;

    //
    // parallel methods
    //

    //! sum up nparts parallel copies
    //! (if bs != NULL it will be used)
    virtual void sum ( const TProcSet & p,
                       const uint       pid,
                       const uint       nparts,
                       TByteStream *    bs = NULL );

    //! same as \see sum but sum up between all processors in \a p
    virtual void sum ( const TProcSet & p );
    
    //
    // output
    //

    //! print vector to stdout
    virtual void print ( const uint ofs = 0 ) const;

    HLIB_RTTI_BASE( TVector )
    
    //! \}
};

//////////////////////////////////////////////////////////
//
// wrappers for vector functions
//

//! return dot product <x,y> = x^H · y
inline complex dot      ( const TVector * x, const TVector * y ) { return x->dot( y );   }

//! return dot product <x,y> = x^T · y
inline complex dotu     ( const TVector * x, const TVector * y ) { return x->dotu( y );  }

//! return euclidean norm of \a x
inline real    norm_2   ( const TVector * x )                    { return x->norm2();    }

//! return infimum norm of \a x
inline real    norm_inf ( const TVector * x )                    { return x->norm_inf(); }

//////////////////////////////////////////////////////////
//
// template wrappers
//

//
// scale
//
template <typename T> void
scale ( const T alpha, TVector * x );

template <> inline void
scale<real>    ( const real     alpha, TVector * x )
{ x->scale( alpha ); }

template <> inline void
scale<complex> ( const complex  alpha, TVector * x )
{ x->cscale( alpha ); }
    

//
// axpy
//
template <typename T> void
axpy ( const T alpha, const TVector * x, TVector * y );

template <> inline void
axpy<real>    ( const real     alpha, const TVector * x, TVector * y )
{ y->axpy( alpha, x ); }

template <> inline void
axpy<complex> ( const complex  alpha, const TVector * x, TVector * y )
{ y->caxpy( alpha, x ); }


//
// dot product
//
template <typename T>  T
tdot ( const TVector * x, const TVector * y );

template <> inline real
tdot<real>    ( const TVector * x, const TVector * y )
{ return std::real( dot( x, y ) ); }

template <> inline complex
tdot<complex> ( const TVector * x, const TVector * y )
{ return dot( x, y ); }

template <typename T>  T
tdotu ( const TVector * x, const TVector * y );

template <> inline real
tdotu<real>    ( const TVector * x, const TVector * y )
{ return std::real( dotu( x, y ) ); }

template <> inline complex
tdotu<complex> ( const TVector * x, const TVector * y )
{ return dotu( x, y ); }

//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//
// write vector to file
//
void write ( const TVector *      v,
             const std::string &  filename,
             const std::string &  vecname );

}// namespace DBG

}// namespace HLIB

#endif  // __HLIB_TVECTOR_HH
