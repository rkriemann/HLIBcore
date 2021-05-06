//
// Project     : HLib
// File        : TScalarVector.cc
// Description : class for a vector of scalar type
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <iostream>

#include "hpro/blas/Algebra.hh"
#include "hpro/parallel/NET.hh"

#include "TRNG.hh"

#include "hpro/vector/TScalarVector.hh"

namespace HLIB
{

// namespace abbr.
namespace B = BLAS;

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//
// TScalarVector
//
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////
//
// misc. methods
//

//
// set size of vector
//
void
TScalarVector::set_size ( const size_t  n )
{
    if ( n != _size )
    {
        if ( is_complex() )
        {
            if ( n == 0 )
            {
                _cvec = B::Vector< complex >();
            }// if
            else
            {
                _cvec = B::Vector< complex >( n );
            }// else
        }// if
        else
        {
            if ( n == 0 )
            {
                _rvec = B::Vector< real >();
            }// if
            else
            {
                _rvec = B::Vector< real >( n );
            }// else
        }// else

        _size = n;

        init_chunk_mutices();        
    }// if
}

//
// switch between real and complex values
//
void
TScalarVector::to_real ()
{
    if ( ! is_complex() || ( _size == 0 ))
        return;

    for ( idx_t  i = 0; i < idx_t(_size); i++ )
    {
        if ( std::imag( _cvec(i) ) != 0.0 )
            HERROR( ERR_COMPLEX, "(TScalarVector) to_real", "vector has imaginary part" );
    }// for

    _rvec = B::Vector< real >( _size );

    for ( idx_t  i = 0; i < idx_t(_size); i++ )
        _rvec(i) = std::real( _cvec(i) );

    _cvec = B::Vector< complex >();
}

void
TScalarVector::to_complex ()
{
    if ( is_complex() || ( _size == 0 ))
        return;

    _cvec = B::Vector< complex >( _size );

    for ( idx_t  i = 0; i < idx_t(_size); i++ )
        _cvec(i) = _rvec(i);

    _rvec = B::Vector< real >();
}
    
//
// return size in bytes used by this object
//
size_t
TScalarVector::byte_size () const
{
    return TVector::byte_size() + sizeof(size_t) + sizeof(B::Vector<real>) + sizeof(B::Vector<complex>) +
        (is_complex() ? sizeof(complex) * size() : sizeof(real) * size());
}

//////////////////////////////////////////////////
//
// BLAS-routines (real valued)
//

//
// fill with constant
//
void
TScalarVector::fill ( const real f )
{
    if ( is_complex() ) B::fill( complex(f), _cvec );
    else                B::fill( f, _rvec );
}

//
// fill with random numbers (but prevent zero entries)
//
void
TScalarVector::fill_rand ( const uint seed )
{
    TRNG  rng( seed );

    if ( is_complex() )
    {
        for ( idx_t  i = 0; i < idx_t(_size); i++ )
        {
            complex  f;
            
            do
            {
                f = complex( real(rng( 1000.0 ) - 500.0),
                             real(rng( 1000.0 ) - 500.0) );
            } while ( f == complex(0) );

            _cvec(i) = f;
        }// for
    }// if
    else
    {
        for ( idx_t  i = 0; i < idx_t(_size); i++ )
        {
            real  f;
            
            do
            {
                f = real(rng( 1000.0 ) - 500.0);
            } while ( f == real(0) );

            _rvec(i) = f;
        }// for
    }// else
}

//
// this = a * this
//
void
TScalarVector::scale ( const real f )
{
    if ( f == real(0) )
        fill( real(0) );
    else
    {
        if ( is_complex() ) B::scale( complex(f), _cvec );
        else                B::scale( f, _rvec );
    }// else
}

//
// this = x
//
void
TScalarVector::assign ( const real f, const TVector * x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) assign", "vector is nullptr" );
    
    if ( IS_TYPE( x, TScalarVector ) )
    {
        const TScalarVector *  sx = cptrcast( x, TScalarVector );
        
        set_complex( x->is_complex() );
        set_is( sx->is() );

        if ( is_complex() )
        {
            B::copy( sx->_cvec, _cvec );

            if ( f != complex(1) )
                scale( f );
        }// if
        else
        {
            B::copy( sx->_rvec, _rvec );

            if ( f != real(1) )
                scale( f );
        }// else
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TScalarVector) assign", x->typestr() );
}

//
// return euclidean norm
//
real
TScalarVector::norm2 () const
{
    if ( is_complex() ) return real(B::norm2( _cvec ));
    else                return real(B::norm2( _rvec ));
}

//
// return infimum norm
//
real
TScalarVector::norm_inf () const
{
    if ( is_complex() )
        return Math::abs( _cvec( B::max_idx( _cvec ) ) );
    else
        return Math::abs( _rvec( B::max_idx( _rvec ) ) );
}   

//
// this = this + a * x
//
void
TScalarVector::axpy ( const real f, const TVector * x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) axpy", "vector is nullptr" );
    
    if ( ! IS_TYPE( x, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TScalarVector) axpy", x->typestr() );

    const TScalarVector * sx = cptrcast( x, TScalarVector );

    if ( size() != sx->size() )
        HERROR( ERR_VEC_SIZE, "(TScalarVector) axpy", "x has wrong size" );
    
    if ( is_complex() != x->is_complex() )
        HERROR( ERR_COMPLEX, "(TScalarVector) axpy", "can not mix real and complex" );

    if ( x->is_complex() ) B::add( complex(f), sx->_cvec, _cvec );
    else                   B::add( f, sx->_rvec, _rvec );
}

//
// this = this + x (multi-thread safe version, x may be sub vector)
//
void
TScalarVector::add_sub_mt ( const TScalarVector &  x )
{
    HERROR( ERR_NOT_IMPL, "", "" );
    
    if ( ! x.is().is_subset_of( is() ) )
        HERROR( ERR_INDEXSET, "(TScalarVector) add_sub_mt", "given vector is not a sub-vector" );
    
    if ( is_complex() != x.is_complex() )
        HERROR( ERR_COMPLEX, "(TScalarVector) add_sub_mt", "can not mix real and complex" );

    const TIndexSet  is_x( x.is() );
    const TIndexSet  is_y( is() );
    const idx_t      ofs_loc_glo = is_x.first() - is_y.first();
    idx_t            start_idx   = ofs_loc_glo;
    idx_t            chunk       = start_idx / SCALAR_CHUNK_SIZE;
    const idx_t      last_idx    = is_x.last() - is_y.first();
    idx_t            end_idx     = std::min< idx_t >( (chunk+1) * SCALAR_CHUNK_SIZE - 1, last_idx );
    
    if ( is_complex() )
    {
        while ( start_idx <= end_idx )
        {
            const B::Range        is_chunk( start_idx, end_idx );
            B::Vector< complex >  y_chunk( blas_cvec(), is_chunk );
            B::Vector< complex >  x_chunk( x.blas_cvec(), is_chunk - ofs_loc_glo );

            {
                // lock_t  lock( _mutices[ chunk ] );
                
                B::add( complex(1), x_chunk, y_chunk );
            }

            ++chunk;
            start_idx = end_idx + 1;
            end_idx   = std::min< idx_t >( end_idx + SCALAR_CHUNK_SIZE, last_idx );
        }// while
    }// if
    else
    {
        while ( start_idx <= end_idx )
        {
            const B::Range     is_chunk( start_idx, end_idx );
            B::Vector< real >  y_chunk( blas_rvec(), is_chunk );
            B::Vector< real >  x_chunk( x.blas_rvec(), is_chunk - ofs_loc_glo );

            {
                // lock_t  lock( _mutices[ chunk ] );
                
                B::add( real(1), x_chunk, y_chunk );
            }

            ++chunk;
            start_idx = end_idx + 1;
            end_idx   = std::min< idx_t >( end_idx + SCALAR_CHUNK_SIZE, last_idx );
        }// while
    }// else
}

//////////////////////////////////////////////////
//
// BLAS-routines (complex valued)
//

//
// fill with constant
//
void
TScalarVector::cfill ( const complex & f )
{
    if ( std::imag(f) != real(0) )
        set_complex( true );
    
    if ( is_complex() ) B::fill( f, _cvec );
    else                B::fill( std::real(f), _rvec );
}

//
// this = a * this
//
void
TScalarVector::cscale ( const complex & f )
{
    if ( f == complex(0) )
        fill( real(0) );
    else
    {
        if ( is_complex() )
            B::scale( f, _cvec );
        else
        {
            if ( f.imag() == real(0) ) 
                B::scale( std::real(f), _rvec );
            else
            {
                set_complex( true );
                B::scale( f, _cvec );
            }// else
        }// else
    }// else
}

//
// this = f * x
//
void
TScalarVector::cassign ( const complex & f, const TVector * x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) cassign", "vector is nullptr" );
    
    if ( IS_TYPE( x, TScalarVector ) )
    {
        const TScalarVector *  sx = cptrcast( x, TScalarVector );
        
        set_complex( x->is_complex() || ( f.imag() != real(0) ));
        set_is( sx->is() );

        if ( x->is_complex() )
        {
            B::copy( sx->_cvec, _cvec );
            
            if ( f != complex(1) )
                cscale( f );
        }// if
        else
        {
            if ( is_complex() )
                HERROR( ERR_COMPLEX, "(TScalarVector) cassign", "can not mix real and complex" );
            else
            {
                B::copy( sx->_rvec, _rvec );

                if ( f != real(1) )
                    scale( std::real(f) );
            }// else
        }// else
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TScalarVector) cassign", x->typestr() );
}

//
// return inner product ( this^H * x )
//
complex
TScalarVector::dot ( const TVector * x ) const
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) dot", "vector is nullptr" );
    
    if ( IS_TYPE( x, TScalarVector ) )
    {
        if ( x->is_complex() )
        {
            if ( is_complex() )
                return B::dot( _cvec, cptrcast( x, TScalarVector )->_cvec );
            else
                HERROR( ERR_COMPLEX, "(TScalarVector) dot", "can not mix real and complex" );
        }// if
        else
        {
            if ( is_complex() )
                HERROR( ERR_COMPLEX, "(TScalarVector) dot", "can not mix real and complex" );
            else
                return real(B::dot( _rvec, cptrcast( x, TScalarVector )->_rvec ));
        }// else
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TScalarVector) dot", x->typestr() );
    
    return complex(0);
}

//
// return inner product ( this^T * x )
//
complex
TScalarVector::dotu ( const TVector * x ) const
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) dotu", "vector is nullptr" );
    
    if ( IS_TYPE( x, TScalarVector ) )
    {
        if ( x->is_complex() )
        {
            if ( is_complex() )
                return B::dotu( _cvec, cptrcast( x, TScalarVector )->_cvec );
            else
                HERROR( ERR_COMPLEX, "(TScalarVector) dot", "can not mix real and complex" );
        }// if
        else
        {
            if ( is_complex() )
                HERROR( ERR_COMPLEX, "(TScalarVector) dot", "can not mix real and complex" );
            else
                return real(B::dotu( _rvec, cptrcast( x, TScalarVector )->_rvec ));
        }// else
    }// if
    else
        HERROR( ERR_VEC_TYPE, "(TScalarVector) dotu", x->typestr() );
    
    return complex(0);
}

//
// this = this + a * x
//
void
TScalarVector::caxpy ( const complex & f, const TVector * x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) caxpy", "vector is nullptr" );
    
    if ( ! IS_TYPE( x, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TScalarVector) caxpy", x->typestr() );

    const TScalarVector * sx = cptrcast( x, TScalarVector );

    if ( size() != sx->size() )
        HERROR( ERR_VEC_SIZE, "(TScalarVector) caxpy", "x has wrong size" );
    
    if ( x->is_complex() || ( f.imag() != real(0) ))
        set_complex( true );
    
    if ( x->is_complex() )
        B::add( f, sx->_cvec, _cvec );
    else
    {
        if ( is_complex() )
            HERROR( ERR_COMPLEX, "(TScalarVector) cassign", "can not mix real and complex" );
        else
            B::add( std::real(f), sx->_rvec, _rvec );
    }// else
}

//////////////////////////////////////////////////
//
// misc. methods
//

//
// create vector restricted to real/imaginary part of coefficients
//
auto
TScalarVector::restrict_re () const -> std::unique_ptr< TVector >
{
    auto  v = create();
    auto  t = ptrcast( v.get(), TScalarVector );

    t->set_is( is() );

    if ( is_complex() )
    {
        for ( idx_t  i = 0; i < idx_t(_size); i++ )
            t->_rvec(i) = std::real( _cvec(i) );
    }// if
    else
    {
        for ( idx_t  i = 0; i < idx_t(_size); i++ )
            t->_rvec(i) = _rvec(i);
    }// else

    return v;
}

auto
TScalarVector::restrict_im () const -> std::unique_ptr< TVector >
{
    auto  v = create();
    auto  t = ptrcast( v.get(), TScalarVector );

    t->set_is( is() );

    if ( is_complex() )
    {
        for ( idx_t  i = 0; i < idx_t(_size); i++ )
            t->_rvec(i) = imag( _cvec(i) );
    }// if
    else
    {
        for ( idx_t  i = 0; i < idx_t(_size); i++ )
            t->_rvec(i) = real(0);
    }// else

    return v;
}

auto
TScalarVector::restrict_abs () const -> std::unique_ptr< TVector >
{
    auto  v = create();
    auto  t = ptrcast( v.get(), TScalarVector );

    t->set_is( is() );

    if ( is_complex() )
    {
        for ( idx_t  i = 0; i < idx_t(_size); i++ )
            t->_rvec(i) = abs( _cvec(i) );
    }// if
    else
    {
        for ( idx_t  i = 0; i < idx_t(_size); i++ )
            t->_rvec(i) = real(0);
    }// else

    return v;
}

//
// copy from C array \a v
//
void 
TScalarVector::copy_from ( const real *  v )
{
    if ( v == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) copy_from", "array is nullptr" );

    if ( is_complex() )
    { 
        for ( idx_t  i = 0; i < idx_t(_cvec.length()); ++i )
        {
            _cvec(i) = v[i];
        }// for
    }// if
    else
    {
        for ( idx_t  i = 0; i < idx_t(_rvec.length()); ++i )
        {
            _rvec(i) = v[i];
        }// for
    }// else
}

//
// copy to C array \a v
//
void 
TScalarVector::copy_to ( real *  v )
{
    if ( v == nullptr )
        HERROR( ERR_ARG, "(TScalarVector) copy_to", "array is nullptr" );

    if ( is_complex() )
    { 
        HERROR( ERR_REAL_CMPLX, "(TScalarVector) copy_to", "vector is complex valued" );
    }// if
    else
    {
        for ( idx_t  i = 0; i < idx_t(_rvec.length()); ++i )
        {
            v[i] = _rvec(i);
        }// for
    }// else
}

//
// permute entries according to \a perm
//
void
TScalarVector::permute ( const TPermutation &  perm )
{
    perm.permute( this );
}

//
// stream output
//
void
TScalarVector::print ( const uint offset ) const
{
    for ( uint i = 0; i < offset; i++ )
        std::cout << ' ';

    std::cout << typestr() << " ( ";

    if ( is_complex() )
    {
        for ( idx_t  i = 0; i < idx_t(_size); i++ )
        {
            if ( i > 0 ) std::cout << ", ";
            std::cout << this->centry(i);
        }// for
    }// if
    else
    {
        for ( idx_t  i = 0; i < idx_t(_size); i++ )
        {
            if ( i > 0 ) std::cout << ", ";
            std::cout << this->entry(i);
        }// for
    }// else

    std::cout << " )" << std::endl;
}

//
// serialisation
//

void
TScalarVector::read  ( TByteStream & s )
{
    TVector::read( s );

    size_t  vec_size;
    
    s.get( vec_size );

    set_size( vec_size );

    if ( is_complex() )
        s.get( blas_cvec().data(), sizeof(complex) * _size );
    else
        s.get( blas_rvec().data(), sizeof(real) * _size );
}

void
TScalarVector::write ( TByteStream & s ) const
{
    TVector::write( s );
    
    s.put( _size );
    
    if ( is_complex() )
        s.put( blas_cvec().data(), sizeof(complex) * _size );
    else
        s.put( blas_rvec().data(), sizeof(real) * _size );
}

//
// returns size of object in bytestream
//
size_t
TScalarVector::bs_size () const
{
    return TVector::bs_size() + sizeof(_size) +
        (is_complex() ? (sizeof(complex) * _size) : (sizeof(real) * _size));
}

//
// pointwise summation between all vectors in \a ps
//
void
TScalarVector::sum ( const TProcSet & ps )
{
    if ( is_complex() )
    {
        B::Vector< complex >  t( _cvec.length() );

        NET::reduce_all( ps, _cvec.data(), t.data(), _cvec.length(), NET::OP_SUM );
        B::copy( t, _cvec );
    }// if
    else
    {
        B::Vector< real >  t( _rvec.length() );

        NET::reduce_all( ps, _rvec.data(), t.data(), _rvec.length(), NET::OP_SUM );
        B::copy( t, _rvec );
    }// else
}

}// namespace
