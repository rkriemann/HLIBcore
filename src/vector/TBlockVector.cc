//
// Project     : HLib
// File        : TBlockVector.cc
// Description : class for a vector build out of other vectors
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <iostream>

#include "hpro/vector/TScalarVector.hh"

#include "hpro/vector/TBlockVector.hh"

namespace HLIB
{

// namespace abbr.
namespace B = BLAS;

///////////////////////////////////////////////
//
// constructor and destructor
//

TBlockVector::~TBlockVector ()
{
    for ( uint i = 0; i < n_blocks(); i++ )
        delete block(i);
}

///////////////////////////////////////////////
//
// access internal data
//

//
// return size of vector
//
size_t
TBlockVector::size () const
{
    size_t  s = 0;
    
    for ( uint i = 0; i < n_blocks(); i++ )
    {
        if ( block(i) != nullptr )
            s += block(i)->size();
    }// for

    return s;
}

//
// set single vector block
//
void
TBlockVector::set_block ( const uint i, TVector * v )
{
    if ( i >= n_blocks() )
        HERROR( ERR_VEC_STRUCT, "(TBlockVector) set_block", "index too high" );

    if ( block(i) != nullptr )
        delete block(i);

    _blocks[i] = v;
}
    
//
// setup blockstructure of vector
//
void
TBlockVector::set_block_struct ( const uint n )
{
    uint  old_n = n_blocks();

    for ( uint i = n; i < old_n; i++ )
        if ( _blocks[i] != nullptr )
            delete _blocks[i];
    
    _blocks.resize( n, nullptr );
}

//
// switch between complex and real format
//
void
TBlockVector::to_real ()
{
    if ( ! is_complex() )
        return;

    for ( uint i = 0; i < n_blocks(); i++ )
        if ( _blocks[i] != nullptr )
            _blocks[i]->set_complex( false );
}

void
TBlockVector::to_complex ()
{
    if ( is_complex() )
        return;


    for ( uint i = 0; i < n_blocks(); i++ )
        if ( _blocks[i] != nullptr )
            _blocks[i]->set_complex( true );
}
    
////////////////////////////////////////////////
//
// access entries
//

//
// return i'th entry
//
real
TBlockVector::entry  ( const idx_t  row ) const
{
    const idx_t  vofs = ofs();
    
    for ( uint i = 0; i < n_blocks(); i++ )
    {
        const TVector * v_i = block(i);
        
        if ( v_i == nullptr )
            continue;
        
        if (( v_i->ofs() - vofs <= row ) && ( row < v_i->ofs() - vofs + idx_t( v_i->size() ) ))
            return v_i->entry( row - v_i->ofs() + vofs );
    }// for

    return 0.0;
}

const complex
TBlockVector::centry ( const idx_t  row ) const
{
    const idx_t  vofs = ofs();
    
    for ( uint i = 0; i < n_blocks(); i++ )
    {
        const TVector * v_i = block(i);
        
        if ( v_i == nullptr )
            continue;
        
        if (( v_i->ofs() - vofs <= row ) && ( row < v_i->ofs() - vofs + idx_t( v_i->size() ) ))
            return v_i->centry( row - v_i->ofs() + vofs );
    }// for

    return 0.0;
}

//
// set i'th entry
//
void
TBlockVector::set_entry  ( const idx_t  row,
                           const real   f )
{
    const idx_t  vofs = ofs();
    
    for ( uint i = 0; i < n_blocks(); i++ )
    {
        TVector * v_i = block(i);
        
        if ( v_i == nullptr )
            continue;
        
        if (( v_i->ofs() - vofs <= row ) && ( row < v_i->ofs() - vofs + idx_t( v_i->size() ) ))
        {
            v_i->set_entry( row - v_i->ofs() + vofs, f );
            return;
        }// if
    }// for
}

void
TBlockVector::set_centry ( const idx_t    row,
                           const complex  f )
{
    const idx_t  vofs = ofs();
    
    for ( uint i = 0; i < n_blocks(); i++ )
    {
        TVector * v_i = block(i);
        
        if ( v_i == nullptr )
            continue;
        
        if (( v_i->ofs() - vofs <= row ) && ( row < v_i->ofs() - vofs + idx_t( v_i->size() ) ))
        {
            v_i->set_centry( row - v_i->ofs() + vofs, f );
            return;
        }// if
    }// for
}

//////////////////////////////////////////////////
//
// BLAS-routines (real valued)
//

//
// fill with constant
//
void
TBlockVector::fill ( const real f )
{
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->fill( f );
}
    
//
// fill with random numbers
//
void
TBlockVector::fill_rand ( const uint seed )
{
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->fill_rand( seed );
}

//
// this = this + a * x
//
void
TBlockVector::axpy ( const real f, const TVector * x )
{
    if ( f == 0.0 )
        return;
    
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TBlockVector) axpy", "given vector is nullptr" );

    if ( is_complex() != x->is_complex() )
        HERROR( ERR_REAL_CMPLX, "(TBlockVector) axpy", "" );
    
    if ( IS_TYPE( x, TBlockVector ) )
    {
        const TBlockVector * b = cptrcast( x, TBlockVector );

        if ( b->n_blocks() != n_blocks() )
            HERROR( ERR_VEC_STRUCT, "(TBlockVector) axpy",
                    "given vector has different number of blocks" );
        
        for ( uint i = 0; i < n_blocks(); i++ )
            if ( block(i) != nullptr )
                block(i)->axpy( f, b->block(i) );
    }// if
    else if ( IS_TYPE( x, TScalarVector ) )
    {
        const TScalarVector * s = cptrcast( x, TScalarVector );

        for ( uint i = 0; i < n_blocks(); i++ )
        {
            TVector *  bi = block(i);
                    
            if ( bi != nullptr )
            {
                if ( is_complex() )
                {
                    B::Vector< complex >  t( s->blas_cvec(), bi->is() - s->ofs() );
                    TScalarVector         st( bi->is(), t );
                                              
                    bi->axpy( f, & st );
                }// if
                else
                {
                    B::Vector< real >  t( s->blas_rvec(), bi->is() - s->ofs() );
                    TScalarVector      st( bi->is(), t );
                                              
                    bi->axpy( f, & st );
                }// else
            }// if
        }// for
    }// if
    else
    {
        HERROR( ERR_VEC_TYPE, "(TBlockVector) axpy", x->typestr() );
    }// else

}

//
// this = a * this
//
void
TBlockVector::scale ( const real f )
{
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->scale( f );
}

//
// this = f * x
//
void
TBlockVector::assign ( const real f, const TVector * x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TBlockVector) assign", "given vector is nullptr" );
    
    if ( ! IS_TYPE( x, TBlockVector ) )
        HERROR( ERR_VEC_TYPE, "(TBlockVector) assign",
               "expected TBlockVector, got " + x->typestr() );

    const TBlockVector * b = cptrcast( x, TBlockVector );

    set_ofs( b->ofs() );
    
    if ( b->n_blocks() != n_blocks() )
        HERROR( ERR_VEC_STRUCT, "(TBlockVector) assign",
               "given vector has different number of blocks" );
    
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->assign( f, b->block(i) );
}

//
// return infimum norm
//
real
TBlockVector::norm_inf () const
{
    real  f = 0.0;

    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            f = std::max( f, block(i)->norm_inf() );

    return f;
}

//////////////////////////////////////////////////
//
// BLAS-routines (complex valued)
//

//
// conjugate entries
//
void
TBlockVector::conjugate ()
{
    for ( uint  i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->conjugate();
}

//
// fill with constant
//
void
TBlockVector::cfill ( const complex & f )
{
    for ( uint  i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->cfill( f );
}
    
//
// this = a * this
//
void
TBlockVector::cscale ( const complex & f )
{
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->cscale( f );
}

//
// this = f * x
//
void
TBlockVector::cassign ( const complex & f, const TVector * x )
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TBlockVector) cassign", "given vector is nullptr" );
    
    if ( ! IS_TYPE( x, TBlockVector ) )
        HERROR( ERR_VEC_TYPE, "(TBlockVector) cassign",
               "expected TBlockVector, got " + x->typestr() );

    const TBlockVector * b = cptrcast( x, TBlockVector );

    set_ofs( b->ofs() );
    
    if ( b->n_blocks() != n_blocks() )
        HERROR( ERR_VEC_STRUCT, "(TBlockVector) cassign",
               "given vector has different number of blocks" );
    
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->cassign( f, b->block(i) );
}

//
// this = this + a * x
//
void
TBlockVector::caxpy ( const complex & f, const TVector * x )
{
    if ( f == complex(0) )
        return;
    
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TBlockVector) caxpy", "given vector is nullptr" );
    
    if ( ! IS_TYPE( x, TBlockVector ) )
        HERROR( ERR_VEC_TYPE, "(TBlockVector) caxpy",
               "expected TBlockVector, got " + x->typestr() );

    const TBlockVector * b = cptrcast( x, TBlockVector );

    if ( b->n_blocks() != n_blocks() )
        HERROR( ERR_VEC_STRUCT, "(TBlockVector) caxpy",
               "given vector has different number of blocks" );
    
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            block(i)->caxpy( f, b->block(i) );
}

//
// return inner product ( this^H * x )
//
complex
TBlockVector::dot ( const TVector * x ) const
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TBlockVector) dot", "given vector is nullptr" );
    
    if ( ! IS_TYPE( x, TBlockVector ) )
        HERROR( ERR_VEC_TYPE, "(TBlockVector) dot",
               "expected TBlockVector, got " + x->typestr() );

    const TBlockVector * b = cptrcast( x, TBlockVector );

    if ( b->n_blocks() != n_blocks() )
        HERROR( ERR_VEC_STRUCT, "(TBlockVector) dot", "given vector has different number of blocks" );

    complex  f = 0.0;
    
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            f += block(i)->dot( b->block(i) );

    return f;
}

//
// return inner product ( this^T * x )
//
complex
TBlockVector::dotu ( const TVector * x ) const
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TBlockVector) dotu", "given vector is nullptr" );
    
    if ( ! IS_TYPE( x, TBlockVector ) )
        HERROR( ERR_VEC_TYPE, "(TBlockVector) dotu",
               "expected TBlockVector, got " + x->typestr() );

    const TBlockVector * b = cptrcast( x, TBlockVector );

    if ( b->n_blocks() != n_blocks() )
        HERROR( ERR_VEC_STRUCT, "(TBlockVector) dotu", "given vector has different number of blocks" );

    complex  f = 0.0;
    
    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            f += block(i)->dotu( b->block(i) );

    return f;
}

///////////////////////////////////////////////
//
// misc.
//

//
// return size in bytes used by this object
//
size_t
TBlockVector::byte_size () const
{
    size_t  s = TVector::byte_size() + sizeof(_blocks);

    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            s += block(i)->byte_size();

    return s;
}

//
// return copy of matrix
//
std::unique_ptr< TVector >
TBlockVector::copy  () const
{
    auto  v = create();
    auto  x = ptrcast( v.get(), TBlockVector );

    x->set_ofs( ofs() );
    x->set_block_struct( n_blocks() );

    for ( uint i = 0; i < n_blocks(); i++ )
        if ( block(i) != nullptr )
            x->set_block( i, block(i)->copy().release() );

    return v;
}
    
//
// stream output
//
void
TBlockVector::print ( const uint offset ) const
{
    for ( uint i = 0; i < offset; i++ )
        std::cout << ' ';

    std::cout << "TBlockVector [ ";
    
    for ( uint i = 0; i < n_blocks(); i++ )
    {
        if ( block(i) != nullptr )
            block(i)->print( offset+2 );
        else
        {
            for ( uint j = 0; j < offset+2; j++ )
                std::cout << ' ';
            
            std::cout << "nullptr";
        }// else

        std::cout << std::endl;
    }// if

    std::cout << " ]" << std::endl;
}

}// namespace
