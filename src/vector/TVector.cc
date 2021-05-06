//
// Project     : HLib
// File        : TVector.cc
// Description : baseclass for all vector-classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <iostream>

#include "hpro/base/error.hh"
#include "hpro/parallel/NET.hh"
#include "hpro/io/TVectorIO.hh"

#include "hpro/vector/TVector.hh"

namespace HLIB
{

using std::unique_ptr;
using std::make_unique;

////////////////////////////////////////////////
//
// access entries
//

real
TVector::entry  ( const idx_t ) const
{
    HERROR( ERR_NOT_IMPL, "(TVector) entry", "" );
}

const complex
TVector::centry ( const idx_t ) const
{
    HERROR( ERR_NOT_IMPL, "(TVector) entry", "" );
}

void
TVector::set_entry  ( const idx_t, const real )
{
    HERROR( ERR_NOT_IMPL, "(TVector) set_entry", "" );
}

void
TVector::set_centry ( const idx_t, const complex  )
{
    HERROR( ERR_NOT_IMPL, "(TVector) set_entry", "" );
}

////////////////////////////////////////////////
//
// memory consumption
//

//
// return size in bytes used by this distributed object,
// i.e. of all distributed sub matrices
//
size_t
TVector::global_byte_size () const
{
    return byte_size();
}

////////////////////////////////////////////////
//
// restriction
//

//
// create vector restricted to real/imaginary part of coefficients
//
auto
TVector::restrict_re () const -> std::unique_ptr< TVector >
{
    HERROR( ERR_NOT_IMPL, "(TVector) restrict_re", "" );
}

auto
TVector::restrict_im () const -> std::unique_ptr< TVector >
{
    HERROR( ERR_NOT_IMPL, "(TVector) restrict_im", "" );
}
    
auto
TVector::restrict_abs () const -> std::unique_ptr< TVector >
{
    HERROR( ERR_NOT_IMPL, "(TVector) restrict_abs", "" );
}
    
////////////////////////////////////////////////
//
// serialisation
//

void
TVector::read  ( TByteStream & s )
{
    typeid_t  t;

    TStreamable::read( s );
    
    s.get( t );

    if ( t != type() )
        HERROR( ERR_BS_TYPE, "(TVector) read", "found " + RTTI::id_to_type( t ) + ","
               + ", expected " + typestr() );
    
    s.get( _ofs );
    s.get( _complex );
}

void
TVector::write ( TByteStream & s ) const
{
    typeid_t  t = type();

    TStreamable::write( s );
    
    s.put( t );
    s.put( _ofs );
    s.put( _complex );
}

//
// returns size of object in bytestream
//
size_t
TVector::bs_size () const
{
    return TStreamable::bs_size() + sizeof(typeid_t) + sizeof(_ofs) + sizeof(bool);
}

//
// parallel methods
//

//
// sum up nparts parallel copies
// (if bs != nullptr it will be used)
//
void
TVector::sum ( const TProcSet & ps,
               const uint       pid,
               const uint       parts,
               TByteStream    * bs )
{
    unique_ptr< TVector >      v;
    unique_ptr< TByteStream >  tbs;
    uint                       nparts = parts;
    
    if ( bs == nullptr )
    {
        tbs = make_unique< TByteStream >();
        bs  = tbs.get();
    }// if

    while ( nparts > 1 )
    {
        std::vector< TProcSet >  psets;
        int                      idx = -1;

        // create processor sets for 
        ps.split( nparts, psets );

        for ( uint i = 0; i < nparts; i++ )
            if ( pid == psets[i].master() )
            {
                idx = i;
                break;
            }// if

        nparts /= 2;
                
        // if not involved in communication, continue
        if ( idx == -1 )
            continue;

        // send data at all "odd" processors
        if ( idx % 2 == 1 )
        {
            bs->set_size( bs_size() );
            
            write( * bs );
            
            NET::dsend( psets[idx-1].master(), bs->data(), bs->size() );
        }// if

        // and receive it at all "even" procs.
        if ( idx % 2 == 0 )
        {
            bs->set_size( NET::dprobe( psets[idx+1].master() ) );
            
            NET::drecv( psets[idx+1].master(), bs->data(), bs->size() );

            if ( v.get() == nullptr )
                v = unique_ptr< TVector >( create() );
            
            v->read( * bs );

            axpy( 1.0, v.get() );
        }// if
    }// while
}

//
// sum up all vectors in parallel machine defined by <procs>
//
void
TVector::sum ( const TProcSet & procs )
{
    unique_ptr< TVector >  v;
    uint                   nparts = procs.size();
    TByteStream            bs;
    const uint             pid    = NET::pid();

    LOG::print( "sum()" );
    
    while ( nparts > 1 )
    {
        std::vector< TProcSet >  psets;
        int                      idx = -1;

        // create processor sets for 
        procs.split( nparts, psets );

        for ( uint i = 0; i < nparts; i++ )
            if ( pid == psets[i].master() )
            {
                idx = i;
                break;
            }// if

        nparts /= 2;
                
        // if not involved in communication, continue
        if ( idx == -1 )
            continue;

        if ( idx % 2 == 1 )
        {
            //
            // send data at all "odd" processors
            //
            
            bs.set_size( bs_size() );
            
            write( bs );
            
            NET::dsend( psets[idx-1].master(), bs.data(), bs.size() );
        }// if
        else // if ( idx % 2 == 0 )
        {
            //
            // and receive it at all "even" procs.
            //
            
            bs.set_size( NET::dprobe( psets[idx+1].master() ) );
            
            NET::drecv( psets[idx+1].master(), bs.data(), bs.size() );

            if ( v.get() == nullptr )
                v = unique_ptr< TVector >( create() );
            
            v->read( bs );

            axpy( 1.0, v.get() );
        }// if
    }// while

    //
    // finally, scatter result back from master to all other
    // nodes in local processor group
    //

    LOG::print( "scatter" );
    
    // scatter( procs );
}

//
// stream output
//

void
TVector::print ( const uint offset ) const
{
    for ( uint i = 0; i < offset; i++ )
        std::cout << ' ';

    std::cout << typestr() << std::endl;
}

//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//
// write vector to file
//
void
write ( const TVector *      v,
        const std::string &  filename,
        const std::string &  vecname )
{
    TMatlabVectorIO  vio;

    vio.write( v, filename, vecname );
}

}// namespace

}// namespace
