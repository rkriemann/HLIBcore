//
// Project     : HLib
// File        : TFnSpace.cc
// Description : implements function spaces over grids
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/bem/TFnSpace.hh"

namespace HLIB
{

using std::vector;
using std::unique_ptr;
using std::make_unique;

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// TFnSpace
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
//
// constructor and destructor
//

TFnSpace::TFnSpace ( const TGrid * agrid )
        : _grid(agrid)
{
    if ( agrid == nullptr )
        HERROR( ERR_ARG, "(TFnSpace)", "given grid is nullptr" );
}

TFnSpace::~TFnSpace ()
{
}

///////////////////////////////////////////////////////////////////////
//
// function evaluation
//
///////////////////////////////////////////////////////////////////////

//
// evaluate function <fn> at index positions on grid
// and build corresponding vector
//
template < typename T_val >
unique_ptr< TScalarVector >
TFnSpace::eval ( const TBEMFunction< T_val > * fn ) const
{
    if ( fn == nullptr )
        HERROR( ERR_ARG, "(TFnSpace) eval", "given function is nullptr" );

    //
    // build vector
    //

    auto  v = make_unique< TScalarVector >( this->n_indices(), 0, fn->is_complex() );
    
    for ( uint i = 0; i < this->n_indices(); i++ )
    {
        // get index position
        const T3Point  t = index_coord( i );

        // get normal at index (only valid for constant ansatz)
        const T3Point  n( _grid->tri_normal( _supp_list_ptr[i] ) );

        // evaluate
        if ( fn->is_complex() ) v->set_centry( i, complex( fn->eval( t, n ) ) );
        else                    v->set_entry(  i, std::real( fn->eval( t, n ) ) );
    }// for

    return v;
}
    
////////////////////////////////////////////////////////
//
// misc.
//

//
// return coordinate info for grid
//
unique_ptr< TCoordinate >
TFnSpace::build_coord ( const bool  with_bbox ) const
{
    const size_t       nindices = _indices.size();
    vector< T3Point >  coord( nindices );

    if ( with_bbox )
    {
        vector< TBBox >  bbox( nindices );

        for ( size_t  i = 0; i < nindices; ++i )
        {
            //
            // start with position of index
            //
                
            coord[i] = _indices[i];

            TBBox  bbox_i( coord[i], coord[i] );
                
            //
            // and enhance by support, e.g. coordinates of
            // triangles in support
            //
                
            for ( auto  tri : support( idx_t(i) ) )
            {
                // loop over the three triangle vertices
                for ( uint k = 0; k < 3; k++ )
                {
                    const T3Point  pos( _grid->vertex( _grid->triangle( tri ).vtx[k] ) );
                    const TBBox    pos_bbox( pos, pos );
                        
                    bbox_i.join( pos_bbox );
                }// for
            }// for

            bbox[i] = bbox_i;
        }// for

        return make_unique< TCoordinate >( coord, bbox );
    }// if
    else 
    {
        for ( size_t  i = 0; i < nindices; i++ )
        {
            coord[i] = _indices[i];
        }// for

        return make_unique< TCoordinate >( coord );
    }// if
}

//
// return size of grid in bytes
//
size_t
TFnSpace::byte_size () const
{
    size_t  size = sizeof(TGrid*);
    
    size += sizeof( _supp_list_ptr ) + _supp_list_ptr.size() * sizeof( idx_t );
    size += sizeof( _supp_list ) + _supp_list.size() * sizeof( idx_t );
    size += sizeof( _indices ) + _indices.size() * sizeof( T3Point );

    return size;
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// TConstFnSpace
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
//
// constructor and destructor
//

TConstFnSpace::TConstFnSpace ( const TGrid * agrid )
        : TFnSpace( agrid )
{
    construct();
}

TConstFnSpace::~TConstFnSpace ()
{
}
    
////////////////////////////////////////////////////////
//
// private functions
//

//
// construct function space by building indices
// and their support
//
void
TConstFnSpace::construct ()
{
    //
    // indices are built at the midpoints of each triangle
    //

    const size_t  n_triangles = _grid->n_triangles();

    _indices.resize( n_triangles );

    for ( size_t  i = 0; i < n_triangles; i++ )
    {
        const idx_t  vid[3] = { _grid->triangle( idx_t(i) ).vtx[0],
                                _grid->triangle( idx_t(i) ).vtx[1],
                                _grid->triangle( idx_t(i) ).vtx[2] };
            
        _indices[i]  = T3Point( _grid->vertex( vid[0] ) +
                                _grid->vertex( vid[1] ) +
                                _grid->vertex( vid[2] ) );
        _indices[i] *= 1.0 / 3.0;
    }// for

    //
    // support is defined by the single triangle enclosing index i
    // and since they are numbered in the same order, the mapping
    // is just the identity
    //

    _supp_list_ptr.resize( n_triangles+1, 0 );
    _supp_list.resize( n_triangles, 0 );

    for ( uint i = 0; i < n_triangles; i++ )
    {
        _supp_list_ptr[i] = i;
        _supp_list[i]     = i;
    }// for
        
    _supp_list_ptr[n_triangles] = uint( n_triangles );

    //
    // build mapping of triangle to triangle local indices
    // - one index per triangle
    //

    _tri_idx_ptr.resize( n_triangles+1, 0 );
    _tri_idx.resize( n_triangles+1, 0 );

    for ( uint i = 0; i < n_triangles; i++ )
    {
        _tri_idx_ptr[i] = i;
        _tri_idx[i]     = i;
    }// for
        
    _tri_idx_ptr[n_triangles] = uint( n_triangles );
}

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// TLinFnSpace
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////
//
// constructor and destructor
//

TLinearFnSpace::TLinearFnSpace ( const TGrid * agrid )
        : TFnSpace( agrid )
{
    construct();
}

TLinearFnSpace::~TLinearFnSpace ()
{
}
    
////////////////////////////////////////////////////////
//
// private functions
//

//
// construct function space by building indices
// and their support
//
void
TLinearFnSpace::construct ()
{
    //
    // just copy position of vertices to index-array
    //

    const size_t  n_vertices = _grid->n_vertices();
        
    _indices.resize( n_vertices );

    for ( idx_t  i = 0; i < idx_t( n_vertices ); ++i )
    {
        _indices[i] = _grid->vertex( i );
    }// for
        
    //
    // count number of triangles in support per index and
    // built start-end-array of positions in the support-list array
    //

    const size_t    n_triangles = _grid->n_triangles();
    vector< uint >  count( n_vertices, 0 );
    
    for ( idx_t  i = 0; i < idx_t( n_triangles ); ++i )
    {
        count[ _grid->triangle(i).vtx[0] ]++;
        count[ _grid->triangle(i).vtx[1] ]++;
        count[ _grid->triangle(i).vtx[2] ]++;
    }// for

    idx_t  pos = 0;

    _supp_list_ptr.resize( n_vertices+1, 0 );

    for ( size_t  i = 0; i < n_vertices; ++i )
    {
        _supp_list_ptr[i] = pos;
        pos              += idx_t( count[i] );
    }// for

    _supp_list_ptr[n_vertices] = pos;

    //
    // fill support-list per index with the corresponding triangle indices
    // - pos now holds the total number of references between indices
    //   and triangles
    //

    _supp_list.resize( pos, 0 );
    
    for ( size_t  i = 0; i < n_vertices; ++i )
        count[i] = 0;

    for ( size_t  i = 0; i < n_triangles; ++i )
    {
        for ( uint  j = 0; j < 3; j++ )
        {
            const idx_t  vid = _grid->triangle(  idx_t(i) ).vtx[j];
            
            _supp_list[ _supp_list_ptr[ vid ] + count[vid] ] = idx_t( i );
            count[vid]++;
        }// for
    }// for

    //
    // build mapping of triangle to triangle local indices
    // - three indices per triangle
    //

    _tri_idx_ptr.resize( n_triangles+1, 0 );
    _tri_idx.resize( 3 * n_triangles, 0 );

    pos = 0;
    
    for ( size_t  i = 0; i < n_triangles; ++i )
    {
        _tri_idx_ptr[ i ] = pos;
        
        for ( uint j = 0; j < 3; j++ )
        {
            _tri_idx[ pos ] = _grid->triangle( idx_t(i) ).vtx[j];
            ++pos;
        }// for
    }// for

    _tri_idx_ptr[ n_triangles ] = pos;
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template unique_ptr< TScalarVector >  TFnSpace::eval< real >     ( const TBEMFunction< real > *    ) const;
template unique_ptr< TScalarVector >  TFnSpace::eval< complex >  ( const TBEMFunction< complex > * ) const;

}// namespace
