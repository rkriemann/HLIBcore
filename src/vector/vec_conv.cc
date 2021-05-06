//
// Project     : HLib
// File        : vec_conv.cc
// Description : vector conversion functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "list.hh"

#include "hpro/vector/TScalarVector.hh"
#include "hpro/vector/TBlockVector.hh"

#include "hpro/vector/vec_conv.hh"

namespace HLIB
{

using namespace std;

// namespace abbr.
namespace B = BLAS;

namespace
{

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// local functions
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
//
// conversion to scalar
//
///////////////////////////////////////////////////////////////////////

TVector *
block_to_scalar ( const TBlockVector *  v )
{
    auto                      s = make_unique< TScalarVector >( v->is(), v->is_complex() );
    list< const TVector * >   subvec;

    subvec.push_back( v );
    
    while ( ! subvec.empty() )
    {
        const TVector *  w = behead( subvec );

        if ( IS_TYPE( w, TBlockVector ) )
        {
            const TBlockVector *  bw = cptrcast( w, TBlockVector );

            for ( uint  i = 0; i < bw->n_blocks(); ++i )
                if ( bw->block(i) != nullptr )
                    subvec.push_back( bw->block(i) );
        }// if
        else if ( IS_TYPE( w, TScalarVector ) )
        {
            const TScalarVector *  sw = cptrcast( w, TScalarVector );

            if ( v->is_complex() )
            {
                B::Vector< complex >  s_part( s->blas_cvec(), sw->is() - s->ofs() );

                B::copy( sw->blas_cvec(), s_part );
            }// if
            else
            {
                B::Vector< real >  s_part( s->blas_rvec(), sw->is() - s->ofs() );

                B::copy( sw->blas_rvec(), s_part );
            }// else
        }// if
        else
            HERROR( ERR_VEC_TYPE, "block_to_scalar", w->typestr() );
    }// while

    return s.release();
}

///////////////////////////////////////////////////////////////////////
//
// conversion to blocked
//
///////////////////////////////////////////////////////////////////////

TVector *
scalar_to_blocked ( const TScalarVector *  v,
                    const TCluster *       ct )
{
    unique_ptr< TVector >  w;
        
    if ( ct->is_leaf() )
    {
        if ( v->is_complex() )
        {
            B::Vector< complex >  v_ct( v->blas_cvec(), *ct - v->ofs(), copy_value );

            w = make_unique< TScalarVector >( *ct, v_ct );
        }// if
        else
        {
            B::Vector< real >     v_ct( v->blas_rvec(), *ct - v->ofs(), copy_value );

            w = make_unique< TScalarVector >( *ct, v_ct );
        }// else
    }// if
    else
    {
        vector< TVector * >  subblocks( ct->nsons() );
        
        for ( uint i = 0; i < ct->nsons(); ++i )
        {
            if ( ct->son(i) != nullptr )
                subblocks[i] = scalar_to_blocked( v, ct->son(i) );
        }// for

        w = make_unique< TBlockVector >( *ct, subblocks );
    }// else

    return w.release();
}

}// namespace anonymous

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// global functions
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//
// convert \a v into scalar vector
//
TVector *
to_scalar ( const TVector *  v )
{
    if ( v == nullptr )
        return nullptr;

    unique_ptr< TVector >  res;
    
    if ( IS_TYPE( v, TBlockVector ) )
    {
        res = unique_ptr< TVector >( block_to_scalar( cptrcast( v, TBlockVector ) ) );
    }// if
    else
        HERROR( ERR_VEC_TYPE, "to_scalar", v->typestr() );
        
    return res.release();
}

//
// convert \a v into blocked vector as defined by
// to cluster tree \a ct
//
TVector *
to_blocked ( const TVector *   v,
             const TCluster *  ct )
{
    if ( v == nullptr )
        return nullptr;

    if ( ct == nullptr )
        HERROR( ERR_ARG, "to_blocked", "cluster tree is nullptr" );

    unique_ptr< TVector >  res;
    
    if ( IS_TYPE( v, TScalarVector ) )
        res = unique_ptr< TVector >( scalar_to_blocked( cptrcast( v, TScalarVector ), ct ) );
    else
        HERROR( ERR_VEC_TYPE, "to_blocked", v->typestr() );
        
    return res.release();
}

}// namespace
