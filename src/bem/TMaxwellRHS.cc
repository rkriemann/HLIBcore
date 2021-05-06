//
// Project     : HLib
// File        : TMaxwellRHS.cc
// Description : classes for building RHS-vectors in Maxwell applications
// Author      : Ronald Kriemann, Jonas Ballani
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/base/error.hh"
#include "hpro/bem/TGaussQuad.hh"
#include "hpro/bem/TConstEdgeFnSpace.hh"

#include "hpro/bem/TMaxwellRHS.hh"

namespace HLIB
{

using std::vector;
using std::unique_ptr;
using std::make_unique;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TMaxwellEFIERHS
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

template < typename T_fnspace >
TMaxwellEFIERHS< T_fnspace >::TMaxwellEFIERHS ( const uint  aquad_order )
        : TBEMRHS< fnspace_t, value_t >(),
          _quad_order( aquad_order )
{}

///////////////////////////////////////////////////////////////
//
// build RHS-vector for given <cluster> over <grid> and RHS
// function <rhs>; by <perm_i2e> one can set a permutation from
// internal to external numbering
//
template < typename T_fnspace >
unique_ptr< TVector >
TMaxwellEFIERHS< T_fnspace >::build ( const T_fnspace *                fnspace,
                                      const TBEMFunction< value_t > *  rhs,
                                      const TPermutation *             perm ) const
{
    using  real_t = typename real_type< value_t >::type_t;

    if ( fnspace == nullptr )
        HERROR( ERR_ARG, "(TMaxwellEFIERHS) build", "function space is nullptr" );

    if ( rhs     == nullptr )
        HERROR( ERR_ARG, "(TMaxwellEFIERHS) build", "RHS is nullptr" );

    //
    // build vector
    //

    unique_ptr< TScalarVector >  v;           // resulting vector
    vector< T3Point >            t( 3 );      // holds triangle coordinates
    vector< T3Point >            qpts;        // quadrature points for each triangle
    vector< double >             qwghts;      // quadrature weights for each triangle
    TTriGaussQuad                quad;        // Gauss quadrature rules for triangles
    vector< T3Point >            phi_values;  // basis function values
    TGrid::triangle_t            vid;         // vertex indices of triangle

    v = make_unique< TScalarVector >( fnspace->n_indices(), 0, true );

    for ( idx_t  idx = 0; idx < idx_t(fnspace->n_indices()); idx++ )
    {
        const idx_t  pidx = ( perm != nullptr ? perm->permute( idx ) : idx );

        //
        // integrate rhs function over support of index
        //

        for ( auto  tri : fnspace->support( idx ) )
        {
            const real_t  J = real( fnspace->grid()->tri_size( tri ) );

            vid = fnspace->grid()->triangle( tri );
            
            for ( uint vtx = 0; vtx < 3; vtx++ )
                t[vtx] = fnspace->grid()->vertex( vid.vtx[vtx] );

            const T3Point  n = fnspace->grid()->tri_normal( tri );

            // get quadrature points
            quad.build( _quad_order, t, qpts, qwghts );

            const size_t  npts = qpts.size();

            if ( phi_values.size() < npts )
                phi_values.resize( npts );

            // evaluate basis function
            fnspace->eval_basis( idx, vid, qpts, phi_values );

            const real_t  scale_phi = fnspace->scaling_factor_basis( idx, tri );

            // integrate
            complex  f = complex(0);

            for ( size_t  j = 0; j < npts; j++ )
            {
                value_t  fn_val  = rhs->eval( qpts[j], n );
                complex  dot_val = ( fn_val.val[0] * real_t( phi_values[j][0] ) +
                                     fn_val.val[1] * real_t( phi_values[j][1] ) +
                                     fn_val.val[2] * real_t( phi_values[j][2] ) );
                                     
                f += real_t(qwghts[j]) * dot_val;
            }// for
            
            v->add_centry( pidx, J * f * scale_phi );
        }// for
    }// for

    return std::unique_ptr< TVector >( v.release() );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TMaxwellMFIERHS
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

template < typename T_fnspace >
TMaxwellMFIERHS< T_fnspace >::TMaxwellMFIERHS ( const uint  aquad_order )
        : TBEMRHS< fnspace_t, value_t >(),
          _quad_order( aquad_order )
{}  

///////////////////////////////////////////////////////////////
//
// build RHS-vector for given <cluster> over <grid> and RHS
// function <rhs>; by <perm_i2e> one can set a permutation from
// internal to external numbering
//
template < typename T_fnspace >
unique_ptr< TVector >
TMaxwellMFIERHS< T_fnspace >::build ( const T_fnspace *                fnspace,
                                      const TBEMFunction< value_t > *  rhs,
                                      const TPermutation *             perm ) const
{
    using  real_t = typename real_type< value_t >::type_t;

    if ( fnspace == nullptr )
        HERROR( ERR_ARG, "(TMaxwellMFIERHS) build", "function space is nullptr" );

    if ( rhs     == nullptr )
        HERROR( ERR_ARG, "(TMaxwellMFIERHS) build", "RHS is nullptr" );

    //
    // build vector
    //

    unique_ptr< TScalarVector >  v;           // resulting vector
    vector< T3Point >            t( 3 );      // holds triangle coordinates
    vector< T3Point >            qpts;        // quadrature points for each triangle
    vector< double >             qwghts;      // quadrature weights for each triangle
    TTriGaussQuad                quad;        // Gauss quadrature rules for triangles
    vector< T3Point >            phi_values;  // basis function values
    TGrid::triangle_t            vid;         // vertex indices of triangle

    v = make_unique< TScalarVector >( fnspace->n_indices(), 0, true );

    for ( idx_t  idx = 0; idx < idx_t(fnspace->n_indices()); idx++ )
    {
        const idx_t  pidx = ( perm != nullptr ? perm->permute( idx ) : idx );

        //
        // integrate rhs function over support of index
        //

        for ( auto  tri : fnspace->support( idx ) )
        {
            const real_t  J = real( fnspace->grid()->tri_size( tri ) );

            vid = fnspace->grid()->triangle( tri );
            
            for ( uint vtx = 0; vtx < 3; vtx++ )
                t[vtx] = fnspace->grid()->vertex( vid.vtx[vtx] );

            const T3Point  n = fnspace->grid()->tri_normal( tri );

            // get quadrature points
            quad.build( _quad_order, t, qpts, qwghts );

            const size_t  npts = qpts.size();

            if ( phi_values.size() < npts )
                phi_values.resize( npts );

            // evaluate basis function
            fnspace->eval_basis( idx, vid, qpts, phi_values );

            const real_t  scale_phi = fnspace->scaling_factor_basis( idx, tri );

            // integrate
            complex  f = complex(0);

            for ( size_t  j = 0; j < npts; j++ )
            {
                value_t        fn_val  = rhs->eval( qpts[j], n );
                const T3Point  n_x_phi = cross( n, phi_values[j] );
                complex        dot_val = ( fn_val.val[0] * real_t( n_x_phi[0] ) +
                                           fn_val.val[1] * real_t( n_x_phi[1] ) +
                                           fn_val.val[2] * real_t( n_x_phi[2] ) );

                f += real(qwghts[j]) * dot_val;
            }// for

            v->add_centry( pidx, J * f * (-scale_phi) );
        }// for
    }// for

    return std::unique_ptr< TVector >( v.release() );
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instantiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template class TMaxwellEFIERHS< TConstEdgeFnSpace >;
template class TMaxwellMFIERHS< TConstEdgeFnSpace >;

}// namespace
