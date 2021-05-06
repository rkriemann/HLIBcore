//
// Project     : HLib
// File        : TBEMRHS.cc
// Description : classes for building RHS-vectors in BEM-applications
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/base/error.hh"
#include "hpro/bem/TGaussQuad.hh"

#include "hpro/bem/TBEMRHS.hh"

namespace HLIB
{

using std::vector;
using std::unique_ptr;
using std::make_unique;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TQuadBEMRHS
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
//
// constructur and destructor
//

template < typename  T_fnspace,
           typename  T_val >
TQuadBEMRHS< T_fnspace, T_val >::TQuadBEMRHS ( const uint  aquad_order )
        : _quad_order( std::max( 1u, aquad_order ) )
{
}

///////////////////////////////////////////////////////////////
//
// build RHS-vector for given <cluster> over <grid> and RHS
// function <rhs>; by <perm_i2e> one can set a permutation from
// internal to external numbering
//
template < typename  T_fnspace,
           typename  T_val >
unique_ptr< TVector >
TQuadBEMRHS< T_fnspace, T_val >::build ( const T_fnspace *              fnspace,
                                         const TBEMFunction< T_val > *  rhs,
                                         const TPermutation *           perm ) const
{
    HASSERT( fnspace != nullptr, ERR_ARG, "(TQuadBEMRHS) build", "function space is nullptr" );
    HASSERT( rhs     != nullptr, ERR_ARG, "(TQuadBEMRHS) build", "RHS is nullptr" );

    //
    // build vector
    //

    vector< T2Point >            qpts;   // quadrature points for each triangle
    vector< double >             qwghts; // quadrature weights for each triangle
    TTriGaussQuad                quad;   // Gauss quadrature rules for triangles
    const TGrid *                grid = fnspace->grid();
    auto                         v    = make_unique< TScalarVector >( fnspace->n_indices(), 0, rhs->is_complex() );
    
    // get quadrature points
    quad.build( _quad_order, qpts, qwghts );

    const size_t  npts = qpts.size();
    
    for ( idx_t  idx = 0; idx < idx_t(fnspace->n_indices()); idx++ )
    {
        const idx_t  pidx = ( perm != nullptr ? perm->permute( idx ) : idx );
        value_t      val  = value_t(0);

        //
        // integrate rhs function over support of index
        //

        for ( auto  tri_idx : fnspace->support( idx ) )
        {
            const TGrid::triangle_t  tri = grid->triangle( tri_idx );
            const T3Point            v0  = grid->vertex( tri.vtx[0] );
            const T3Point            v1  = grid->vertex( tri.vtx[1] );
            const T3Point            v2  = grid->vertex( tri.vtx[2] );
            const real               J   = real( grid->tri_size( tri_idx ) );
            const T3Point            n   = grid->tri_normal( tri_idx );

            //
            // evaluate at quadrature points
            //
            
            value_t  f = value_t(0);
            
            for ( size_t  j = 0; j < npts; j++ )
            {
                // transform to triangle coordinates
                const double   a1 = qpts[j][0];
                const double   a2 = qpts[j][1];
                const double   a0 = 1.0 - a1 - a2;
                const T3Point  x  = a0 * v0 + a1 * v1 + a2 * v2;
                
                f += real( qwghts[j] ) * rhs->eval( x, n ) * fnspace->eval_basis_unit( idx, a1, a2, tri.vtx );
            }// for

            val += J * f;
        }// for
            
        if ( rhs->is_complex() ) v->set_centry( pidx, complex(val) );
        else                     v->set_entry(  pidx, std::real(val) );
    }// for

    return std::unique_ptr< TVector >( v.release() );
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template class TQuadBEMRHS< TConstFnSpace,  real >;
template class TQuadBEMRHS< TLinearFnSpace, real >;
template class TQuadBEMRHS< TConstFnSpace,  complex >;
template class TQuadBEMRHS< TLinearFnSpace, complex >;

}// namespace
