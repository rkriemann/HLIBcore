//
// Project     : HLib
// File        : TBEMBF.cc
// Description : classes for bilinearforms in BEM-applications
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <set>

#include "unordered_map.hh"

#include "hpro/base/error.hh"
#include "hpro/base/config.hh"
#include "hpro/bem/TConstEdgeFnSpace.hh"
#include "hpro/bem/TSauterTriQuad.hh"

#include "hpro/bem/TQuadBEMBF.hh"

namespace HLIB
{

using namespace std;

namespace B = BLAS;

//
// ctor
//
template < typename T_ansatzsp, typename T_testsp, typename T_val >
TQuadBEMBF< T_ansatzsp, T_testsp, T_val >::TQuadBEMBF ( const T_ansatzsp *  aansatzsp,
                                                        const T_testsp *    atestsp,
                                                        const uint          aorder,
                                                        const bool          dist_ada )
        : TBEMBF< T_ansatzsp, T_testsp, T_val >( aansatzsp, atestsp )
        , _quad_order( aorder )
        , _quad_dist_adaptive( dist_ada )
{
    //
    // build quad-points and with up to max-order
    //

    TSauterTriQuad  quad;

    _quad_rules.resize( 4 );
    
    for ( uint  ncommon = 0; ncommon < 4; ++ncommon )
    {
        _quad_rules[ ncommon ].resize(  _quad_order + 1 );
        
        for ( uint  order = 1; order < _quad_order + 1; ++order )
        {
            std::vector< TSauterTriQuad::rule_t >  rule;
            
            quad.build( ncommon, order, rule  );

            //
            // copy data into decoupled form
            //

            const size_t  npts        = rule.size();
            const size_t  padded_npts = CFG::Mach::simd_padded_size< real >( npts );

            _quad_rules[ ncommon ][ order ].npts = npts;
            
            _quad_rules[ ncommon ][ order ].pts1.resize( padded_npts );
            _quad_rules[ ncommon ][ order ].pts2.resize( padded_npts );
            _quad_rules[ ncommon ][ order ].w.resize( padded_npts );

            _quad_rules[ ncommon ][ order ].x1.resize( padded_npts );
            _quad_rules[ ncommon ][ order ].x2.resize( padded_npts );
            _quad_rules[ ncommon ][ order ].y1.resize( padded_npts );
            _quad_rules[ ncommon ][ order ].y2.resize( padded_npts );

            size_t  i = 0;

            // copy "real" quadrature data
            for ( ; i < npts; ++i )
            {
                _quad_rules[ ncommon ][ order ].pts1[i] = rule[i].pttri1;
                _quad_rules[ ncommon ][ order ].pts2[i] = rule[i].pttri2;
                _quad_rules[ ncommon ][ order ].w[i]    = rule[i].wght;

                _quad_rules[ ncommon ][ order ].x1[i]   = rule[i].pttri1[0];
                _quad_rules[ ncommon ][ order ].y1[i]   = rule[i].pttri1[1];
                _quad_rules[ ncommon ][ order ].x2[i]   = rule[i].pttri2[0];
                _quad_rules[ ncommon ][ order ].y2[i]   = rule[i].pttri2[1];
            }// for

            // fill rest with zero (only weights to avoid division by zero in coordinate computations)
            for ( ; i < padded_npts; ++i )
            {
                _quad_rules[ ncommon ][ order ].pts1[i] = rule[0].pttri1;
                _quad_rules[ ncommon ][ order ].pts2[i] = rule[0].pttri2;
                _quad_rules[ ncommon ][ order ].w[i]    = 0.0;

                _quad_rules[ ncommon ][ order ].x1[i]   = rule[0].pttri1[0];
                _quad_rules[ ncommon ][ order ].y1[i]   = rule[0].pttri1[1];
                _quad_rules[ ncommon ][ order ].x2[i]   = rule[0].pttri2[0];
                _quad_rules[ ncommon ][ order ].y2[i]   = rule[0].pttri2[1];
            }// for
        }// for
    }// for
}

//////////////////////////////////////
//
// evaluate bilinearform
//

//
// evaluate subblock defined by \a row_ind × \a col_ind; the indices
// in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
// contiguous
//
template < typename T_ansatzsp, typename T_testsp, typename T_val >
void
TQuadBEMBF< T_ansatzsp, T_testsp, T_val >::eval  ( const vector< idx_t > &,
                                                   const vector< idx_t > &,
                                                   BLAS::Matrix< T_val > & ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// reorder triangle vertices such that common vertices are ordered
// from 0,…,#common-1; number of common vertices is returned
//
template < typename T_ansatzsp, typename T_testsp, typename T_val >
uint
TQuadBEMBF< T_ansatzsp, T_testsp, T_val >::reorder_common ( idx_t *  vtx0idxs,
                                                            idx_t *  vtx1idxs ) const
{
    //
    // count number of common vertices
    //
    
    uint  ncommon = 0;

    for ( uint i = 0; i < 3; i++ )
    {
        for ( uint j = 0; j < 3; j++ )
        {
            if ( vtx0idxs[i] == vtx1idxs[j] )
            {
                swap( vtx0idxs[ncommon], vtx0idxs[i] );
                swap( vtx1idxs[ncommon], vtx1idxs[j] );
                ncommon++;
                break;
            }// if
        }// for
    }// for
    
    if ( ncommon > 3 )
        HERROR( ERR_CONSISTENCY, "(TQuadBEMBF) reorder_common",
               "triangles have more than 3 common vertices" );

    return ncommon;
}

//
// adjust quadrature order depending on diameter and distance of triangles
//
template < typename T_ansatzsp, typename T_testsp, typename T_val >
uint
TQuadBEMBF< T_ansatzsp, T_testsp, T_val >::adjust_order ( const idx_t *  vtx0idxs,
                                                          const idx_t *  vtx1idxs,
                                                          const uint     order ) const
{
    const TGrid *  ansatz_grid = this->ansatz_space()->grid();
    const TGrid *  test_grid   = this->test_space()->grid();
    
    //
    // compute bounding box of vertex set
    //
    
    T3Point  bbmin0( ansatz_grid->vertex( vtx0idxs[0] ) );
    T3Point  bbmax0( bbmin0 );
    T3Point  bbmin1( test_grid->vertex( vtx1idxs[0] ) );
    T3Point  bbmax1( bbmin1 );

    for ( uint  j = 1; j < 3; j++ )
    {
        for ( uint  i = 0; i < 3; i++ )
        {
            bbmin0[i] = min( bbmin0[i], ansatz_grid->vertex( vtx0idxs[j] )[i] );
            bbmax0[i] = max( bbmax0[i], ansatz_grid->vertex( vtx0idxs[j] )[i] );
            bbmin1[i] = min( bbmin1[i], test_grid->vertex( vtx1idxs[j] )[i] );
            bbmax1[i] = max( bbmax1[i], test_grid->vertex( vtx1idxs[j] )[i] );
        }// for
    }// for
    
    const double  diam_t0 = norm2( bbmax0 - bbmin0 );
    const double  diam_t1 = norm2( bbmax1 - bbmin1 );
    double        dist    = 0.0;
    
#if 1

    //
    // determine distance between bounding boxes
    //
    
    for ( uint  i = 0; i < 3; i++ )
    {
        const double min0 = bbmin0[i];
        const double max0 = bbmax0[i];
        const double min1 = bbmin1[i];
        const double max1 = bbmax1[i];
        
        if      (min0 > max1) dist += Math::square( min0 - max1 );
        else if (max0 < min1) dist += Math::square( min1 - max0 );
    }// for
    
    dist = Math::sqrt( dist );
    
#else

    //
    // distance = ∥ dist( center0 - center1 ) ∥₂
    //
    
    dist  = norm2( 0.5 * ((bbmax0 + bbmin0) - (bbmax1 + bbmin1)) );
    dist -= 0.5 * ( diam_t0 + diam_t1 );
    
#endif

    //
    // inflate diameter and reduce order as long as distance is still
    // large compared to diameter
    //
    
    double  diam   = 5.0 * max( diam_t0, diam_t1 );
    uint    torder = order;
    
    while (( dist > diam ) && ( torder > 1 ))
    {
        diam *= 5.0;
        torder--;
    }// while

    return torder;
}    

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instantiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template class TQuadBEMBF< TConstFnSpace,  TConstFnSpace,  real >;
template class TQuadBEMBF< TConstFnSpace,  TConstFnSpace,  complex >;
template class TQuadBEMBF< TConstFnSpace,  TLinearFnSpace, real >;
template class TQuadBEMBF< TConstFnSpace,  TLinearFnSpace, complex >;
template class TQuadBEMBF< TLinearFnSpace, TConstFnSpace,  real >;
template class TQuadBEMBF< TLinearFnSpace, TConstFnSpace,  complex >;
template class TQuadBEMBF< TLinearFnSpace, TLinearFnSpace, real >;
template class TQuadBEMBF< TLinearFnSpace, TLinearFnSpace, complex >;

template class TQuadBEMBF< TConstEdgeFnSpace, TConstEdgeFnSpace, complex >;

}// namespace
