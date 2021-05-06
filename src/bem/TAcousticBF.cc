//
// Project     : HLib
// File        : TAcousticBF.cc
// Description : bilinearform for acoustic scattering
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/base/packed.hh"

#include "hpro/bem/TAcousticBF.hh"

namespace HLIB
{
    
using namespace std;

namespace B = BLAS;

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// external function definitions
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template < typename  T_ansatzsp,
           typename  T_testsp,
           typename  T_packed >
void
acoustic_simd ( const complex                ikappa,
                const idx_t                  tri0idx,
                const TGrid::triangle_t &    tri0,
                const TGrid::triangle_t &    tri1,
                const tripair_quad_rule_t *  rule,
                const T_ansatzsp *           ansatz_sp,
                const T_testsp *             test_sp,
                std::vector< complex > &     values );

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// kernel function: default implementation
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename  T_ansatzsp,
           typename  T_testsp >
void
acoustic_flt ( const complex                ikappa,
               const idx_t                  tri0idx,
               const TGrid::triangle_t &    tri0,
               const TGrid::triangle_t &    tri1,
               const tripair_quad_rule_t *  rule,
               const T_ansatzsp *           ansatz_sp,
               const T_testsp *             test_sp,
               std::vector< complex > &     values )
{
    const real     ONE_OVER_2PI = 1.0 / (2.0 * Math::pi< double >());
    const size_t   n_pts        = rule->npts;
    const T3Point  v00( ansatz_sp->grid()->vertex( tri0.vtx[0] ) );
    const T3Point  v01( ansatz_sp->grid()->vertex( tri0.vtx[1] ) );
    const T3Point  v02( ansatz_sp->grid()->vertex( tri0.vtx[2] ) );
    const T3Point  v10( test_sp->grid()->vertex( tri1.vtx[0] ) );
    const T3Point  v11( test_sp->grid()->vertex( tri1.vtx[1] ) );
    const T3Point  v12( test_sp->grid()->vertex( tri1.vtx[2] ) );
    
    #pragma omp simd
    for ( size_t  i = 0; i < n_pts; i++ )
    {
        // get quadrature points and
        // transform to local triangles
        const auto     a1 = rule->pts1[i].x();
        const auto     a2 = rule->pts1[i].y();
        const auto     a0 = 1.0 - a1 - a2;
        
        const auto     b1 = rule->pts2[i].x();
        const auto     b2 = rule->pts2[i].y();
        const auto     b0 = 1.0 - b1 - b2;

        const T3Point  x = a0 * v00 + a1 * v01 + a2 * v02;
        const T3Point  y = b0 * v10 + b1 * v11 + b2 * v12;
        const T3Point  n( ansatz_sp->grid()->tri_normal( tri0idx, tri0, a1, a2 ) );

        //
        //   (i·κ·|x-y|)
        // e            (1 - i·κ·|x-y|) <x-y,n>
        // ────────────────────────────────────
        //                |x-y|³
        //
        //
        
        const T3Point  xmy( x - y );
        const auto     norm_xmy2 = real( dot( xmy, xmy ));   // |x-y|²
        const auto     norm_xmy  = Math::sqrt( norm_xmy2 );  // |x-y|
        const auto     norm_xmy3 = norm_xmy2 * norm_xmy;     // |x-y|³
        const auto     xmydotn   = real( dot( xmy, n ) );    // <x-y,n(x)>
        const auto     ikxmy     = ikappa * norm_xmy;

        values[ i ] = ONE_OVER_2PI * Math::exp( ikxmy ) * ( real(1) - ikxmy ) * xmydotn / norm_xmy3;
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TAcousticScatterBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template < typename  T_ansatzsp,
           typename  T_testsp >
TAcousticScatterBF< T_ansatzsp, T_testsp >::TAcousticScatterBF ( const complex       kappa,
                                                                 const ansatzsp_t *  aansatzsp,
                                                                 const testsp_t *    atestsp,
                                                                 const uint          quad_order )
  : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, complex >( aansatzsp, atestsp, quad_order )
  , _ikappa( complex(0,1) * kappa )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::BEM::use_simd_avx512f && CFG::Mach::has_avx512f() )
        {
            HINFO( "(TAcousticScatterBF) using AVX512F kernel" );
            
            _kernel_fn = acoustic_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX512F > >;
        }// if
        else if ( CFG::BEM::use_simd_mic && CFG::Mach::has_mic() )
        {
            HINFO( "(TAcousticScatterBF) using MIC kernel" );
            
            _kernel_fn = acoustic_simd< T_ansatzsp, T_testsp, packed< real, ISA_MIC > >;
        }// if
        else if ( CFG::BEM::use_simd_avx2 && CFG::Mach::has_avx2() )
        {
            HINFO( "(TAcousticScatterBF) using AVX2 kernel" );
            
            _kernel_fn = acoustic_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX2 > >;
        }// if
        else if ( CFG::BEM::use_simd_avx && CFG::Mach::has_avx() )
        {
            HINFO( "(TAcousticScatterBF) using AVX kernel" );
            
            _kernel_fn = acoustic_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX > >;
        }// if
        else if ( CFG::BEM::use_simd_sse3 && CFG::Mach::has_sse3() )
        {
            HINFO( "(TAcousticScatterBF) using SSE3 kernel" );
            
            _kernel_fn = acoustic_simd< T_ansatzsp, T_testsp, packed< real, ISA_SSE3 > >;
        }// if
        else if ( CFG::BEM::use_simd_vsx && CFG::Mach::has_vsx() )
        {
            HINFO( "(TAcousticScatterBF) using VSX kernel" );
            
            _kernel_fn = acoustic_simd< T_ansatzsp, T_testsp, packed< real, ISA_VSX > >;
        }// if
        else if ( CFG::BEM::use_simd_neon && CFG::Mach::has_neon() )
        {
            HINFO( "(TAcousticScatterBF) using NEON kernel" );
            
            _kernel_fn = acoustic_simd< T_ansatzsp, T_testsp, packed< real, ISA_NEON > >;
        }// if
        else if ( CFG::BEM::use_simd_neon && CFG::Mach::has_neon() )
        {
            HINFO( "(TAcousticScatterBF) using NEON kernel" );
            
            _kernel_fn = acoustic_simd< T_ansatzsp, T_testsp, packed< real, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TAcousticScatterBF) using standard kernel" );
            
            _kernel_fn = acoustic_flt;
        }// else
    }// if
    else
    {
        HINFO( "(TAcousticScatterBF) using standard kernel" );
            
        _kernel_fn = acoustic_flt;
    }// else
}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
TAcousticScatterBF< T_ansatzsp, T_testsp >::eval_kernel ( const idx_t                  tri0idx,
                                                          const idx_t,
                                                          const TGrid::triangle_t &    tri0,
                                                          const TGrid::triangle_t &    tri1,
                                                          const tripair_quad_rule_t *  rule,
                                                          std::vector< complex > &     values ) const
{
    _kernel_fn( _ikappa, tri0idx, tri0, tri1, rule, this->ansatz_space(), this->test_space(), values );
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template class TAcousticScatterBF< TConstFnSpace,  TConstFnSpace >;
template class TAcousticScatterBF< TConstFnSpace,  TLinearFnSpace >;
template class TAcousticScatterBF< TLinearFnSpace, TConstFnSpace >;
template class TAcousticScatterBF< TLinearFnSpace, TLinearFnSpace >;

}// namespace
