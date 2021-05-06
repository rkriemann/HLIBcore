//
// Project     : HLib
// File        : THelmholtzBF.cc
// Description : bilinear forms for Helmholtz operator
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/base/packed.hh"

#include "hpro/bem/TConstEdgeFnSpace.hh"

#include "hpro/bem/THelmholtzBF.hh"

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

template < typename T_ansatzsp,
           typename T_testsp >
void
helmholtz_slp_omp ( const TGrid::triangle_t &   tri0,
                    const TGrid::triangle_t &   tri1,
                    const tripair_quad_rule_t * rule,
                    const complex               ikappa,
                    const T_ansatzsp *          ansatz_sp,
                    const T_testsp *            test_sp,
                    std::vector< complex > &    values );

template < typename T_ansatzsp,
           typename T_testsp >
void
helmholtz_dlp_omp ( const idx_t                 tri1_id,
                    const TGrid::triangle_t &   tri0,
                    const TGrid::triangle_t &   tri1,
                    const tripair_quad_rule_t * rule,
                    const complex               ikappa,
                    const T_ansatzsp *          ansatz_sp,
                    const T_testsp *            test_sp,
                    std::vector< complex > &    values );


#define  HELMHOLTZ_SLP( suffix )                                        \
    template < typename  T_ansatzsp,                                    \
               typename  T_testsp,                                      \
               typename  T_packed >                                     \
    void                                                                \
    helmholtz_slp_##suffix ( const TGrid::triangle_t &  tri0,           \
                             const TGrid::triangle_t &  tri1,           \
                             const tripair_quad_rule_t *  rule,         \
                             const complex              ikappa,         \
                             const T_ansatzsp *         ansatz_sp,      \
                             const T_testsp *           test_sp,        \
                             vector< complex > &        values )

HELMHOLTZ_SLP( simd );
HELMHOLTZ_SLP( re_simd );
HELMHOLTZ_SLP( im_simd );

#define  HELMHOLTZ_DLP( suffix )                                        \
    template < typename  T_ansatzsp,                                    \
               typename  T_testsp,                                      \
               typename  T_packed >                                     \
    void                                                                \
    helmholtz_dlp_##suffix ( const idx_t                tri1_id,        \
                             const TGrid::triangle_t &  tri0,           \
                             const TGrid::triangle_t &  tri1,           \
                             const tripair_quad_rule_t *  rule,         \
                             const complex              ikappa,         \
                             const T_ansatzsp *         ansatz_sp,      \
                             const T_testsp *           test_sp,        \
                             vector< complex > &        values )

HELMHOLTZ_DLP( simd );
HELMHOLTZ_DLP( re_simd );
HELMHOLTZ_DLP( im_simd );

#define  HELMHOLTZ_SLP_HCA( suffix )                                    \
    template < typename  T_packed >                                     \
    void                                                                \
    helmholtz_slp_eval_dx_##suffix ( const complex &            ikappa, \
                                     const tri_quad_rule_t &    quad_rule, \
                                     const T3Point              x[3],   \
                                     const T3Point &            y,      \
                                     vector< complex > &        values )

HELMHOLTZ_SLP_HCA( simd );
HELMHOLTZ_SLP_HCA( re_simd );
HELMHOLTZ_SLP_HCA( im_simd );

#define  HELMHOLTZ_DLP_HCA( suffix )                                    \
    template < typename  T_packed >                                     \
    void                                                                \
    helmholtz_dlp_eval_dy_##suffix ( const complex &            ikappa, \
                                     const tri_quad_rule_t &    quad_rule, \
                                     const T3Point &            x,      \
                                     const T3Point              y[3],   \
                                     const T3Point &            normal, \
                                     vector< complex > &        values )

HELMHOLTZ_DLP_HCA( simd );
HELMHOLTZ_DLP_HCA( re_simd );
HELMHOLTZ_DLP_HCA( im_simd );

namespace
{

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// local constants
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

const real ONE_OVER_4PI = real(1) / (real(4) * Math::pi< real >());

}// namespace anonymous

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// THelmholtzSLPBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//
// standard evaluation of Helmholtz SLP kernel
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
helmholtz_slp_flt ( const TGrid::triangle_t &   tri0,
                    const TGrid::triangle_t &   tri1,
                    const tripair_quad_rule_t * rule,
                    const complex               ikappa,
                    const T_ansatzsp *          ansatz_sp,
                    const T_testsp *            test_sp,
                    vector< complex > &         values )
{
    const size_t   n_pts = rule->npts;
    const T3Point  v00( ansatz_sp->grid()->vertex( tri0.vtx[0] ) );
    const T3Point  v01( ansatz_sp->grid()->vertex( tri0.vtx[1] ) );
    const T3Point  v02( ansatz_sp->grid()->vertex( tri0.vtx[2] ) );
    const T3Point  v10( test_sp->grid()->vertex( tri1.vtx[0] ) );
    const T3Point  v11( test_sp->grid()->vertex( tri1.vtx[1] ) );
    const T3Point  v12( test_sp->grid()->vertex( tri1.vtx[2] ) );

    #pragma omp simd
    for ( size_t  i = 0; i < n_pts; i++ )
    {
        const double  a1 = rule->x1[i];
        const double  a2 = rule->y1[i];
        const double  a0 = 1.0 - a1 - a2;
        
        const double  b1 = rule->x2[i];
        const double  b2 = rule->y2[i];
        const double  b0 = 1.0 - b1 - b2;

        // transform coordinates
        const T3Point  x    = a0 * v00 + a1 * v01 + a2 * v02;
        const T3Point  y    = b0 * v10 + b1 * v11 + b2 * v12;
        const auto     dist = real( dot( x-y ) );

        // e^( i * kappa * |x-y| ) / |x-y|
        values[ i ] = ( ONE_OVER_4PI * Math::exp( ikappa * Math::sqrt( dist ) ) * Math::rsqrt( dist ) );
    }// for
}

//
// ctor
//
template < typename  T_ansatzsp,
           typename  T_testsp >
THelmholtzSLPBF< T_ansatzsp, T_testsp >::THelmholtzSLPBF ( const complex       kappa,
                                                           const ansatzsp_t *  aansatzsp,
                                                           const testsp_t *    atestsp,
                                                           const uint          quad_order )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, complex >( aansatzsp, atestsp, quad_order ),
          _ikappa( complex( 0, 1 ) * kappa )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            using  packed_t = packed< real, ISA_AVX512F >;
            
            HINFO( "(THelmholtzSLPBF) using AVX512F kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_slp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            using  packed_t = packed< real, ISA_MIC >;
            
            HINFO( "(THelmholtzSLPBF) using MIC kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_slp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            using  packed_t = packed< real, ISA_AVX2 >;
            
            HINFO( "(THelmholtzSLPBF) using AVX2 kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_slp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            using  packed_t = packed< real, ISA_AVX >;
            
            HINFO( "(THelmholtzSLPBF) using AVX kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_slp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            using  packed_t = packed< real, ISA_SSE3 >;
            
            HINFO( "(THelmholtzSLPBF) using SSE3 kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_slp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            using  packed_t = packed< real, ISA_VSX >;
            
            HINFO( "(THelmholtzSLPBF) using VSX kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_slp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            using  packed_t = packed< real, ISA_NEON >;
            
            HINFO( "(THelmholtzSLPBF) using NEON kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_slp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_slp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else
        {
            HINFO( "(THelmholtzSLPBF) using OpenMP kernel" );
            
            _kernel_fn = helmholtz_slp_omp;
        }// else
    }// if
    else
    {
        HINFO( "(THelmholtzSLPBF) using standard kernel" );
            
        _kernel_fn = helmholtz_slp_flt;
    }// else
}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
THelmholtzSLPBF< T_ansatzsp, T_testsp >::eval_kernel ( const                       idx_t,
                                                       const                       idx_t,
                                                       const TGrid::triangle_t &   tri0,
                                                       const TGrid::triangle_t &   tri1,
                                                       const tripair_quad_rule_t * rule,
                                                       std::vector< complex > &    values ) const
{
    _kernel_fn( tri0, tri1, rule, _ikappa, this->ansatz_space(), this->test_space(), values );
}

//
// return format of bilinear form, e.g. symmetric
//
template < typename  T_ansatzsp,
           typename  T_testsp >
matform_t
THelmholtzSLPBF< T_ansatzsp, T_testsp >::format () const
{
    if ( ( cptrcast( this->ansatz_space(), TFnSpace ) == cptrcast( this->test_space(), TFnSpace ) ) ||
         (( this->ansatz_space()->type() == this->test_space()->type() ) &&
          ( this->ansatz_space()->grid() == this->test_space()->grid() )) )
        return symmetric;
    else
        return unsymmetric;
}




/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// THelmholtzDLPBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//
// standard evaluation of Helmholtz DLP kernel
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
helmholtz_dlp_flt ( const idx_t                 tri1_id,
                    const TGrid::triangle_t &   tri0,
                    const TGrid::triangle_t &   tri1,
                    const tripair_quad_rule_t * rule,
                    const complex               ikappa,
                    const T_ansatzsp *          ansatz_sp,
                    const T_testsp *            test_sp,
                    vector< complex > &         values )
{
    const size_t   n_pts = rule->npts;
    const T3Point  v00( ansatz_sp->grid()->vertex( tri0.vtx[0] ) );
    const T3Point  v01( ansatz_sp->grid()->vertex( tri0.vtx[1] ) );
    const T3Point  v02( ansatz_sp->grid()->vertex( tri0.vtx[2] ) );
    const T3Point  v10( test_sp->grid()->vertex( tri1.vtx[0] ) );
    const T3Point  v11( test_sp->grid()->vertex( tri1.vtx[1] ) );
    const T3Point  v12( test_sp->grid()->vertex( tri1.vtx[2] ) );
    const T3Point  n( test_sp->grid()->tri_normal( tri1_id ) );

    #pragma omp simd
    for ( size_t  i = 0; i < n_pts; i++ )
    {
        const double  a1 = rule->x1[i];
        const double  a2 = rule->y1[i];
        const double  a0 = 1.0 - a1 - a2;
        
        const double  b1 = rule->x2[i];
        const double  b2 = rule->y2[i];
        const double  b0 = 1.0 - b1 - b2;

        // u = y - x
        const T3Point  x      = a0 * v00 + a1 * v01 + a2 * v02;
        const T3Point  y      = b0 * v10 + b1 * v11 + b2 * v12;
        const T3Point  ymx    = y - x;
        const double   sqdist = real( dot( ymx ) );       // |y-x|²

        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //

        const auto  dist      = Math::sqrt( sqdist );     // |y-x|
        const auto  cudist    = real( sqdist * dist );    // |y-x|³
        const auto  n_dot_ymx = real( dot( n, ymx ) );    // <n(y),y-x>
        const auto  ikymx     = ikappa * real( dist );    // i · κ · |x-y|
        
        values[ i ] = ONE_OVER_4PI * Math::exp( ikymx ) * ( ikymx - real(1) ) * n_dot_ymx / cudist;
    }// for
}

//
// ctor
//
template < typename  T_ansatzsp,
           typename  T_testsp >
THelmholtzDLPBF< T_ansatzsp, T_testsp >::THelmholtzDLPBF ( const complex       kappa,
                                                           const ansatzsp_t *  aansatzsp,
                                                           const testsp_t *    atestsp,
                                                           const uint          quad_order )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, complex >( aansatzsp, atestsp, quad_order ),
          _ikappa( complex( 0, 1 ) * kappa )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            using  packed_t = packed< real, ISA_AVX512F >;
            
            HINFO( "(THelmholtzDLPBF) using AVX512F kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_dlp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic  )
        {
            using  packed_t = packed< real, ISA_MIC >;
            
            HINFO( "(THelmholtzDLPBF) using MIC kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_dlp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            using  packed_t = packed< real, ISA_AVX2 >;
            
            HINFO( "(THelmholtzDLPBF) using AVX2 kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_dlp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            using  packed_t = packed< real, ISA_AVX >;
            
            HINFO( "(THelmholtzDLPBF) using AVX kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_dlp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            using  packed_t = packed< real, ISA_SSE3 >;
            
            HINFO( "(THelmholtzDLPBF) using SSE3 kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_dlp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            using  packed_t = packed< real, ISA_VSX >;
            
            HINFO( "(THelmholtzDLPBF) using VSX kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_dlp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            using  packed_t = packed< real, ISA_NEON >;
            
            HINFO( "(THelmholtzDLPBF) using NEON kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_re_simd< T_ansatzsp, T_testsp, packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _kernel_fn = helmholtz_dlp_im_simd< T_ansatzsp, T_testsp, packed_t >;
            else                                        _kernel_fn = helmholtz_dlp_simd< T_ansatzsp, T_testsp, packed_t >;
        }// if
        else
        {
            HINFO( "(THelmholtzDLPBF) using OpenMP kernel" );
            
            _kernel_fn = helmholtz_dlp_omp;
        }// else
    }// if
    else
    {
        HINFO( "(THelmholtzDLPBF) using standard kernel" );
            
        _kernel_fn = helmholtz_dlp_flt;
    }// else
}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
THelmholtzDLPBF< T_ansatzsp, T_testsp >::eval_kernel ( const                       idx_t,
                                                       const idx_t                 tri1idx,
                                                       const TGrid::triangle_t &   tri0,
                                                       const TGrid::triangle_t &   tri1,
                                                       const tripair_quad_rule_t * rule,
                                                       std::vector< complex > &   values ) const
{
    _kernel_fn( tri1idx, tri0, tri1, rule, _ikappa, this->ansatz_space(), this->test_space(), values );
}

//
// return format of bilinear form, e.g. symmetric
//
template < typename  T_ansatzsp,
           typename  T_testsp >
matform_t
THelmholtzDLPBF< T_ansatzsp, T_testsp >::format () const
{
    return unsymmetric;
}



///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// auxiliary functions for SLP and DLP generators
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// computes x-integral of SLP
//
void
helmholtz_slp_eval_dx_flt ( const complex &          ikappa,
                            const tri_quad_rule_t &  quad_rule,
                            const T3Point            vx[3],
                            const T3Point &          vy,
                            vector< complex > &      values )
{
    const size_t  npts = quad_rule.npts;
    
    #pragma omp simd
    for ( size_t  k = 0; k < npts; ++k )
    {
        const double  a2 = quad_rule.y[k];
        const double  a1 = quad_rule.x[k];
        const double  a0 = 1.0 - a1 - a2;

        //
        // evaluate SLP kernel
        //

        // d = y - x
        const T3Point  xmy( vy - ( a0 * vx[0] + a1 * vx[1] + a2 * vx[2] ) );
        const auto     ddotd = real( Math::sqrt( dot( xmy, xmy ) ) );

        // e^( i * kappa * | x - y | ) / | x - y |
        values[k] = ONE_OVER_4PI * Math::exp( ikappa * ddotd ) / ddotd;
    }// for
}

//
// computes y-integral of DLP
//
void
helmholtz_dlp_eval_dy_flt ( const complex &          ikappa,
                            const tri_quad_rule_t &  quad_rule,
                            const T3Point &          vx,
                            const T3Point            vy[3],
                            const T3Point &          normal,
                            vector< complex > &      values )
{
    const size_t  npts = quad_rule.npts;
    
    #pragma omp simd
    for ( size_t  k = 0; k < npts; ++k )
    {
        const double  a2 = quad_rule.y[k];
        const double  a1 = quad_rule.x[k];
        const double  a0 = 1.0 - a1 - a2;

        //
        // evaluate DLP kernel
        //
                
        // d = y - x
        const T3Point  ymx( ( a0 * vy[0] + a1 * vy[1] + a2 * vy[2] ) - vx );
        const auto     dist2  = real( dot( ymx, ymx ) );
        const auto     dist   = Math::sqrt( dist2 );
        const auto     dist3  = dist * dist2;
        const auto     ikdist = ikappa * dist;
                    
        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //

        values[k] = ( ONE_OVER_4PI * ( Math::exp( ikdist ) * ( ikdist - real(1) ) * real(dot( normal, ymx )) )
                      / dist3 );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// THelmholtzSLPGenFn
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// constructor
//
template < typename  T_ansatzsp,
           typename  T_testsp >
THelmholtzSLPGenFn< T_ansatzsp, T_testsp >::THelmholtzSLPGenFn ( const complex         kappa,
                                                                 const ansatzsp_t *    ansatzsp,
                                                                 const testsp_t *      testsp,
                                                                 const TPermutation *  row_perm_i2e,
                                                                 const TPermutation *  col_perm_i2e,
                                                                 const uint            quad_order )
        : TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, complex >( ansatzsp, testsp, quad_order,
                                                                    row_perm_i2e, col_perm_i2e ),
          _ikappa( complex( 0, 1 ) * kappa )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            using  packed_t = packed< real, ISA_AVX512F >;
            
            HINFO( "(THelmholtzSLPGenFn) using AVX512F kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
            else                                        _eval_dxy_impl = helmholtz_slp_eval_dx_simd< packed_t >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            using  packed_t = packed< real, ISA_MIC >;
            
            HINFO( "(THelmholtzSLPGenFn) using MIC kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
            else                                        _eval_dxy_impl = helmholtz_slp_eval_dx_simd< packed_t >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            using  packed_t = packed< real, ISA_AVX2 >;
            
            HINFO( "(THelmholtzSLPGenFn) using AVX2 kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
            else                                        _eval_dxy_impl = helmholtz_slp_eval_dx_simd< packed_t >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            using  packed_t = packed< real, ISA_AVX >;
            
            HINFO( "(THelmholtzSLPGenFn) using AVX kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
            else                                        _eval_dxy_impl = helmholtz_slp_eval_dx_simd< packed_t >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            using  packed_t = packed< real, ISA_SSE3 >;
            
            HINFO( "(THelmholtzSLPGenFn) using SSE3 kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
            else                                        _eval_dxy_impl = helmholtz_slp_eval_dx_simd< packed_t >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            using  packed_t = packed< real, ISA_VSX >;
            
            HINFO( "(THelmholtzSLPGenFn) using VSX kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
            else                                        _eval_dxy_impl = helmholtz_slp_eval_dx_simd< packed_t >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            using  packed_t = packed< real, ISA_NEON >;
            
            HINFO( "(THelmholtzSLPGenFn) using NEON kernel" );
            
            if      ( std::real( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
            else if ( std::imag( _ikappa ) == real(0) ) _eval_dxy_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
            else                                        _eval_dxy_impl = helmholtz_slp_eval_dx_simd< packed_t >;
        }// if
        else
        {
            HINFO( "(THelmholtzSLPGenFn) using standard kernel" );
            
            _eval_dxy_impl = helmholtz_slp_eval_dx_flt;
        }// else
    }// if
    else
    {
        HINFO( "(THelmholtzSLPGenFn) using standard kernel" );
            
        _eval_dxy_impl = helmholtz_slp_eval_dx_flt;
    }// else
}
        
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// THelmholtzDLPGenFn
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// constructor
//
template < typename  T_ansatzsp,
           typename  T_testsp >
THelmholtzDLPGenFn< T_ansatzsp, T_testsp >::THelmholtzDLPGenFn ( const complex         kappa,
                                                                 const ansatzsp_t *    ansatzsp,
                                                                 const testsp_t *      testsp,
                                                                 const TPermutation *  row_perm_i2e,
                                                                 const TPermutation *  col_perm_i2e,
                                                                 const uint            quad_order )
        : TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, complex >( ansatzsp, testsp, quad_order,
                                                                    row_perm_i2e, col_perm_i2e ),
          _ikappa( complex( 0, 1 ) * kappa )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            using  packed_t = packed< real, ISA_AVX512F >;
            
            HINFO( "(THelmholtzDLPGenFn) using AVX512F kernel" );
            
            if ( std::real( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            using  packed_t = packed< real, ISA_MIC >;
            
            HINFO( "(THelmholtzDLPGenFn) using MIC kernel" );
            
            if ( std::real( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            using  packed_t = packed< real, ISA_AVX2 >;
            
            HINFO( "(THelmholtzDLPGenFn) using AVX2 kernel" );
            
            if ( std::real( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            using  packed_t = packed< real, ISA_AVX >;
            
            HINFO( "(THelmholtzDLPGenFn) using AVX kernel" );
            
            if ( std::real( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            using  packed_t = packed< real, ISA_SSE3 >;
            
            HINFO( "(THelmholtzDLPGenFn) using SSE3 kernel" );
            
            if ( std::real( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            using  packed_t = packed< real, ISA_VSX >;
            
            HINFO( "(THelmholtzDLPGenFn) using VSX kernel" );
            
            if ( std::real( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< packed_t >;
            }// else
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            using  packed_t = packed< real, ISA_NEON >;
            
            HINFO( "(THelmholtzDLPGenFn) using NEON kernel" );
            
            if ( std::real( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_re_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_re_simd< packed_t >;
            }// if
            else if ( std::imag( _ikappa ) == real(0) )
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_im_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_im_simd< packed_t >;
            }// if
            else
            {
                _eval_dx_impl = helmholtz_slp_eval_dx_simd< packed_t >;
                _eval_dy_impl = helmholtz_dlp_eval_dy_simd< packed_t >;
            }// else
        }// if
        else
        {
            HINFO( "(THelmholtzDLPGenFn) using standard kernel" );
            
            _eval_dx_impl = helmholtz_slp_eval_dx_flt;
            _eval_dy_impl = helmholtz_dlp_eval_dy_flt;
        }// else
    }// if
    else
    {
        HINFO( "(THelmholtzDLPGenFn) using standard kernel" );
            
        _eval_dx_impl = helmholtz_slp_eval_dx_flt;
        _eval_dy_impl = helmholtz_dlp_eval_dy_flt;
    }// if
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template class THelmholtzSLPBF< TConstFnSpace,  TConstFnSpace >;
template class THelmholtzSLPBF< TConstFnSpace,  TLinearFnSpace >;
template class THelmholtzSLPBF< TLinearFnSpace, TConstFnSpace >;
template class THelmholtzSLPBF< TLinearFnSpace, TLinearFnSpace >;

template class THelmholtzDLPBF< TConstFnSpace,  TConstFnSpace >;
template class THelmholtzDLPBF< TConstFnSpace,  TLinearFnSpace >;
template class THelmholtzDLPBF< TLinearFnSpace, TConstFnSpace >;
template class THelmholtzDLPBF< TLinearFnSpace, TLinearFnSpace >;

template class THelmholtzSLPGenFn< TConstFnSpace,  TConstFnSpace >;
template class THelmholtzSLPGenFn< TConstFnSpace,  TLinearFnSpace >;
template class THelmholtzSLPGenFn< TLinearFnSpace, TConstFnSpace >;
template class THelmholtzSLPGenFn< TLinearFnSpace, TLinearFnSpace >;

template class THelmholtzDLPGenFn< TConstFnSpace,  TConstFnSpace >;
template class THelmholtzDLPGenFn< TConstFnSpace,  TLinearFnSpace >;
template class THelmholtzDLPGenFn< TLinearFnSpace, TConstFnSpace >;
template class THelmholtzDLPGenFn< TLinearFnSpace, TLinearFnSpace >;

//
// instantiate Helmholtz for TConstEdgeFnSpace for MaxwellEFIE
//
template
void
helmholtz_slp_flt < TConstEdgeFnSpace, TConstEdgeFnSpace > (
    const TGrid::triangle_t &   tri0,
    const TGrid::triangle_t &   tri1,
    const tripair_quad_rule_t * rule,
    const complex               ikappa,
    const TConstEdgeFnSpace *   ansatz_sp,
    const TConstEdgeFnSpace *   test_sp,
    vector< complex > &         values );

}// namespace
