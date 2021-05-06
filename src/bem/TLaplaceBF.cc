//
// Project     : HLib
// File        : TLaplaceBF.cc
// Description : bilinearforms for Laplace operator
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/base/packed.hh"
#include "hpro/blas/lapack.hh"

#include "hpro/bem/TLaplaceBF.hh"

/////////////////////////////////////////////////////////////

namespace HLIB
{

using namespace std;

namespace B = BLAS;

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// external functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
laplace_slp_simd  ( const TGrid::triangle_t &   tri0,
                    const TGrid::triangle_t &   tri1,
                    const tripair_quad_rule_t * rule,
                    std::vector< real > &       values,
                    const T_ansatzsp *          ansatz_sp,
                    const T_testsp *            test_sp );

template < typename T_ansatzsp,
           typename T_testsp >
void
laplace_slp_omp   ( const TGrid::triangle_t &   tri0,
                    const TGrid::triangle_t &   tri1,
                    const tripair_quad_rule_t * rule,
                    std::vector< real > &       values,
                    const T_ansatzsp *          ansatz_sp,
                    const T_testsp *            test_sp );

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
laplace_dlp_simd  ( const TGrid::triangle_t &   tri0,
                    const TGrid::triangle_t &   tri1,
                    const tripair_quad_rule_t * rule,
                    std::vector< real > &       values,
                    const T_ansatzsp *          ansatz_sp,
                    const T_testsp *            test_sp,
                    const T3Point &             n_vec );

template < typename T_packed >
void
laplace_slp_eval_dx_simd  ( const tri_quad_rule_t &   quad_rule,
                            const T3Point             x[3],
                            const T3Point &           y,
                            vector< real > &          values );

template < typename T_packed >
void
laplace_dlp_eval_dy_simd  ( const tri_quad_rule_t &   quad_rule,
                            const T3Point &           x,
                            const T3Point             y[3],
                            const T3Point &           normal,
                            vector< real > &          values );

namespace
{

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// local constants
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

const real ONE_OVER_4PI = real( 1.0 / (4.0 * Math::pi< double >()) );

// define to disable all accelleration techniques, e.g. SSE3, AVX, MIC
//#define NO_ACCEL

}// namespace anonymous

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// Laplace SLP
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// machine dependent implementations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_ansatzsp,
           typename T_testsp >
void
laplace_slp_flt ( const TGrid::triangle_t &   tri0,
                  const TGrid::triangle_t &   tri1,
                  const tripair_quad_rule_t * rule,
                  std::vector< real > &       values,
                  const T_ansatzsp *          ansatz_sp,
                  const T_testsp *            test_sp )
{
    const size_t   n_pts = rule->npts;
    const T3Point  v00( ansatz_sp->grid()->vertex( tri0.vtx[0] ) );
    const T3Point  v01( ansatz_sp->grid()->vertex( tri0.vtx[1] ) );
    const T3Point  v02( ansatz_sp->grid()->vertex( tri0.vtx[2] ) );
    const T3Point  v10( test_sp->grid()->vertex( tri1.vtx[0] ) );
    const T3Point  v11( test_sp->grid()->vertex( tri1.vtx[1] ) );
    const T3Point  v12( test_sp->grid()->vertex( tri1.vtx[2] ) );
    T3Point        x, y, u;

    #pragma omp simd
    for ( size_t  i = 0; i < n_pts; ++i )
    {
        const double  a1 = rule->x1[i];
        const double  a2 = rule->y1[i];
        const double  a0 = real(1) - a1 - a2;
        
        const double  b1 = rule->x2[i];
        const double  b2 = rule->y2[i];
        const double  b0 = real(1) - b1 - b2;

        // u = x - y
        x = a0 * v00 + a1 * v01 + a2 * v02;
        y = b0 * v10 + b1 * v11 + b2 * v12;
        u = x - y;

        //
        //   1
        // ─────
        // |x-y|
        //
        
        values[ i ] = ( ONE_OVER_4PI * Math::rsqrt( real( dot( u, u ) ) ) );
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TLaplaceSLPBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template < typename  T_ansatzsp,
           typename  T_testsp >
TLaplaceSLPBF< T_ansatzsp, T_testsp >::TLaplaceSLPBF ( const ansatzsp_t *  aansatzsp,
                                                       const testsp_t *    atestsp,
                                                       const uint          quad_order )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, real >( aansatzsp, atestsp, quad_order )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f()  && CFG::BEM::use_simd_avx512f  )
        {
            HINFO( "(TLaplaceSLPBF) using AVX512F kernel" );
            
            _kernel_fn = laplace_slp_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic()  && CFG::BEM::use_simd_mic  )
        {
            HINFO( "(TLaplaceSLPBF) using MIC kernel" );
            
            _kernel_fn = laplace_slp_simd< T_ansatzsp, T_testsp, packed< real, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2()  && CFG::BEM::use_simd_avx2  )
        {
            HINFO( "(TLaplaceSLPBF) using AVX2 kernel" );
            
            _kernel_fn = laplace_slp_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx()  && CFG::BEM::use_simd_avx  )
        {
            HINFO( "(TLaplaceSLPBF) using AVX kernel" );
            
            _kernel_fn = laplace_slp_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TLaplaceSLPBF) using SSE3 kernel" );
            
            _kernel_fn = laplace_slp_simd< T_ansatzsp, T_testsp, packed< real, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TLaplaceSLPBF) using VSX kernel" );
            
            _kernel_fn = laplace_slp_simd< T_ansatzsp, T_testsp, packed< real, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TLaplaceSLPBF) using NEON kernel" );
            
            _kernel_fn = laplace_slp_simd< T_ansatzsp, T_testsp, packed< real, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TLaplaceSLPBF) using OpenMP kernel" );
            
            _kernel_fn = laplace_slp_omp;
        }// else
    }// if
    else
    {
        HINFO( "(TLaplaceSLPBF) using FPU kernel" );
            
        _kernel_fn = laplace_slp_flt;
    }// else
}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
TLaplaceSLPBF< T_ansatzsp, T_testsp >::eval_kernel ( const                       idx_t,
                                                     const                       idx_t,
                                                     const TGrid::triangle_t &   tri0,
                                                     const TGrid::triangle_t &   tri1,
                                                     const tripair_quad_rule_t * rule,
                                                     std::vector< real > &       values ) const
{
    _kernel_fn( tri0, tri1, rule, values, this->ansatz_space(), this->test_space() );
    ADD_FLOPS( rule->npts * ( 15 + 15 + 3 + 8 ) );
}

//
// return format of bilinear form, e.g. symmetric
//
template < typename  T_ansatzsp,
           typename  T_testsp >
matform_t
TLaplaceSLPBF< T_ansatzsp, T_testsp >::format () const
{
    if ( ( cptrcast( this->ansatz_space(), TFnSpace ) == cptrcast( this->test_space(), TFnSpace ) ) ||
         (( this->ansatz_space()->type() == this->test_space()->type() ) &&
          ( this->ansatz_space()->grid() == this->test_space()->grid() )) )
        return symmetric;
    else
        return unsymmetric;
}



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
// Laplace DLP
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// machine dependent implementations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_ansatzsp,
           typename T_testsp >
void
laplace_dlp_flt ( const TGrid::triangle_t &   tri0,
                  const TGrid::triangle_t &   tri1,
                  const tripair_quad_rule_t * rule,
                  std::vector< real > &       values,
                  const T_ansatzsp *          ansatz_sp,
                  const T_testsp *            test_sp,
                  const T3Point &             n )
{
    const size_t   n_pts = rule->npts;
    const T3Point  v00( ansatz_sp->grid()->vertex( tri0.vtx[0] ) );
    const T3Point  v01( ansatz_sp->grid()->vertex( tri0.vtx[1] ) );
    const T3Point  v02( ansatz_sp->grid()->vertex( tri0.vtx[2] ) );
    const T3Point  v10( test_sp->grid()->vertex( tri1.vtx[0] ) );
    const T3Point  v11( test_sp->grid()->vertex( tri1.vtx[1] ) );
    const T3Point  v12( test_sp->grid()->vertex( tri1.vtx[2] ) );
    T3Point        x, y, u;

    #pragma omp simd
    for ( size_t  i = 0; i < n_pts; i++ )
    {
        const double  a1 = rule->x1[i];
        const double  a2 = rule->y1[i];
        const double  a0 = 1.0 - a1 - a2;
        
        const double  b1 = rule->x2[i];
        const double  b2 = rule->y2[i];
        const double  b0 = 1.0 - b1 - b2;

        // u = x - y
        x = a0 * v00 + a1 * v01 + a2 * v02;
        y = b0 * v10 + b1 * v11 + b2 * v12;
        u = y - x;

        //
        // ⟨n,x-y⟩
        // ───────
        //  |x-y|³
        //

        const real  inxy = Math::rsqrt( real( dot( u, u ) ) );
        
        values[ i ] = ONE_OVER_4PI * real( dot( n, u ) ) * inxy * inxy * inxy;
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TLaplaceDLPBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

//
// ctor
//
template < typename  T_ansatzsp,
           typename  T_testsp >
TLaplaceDLPBF< T_ansatzsp, T_testsp >::TLaplaceDLPBF ( const ansatzsp_t *  aansatzsp,
                                                       const testsp_t *    atestsp,
                                                       const uint          quad_order )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, real >( aansatzsp, atestsp, quad_order )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f()  && CFG::BEM::use_simd_avx512f  )
        {
            HINFO( "(TLaplaceDLPBF) using AVX512F kernel" );
            
            _kernel_fn = laplace_dlp_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic()  && CFG::BEM::use_simd_mic  )
        {
            HINFO( "(TLaplaceDLPBF) using MIC kernel" );
            
            _kernel_fn = laplace_dlp_simd< T_ansatzsp, T_testsp, packed< real, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2()  && CFG::BEM::use_simd_avx2  )
        {
            HINFO( "(TLaplaceDLPBF) using AVX2 kernel" );
            
            _kernel_fn = laplace_dlp_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx()  && CFG::BEM::use_simd_avx  )
        {
            HINFO( "(TLaplaceDLPBF) using AVX kernel" );
            
            _kernel_fn = laplace_dlp_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TLaplaceDLPBF) using SSE3 kernel" );
            
            _kernel_fn = laplace_dlp_simd< T_ansatzsp, T_testsp, packed< real, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TLaplaceDLPBF) using VSX kernel" );
            
            _kernel_fn = laplace_dlp_simd< T_ansatzsp, T_testsp, packed< real, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TLaplaceDLPBF) using NEON kernel" );
            
            _kernel_fn = laplace_dlp_simd< T_ansatzsp, T_testsp, packed< real, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TLaplaceDLPBF) using FPU kernel" );
            
            _kernel_fn = laplace_dlp_flt;
        }// else
    }// if
    else
    {
        HINFO( "(TLaplaceDLPBF) using FPU kernel" );
            
        _kernel_fn = laplace_dlp_flt;
    }// else
}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
TLaplaceDLPBF< T_ansatzsp, T_testsp >::eval_kernel ( const idx_t ,
                                                     const idx_t                 tri1idx,
                                                     const TGrid::triangle_t &   tri0,
                                                     const TGrid::triangle_t &   tri1,
                                                     const tripair_quad_rule_t * rule,
                                                     std::vector< real > &       values ) const
{
    _kernel_fn( tri0, tri1, rule, values, this->ansatz_space(), this->test_space(),
                this->test_space()->grid()->tri_normal( tri1idx ) );
    ADD_FLOPS( rule->npts * ( 15 + 15 + 3 + 7 + 9 ) );
}


//
// return format of bilinear form, e.g. symmetric
//
template < typename  T_ansatzsp,
           typename  T_testsp >
matform_t
TLaplaceDLPBF< T_ansatzsp, T_testsp >::format () const
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
// evaluate dx-integral of SLP
//
void
laplace_slp_eval_dx_flt ( const tri_quad_rule_t &  quad_rule,
                          const T3Point            vx[3],
                          const T3Point &          vy,
                          vector< real > &         values )
{
    const size_t  npts  = quad_rule.npts;
    
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
        const T3Point  x = a0 * vx[0]  +  a1 * vx[1]  +  a2 * vx[2];
        const T3Point  d = vy - x;

        values[k] = ONE_OVER_4PI * Math::rsqrt( dot( d, d ) );
    }// for
}

//
// evaluate dy-integral of DLP
//
void
laplace_dlp_eval_dy_flt ( const tri_quad_rule_t &  quad_rule,
                          const T3Point &          vx,
                          const T3Point            vy[3],
                          const T3Point &          normal,
                          vector< real > &         values )
{
    const size_t  npts  = quad_rule.npts;
    
    //
    // fall back or handle remaining
    //

    #pragma omp simd
    for ( size_t  k = 0; k < npts; k++ )
    {
        const double  a2 = quad_rule.y[k];
        const double  a1 = quad_rule.x[k];
        const double  a0 = 1.0 - a1 - a2;

        //
        // evaluate DLP kernel
        //
                
        // d = y - x
        const T3Point  y = a0 * vy[0]  +  a1 * vy[1]  +  a2 * vy[2];
        const T3Point  d = y - vx;

        // r = 1 / |d| = 1 / sqrt( <d,d> )
        const real     r = Math::rsqrt( d.dot( d ) );
                    
        // store < n, d > / |d|^3
        values[k] = ONE_OVER_4PI * normal.dot( d ) * (r*r*r);
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// TLaplaceSLPGenFn
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// constructor
//
template < typename T_ansatzsp, typename T_testsp >
TLaplaceSLPGenFn< T_ansatzsp, T_testsp >::TLaplaceSLPGenFn ( const ansatzsp_t *    ansatzsp,
                                                             const testsp_t *      testsp,
                                                             const TPermutation *  row_perm_i2e,
                                                             const TPermutation *  col_perm_i2e,
                                                             const uint            quad_order )
        : TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, real >( ansatzsp,
                                                                 testsp,
                                                                 quad_order,
                                                                 row_perm_i2e,
                                                                 col_perm_i2e )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            HINFO( "(TLaplaceSLPGenFn) using AVX512F kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic()  && CFG::BEM::use_simd_mic  )
        {
            HINFO( "(TLaplaceSLPGenFn) using MIC kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            HINFO( "(TLaplaceSLPGenFn) using AVX2 kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx()  && CFG::BEM::use_simd_avx  )
        {
            HINFO( "(TLaplaceSLPGenFn) using AVX kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TLaplaceSLPGenFn) using SSE3 kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TLaplaceSLPGenFn) using VSX kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TLaplaceSLPGenFn) using NEON kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_simd< packed< real, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TLaplaceSLPGenFn) using FPU kernel" );
            
            _eval_dxy_impl = laplace_slp_eval_dx_flt;
        }// else
    }// if
    else
    {
        HINFO( "(TLaplaceSLPGenFn) using FPU kernel" );
            
        _eval_dxy_impl = laplace_slp_eval_dx_flt;
    }// else
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// TLaplaceDLPGenFn
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// constructor
//
template < typename T_ansatzsp, typename T_testsp >
TLaplaceDLPGenFn< T_ansatzsp, T_testsp >::TLaplaceDLPGenFn ( const ansatzsp_t *    ansatzsp,
                                                             const testsp_t *      testsp,
                                                             const TPermutation *  row_perm_i2e,
                                                             const TPermutation *  col_perm_i2e,
                                                             const uint            quad_order )
        : TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, real >( ansatzsp,
                                                                 testsp,
                                                                 quad_order,
                                                                 row_perm_i2e,
                                                                 col_perm_i2e )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            HINFO( "(TLaplaceDLPGenFn) using AVX512F kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real, ISA_AVX512F > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            HINFO( "(TLaplaceDLPGenFn) using MIC kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real, ISA_MIC > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            HINFO( "(TLaplaceDLPGenFn) using AVX2 kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real, ISA_AVX2 > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            HINFO( "(TLaplaceDLPGenFn) using AVX kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real, ISA_AVX > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TLaplaceDLPGenFn) using SSE3 kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real, ISA_SSE3 > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TLaplaceDLPGenFn) using VSX kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real, ISA_VSX > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TLaplaceDLPGenFn) using NEON kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_simd< packed< real, ISA_NEON > >;
            _eval_dy_impl = laplace_dlp_eval_dy_simd< packed< real, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TLaplaceDLPGenFn) using standard kernel" );
            
            _eval_dx_impl = laplace_slp_eval_dx_flt;
            _eval_dy_impl = laplace_dlp_eval_dy_flt;
        }// else
    }// if
    else
    {
        HINFO( "(TLaplaceDLPGenFn) using standard kernel" );
            
        _eval_dx_impl = laplace_slp_eval_dx_flt;
        _eval_dy_impl = laplace_dlp_eval_dy_flt;
    }// if
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template class TLaplaceSLPBF< TConstFnSpace,  TConstFnSpace >;
template class TLaplaceSLPBF< TConstFnSpace,  TLinearFnSpace >;
template class TLaplaceSLPBF< TLinearFnSpace, TConstFnSpace >;
template class TLaplaceSLPBF< TLinearFnSpace, TLinearFnSpace >;

template class TLaplaceDLPBF< TConstFnSpace,  TConstFnSpace >;
template class TLaplaceDLPBF< TConstFnSpace,  TLinearFnSpace >;
template class TLaplaceDLPBF< TLinearFnSpace, TConstFnSpace >;
template class TLaplaceDLPBF< TLinearFnSpace, TLinearFnSpace >;

template class TLaplaceSLPGenFn< TConstFnSpace,  TConstFnSpace >;
template class TLaplaceSLPGenFn< TConstFnSpace,  TLinearFnSpace >;
template class TLaplaceSLPGenFn< TLinearFnSpace, TConstFnSpace >;
template class TLaplaceSLPGenFn< TLinearFnSpace, TLinearFnSpace >;

template class TLaplaceDLPGenFn< TConstFnSpace,  TConstFnSpace >;
template class TLaplaceDLPGenFn< TConstFnSpace,  TLinearFnSpace >;
template class TLaplaceDLPGenFn< TLinearFnSpace, TConstFnSpace >;
template class TLaplaceDLPGenFn< TLinearFnSpace, TLinearFnSpace >;

}// namespace HLIB
