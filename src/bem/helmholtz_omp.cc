//
// Project     : HLib
// File        : helmholtz_omp.cc
// Description : Helmholtz kernels using OpenMP SIMD functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/bem/THelmholtzBF.hh"

namespace HLIB
{

namespace
{

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// local constants
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

const real  ONE_OVER_4PI = real(1) / (real(4) * Math::pi< real >());

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// auxiliary functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#pragma omp declare simd
double
norm ( const double x0,
       const double x1,
       const double x2 )
{
    return std::sqrt( x0*x0 + x1*x1 + x2 * x2 );
}

#pragma omp declare simd
double
normsqr ( const double x0,
          const double x1,
          const double x2 )
{
    return x0*x0 + x1*x1 + x2 * x2;
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Helmholtz SLP kernel
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#pragma omp declare simd uniform(x1,y1,x2,y2,v00,v01,v02,v10,v11,v12,ikappa_re,ikappa_im) linear(i:1)
void
helmholtz_slp_kernel ( const size_t    i,
                       const real *    x1,
                       const real *    y1,
                       const real *    x2,
                       const real *    y2,
                       const double *  v00,
                       const double *  v01,
                       const double *  v02,
                       const double *  v10,
                       const double *  v11,
                       const double *  v12,
                       const real      ikappa_re,
                       const real      ikappa_im,
                       complex *       values )
{
    const real  a1 = x1[i];
    const real  a2 = y1[i];
    const real  a0 = 1.0 - a1 - a2;
        
    const real  b1 = x2[i];
    const real  b2 = y2[i];
    const real  b0 = 1.0 - b1 - b2;

    // u = x - y
    const real  X0 = a0 * v00[0] + a1 * v01[0] + a2 * v02[0];
    const real  X1 = a0 * v00[1] + a1 * v01[1] + a2 * v02[1];
    const real  X2 = a0 * v00[2] + a1 * v01[2] + a2 * v02[2];
    
    const real  Y0 = b0 * v10[0] + b1 * v11[0] + b2 * v12[0];
    const real  Y1 = b0 * v10[1] + b1 * v11[1] + b2 * v12[1];
    const real  Y2 = b0 * v10[2] + b1 * v11[2] + b2 * v12[2];

    const real  u0 = X0 - Y0;
    const real  u1 = X1 - Y1;
    const real  u2 = X2 - Y2;

    //
    //   ( i·κ·|x-y| )
    // e
    // ───────────────
    //      |x-y|
    //

    //
    // compute exp( i·κ·|x-y| )
    //         = polar( exp( re ), im )
    //         = ( exp(re) · cos( im ), exp(re) · sin( im ) )
    //

    const real   dist       = norm( u0, u1, u2 );
    
    // exp( re( i·κ·d ) ); (also: multiply with 1/(4·π) / |x-y| as this would be done later anyway)
    const real   exp_ikd_re = ( ONE_OVER_4PI / dist ) * std::exp( ikappa_re * dist ); 

    // im( i·κ·d )
    const real   ikd_im     = ikappa_im * dist; 

    // cos( i·κ·d ), sin( i·κ·d )
    real         c, s;

    #if HAS_SINCOS == 1
        
    Math::sincos( ikd_im, s, c );
        
    #else
    
    s = std::sin( ikd_im );
    c = std::cos( ikd_im );
    
    #endif

    // exp(re)·cos(), exp(re)·sin()
    c = exp_ikd_re * c;
    s = exp_ikd_re * s;

    // combine real and imaginary parts
    values[ i ] = { real(c), real(s) };
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Helmholtz DLP kernel
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#pragma omp declare simd uniform(v00,v01,v02,v10,v11,v12,ikappa_re,ikappa_im,n) linear(i:1)
void
helmholtz_dlp_kernel ( size_t          i,
                       const real *    x1,
                       const real *    y1,
                       const real *    x2,
                       const real *    y2,
                       const double *  v00,
                       const double *  v01,
                       const double *  v02,
                       const double *  v10,
                       const double *  v11,
                       const double *  v12,
                       const double    ikappa_re,
                       const double    ikappa_im,
                       const double *  n,
                       complex *       values )
{
    const real  a1 = x1[i];
    const real  a2 = y1[i];
    const real  a0 = 1.0 - a1 - a2;
        
    const real  b1 = x2[i];
    const real  b2 = y2[i];
    const real  b0 = 1.0 - b1 - b2;

    // u = x - y
    const real  X0 = a0 * v00[0] + a1 * v01[0] + a2 * v02[0];
    const real  X1 = a0 * v00[1] + a1 * v01[1] + a2 * v02[1];
    const real  X2 = a0 * v00[2] + a1 * v01[2] + a2 * v02[2];
    
    const real  Y0 = b0 * v10[0] + b1 * v11[0] + b2 * v12[0];
    const real  Y1 = b0 * v10[1] + b1 * v11[1] + b2 * v12[1];
    const real  Y2 = b0 * v10[2] + b1 * v11[2] + b2 * v12[2];

    const real  u0 = X0 - Y0;
    const real  u1 = X1 - Y1;
    const real  u2 = X2 - Y2;

    // |y-x|²
    const real  dist2 = normsqr( u0, u1, u2 );
    // |y-x|
    const real  dist  = std::sqrt( dist2 );

    // |y-x|³
    const real  dist3 = dist * dist2;
        
    // < n, y-x >
    const real  ndu   = normsqr( n[0] - u0, n[1] - u1, n[2] - u2 );
    
    //
    //   (i·κ·|x-y|)
    // e            (i·κ·|x-y| - 1) <n,y-x>
    // ────────────────────────────────────
    //                |x-y|³
    //

    //
    // compute exp( i·κ·|x-y| )
    //         = polar( exp( re ), im )
    //         = ( exp(re) · cos( im ), exp(re) · sin( im ) )
    //

    // re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
    real        ikd_re = ikappa_re * dist;
        
    // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
    const real  ikd_im = ikappa_im * dist;
        
    // exp( re(·) )
    // (also: multiply with 1/(4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
    const real  exp_ikd_re = ( ONE_OVER_4PI * ndu / dist3 ) * std::exp( ikd_re ); 

    // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
    real  c, s;

    #if HAS_SINCOS == 1
        
    Math::sincos( ikd_im, s, c );
        
    #else
    
    s = std::sin( ikd_im );
    c = std::cos( ikd_im );
    
    #endif

    // exp(re)·cos(), exp(re)·sin()   (here we would multiply with 1/(4·π) <n,y-x> / |x-y|³)
    c = exp_ikd_re * c;
    s = exp_ikd_re * s;

    // i·κ·|x-y| - 1
    ikd_re = ikd_re - 1.0;
        
    // compute exp() · (i·κ·|x-y| - 1) and store final values
    values[i] = { real( c * ikd_re - s * ikd_im ),
                  real( c * ikd_im - s * ikd_re ) };
}

}// namespace anonymous

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// public functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_ansatzsp,
           typename T_testsp >
void
helmholtz_slp_omp ( const TGrid::triangle_t &   tri0,
                    const TGrid::triangle_t &   tri1,
                    const tripair_quad_rule_t * rule,
                    const complex               ikappa,
                    const T_ansatzsp *          ansatz_sp,
                    const T_testsp *            test_sp,
                    std::vector< complex > &    values )
{
    const size_t   n_pts = rule->npts;
    const T3Point  v00( ansatz_sp->grid()->vertex( tri0.vtx[0] ) );
    const T3Point  v01( ansatz_sp->grid()->vertex( tri0.vtx[1] ) );
    const T3Point  v02( ansatz_sp->grid()->vertex( tri0.vtx[2] ) );
    const T3Point  v10( test_sp->grid()->vertex( tri1.vtx[0] ) );
    const T3Point  v11( test_sp->grid()->vertex( tri1.vtx[1] ) );
    const T3Point  v12( test_sp->grid()->vertex( tri1.vtx[2] ) );

    const auto *   cx1 = rule->x1.data();
    const auto *   cy1 = rule->y1.data();
    const auto *   cx2 = rule->x2.data();
    const auto *   cy2 = rule->y2.data();

    const double *  cv00 = v00.vector();
    const double *  cv01 = v01.vector();
    const double *  cv02 = v02.vector();
    const double *  cv10 = v10.vector();
    const double *  cv11 = v11.vector();
    const double *  cv12 = v12.vector();

    complex *       cvalues = values.data();
        
    #pragma omp simd aligned(cx1,cy1,cx2,cy2:64) linear(cx1,cy1,cx2,cy2:1)
    for ( size_t  i = 0; i < n_pts; ++i )
    {
        helmholtz_slp_kernel( i, cx1, cy1, cx2, cy2,
                              cv00, cv01, cv02, cv10, cv11, cv12,
                              std::real( ikappa ), std::imag( ikappa ), cvalues );
    }// for
}

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
                    std::vector< complex > &    values )
{
    //
    // get normal direction
    //
    
    const T3Point  sn( test_sp->grid()->tri_normal( tri1_id ) );
    
    const size_t   n_pts = rule->npts;
    const T3Point  v00( ansatz_sp->grid()->vertex( tri0.vtx[0] ) );
    const T3Point  v01( ansatz_sp->grid()->vertex( tri0.vtx[1] ) );
    const T3Point  v02( ansatz_sp->grid()->vertex( tri0.vtx[2] ) );
    const T3Point  v10( test_sp->grid()->vertex( tri1.vtx[0] ) );
    const T3Point  v11( test_sp->grid()->vertex( tri1.vtx[1] ) );
    const T3Point  v12( test_sp->grid()->vertex( tri1.vtx[2] ) );

    const auto *   cx1 = & rule->x1[0];
    const auto *   cy1 = & rule->y1[0];
    const auto *   cx2 = & rule->x2[0];
    const auto *   cy2 = & rule->y2[0];

    const double *  cv00 = v00.vector();
    const double *  cv01 = v01.vector();
    const double *  cv02 = v02.vector();
    const double *  cv10 = v10.vector();
    const double *  cv11 = v11.vector();
    const double *  cv12 = v12.vector();

    const double *  cn   = sn.vector();
    
    complex *       cvalues = values.data();
        
    #pragma omp simd aligned(cx1,cy1,cx2,cy2:64) linear(cx1,cy1,cx2,cy2:1)
    for ( size_t  i = 0; i < n_pts; ++i )
    {
        helmholtz_dlp_kernel( i, cx1, cy1, cx2, cy2,
                              cv00, cv01, cv02, cv10, cv11, cv12,
                              std::real( ikappa ), std::imag( ikappa ), cn, cvalues );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// explicit template instantiations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#define  INST_HELMHOLTZ_SLP_OMP( T_ansatzsp, T_testsp )                 \
    template void                                                       \
    helmholtz_slp_omp< T_ansatzsp, T_testsp > (                         \
        const TGrid::triangle_t &   tri0,                               \
        const TGrid::triangle_t &   tri1,                               \
        const tripair_quad_rule_t * rule,                               \
        const complex               ikappa,                             \
        const T_ansatzsp *          ansatz_sp,                          \
        const T_testsp *            test_sp,                            \
        std::vector< complex > &    values );

INST_HELMHOLTZ_SLP_OMP( TConstFnSpace,  TConstFnSpace )
INST_HELMHOLTZ_SLP_OMP( TLinearFnSpace, TConstFnSpace )
INST_HELMHOLTZ_SLP_OMP( TConstFnSpace,  TLinearFnSpace )
INST_HELMHOLTZ_SLP_OMP( TLinearFnSpace, TLinearFnSpace )

#define  INST_HELMHOLTZ_DLP_OMP( T_ansatzsp, T_testsp )                 \
    template void                                                       \
    helmholtz_dlp_omp< T_ansatzsp, T_testsp > (                         \
        const idx_t                 tri1_id,                            \
        const TGrid::triangle_t &   tri0,                               \
        const TGrid::triangle_t &   tri1,                               \
        const tripair_quad_rule_t * rule,                               \
        const complex               ikappa,                             \
        const T_ansatzsp *          ansatz_sp,                          \
        const T_testsp *            test_sp,                            \
        std::vector< complex > &    values );

INST_HELMHOLTZ_DLP_OMP( TConstFnSpace,  TConstFnSpace )
INST_HELMHOLTZ_DLP_OMP( TLinearFnSpace, TConstFnSpace )
INST_HELMHOLTZ_DLP_OMP( TConstFnSpace,  TLinearFnSpace )
INST_HELMHOLTZ_DLP_OMP( TLinearFnSpace, TLinearFnSpace )

}// namespace HLIB
