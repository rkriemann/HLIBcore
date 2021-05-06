//
// Project     : HLib
// File        : helmholtz_simd.inc
// Description : Helmholtz kernels using SIMD functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#if !defined(SIMD_ISA)
#  error "no SIMD ISA defined"
#endif

#include <vector>

#include "hpro/base/packed.hh"

#include "hpro/bem/TConstEdgeFnSpace.hh"
#include "hpro/bem/THelmholtzBF.hh"

namespace HLIB
{

using std::vector;

namespace
{

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// local constants
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

constexpr real  ONE_OVER_4PI = real(1) / (real(4) * Math::pi< real >());

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// local functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// form local quadrature vectors x and y and compute  y-x (Attention: y-x !!!)
//
template < typename T_packed >
inline
void
comp_diff_simd ( const real *    x1,
                 const real *    y1,
                 const real *    x2,
                 const real *    y2,
                 const T_packed  tri0coo[3][3],
                 const T_packed  tri1coo[3][3],
                 T_packed        diff[3] )
{
    using vpacked = T_packed;
    using value_t = typename vpacked::value_t;

    const vpacked  vONE( value_t(1) );

    // load quadrature points for first triangle
    const vpacked  a1 = load< vpacked >( x1 );
    const vpacked  a2 = load< vpacked >( y1 );
    const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );

    // load quadrature points for second triangle
    const vpacked  b1 = load< vpacked >( x2 );
    const vpacked  b2 = load< vpacked >( y2 );
    const vpacked  b0 = sub( sub( vONE, b1 ), b2 );

    //
    // compute y-x
    //

    vpacked  tmp1;
    
    diff[0] = sub( muladd( b2, tri1coo[0][2], muladd( b1, tri1coo[0][1], mul( b0, tri1coo[0][0] ) ) ),
                   muladd( a2, tri0coo[0][2], muladd( a1, tri0coo[0][1], mul( a0, tri0coo[0][0] ) ) ) );
                                                                                                     
    diff[1] = sub( muladd( b2, tri1coo[1][2], muladd( b1, tri1coo[1][1], mul( b0, tri1coo[1][0] ) ) ),
                   muladd( a2, tri0coo[1][2], muladd( a1, tri0coo[1][1], mul( a0, tri0coo[1][0] ) ) ) );
                                                                                                     
    diff[2] = sub( muladd( b2, tri1coo[2][2], muladd( b1, tri1coo[2][1], mul( b0, tri1coo[2][0] ) ) ),
                   muladd( a2, tri0coo[2][2], muladd( a1, tri0coo[2][1], mul( a0, tri0coo[2][0] ) ) ) );
}

//
// form local quadrature vectors x and y and compute  |x-y|²
//
template < typename T_packed >
inline
T_packed
comp_sqdist_simd ( const real *    x1,
                   const real *    y1,
                   const real *    x2,
                   const real *    y2,
                   const T_packed  tri0coo[3][3],
                   const T_packed  tri1coo[3][3] )
{
    using vpacked = T_packed;
    using value_t = typename vpacked::value_t;

    const vpacked  vONE( value_t(1) );

    // load quadrature points for first triangle
    const vpacked  a1 = load< vpacked >( x1 );
    const vpacked  a2 = load< vpacked >( y1 );
    const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );

    // load quadrature points for second triangle
    const vpacked  b1 = load< vpacked >( x2 );
    const vpacked  b2 = load< vpacked >( y2 );
    const vpacked  b0 = sub( sub( vONE, b1 ), b2 );

    //
    // compute |y-x|²
    //

    vpacked  dot, tmp1;
    
    // u0-v0
    tmp1 = sub( muladd( b2, tri1coo[0][2], muladd( b1, tri1coo[0][1], mul( b0, tri1coo[0][0] ) ) ),
                muladd( a2, tri0coo[0][2], muladd( a1, tri0coo[0][1], mul( a0, tri0coo[0][0] ) ) ) );
    dot  = mul( tmp1, tmp1 );                                                              
                                                                                                     
    // u1-v1                                                                               
    tmp1 = sub( muladd( b2, tri1coo[1][2], muladd( b1, tri1coo[1][1], mul( b0, tri1coo[1][0] ) ) ),
                muladd( a2, tri0coo[1][2], muladd( a1, tri0coo[1][1], mul( a0, tri0coo[1][0] ) ) ) );
    dot  = muladd( tmp1, tmp1, dot )   ;                                                   
                                                                                                     
    // u2-v2                                                                               
    tmp1 = sub( muladd( b2, tri1coo[2][2], muladd( b1, tri1coo[2][1], mul( b0, tri1coo[2][0] ) ) ),
                muladd( a2, tri0coo[2][2], muladd( a1, tri0coo[2][1], mul( a0, tri0coo[2][0] ) ) ) );
    dot  = muladd( tmp1, tmp1, dot );

    return dot;
}

}// namespace anonymous

///////////////////////////////////////////////////////////////
//
// some standard definitions used in all functions
//
///////////////////////////////////////////////////////////////

#define DEFINE_LOCAL_TYPE                               \
    using vpacked = T_packed;                           \
    using value_t = typename vpacked::value_t;          \
    const size_t   VECTOR_SIZE = vpacked::vector_size

#define LOAD_TRIANGLE_COORD                                             \
    const vpacked  t0[3][3] = { \
        { ansatz_sp->grid()->vertex( tri0.vtx[0] )[0],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[1] )[0],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[2] )[0] },                \
        { ansatz_sp->grid()->vertex( tri0.vtx[0] )[1],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[1] )[1],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[2] )[1] },                \
        { ansatz_sp->grid()->vertex( tri0.vtx[0] )[2],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[1] )[2],                  \
          ansatz_sp->grid()->vertex( tri0.vtx[2] )[2] } };              \
                                                                        \
    const vpacked  t1[3][3] = { \
        { test_sp->grid()->vertex( tri1.vtx[0] )[0],                    \
          test_sp->grid()->vertex( tri1.vtx[1] )[0],                    \
          test_sp->grid()->vertex( tri1.vtx[2] )[0] },                  \
        { test_sp->grid()->vertex( tri1.vtx[0] )[1],                    \
          test_sp->grid()->vertex( tri1.vtx[1] )[1],                    \
          test_sp->grid()->vertex( tri1.vtx[2] )[1] },                  \
        { test_sp->grid()->vertex( tri1.vtx[0] )[2],                    \
          test_sp->grid()->vertex( tri1.vtx[1] )[2],                    \
          test_sp->grid()->vertex( tri1.vtx[2] )[2] } }

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Helmholtz SLP kernel
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
helmholtz_slp_simd ( const TGrid::triangle_t &   tri0,
                     const TGrid::triangle_t &   tri1,
                     const tripair_quad_rule_t * rule,
                     const complex               ikappa,
                     const T_ansatzsp *          ansatz_sp,
                     const T_testsp *            test_sp,
                     vector< complex > &         values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;
    
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    const vpacked  vikappa_re( ikappa.real() );
    const vpacked  vikappa_im( ikappa.imag() );
    const size_t   n_pts = rule->npts;

    //
    // eval VECTOR_SIZE add the same time
    //

    value_t  res_re[ VECTOR_SIZE ];
    value_t  res_im[ VECTOR_SIZE ];

    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute and store |x-y|
        //
        
        const vpacked  vdist = sqrt( comp_sqdist_simd< vpacked >( & rule->x1[i], & rule->y1[i],
                                                                  & rule->x2[i], & rule->y2[i],
                                                                  t0, t1 ) );

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
        
        // exp( re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) ) )
        // (also: multiply with 1/(4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( vONE_OVER_4PI, vdist ), exp( mul( vikappa_re, vdist ) ) ); 

        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist ); 

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        vpacked        vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // combine real and imaginary parts
        store( vc, res_re );
        store( vs, res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; j++ )
            values[ i+j ] = complex( res_re[j], res_im[j] );
    }// for
}

///////////////////////////////////////////////////////////////
//
// Helmholtz SLP kernel; real wave number
//

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
helmholtz_slp_re_simd ( const TGrid::triangle_t &   tri0,
                        const TGrid::triangle_t &   tri1,
                        const tripair_quad_rule_t * rule,
                        const complex               ikappa,
                        const T_ansatzsp *          ansatz_sp,
                        const T_testsp *            test_sp,
                        vector< complex > &         values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;
    
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    const size_t   n_pts = rule->npts;
    const vpacked  vikappa_im( ikappa.imag() );

    //
    // eval VECTOR_SIZE add the same time
    //

    value_t  res_re[ VECTOR_SIZE ];
    value_t  res_im[ VECTOR_SIZE ];

    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute and store |x-y|
        //
        
        const vpacked  vdist = sqrt( comp_sqdist_simd< vpacked >( & rule->x1[i], & rule->y1[i],
                                                                  & rule->x2[i], & rule->y2[i],
                                                                  t0, t1 ) );

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
        
        // exp( re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) ) ) = 1 since re = 0
        // (also: multiply with 1/(4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = div( vONE_OVER_4PI, vdist ); 

        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im     = mul( vikappa_im, vdist ); 
        vpacked        vc, vs;

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // combine real and imaginary parts
        store( vc, res_re );
        store( vs, res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; j++ )
            values[ i+j ] = complex( res_re[j], res_im[j] );
    }// for
}

///////////////////////////////////////////////////////////////
//
// Helmholtz SLP kernel; imaginary wave number
//

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
helmholtz_slp_im_simd ( const TGrid::triangle_t &   tri0,
                        const TGrid::triangle_t &   tri1,
                        const tripair_quad_rule_t * rule,
                        const complex               ikappa,
                        const T_ansatzsp *          ansatz_sp,
                        const T_testsp *            test_sp,
                        vector< complex > &         values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;
    
    const vpacked   vONE_OVER_4PI( ONE_OVER_4PI );
    const size_t    n_pts   = rule->npts;
    const vpacked   vikappa_re( ikappa.real() );

    //
    // eval VECTOR_SIZE add the same time
    //

    value_t  res_re[ VECTOR_SIZE ];
    
    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute and store |x-y|
        //
        
        const vpacked  vdist = sqrt( comp_sqdist_simd< vpacked >( & rule->x1[i], & rule->y1[i],
                                                                  & rule->x2[i], & rule->y2[i],
                                                                  t0, t1 ) );

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
        
        // exp( re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) ) )
        // (also: multiply with 1/(4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( vONE_OVER_4PI, vdist ), exp( mul( vikappa_re, vdist ) ) ); 

        // exp(re)·cos(im), exp(re)·sin(im) = ( exp(re), 0 ) since im = 0
        // combine real and imaginary parts
        store( vexp_ikd_re, res_re );

        for ( size_t  j = 0; j < VECTOR_SIZE; j++ )
            values[ i+j ] = complex( res_re[j], real(0) );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Helmholtz DLP kernel
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
helmholtz_dlp_simd ( const idx_t                 tri1_id,
                     const TGrid::triangle_t &   tri0,
                     const TGrid::triangle_t &   tri1,
                     const tripair_quad_rule_t * rule,
                     const complex               ikappa,
                     const T_ansatzsp *          ansatz_sp,
                     const T_testsp *            test_sp,
                     vector< complex > &         values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;

    const vpacked  vONE( value_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    //
    // get normal direction
    //
    
    const T3Point  sn( test_sp->grid()->tri_normal( tri1_id ) );
    const vpacked  n[3] = { sn[0], sn[1], sn[2] };
    
    //
    // eval VECTOR_SIZE  add the same time
    //

    const vpacked  vikappa_re( ikappa.real() );
    const vpacked  vikappa_im( ikappa.imag() );
    const size_t   n_pts = rule->npts;
    value_t        res_re[ VECTOR_SIZE ];
    value_t        res_im[ VECTOR_SIZE ];
    
    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute y-x
        //
        
        vpacked  diff[3];
        
        comp_diff_simd( & rule->x1[i], & rule->y1[i],
                        & rule->x2[i], & rule->y2[i],
                        t0, t1, diff );

        // |y-x|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |y-x|
        const vpacked  vdist  = sqrt( vdist2 );

        // |y-x|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // < n, y-x >
        const vpacked  vndu   = muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) );

        
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
        vpacked        vikd_re = mul( vikappa_re, vdist );
        
        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist );
        
        // exp( re(·) )
        // (also: multiply with 1/(4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( mul( vONE_OVER_4PI, vndu ), vdist3 ), exp( vikd_re ) ); 

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        vpacked  vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()   (here we would multiply with 1/(4·π) <n,y-x> / |x-y|³)
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // i·κ·|x-y| - 1
        vikd_re = sub( vikd_re, vONE );
        
        // compute exp() · (i·κ·|x-y| - 1) and store final values
        store( mulsub( vc, vikd_re, mul( vs, vikd_im ) ), res_re );
        store( muladd( vc, vikd_im, mul( vs, vikd_re ) ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ i+j ] = complex( res_re[j], res_im[j] );
    }// for
}

///////////////////////////////////////////////////////////////
//
// Helmholtz DLP kernel; real wavenumber
//

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
helmholtz_dlp_re_simd ( const idx_t                 tri1_id,
                        const TGrid::triangle_t &   tri0,
                        const TGrid::triangle_t &   tri1,
                        const tripair_quad_rule_t * rule,
                        const complex               ikappa,
                        const T_ansatzsp *          ansatz_sp,
                        const T_testsp *            test_sp,
                        vector< complex > &         values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;

    const vpacked  vMINUS_ONE( value_t(-1) );  // used for: re(i·κ·|x-y|)-1 since re() = 0
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    //
    // get normal direction
    //
    
    const T3Point  sn( test_sp->grid()->tri_normal( tri1_id ) );
    const vpacked  n[3] = { sn[0], sn[1], sn[2] };
    
    //
    // eval VECTOR_SIZE  add the same time
    //

    const vpacked  vikappa_im( ikappa.imag() );
    const size_t   n_pts = rule->npts;
    value_t        res_re[ VECTOR_SIZE ];
    value_t        res_im[ VECTOR_SIZE ];
    
    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute y-x
        //
        
        vpacked  diff[3];
        
        comp_diff_simd( & rule->x1[i], & rule->y1[i],
                        & rule->x2[i], & rule->y2[i],
                        t0, t1, diff );

        // |y-x|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |y-x|
        const vpacked  vdist  = sqrt( vdist2 );

        // |y-x|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // < n, y-x >
        const vpacked  vndu   = muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) );

        
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

        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist );
        
        // exp( re(·) ) = 1 since re() = 0
        // (also: multiply with 1/(4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = div( mul( vONE_OVER_4PI, vndu ), vdist3 ); 

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        vpacked  vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()   (here we would multiply with 1/(4·π) <n,y-x> / |x-y|³)
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // compute exp() · (i·κ·|x-y| - 1) and store final values
        store( mulsub( vc, vMINUS_ONE, mul( vs, vikd_im    ) ), res_re );
        store( muladd( vc, vikd_im,    mul( vs, vMINUS_ONE ) ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ i+j ] = complex( res_re[j], res_im[j] );
    }// for
}

///////////////////////////////////////////////////////////////
//
// Helmholtz DLP kernel; imaginary wavenumber
//

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
helmholtz_dlp_im_simd ( const idx_t                 tri1_id,
                        const TGrid::triangle_t &   tri0,
                        const TGrid::triangle_t &   tri1,
                        const tripair_quad_rule_t * rule,
                        const complex               ikappa,
                        const T_ansatzsp *          ansatz_sp,
                        const T_testsp *            test_sp,
                        vector< complex > &         values )
{
    DEFINE_LOCAL_TYPE;
    LOAD_TRIANGLE_COORD;

    const vpacked  vONE( value_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    //
    // get normal direction
    //
    
    const T3Point  sn( test_sp->grid()->tri_normal( tri1_id ) );
    const vpacked  n[3] = { sn[0], sn[1], sn[2] };
    
    //
    // eval VECTOR_SIZE  add the same time
    //

    const vpacked  vikappa_re( ikappa.real() );
    const size_t   n_pts = rule->npts;
    value_t        res_re[ VECTOR_SIZE ];
    
    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute y-x
        //
        
        vpacked  diff[3];
        
        comp_diff_simd( & rule->x1[i], & rule->y1[i],
                        & rule->x2[i], & rule->y2[i],
                        t0, t1, diff );

        // |y-x|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |y-x|
        const vpacked  vdist  = sqrt( vdist2 );

        // |y-x|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // < n, y-x >
        const vpacked  vndu   = muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) );

        
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
        vpacked        vikd_re = mul( vikappa_re, vdist );
        
        // exp( re(·) )
        // (also: multiply with 1/(4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( mul( vONE_OVER_4PI, vndu ), vdist3 ), exp( vikd_re ) ); 

        // exp(re)·cos(), exp(re)·sin() = ( exp(re), 0 ) since im = 0

        // i·κ·|x-y| - 1
        vikd_re = sub( vikd_re, vONE );
        
        // compute exp() · (i·κ·|x-y| - 1) ( only real numbers ! )
        vikd_re = mul( vexp_ikd_re, vikd_re );
        
        // store final values
        store( vikd_re, res_re );
        
        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ i+j ] = complex( res_re[j], real(0) );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Helmholtz HCA functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_packed >
void
helmholtz_slp_eval_dx_simd ( const complex &          ikappa,
                             const tri_quad_rule_t &  quad_rule,
                             const T3Point            vx[3],
                             const T3Point &          vy,
                             vector< complex > &      values )
{
    DEFINE_LOCAL_TYPE;

    const vpacked  vONE( value_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_re( ikappa.real() );
    const vpacked  vikappa_im( ikappa.imag() );
    const size_t   npts  = quad_rule.npts;
    const vpacked  y[3]  = { vy.x(),    vy.y(),    vy.z() };
    const vpacked  x0[3] = { vx[0].x(), vx[0].y(), vx[0].z() };
    const vpacked  x1[3] = { vx[1].x(), vx[1].y(), vx[1].z() };
    const vpacked  x2[3] = { vx[2].x(), vx[2].y(), vx[2].z() };
    value_t        res_re[ VECTOR_SIZE ];
    value_t        res_im[ VECTOR_SIZE ];

    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>
        //

        vpacked  diff[3];

        diff[0] = sub( y[0], muladd( a0, x0[0], muladd( a1, x1[0], mul( a2, x2[0] ) ) ) );
        diff[1] = sub( y[1], muladd( a0, x0[1], muladd( a1, x1[1], mul( a2, x2[1] ) ) ) );
        diff[2] = sub( y[2], muladd( a0, x0[2], muladd( a1, x1[2], mul( a2, x2[2] ) ) ) );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        //
        //   ( i·κ·|x-y| )
        // e
        // ───────────────
        //      |x-y|
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // exp( re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) ) )
        // (also: multiply with (1/4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( vONE_OVER_4PI, vdist ), exp( mul( vikappa_re, vdist ) ) );

        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist ); 
        vpacked        vc, vs;

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()
        // combine real and imaginary parts
        store( mul( vexp_ikd_re, vc ), res_re );
        store( mul( vexp_ikd_re, vs ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = complex( res_re[j], res_im[j] );
    }// for
}

template < typename T_packed >
void
helmholtz_slp_eval_dx_re_simd ( const complex &          ikappa,
                                const tri_quad_rule_t &  quad_rule,
                                const T3Point            vx[3],
                                const T3Point &          vy,
                                vector< complex > &      values )
{
    DEFINE_LOCAL_TYPE;

    const vpacked  vONE( value_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_im( ikappa.imag() );  // κ real -> i·κ imaginary
    const size_t   npts  = quad_rule.npts;
    const vpacked  y[3]  = { vy.x(),    vy.y(),    vy.z() };
    const vpacked  x0[3] = { vx[0].x(), vx[0].y(), vx[0].z() };
    const vpacked  x1[3] = { vx[1].x(), vx[1].y(), vx[1].z() };
    const vpacked  x2[3] = { vx[2].x(), vx[2].y(), vx[2].z() };
    value_t        res_re[ VECTOR_SIZE ];
    value_t        res_im[ VECTOR_SIZE ];
    
    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>
        //

        vpacked  diff[3];

        diff[0] = sub( y[0], muladd( a0, x0[0], muladd( a1, x1[0], mul( a2, x2[0] ) ) ) );
        diff[1] = sub( y[1], muladd( a0, x0[1], muladd( a1, x1[1], mul( a2, x2[1] ) ) ) );
        diff[2] = sub( y[2], muladd( a0, x0[2], muladd( a1, x1[2], mul( a2, x2[2] ) ) ) );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        //
        //   ( i·κ·|x-y| )
        // e
        // ───────────────
        //      |x-y|
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // exp( re( ( i·κ·d[0], i·κ·d[1] ) ) ) = 1 since re() = 0
        // (also: multiply with 1/(4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = div( vONE_OVER_4PI, vdist ); 

        // im( ( i·κ·d[0], i·κ·d[1] ) )
        const vpacked  vikd_im     = mul( vikappa_im, vdist ); 

        // cos( i·κ·d[0..1] ), sin( i·κ·d[0..1] )
        vpacked        vc, vs;
        
        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()
        // combine real and imaginary parts
        store( mul( vexp_ikd_re, vc ), res_re );
        store( mul( vexp_ikd_re, vs ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = complex( res_re[j], res_im[j] );
    }// for
}

template < typename T_packed >
void
helmholtz_slp_eval_dx_im_simd ( const complex &         ikappa,
                               const tri_quad_rule_t &  quad_rule,
                               const T3Point            vx[3],
                               const T3Point &          vy,
                               vector< complex > &      values )
{
    DEFINE_LOCAL_TYPE;

    const vpacked  vONE( value_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_re( ikappa.real() ); // κ imaginary -> i·κ real
    const size_t   npts  = quad_rule.npts;
    const vpacked  y[3]  = { vy.x(),    vy.y(),    vy.z() };
    const vpacked  x0[3] = { vx[0].x(), vx[0].y(), vx[0].z() };
    const vpacked  x1[3] = { vx[1].x(), vx[1].y(), vx[1].z() };
    const vpacked  x2[3] = { vx[2].x(), vx[2].y(), vx[2].z() };
    value_t        res_re[ VECTOR_SIZE ];
    
    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>
        //

        vpacked  diff[3];

        diff[0] = sub( y[0], muladd( a0, x0[0], muladd( a1, x1[0], mul( a2, x2[0] ) ) ) );
        diff[1] = sub( y[1], muladd( a0, x0[1], muladd( a1, x1[1], mul( a2, x2[1] ) ) ) );
        diff[2] = sub( y[2], muladd( a0, x0[2], muladd( a1, x1[2], mul( a2, x2[2] ) ) ) );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        //
        //   ( i·κ·|x-y| )
        // e
        // ───────────────
        //      |x-y|
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // exp( re( ( i·κ·d[0], i·κ·d[1] ) ) )
        // (also: multiply with 1/(4·π) / |x-y| as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( vONE_OVER_4PI, vdist ), exp( mul( vikappa_re, vdist ) ) ); 

        // exp(re)·cos(im), exp(re)·sin(im) = ( exp(re), 0 ) since im = 0
        // (here we would multiply with 1/(4·π) / |x-y|)

        // combine real and imaginary parts
        store( vexp_ikd_re, res_re );
        
        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = complex( res_re[j], real(0) );
    }// for
}

template < typename T_packed >
void
helmholtz_dlp_eval_dy_simd ( const complex &          ikappa,
                             const tri_quad_rule_t &  quad_rule,
                             const T3Point &          vx,
                             const T3Point            vy[3],
                             const T3Point &          normal,
                             vector< complex > &      values )
{
    DEFINE_LOCAL_TYPE;

    const vpacked  vONE( value_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_re( ikappa.real() );
    const vpacked  vikappa_im( ikappa.imag() );
    const size_t   npts  = quad_rule.npts;
    const vpacked  x[3]  = { vx.x(),     vx.y(),     vx.z() };
    const vpacked  y0[3] = { vy[0].x(),  vy[0].y(),  vy[0].z() };
    const vpacked  y1[3] = { vy[1].x(),  vy[1].y(),  vy[1].z() };
    const vpacked  y2[3] = { vy[2].x(),  vy[2].y(),  vy[2].z() };
    const vpacked  n[3]  = { normal.x(), normal.y(), normal.z() };
    value_t        res_re[ VECTOR_SIZE ];
    value_t        res_im[ VECTOR_SIZE ];

    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>; <n,d>
        //

        vpacked  diff[3];

        diff[0] = sub( muladd( a0, y0[0], muladd( a1, y1[0], mul( a2, y2[0] ) ) ), x[0] );
        diff[1] = sub( muladd( a0, y0[1], muladd( a1, y1[1], mul( a2, y2[1] ) ) ), x[1] );
        diff[2] = sub( muladd( a0, y0[2], muladd( a1, y1[2], mul( a2, y2[2] ) ) ), x[2] );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        // |x-y|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // < n, y-x >
        const vpacked  vndu   = muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) );

        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // re( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        vpacked        vikd_re = mul( vikappa_re, vdist );
        
        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist );
        
        // exp( re(·) )
        // (also: multiply with (1/4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( mul( vONE_OVER_4PI, vndu ), vdist3 ), exp( vikd_re ) ); 

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        vpacked  vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()   (here we would multiply with <n,y-x> / |x-y|³)
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // i·κ·|x-y| - 1
        vikd_re = sub( vikd_re, vONE );
        
        // compute exp() · (i·κ·|x-y| - 1) and store final values
        store( mulsub( vc, vikd_re, mul( vs, vikd_im ) ), res_re );
        store( muladd( vc, vikd_im, mul( vs, vikd_re ) ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = complex( res_re[j], res_im[j] );
    }// for
}

template < typename T_packed >
void
helmholtz_dlp_eval_dy_re_simd ( const complex &          ikappa,
                                const tri_quad_rule_t &  quad_rule,
                                const T3Point &          vx,
                                const T3Point            vy[3],
                                const T3Point &          normal,
                                vector< complex > &      values )
{
    DEFINE_LOCAL_TYPE;

    const vpacked  vONE( value_t(1) );
    const vpacked  vMINUS_ONE( value_t(-1) );  // used for: re(i·κ·|x-y|)-1 since re() = 0
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_im( ikappa.imag() );  // κ real -> i·κ imaginary
    const size_t   npts  = quad_rule.npts;
    const vpacked  x[3]  = { vx.x(),     vx.y(),     vx.z() };
    const vpacked  y0[3] = { vy[0].x(),  vy[0].y(),  vy[0].z() };
    const vpacked  y1[3] = { vy[1].x(),  vy[1].y(),  vy[1].z() };
    const vpacked  y2[3] = { vy[2].x(),  vy[2].y(),  vy[2].z() };
    const vpacked  n[3]  = { normal.x(), normal.y(), normal.z() };
    value_t        res_re[ VECTOR_SIZE ];
    value_t        res_im[ VECTOR_SIZE ];

    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>; <n,d>
        //

        vpacked  diff[3];

        diff[0] = sub( muladd( a0, y0[0], muladd( a1, y1[0], mul( a2, y2[0] ) ) ), x[0] );
        diff[1] = sub( muladd( a0, y0[1], muladd( a1, y1[1], mul( a2, y2[1] ) ) ), x[1] );
        diff[2] = sub( muladd( a0, y0[2], muladd( a1, y1[2], mul( a2, y2[2] ) ) ), x[2] );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        // |x-y|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // < n, y-x >
        const vpacked  vndu   = muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) );

        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // im( ( i·κ·d[0], i·κ·d[1] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist );
        
        // exp( re(·) ) = 1 since re() = 0
        // (also: multiply with 1/(4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = div( mul( vONE_OVER_4PI, vndu ), vdist3 ); 

        // cos( i·κ·d[0..1] ), sin( i·κ·d[0..1] )
        vpacked  vc, vs;

        sincos( vikd_im, vs, vc );

        // exp(re)·cos(), exp(re)·sin()   (here we would multiply with 1/(4·π) <n,y-x> / |x-y|³)
        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // compute exp() · (i·κ·|x-y| - 1) and store final values
        store( mulsub( vc, vMINUS_ONE, mul( vs, vikd_im    ) ), res_re );
        store( muladd( vc, vikd_im,    mul( vs, vMINUS_ONE ) ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = complex( res_re[j], res_im[j] );
    }// for
}

template < typename T_packed >
void
helmholtz_dlp_eval_dy_im_simd ( const complex &          ikappa,
                                const tri_quad_rule_t &  quad_rule,
                                const T3Point &          vx,
                                const T3Point            vy[3],
                                const T3Point &          normal,
                                vector< complex > &      values )
{
    DEFINE_LOCAL_TYPE;

    const vpacked  vONE( value_t(1) );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    const vpacked  vikappa_re( ikappa.real() ); // κ imaginary -> i·κ real
    const size_t   npts  = quad_rule.npts;
    const vpacked  x[3]  = { vx.x(),     vx.y(),     vx.z() };
    const vpacked  y0[3] = { vy[0].x(),  vy[0].y(),  vy[0].z() };
    const vpacked  y1[3] = { vy[1].x(),  vy[1].y(),  vy[1].z() };
    const vpacked  y2[3] = { vy[2].x(),  vy[2].y(),  vy[2].z() };
    const vpacked  n[3]  = { normal.x(), normal.y(), normal.z() };
    value_t        res_re[ VECTOR_SIZE ];

    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        // load quadrature points and weight for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // compute d = y - x; <d,d>; <n,d>
        //

        vpacked  diff[3];

        diff[0] = sub( muladd( a0, y0[0], muladd( a1, y1[0], mul( a2, y2[0] ) ) ), x[0] );
        diff[1] = sub( muladd( a0, y0[1], muladd( a1, y1[1], mul( a2, y2[1] ) ) ), x[1] );
        diff[2] = sub( muladd( a0, y0[2], muladd( a1, y1[2], mul( a2, y2[2] ) ) ), x[2] );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        // |x-y|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        // < n, y-x >
        const vpacked  vndu   = muladd( n[0], diff[0], muladd( n[1], diff[1], mul( n[2], diff[2] ) ) );

        //
        //   (i·κ·|x-y|)
        // e            (i·κ·|x-y| - 1) <n,y-x>
        // ────────────────────────────────────
        //                |x-y|³
        //
        // compute exp( i·κ·|x-y| ) = polar( exp( re ), im )
        //                          = ( exp(re) · cos( im ), exp(re) · sin( im ) )
        //

        // re( ( i·κ·d[0], i·κ·d[1] ) )
        vpacked        vikd_re = mul( vikappa_re, vdist );
        
        // exp( re(·) )
        // (also: multiply with <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( mul( vONE_OVER_4PI, vndu ), vdist3 ),
                                          exp( vikd_re ) );

        // exp(re)·cos(im), exp(re)·sin(im) = ( exp(re), 0 ) since im = 0 
        // (here we would multiply with 1/(4·π) <n,y-x> / |x-y|³)

        // i·κ·|x-y| - 1
        vikd_re = sub( vikd_re, vONE );
        
        // compute exp() · (i·κ·|x-y| - 1)  ( only real numbers ! )
        vikd_re = mul( vexp_ikd_re, vikd_re );

        // compute exp() · (i·κ·|x-y| - 1) and store final values
        store( vikd_re, res_re );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ k+j ] = complex( res_re[j], real(0) );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// explicit template instantiations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#define  INST_HELMHOLTZ_SLP_SIMD( T_value, T_ansatzsp, T_testsp, suffix ) \
    template void                                                       \
    helmholtz_slp_##suffix< T_ansatzsp, T_testsp, packed< T_value, SIMD_ISA > > ( \
        const TGrid::triangle_t &   tri0,                               \
        const TGrid::triangle_t &   tri1,                               \
        const tripair_quad_rule_t * rule,                               \
        const complex               ikappa,                             \
        const T_ansatzsp *          ansatz_sp,                          \
        const T_testsp *            test_sp,                            \
        vector< complex > &         values )

INST_HELMHOLTZ_SLP_SIMD( real, TConstFnSpace,     TConstFnSpace,     simd );
INST_HELMHOLTZ_SLP_SIMD( real, TLinearFnSpace,    TConstFnSpace,     simd );
INST_HELMHOLTZ_SLP_SIMD( real, TConstFnSpace,     TLinearFnSpace,    simd );
INST_HELMHOLTZ_SLP_SIMD( real, TLinearFnSpace,    TLinearFnSpace,    simd );
INST_HELMHOLTZ_SLP_SIMD( real, TConstEdgeFnSpace, TConstEdgeFnSpace, simd );

INST_HELMHOLTZ_SLP_SIMD( real, TConstFnSpace,     TConstFnSpace,     re_simd );
INST_HELMHOLTZ_SLP_SIMD( real, TLinearFnSpace,    TConstFnSpace,     re_simd );
INST_HELMHOLTZ_SLP_SIMD( real, TConstFnSpace,     TLinearFnSpace,    re_simd );
INST_HELMHOLTZ_SLP_SIMD( real, TLinearFnSpace,    TLinearFnSpace,    re_simd );
INST_HELMHOLTZ_SLP_SIMD( real, TConstEdgeFnSpace, TConstEdgeFnSpace, re_simd );

INST_HELMHOLTZ_SLP_SIMD( real, TConstFnSpace,     TConstFnSpace,     im_simd );
INST_HELMHOLTZ_SLP_SIMD( real, TLinearFnSpace,    TConstFnSpace,     im_simd );
INST_HELMHOLTZ_SLP_SIMD( real, TConstFnSpace,     TLinearFnSpace,    im_simd );
INST_HELMHOLTZ_SLP_SIMD( real, TLinearFnSpace,    TLinearFnSpace,    im_simd );
INST_HELMHOLTZ_SLP_SIMD( real, TConstEdgeFnSpace, TConstEdgeFnSpace, im_simd );

#define  INST_HELMHOLTZ_DLP_SIMD( T_value, T_ansatzsp, T_testsp, suffix ) \
    template void                                                       \
    helmholtz_dlp_##suffix< T_ansatzsp, T_testsp, packed< T_value, SIMD_ISA > > ( \
        const idx_t                 tri1_id,                            \
        const TGrid::triangle_t &   tri0,                               \
        const TGrid::triangle_t &   tri1,                               \
        const tripair_quad_rule_t * rule,                               \
        const complex               ikappa,                             \
        const T_ansatzsp *          ansatz_sp,                          \
        const T_testsp *            test_sp,                            \
        vector< complex > &         values )

INST_HELMHOLTZ_DLP_SIMD( real, TConstFnSpace,  TConstFnSpace,  simd );
INST_HELMHOLTZ_DLP_SIMD( real, TLinearFnSpace, TConstFnSpace,  simd );
INST_HELMHOLTZ_DLP_SIMD( real, TConstFnSpace,  TLinearFnSpace, simd );
INST_HELMHOLTZ_DLP_SIMD( real, TLinearFnSpace, TLinearFnSpace, simd );

INST_HELMHOLTZ_DLP_SIMD( real, TConstFnSpace,  TConstFnSpace,  re_simd );
INST_HELMHOLTZ_DLP_SIMD( real, TLinearFnSpace, TConstFnSpace,  re_simd );
INST_HELMHOLTZ_DLP_SIMD( real, TConstFnSpace,  TLinearFnSpace, re_simd );
INST_HELMHOLTZ_DLP_SIMD( real, TLinearFnSpace, TLinearFnSpace, re_simd );

INST_HELMHOLTZ_DLP_SIMD( real, TConstFnSpace,  TConstFnSpace,  im_simd );
INST_HELMHOLTZ_DLP_SIMD( real, TLinearFnSpace, TConstFnSpace,  im_simd );
INST_HELMHOLTZ_DLP_SIMD( real, TConstFnSpace,  TLinearFnSpace, im_simd );
INST_HELMHOLTZ_DLP_SIMD( real, TLinearFnSpace, TLinearFnSpace, im_simd );


#define  INST_HELMHOLTZ_SLP_HCA_SIMD( T_value, suffix )             \
    template void                                                   \
    helmholtz_slp_eval_dx_##suffix< packed< T_value, SIMD_ISA > > ( \
        const complex &            ikappa,                          \
        const tri_quad_rule_t &    quad_rule,                       \
        const T3Point              vx[3],                           \
        const T3Point &            vy,                              \
        vector< complex > &        values )

INST_HELMHOLTZ_SLP_HCA_SIMD( real, simd );
INST_HELMHOLTZ_SLP_HCA_SIMD( real, re_simd );
INST_HELMHOLTZ_SLP_HCA_SIMD( real, im_simd );


#define  INST_HELMHOLTZ_DLP_HCA_SIMD( T_value, suffix )             \
    template void                                                   \
    helmholtz_dlp_eval_dy_##suffix< packed< T_value, SIMD_ISA > > (  \
        const complex &            ikappa,                          \
        const tri_quad_rule_t &    quad_rule,                       \
        const T3Point &            vx,                              \
        const T3Point              vy[3],                           \
        const T3Point &            normal,                          \
        vector< complex > &        values )

INST_HELMHOLTZ_DLP_HCA_SIMD( real, simd );
INST_HELMHOLTZ_DLP_HCA_SIMD( real, re_simd );
INST_HELMHOLTZ_DLP_HCA_SIMD( real, im_simd );

}// namespace HLIB

// Local Variables:
// mode: c++
// End:
