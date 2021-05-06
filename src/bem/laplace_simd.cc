//
// Project     : HLib
// File        : laplace_simd.cc
// Description : Laplace kernels using SIMD functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#if !defined(SIMD_ISA)
#  error "no SIMD ISA defined"
#endif

#include <vector>

#include "hpro/base/packed.hh"

#include "hpro/bem/TLaplaceBF.hh"

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

const real  ONE_OVER_4PI = real(1) / (real(4) * Math::pi< real >());

}// namespace anonymous

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Laplace SLP kernel
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
laplace_slp_simd ( const TGrid::triangle_t &    tri0,
                   const TGrid::triangle_t &    tri1,
                   const tripair_quad_rule_t *  rule,
                   vector< real > &             values,
                   const T_ansatzsp *           ansatz_sp,
                   const T_testsp *             test_sp )
{
    using vpacked = T_packed;
    using value_t = typename vpacked::value_t;

    const size_t   VECTOR_SIZE = vpacked::vector_size;
    
    const vpacked  vONE( 1.0 );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );

    // load triangle coordinates, all coefficients in t0[i][j]
    // contain coord j (x, y, or z) of vertex i
    const vpacked  t0[3][3] = { { ansatz_sp->grid()->vertex( tri0.vtx[0] )[0],
                                  ansatz_sp->grid()->vertex( tri0.vtx[0] )[1],
                                  ansatz_sp->grid()->vertex( tri0.vtx[0] )[2] },
                                { ansatz_sp->grid()->vertex( tri0.vtx[1] )[0],
                                  ansatz_sp->grid()->vertex( tri0.vtx[1] )[1],
                                  ansatz_sp->grid()->vertex( tri0.vtx[1] )[2] },
                                { ansatz_sp->grid()->vertex( tri0.vtx[2] )[0],
                                  ansatz_sp->grid()->vertex( tri0.vtx[2] )[1],
                                  ansatz_sp->grid()->vertex( tri0.vtx[2] )[2] } };
    
    const vpacked  t1[3][3] = { { test_sp->grid()->vertex( tri1.vtx[0] )[0],
                                  test_sp->grid()->vertex( tri1.vtx[0] )[1],
                                  test_sp->grid()->vertex( tri1.vtx[0] )[2] },
                                { test_sp->grid()->vertex( tri1.vtx[1] )[0],
                                  test_sp->grid()->vertex( tri1.vtx[1] )[1],
                                  test_sp->grid()->vertex( tri1.vtx[1] )[2] },
                                { test_sp->grid()->vertex( tri1.vtx[2] )[0],
                                  test_sp->grid()->vertex( tri1.vtx[2] )[1],
                                  test_sp->grid()->vertex( tri1.vtx[2] )[2] } };

    //
    // eval VECTOR_SIZE add the same time
    //

    const size_t  n_pts = rule->npts;
    value_t       res[ VECTOR_SIZE ];

    #pragma ivdep
    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        //
        // compute 1/|x-y|
        //
        
        // load quadrature points for first triangle
        const vpacked  a1 = load< vpacked >( & rule->x1[i] );  // ( x1[0], x1[1], ... x1[N] )
        const vpacked  a2 = load< vpacked >( & rule->y1[i] );  // ( y1[0], y1[1], ... y1[N] )
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );

        // load quadrature points for second triangle
        const vpacked  b1 = load< vpacked >( & rule->x2[i] );  // ( x2[0], x2[1], ... x2[N] )
        const vpacked  b2 = load< vpacked >( & rule->y2[i] );  // ( y2[0], y2[1], ... y2[N] )
        const vpacked  b0 = sub( sub( vONE, b1 ), b2 );

        //
        // compute |x-y|²
        //

        vpacked  dot( 0 );
        vpacked  tmp1;

        // u0-v0
        tmp1 =            muladd( a0, t0[0][0], muladd( a1, t0[1][0], mul( a2, t0[2][0] ) ) );
        tmp1 = sub( tmp1, muladd( b0, t1[0][0], muladd( b1, t1[1][0], mul( b2, t1[2][0] ) ) ) );
        dot  = mul( tmp1, tmp1 );

        // u1-v1
        tmp1 =            muladd( a0, t0[0][1], muladd( a1, t0[1][1], mul( a2, t0[2][1] ) ) );
        tmp1 = sub( tmp1, muladd( b0, t1[0][1], muladd( b1, t1[1][1], mul( b2, t1[2][1] ) ) ) );
        dot  = muladd( tmp1, tmp1, dot );

        // u2-v2
        tmp1 =            muladd( a0, t0[0][2], muladd( a1, t0[1][2], mul( a2, t0[2][2] ) ) );
        tmp1 = sub( tmp1, muladd( b0, t1[0][2], muladd( b1, t1[1][2], mul( b2, t1[2][2] ) ) ) );
        dot  = muladd( tmp1, tmp1, dot );

        // compute (4π / |x-y|) and store in <values>
        store( mul( vONE_OVER_4PI, rsqrt( dot ) ), res );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[i+j] = real(res[j]);
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Laplace DLP kernel
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_ansatzsp,
           typename T_testsp,
           typename T_packed >
void
laplace_dlp_simd ( const TGrid::triangle_t &  tri0,
                   const TGrid::triangle_t &  tri1,
                   const tripair_quad_rule_t * rule,
                   vector< real > &           values,
                   const T_ansatzsp *         ansatz_sp,
                   const T_testsp *           test_sp,
                   const T3Point &            n_vec )
{
    using vpacked = T_packed;
    using value_t = typename vpacked::value_t;

    const size_t   VECTOR_SIZE = vpacked::vector_size;
    
    const vpacked  vZERO( 0.0 );
    const vpacked  vONE(  1.0 );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );

    // set normal direction
    const vpacked  n[3]     = { n_vec[0], n_vec[1], n_vec[2] };

    // load triangle coordinates
    const vpacked  t0[3][3] = { { ansatz_sp->grid()->vertex( tri0.vtx[0] )[0],
                                  ansatz_sp->grid()->vertex( tri0.vtx[0] )[1],
                                  ansatz_sp->grid()->vertex( tri0.vtx[0] )[2] },
                                { ansatz_sp->grid()->vertex( tri0.vtx[1] )[0],
                                  ansatz_sp->grid()->vertex( tri0.vtx[1] )[1],
                                  ansatz_sp->grid()->vertex( tri0.vtx[1] )[2] },
                                { ansatz_sp->grid()->vertex( tri0.vtx[2] )[0],
                                  ansatz_sp->grid()->vertex( tri0.vtx[2] )[1],
                                  ansatz_sp->grid()->vertex( tri0.vtx[2] )[2] } };
    
    const vpacked  t1[3][3] = { { test_sp->grid()->vertex( tri1.vtx[0] )[0],
                                  test_sp->grid()->vertex( tri1.vtx[0] )[1],
                                  test_sp->grid()->vertex( tri1.vtx[0] )[2] },
                                { test_sp->grid()->vertex( tri1.vtx[1] )[0],
                                  test_sp->grid()->vertex( tri1.vtx[1] )[1],
                                  test_sp->grid()->vertex( tri1.vtx[1] )[2] },
                                { test_sp->grid()->vertex( tri1.vtx[2] )[0],
                                  test_sp->grid()->vertex( tri1.vtx[2] )[1],
                                  test_sp->grid()->vertex( tri1.vtx[2] )[2] } };
    
    //
    // eval VECTOR_SIZE add the same time
    //
    
    const size_t  n_pts = rule->npts;
    value_t       res[ VECTOR_SIZE ];
    
    #pragma ivdep
    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE )
    {
        // load quadrature points for first triangle
        const vpacked  a1 = load< vpacked >( & rule->x1[i] );  // ( x1[0], x1[1], ... x1[N] )
        const vpacked  a2 = load< vpacked >( & rule->y1[i] );  // ( y1[0], y1[1], ... y1[N] )
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );

        // load quadrature points for second triangle
        const vpacked  b1 = load< vpacked >( & rule->x2[i] );  // ( x2[0], x2[1], ... x2[N] )
        const vpacked  b2 = load< vpacked >( & rule->y2[i] );  // ( y2[0], y2[1], ... y2[N] )
        const vpacked  b0 = sub( sub( vONE, b1 ), b2 );

        //
        // compute || y-x ||^2 and <n,y-x>
        //

        vpacked  dot( 0.0 );
        vpacked  ndu( 0.0 );
        vpacked  tmp1;
        
        // u0-v0
        tmp1 =      muladd( a0, t0[0][0], muladd( a1, t0[1][0], mul( a2, t0[2][0] ) ) );
        tmp1 = sub( muladd( b0, t1[0][0], muladd( b1, t1[1][0], mul( b2, t1[2][0] ) ) ), tmp1 );
        
        // n_x * u_x
        ndu  = mul( n[0], tmp1 );
        dot  = mul( tmp1, tmp1 );

        // u1-v1
        tmp1 =      muladd( a0, t0[0][1], muladd( a1, t0[1][1], mul( a2, t0[2][1] ) ) );
        tmp1 = sub( muladd( b0, t1[0][1], muladd( b1, t1[1][1], mul( b2, t1[2][1] ) ) ), tmp1 );

        // n_y * u_y
        ndu  = muladd( n[1], tmp1, ndu );
        dot  = muladd( tmp1, tmp1, dot );

        // u2-v2
        tmp1 =      muladd( a0, t0[0][2], muladd( a1, t0[1][2], mul( a2, t0[2][2] ) ) );
        tmp1 = sub( muladd( b0, t1[0][2], muladd( b1, t1[1][2], mul( b2, t1[2][2] ) ) ), tmp1 );

        // n_z * u_z
        ndu  = muladd( n[2], tmp1, ndu );
        dot  = muladd( tmp1, tmp1, dot );

        //
        // compute (1 / sqrt(dot))³
        //

        vpacked  tmp0 = rsqrt( dot );
        
        tmp0 = mul( tmp0, mul( tmp0, tmp0 ) );

        // tmp0 contains (1 / sqrt(dot))³; multiply <n,u> and store result
        store( mul( ndu, mul( tmp0, vONE_OVER_4PI ) ), res );
        
        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[i+j] = real(res[j]);
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// Laplace HCA functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

template < typename T_packed >
void
laplace_slp_eval_dx_simd ( const tri_quad_rule_t &   quad_rule,
                           const T3Point             vx[3],
                           const T3Point &           vy,
                           vector< real > &          values )
{
    using vpacked = T_packed;
    using value_t = typename vpacked::value_t;

    const size_t   VECTOR_SIZE = vpacked::vector_size;
    
    const vpacked  vZERO( 0.0 );
    const vpacked  vONE( 1.0 );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    // store given points in AVX variables
    const vpacked  y[3]  = { vy.x(),    vy.y(),    vy.z()    };
    const vpacked  x0[3] = { vx[0].x(), vx[0].y(), vx[0].z() };
    const vpacked  x1[3] = { vx[1].x(), vx[1].y(), vx[1].z() };
    const vpacked  x2[3] = { vx[2].x(), vx[2].y(), vx[2].z() };

    const size_t   npts  = quad_rule.npts;
    value_t        res[ VECTOR_SIZE ];
                
    #pragma ivdep
    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        vpacked  tmp1;
        vpacked  d_dot_d = vZERO;

        // load quadrature points for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // evaluate SLP kernel
        //

        // d0 = y0 - x0
        tmp1    = sub( y[0], muladd( a0, x0[0], muladd( a1, x1[0], mul( a2, x2[0] ) ) ) );
        d_dot_d = mul( tmp1, tmp1 );

        // d1 = y1 - x1
        tmp1    = sub( y[1], muladd( a0, x0[1], muladd( a1, x1[1], mul( a2, x2[1] ) ) ) );
        d_dot_d = muladd( tmp1, tmp1, d_dot_d );

        // d2 = y2 - x2
        tmp1    = sub( y[2], muladd( a0, x0[2], muladd( a1, x1[2], mul( a2, x2[2] ) ) ) );
        d_dot_d = muladd( tmp1, tmp1, d_dot_d );

        // store (1 / sqrt(dot)) in output array
        store< vpacked >( mul( vONE_OVER_4PI, rsqrt( d_dot_d ) ), res );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[k+j] = real(res[j]);
    }// for
}

template < typename T_packed >
void
laplace_dlp_eval_dy_simd ( const tri_quad_rule_t &   quad_rule,
                           const T3Point &           vx,
                           const T3Point             vy[3],
                           const T3Point &           normal,
                           vector< real > &          values )
{
    using vpacked = T_packed;
    using value_t = typename vpacked::value_t;

    const size_t   VECTOR_SIZE = vpacked::vector_size;
    
    const vpacked  vZERO( 0.0 );
    const vpacked  vONE( 1.0 );
    const vpacked  vONE_OVER_4PI( ONE_OVER_4PI );
    
    // store given points in AVX variables
    const vpacked  x[3]  = { vx.x()   , vx.y()   , vx.z()    };
    const vpacked  y0[3] = { vy[0].x(), vy[0].y(), vy[0].z() };
    const vpacked  y1[3] = { vy[1].x(), vy[1].y(), vy[1].z() };
    const vpacked  y2[3] = { vy[2].x(), vy[2].y(), vy[2].z() };
    const vpacked  n[3]  = { normal.x(),
                             normal.y(),
                             normal.z() };

    const size_t    npts  = quad_rule.npts;
    value_t         res[ VECTOR_SIZE ];
                
    #pragma ivdep
    for ( size_t  k = 0; k < npts; k += VECTOR_SIZE )
    {
        vpacked  tmp1;
        vpacked  d_dot_d = vZERO;
        vpacked  n_dot_d = vZERO;
                    
        // load quadrature points for triangle
        const vpacked  a1 = load< vpacked >( & quad_rule.x[k] );
        const vpacked  a2 = load< vpacked >( & quad_rule.y[k] );
        const vpacked  a0 = sub( sub( vONE, a1 ),  a2 );
        
        //
        // evaluate DLP kernel
        //

        // d0 = y0 - x0
        tmp1    = sub( muladd( a0, y0[0], muladd( a1, y1[0], mul( a2, y2[0] ) ) ), x[0] );
        d_dot_d = mul( tmp1, tmp1 );
        n_dot_d = mul( n[0], tmp1 );

        // d1 = y1 - x1
        tmp1    = sub( muladd( a0, y0[1], muladd( a1, y1[1], mul( a2, y2[1] ) ) ), x[1] );
        d_dot_d = muladd( tmp1, tmp1, d_dot_d );
        n_dot_d = muladd( n[1], tmp1, n_dot_d );

        // d2 = y2 - x2
        tmp1    = sub( add( add( mul( a0, y0[2] ), mul( a1, y1[2] ) ), mul( a2, y2[2] ) ), x[2] );
        d_dot_d = muladd( tmp1, tmp1, d_dot_d );
        n_dot_d = muladd( n[2], tmp1, n_dot_d );

        // tmp1 = 1 / |d|³ = 1 / sqrt( d_dot_d )³
        tmp1 = rsqrt( d_dot_d );
        tmp1 = mul( mul( tmp1, tmp1 ), tmp1 );
        
        // multiply (1 / sqrt(dot)³) with <n,d> and store values
        store( mul( vONE_OVER_4PI, mul( n_dot_d, tmp1 ) ), res );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[k+j] = real(res[j]);
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// explicit template instantiations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#define  INST_LAPLACE_SLP_SIMD( T_value, T_ansatzsp, T_testsp )         \
    template void                                                       \
    laplace_slp_simd< T_ansatzsp, T_testsp, packed< T_value, SIMD_ISA > > (     \
        const TGrid::triangle_t &   tri0,                               \
        const TGrid::triangle_t &   tri1,                               \
        const tripair_quad_rule_t * rule,                               \
        vector< real > &            values,                             \
        const T_ansatzsp *          ansatz_sp,                          \
        const T_testsp *            test_sp );

INST_LAPLACE_SLP_SIMD( real, TConstFnSpace,  TConstFnSpace )
INST_LAPLACE_SLP_SIMD( real, TLinearFnSpace, TConstFnSpace )
INST_LAPLACE_SLP_SIMD( real, TConstFnSpace,  TLinearFnSpace )
INST_LAPLACE_SLP_SIMD( real, TLinearFnSpace, TLinearFnSpace )

#define  INST_LAPLACE_DLP_SIMD( T_value, T_ansatzsp, T_testsp )         \
    template void                                                       \
    laplace_dlp_simd< T_ansatzsp, T_testsp, packed< T_value, SIMD_ISA > > ( \
        const TGrid::triangle_t &   tri0,                               \
        const TGrid::triangle_t &   tri1,                               \
        const tripair_quad_rule_t * rule,                               \
        vector< real > &            values,                             \
        const T_ansatzsp *          ansatz_sp,                          \
        const T_testsp *            test_sp,                            \
        const T3Point &             n_vec );

INST_LAPLACE_DLP_SIMD( real, TConstFnSpace,  TConstFnSpace )
INST_LAPLACE_DLP_SIMD( real, TLinearFnSpace, TConstFnSpace )
INST_LAPLACE_DLP_SIMD( real, TConstFnSpace,  TLinearFnSpace )
INST_LAPLACE_DLP_SIMD( real, TLinearFnSpace, TLinearFnSpace )

template
void
laplace_slp_eval_dx_simd< packed< real, SIMD_ISA > >  ( const tri_quad_rule_t &   quad_rule,
                                                        const T3Point             vx[3],
                                                        const T3Point &           vy,
                                                        vector< real > &          values );
template
void
laplace_dlp_eval_dy_simd< packed< real, SIMD_ISA > >  ( const tri_quad_rule_t &   quad_rule,
                                                        const T3Point &           vx,
                                                        const T3Point             vy[3],
                                                        const T3Point &           normal,
                                                        vector< real > &          values );

}// namespace HLIB

// Local Variables:
// mode: c++
// End:
