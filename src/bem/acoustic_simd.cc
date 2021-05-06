//
// Project     : HLib
// File        : acoustic_simd.inc
// Description : Acoustic Scattering kernel using SIMD instructions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#if !defined(SIMD_ISA)
#  error "no SIMD ISA defined"
#endif

#include <vector>

#include "hpro/base/packed.hh"

#include "hpro/bem/TAcousticBF.hh"

namespace HLIB
{

using std::vector;

//
// local functions and constants
//

namespace
{

const real  ONE_OVER_2PI = real(1) / (real(2) * Math::pi< real >());

//
// get normal direction for all quadrature points in vector data
//
template < typename  T_ansatzsp,
           typename  T_packed >
inline
void
get_normal ( const T_packed             va1,
             const T_packed             va2,
             const idx_t                tri0idx,
             const TGrid::triangle_t &  tri0,
             const T_ansatzsp *         ansatz_sp,
             T_packed                   vn[3] )
{
    using vpacked = T_packed;
    using value_t = typename vpacked::value_t;

    const size_t   VECTOR_SIZE = vpacked::vector_size;
    
    //
    // convert vector data to array
    //
    
    value_t  a1[ VECTOR_SIZE ], a2[ VECTOR_SIZE ];
    value_t  x[ VECTOR_SIZE ],  y[ VECTOR_SIZE ],  z[ VECTOR_SIZE ];

    store( va1, a1 );
    store( va2, a2 );
    
    for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
    {
        const T3Point  n( ansatz_sp->grid()->tri_normal( tri0idx, tri0, a1[j], a2[j] ) );

        x[j] = n.x();
        y[j] = n.y();
        z[j] = n.z();
    }// for

    vn[0] = load< vpacked >( x );
    vn[1] = load< vpacked >( y );
    vn[2] = load< vpacked >( z );
}

}// namespace anonymous

//
// kernel function: mic implementation
//
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
                std::vector< complex > &     values )
{
    using vpacked = T_packed;
    using value_t = typename vpacked::value_t;
    
    const size_t   VECTOR_SIZE = vpacked::vector_size;
        
    const vpacked  vONE_OVER_2PI( ONE_OVER_2PI );
    const vpacked  vONE( 1.0 );
    const vpacked  vZERO( 0.0 );

    //
    // load triangle coordinates
    //

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
    // eval 4 add the same time
    //

    const vpacked     vikappa_re( ikappa.real() );
    const vpacked     vikappa_im( ikappa.imag() );
    const size_t      n_pts = rule->npts;
    value_t           res_re[ VECTOR_SIZE ];
    value_t           res_im[ VECTOR_SIZE ];
    
    #pragma ivdep
    for ( size_t  i = 0; i < n_pts; i += VECTOR_SIZE  )
    {
        //
        // load quadrature points
        //
        
        // for first triangle
        const vpacked  a1  = load< vpacked >( & rule->x1[i] );
        const vpacked  a2  = load< vpacked >( & rule->y1[i] );
        const vpacked  a0  = sub( sub( vONE, a1 ),  a2 );

        // for second triangle
        const vpacked  b1  = load< vpacked >( & rule->x2[i] );
        const vpacked  b2  = load< vpacked >( & rule->y2[i] );
        const vpacked  b0  = sub( sub( vONE, b1 ), b2 );

        //
        // compute x-y
        //

        vpacked  diff[3];
        
        diff[0] = sub( muladd( a0, t0[0][0], muladd( a1, t0[1][0], mul( a2, t0[2][0] ) ) ),
                       muladd( b0, t1[0][0], muladd( b1, t1[1][0], mul( b2, t1[2][0] ) ) ) );

        diff[1] = sub( muladd( a0, t0[0][1], muladd( a1, t0[1][1], mul( a2, t0[2][1] ) ) ),
                       muladd( b0, t1[0][1], muladd( b1, t1[1][1], mul( b2, t1[2][1] ) ) ) );
        
        diff[2] = sub( muladd( a0, t0[0][2], muladd( a1, t0[1][2], mul( a2, t0[2][2] ) ) ),
                       muladd( b0, t1[0][2], muladd( b1, t1[1][2], mul( b2, t1[2][2] ) ) ) );

        // |x-y|²
        const vpacked  vdist2 = muladd( diff[0], diff[0], muladd( diff[1], diff[1], mul( diff[2], diff[2] ) ) );
        
        // |x-y|
        const vpacked  vdist  = sqrt( vdist2 );

        // |x-y|³
        const vpacked  vdist3 = mul( vdist, vdist2 );
        
        //
        // get normal directions
        //

        vpacked  n[3];

        get_normal( a1, a2, tri0idx, tri0, ansatz_sp, n );
        
        // < n, y-x >
        const vpacked  vndu = muladd( diff[0], n[0], muladd( diff[1], n[1], mul( diff[2], n[2] ) ) );
        
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
        const vpacked  vikd_re = mul( vikappa_re, vdist );
        
        // im( ( i·κ·d[0], i·κ·d[1], i·κ·d[2], i·κ·d[3] ) )
        const vpacked  vikd_im = mul( vikappa_im, vdist );
        
        // exp( re(·) )
        // (also: multiply with 1/(4·π) <n,y-x> / |x-y|³ as this would be done later anyway)
        const vpacked  vexp_ikd_re = mul( div( mul( vONE_OVER_2PI, vndu ), vdist3 ), exp( vikd_re ) ); 

        // cos( i·κ·d[0..3] ), sin( i·κ·d[0..3] )
        // exp(re)·cos(), exp(re)·sin()   (here we would multiply with 1/(2·π) <n,y-x> / |x-y|³)
        vpacked  vc, vs;

        sincos( vikd_im, vs, vc );

        vc = mul( vexp_ikd_re, vc );
        vs = mul( vexp_ikd_re, vs );

        // 1 - i·κ·|x-y| 
        const vpacked  vone_m_ikd_re = sub( vONE,  vikd_re );
        const vpacked  vone_m_ikd_im = sub( vZERO, vikd_im );
        
        // compute exp() · (1 - i·κ·|x-y|) and store final values
        store( sub( mul( vc, vone_m_ikd_re ), mul( vs, vone_m_ikd_im ) ), res_re );
        store( add( mul( vc, vone_m_ikd_im ), mul( vs, vone_m_ikd_re ) ), res_im );

        for ( size_t  j = 0; j < VECTOR_SIZE; ++j )
            values[ i+j ] = complex( res_re[j], res_im[j] );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// explicit template instantiations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#define  INST_ACOUSTIC( T_ansatzsp, T_testsp )                   \
    template void                                                \
    acoustic_simd< T_ansatzsp, T_testsp, packed< real, SIMD_ISA > > (   \
        const complex               ikappa,                      \
        const idx_t                 tri0idx,                     \
        const TGrid::triangle_t &   tri0,                        \
        const TGrid::triangle_t &   tri1,                        \
        const tripair_quad_rule_t * rule,                        \
        const T_ansatzsp *          ansatz_sp,                   \
        const T_testsp *            test_sp,                     \
        std::vector< complex > &    values )

INST_ACOUSTIC( TConstFnSpace,  TConstFnSpace );
INST_ACOUSTIC( TLinearFnSpace, TConstFnSpace );
INST_ACOUSTIC( TConstFnSpace,  TLinearFnSpace );
INST_ACOUSTIC( TLinearFnSpace, TLinearFnSpace );

}// namespace HLIB

// Local Variables:
// mode: c++
// End:
