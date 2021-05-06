//
// Project     : HLib
// File        : TExpBF.cc
// Description : bilinearforms for exponential kernel
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/base/packed.hh"

#include "hpro/bem/TExpBF.hh"

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
exp_simd  ( const TGrid::triangle_t &   tri0,
            const TGrid::triangle_t &   tri1,
            const tripair_quad_rule_t * rule,
            std::vector< real > &       values,
            const T_ansatzsp *          ansatz_sp,
            const T_testsp *            test_sp );

template < typename T_packed >
void
exp_eval_dx_simd  ( const tri_quad_rule_t &   quad_rule,
                    const T3Point             x[3],
                    const T3Point &           y,
                    vector< real > &          values );

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
exp_flt ( const TGrid::triangle_t &   tri0,
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
        const double  a0 = 1.0 - a1 - a2;
        
        const double  b1 = rule->x2[i];
        const double  b2 = rule->y2[i];
        const double  b0 = 1.0 - b1 - b2;

        // u = x - y
        x = a0 * v00 + a1 * v01 + a2 * v02;
        y = b0 * v10 + b1 * v11 + b2 * v12;
        u = x - y;

        // exp( - |x-y| )
        values[ i ] = Math::exp( - norm2( u ) );

        // x_0 * exp( - |x-y| )
        // values[ i ] = x[0] * Math::exp( - norm2( u ) );
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// TExpBF
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template < typename  T_ansatzsp,
           typename  T_testsp >
TExpBF< T_ansatzsp, T_testsp >::TExpBF ( const ansatzsp_t *  aansatzsp,
                                         const testsp_t *    atestsp,
                                         const uint          quad_order )
        : TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, real >( aansatzsp, atestsp, quad_order )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f()  && CFG::BEM::use_simd_avx512f  )
        {
            HINFO( "(TExpBF) using AVX512F kernel" );
            
            _kernel_fn = exp_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic()  && CFG::BEM::use_simd_mic  )
        {
            HINFO( "(TExpBF) using MIC kernel" );
            
            _kernel_fn = exp_simd< T_ansatzsp, T_testsp, packed< real, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2()  && CFG::BEM::use_simd_avx2  )
        {
            HINFO( "(TExpBF) using AVX2 kernel" );
            
            _kernel_fn = exp_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx()  && CFG::BEM::use_simd_avx  )
        {
            HINFO( "(TExpBF) using AVX kernel" );
            
            _kernel_fn = exp_simd< T_ansatzsp, T_testsp, packed< real, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TExpBF) using SSE3 kernel" );
            
            _kernel_fn = exp_simd< T_ansatzsp, T_testsp, packed< real, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TExpBF) using VSX kernel" );
            
            _kernel_fn = exp_simd< T_ansatzsp, T_testsp, packed< real, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TExpBF) using NEON kernel" );
            
            _kernel_fn = exp_simd< T_ansatzsp, T_testsp, packed< real, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TExpBF) using FPU kernel" );
            
            _kernel_fn = exp_flt;
        }// else
    }// if
    else
    {
        HINFO( "(TExpBF) using FPU kernel" );
            
        _kernel_fn = exp_flt;
    }// else
}

//
// eval kernel function at quadrature points
//
template < typename  T_ansatzsp,
           typename  T_testsp >
void
TExpBF< T_ansatzsp, T_testsp >::eval_kernel ( const                       idx_t,
                                              const                       idx_t,
                                              const TGrid::triangle_t &   tri0,
                                              const TGrid::triangle_t &   tri1,
                                              const tripair_quad_rule_t * rule,
                                              std::vector< real > &       values ) const
{
    _kernel_fn( tri0, tri1, rule, values, this->ansatz_space(), this->test_space() );
}

//
// return format of bilinear form, e.g. symmetric
//
template < typename  T_ansatzsp,
           typename  T_testsp >
matform_t
TExpBF< T_ansatzsp, T_testsp >::format () const
{
    if ( ( cptrcast( this->ansatz_space(), TFnSpace ) == cptrcast( this->test_space(), TFnSpace ) ) ||
         (( this->ansatz_space()->type() == this->test_space()->type() ) &&
          ( this->ansatz_space()->grid() == this->test_space()->grid() )) )
        return symmetric;
    else
        return unsymmetric;
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// TExpGenFn
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

//
// computes x-integral of SLP
//
void
exp_eval_dx_flt ( const tri_quad_rule_t &  quad_rule,
                  const T3Point            vx[3],
                  const T3Point &          vy,
                  vector< real > &         values )
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
        const double   ddotd = double( Math::sqrt( dot( xmy, xmy ) ) );

        // e^( - | x - y | )
        values[k] = Math::exp( - ddotd );
    }// for
}

//
// constructor
//
template < typename  T_ansatzsp,
           typename  T_testsp >
TExpGenFn< T_ansatzsp, T_testsp >::TExpGenFn ( const ansatzsp_t *    ansatzsp,
                                               const testsp_t *      testsp,
                                               const TPermutation *  row_perm_i2e,
                                               const TPermutation *  col_perm_i2e,
                                               const uint            quad_order )
        : TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, real >( ansatzsp, testsp, quad_order,
                                                                 row_perm_i2e, col_perm_i2e )
{
    if ( CFG::BEM::use_simd )
    {
        if ( CFG::Mach::has_avx512f() && CFG::BEM::use_simd_avx512f )
        {
            HINFO( "(TExpGenFn) using AVX512F kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< packed< real, ISA_AVX512F > >;
        }// if
        else if ( CFG::Mach::has_mic() && CFG::BEM::use_simd_mic )
        {
            HINFO( "(TExpGenFn) using MIC kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< packed< real, ISA_MIC > >;
        }// if
        else if ( CFG::Mach::has_avx2() && CFG::BEM::use_simd_avx2 )
        {
            HINFO( "(TExpGenFn) using AVX2 kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< packed< real, ISA_AVX2 > >;
        }// if
        else if ( CFG::Mach::has_avx() && CFG::BEM::use_simd_avx )
        {
            HINFO( "(TExpGenFn) using AVX kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< packed< real, ISA_AVX > >;
        }// if
        else if ( CFG::Mach::has_sse3() && CFG::BEM::use_simd_sse3 )
        {
            HINFO( "(TExpGenFn) using SSE3 kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< packed< real, ISA_SSE3 > >;
        }// if
        else if ( CFG::Mach::has_vsx() && CFG::BEM::use_simd_vsx )
        {
            HINFO( "(TExpGenFn) using VSX kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< packed< real, ISA_VSX > >;
        }// if
        else if ( CFG::Mach::has_neon() && CFG::BEM::use_simd_neon )
        {
            HINFO( "(TExpGenFn) using NEON kernel" );
            
            _eval_dxy_impl = exp_eval_dx_simd< packed< real, ISA_NEON > >;
        }// if
        else
        {
            HINFO( "(TExpGenFn) using standard kernel" );
            
            _eval_dxy_impl = exp_eval_dx_flt;
        }// else
    }// if
    else
    {
        HINFO( "(TExpGenFn) using standard kernel" );
            
        _eval_dxy_impl = exp_eval_dx_flt;
    }// else
}
        
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template class TExpBF< TConstFnSpace,  TConstFnSpace >;
template class TExpBF< TConstFnSpace,  TLinearFnSpace >;
template class TExpBF< TLinearFnSpace, TConstFnSpace >;
template class TExpBF< TLinearFnSpace, TLinearFnSpace >;

template class TExpGenFn< TConstFnSpace,  TConstFnSpace >;
template class TExpGenFn< TConstFnSpace,  TLinearFnSpace >;
template class TExpGenFn< TLinearFnSpace, TConstFnSpace >;
template class TExpGenFn< TLinearFnSpace, TLinearFnSpace >;

}// namespace HLIB
