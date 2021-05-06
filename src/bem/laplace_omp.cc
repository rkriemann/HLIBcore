//
// Project     : HLib
// File        : laplace_omp.cc
// Description : Laplace kernels using OpenMP SIMD functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>

#include "hpro/bem/TLaplaceBF.hh"

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

const double  ONE_OVER_4PI = 1.0 / (4.0 * Math::pi< double >());

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// auxiliary functions
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#pragma omp declare simd
double
rsqrt ( const double x )
{
    return 1.0 / std::sqrt( x );
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
// Laplace SLP kernel
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#pragma omp declare simd uniform(v00,v01,v02,v10,v11,v12) linear(i:1)
template < typename value_t >
void
laplace_slp_kernel ( size_t          i,
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
                     value_t *       values )
{
    const auto  a1 = x1[i];
    const auto  a2 = y1[i];
    const auto  a0 = 1.0 - a1 - a2;
        
    const auto  b1 = x2[i];
    const auto  b2 = y2[i];
    const auto  b0 = 1.0 - b1 - b2;

    // u = x - y
    const double  X0 = a0 * v00[0] + a1 * v01[0] + a2 * v02[0];
    const double  X1 = a0 * v00[1] + a1 * v01[1] + a2 * v02[1];
    const double  X2 = a0 * v00[2] + a1 * v01[2] + a2 * v02[2];
    
    const double  Y0 = b0 * v10[0] + b1 * v11[0] + b2 * v12[0];
    const double  Y1 = b0 * v10[1] + b1 * v11[1] + b2 * v12[1];
    const double  Y2 = b0 * v10[2] + b1 * v11[2] + b2 * v12[2];

    const double  u0 = X0 - Y0;
    const double  u1 = X1 - Y1;
    const double  u2 = X2 - Y2;

    //
    //   1
    // ─────
    // |x-y|
    //
        
    values[i] = ( ONE_OVER_4PI * rsqrt( normsqr( u0, u1, u2 ) ) );
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
laplace_slp_omp ( const TGrid::triangle_t &   tri0,
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

    const auto *   x1 = rule->x1.data();
    const auto *   y1 = rule->y1.data();
    const auto *   x2 = rule->x2.data();
    const auto *   y2 = rule->y2.data();

    const double *  cv00 = v00.vector();
    const double *  cv01 = v01.vector();
    const double *  cv02 = v02.vector();
    const double *  cv10 = v10.vector();
    const double *  cv11 = v11.vector();
    const double *  cv12 = v12.vector();

    real *          cvalues = values.data();
        
    #pragma omp simd linear(x1,y1,x2,y2:1) aligned(x1,y1,x2,y2:64)
    for ( size_t  i = 0; i < n_pts; ++i )
    {
        laplace_slp_kernel( i, x1, y1, x2, y2,
                            cv00, cv01, cv02, cv10, cv11, cv12,
                            cvalues );
    }// for
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
//
// explicit template instantiations
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

#define  INST_LAPLACE_SLP_OMP( T_ansatzsp, T_testsp )                   \
    template void                                                       \
    laplace_slp_omp< T_ansatzsp, T_testsp > (                           \
        const TGrid::triangle_t &   tri0,                               \
        const TGrid::triangle_t &   tri1,                               \
        const tripair_quad_rule_t * rule,                               \
        std::vector< real > &       values,                             \
        const T_ansatzsp *          ansatz_sp,                          \
        const T_testsp *            test_sp );

INST_LAPLACE_SLP_OMP( TConstFnSpace,  TConstFnSpace )
INST_LAPLACE_SLP_OMP( TLinearFnSpace, TConstFnSpace )
INST_LAPLACE_SLP_OMP( TConstFnSpace,  TLinearFnSpace )
INST_LAPLACE_SLP_OMP( TLinearFnSpace, TLinearFnSpace )

}// namespace HLIB
