//
// Project     : HLib
// File        : TMaternCovCoeffFn.cc
// Description : matrix coefficients for matern covariance kernel
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hlib-config.h"

#if USE_GSL == 1
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#else
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>
#endif

#include "hpro/misc/TMaternCovCoeffFn.hh"

namespace HLIB
{

namespace
{

//
// some constants
//
const double  SQRT_3 = std::sqrt( 3 );
const double  SQRT_5 = std::sqrt( 5 );

//
// compute modified bessel function
//
real
compute_Bessel ( const real  d,
                 const real  rho,
                 const real  nu,
                 const real  sigmasq,
                 const real  scale_fac ) // constant scaling factor for optimization
{
    if ( d < Limits::epsilon< real >() )
        return sigmasq;
    else
    {
        const double temp =  d / rho;

        #if USE_GSL == 1
        return scale_fac * gsl_sf_bessel_Knu( nu, temp ) * std::pow( temp, nu );
        #else
        return scale_fac * boost::math::cyl_bessel_k( nu, temp ) * std::pow( temp, nu );
        #endif
    }// else
}

}// namespace anonymous

//
// constructor
//
template < typename T_point >
TMaternCovCoeffFn< T_point >::TMaternCovCoeffFn ( const real                      sigma,
                                                  const real                      length,
                                                  const real                      nu,
                                                  const std::vector< T_point > &  vertices )
    : TCoeffFn< real >()
    , _length( length )
    , _nu( nu )
    , _sigmasq( sigma*sigma )
    #if USE_GSL == 1
    , _scale_fac( _sigmasq / ( std::pow(2.0, nu-1) * gsl_sf_gamma(nu) ) )
    #else
    , _scale_fac( _sigmasq / ( std::pow(2.0, nu-1) * boost::math::tgamma( nu ) ) )
    #endif
    , _vertices( vertices )
    , _matern_type( general )
{}

//
// coefficient evaluation
//
template < typename T_point >
void
TMaternCovCoeffFn< T_point >::eval  ( const std::vector< idx_t > &  rowidxs,
                                      const std::vector< idx_t > &  colidxs,
                                      real *                        matrix ) const
{
    const size_t  n = rowidxs.size();
    const size_t  m = colidxs.size();

    for ( size_t  j = 0; j < m; ++j )
    {
        const idx_t      idx1 = colidxs[ j ];
        const T_point &  y    = _vertices[ idx1 ];
            
        for ( size_t  i = 0; i < n; ++i )
        {
            const idx_t      idx0 = rowidxs[ i ];
            const T_point &  x    = _vertices[ idx0 ];
            const real       dist = norm2( x - y );

            switch ( _matern_type )
            {
                case one_half :
                {
                    matrix[ j*n + i ] = _sigmasq * std::exp( - dist / _length );
                }
                break;

                case three_half :
                {
                    const auto  sqrt3drho = SQRT_3 * dist / _length;
                    
                    matrix[ j*n + i ] = _sigmasq * ( 1.0 + sqrt3drho ) * std::exp( - sqrt3drho );
                }
                break;

                case five_half :
                {
                    const auto  sqrt5drho = SQRT_5 * dist / _length;
                    
                    matrix[ j*n + i ] = _sigmasq * ( 1.0 + sqrt5drho + ( 5 * dist * dist ) / ( 3 * _length * _length ) ) * std::exp( - sqrt5drho );
                }
                break;

                case general :
                default :
                {
                    matrix[ j*n + i ] = compute_Bessel( dist, _length, _nu, _sigmasq, _scale_fac );
                }
            }// switch
        }// for
    }// for
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
//
// explicit template instatiation
//
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

template class TMaternCovCoeffFn< TPoint >;
template class TMaternCovCoeffFn< T2Point >;
template class TMaternCovCoeffFn< T3Point >;

}// namespace HLIB
