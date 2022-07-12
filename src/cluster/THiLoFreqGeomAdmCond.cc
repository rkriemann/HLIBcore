//
// Project     : HLIBpro
// File        : THiLoFreqGeomAdmCond.cc
// Description : 
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/error.hh"
#include "hpro/cluster/TGeomCluster.hh"

#include "hpro/cluster/TGeomAdmCond.hh"

namespace Hpro
{

THiLoFreqGeomAdmCond::THiLoFreqGeomAdmCond ( const std::complex< double >  kappa,
                                             const double                  max_waves,
                                             const double                  eta )
        : _kappa( Math::abs(kappa) )
        , _eta( eta )
        , _max_waves( max_waves )
{}

//
// check admissibility
//
bool
THiLoFreqGeomAdmCond::is_adm ( const TBlockCluster * c ) const
{
    if ( c == nullptr )
        HERROR( ERR_ARG, "(THiLoFreqGeomAdmCond) is_adm", "cluster is NULL" );
    
    const TGeomCluster  * rowcl, * colcl;

    if ( ! IS_TYPE( c->rowcl(), TGeomCluster ) ||
         ! IS_TYPE( c->colcl(), TGeomCluster ) )
        HERROR( ERR_CT_TYPE, "(THiLoFreqGeomAdmCond) is_adm", 
               c->rowcl()->typestr() + " x " + c->colcl()->typestr() );
    
    rowcl = cptrcast( c->rowcl(), TGeomCluster );
    colcl = cptrcast( c->colcl(), TGeomCluster );

    if ( rowcl == colcl )
        return false;

    //
    // build balls around center with bb-radius as radius of ball
    //

    const double  diam_rowcl = rowcl->diameter();
    const double  diam_colcl = colcl->diameter();

    #if 1

    //
    // determine minimal axis in bboxes of cluster
    //
    
    const TBBox   row_bbox   = rowcl->bbox();
    const TBBox   col_bbox   = colcl->bbox();
    const uint    dim        = row_bbox.min().dim();
    double        row_min    = row_bbox.max()[0] - row_bbox.min()[0];
    double        col_min    = col_bbox.max()[0] - col_bbox.min()[0];
    
    for ( uint  i = 1; i < dim; ++i )
    {
        row_min = std::min( row_min, row_bbox.max()[i] - row_bbox.min()[i] );
        col_min = std::min( col_min, col_bbox.max()[i] - col_bbox.min()[i] );
    }// for
    
    const double  r_kappa    = std::min( row_min, col_min ) * _kappa;

    #else

    //
    // use minimal diameter
    //
    
    const double  r_kappa    = std::min( diam_rowcl, diam_colcl ) * _kappa;

    #endif

    if ( r_kappa <= _max_waves )
    {
        //
        // low frequency regime: use standard admissibility
        //

        return Math::max( diam_rowcl, diam_colcl ) <= ( _eta * rowcl->distance( colcl ) );
    }// if
    else
    {
        return false;
        
        // //
        // // high frequency regime: use directional parabolic separation condition
        // //

        // const double  alpha     = 1.0 / r_kappa;
        // const double  tan_alpha = std::tan( alpha );
        // const double  dt        = r_kappa / tan_alpha;
        // const double  ds        = dt + distance;
        
        // // return _kappa * Math::square( Math::max( diam_rowcl, diam_colcl ) ) <= ( _eta_high * distance );
        // return (distance >= Math::square( diam_rowcl )) && ( _kappa * diam_colcl <= ds * tan_alpha);







        
        // //
        // // high frequency regime: use directional parabolic separation condition
        // //

        // const double  alpha        = 1.0 / r_kappa;
        // const double  tan_alpha    = std::tan( alpha );
        // const double  center_dist  = distance + radius_rowcl;  // distance from center of rowcl to colcl
        // const double  dt           = r_kappa / tan_alpha;
        // const double  radius_colcl = diam_colcl / 2.0;
        // const double  ds           = dt + center_dist + radius_colcl;

        // if ( ( _kappa * center_dist >= Math::square( r_kappa )) &&
        //      ( _kappa * radius_colcl <= ds * tan_alpha) )
        // {
        //     DBG::noop();
        // }
        
        // // return _kappa * Math::square( Math::max( diam_rowcl, diam_colcl ) ) <= ( _eta_high * distance );
        // return ( _kappa * center_dist >= Math::square( r_kappa )) &&
        //     ( _kappa * radius_colcl <= ds * tan_alpha);
        
    }// switch

    return false; 
}

}// namespace Hpro
