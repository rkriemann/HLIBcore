//
// Project     : HLIBpro
// File        : TGeomAdmCond.cc
// Description : classes for geometric admissibility conditions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2024. All Rights Reserved.
//

#include <hpro/base/error.hh>
#include <hpro/cluster/TGeomCluster.hh>

#include <hpro/cluster/TGeomAdmCond.hh>

namespace Hpro
{

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TStrongGeomAdmCond (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

///////////////////////////////////////////
//
// check block-cluster if admissible
//

//
// implement basic rule for admissible clusters :
// min{ diam(t_1) , diam(t_2) } <= 2 eta dist(t_1, t_2)
//
bool
TStrongGeomAdmCond::is_adm ( const TBlockCluster * c ) const
{
    if ( ! ( is_geom_cluster( c->rowcl() ) && is_geom_cluster( c->colcl() ) ) )
        HERROR( ERR_CT_TYPE, "(TStrongGeomAdmCond) is_adm", c->rowcl()->typestr() + " x " + c->colcl()->typestr() );
    
    auto  rowcl = cptrcast( c->rowcl(), TGeomCluster );
    auto  colcl = cptrcast( c->colcl(), TGeomCluster );

    if ( rowcl == colcl )
        return false;

    //
    // if both clusters are domain clusters it is admissible (since rowcl != colcl)
    //

    if ( rowcl->is_domain() && colcl->is_domain() )
        return true;
    
    //
    // build balls around center with bb-radius as radius of ball
    //

    double  distance;

    if ( _period.dim() != 0 )
        distance = rowcl->distance( colcl, _period );
    else
        distance = rowcl->distance( colcl );

    switch ( _diam_mode )
    {
        case  use_max_diam :
            return ( std::max( rowcl->diameter(), colcl->diameter() ) <= ( _eta * distance ) );

        default :
        case  use_min_diam :
            return ( std::min( rowcl->diameter(), colcl->diameter() ) <= ( _eta * distance ) );
    }// switch

    return false; 
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TVertexGeomAdmCond (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//
// apply weak admissibility for small clusters and
// standard admissibility for large clusters
//
bool
TVertexGeomAdmCond::is_adm ( const TBlockCluster * c ) const
{
    if ( ! ( is_geom_cluster( c->rowcl() ) && is_geom_cluster( c->colcl() ) ) )
        HERROR( ERR_CT_TYPE, "(TVertexGeomAdmCond) is_adm", c->rowcl()->typestr() + " x " + c->colcl()->typestr() );
    
    auto  rowcl = cptrcast( c->rowcl(), TGeomCluster );
    auto  colcl = cptrcast( c->colcl(), TGeomCluster );

    if ( rowcl == colcl )
        return false;

    //
    // if both clusters are domain clusters it is admissible (since rowcl != colcl)
    //

    if ( rowcl->is_domain() && colcl->is_domain() )
        return true;

    //
    // count number of empty intersections per axis (up to single point)
    //

    // test real weak admissibility, i.e., no overlap
    if ( rowcl->bvol().overlap_dim( colcl->bvol() ) == 0 )
        return true;
    
    //
    // test standard adm.
    //
    
    return TStrongGeomAdmCond::is_adm( c );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TWeakGeomAdmCond (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//
// test positive distance for admissibility
//
bool
TWeakGeomAdmCond::is_adm ( const TBlockCluster * c ) const
{
    if ( ! ( is_geom_cluster( c->rowcl() ) && is_geom_cluster( c->colcl() ) ) )
        HERROR( ERR_CT_TYPE, "(TWeakGeomAdmCond) is_adm", c->rowcl()->typestr() + " x " + c->colcl()->typestr() );
    
    auto  rowcl = cptrcast( c->rowcl(), TGeomCluster );
    auto  colcl = cptrcast( c->colcl(), TGeomCluster );

    if ( rowcl == colcl )
        return false;

    //
    // if both clusters are domain clusters it is admissible (since rowcl != colcl)
    //

    if ( rowcl->is_domain() && colcl->is_domain() )
        return true;

    //
    // test distance between clusters
    //
    
    return rowcl->distance( colcl ) > 0.0;
}

}// namespace Hpro
