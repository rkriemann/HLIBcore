//
// Project     : HLib
// File        : TStdGeomAdmCond.cc
// Description : standard admissibility classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/base/error.hh"
#include "hpro/cluster/TGeomCluster.hh"

#include "hpro/cluster/TGeomAdmCond.hh"

namespace HLIB
{

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TStdGeomAdmCond (implementation)
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
TStdGeomAdmCond::is_adm ( const TBlockCluster * c ) const
{
    const TGeomCluster  * rowcl, * colcl;

    if ( ! IS_TYPE( c->rowcl(), TGeomCluster ) ||
         ! IS_TYPE( c->colcl(), TGeomCluster ) )
        HERROR( ERR_CT_TYPE, "(TStdGeomAdmCond) is_adm", 
               c->rowcl()->typestr() + " x " + c->colcl()->typestr() );
    
    rowcl = cptrcast( c->rowcl(), TGeomCluster );
    colcl = cptrcast( c->colcl(), TGeomCluster );

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
// TWeakStdGeomAdmCond (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//
// apply weak admissibility for small clusters and
// standard admissibility for large clusters
//
bool
TWeakStdGeomAdmCond::is_adm ( const TBlockCluster * c ) const
{
    const TGeomCluster  * rowcl, * colcl;

    if ( ! IS_TYPE( c->rowcl(), TGeomCluster ) ||
         ! IS_TYPE( c->colcl(), TGeomCluster ) )
        HERROR( ERR_CT_TYPE, "(TWeakStdGeomAdmCond) is_adm", 
               c->rowcl()->typestr() + " x " + c->colcl()->typestr() );
    
    rowcl = cptrcast( c->rowcl(), TGeomCluster );
    colcl = cptrcast( c->colcl(), TGeomCluster );

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

    const uint  dim     = rowcl->bbox().min().dim();
    uint        n_empty = 0;
    
    // in 1D, different clusters share at most one vertex => admissible
    if ( dim == 1 )
        return true;
    
    if (( rowcl->bbox().max().dim() < dim ) ||
        ( colcl->bbox().min().dim() < dim ) ||
        ( colcl->bbox().max().dim() < dim ))
        HERROR( ERR_ARG, "(TWeakStdGeomAdmCond) is_adm",
               "dimension of vectors in rowcl and colcl differs" );

    for ( uint i = 0; i < dim; i++ )
    {
        if (( rowcl->bbox().max()[i] <= colcl->bbox().min()[i] ) ||
            ( colcl->bbox().max()[i] <= rowcl->bbox().min()[i] ))
            n_empty++;
    }// for

    // test real weak admissibility, i.e., no overlap
    if ( n_empty == dim )
        return true;
    
    // if "small" use weak admissibility for lower dimension
    if (( std::max( rowcl->size(), colcl->size() ) < 1000 ) && ( n_empty == dim-1 ))
        return true;
    
    //
    // test standard adm.
    //
    
    return TStdGeomAdmCond::is_adm( c );
}

}// namespace HLIB
