#ifndef __HPRO_TGEOMCLUSTER_HH
#define __HPRO_TGEOMCLUSTER_HH
//
// Project     : HLIBpro
// File        : TGeomCluster.hh
// Description : class for a cluster with geometrical information
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/types.hh"
#include "hpro/cluster/TBBox.hh"
#include "hpro/cluster/TCluster.hh"

namespace Hpro
{

//
// local type
//
DECLARE_TYPE( TGeomCluster );

//!
//! \ingroup  Cluster_Module
//! \class    TGeomCluster
//! \brief    Extend standard cluster by bounding box.
//!
class TGeomCluster : public TCluster
{
private:
    //! the bounding box of the cluster
    TBBox  _bbox;

public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct empty cluster
    TGeomCluster () {}

    //! construct cluster with index set [\a first_idx, \a last_idx]
    TGeomCluster ( const idx_t  first_idx,
                   const idx_t  last_idx )
            : TCluster( first_idx, last_idx )
    {}
    
    //! construct cluster with index set [\a first_idx, \a last_idx]
    TGeomCluster ( const idx_t    first_idx,
                   const idx_t    last_idx,
                   const TBBox &  abbox )
            : TCluster( first_idx, last_idx ), _bbox( abbox )
    {}
    
    ///////////////////////////////////////////////
    //
    // access local variables
    //

    //! return bounding box
    TBBox &       bbox ()       { return _bbox; }
    const TBBox & bbox () const { return _bbox; }

    //! set bounding box
    void  set_bbox ( const TBBox &  abbox )
    {
        _bbox = abbox;
    }
    
    ///////////////////////////////////////////////
    //
    // geometrical cluster properties
    //

    //! return diameter of cluster
    double diameter () const
    {
        return bbox().diameter();
    }
    
    //! return distance to cluster \a cl
    double distance ( const TGeomCluster * cl ) const
    {
        return bbox().distance( cl->bbox() );
    }
    
    //! return distance to cluster \a cl where \a period defines periodicity of coordinates
    double distance ( const TGeomCluster * cl,
                      const TPoint &       period ) const
    {
        return bbox().distance( cl->bbox(), period );
    }
    
    ////////////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! return object of same type
    virtual TGeomCluster * create () const { return new TGeomCluster; }

    //! return size in bytes used by this object
    virtual size_t byte_size () const
    {
        return TCluster::byte_size() + _bbox.byte_size();
    }
    
    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //

    HPRO_RTTI_DERIVED( TGeomCluster, TCluster );
};

//
// return true of given cluster is of type TGeomCluster
//
inline bool is_geom_cluster  ( const TCluster *  cl ) noexcept { return IS_TYPE( cl, TGeomCluster ); }
inline bool is_geom_cluster  ( const TCluster &  cl ) noexcept { return is_geom_cluster( & cl ); }

}// namespace Hpro

#endif  // __HPRO_TGEOMCLUSTER_HH
