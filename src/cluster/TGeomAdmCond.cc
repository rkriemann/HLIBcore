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
// TWeakGeomAdmCond (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//
// test positive distance for admissibility
//
bool
TWeakGeomAdmCond::is_adm ( const TBlockCluster *  bcl ) const
{
    if ( ! ( is_geom_cluster( bcl->rowcl() ) && is_geom_cluster( bcl->colcl() ) ) )
        HERROR( ERR_CT_TYPE, "(TWeakGeomAdmCond) is_adm", bcl->rowcl()->typestr() + " x " + bcl->colcl()->typestr() );
    
    auto  rowcl = cptrcast( bcl->rowcl(), TGeomCluster );
    auto  colcl = cptrcast( bcl->colcl(), TGeomCluster );

    // identical clusters are always inadmissibile
    if ( rowcl == colcl )
        return false;

    //
    // nested dissection case: domain-domain clusters are admissible
    //

    if ( rowcl->is_domain() && colcl->is_domain() )
        return true;

    //
    // compute number of axes without overlap
    //

    const uint  dim        = rowcl->bbox().min().dim();
    uint        nnooverlap = 0;
    
    const auto  rbox = rowcl->bbox();
    const auto  cbox = colcl->bbox();

    if ( ! (( rbox.max().dim() == dim ) && 
            ( cbox.min().dim() == dim ) &&
            ( cbox.max().dim() == dim ) ) )
        HERROR( ERR_DIM, "(TWeakGeomAdmCond) is_adm", "cluster dimensions differ" );

    for ( uint  i = 0; i < dim; ++i )
    {
        const auto  rmin   = rbox.min()[i];
        const auto  rmax   = rbox.max()[i];
            
        const auto  cmin   = cbox.min()[i];
        const auto  cmax   = cbox.max()[i];
            
        const auto  rlen   = rmax - rmin;
        const auto  clen   = cmax - cmin;
        const auto  minlen = std::min( rlen, clen );
            
        if (( rmax <= cmin ) ||   // ├── τ ──┤├── σ ──┤
            ( cmax <= rmin ))     // ├── σ ──┤├── τ ──┤
        {
            nnooverlap++;
        }// if
        else if ( _hoverlap > 0 )
        {
            //
            // test relative overlap size
            //

            // test ├── τ ──┼h┼── σ ──┤
            if (( rmax >= cmin ) && ( rmax <= cmax ) && (( rmax - cmin ) < _hoverlap ))
                nnooverlap++;
                    
            // test ├── σ ──┼h┼── τ ──┤
            else if (( cmax >= rmin ) && ( cmax <= rmax ) && (( cmax - rmin ) < _hoverlap ))
                nnooverlap++;
        }// if
    }// for

    // auto  inter = intersection( rbox, cbox );
    // auto  l0    = inter.max()[0] - inter.min()[0];
    // auto  l1    = inter.max()[1] - inter.min()[1];
    // auto  l2    = 0; // inter.max()[2] - inter.min()[2];
    // double  ol  = 1e10;

    // for ( uint  i = 0; i < 2; ++i )
    //     if ( inter.max()[i] - inter.min()[i] > 0 )
    //         ol = std::min( ol, inter.max()[i] - inter.min()[i] );
        
    // if ( nnooverlap <= 1 )
    //     std::cout << rbox.to_string() << " x " << cbox.to_string() << " : "
    //               << term::bold() << nnooverlap << ", " << ol << term::reset()
    //               << " / " << l0 << " / " << l1 << " / " << l2 << std::endl;
    // else if ( nnooverlap == 2 )
    //     std::cout << term::red() << rbox.to_string() << " x " << cbox.to_string() << " : "
    //               << term::bold() << nnooverlap << ", " << ol << term::reset()
    //               << " / " << l0 << " / " << l1 << " / " << l2 << std::endl;
    // else if ( nnooverlap == 3 )
    //     std::cout << term::green() << rbox.to_string() << " x " << cbox.to_string() << " : "
    //               << term::bold() << nnooverlap << ", " << ol << term::reset()
    //               << " / " << l0 << " / " << l1 << " / " << l2 << std::endl;
        
    if ( nnooverlap >= _noverlap )
        return true;
    
    //
    // test standard admissibility
    //

    return std::min( rowcl->diameter(), colcl->diameter() ) <= ( _eta * rowcl->distance( colcl ) );
}

}// namespace Hpro
