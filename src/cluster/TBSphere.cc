//
// Project     : HLIBpro
// File        : TBSphere.cc
// Description : class for a bounding sphere
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2024. All Rights Reserved.
//

#include "hpro/cluster/TBSphere.hh"

namespace Hpro
{

///////////////////////////////////////////////
//
// bounding bsphere properties
//

//
// return true if given point is inside bsphere
//
bool
TBSphere::is_inside ( const TPoint & x ) const
{
    if ( x.dim() != center().dim() )
        HERROR( ERR_DIM, "(TBSphere) is_inside", "given point has wrong dimension" );

    return ( _center - x ).norm2() <= _radius;
}
    
//
// return volume to given bsphere
//
double
TBSphere::volume () const
{
    switch ( _center.dim() )
    {
        case  1 : return 2.0                              * _radius;
        case  2 : return             Math::pi< double >() * _radius * _radius;
        case  3 : return 4.0 / 3.0 * Math::pi< double >() * _radius * _radius * _radius;
        default :
            HERROR( ERR_DIM, "(TBSphere) volume", "unsupported dimension" );
    }// switch
}

//
// return distance to given bsphere
//
double
TBSphere::distance ( const TBSphere &  bsphere ) const
{
    const auto  dist = ( _center - bsphere._center ).norm2() - _radius - bsphere._radius;
    
    return  std::max( dist, 0.0 );
}
    
///////////////////////////////////////////////
//
// misc.
//

//
// copy operator
//
TBSphere &
TBSphere::operator = ( const TBSphere & bsphere )
{
    _center = bsphere._center;
    _radius = bsphere._radius;

    return *this;
}

//
// return string representation
//
std::string
TBSphere::to_string () const
{
    return "( " + _center.to_string() + ", " + Hpro::to_string( "%f", _radius ) + " )";
}
    
}// namespace Hpro
