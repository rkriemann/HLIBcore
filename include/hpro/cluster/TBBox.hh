#ifndef __HPRO_TBBOX_HH
#define __HPRO_TBBOX_HH
//
// Project     : HLIBpro
// File        : TBBox.hh
// Description : class for a axis aligned bounding box 
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/types.hh"
#include "hpro/base/TPoint.hh"

namespace Hpro
{

/////////////////////////////////////////////////////////////////
//
// type for a bounding box
//
class TBBox
{
private:
    // minimal and maximal point of box
    TPoint  _bb_min, _bb_max;

public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor for empty bbox
    TBBox ()
    {}
    
    //! ctor for bbox [\a bbmin, \a bbmax ]
    TBBox ( const TPoint &  bbmin,
            const TPoint &  bbmax )
            : _bb_min( bbmin )
            , _bb_max( bbmax )
    {
        if ( bbmin.dim() != bbmax.dim() )
            HERROR( ERR_ARG, "(TBBox) ctor", "different spatial dimension in bbox coordinates" );
    }

    //! ctor for bbox [\a bbmin, \a bbmax ] (special version for T2Point)
    TBBox ( const T2Point &  bbmin,
            const T2Point &  bbmax )
            : _bb_min( 2, bbmin.vector() )
            , _bb_max( 2, bbmax.vector() )
    {}

    //! ctor for bbox [\a bbmin, \a bbmax ] (special version for T3Point)
    TBBox ( const T3Point &  bbmin,
            const T3Point &  bbmax )
            : _bb_min( 3, bbmin.vector() )
            , _bb_max( 3, bbmax.vector() )
    {}

    //! copy ctor
    TBBox ( const TBBox & box )
    {
        *this = box;
    }
    
    ///////////////////////////////////////////////
    //
    // access local variables
    //

    //! return minimal coordinate of bbox
    TPoint &       min ()       { return _bb_min; }
    const TPoint & min () const { return _bb_min; }

    //! return maximal coordinate of bbox
    TPoint &       max ()       { return _bb_max; }
    const TPoint & max () const { return _bb_max; }

    //! return spatial dimension of bbox
    uint           dim () const { return _bb_min.dim(); }
    
    ///////////////////////////////////////////////
    //
    // bounding box properties
    //

    //! return true if point \a x is inside box
    bool is_inside ( const TPoint & x ) const;
    
    //! return diameter of box
    double diameter () const;

    //! return distance to \a box
    double distance ( const TBBox & box ) const;

    //! return distance to \a box but coordinates have
    //! periodicity defined by \a period
    double distance ( const TBBox &   box,
                      const TPoint &  period ) const;

    //! join local bbox with \a box
    void join ( const TBBox & box );
    
    ///////////////////////////////////////////////
    //
    // misc.
    //

    //! return size in bytes used by this object
    size_t byte_size () const
    { return _bb_min.byte_size() + _bb_max.byte_size(); }

    //! copy operator
    TBBox & operator = ( const TBBox & box );

    //! return string representation
    std::string  to_string () const;
};

}// namespace

#endif  // __HPRO_TBBOX_HH
