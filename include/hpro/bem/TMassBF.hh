#ifndef __HLIB_TMASSBF_HH
#define __HLIB_TMASSBF_HH
//
// Project     : HLib
// File        : TMassBF.hh
// Description : bilinear form for mass matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/bem/TBEMBF.hh"

namespace HLIB
{
    
//!
//! \class  TMassBF
//! \brief  bilinear form for mass matrix
//!
template < class T_ansatzsp,
           class T_testsp >
class TMassBF : public TBEMBF< T_ansatzsp,
                               T_testsp,
                               real >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = real;

protected:
    
    // quadrature rule
    std::vector< T2Point >  _quad_pts;
    std::vector< double >   _quad_wghts;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TMassBF ( const ansatzsp_t *  aansatzsp,
              const testsp_t *    atestsp,
              const uint          order = CFG::BEM::quad_order );

    virtual ~TMassBF () {}

    // return format of bilinear form, e.g. symmetric
    matform_t  format () const;

    //////////////////////////////////////
    //
    // evaluate bilinearform
    //

    //! evaluate subblock defined by \a row_ind × \a col_ind; the indices
    //! in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
    //! contiguous
    virtual void  eval  ( const std::vector< idx_t > &  row_ind,
                          const std::vector< idx_t > &  col_ind,
                          BLAS::Matrix< value_t > &     values ) const;
};

}// namespace HLIB

#endif  // __HLIB_TMASSBF_HH
