#ifndef __HLIB_TMAXWELLBF_HH
#define __HLIB_TMAXWELLBF_HH
//
// Project     : HLib
// File        : TMaxwellBF.hh
// Description : bilinear forms for Maxwell operator
// Author      : Ronald Kriemann, Jonas Ballani
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/bem/TQuadBEMBF.hh"

namespace HLIB
{

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// standard bilinear forms for Maxwell
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
//
// Electric Field Integral Equation for Maxwell
//
////////////////////////////////////////////////////

template < typename  T_ansatzsp,
           typename  T_testsp >
class TMaxwellEFIEBF : public TQuadBEMBF< T_ansatzsp, T_testsp, complex >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t     = T_ansatzsp;
    using  testsp_t       = T_testsp;
    using  value_t        = complex;

    using  ansatz_value_t = typename ansatzsp_t::value_t;
    using  test_value_t   = typename testsp_t::value_t;

private:

    //! @cond
    
    // i · wave number
    const complex  _ikappa;

    // √( mu / epsilon )
    const real     _eta;

    // constants used in bilinear form
    const complex  _ikappa_times_eta;
    const complex  _eta_over_ikappa;

    // kernel function implementation
    void ( * _kernel_fn ) ( const TGrid::triangle_t &    tri0,
                            const TGrid::triangle_t &    tri1,
                            const tripair_quad_rule_t *  rule,
                            const complex                ikappa,
                            const ansatzsp_t *           ansatz_sp,
                            const testsp_t *             test_sp,
                            std::vector< complex > &     values );
    
    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TMaxwellEFIEBF ( const complex       kappa,
                     const real          eta,
                     const ansatzsp_t *  aansatzsp,
                     const testsp_t *    atestsp,
                     const uint          quad_order = CFG::BEM::quad_order );

    virtual ~TMaxwellEFIEBF () {}

    // return format of bilinear form, e.g. symmetric
    matform_t  format () const;

    //! evaluate subblock defined by \a row_ind × \a col_ind; the indices
    //! in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
    //! contiguous
    virtual void  eval  ( const std::vector< idx_t > &  row_ind,
                          const std::vector< idx_t > &  col_ind,
                          BLAS::Matrix< complex > &     values ) const;



protected:
    // eval kernel function at quadrature points
    virtual void  eval_kernel  ( const idx_t                  tri0idx,
                                 const idx_t                  tri1idx,
                                 const TGrid::triangle_t &    tri0,
                                 const TGrid::triangle_t &    tri1,
                                 const tripair_quad_rule_t *  quad_rule,
                                 std::vector< complex > &     values ) const;
};

//!
//! \class  TMaxwellEFIEMassBF
//! \brief  bilinear form for Maxwell EFIE mass matrix
//!
template < class T_ansatzsp,
           class T_testsp >
class TMaxwellEFIEMassBF : public TBEMBF< T_ansatzsp,
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

    TMaxwellEFIEMassBF ( const ansatzsp_t *  ansatzsp,
                         const testsp_t *    testsp,
                         const uint          quadorder = CFG::BEM::quad_order );

    virtual ~TMaxwellEFIEMassBF () {}

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


////////////////////////////////////////////////////
//
// Magnetic Field Integral Equation for Maxwell
//
////////////////////////////////////////////////////

template < typename  T_ansatzsp,
           typename  T_testsp >
class TMaxwellMFIEBF : public TQuadBEMBF< T_ansatzsp, T_testsp, complex >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t     = T_ansatzsp;
    using  testsp_t       = T_testsp;
    using  value_t        = complex;

    using  ansatz_value_t = typename ansatzsp_t::value_t;
    using  test_value_t   = typename testsp_t::value_t;

private:

    //! @cond
    
    // i · wave number
    const complex  _ikappa;

    // kernel function implementation
    void ( * _kernel_fn ) ( const TGrid::triangle_t &    tri0,
                            const TGrid::triangle_t &    tri1,
                            const tripair_quad_rule_t *  rule,
                            const complex                ikappa,
                            const ansatzsp_t *           ansatz_sp,
                            const testsp_t *             test_sp,
                            std::vector< complex > &     values );
    
    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TMaxwellMFIEBF ( const complex       kappa,
                     const ansatzsp_t *  aansatzsp,
                     const testsp_t *    atestsp,
                     const uint          quad_order = CFG::BEM::quad_order );

    virtual ~TMaxwellMFIEBF () {}

    // return format of bilinear form, e.g. symmetric
    matform_t  format () const;

    //! evaluate subblock defined by \a row_ind × \a col_ind; the indices
    //! in \a row_ind and \a col_ind can be arbitrary, e.g. must not be
    //! contiguous
    virtual void  eval  ( const std::vector< idx_t > &  row_ind,
                          const std::vector< idx_t > &  col_ind,
                          BLAS::Matrix< value_t > &     values ) const;

protected:
    // eval kernel function at quadrature points
    virtual void  eval_kernel  ( const idx_t                  tri0idx,
                                 const idx_t                  tri1idx,
                                 const TGrid::triangle_t &    tri0,
                                 const TGrid::triangle_t &    tri1,
                                 const tripair_quad_rule_t *  quad_rule,
                                 std::vector< complex > &     values ) const;
};

//!
//! \class  TMaxwellMFIEMassBF
//! \brief  bilinear form for Maxwell MFIE mass matrix
//!
template < class T_ansatzsp,
           class T_testsp >
class TMaxwellMFIEMassBF : public TBEMBF< T_ansatzsp,
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

    TMaxwellMFIEMassBF ( const ansatzsp_t *  ansatzsp,
                         const testsp_t *    testsp,
                         const uint          order = CFG::BEM::quad_order );

    virtual ~TMaxwellMFIEMassBF () {}

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


}// namespace

#endif  // __TMAXWELLBF_HH
