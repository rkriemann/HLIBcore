#ifndef __HLIB_TACOUSTICBF_HH
#define __HLIB_TACOUSTICBF_HH
//
// Project     : HLib
// File        : TAcousticBF.hh
// Description : bilinear form for acoustic scattering
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/bem/TQuadBEMBF.hh"

namespace HLIB
{
    
////////////////////////////////////////////////////
//
// Accoustic Scattering operator
//
////////////////////////////////////////////////////

//!
//! \ingroup BEM_Module
//! \class   TAcousticScatterBF
//! \brief   Bilinear form for acoustic scattering
//!
template < typename  T_ansatzsp,
           typename  T_testsp >
class TAcousticScatterBF : public TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, complex >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = complex;

private:
    //! @cond

    // i · wave number
    const complex  _ikappa;
    
    // kernel function implementation
    void ( * _kernel_fn ) ( const complex                ikappa,
                            const idx_t                  tri0idx,
                            const TGrid::triangle_t &    tri0,
                            const TGrid::triangle_t &    tri1,
                            const tripair_quad_rule_t *  rule,
                            const ansatzsp_t *           ansatz_sp,
                            const testsp_t *             test_sp,
                            std::vector< complex > &     values );
    
    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //
    
    TAcousticScatterBF ( const complex       kappa,
                         const ansatzsp_t *  aansatzsp,
                         const testsp_t *    atestsp,
                         const uint          quad_order = CFG::BEM::quad_order );

    virtual ~TAcousticScatterBF () {}

protected:
    // eval kernel function at quadrature points
    virtual void  eval_kernel  ( const idx_t                  tri0idx,
                                 const idx_t                  tri1idx,
                                 const TGrid::triangle_t &    tri0,
                                 const TGrid::triangle_t &    tri1,
                                 const tripair_quad_rule_t *  quad_rule,
                                 std::vector< complex > &     values ) const;
};

}// namespace

#endif  // __HLIB_TACOUSTICBF_HH
