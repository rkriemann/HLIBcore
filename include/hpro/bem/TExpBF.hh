#ifndef __HLIB_TEXPBF_HH
#define __HLIB_TEXPBF_HH
//
// Project     : HLib
// File        : TEXPBF.hh
// Description : bilinear forms for exponential kernel
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/bem/TQuadBEMBF.hh"
#include "hpro/bem/TQuadHCAGenFn.hh"

namespace HLIB
{
    
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// standard bilinear forms for exp kernel
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
//
// single layer potential of Exp operator
//
////////////////////////////////////////////////////

//!
//! \ingroup BEM_Module
//! \class   TExpBF
//! \brief   Bilinear form for expontential kernel
//!
//!          TExpBF implements the bilinear form for the expontential
//!          kernel \f[ e^{ - \|x-y\|_2} \f].
//!
template < typename  T_ansatzsp,
           typename  T_testsp >
class TExpBF : public TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, real >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = real;

public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TExpBF ( const ansatzsp_t *  aansatzsp,
             const testsp_t *    atestsp,
             const uint          quad_order = CFG::BEM::quad_order );

    virtual ~TExpBF () {}

    // return format of bilinear form, e.g. symmetric
    matform_t  format () const;
            
protected:
    //
    // eval kernel function at quadrature points
    //
    virtual void  eval_kernel  ( const idx_t                  tri0idx,
                                 const idx_t                  tri1idx,
                                 const TGrid::triangle_t &    tri0,
                                 const TGrid::triangle_t &    tri1,
                                 const tripair_quad_rule_t *  quad_rule,
                                 std::vector< real > &        values ) const;

private:

    //
    // kernel function impl.
    //

    void  ( * _kernel_fn ) ( const TGrid::triangle_t &    tri0,
                             const TGrid::triangle_t &    tri1,
                             const tripair_quad_rule_t *  rule,
                             std::vector< real > &        values,
                             const ansatzsp_t *           ansatz_sp,
                             const testsp_t *             test_sp );

};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// HCA functions for Exp Kernel
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

//!
//! \ingroup BEM_Module
//! \class   TExpGenFn
//! \brief   kernel generator function for Exp Kernel
//!
template < typename T_ansatzsp,
           typename T_testsp >
class TExpGenFn : public TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, real >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = real;
    
    //
    // inherit from base class
    //
    using  stat_t     = typename TQuadHCAGenFn< ansatzsp_t, testsp_t, value_t >::stat_t;
    
private:

    //
    // HCA function impl.
    //

    void ( * _eval_dxy_impl ) ( const tri_quad_rule_t &   quad_rule,
                                const T3Point             x[3],
                                const T3Point &           y,
                                std::vector< value_t > &  values );
    
public:
    //
    // constructor
    //
    TExpGenFn ( const ansatzsp_t *    ansatzsp,
                const testsp_t *      testsp,
                const TPermutation *  row_perm_i2e,
                const TPermutation *  col_perm_i2e,
                const uint            quad_order = CFG::BEM::quad_order );
    
    //!
    //! evaluate generator function at (\a x, \a y)
    //!
    value_t  eval  ( const T3Point &  x,
                     const T3Point &  y ) const
    {
        return std::exp( - norm2( x - y ) );
    }

protected:
    //! Evaluate \f$ D_x \gamma(x, y) \f$ on with \f$x\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dx  ( const idx_t               tri_idx,
                             const T3Point &           y,
                             const tri_quad_rule_t &   quad_rule,
                             std::vector< value_t > &  values ) const
    {
        const TGrid *            grid = this->ansatz_space()->grid();
        const TGrid::triangle_t  tri  = grid->triangle( tri_idx );
        const T3Point            x[3] = { grid->vertex( tri.vtx[0] ),
                                          grid->vertex( tri.vtx[1] ),
                                          grid->vertex( tri.vtx[2] ) };
    
        _eval_dxy_impl( quad_rule, x, y, values );
    }
    
    //! Evaluate \f$ D_y \gamma(x, y) \f$ on with \f$y\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dy  ( const T3Point &           x,
                             const idx_t               tri_idx,
                             const tri_quad_rule_t &   quad_rule,
                             std::vector< value_t > &  values ) const
    {
        const TGrid *            grid = this->test_space()->grid();
        const TGrid::triangle_t  tri  = grid->triangle( tri_idx );
        const T3Point            y[3] = { grid->vertex( tri.vtx[0] ),
                                          grid->vertex( tri.vtx[1] ),
                                          grid->vertex( tri.vtx[2] ) };
    
        _eval_dxy_impl( quad_rule, y, x, values );
    }

    DISABLE_COPY_OP( TExpGenFn );
};

}// namespace HLIB

#endif  // __HLIB_TEXPBF_HH
