#ifndef __HPRO_TSTDADMCOND_HH
#define __HPRO_TSTDADMCOND_HH
//
// Project     : HLIBpro
// File        : TGeomAdmCond.hh
// Description : baseclass for all admissibility-conditions
//               for block-clusters
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include "hpro/base/TPoint.hh"
#include "hpro/cluster/types.hh"

#include "hpro/cluster/TAdmCondition.hh"

namespace Hpro
{

//!
//! \ingroup  Cluster_Module
//! \class    TStdGeomAdmCond
//! \brief    Standard admissibility for FEM/BEM applications
//!           normal :  adm  iff  min( diam(τ), diam(σ) ) ≤ η·dist(τ,σ)
//!           use_max:  adm  iff  max( diam(τ), diam(σ) ) ≤ η·dist(τ,σ)
//!
class TStdGeomAdmCond : public TAdmCondition
{
protected:
    //! parameter for ratio between diameter and distance
    const double       _eta;

    //! choose between min/max cluster diameters
    const diam_mode_t  _diam_mode;

    //! defines periodicity of coordinates
    const TPoint       _period;

public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct standard admissiblity 
    TStdGeomAdmCond ( const double       eta       = 2.0,
                      const diam_mode_t  diam_mode = use_min_diam )
            : _eta(eta)
            , _diam_mode(diam_mode)
    {}

    //! construct standard admissiblity with periodic geometry
    TStdGeomAdmCond ( const TPoint &     period,
                      const double       eta       = 2.0,
                      const diam_mode_t  diam_mode = use_min_diam )
            : _eta(eta)
            , _diam_mode(diam_mode)
            , _period(period)
    {}

    //! dtor
    virtual ~TStdGeomAdmCond () {}

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if cluster \a cl is admissible
    virtual bool is_adm ( const TBlockCluster * cl ) const;

    DISABLE_COPY_OP( TStdGeomAdmCond );
};

//!
//! \ingroup  Cluster_Module
//! \class    TWeakStdGeomAdmCond
//! \brief    Combination of standard and weak (vertex) admissibility
//!
class TWeakStdGeomAdmCond : public TStdGeomAdmCond
{
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    TWeakStdGeomAdmCond ( const double  eta = 2.0 )
            : TStdGeomAdmCond( eta, use_min_diam )
    {}

    TWeakStdGeomAdmCond ( const TPoint &  period,
                          const double    eta = 2.0 )
            : TStdGeomAdmCond( period, eta, use_min_diam )
    {}

    virtual ~TWeakStdGeomAdmCond () {}

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if \a cl is either weakly or strongly admissible
    virtual bool is_adm ( const TBlockCluster * cl ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TWeakGeomAdmCond
//! \brief    blocks are admissible if clusters have positive distance
//!
class TWeakGeomAdmCond : public TAdmCondition
{
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    TWeakGeomAdmCond ()
            : TAdmCondition()
    {}

    virtual ~TWeakGeomAdmCond () {}

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if \a cl is either weakly or strongly admissible
    virtual bool is_adm ( const TBlockCluster * cl ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    THiLoFreqGeomAdmCond
//! \brief    Admissibility for high and low frequency regimes
//!
class THiLoFreqGeomAdmCond : public TAdmCondition
{
protected:
    // (modulus of) wave number 
    const double  _kappa;
    
    // parameter η of standard admissibility
    const double  _eta;

    // maximal number of wave lengths per cluster
    const double  _max_waves;

public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    THiLoFreqGeomAdmCond ( const std::complex< double >  kappa,
                           const double                  max_waves,
                           const double                  eta = 2.0 );

    virtual ~THiLoFreqGeomAdmCond () {}

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if \a cl is admissible
    virtual bool is_adm ( const TBlockCluster * cl ) const;
};

}// namespace Hpro

#endif  // __HPRO_TGEOMADMCOND_HH
