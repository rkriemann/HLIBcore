#ifndef __HLIB_INIT_HH
#define __HLIB_INIT_HH
//
// Project     : HLib
// File        : init.hh
// Description : initialisation and finalisation
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

namespace HLIB
{

///////////////////////////////////////////////////////////////////
//
// initialisation and finalisation of HLIBpro
//
///////////////////////////////////////////////////////////////////

//!
//! global initialisation routine for HLib
//!
void INIT ();

//!
//! finalisation routing for HLib
//!
void DONE ();

//!
//! return true, if HLIBpro is initialised
//!
bool is_init ();

}// namespace

#endif // __INIT_HH
