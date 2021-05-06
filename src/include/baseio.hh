#ifndef __HLIB_BASEIO_HH
#define __HLIB_BASEIO_HH
//
// Project     : HLib
// File        : baseio.hh
// Description : basic IO related functions and classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <istream>
#include <ostream>
#include <string>

#include <boost/cstdint.hpp>

#include "hpro/base/types.hh"

namespace HLIB
{

//////////////////////////////////////////////////////////////
//
// system independent types
//
//////////////////////////////////////////////////////////////

using  int8_t   = boost::int8_t;
using  uint8_t  = boost::uint8_t;
                
using  int16_t  = boost::int16_t;
using  uint16_t = boost::uint16_t;
                
using  int32_t  = boost::int32_t;
using  uint32_t = boost::uint32_t;
                
using  int64_t  = boost::int64_t;
using  uint64_t = boost::uint64_t;

///////////////////////////////////////////////////
// 
// byte swapping (little <-> big endian)
//

template <class T>
void
swap_bytes  ( T &  data );

///////////////////////////////////////////////////
// 
// open file with suitable stream, e.g. for
// compression
//

std::istream *
open_read  ( const std::string &  filename );

std::ostream *
open_write ( const std::string &  filename );

///////////////////////////////////////////////////
// 
// return extension without .gz etc.
//

std::string
extension ( const std::string &  filename );

///////////////////////////////////////////////////
// 
// file format detection
//

enum fmt_t {
        FMT_UNKNOWN,
        FMT_HLIB,
        FMT_HLIB_GRID,
        FMT_OCTAVE,
        FMT_SAMG,
        FMT_MATLAB,
        FMT_HDF5,
        FMT_HB,
        FMT_MTX,
        FMT_PLTMG,
        FMT_PLY,
        FMT_SURFMESH,
        FMT_GMSH
};

//
// guess file format of file \a filename
// - if \a filename exists, content of it
//   may be read
//
fmt_t
guess_format ( const std::string &  filename );


///////////////////////////////////////////////////
// 
// lexical cast with default value
//

long
str_to_int ( const std::string &  str,
             const long           def_val = 0 );

long
str_to_int ( const char *         begin,
             const char *         end,
             const long           def_val = 0 );

double
str_to_dbl ( const std::string &  str,
             const double         def_val = 0 );

double
str_to_dbl ( const char *         begin,
             const char *         end,
             const double         def_val = 0 );

///////////////////////////////////////////////////
//
// append file suffix if not already present
//

std::string
add_extension ( const std::string &  filename,
                const std::string &  suffix );

}// namespace HLIB

#endif  // __HLIB_BASEIO_HH
