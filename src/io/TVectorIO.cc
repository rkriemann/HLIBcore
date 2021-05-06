//
// Project     : HLib
// File        : TVectorIO.cc
// Description : classes for matrix input/output
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include "hpro/vector/TScalarVector.hh"

#include "baseio.hh"

#include "hpro/io/TVectorIO.hh"

namespace HLIB
{

namespace fs = boost::filesystem;

using std::string;
using std::vector;
using std::unique_ptr;
using std::make_unique;

using boost::to_lower_copy;
using boost::spirit::qi::long_;
using boost::spirit::qi::double_;
using boost::spirit::qi::_1;
using boost::spirit::qi::phrase_parse;
using boost::spirit::ascii::space;
using boost::phoenix::ref;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// class TVectorIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

TVectorIO::TVectorIO ()
{
}

void
TVectorIO::write ( const TVector *, const string & ) const
{
    HERROR( ERR_NOT_IMPL, "(TVectorIO) write", "" );
}

unique_ptr< TVector >
TVectorIO::read ( const string & ) const
{
    HERROR( ERR_NOT_IMPL, "(TVectorIO) read", "" );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// class TAutoVectorIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void
TAutoVectorIO::write ( const TVector *  A,
                       const string &   filename ) const
{
    if ( filename == "" )
        HERROR( ERR_ARG, "(TAutoVectorIO) write", "empty filename" );
    
    unique_ptr< TVectorIO >  mio;
        
    switch ( guess_format( filename ) )
    {
        case FMT_HLIB    : mio = make_unique< THLibVectorIO >();   break;
        case FMT_MATLAB  : mio = make_unique< TMatlabVectorIO >(); break;
        case FMT_SAMG    : mio = make_unique< TSAMGVectorIO >();   break;
        case FMT_HB      : mio = make_unique< THBVectorIO >();     break;
        case FMT_MTX     : mio = make_unique< TMMVectorIO >();     break;

        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoVectorIO) read", "in file " + filename );
    }// switch

    mio->write( A, filename );
}

void
TAutoVectorIO::write ( const TVector *  A,
                       const string &   filename,
                       const string &   vecname ) const
{
    if ( filename == "" )
        HERROR( ERR_ARG, "(TAutoVectorIO) write", "empty filename" );
    
    unique_ptr< TVectorIO >  mio;
    bool                  support_named = false;
        
    switch ( guess_format( filename ) )
    {
        case FMT_HLIB    : mio = make_unique< THLibVectorIO >();   break;
        case FMT_MATLAB  : mio = make_unique< TMatlabVectorIO >(); support_named = true; break;
        case FMT_SAMG    : mio = make_unique< TSAMGVectorIO >();   break;
        case FMT_HB      : mio = make_unique< THBVectorIO >();     break;
        case FMT_MTX     : mio = make_unique< TMMVectorIO >();     break;
            
        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoVectorIO) read", "in file " + filename );
    }// switch

    if ( support_named )
        ptrcast( mio.get(), TMatlabVectorIO )->write( A, filename, vecname );
    else
        mio->write( A, filename );
}

unique_ptr< TVector >
TAutoVectorIO::read  ( const string &  filename ) const
{
    if ( filename == "" )
        HERROR( ERR_ARG, "(TAutoVectorIO) read", "empty filename" );

    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TAutoVectorIO) read", filename );
    
    unique_ptr< TVectorIO >  mio;
        
    switch ( guess_format( filename ) )
    {
        case FMT_HLIB    : mio = make_unique< THLibVectorIO >();   break;
        case FMT_MATLAB  : mio = make_unique< TMatlabVectorIO >(); break;
        case FMT_SAMG    : mio = make_unique< TSAMGVectorIO >();   break;
        case FMT_HB      : mio = make_unique< THBVectorIO >();     break;
        case FMT_MTX     : mio = make_unique< TMMVectorIO >();     break;

        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoVectorIO) read", "in file " + filename );
    }// switch

    return mio->read( filename );
}

unique_ptr< TVector >
TAutoVectorIO::read  ( const string &  filename,
                       const string &  vecname ) const
{
    if ( filename == "" )
        HERROR( ERR_ARG, "(TAutoVectorIO) read", "empty filename" );

    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TAutoVectorIO) read", filename );
    
    unique_ptr< TVectorIO >  mio;
    bool                  support_named = false;
        
    switch ( guess_format( filename ) )
    {
        case FMT_HLIB    : mio = make_unique< THLibVectorIO >();   break;
        case FMT_MATLAB  : mio = make_unique< TMatlabVectorIO >(); support_named = true; break;
        case FMT_SAMG    : mio = make_unique< TSAMGVectorIO >();   break;
        case FMT_HB      : mio = make_unique< THBVectorIO >();     break;
        case FMT_MTX     : mio = make_unique< TMMVectorIO >();     break;

        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoVectorIO) read", "in file " + filename );
    }// switch

    if ( support_named )
        return ptrcast( mio.get(), TMatlabVectorIO )->read( filename, vecname );
    else
        return mio->read( filename );
}

unique_ptr< TVector >
read_vector  ( const char *  filename )
{
    TAutoVectorIO  mio;

    return mio.read( filename );
}

unique_ptr< TVector >
read_vector  ( const char *  filename,
               const char *  vecname )
{
    TAutoVectorIO  mio;

    return mio.read( filename, vecname );
}

void
write_vector  ( const TVector *  A,
                const char *     filename )
{
    TAutoVectorIO  mio;

    mio.write( A, filename );
}

void
write_vector  ( const TVector *  A,
                const char *     filename,
                const char *     vecname )
{
    TAutoVectorIO  mio;

    mio.write( A, filename, vecname );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class TSAMGVectorIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// write vector in SAMG format
//
void
TSAMGVectorIO::write ( const TVector *  v,
                       const string &   filename ) const
{
    if ( v == nullptr )
        HERROR( ERR_ARG, "(TSAMGVectorIO) write", "vector is nullptr" );

    if ( ! IS_TYPE( v, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TSAMGVectorIO) write", v->typestr() );

    if ( v->is_complex() )
        HERROR( ERR_NREAL, "(TSAMGVectorIO) write", "SAMG only supports real valued data" );
    
    //
    // write vector
    //

    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();
    const TScalarVector *       sv  = cptrcast( v, TScalarVector );
    
    for ( uint i = 0; i < v->size(); i++ )
        out << boost::format( "%.8g" ) % sv->blas_rvec()(i) << std::endl;
}

//
// read vector in SAMG format
//
unique_ptr< TVector >
TSAMGVectorIO::read ( const string & filename ) const
{
    //
    // open format file to read further information
    //

    fs::path  filepath( filename );
    fs::path  ext = filepath.extension();
    string    compress_ext;
    fs::path  frmpath;

    while (( to_lower_copy( ext.string() ) == ".gz"    ) ||
           ( to_lower_copy( ext.string() ) == ".bzip2" ) ||
           ( to_lower_copy( ext.string() ) == ".zlib"  ))
    {
        compress_ext = compress_ext + ext.string();
        filepath.replace_extension( "" );
        ext = filepath.extension();
    }// while
    
    if (( ext == ".rhs" ) || ( ext == ".sol" ))
        frmpath = filepath.replace_extension( ".frm" );
    else
        frmpath = filepath.string() + ".frm";

    unique_ptr< std::istream >  in_ptr;
    
    if ( fs::exists( frmpath ) )
        in_ptr.reset( open_read( frmpath.string() ) );
    else
    {
        if ( fs::exists( frmpath.string() + compress_ext ) )
            in_ptr.reset( open_read( frmpath.string() + compress_ext ) );
        else
            HERROR( ERR_FNEXISTS, "(TSAMGMatrixIO) read", frmpath.string() + "(" + compress_ext + ")" );
    }// else

    std::istream &  frm_in = * in_ptr.get();
    
    ulong             n;
    string            line;
    vector< string >  parts;

    std::getline( frm_in, line );
    std::getline( frm_in, line );
    split( line, " \t\r\n", parts );

    if ( parts.size() < 3 )
        HERROR( ERR_FMT_SAMG, "(TSAMGMatrixIO) read", "expected <nnz> <nrows> <type> in format file" );
    
    n = str_to_int( parts[1] );

    //
    // read vector
    //

    in_ptr.reset( open_read( filename ) );

    std::istream &  amg_in = * in_ptr.get();

    auto  v = make_unique< TScalarVector >();
    
    v->set_size( n );

    for ( uint i = 0; i < n; i++ )
    {
        double  val;
            
        std::getline( amg_in, line );
        val = str_to_dbl( line );

        v->set_entry( i, real( val ) );
    }// for
    
    return unique_ptr< TVector >( v.release() );
}

///////////////////////////////////////////////////
// 
// input and output in Harwell-Boeing and the
// Rutherford-Boeing format
//

///////////////////////////////////////////////////
// 
// input and output in MatrixMarket format
//

enum mtxformat_t { MTX_ARRAY, MTX_COORD };
enum mtxfield_t  { MTX_REAL, MTX_INTEGER, MTX_COMPLEX, MTX_PATTERN };
enum mtxsym_t    { MTX_GENERAL, MTX_SYM, MTX_SKEWSYM, MTX_HERM };

//
// write matrix to <fname>
//
void
TMMVectorIO::write ( const TVector *  x,
                     const string &   filename ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();

    //
    // file header defining data type and size
    //
    
    out << "%%MatrixMarket matrix array";

    if ( x->is_complex() ) out << " complex";
    else                   out << " real";

    out << " general" << std::endl;

    out << "% vector generated by HLIBpro" << std::endl;
    out << x->size() << " 1 " << std::endl;  // n×1 matrix

    //
    // write vector coefficients
    //

    if ( x->is_complex() )
    {
        for ( auto  i : x->is() )
        {
            const auto  f = x->centry( i );
            
            out << boost::format( "%.16e %.16e" ) % std::real( f ) % std::imag( f ) << std::endl;
        }// for
    }// if
    else
    {
        for ( auto  i : x->is() )
        {
            const auto  f = x->entry( i );
            
            out << boost::format( "%.16e" ) % f << std::endl;
        }// for
    }// else
}

//
// read matrix from <fname>
//
unique_ptr< TVector >
TMMVectorIO::read  ( const string & filename ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TMMVectorIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    ///////////////////////////////////////////////////
    //
    // read header and determine format of file
    //

    string            line;
    vector< string >  parts;
    mtxformat_t       mat_format = MTX_ARRAY;
    mtxfield_t        mat_field  = MTX_REAL;
    mtxsym_t          mat_sym    = MTX_GENERAL;

    while ( in.good() )
    {
        std::getline( in, line );
        split( line, " \t\r\n", parts );

        if ( parts[0] == "%%MatrixMarket" )
        {
            if ( parts.size() < 3 )
                HERROR( ERR_FMT_MTX, "(TMMVectorIO) read",
                        "missing header data, expected <object> <format>" );
            
            if ( parts[1] != "matrix" )
                HERROR( ERR_FMT_MTX, "(TMMVectorIO) read", "file does not contain matrix" );
            
            if      ( parts[2] == "coordinate" ) mat_format = MTX_COORD;
            else if ( parts[2] == "array" )      mat_format = MTX_ARRAY;
            else
                HERROR( ERR_FMT_MTX, "(TMMVectorIO) read", "unknown matrix format \"" + parts[2] + "\"" );

            // look for qualifiers
            for ( uint i = 3; i < parts.size(); i++ )
            {
                if      ( parts[i] == "real"           ) mat_field = MTX_REAL;
                else if ( parts[i] == "integer"        ) mat_field = MTX_INTEGER;
                else if ( parts[i] == "complex"        ) mat_field = MTX_COMPLEX;
                else if ( parts[i] == "pattern"        ) mat_field = MTX_PATTERN;
                else if ( parts[i] == "general"        ) mat_sym   = MTX_GENERAL;
                else if ( parts[i] == "symmetric"      ) mat_sym   = MTX_SYM;
                else if ( parts[i] == "skew-symmetric" ) mat_sym   = MTX_SKEWSYM;
                else if ( parts[i] == "hermitian"      ) mat_sym   = MTX_HERM;
                else
                    HERROR( ERR_FMT_MTX, "(TMMVectorIO) read", "unknown qualifier \"" + parts[i] + "\"" );
            }// for
            
            break;
        }// if
    }// while

    if ( mat_sym != MTX_GENERAL )
        HERROR( ERR_NOT_IMPL, "(TMMVectorIO) read", "only unsymmetric matrices supported" );
    
    ///////////////////////////////////////////////////
    //
    // read content depending on format
    //

    if ( mat_format == MTX_ARRAY )
    {
        //
        // read size of vector
        //

        uint  nrows = 0, ncols = 0;

        while ( in.good() )
        {
            std::getline( in, line );

            if ( line[0] == '%' )
                continue;

            const bool  r = phrase_parse( line.begin(), line.end(),
                                          ( long_[ ref(nrows) = _1 ] >>
                                            long_[ ref(ncols) = _1 ] ),
                                          space );

            if ( ! r )
                HERROR( ERR_FMT_MTX, "(TMMVectorIO) read", "missing data, expected <nrows> <ncols>" );

            break;
        }// while

        HINFO( to_string( "(TMMVectorIO) read : reading vector of size %d", nrows ) );

        if ( ncols > 1 )
            HWARNING( "(TMMVectorIO) read : only reading first vector in file" );
        
        //
        // now, first go through file and count entries per row
        //

        const bool  is_complex = ( mat_field == MTX_COMPLEX );
        auto        v          = make_unique< TScalarVector >( nrows, 0, is_complex );

        if ( is_complex )
        {
            for ( uint  i = 0; i < nrows; ++i )
            {
                std::getline( in, line );

                double      re = 0, im = 0;
                const bool  r  = phrase_parse( line.begin(), line.end(),
                                               ( double_[ ref(re) = _1 ] >>
                                                 double_[ ref(im) = _1 ] ),
                                               space );
                
                if ( ! r )
                    HERROR( ERR_FMT_MTX, "(TMMVectorIO) read",
                            "missing data, expected <real> <imag>" );
                
                const complex  val( re, im );

                v->set_centry( i, val );
            }// for
        }// if
        else
        {
            for ( uint  i = 0; i < nrows; ++i )
            {
                std::getline( in, line );

                double      val = 0;
                const bool  r   = phrase_parse( line.begin(), line.end(),
                                                ( double_[ ref(val) = _1 ] ),
                                                space );
                
                if ( ! r )
                    HERROR( ERR_FMT_MTX, "(TMMVectorIO) read",
                            "missing data, expected <float>" );
                
                v->set_entry( i, real(val) );
            }// for
        }// else

        return unique_ptr< TVector >( v.release() );
    }// if
    else
    {
        //
        // read size of vector
        //

        uint  nrows = 0, ncols = 0, nnz = 0;

        while ( in.good() )
        {
            std::getline( in, line );

            if ( line[0] == '%' )
                continue;

            const bool  r = phrase_parse( line.begin(), line.end(),
                                          ( long_[ ref(nrows) = _1 ] >>
                                            long_[ ref(ncols) = _1 ] >>
                                            long_[ ref(nnz)   = _1 ] ),
                                          space );

            if ( ! r )
                HERROR( ERR_FMT_MTX, "(TMMVectorIO) read", "missing data, expected <nrows> <ncols> <nnz>" );

            break;
        }// while

        HINFO( to_string( "(TMMVectorIO) read : reading vector of size %d", nrows ) );

        if ( ncols > 1 )
            HWARNING( "(TMMVectorIO) read : only reading first vector in file" );
        
        //
        // read coefficients
        //

        const bool  is_complex = ( mat_field == MTX_COMPLEX );
        auto        v          = make_unique< TScalarVector >( nrows, 0, is_complex );

        if ( is_complex )
        {
            for ( uint  i = 0; i < nrows; ++i )
            {
                std::getline( in, line );

                long        row = 0, col = 0;
                double      re = 0, im = 0;
                const bool  r  = phrase_parse( line.begin(), line.end(),
                                               ( long_[ ref(row)  = _1 ] >>
                                                 long_[ ref(col)  = _1 ] >>
                                                 double_[ ref(re) = _1 ] >>
                                                 double_[ ref(im) = _1 ] ),
                                               space );
                
                if ( ! r )
                    HERROR( ERR_FMT_MTX, "(TMMVectorIO) read",
                            "missing data, expected <row> <col> <real> <imag>" );

                if ( col != 1 )
                    HERROR( ERR_FMT_MTX, "(TMMVectorIO) read",
                            "vector can have only 1 column" );
                    
                row--;
                col--;

                const complex  val( re, im );

                v->set_centry( row, val );
            }// for
        }// if
        else
        {
            for ( uint  i = 0; i < nrows; ++i )
            {
                std::getline( in, line );

                long        row = 0, col = 0;
                double      val = 0;
                const bool  r   = phrase_parse( line.begin(), line.end(),
                                                ( long_[ ref(row)  = _1 ] >>
                                                  long_[ ref(col)  = _1 ] >>
                                                  double_[ ref(val) = _1 ] ),
                                                space );
                
                if ( ! r )
                    HERROR( ERR_FMT_MTX, "(TMMVectorIO) read",
                            "missing data, expected <float>" );
                
                if ( col != 1 )
                    HERROR( ERR_FMT_MTX, "(TMMVectorIO) read",
                            "vector can have only 1 column" );
                    
                row--;
                col--;

                v->set_entry( row, real(val) );
            }// for
        }// else

        return unique_ptr< TVector >( v.release() );
    }// else
}

}// namespace HLIB
