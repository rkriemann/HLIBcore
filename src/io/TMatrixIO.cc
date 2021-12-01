//
// Project     : HLib
// File        : TMatrixIO.cc
// Description : classes for matrix input/output
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <string>
#include <sstream>
#include <vector>
#include <memory>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include "hlib-config.h"

#if USE_HDF5 == 1

#  if (defined(__GNUC__) || defined(__clang__)) && !defined(__ICC)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wold-style-cast"
#  endif

#include "H5Cpp.h"

#  if (defined(__GNUC__) || defined(__clang__)) && !defined(__ICC)
#    pragma GCC diagnostic pop
#  endif

#endif

#include "hpro/matrix/TSparseMatrix.hh"
#include "hpro/matrix/THMatrix.hh"
#include "hpro/matrix/structure.hh"

#include "baseio.hh"

#include "hpro/io/TMatrixIO.hh"

namespace HLIB
{

using std::string;
using std::vector;
using std::unique_ptr;
using std::make_unique;

namespace fs = boost::filesystem;

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
// class TMatrixIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

TMatrixIO::TMatrixIO ()
{
}

void
TMatrixIO::write ( const TMatrix *, const string & ) const
{
    HERROR( ERR_NOT_IMPL, "(TMatrixIO) write", "" );
}

unique_ptr< TMatrix >
TMatrixIO::read  ( const string & ) const
{
    HERROR( ERR_NOT_IMPL, "(TMatrixIO) write", "" );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// class TAutoMatrixIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void
TAutoMatrixIO::write ( const TMatrix *  A,
                       const string &   filename ) const
{
    unique_ptr< TMatrixIO >  mio;
        
    switch ( guess_format( filename ) )
    {
        case FMT_HLIB    : mio = make_unique< THLibMatrixIO >();   break;
        case FMT_MATLAB  : mio = make_unique< TMatlabMatrixIO >(); break;
        case FMT_SAMG    : mio = make_unique< TSAMGMatrixIO >();   break;
        case FMT_HB      : mio = make_unique< THBMatrixIO >();     break;
        case FMT_MTX     : mio = make_unique< TMMMatrixIO >();     break;
        case FMT_HDF5    : mio = make_unique< THDF5MatrixIO >();   break;

        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoMatrixIO) read", "" );
            return;
    }// switch

    if ( mio.get() == nullptr )
        return;

    mio->write( A, filename );
}

void
TAutoMatrixIO::write ( const TMatrix *  A,
                       const string &   filename,
                       const string &   matname ) const
{
    unique_ptr< TMatrixIO >  mio;
    bool                  support_named = false;
        
    switch ( guess_format( filename ) )
    {
        case FMT_HLIB    : mio = make_unique< THLibMatrixIO >();   break;
        case FMT_MATLAB  : mio = make_unique< TMatlabMatrixIO >(); support_named = true; break;
        case FMT_SAMG    : mio = make_unique< TSAMGMatrixIO >();   break;
        case FMT_HB      : mio = make_unique< THBMatrixIO >();     break;
        case FMT_MTX     : mio = make_unique< TMMMatrixIO >();     break;
        case FMT_HDF5    : mio = make_unique< THDF5MatrixIO >();   support_named = true; break;

        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoMatrixIO) read", "" );
            return;
    }// switch

    if ( mio.get() == nullptr )
        return;

    if ( support_named )
        ptrcast( mio.get(), TMatlabMatrixIO )->write( A, filename, matname );
    else
        mio->write( A, filename );
}

unique_ptr< TMatrix >
TAutoMatrixIO::read  ( const string &  filename ) const
{
    if ( filename == "" )
        HERROR( ERR_ARG, "(TAutoMatrixIO) read", "empty filename" );

    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TAutoMatrixIO) read", filename );
    
    unique_ptr< TMatrixIO >  mio;
        
    switch ( guess_format( filename ) )
    {
        case FMT_HLIB    : mio = make_unique< THLibMatrixIO >();   break;
        case FMT_MATLAB  : mio = make_unique< TMatlabMatrixIO >(); break;
        case FMT_SAMG    : mio = make_unique< TSAMGMatrixIO >();   break;
        case FMT_HB      : mio = make_unique< THBMatrixIO >();     break;
        case FMT_MTX     : mio = make_unique< TMMMatrixIO >();     break;
        case FMT_HDF5    : mio = make_unique< THDF5MatrixIO >();   break;

        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoMatrixIO) read", "" );
    }// switch

    return mio->read( filename );
}

unique_ptr< TMatrix >
TAutoMatrixIO::read  ( const string &  filename,
                       const string &  matname ) const
{
    if ( filename == "" )
        HERROR( ERR_ARG, "(TAutoMatrixIO) read", "empty filename" );

    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TAutoMatrixIO) read", filename );
    
    unique_ptr< TMatrixIO >  mio;
    bool                     support_named = false;
        
    switch ( guess_format( filename ) )
    {
        case FMT_HLIB    : mio = make_unique< THLibMatrixIO >();   break;
        case FMT_MATLAB  : mio = make_unique< TMatlabMatrixIO >(); support_named = true; break;
        case FMT_SAMG    : mio = make_unique< TSAMGMatrixIO >();   break;
        case FMT_HB      : mio = make_unique< THBMatrixIO >();     break;
        case FMT_MTX     : mio = make_unique< TMMMatrixIO >();     break;
        case FMT_HDF5    : mio = make_unique< THDF5MatrixIO >();   break;

        case FMT_UNKNOWN :
        default:
            HERROR( ERR_FMT_UNKNOWN, "(TAutoMatrixIO) read", "" );
    }// switch

    if ( support_named )
        return ptrcast( mio.get(), TMatlabMatrixIO )->read( filename, matname );
    else
        return mio->read( filename );
}

unique_ptr< TMatrix >
read_matrix  ( const std::string &  filename )
{
    TAutoMatrixIO  mio;

    return mio.read( filename );
}

unique_ptr< TMatrix >
read_matrix  ( const std::string &  filename,
               const std::string &  matname )
{
    TAutoMatrixIO  mio;

    return mio.read( filename, matname );
}

void
write_matrix  ( const TMatrix *      A,
                const std::string &  filename )
{
    TAutoMatrixIO  mio;

    mio.write( A, filename );
}

void
write_matrix  ( const TMatrix *      A,
                const std::string &  filename,
                const std::string &  matname )
{
    TAutoMatrixIO  mio;

    mio.write( A, filename, matname );
}

unique_ptr< TLinearOperator >
read_linop  ( const std::string &  filename )
{
    THLibMatrixIO  mio;

    return mio.read_linop( filename );
}

void
write_linop ( const TLinearOperator *  A,
              const std::string &      filename )
{
    THLibMatrixIO  mio;

    mio.write_linop( A, filename );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// class TOctaveMatrixIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void
TOctaveMatrixIO::write ( const TMatrix * A, const string & name ) const
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TOctaveMatrixIO) write", "argument is nullptr" );
    
    unique_ptr< std::ostream >  out_ptr( open_write( name ) );
    std::ostream &              out = * out_ptr.get();

    if ( ! is_dense( A ) )
        HERROR( ERR_MAT_TYPE, "(TOctaveMatrixIO) write", "only dense matrices supported" );
    
    //
    // convert matrix to dense and then write dense matrix
    //

    auto  D = cptrcast( A, TDenseMatrix );

    //
    // write header
    //

    const size_t  n = D->rows();
    const size_t  m = D->cols();

    out << "# Created by HLIBpro" << std::endl
        << "# name: A" << std::endl;
    if ( D->is_complex() ) out << "# type: complex matrix" << std::endl;
    else                   out << "# type: matrix" << std::endl;
    out << "# rows: " << n << std::endl
        << "# columns: " << m << std::endl;

    //
    // write matrix
    //

    if ( D->is_complex() )
    {
        for ( idx_t  i = 0; i < idx_t(n); i++ )
        {
            for ( idx_t  j = 0; j < idx_t(m); j++ )
            {
                const complex f = D->centry(i,j);
                
                out << boost::format( "(%.16e,%.16e) " ) % std::real(f) % std::imag(f);
            }// for
            
            out << std::endl;
        }// for
    }// if
    else
    {
        for ( idx_t  i = 0; i < idx_t(n); i++ )
        {
            for ( idx_t  j = 0; j < idx_t(m); j++ )
                out << boost::format( "%.16e " ) % D->entry(i,j);
            
            out << std::endl;
        }// for
    }// else
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class TSAMGMatrixIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// write matrix to file
//
void
TSAMGMatrixIO::write ( const TMatrix *  A,
                       const string &   filename ) const
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TSAMGMatrixIO) write", "matrix is nullptr" );

    // can not write complex matrices
    if ( A->is_complex() )
        HERROR( ERR_NREAL, "(TSAMGMatrixIO) write", "SAMG only supports real valued data" );

    //
    // set up name for format file
    //

    fs::path  filepath( filename );
    fs::path  ext = filepath.extension();
    string    compress_ext;
    string    frm_filename;

    while (( to_lower_copy( ext.string() ) == ".gz"    ) ||
           ( to_lower_copy( ext.string() ) == ".bzip2" ) ||
           ( to_lower_copy( ext.string() ) == ".zlib"  ))
    {
        compress_ext = compress_ext + ext.string();
        filepath.replace_extension( "" );
        ext = filepath.extension();
    }// while

    if ( to_lower_copy( ext.string() ) == ".amg" )
        filepath.replace_extension( "" );
    
    frm_filename = filepath.string() + ".frm" + compress_ext;
            
    //
    // write files
    //
    
    if ( IS_TYPE( A, TSparseMatrix ) )
    {
        const TSparseMatrix * S = cptrcast( A, TSparseMatrix );

        //
        // open format file to write basic information
        //
        
        const size_t    nnz = S->n_non_zero();
        const size_t    n   = S->rows();
        ulong           type;
        bool            row_sum_zero = true;
        vector< real >  diag_values( S->rows(), 0.0 );
        
        // determine, if row-sum is zero in matrix
        for ( idx_t  row = 0; row < idx_t(S->rows()); row++ )
        {
            const idx_t  lb = S->rowptr(row);
            const idx_t  ub = S->rowptr(row+1);
            real         f  = 0.0;
            bool         found_diag = false;
            
            for ( idx_t j = lb; j < ub; j++ )
            {
                if ( S->colind(j) == idx_t(row) )
                {
                    diag_values[row] = S->rcoeff(j);
                    found_diag       = true;
                }// if
                
                f += S->rcoeff(j);
            }// for

            if ( ! found_diag )
                HERROR( ERR_NOT_IMPL, "", "add zero diagonal elements" );
            
            if ( Math::abs( f ) > 1e-32 )
                row_sum_zero = false;
        }// for
        
        if ( row_sum_zero ) type = (S->is_symmetric() ? 11 : 21);
        else                type = (S->is_symmetric() ? 12 : 22);

        
        //
        // write format to frm file
        //
        
        unique_ptr< std::ostream >  out_ptr( open_write( frm_filename ) );
        std::ostream &              frm_out = * out_ptr.get();
        
        frm_out << "f  4" << std::endl
                << nnz << " " << n << " " << type << " 0 0" << std::endl;
    
        //
        // write matrix to amg file
        //
        
        out_ptr.reset( open_write( filename ) );

        std::ostream &  amg_out = * out_ptr.get();
        
        // first the row-index
        for ( uint row = 0; row <= n; row++ )
            amg_out << S->rowptr(row)+1 << std::endl;
        
        // then the column-index, but make sure that the diagonal entry comes first
        for ( uint row = 0; row < n; row++ )
        {
            const idx_t  lb = S->rowptr(row);
            const idx_t  ub = S->rowptr(row+1);
            
            // first the diagonal entry
            amg_out << row+1 << std::endl;
            
            // then the rest
            for ( idx_t  j = lb; j < ub; j++ )
            {
                if ( S->colind(j) != idx_t(row) )
                    amg_out << S->colind(j)+1 << std::endl;
            }// for
        }// for
        
        // finally the coefficients, again with the diag-entry first
        for ( uint row = 0; row < n; row++ )
        {
            const idx_t  lb = S->rowptr(row);
            const idx_t  ub = S->rowptr(row+1);
            
            amg_out << boost::format( "%.16E" ) % double(diag_values[row]) << std::endl;
            
            for ( idx_t  j = lb; j < ub; j++ )
            {
                if ( S->colind(j) != idx_t(row) )
                    amg_out << boost::format( "%.16E" ) % double(S->rcoeff(j)) << std::endl;
            }// for
        }// for
    }// if
    else if ( is_dense( A ) )
    {
        //
        // the matrix is converted to a dense matrix and then written
        // in sparse format
        //

        auto            D = cptrcast( A, TDenseMatrix );
        ulong           nnz;
        const size_t    n = D->rows();
        const size_t    m = D->cols();
        ulong           type;
        bool            row_sum_zero = true;
        vector< uint >  row_ptr( n+1, 0 );
        uint            pos  = 0;
        
        // determine, if row-sum is zero in matrix
        for ( idx_t  i = 0; i < idx_t(n); i++ )
        {
            real  f = 0.0;

            row_ptr[i] = pos;
            
            for ( idx_t  j = 0; j < idx_t(m); j++ )
            {
                const real val = D->entry( i, j );

                if (( i == j ) || ( val != 0.0 ))
                    pos++;
                
                f += val;
            }// for
            
            if ( Math::abs( f ) > 1e-32 )
                row_sum_zero = false;
        }// for

        row_ptr[n] = pos;
        nnz        = pos;
        
        if ( row_sum_zero ) type = (D->is_symmetric() ? 11 : 21);
        else                type = (D->is_symmetric() ? 12 : 22);
        
        //
        // write format to frm file
        //
        
        unique_ptr< std::ostream >  out_ptr( open_write( frm_filename ) );
        std::ostream &              frm_out = * out_ptr.get();
        
        frm_out << "f  4" << std::endl
                << nnz << " " << n << " " << type << " 0 0" << std::endl;
    
        //
        // write matrix to amg file
        //
        
        out_ptr.reset( open_write( filename ) );

        std::ostream &  amg_out = * out_ptr.get();
        
        // first the row-index
        for ( uint i = 0; i <= n; i++ )
            amg_out << row_ptr[i] + 1 << std::endl;
        
        // then the column-index, but make sure that the diagonal entry comes first
        for ( uint i = 0; i < n; i++ )
        {
            // first the diagonal entry
            amg_out << i+1 << std::endl;
            
            // then the rest
            for ( uint j = 0; j < m; j++ )
            {
                if (( j != i ) && ( D->entry( i, j ) != 0.0 ))
                    amg_out << j+1 << std::endl;
            }// for
        }// for
        
        // finally the coefficients, again with the diag-entry first
        for ( uint i = 0; i < n; i++ )
        {
            amg_out << boost::format( "%.16E" ) % double(D->entry( i, i )) << std::endl;
            
            for ( uint j = 0; j < m; j++ )
            {
                const double val = D->entry( i, j );
                
                if (( j != i ) && ( val != 0.0 ))
                    amg_out << boost::format( "%.16E" ) % val << std::endl;
            }// for
        }// for
    }// if
    else
        HERROR( ERR_MAT_TYPE, "(TSAMGMatrixIO) write", "unsupported matrix type " + A->typestr() );
}

//
// read matrix from file
//
unique_ptr< TMatrix >
TSAMGMatrixIO::read ( const string &  filename ) const
{
    //
    // open format file to read further information
    //

    fs::path  amgpath( filename );
    fs::path  ext = amgpath.extension();
    string    compress_ext;
    fs::path  frmpath;

    while (( to_lower_copy( ext.string() ) == ".gz"    ) ||
           ( to_lower_copy( ext.string() ) == ".bzip2" ) ||
           ( to_lower_copy( ext.string() ) == ".zlib"  ))
    {
        compress_ext = compress_ext + ext.string();
        amgpath.replace_extension( "" );
        ext = amgpath.extension();
    }// while
    
    if ( ext == ".amg" )
        frmpath = amgpath.replace_extension( ".frm" );
    else
        frmpath = amgpath.string() + ".frm";

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
    
    ulong             nnz, n, type;
    string            line;
    vector< string >  parts;

    std::getline( frm_in, line );
    std::getline( frm_in, line );
    split( line, " \t\r\n", parts );

    if ( parts.size() < 3 )
        HERROR( ERR_FMT_SAMG, "(TSAMGMatrixIO) read", "expected <nnz> <nrows> <type> in format file" );
    
    nnz  = str_to_int( parts[0] );
    n    = str_to_int( parts[1] );
    type = str_to_int( parts[2] );

    //
    // read matrix
    //

    in_ptr.reset( open_read( filename ) );

    std::istream &  amg_in = * in_ptr.get();

    auto  S = make_unique< TSparseMatrix >();

    S->set_size( n, n );
    S->init( nnz );

    // first fill row-indices
    for ( uint i = 0; i <= n; i++ )
    {
        ulong  idx;
            
        std::getline( amg_in, line );
        idx = str_to_int( line )-1;

        if ( idx > nnz )
            HERROR( ERR_FMT_SAMG, "(TSAMGMatrixIO) read", "row pointer index > number of non-zeroes" );
        
        S->rowptr(i) = idx;
    }// for

    // now fill column indices
    for ( uint i = 0; i < nnz; i++ )
    {
        ulong  idx;
            
        std::getline( amg_in, line );
        idx = str_to_int( line )-1;

        if ( idx >= n )
            HERROR( ERR_FMT_SAMG, "(TSAMGMatrixIO) read", "column index > number of columns" );
        
        S->colind(i) = idx;
    }// for

    // and finally the coefficients
    for ( uint i = 0; i < nnz; i++ )
    {
        double  val;
            
        std::getline( amg_in, line );
        val = str_to_dbl( line );

        S->rcoeff(i) = real(val);
    }// for

    S->sort_entries();

    if (( type == 11 ) || ( type == 12 )) S->set_symmetric();
    else                                  S->set_nonsym();

    return unique_ptr< TMatrix >( S.release() );
}

///////////////////////////////////////////////////
// 
// input and output in PLTMG format
//

//
// write matrix with name
//
void
TPLTMGMatrixIO::write ( const TMatrix *, const string & ) const
{
    HERROR( ERR_NOT_IMPL, "(TPLTMGMatrixIO) write", "" );
}

//
// read matrix from file
//
unique_ptr< TMatrix >
TPLTMGMatrixIO::read  ( const string & filename ) const
{
    //
    // read matrix file
    //

    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TPLTMGMatrixIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();
    string                      line;
    uint                        n;

    std::getline( in, line );
    n = str_to_int( line );
    std::getline( in, line );
    std::getline( in, line );

    //
    // read rhs
    //

    for ( uint i = 0; i < n; i++ )
        std::getline( in, line );

    //
    // read sparse matrix
    //

    auto  pos = in.tellg();
    uint  nnz = 0;
    
    vector< uint >    rows( n, 0 );
    vector< string >  parts;

    // first count number of entries per row to set up CRS structure
    while ( in.good() )
    {
        uint  i, j;
        
        std::getline( in, line );
        split( line, " \t\n", parts );
        if ( parts.size() < 2 )
            HERROR( ERR_FMT_PLTMG, "(TPLTMGMatrixIO) read", "expecting <row> <col>" );

        i = str_to_int( parts[0] );
        j = str_to_int( parts[1] );

        --i; --j;

        rows[i]++;
        nnz++;
    }// while

    // now reread file and fill CRS data
    auto   S   = make_unique< TSparseMatrix >( n, n );
    idx_t  ofs = 0;

    S->init( nnz );

    in.seekg( pos );

    for ( uint i = 0; i < n; i++ )
    {
        S->rowptr(i) = ofs;
        ofs += rows[i];
        rows[i] = 0; // reset for later usage
    }// for

    S->rowptr(n) = nnz;

    while ( in.good() )
    {
        uint    i, j;
        double  entry;
        
        std::getline( in, line );
        split( line, " \t\n", parts );
        if ( parts.size() < 3 )
            HERROR( ERR_FMT_PLTMG, "(TPLTMGMatrixIO) read", "expecting <row> <col> <coeff>" );

        i     = str_to_int( parts[0] );
        j     = str_to_int( parts[1] );
        entry = str_to_dbl( parts[2] );

        --i; --j;

        ofs = S->rowptr(i) + rows[i];

        S->colind(ofs) = j;
        S->rcoeff(ofs) = real(entry);
        
        rows[i]++;
    }// while
    
    return unique_ptr< TMatrix >( S.release() );
}

///////////////////////////////////////////////////
// 
// input and output in MatrixMarket format
//

namespace
{

enum mtxformat_t { MTX_ARRAY, MTX_COORD };
enum mtxfield_t  { MTX_REAL, MTX_INTEGER, MTX_COMPLEX, MTX_PATTERN };
enum mtxsym_t    { MTX_GENERAL, MTX_SYM, MTX_SKEWSYM, MTX_HERM };

}// namespace anonymous

//
// write matrix to <fname>
//
void
TMMMatrixIO::write ( const TMatrix *, const string & ) const
{
    HERROR( ERR_NOT_IMPL, "", "" );
}

//
// read matrix from <fname>
//
unique_ptr< TMatrix >
TMMMatrixIO::read  ( const string &  filename ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TMMMatrixIO) read", filename );
    
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

        if (( parts[0] == "%%MatrixMarket" ) || ( parts[0] == "%MatrixMarket" ))
        {
            if ( parts.size() < 3 )
                HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read",
                        "missing header data, expected <object> <format>" );
            
            if ( parts[1] != "matrix" )
                HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read", "file does not contain matrix" );
            
            if      ( parts[2] == "coordinate" ) mat_format = MTX_COORD;
            else if ( parts[2] == "array" )      mat_format = MTX_ARRAY;
            else
                HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read", "unknown matrix format \"" + parts[2] + "\"" );

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
                else if ( parts[i] != ""               )
                    HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read", "unknown qualifier \"" + parts[i] + "\"" );
            }// for
            
            break;
        }// if
    }// while

    if ( mat_sym == MTX_SKEWSYM )
        HERROR( ERR_NOT_IMPL, "(TMMMatrixIO) read", "skewsym matrices" );
    
    ///////////////////////////////////////////////////
    //
    // read content depending on format
    //

    if ( mat_format == MTX_COORD )
    {
        //
        // read number of rows, columns and non-zero entries
        //

        long  nrows = 0, ncols = 0, nnz = 0;

        while ( in.good() )
        {
            std::getline( in, line );

            if ( line[0] == '%' )
                continue;

            const bool  r = phrase_parse( line.begin(), line.end(),
                                          ( long_[ ref(nrows) = _1 ] >>
                                            long_[ ref(ncols) = _1 ] >>
                                            long_[ ref(nnz) = _1 ] ),
                                          space );

            if ( ! r )
                HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read", "missing data, expected <nrows> <ncols> <nnz>" );
            
            break;
        }// while

        HINFO( to_string( "(TMMMatrixIO) read : reading sparse %dx%d matrix with %d non-zero entries",
                          nrows, ncols, nnz ) );

        //
        // now, first go through file and count entries per row
        //

        const bool      is_sym   = (( mat_sym == MTX_SYM ) || ( mat_sym == MTX_HERM ));
        auto            file_pos = in.tellg();
        vector< uint >  row_count( nrows, 0 );
        size_t          real_nnz = 0;
        
        for ( long  i = 0; i < nnz; i++ )
        {
            std::getline( in, line );

            uint        row, col;
            const bool  r = phrase_parse( line.begin(), line.end(),
                                          ( long_[ ref(row) = _1 ] >>
                                            long_[ ref(col) = _1 ] ),
                                          space );

            if ( ! r )
                HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read", "missing data, expected <nrows> <ncols>" );
            
            row--;
            col--;
            
            if ( is_sym && ( row != col ))
            {
                row_count[col]++;
                real_nnz++;
            }// if
            
            row_count[row]++;
            real_nnz++;
        }// for

        //
        // finally, reread the file and construct CRS data for the
        // sparse matrix with the previously collected data
        //

        auto  S   = make_unique< TSparseMatrix >( nrows, ncols );
        uint  pos = 0;

        if ( mat_field == MTX_COMPLEX ) S->set_complex( true );

        S->init( real_nnz );

        // first the row pointers
        for ( long  i = 0; i < nrows; i++ )
        {
            S->rowptr(i) = pos;
            pos         += row_count[i];
            row_count[i] = 0; // reset for later usage
        }// for
        S->rowptr(nrows) = pos;

        if ( pos != real_nnz )
            HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read",
                    "number of non-zero entries and counted entries differs" );
        
        in.seekg( file_pos );
        
        for ( long  i = 0; i < nnz; i++ )
        {
            long  row = 0, col = 0;
            
            std::getline( in, line );

            if ( mat_field == MTX_COMPLEX )
            {
                double      re = 0, im = 0;
                const bool  r  = phrase_parse( line.begin(), line.end(),
                                               ( long_[ ref(row) = _1 ] >>
                                                 long_[ ref(col) = _1 ] >>
                                                 double_[ ref(re) = _1 ] >>
                                                 double_[ ref(im) = _1 ] ),
                                               space );
                
                if ( ! r )
                    HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read",
                            "missing data, expected <row> <cols> <real> <imag>" );

                // adjust since starts counting at 1
                row--;
                col--;
                
                const complex  val( re, im );

                {
                    const idx_t    idx = S->rowptr(row) + row_count[row];
                
                    S->colind(idx) = col;
                    S->ccoeff(idx) = val;
                }

                if ( is_sym && ( row != col ))
                {
                    const idx_t  idx = S->rowptr(col) + row_count[col];
                
                    S->colind(idx) = row;
                    S->ccoeff(idx) = val;
                }// if
            }// if
            else if ( mat_field == MTX_REAL )
            {
                double      val = 0;
                const bool  r   = phrase_parse( line.begin(), line.end(),
                                                ( long_[ ref(row) = _1 ] >>
                                                  long_[ ref(col) = _1 ] >>
                                                  double_[ ref(val) = _1 ] ),
                                                space );
                
                if ( ! r )
                    HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read",
                            "missing data, expected <row> <cols> <coeff>" );
                
                // adjust since starts counting at 1
                row--;
                col--;

                {
                    const idx_t  idx = S->rowptr(row) + row_count[row];
                
                    S->colind(idx) = col;
                    S->rcoeff(idx) = real(val);
                }

                if ( is_sym && ( row != col ))
                {
                    const idx_t  idx = S->rowptr(col) + row_count[col];
                
                    S->colind(idx) = row;
                    S->rcoeff(idx) = real(val);
                }// if
            }// if
            else if ( mat_field == MTX_PATTERN )
            {
                const double  val = 1;
                const bool    r   = phrase_parse( line.begin(), line.end(),
                                                  ( long_[ ref(row) = _1 ] >>
                                                    long_[ ref(col) = _1 ] ),
                                                  space );
                
                if ( ! r )
                    HERROR( ERR_FMT_MTX, "(TMMMatrixIO) read",
                            "missing data, expected <row> <cols> <coeff>" );
                
                // adjust since starts counting at 1
                row--;
                col--;

                {
                    const idx_t  idx = S->rowptr(row) + row_count[row];
                
                    S->colind(idx) = col;
                    S->rcoeff(idx) = real(val);
                }
                
                if ( is_sym && ( row != col ))
                {
                    const idx_t  idx = S->rowptr(col) + row_count[col];
                
                    S->colind(idx) = row;
                    S->rcoeff(idx) = real(val);
                }// if
            }// else
            
            row_count[row]++;
        }// for

        return unique_ptr< TMatrix >( S.release() );
    }// if
    else
        HERROR( ERR_NOT_IMPL, "(TMMMatrixIO) read", "array format" );
}

///////////////////////////////////////////////////
// 
// input and output in HDF5 format
//

#if USE_HDF5 == 1
    
using namespace H5;

namespace
{

template < typename value_t >
void
h5_write_dense ( H5File *                         file,
                 const std::string &              gname,
                 const std::string &              mname,
                 const BLAS::Matrix< value_t > &  A )
{
    if ( is_complex_type< value_t >::value )
        HERROR( ERR_CONSISTENCY, "", "" );
    
    hsize_t    dims[2] = { A.nrows(), A.ncols() };
    DataSpace  dataspace( 2, dims );

    //
    // define datatype
    //

    unique_ptr< FloatType >  datatype;
    
    if ( is_single_prec< value_t >::value )
        datatype = make_unique< FloatType >( PredType::NATIVE_FLOAT );
    else
        datatype = make_unique< FloatType >( PredType::NATIVE_DOUBLE );

    // // little endian
    // datatype->setOrder( H5T_ORDER_LE );
            
    //
    // Create dataset for matrix data
    //
            
    DataSet dataset = file->createDataSet( gname + "/" + mname, *datatype, dataspace );
            
    //
    // write the data to the dataset using default memory space, file
    // space, and transfer properties.
    //

    if ( is_single_prec< value_t >::value )
        dataset.write( A.data(), PredType::NATIVE_FLOAT );
    else
        dataset.write( A.data(), PredType::NATIVE_DOUBLE );
}

template < typename value_t >
void
h5_write_dense ( H5File *                                    file,
                 const std::string &                         gname,
                 const std::string &                         mname,
                 const BLAS::Matrix< Complex< value_t > > &  A )
{
    if ( is_complex_type< value_t >::value )
        HERROR( ERR_CONSISTENCY, "", "" );
    
    hsize_t    dims[2] = { A.nrows(), A.ncols() };
    DataSpace  dataspace( 2, dims );

    //
    // define datatype
    //

    CompType  datatype( sizeof( Complex< value_t > ) );
    
    if ( is_single_prec< value_t >::value )
    {
        datatype.insertMember( "real", 0,               PredType::NATIVE_FLOAT );
        datatype.insertMember( "imag", sizeof(value_t), PredType::NATIVE_FLOAT );
    }// if
    else
    {
        datatype.insertMember( "real", 0,               PredType::NATIVE_DOUBLE );
        datatype.insertMember( "imag", sizeof(value_t), PredType::NATIVE_DOUBLE );
    }// else
    
    //
    // Create dataset for matrix data
    //
            
    DataSet dataset = file->createDataSet( gname + "/" + mname, datatype, dataspace );
            
    //
    // write the data to the dataset using default memory space, file
    // space, and transfer properties.
    //

    if ( is_single_prec< value_t >::value )
        dataset.write( A.data(), datatype );
    else
        dataset.write( A.data(), datatype );
}

//
// write structural part of matrix (type, indexsets, etc.)
//
void
h5_write_struct ( H5File *             file,
                  const std::string &  gname,
                  const TMatrix *      M,
                  const std::string &  type )
{
    // local type for storing data
    struct  structure_t
    {
        char *  type;
        long    id;
        long    row_first, row_last;
        long    col_first, col_last;
    };

    structure_t  data;
    TIndexSet    row_is = M->row_is();
    TIndexSet    col_is = M->col_is();

    data.type      = const_cast< char * >( type.c_str() );
    data.id        = M->id();
    data.row_first = row_is.first();
    data.row_last  = row_is.last();
    data.col_first = col_is.first();
    data.col_last  = col_is.last();

    hid_t  str_type = H5Tcopy( H5T_C_S1 );
    
    H5Tset_size( str_type, H5T_VARIABLE );

    hsize_t    dims[] = { 1 };
    DataSpace  data_space( 1, dims );
    CompType   data_type( sizeof(structure_t) );
        
    data_type.insertMember( "type",      HOFFSET(structure_t, type),      str_type );
    data_type.insertMember( "id",        HOFFSET(structure_t, id),        PredType::NATIVE_LONG );
    data_type.insertMember( "row_first", HOFFSET(structure_t, row_first), PredType::NATIVE_LONG );
    data_type.insertMember( "row_last",  HOFFSET(structure_t, row_last),  PredType::NATIVE_LONG );
    data_type.insertMember( "col_first", HOFFSET(structure_t, col_first), PredType::NATIVE_LONG );
    data_type.insertMember( "col_last",  HOFFSET(structure_t, col_last),  PredType::NATIVE_LONG );

    auto  data_set = file->createDataSet( gname + "/structure", data_type, data_space );

    data_set.write( & data, data_type );

    data_space.close();
    data_type.close();
    data_set.close();
}

void
h5_write ( H5File *             file,
           const std::string &  gname,
           const TMatrix *      M );

void
h5_write ( H5File *              file,
           const std::string &   gname,
           const TDenseMatrix *  M )
{
    h5_write_struct( file, gname, M, "dense" );

    if ( M->is_complex() )
        h5_write_dense( file, gname, "D", M->blas_cmat() );
    else
        h5_write_dense( file, gname, "D", M->blas_rmat() );
}

void
h5_write ( H5File *              file,
           const std::string &   gname,
           const TRkMatrix *     M )
{
    h5_write_struct( file, gname, M, "lowrank" );

    if ( M->is_complex() )
    {
        h5_write_dense( file, gname, "U", M->blas_cmat_A() );
        h5_write_dense( file, gname, "V", M->blas_cmat_B() );
    }// if
    else
    {
        h5_write_dense( file, gname, "U", M->blas_rmat_A() );
        h5_write_dense( file, gname, "V", M->blas_rmat_B() );
    }// else
}

void
h5_write ( H5File *              file,
           const std::string &   gname,
           const TBlockMatrix *  M )
{
    h5_write_struct( file, gname, M, "structured" );

    //
    // write sub block structure and sub block ids
    //

    struct  structure_t
    {
        long    nblock_rows, nblock_cols;
        char *  ids;
    };

    std::ostringstream  ids;

    for ( uint i = 0; i < M->nblock_rows(); ++i )
    {
        for ( uint j = 0; j < M->nblock_cols(); ++j )
        {
            ids << M->block(i,j)->id() << ",";
        }// for
    }// for
                
    structure_t  data;

    data.nblock_rows = M->nblock_rows();
    data.nblock_cols = M->nblock_rows();
    data.ids         = const_cast< char * >( ids.str().c_str() );

    hid_t  str_type = H5Tcopy( H5T_C_S1 );
    
    H5Tset_size( str_type, H5T_VARIABLE );

    hsize_t    dims[] = { 1 };
    DataSpace  data_space( 1, dims );
    CompType   data_type( sizeof(structure_t) );
        
    data_type.insertMember( "nblock_rows", HOFFSET(structure_t, nblock_rows), PredType::NATIVE_LONG );
    data_type.insertMember( "nblock_cols", HOFFSET(structure_t, nblock_cols), PredType::NATIVE_LONG );
    data_type.insertMember( "ids",         HOFFSET(structure_t, ids),         str_type );

    auto  data_set = file->createDataSet( gname + "/blockstruct", data_type, data_space );

    data_set.write( & data, data_type );

    data_space.close();
    data_type.close();
    data_set.close();
    
    //
    // write sub blocks
    //
    
    auto  group = make_unique< Group >( file->createGroup( gname + "/subblocks" ) );

    // std::cout << gname + "/subblocks" << std::endl;
    
    for ( uint i = 0; i < M->nblock_rows(); ++i )
    {
        for ( uint j = 0; j < M->nblock_cols(); ++j )
        {
            auto  sgname = gname + "/subblocks/" + to_string( "%d", M->block(i,j)->id() );

            h5_write( file, sgname, M->block(i,j) );
        }// for
    }// for
}

void
h5_write ( H5File *             file,
           const std::string &  gname,
           const TMatrix *      M )
{
    // std::cout << gname << std::endl;
    
    auto  group = make_unique< Group >( file->createGroup( gname ) );
    
    if ( is_dense( M ) )
    {
        h5_write( file, gname, cptrcast( M, TDenseMatrix ) );
    }// if
    else if ( is_lowrank( M ) )
    {
        h5_write( file, gname, cptrcast( M, TRkMatrix ) );
    }// if
    else if ( is_blocked( M ) )
    {
        h5_write( file, gname, cptrcast( M, TBlockMatrix ) );
    }// if
    else
        HERROR( ERR_MAT_TYPE, "(THDF5MatrixIO) h5_write", M->typestr() );
}

template < typename value_t >
const BLAS::Matrix< value_t >
h5_read_dense ( H5File *             file,
                const std::string &  mname )
{
    if ( is_complex_type< value_t >::value )
        HERROR( ERR_CONSISTENCY, "", "" );

    auto  data_name = std::string( "" );
    auto  iter_op   = [] ( H5::H5Object &      loc,
                           const std::string   attr_name,
                           const H5O_info_t *  oinfo, 
                           void *              operator_data ) -> int
    {
        std::string *  dname = static_cast< std::string * >( operator_data );
        
        // if ( mname != "" )
        // {
        //     // use given matrix name if in file
        //     if ( mname == attr_name )
        //         *dname = attr_name;
        // }// if
        // else
        if ( attr_name != "." )
        {
            // use first name encountered
            if ( *dname == "" )
                *dname = attr_name;
            else
            {
                // reset if not expected value
                if (( attr_name != *dname + "/type" ) &&
                    ( attr_name != *dname + "/value" ))
                    *dname = "";
            }// else
        }// if
            
        // std::cout << attr_name << std::endl;
        
        return 0;
    };
    
    file->visit( H5_INDEX_NAME, H5_ITER_INC, iter_op, & data_name, 0 );

    // std::cout << "data : " << data_name << std::endl;
    
    auto  dataset    = file->openDataSet( data_name + "/value" );
    auto  type_class = dataset.getTypeClass();

    // if ( type_class == H5T_FLOAT )
    //     std::cout << "double" << std::endl;
    // else
    //     std::cout << type_class << std::endl;

    auto  dataspace = dataset.getSpace();
    auto  ndims     = dataspace.getSimpleExtentNdims();
    auto  dims      = std::vector< hsize_t >( ndims );
    
    dataspace.getSimpleExtentDims( dims.data() );

    if ( ndims != 2 )
        std::cout << "not a matrix" << std::endl;
                                                
    // std::cout << dims[0] << " × " << dims[1] << std::endl;

    auto  M = BLAS::Matrix< value_t >( dims[0], dims[1] );

    dataset.read( M.data(), PredType::NATIVE_DOUBLE );
    
    return  M;
}

}// namespace anonymous

#endif

//
// write matrix \a A to file \a fname
//
void
THDF5MatrixIO::write ( const TMatrix *      A,
                       const std::string &  fname ) const
{
    write( A, fname, "M" );
}

//
// write matrix \a A with name \a mname to file \a fname
//
#if USE_HDF5 == 1
void
THDF5MatrixIO::write ( const TMatrix *      A,
                       const std::string &  fname,
                       const std::string &  mname ) const
{
    try
    {
        H5File  file( fname, H5F_ACC_TRUNC );

        if ( A == nullptr )
            return;

        h5_write( & file, "/" + mname, A );
    }// try
    catch( FileIException &      error ) { error.printErrorStack(); }
    catch( DataSetIException &   error ) { error.printErrorStack(); }
    catch( DataSpaceIException & error ) { error.printErrorStack(); }
    catch( DataTypeIException &  error ) { error.printErrorStack(); }
}
#else
void
THDF5MatrixIO::write ( const TMatrix *      ,
                       const std::string &  ,
                       const std::string &   ) const
{
    HERROR( ERR_NOHDF5, "(THDF5MatrixIO) write", "" );
}
#endif


#if USE_HDF5 == 1
void
THDF5MatrixIO::write ( const BLAS::Matrix< real > &  A,
                       const std::string &           fname,
                       const std::string &           mname ) const
{
    try
    {
        H5File  file( fname, H5F_ACC_TRUNC );

        h5_write_dense( & file, "/", mname, A );
    }// try
    catch( FileIException &      error ) { error.printErrorStack(); }
    catch( DataSetIException &   error ) { error.printErrorStack(); }
    catch( DataSpaceIException & error ) { error.printErrorStack(); }
    catch( DataTypeIException &  error ) { error.printErrorStack(); }
}
#else
void
THDF5MatrixIO::write ( const BLAS::Matrix< real > &  ,
                       const std::string &           ,
                       const std::string &            ) const
{
    HERROR( ERR_NOHDF5, "(THDF5MatrixIO) write", "" );
}
#endif

#if USE_HDF5 == 1
void
THDF5MatrixIO::write ( const BLAS::Matrix< complex > &  A,
                       const std::string &              fname,
                       const std::string &              mname ) const
{
    try
    {
        H5File  file( fname, H5F_ACC_TRUNC );

        h5_write_dense( & file, "/", mname, A );
    }// try
    catch( FileIException &      error ) { error.printErrorStack(); }
    catch( DataSetIException &   error ) { error.printErrorStack(); }
    catch( DataSpaceIException & error ) { error.printErrorStack(); }
    catch( DataTypeIException &  error ) { error.printErrorStack(); }
}
#else
void
THDF5MatrixIO::write ( const BLAS::Matrix< complex > &  ,
                       const std::string &              ,
                       const std::string &               ) const
{
    HERROR( ERR_NOHDF5, "(THDF5MatrixIO) write", "" );
}
#endif

//
// read and return matrix from file \a fname
//
#if USE_HDF5 == 1
std::unique_ptr< TMatrix >
THDF5MatrixIO::read  ( const std::string &  filename ) const
{
    try
    {
        H5File  file( filename, H5F_ACC_RDONLY );

        auto  M = h5_read_dense< real >( & file, "" );

        return std::make_unique< TDenseMatrix >( is( 0, M.nrows()-1 ), is( 0, M.ncols()-1 ), std::move( M ) );
    }// try
    catch( FileIException &      error ) { error.printErrorStack(); }
    catch( DataSetIException &   error ) { error.printErrorStack(); }
    catch( DataSpaceIException & error ) { error.printErrorStack(); }
    catch( DataTypeIException &  error ) { error.printErrorStack(); }

    return std::unique_ptr< TMatrix >();
}
#endif

}// namespace HLIB
