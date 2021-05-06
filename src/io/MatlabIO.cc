//
// Project     : HLib
// File        : matlab_io.cc
// Description : Matlab related IO functions and classes
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <cstring>
#include <vector>
#include <sstream>
#include <memory>

#include "hlib-config.h"

#if HAS_BOOST_IOSTREAMS == 1
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#endif

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/matrix/TSparseMatrix.hh"
#include "hpro/matrix/THMatrix.hh"

#include "baseio.hh"

#include "hpro/io/TMatrixIO.hh"
#include "hpro/io/TVectorIO.hh"
#include "hpro/io/TCoordIO.hh"

namespace HLIB
{

namespace fs = boost::filesystem;

#if HAS_BOOST_IOSTREAMS == 1
namespace io = boost::iostreams;
#endif

using std::vector;
using std::unique_ptr;
using std::make_unique;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// local defines
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// padding for file acces in Matlab
#define PAD( s )             ((s) <= 4 ? (4 - (s)) : (((s)%8 == 0 ? 0 : 8 - ((s)%8))))
#define PAD_SEEKG( s, file ) { if (( (s) != 0 ) && ( PAD(s) != 0 )) { (file).seekg( PAD(s), std::ios_base::cur ); } }
#define PAD_SEEKP( s, file ) { if (( (s) != 0 ) && ( PAD(s) != 0 )) { (file).seekp( PAD(s), std::ios_base::cur ); } }

// for automatic swapping of bytes
#define BYTE_SWAP( t )  { if ( byteswap ) swap_bytes( t ); }

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// local types
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

namespace
{

//
// Matlab element and matrix types with their names
//
enum matmi_t
{
    MI_INT8 = 1,
    MI_UINT8,
    MI_INT16,
    MI_UINT16,
    MI_INT32,
    MI_UINT32,
    MI_SINGLE, 
    MI_RESERVED8,
    MI_DOUBLE,
    MI_RESERVED10,
    MI_RESERVED11,
    MI_INT64,
    MI_UINT64,
    MI_MATRIX,
    MI_COMPRESSED,
    MI_UTF8,
    MI_UTF16,
    MI_UTF32
};

enum matmx_t
{
    MX_CELL = 1,
    MX_STRUCTURE,
    MX_OBJECT,
    MX_CHAR_ARRAY,
    MX_SPARSE_ARRAY,
    MX_DBL_ARRAY,
    MX_SGL_ARRAY,
    MX_INT8_ARRAY,
    MX_UINT8_ARRAY,
    MX_INT16_ARRAY,
    MX_UINT16_ARRAY,
    MX_INT32_ARRAY,
    MX_UINT32_ARRAY
};

}// namespace anonymous

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// local functions
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

namespace
{

//
// template functions for reading and writing data
//
template <class T>
T
mat5_read ( std::istream &  in,
            const bool      byteswap )
{
    T  val;
    
    in.read( reinterpret_cast< char * >( & val ), sizeof(T) );
    BYTE_SWAP( val );

    return val;
}

template <class T>
void
mat5_read ( std::istream &  in,
            T *             data,
            const uint      nelem,
            const bool      byteswap )
{
    in.read( reinterpret_cast< char * >(  data ), sizeof(T)*nelem );
    
    if ( byteswap && ( sizeof(T) > 1 ))
    {
        for ( uint i = 0; i < nelem; i++ )
            BYTE_SWAP( data[i] );
    }// if
}

// template <>
// void
// mat5_read ( std::istream &  in, complex * data, const uint nelem, const bool byteswap )
// {
//     in.read( reinterpret_cast< char * >(  data ), sizeof(complex)*nelem );

//     if ( byteswap )
//     {
//         for ( uint i = 0; i < nelem; i++ )
//         {
//             BYTE_SWAP( std::real(data[i]) );
//             BYTE_SWAP( std::imag(data[i]) );
//         }// for
//     }// if
// }

template <class T>
void
mat5_write ( std::ostream &  out,
             T               data )
{
    out.write( reinterpret_cast< const char * >( & data ), sizeof(T) );
}

template <class T>
void
mat5_write ( std::ostream &  out,
             const T *       data,
             const size_t    nelem )
{
    out.write( reinterpret_cast< const char * >(  data ), sizeof(T)*nelem );
}

void
mat5_write ( std::ostream &       out,
             const std::string &  str )
{
    BOOST_FOREACH( int8_t  ch, str )
    {
        mat5_write< int8_t >( out, ch );
    }// for
}

//
// read element tag with type and size
//
void
mat5_readtag ( std::istream & in,
               matmi_t &      type,
               int32_t &      size,
               const bool     byteswap )
{
    int32_t  ltype = 0;

    //
    // read first part of tag with type information
    //
    
    in.read( reinterpret_cast< char * >( & ltype ), sizeof(ltype) );
    BYTE_SWAP( ltype );

    //
    // type might also contain size of "small elements"
    //
    
    if ( (( ltype >> 16 ) & 0xff) != 0 )
    {
        // small data element
        type = matmi_t( ltype & 0xff );
        size = ( ltype >> 16 ) & 0xff;
    }// if
    else
    {
        // normal data element
        type = matmi_t( ltype );
        
        in.read( reinterpret_cast< char * >( & size ), sizeof(size) );
        BYTE_SWAP( size );
    }// else
}

//
// write element tag with type and size
//
void
mat5_writetag ( std::ostream & out,
                const int32_t  atype,
                const int32_t  size )
{
    int32_t  type = atype;
    
    if ( size <= 4 )
    {
        // write small element tag
        type |= (size & 0xff) << 16;
        out.write( reinterpret_cast< const char * >( & type ), sizeof(type) );
    }// if
    else
    {
        // write standard element tag
        out.write( reinterpret_cast< const char * >( & type ), sizeof(type) );
        out.write( reinterpret_cast< const char * >( & size ), sizeof(size) );
    }// else
}

//
// read coefficient of array and return as double
//
real
mat5_read_coeff ( std::istream &  in,
                  const bool      byteswap,
                  const int32_t   type )
{
    real  coeff;

    switch ( type )
    {
        case MI_INT8   : coeff = real(mat5_read<int8_t>(   in, byteswap )); break;
        case MI_UINT8  : coeff = real(mat5_read<uint8_t>(  in, byteswap )); break;
        case MI_INT16  : coeff = real(mat5_read<int16_t>(  in, byteswap )); break;
        case MI_UINT16 : coeff = real(mat5_read<uint16_t>( in, byteswap )); break;
        case MI_INT32  : coeff = real(mat5_read<int32_t>(  in, byteswap )); break;
        case MI_UINT32 : coeff = real(mat5_read<uint32_t>( in, byteswap )); break;
        case MI_SINGLE : coeff = real(mat5_read<float>(    in, byteswap )); break;
        case MI_DOUBLE : coeff = real(mat5_read<double>(   in, byteswap )); break;
                        
    default:
        HERROR( ERR_NOT_IMPL, "mat5_read_coeff",
                to_string( "unsupported Matlab type %d", type ) );
    }// switch

    return coeff;
}

//
// read sparse double array
//
void
mat5_read_sparse ( std::istream &         in,
                   const bool             imag,
                   const bool             byteswap,
                   const vector< uint > & dims,
                   const uint             nonzero,
                   TMatrix **             M,
                   TVector **             v )
{
    unique_ptr< TSparseMatrix >  S;
    unique_ptr< TScalarVector >  x;

    if ( v != nullptr )
    {
        int  len = (dims[0] != 1 ? dims[0] : dims[1]);

        x.reset( new TScalarVector( len, 0, imag ) );
    }// if
    else
    {
        S.reset( new TSparseMatrix(  dims[0], dims[1] ) );
        S->set_complex( imag );
    }// else
    
    //
    // read row indices
    //

    matmi_t            type = matmi_t( 0 );
    int32_t            size = 0;
    size_t             nnz  = nonzero;
    vector< int32_t >  rowptr( nnz );

    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_INT32 )
        HERROR( ERR_FMT_MATLAB, "mat5_read_sparse", "expected INT32" );

    // adjust number of non-zero
    if ( size != int32_t( sizeof(int32_t) * nnz ) )
        nnz = size / sizeof(int32_t);
    
    for ( size_t  i = 0; i < nnz; i++ )
    {
        const int32_t row = mat5_read<int32_t>( in, byteswap );

        if ( row >= int32_t(dims[0]) )
            HERROR( ERR_NOT_IMPL, "mat5_read_sparse", "rectangular sparse matrices" );

        rowptr[i] = row;
    }// for
    
    PAD_SEEKG( size, in );
    
    //
    // read column indices
    //

    vector< int32_t >  colind( dims[1]+1 );

    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_INT32 )
        HERROR( ERR_FMT_MATLAB, "mat5_read_sparse", "expected INT32" );

    if ( size != int32_t( sizeof(int32_t) * (dims[1]+1) ) )
        HERROR( ERR_SIZE, "mat5_read_sparse", "while reading column entries" );
    
    for ( uint i = 0; i <= dims[1]; i++ )
    {
        const int32_t  col = mat5_read<int32_t>( in, byteswap );

        if ( col > int32_t( nnz ) )
            HERROR( ERR_CONSISTENCY, "mat5_read_sparse", "column index too large" );

        colind[i] = col;
    }// for
    
    PAD_SEEKG( size, in );

    // again, adjust number of non-zero with last number in colind
    if ( int32_t( nnz ) != colind[dims[1]] )
        nnz = colind[dims[1]];
    
    // build correct row indices for vectors in case of
    // column-vectors
    if (( x.get() != nullptr ) && ( dims[0] == 1 ))
    {
        uint  pos = 0;
        
        for ( uint i = 0; i < dims[1]; i++ )
            if ( colind[i] != colind[i+1] )
                rowptr[pos++] = i;
    }// if

    // build row-indices in CRS data of sparse matrices
    if ( S.get() != nullptr )
    {
        S->init( nnz );

        //
        // build row-wise numbering for sparse matrix
        //
        
        vector< int32_t >  tmp( dims[0], 0 );
        
        for ( uint j = 0; j < dims[1]; j++ )
        {
            const uint lb = colind[j];
            const uint ub = colind[j+1];

            for ( uint i = lb; i < ub; i++ )
                tmp[rowptr[i]]++;
        }// for

        //
        // build row indices in sparse matrix
        //

        uint  pos = 0;
        
        for ( uint i = 0; i < dims[0]; i++ )
        {
            S->rowptr(i) = pos;

            pos += tmp[i];
        }// for
        
        S->rowptr(dims[0]) = idx_t(nnz);
    }// if        
    
    //
    // read real parts
    //

    mat5_readtag( in, type, size, byteswap );

    if ( x.get() != nullptr )
    {
        for ( size_t  i = 0; i < nnz; i++ )
        {
            const real coeff = mat5_read_coeff( in, byteswap, type );

            if   ( imag ) x->set_centry( rowptr[i], coeff );
            else          x->set_entry(  rowptr[i], coeff );
        }// for
    }// if
    else
    {
        vector< int32_t >  tmp( dims[0], 0 );
        
        for ( uint col = 0; col < dims[1]; col++ )
        {
            const uint lb = colind[col];
            const uint ub = colind[col+1];

            for ( uint i = lb; i < ub; i++ )
            {
                const idx_t  row   = rowptr[i];
                const idx_t  idx   = S->rowptr(row) + tmp[row];
                const real   coeff = mat5_read_coeff( in, byteswap, type );
                
                S->colind(idx) = col;

                if ( imag ) S->ccoeff(idx) = coeff;
                else        S->rcoeff(idx) = coeff;

                tmp[row]++;
            }// for
        }// for
    }// else
    
    PAD_SEEKG( size, in );

    //
    // read imaginary parts
    //

    if ( imag )
    {
        mat5_readtag( in, type, size, byteswap );
        
        if ( x.get() != nullptr )
        {
            for ( size_t  i = 0; i < nnz; i++ )
            {
                const real coeff = mat5_read_coeff( in, byteswap, type );

                x->set_centry( rowptr[i], x->centry( rowptr[i] ) + complex( 0, coeff ) );
            }// for
        }// if
        else
        {
            vector< int32_t >  tmp( dims[0], 0 );
        
            for ( uint col = 0; col < dims[1]; col++ )
            {
                const uint lb = colind[col];
                const uint ub = colind[col+1];

                for ( uint i = lb; i < ub; i++ )
                {
                    const idx_t  row   = rowptr[i];
                    const idx_t  idx   = S->rowptr(row) + tmp[row];
                    const real   coeff = mat5_read_coeff( in, byteswap, type );
            
                    S->ccoeff(idx) += complex( 0, coeff );

                    tmp[row]++;
                }// for
            }// for
        }// else

        PAD_SEEKG( size, in );
    }// if

    //
    // finally sort entries in sparse matrices
    //

    if ( S.get() != nullptr )
    {
        S->sort_entries();
        S->test_symmetry();
    }// if

    //
    // return local objects or delete them
    //
    
    if ( x.get() != nullptr )
    {
        if ( *v == nullptr )
            *v = x.release();
    }// if
    else if ( S.get() != nullptr )
    {
        if ( *M == nullptr )
            *M = S.release();
    }// if
}

//
// read full double matrix
//
void
mat5_read_dblarray ( std::istream &         in,
                     const bool             imag,
                     const bool             byteswap,
                     const vector< uint > & dims,
                     TMatrix **             M,
                     TVector **             v )
{
    unique_ptr< TDenseMatrix >   D;
    unique_ptr< TScalarVector >  x;

    if ( v != nullptr )
    {
        const uint  len = (dims[0] != 1 ? dims[0] : dims[1]);

        x = std::make_unique< TScalarVector >( len, 0, imag );
    }// if
    else
    {
        D = std::make_unique< TDenseMatrix >();

        if ( D == nullptr )
            HERROR( ERR_MEM, "mat5_read_dblarray", "" );
        
        D->set_complex( imag );
        D->set_size( dims[0], dims[1] );
    }// else
    
    //
    // read real parts
    //

    matmi_t  type = matmi_t( 0 );
    int32_t  size = 0;

    mat5_readtag( in, type, size, byteswap );

    if ( x.get() != nullptr )
    {
        const uint  len = (dims[0] != 1 ? dims[0] : dims[1]);

        for ( uint i = 0; i < len; i++ )
        {
            const real  coeff = mat5_read_coeff( in, byteswap, type );
            
            if   ( imag ) x->set_centry( i, coeff );
            else          x->set_entry(  i, coeff );
        }// for
    }// if
    else
    {
        for ( uint j = 0; j < dims[1]; j++ )
            for ( uint i = 0; i < dims[0]; i++ )
            {
                const real  coeff = mat5_read_coeff( in, byteswap, type );
            
                if ( D->is_complex() ) D->set_centry( i, j, coeff );
                else                   D->set_entry(  i, j, coeff );
            }// for
    }// else
    
    PAD_SEEKG( size, in );

    //
    // read imaginary parts
    //

    if ( imag )
    {
        mat5_readtag( in, type, size, byteswap );
        
        if ( x.get() != nullptr )
        {
            const uint  len = (dims[0] != 1 ? dims[0] : dims[1]);
            
            for ( uint i = 0; i < len; i++ )
            {
                const real  coeff = mat5_read_coeff( in, byteswap, type );
            
                x->set_centry( i, x->centry( i ) + complex( 0, coeff ) );
            }// for
        }// if
        else
        {
            for ( uint j = 0; j < dims[1]; j++ )
                for ( uint i = 0; i < dims[0]; i++ )
                {
                    const real  coeff = mat5_read_coeff( in, byteswap, type );
                    
                    D->set_centry( i, j, D->centry( i, j ) + complex( 0, coeff ) );
                }// for
        }// else
        
        PAD_SEEKG( size, in );
    }// if

    //
    // test for symmetry or hermitian
    //

    if ( D.get() != nullptr )
    {
        if ( imag )
        {
            bool  is_sym  = true;
            bool  is_herm = true;
        
            for ( uint i = 0; i < dims[0]; i++ )
                for ( uint j = i+1; j < dims[1]; j++ )
                {
                    if ( D->centry( i, j ) != D->centry( j, i ) )
                        is_sym = false;

                    if ( D->centry( i, j ) != conj( D->centry( j, i ) ) )
                        is_herm = false;

                    if ( ! is_sym && ! is_herm )
                        break;
                }// for

            if ( is_sym )
                D->set_symmetric();
            else if ( is_herm )
                D->set_hermitian();
        }// if
        else
        {
            bool  is_sym = true;
        
            for ( uint i = 0; i < dims[0]; i++ )
                for ( uint j = i+1; j < dims[1]; j++ )
                    if ( D->entry( i, j ) != D->entry( j, i ) )
                    {
                        is_sym = false;
                        break;
                    }// if

            if ( is_sym )
                D->set_symmetric();
        }// else
    }// if

    //
    // return local objects or delete them
    //
    
    if ( x.get() != nullptr )
    {
        if ( *v == nullptr )
            *v = x.release();
    }// if
    else if ( D.get() != nullptr )
    {
        if ( *M == nullptr )
            *M = D.release();
    }// if
}

//
// read a single element from file
//
void
mat5_read_element ( std::istream &      in,
                    const bool          byteswap,
                    const std::string & name,
                    const bool          is_subfield,
                    const std::string & fieldname,
                    TMatrix **          M,
                    TVector **          v )
{
    int32_t  element_size = 0;
    matmi_t  type = matmi_t( 0 );
        
    mat5_readtag( in, type, element_size, byteswap );

    int32_t  size = element_size;
    auto     pos  = in.tellg();
        
    if ( type == MI_COMPRESSED )
    {
#if HAS_BOOST_IOSTREAMS == 1
        std::string  inbuf( size, '\0' );

        in.read( const_cast< char * >( inbuf.data() ), size );
        
        //
        // read first 8 bytes to obtain size
        //
        
        uint32_t  tmp[2] = { 0, 0 };

        {
            io::filtering_istream  zin;
            std::istringstream     strin( inbuf );

            zin.push( io::zlib_decompressor()  );
            zin.push( strin );

            zin.read( reinterpret_cast< char * >( tmp ), sizeof(tmp) );
        }

        //
        // uncompress all and continue reading matlab format
        // - can not uncompress on the fly using zlib_decompressor as it
        //   does not support tell/seek
        // - "+8" gives correct size
        //
        
        std::string  outbuf( tmp[1] + 8, '\0' );

        {
            io::filtering_istream  zin;
            std::istringstream     strin( inbuf );

            zin.push( io::zlib_decompressor()  );
            zin.push( strin );

            zin.read( const_cast< char * >( outbuf.data() ), tmp[1] + 8 );
        }
        
        std::istringstream  strout( outbuf );
        
        mat5_read_element( strout, byteswap, name, is_subfield, fieldname, M, v );

        // go to point after datastructure in case we did not
        // read everything
        in.seekg( pos );
        in.seekg( element_size, std::ios_base::cur );
        
        PAD_SEEKG( element_size, in );
                
        return;

#else
        
        HERROR( ERR_NOT_IMPL, "mat5_read_element", "nompression not supported" );
        
#endif
    }// if
    else if ( type != MI_MATRIX )
        HERROR( ERR_NOT_IMPL, "mat5_read_element", "only miMatrix is supported" );

    //
    // read array flags
    //

    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_UINT32 )
        HERROR( ERR_FMT_MATLAB, "mat5_read_element", "expected UINT32" );

    if ( size != 8 )
        HERROR( ERR_SIZE, "mat5_read_element", "expected 2 * UINT32" );

    const uint32_t flags  = mat5_read<uint32_t>( in, byteswap );
    const bool     imag   = ((flags & 0x0800) != 0);
    const matmx_t  aclass = matmx_t( flags & 0xff );

    PAD_SEEKG( size, in ); 
    
    //
    // read number of non-zeroes in case of sparse matrix
    // or junk in case of a dense matrix
    //
    
    const uint32_t  nnz = mat5_read<uint32_t>( in, byteswap );
    PAD_SEEKG( size, in );
    
    //
    // read dimensions 
    //
    
    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_INT32 )
        HERROR( ERR_FMT_MATLAB, "mat5_read_dblarray", "expected INT32" );

    const uint     ndims    = size / 4;
    vector< uint > dims( ndims );
    uint           nentries = 1;
    bool           skip     = false;  // indicates skipping of data element

    for ( uint i = 0; i < ndims; i++ )
    {
        dims[i]   = mat5_read<int32_t>( in, byteswap );
        nentries *= dims[i];
    }// for
    
    PAD_SEEKG( size, in ); 

    if ( ndims != 2 )
        HERROR( ERR_NOT_IMPL, "mat5_read_dblarray", "only 2-dim. matrices supported" );

    if ( aclass != MX_STRUCTURE )
    {
        // check if vector is sought but matrix found
        if (( M == nullptr ) && ( v != nullptr ) && ( dims[0] > 1 ) && ( dims[1] > 1 ))
            skip = true;
        
        // check if matrix is sought but vector found
        if (( M != nullptr ) && ( v == nullptr ) && (( dims[0] == 1 ) || ( dims[1] == 1 )))
            skip = true;
    }// if
    
    //
    // read name of array
    //
    
    mat5_readtag( in, type, size, byteswap );

    if ( type != MI_INT8 )
        HERROR( ERR_FMT_MATLAB, "mat5_read_dblarray", "expected INT8" );

    std::string  element_name( size_t(size+1), '\0' );

    element_name[size] = '\0';
    
    in.read( const_cast< char * >( element_name.data() ), size );
    PAD_SEEKG( size, in );

    // in case of a subfield: adjust element name
    if (( element_name == "" ) && is_subfield )
        element_name = fieldname;

    // adjust "skipping" if requested name equals element_name
    if (( name != "" ) && ( name == element_name ))
        skip = false;
        
    //
    // read fields based on class value
    //

    if ( ! skip )
    {
        if ( aclass == MX_STRUCTURE )
        {
            mat5_readtag( in, type, size, byteswap );

            if ( type != MI_INT32 )
                HERROR( ERR_FMT_MATLAB, "mat5_read_element", "expected INT32" );

            const int32_t  fnlen = mat5_read<int32_t>( in, byteswap );
        
            mat5_readtag( in, type, size, byteswap );

            if ( type != MI_INT8 )
                HERROR( ERR_FMT_MATLAB, "mat5_read_element", "expected INT8" );

            const uint        fields = size / fnlen;
            vector< std::string >  fieldnames( fields );

            for ( uint i = 0; i < fields; i++ )
            {
                fieldnames[i].resize( fnlen+1 );
                fieldnames[i][fnlen] = '\0';
                in.read( const_cast< char * >( fieldnames[i].data() ), fnlen );
            }// for

            PAD_SEEKG( size, in );

            for ( uint i = 0; i < fields; i++ )
            {
                mat5_read_element( in, byteswap, name, true, fieldnames[i], M, v );

                if (( M != nullptr ) && ( *M != nullptr ))
                    break;
                else if (( v != nullptr ) && ( *v != nullptr ))
                    break;
            }// for
        }// if
        else if ( aclass == MX_SPARSE_ARRAY )
        {
            if (( name == "" ) || ( name == element_name ))
                mat5_read_sparse( in, imag, byteswap, dims, nnz, M, v );
        }// if
        else if ( aclass == MX_DBL_ARRAY )
        {
            if (( name == "" ) || ( name == element_name ))
                mat5_read_dblarray( in, imag, byteswap, dims, M, v );
        }// if
        else if ( aclass == MX_SGL_ARRAY )
        {
            if (( name == "" ) || ( name == element_name ))
                mat5_read_dblarray( in, imag, byteswap, dims, M, v );
        }// if
        // else
        //     HERROR( ERR_TYPE, "mat5_read_element", "only dense and sparse matrices supported" );
    }// if

    // go to point after datastructure in case we did not
    // read everything
    in.seekg( pos );
    in.seekg( element_size, std::ios_base::cur );
        
    PAD_SEEKG( element_size, in ); 
}

//
// special write methods
//
template <typename T>
void
write_dense ( std::ostream &            out,
              const BLAS::Matrix< T > & D,
              const std::string &       matname )
{
    using  real_t = typename real_type< T >::type_t;
    
    // remember current position for size info
    const auto  pos = out.tellp();
        
    // size = 5 ensures standard element tag
    mat5_writetag( out, MI_MATRIX, 5 );

    //
    // write array flags
    // (if complex valued)
    //

    uint32_t  flags = ( is_complex_type<T>::value ? 0x0800 : 0 );
    uint32_t  tmp   = 0;

    flags |= is_single_prec< T >::value ? MX_SGL_ARRAY : MX_DBL_ARRAY;
    
    mat5_writetag( out, MI_UINT32, 8 );
    mat5_write<uint32_t>( out, flags );
    mat5_write<uint32_t>( out, tmp   );

    //
    // write dimensions
    //

    int32_t  rows = int32_t( D.nrows() );
    int32_t  cols = int32_t( D.ncols() );

    mat5_writetag( out, MI_INT32, 8 );
    mat5_write<int32_t>( out, rows );
    mat5_write<int32_t>( out, cols );

    //
    // write name
    //

    int32_t  size = int32_t( matname.length()*sizeof(int8_t) );
    
    mat5_writetag( out, MI_INT8, size );
    mat5_write( out, matname );

    PAD_SEEKP( size, out );
    
    //
    // write real part of array
    // - no padding neccessary since sizeof(double) == 8
    //

    matmi_t  mi_type = is_single_prec< T >::value ? MI_SINGLE : MI_DOUBLE;
    
    mat5_writetag( out, mi_type, int32_t( sizeof(real_t)*rows*cols ) );

    for ( int j = 0; j < cols; j++ )
        for ( int i = 0; i < rows; i++ )
        {
            mat5_write< real_t >( out, std::real( D(i,j) ) );
        }// for
    
    //
    // write imaginary part of array
    //

    if ( is_complex_type<T>::value )
    {
        mat5_writetag( out, mi_type, int32_t( sizeof(real_t)*rows*cols ) );

        for ( int j = 0; j < cols; j++ )
            for ( int i = 0; i < rows; i++ )
            {
                mat5_write< real_t >( out, std::imag( D(i,j) ) );
            }// for
    }// else

    //
    // compute size of object and rewrite array tag
    // (remove size of tag)
    //

    size = int32_t( out.tellp() - pos - 8 );

    out.seekp( pos );
        
    mat5_writetag( out, MI_MATRIX, size );
}

void
write_rank ( std::ostream &       out,
             const TRkMatrix *    R,
             const std::string &  matname )
{
    //
    // split low-rank matrix in A and B and write matname_A and matname_B
    // first : matrix A
    //
    
    // remember current position for size info
    const auto  pos = out.tellp();
        
    // size = 5 ensures standard element tag
    mat5_writetag( out, MI_MATRIX, 5 );

    //
    // write array flags
    // (if complex valued)
    //

    uint32_t  flags = ( R->is_complex() ? 0x0800 : 0 );
    uint32_t  tmp   = 0;

    flags |= MX_DBL_ARRAY;
    
    mat5_writetag( out, MI_UINT32, 8 );
    mat5_write<uint32_t>( out, flags );
    mat5_write<uint32_t>( out, tmp   );

    //
    // write dimensions
    //

    int32_t  rows = int32_t(R->rows());
    int32_t  rank = int32_t(R->rank());

    mat5_writetag( out, MI_INT32, 8 );
    mat5_write<int32_t>( out, rows );
    mat5_write<int32_t>( out, rank );

    //
    // write name
    //

    std::string   tname = matname + "_A";
    int32_t       size  = int(tname.length()*sizeof(int8_t));
    
    mat5_writetag( out, MI_INT8, size );
    mat5_write( out, tname );

    PAD_SEEKP( size, out );
    
    //
    // write real part of array
    // - no padding neccessary since sizeof(double) == 8
    //

    mat5_writetag( out, MI_DOUBLE, int32_t(sizeof(double)*rows*rank) );

    if ( R->is_complex() )
    {
        for ( int k = 0; k < rank; k++ )
            for ( int i = 0; i < rows; i++ )
            {
                mat5_write<double>( out, std::real( R->blas_cmat_A()( i, k ) ) );
            }// for
    }// if
    else
    {
        for ( int k = 0; k < rank; k++ )
            for ( int i = 0; i < rows; i++ )
            {
                mat5_write<double>( out, R->blas_rmat_A()( i, k ) );
            }// for
    }// else
    
    //
    // write imaginary part of array
    //

    if ( R->is_complex() )
    {
        mat5_writetag( out, MI_DOUBLE, int32_t(sizeof(double)*rows*rank) );

        for ( int k = 0; k < rank; k++ )
            for ( int i = 0; i < rows; i++ )
            {
                mat5_write<double>( out, std::imag( R->blas_cmat_A()( i, k ) ) );
            }// for
    }// if

    // compute size of object and rewrite array tag
    // (remove size of tag)
    auto  fpos = out.tellp();

    size = int32_t( fpos - pos - 8 );

    out.seekp( pos );
    mat5_writetag( out, MI_MATRIX, size );
    out.seekp( fpos );
    
    //
    // now write matrix B
    //

    // remember current position for size info
    fpos = out.tellp();
        
    // size = 5 ensures standard element tag
    mat5_writetag( out, MI_MATRIX, 5 );

    //
    // write array flags
    // (if complex valued)
    //

    tmp = 0;

    mat5_writetag( out, MI_UINT32, 8 );
    mat5_write<uint32_t>( out, flags );
    mat5_write<uint32_t>( out, tmp   );

    //
    // write dimensions
    //

    int32_t  cols = int32_t(R->cols());

    mat5_writetag( out, MI_INT32, 8 );
    mat5_write<int32_t>( out, cols );
    mat5_write<int32_t>( out, rank );

    //
    // write name
    //

    tname = matname + "_B";
    size  = int32_t(tname.length()*sizeof(int8_t));
    
    mat5_writetag( out, MI_INT8, size );
    mat5_write( out, tname );

    PAD_SEEKP( size, out );
    
    //
    // write real part of array
    // - no padding neccessary since sizeof(double) == 8
    //

    mat5_writetag( out, MI_DOUBLE, int32_t(sizeof(double)*cols*rank) );

    if ( R->is_complex() )
    {
        for ( int k = 0; k < rank; k++ )
            for ( int i = 0; i < cols; i++ )
            {
                mat5_write<double>( out, std::real( R->blas_cmat_B()( i, k ) ) );
            }// for
    }// if
    else
    {
        for ( int k = 0; k < rank; k++ )
            for ( int i = 0; i < cols; i++ )
            {
                mat5_write<double>( out, R->blas_rmat_B()( i, k ) );
            }// for
    }// else
    
    //
    // write imaginary part of array
    //

    if ( R->is_complex() )
    {
        mat5_writetag( out, MI_DOUBLE, int32_t(sizeof(double)*cols*rank) );

        for ( int k = 0; k < rank; k++ )
            for ( int i = 0; i < cols; i++ )
            {
                mat5_write<double>( out, std::imag( R->blas_cmat_B()( i, k ) ) );
            }// for
    }// if

    // compute size of object and rewrite array tag
    // (remove size of tag)
    size = int32_t( out.tellp() - fpos - 8);

    out.seekp( fpos );
    mat5_writetag( out, MI_MATRIX, size );
}

void
write_sparse ( std::ostream &         out,
               const TSparseMatrix *  S,
               const std::string &    matname )
{
    // remember current position for size info
    auto  fpos = out.tellp();
        
    // size = 5 ensures standard element tag
    mat5_writetag( out, MI_MATRIX, 5 );
    
    //
    // write array flags
    // (if complex valued)
    //

    uint32_t  flags = ( S->is_complex() ? 0x0800 : 0 );
    uint32_t  nnz   = uint32_t(S->n_non_zero());

    flags |= MX_SPARSE_ARRAY;
    
    mat5_writetag( out, MI_UINT32, 8 );
    mat5_write<uint32_t>( out, flags );
    mat5_write<uint32_t>( out, nnz   );

    //
    // write dimensions
    //

    int32_t  rows = int32_t(S->rows());
    int32_t  cols = int32_t(S->cols());

    mat5_writetag( out, MI_INT32, 8 );
    mat5_write<int32_t>( out, rows );
    mat5_write<int32_t>( out, cols );

    //
    // write name
    //

    int32_t  size = int32_t(matname.length()*sizeof(int8_t));
    
    mat5_writetag( out, MI_INT8, size );
    mat5_write( out, matname );

    PAD_SEEKP( size, out );
    
    //
    // build column-wise numbering for sparse matrix
    //

    vector< int32_t >  tmp( S->cols(), 0 );
    
    for ( int i = 0; i < rows; i++ )
    {
        const idx_t  lb = S->rowptr(i);
        const idx_t  ub = S->rowptr(i+1);
        
        for ( idx_t  j = lb; j < ub; j++ )
            tmp[ S->colind(j) ]++;
    }// for

    //
    // build column indices in CCS format
    //

    vector< int32_t >  ccs_colind( cols+1 );
    int                pos = 0;
    
    for ( int i = 0; i < cols; i++ )
    {
        ccs_colind[i] = pos;
        pos           += tmp[i];
    }// for
    
    ccs_colind[cols] = nnz;
        
    //
    // build row indices in CCS format
    //

    vector< int32_t >  ccs_rowptr( nnz );

    for ( int i = 0; i < cols; i++ )
        tmp[i] = 0;

    for ( int i = 0; i < rows; i++ )
    {
        const idx_t  lb = S->rowptr(i);
        const idx_t  ub = S->rowptr(i+1);

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col = S->colind(j);
            const idx_t  idx = ccs_colind[col] + tmp[col];
                
            ccs_rowptr[ idx ] = i;
            tmp[col]++;
        }// for
    }// for
    
    //
    // build new coeff array in CCS format
    //

    vector< double >  ccs_coeff( nnz );

    for ( int i = 0; i < cols; i++ )
        tmp[i] = 0;
    
    for ( int i = 0; i < rows; i++ )
    {
        const idx_t  lb = S->rowptr(i);
        const idx_t  ub = S->rowptr(i+1);

        for ( idx_t  j = lb; j < ub; j++ )
        {
            const idx_t  col = S->colind(j);
            const idx_t  idx = ccs_colind[col] + tmp[col];

            if ( S->is_complex() )
                ccs_coeff[ idx ] = std::real( S->ccoeff(j) );
            else
                ccs_coeff[ idx ] = S->rcoeff(j);
            
            tmp[col]++;
        }// for
    }// for
    
    //
    // write row indices
    //

    mat5_writetag( out, MI_INT32, int32_t( sizeof(int32_t)*nnz ) );
    mat5_write<int32_t>( out, ccs_rowptr.data(), nnz );
    PAD_SEEKP( sizeof(int32_t)*nnz, out );
    
    //
    // write column indices
    //

    mat5_writetag( out, MI_INT32, int32_t( sizeof(int32_t)*(cols+1) ) );
    mat5_write<int32_t>( out, ccs_colind.data(), cols+1 );
    PAD_SEEKP( sizeof(int32_t)*(cols+1), out );

    //
    // write real parts
    //

    mat5_writetag( out, MI_DOUBLE, int32_t( sizeof(double)*nnz ) );

    for ( uint i = 0; i < nnz; i++ )
        mat5_write<double>( out, ccs_coeff[i] );

    //
    // read imaginary parts
    //

    if ( S->is_complex() )
    {
        for ( int i = 0; i < cols; i++ )
            tmp[i] = 0;
        
        for ( uint i = 0; i < nnz; i++ )
            ccs_coeff[i] = 0.0;
        
        for ( int i = 0; i < rows; i++ )
        {
            const idx_t  lb = S->rowptr(i);
            const idx_t  ub = S->rowptr(i+1);
            
            for ( idx_t  j = lb; j < ub; j++ )
            {
                const idx_t  col = S->colind(j);
                const idx_t  idx = ccs_colind[col] + tmp[col];
                
                ccs_coeff[ idx ] = std::imag( S->ccoeff(j) );
                
                tmp[col]++;
            }// for
        }// for
    
        mat5_writetag( out, MI_DOUBLE, int32_t( sizeof(double)*nnz ) );
        mat5_write<double>( out, ccs_coeff.data(), nnz );
    }// if
    
    //
    // compute size of object and rewrite array tag
    // (remove size of tag)
    //

    size = int32_t( long(out.tellp()) - fpos - 8 );

    out.seekp( fpos );
        
    mat5_writetag( out, MI_MATRIX, size );
}

template <typename T>
void
write_scalar ( std::ostream &            out,
               const BLAS::Vector< T > & v,
               const std::string &       vecname )
{
    // remember current position for size info
    auto  pos = out.tellp();
        
    // size = 5 ensures standard element tag
    mat5_writetag( out, MI_MATRIX, 5 );

    //
    // write array flags
    // (if complex valued)
    //

    uint32_t  flags = ( is_complex_type<T>::value ? 0x0800 : 0 );
    uint32_t  tmp   = 0;

    flags |= MX_DBL_ARRAY;
    
    mat5_writetag( out, MI_UINT32, 8 );
    mat5_write<uint32_t>( out, flags );
    mat5_write<uint32_t>( out, tmp   );

    //
    // write dimensions
    //

    int32_t  dim = int32_t( v.length() );

    mat5_writetag( out, MI_INT32, 8 );
    mat5_write<int32_t>( out, dim );
    mat5_write<int32_t>( out, 1 );

    //
    // write name
    //

    int32_t  size = int32_t( vecname.length()*sizeof(int8_t) );
    
    mat5_writetag( out, MI_INT8, size );
    mat5_write( out, vecname );

    PAD_SEEKP( size, out );
    
    //
    // write real part of array
    // - no padding neccessary since sizeof(double) == 8
    //

    mat5_writetag( out, MI_DOUBLE, int32_t( sizeof(double)*dim ) );

    for ( int i = 0; i < dim; i++ )
        mat5_write<double>( out, std::real( v(i) ) );
    
    //
    // write imaginary part of array
    //

    if ( is_complex_type<T>::value )
    {
        mat5_writetag( out, MI_DOUBLE, int32_t( sizeof(double)*dim ) );

        for ( int i = 0; i < dim; i++ )
            mat5_write<double>( out, std::imag( v(i) ) );
    }// else

    //
    // compute size of object and rewrite array tag
    // (remove size of tag)
    //

    size = int32_t( out.tellp() - pos - 8 );

    out.seekp( pos );
        
    mat5_writetag( out, MI_MATRIX, size );
}

}// namespace anonymous

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class TMatlabMatrixIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// write matrix to file
//
void
TMatlabMatrixIO::write ( const TMatrix *      A,
                         const std::string &  filename ) const
{
    return write( A, filename, "M" );
}

//
// write matrix with given name (options are obsolete here)
//
void
TMatlabMatrixIO::write ( const TMatrix *      A,
                         const std::string &  filename,
                         const std::string &  matname ) const
{
    if ( A == nullptr )
        HERROR( ERR_ARG, "(TMatlabMatrixIO) write", "matrix is NULL" );

    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();

    int16_t  version   = 0x0100;
    int16_t  endian    = 0x4d49; // = MI
    char     buf[116]  = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    if ( IS_TYPE( A, TDenseMatrix ) )
    {
        if ( A->is_complex() )
            write_dense<complex>( out, cptrcast( A, TDenseMatrix )->blas_cmat(), matname );
        else
            write_dense<real>(    out, cptrcast( A, TDenseMatrix )->blas_rmat(), matname );
    }// if
    else if ( IS_TYPE( A, TRkMatrix ) )
    {
        write_rank( out, cptrcast( A, TRkMatrix ), matname );
    }// if
    else if ( IS_TYPE( A, TSparseMatrix ) )
    {
        write_sparse( out, cptrcast( A, TSparseMatrix ), matname );
    }// if
    else
        HERROR( ERR_MAT_TYPE, "(TMatlabMatrixIO) write", "unsupported matrix type " + A->typestr() );
}

//
// write linear operator \a A with name \a mname to file \a fname
//
void
TMatlabMatrixIO::write ( const TLinearOperator *  A,
                         const std::string &      /* fname */,
                         const std::string &      /* mname */ ) const
{
    HERROR( ERR_MAT_TYPE, "(TMatlabMatrixIO) write", A->typestr() );
}

//
// write BLAS matrix with given name
//
void
TMatlabMatrixIO::write ( const BLAS::Matrix< float > &  A,
                         const std::string &            filename,
                         const std::string &            matname ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();

    int16_t  version   = 0x0100;
    int16_t  endian    = 0x4d49; // = MI
    char     buf[116]  = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    write_dense< float >( out, A, matname );
}

void
TMatlabMatrixIO::write ( const BLAS::Matrix< double > &  A,
                         const std::string &             filename,
                         const std::string &             matname ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();

    int16_t  version   = 0x0100;
    int16_t  endian    = 0x4d49; // = MI
    char     buf[116]  = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    write_dense< double >( out, A, matname );
}

void
TMatlabMatrixIO::write ( const BLAS::Matrix< Complex< float > > &  A,
                         const std::string &                        filename,
                         const std::string &                        matname ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();
    
    int16_t  version   = 0x0100;
    int16_t  endian    = 0x4d49; // = MI
    char     buf[116]  = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    write_dense< Complex< float > >( out, A, matname );
}

void
TMatlabMatrixIO::write ( const BLAS::Matrix< Complex< double > > &  A,
                         const std::string &                        filename,
                         const std::string &                        matname ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();
    
    int16_t  version   = 0x0100;
    int16_t  endian    = 0x4d49; // = MI
    char     buf[116]  = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    write_dense< Complex< double > >( out, A, matname );
}

//
// read matrix from file (assuming only one entry available)
//
unique_ptr< TMatrix >
TMatlabMatrixIO::read ( const std::string & filename ) const
{
    return read( filename, "" );
}

//
// read matrix with name <matname> from file
//
unique_ptr< TMatrix >
TMatlabMatrixIO::read  ( const std::string & filename,
                         const std::string & matname ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TMatlabMatrixIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    int16_t    version, endian;
    bool       byteswap = false;
    TMatrix *  M = nullptr;

    in.seekg( 124 );
    in.read( reinterpret_cast< char * >( & version ), sizeof(version) );
    in.read( reinterpret_cast< char * >( & endian ),  sizeof(endian)  );
    
    if      ( endian == 0x4d49 )  byteswap = false; // MI
    else if ( endian == 0x494d )  byteswap = true;  // IM
    else
        HERROR( ERR_FMT_MATLAB, "(TMatlabMatrixIO) read",
                "endianess indicator wrong" );

    BYTE_SWAP( version );
        
    while (( M == nullptr ) && in.good() )
        mat5_read_element( in, byteswap, matname, false, "", & M, nullptr );
    
    return unique_ptr< TMatrix >( M );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class TMatlabVectorIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// write vector with name "v" to file <filename>
//
void
TMatlabVectorIO::write ( const TVector *      v,
                         const std::string &  filename ) const
{
    write( v, filename, "v" );
}

//
// write vector with name <vecname> to file <filename>
//
void
TMatlabVectorIO::write ( const TVector *      v,
                         const std::string &  filename,
                         const std::string &  vecname ) const
{
    if ( v == nullptr )
        HERROR( ERR_ARG, "(TMatlabVectorIO) write", "vector is NULL" );

    if ( ! IS_TYPE( v, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TMatlabVectorIO) write", v->typestr() );

    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();
    
    const TScalarVector * x = cptrcast( v, TScalarVector );
    int16_t               version  = 0x0100;
    int16_t               endian   = 0x4d49; // = MI
    char                  buf[116] = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char                  subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    if ( x->is_complex() )
        write_scalar( out, x->blas_cvec(), vecname );
    else
        write_scalar( out, x->blas_rvec(), vecname );
}

void
TMatlabVectorIO::write ( const BLAS::Vector< float > &  v,
                         const std::string &            filename,
                         const std::string &            vecname ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();
    
    int16_t  version  = 0x0100;
    int16_t  endian   = 0x4d49; // = MI
    char     buf[116] = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    write_scalar( out, v, vecname );
}

void
TMatlabVectorIO::write ( const BLAS::Vector< double > &  v,
                         const std::string &             filename,
                         const std::string &             vecname ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();
    
    int16_t  version  = 0x0100;
    int16_t  endian   = 0x4d49; // = MI
    char     buf[116] = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    write_scalar( out, v, vecname );
}

void
TMatlabVectorIO::write ( const BLAS::Vector< Complex< float > > &  v,
                         const std::string &                       filename,
                         const std::string &                       vecname ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();
    
    int16_t  version  = 0x0100;
    int16_t  endian   = 0x4d49; // = MI
    char     buf[116] = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    write_scalar( out, v, vecname );
}

void
TMatlabVectorIO::write ( const BLAS::Vector< Complex< double > > &  v,
                         const std::string &                        filename,
                         const std::string &                        vecname ) const
{
    unique_ptr< std::ostream >  out_ptr( open_write( filename ) );
    std::ostream &              out = * out_ptr.get();
    
    int16_t  version  = 0x0100;
    int16_t  endian   = 0x4d49; // = MI
    char     buf[116] = "MATLAB 5.0 MAT-file, Created by HLibPro";
    char     subsys[8];

    for ( size_t i = strlen(buf); i < 116; i++ ) buf[i]    = ' ';
    for ( size_t i = 0;           i < 8;   i++ ) subsys[i] = '\0';
    
    out.write( buf, 116 );
    out.write( subsys, 8 );
    out.write( reinterpret_cast< const char * >( & version ), sizeof(version) );
    out.write( reinterpret_cast< const char * >( & endian ),  sizeof(endian)  );

    write_scalar( out, v, vecname );
}

//
// read first vector from file <filename>
//
unique_ptr< TVector >
TMatlabVectorIO::read ( const std::string &  filename ) const
{
    return read( filename, "" );
}

//
// read vector named <vecname> from file (vecname = "" means first vector avail.)
//
unique_ptr< TVector >
TMatlabVectorIO::read ( const std::string &  filename,
                        const std::string &  vecname ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TMatlabVectorIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    int16_t    version, endian;
    bool       byteswap = false;
    TVector *  v = nullptr;

    in.seekg( 124 );
    in.read( reinterpret_cast< char * >( & version ), sizeof(version) );
    in.read( reinterpret_cast< char * >( & endian ),  sizeof(endian)  );
    
    if      ( endian == 0x4d49 )  byteswap = false; // MI
    else if ( endian == 0x494d )  byteswap = true;  // IM
    else
        HERROR( ERR_FMT_MATLAB, "(TMatlabVectorIO) read", "endianess indicator wrong" );

    BYTE_SWAP( version );
        
    while (( v == nullptr ) && in.good() )
        mat5_read_element( in, byteswap, vecname, false, "", nullptr, & v );
    
    return unique_ptr< TVector >( v );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
// 
// class TMatlabCoordIO
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void
TMatlabCoordIO::write ( const TCoordinate *  coord,
                        const std::string &  filename ) const
{
    write( coord, filename, "coord" );
}

void
TMatlabCoordIO::write ( const TCoordinate *  coord,
                        const std::string &  filename,
                        const std::string &  coordname ) const
{
    if ( coord == nullptr )
        HERROR( ERR_ARG, "(TMatlabCoordIO) write", "NULL coordinates" );
    
    //
    // convert coordinates to dense matrix and write that
    //

    const auto  ncoord = idx_t(coord->ncoord());
    const auto  ndim   = idx_t(coord->dim());
    
    BLAS::Matrix< double >  D( ncoord, ndim );

    for ( idx_t  i = 0; i < ncoord; ++i )
    {
        auto  p = coord->coord( i );
        
        for ( idx_t  j = 0; j < ndim; ++j )
        {
            D( i, j ) = p[j];
        }// for
    }// for

    TMatlabMatrixIO  mio;

    mio.write( D, filename, coordname );
}

//
// read first vertices in file <filename>
//
unique_ptr< TCoordinate >
TMatlabCoordIO::read ( const std::string &  filename ) const
{
    return read( filename, "" );
}

//
// read vertices named <cooname> from file <filename>
//
unique_ptr< TCoordinate >
TMatlabCoordIO::read ( const std::string &  filename,
                       const std::string &  cooname ) const
{
    if ( ! fs::exists( filename ) )
        HERROR( ERR_FNEXISTS, "(TMatlabCoordIO) read", filename );
    
    unique_ptr< std::istream >  in_ptr( open_read( filename ) );
    std::istream &              in = * in_ptr.get();

    int16_t                version, endian;
    bool                   byteswap = false;
    unique_ptr< TMatrix >  M;

    in.seekg( 124 );
    in.read( reinterpret_cast< char * >( & version ), sizeof(version) );
    in.read( reinterpret_cast< char * >( & endian ),  sizeof(endian)  );
    
    if      ( endian == 0x4d49 )  byteswap = false; // MI
    else if ( endian == 0x494d )  byteswap = true;  // IM
    else
        HERROR( ERR_FMT_MATLAB, "(TMatlabCoordIO) read",
                "endianess indicator wrong" );

    BYTE_SWAP( version );
        
    while (( M.get() == nullptr ) && in.good() )
    {
        TMatrix *  MT = nullptr;
        
        mat5_read_element( in, byteswap, cooname, false, "", & MT, nullptr );

        M = unique_ptr< TMatrix >( MT );
        
        if ( ! IS_TYPE( M, TDenseMatrix ) || M->is_complex() )
        {
            M.reset( nullptr );
        }// if
    }// while

    if ( M.get() != nullptr )
    {
        //
        // convert dense matrix into coordinates
        //

        TDenseMatrix *           D      = ptrcast( M.get(), TDenseMatrix );
        const idx_t              ncoord = idx_t( D->rows() );
        const idx_t              dim    = idx_t( D->cols() );
        std::vector< double * >  vertices( ncoord );

        for ( idx_t  i = 0; i < ncoord; i++ )
        {
            vertices[i] = new double[ dim ];

            for ( idx_t  j = 0; j < dim; j++ )
                vertices[i][j] = D->entry( i, j );
        }// for

        return make_unique< TCoordinate >( vertices, uint( dim ) );
    }// if

    return nullptr;
}

}// namespace HLIB
