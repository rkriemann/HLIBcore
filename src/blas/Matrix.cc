//
// Project     : HLib
// File        : Matrix.cc
// Description : implements dense matrix class for BLAS operations
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>

#include "hpro/base/config.hh"
#include "hpro/base/System.hh"

#include "TRNG.hh"

#include "hpro/io/TMatrixIO.hh"
#include "hpro/io/TVectorIO.hh"

#include "hpro/blas/Vector.hh"

#include "hpro/blas/Matrix.hh"

namespace HLIB
{

namespace BLAS
{

// special Range type for "whole" indexset
Range  Range::all( -1, -1 );

//
// test data for invalid values, e.g. INF and NAN
//
template < typename value_t >
void
Matrix< value_t >::check_data  () const
{
    for ( idx_t  j = 0; j < idx_t( ncols() ); ++j )
    {
        for ( idx_t  i = 0; i < idx_t( nrows() ); ++i )
        {
            const value_t  val = (*this)(i,j);
            
            if ( Math::is_inf( val ) )
                HERROR( ERR_INF, "(Matrix) check_data", "" );
            
            if ( Math::is_nan( val ) )
                HERROR( ERR_NAN, "(Matrix) check_data", "" );
        }// for
    }// for
}

//
// stream output
//
template < typename T >
std::ostream &
operator << ( std::ostream &           os,
              const MatrixBase< T > &  M )
{
    if ( CFG::IO::use_matlab_syntax )
        os << '[';
    
    for ( idx_t  i = 0; i < idx_t( M.nrows() ); ++i )
    {
        for ( idx_t  j = 0; j < idx_t( M.ncols() )-1; ++j )
        {
            os << M(i,j) << ", ";
        }// for

        if ( M.ncols() > 0 )
        {
            if ( CFG::IO::use_matlab_syntax )
                os << M(i,idx_t(M.ncols())-1) << ';';
            os << std::endl;
        }// if
    }// for
        
    if ( CFG::IO::use_matlab_syntax )
        os << ']';
    
    return os;
}

template < typename T >
std::ostream &
operator << ( std::ostream &           os,
              const VectorBase< T > &  v )
{
    if ( CFG::IO::use_matlab_syntax )
        os << '[';
    
    for ( idx_t  i = 0; i < idx_t(v.length())-1; i++ )
        os << v(i) << ", ";
    
    if ( v.length() > 0 )
        os << v( idx_t(v.length())-1 );
    
    if ( CFG::IO::use_matlab_syntax )
        os << "]'";
    
    return os;
}

//
// create identity matrix
//
template < typename T >
Matrix< T >
identity ( const size_t  n )
{
    Matrix< T >  M( n, n );

    for ( idx_t  i = 0; i < idx_t( n ); i++ )
        M( i, i ) = T(1);

    return M;
}

//
// create random matrix
//
template < typename T >
Matrix< T >
random ( const size_t  nrows,
         const size_t  ncols )
{
    TRNG         rng;
    Matrix< T >  M( nrows, ncols );

    for ( idx_t  j = 0; j < idx_t( ncols ); j++ )
        for ( idx_t  i = 0; i < idx_t( nrows ); i++ )
            M( i, j ) = rng.rand< T >( typename real_type< T >::type_t( 1 ) );

    return M;
}

//
// create random vector
//
template < typename T >
Vector< T >
random ( const size_t  length )
{
    TRNG         rng;
    Vector< T >  v( length );

    for ( idx_t  i = 0; i < idx_t( length ); ++i )
        v( i ) = rng.rand< T >( typename real_type< T >::type_t( 1 ) );

    return v;
}

//
// explicit template instantiation
//

template class Matrix< float >;
template class Matrix< double >;
template class Matrix< Complex< float > >;
template class Matrix< Complex< double > >;


template std::ostream &
operator << < Matrix< float > >             ( std::ostream &,
                                              const MatrixBase< Matrix< float > > & );
template std::ostream &
operator << < Matrix< double > >            ( std::ostream &,
                                              const MatrixBase< Matrix< double > > & );

template std::ostream &
operator << < Matrix< Complex< float > > >  ( std::ostream &,
                                              const MatrixBase< Matrix< Complex< float > > > & );
template std::ostream &
operator << < Matrix< Complex< double > > > ( std::ostream &,
                                              const MatrixBase< Matrix< Complex< double > > > & );


template std::ostream &
operator << < Vector< float > >             ( std::ostream &,
                                              const VectorBase< Vector< float > > & );
template std::ostream &
operator << < Vector< double > >            ( std::ostream &,
                                              const VectorBase< Vector< double > > & );

template std::ostream &
operator << < Vector< Complex< float > > >  ( std::ostream &,
                                              const VectorBase< Vector< Complex< float > > > & );
template std::ostream &
operator << < Vector< Complex< double > > > ( std::ostream &,
                                              const VectorBase< Vector< Complex< double > > > & );


template Matrix< float >
identity< float >             ( const size_t  n );
    
template Matrix< double >
identity< double >            ( const size_t  n );

template Matrix< Complex< float > >
identity< Complex< float > >  ( const size_t  n );
    
template Matrix< Complex< double > >
identity< Complex< double > > ( const size_t  n );


template Matrix< float >
random< float >             ( const size_t  nrows,
                              const size_t  ncols );
    
template Matrix< double >
random< double >            ( const size_t  nrows,
                              const size_t  ncols );
    
template Matrix< Complex< float > >
random< Complex< float > >  ( const size_t  nrows,
                              const size_t  ncols );
    
template Matrix< Complex< double > >
random< Complex< double > > ( const size_t  nrows,
                              const size_t  ncols );
    
template Vector< float >
random< float >             ( const size_t  length );
    
template Vector< double >
random< double >            ( const size_t  length );
    
template Vector< Complex< float > >
random< Complex< float > >  ( const size_t  length );
    
template Vector< Complex< double > >
random< Complex< double > > ( const size_t  length );
    
}// namespace BLAS

//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//!
//! write matrix to file
//!
template < typename T >
void
write ( const BLAS::Matrix< T > & M,
        const std::string &       filename,
        const std::string &       matname )
{
    TMatlabMatrixIO  mio;
    
    mio.write( M, filename, matname );
}

// instantiate
template void write< float >   ( const BLAS::Matrix< float > &,  const std::string &, const std::string & );
template void write< double >  ( const BLAS::Matrix< double > &, const std::string &, const std::string & );

template void write< Complex< float > >   ( const BLAS::Matrix< Complex< float > > &,  const std::string &, const std::string & );
template void write< Complex< double > >  ( const BLAS::Matrix< Complex< double > > &, const std::string &, const std::string & );

//!
//! write vector to file
//!
template < typename T >
void
write ( const BLAS::Vector< T > & v,
        const std::string &       filename,
        const std::string &       matname )
{
    TMatlabVectorIO  vio;
    
    vio.write( v, filename, matname );
}

// instantiate
template void write< float >    ( const BLAS::Vector< float > &,   const std::string &, const std::string & );
template void write< double >   ( const BLAS::Vector< double > &,  const std::string &, const std::string & );

template void write< Complex< float > >   ( const BLAS::Vector< Complex< float > > &,  const std::string &, const std::string & );
template void write< Complex< double > >  ( const BLAS::Vector< Complex< double > > &, const std::string &, const std::string & );

}// namespace DBG

}// namespace HLIB
