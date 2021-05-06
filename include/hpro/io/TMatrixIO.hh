#ifndef __HLIB_TMATRIXIO_HH
#define __HLIB_TMATRIXIO_HH
//
// Project     : HLib
// File        : TMatrixIO.hh
// Description : class for matrix input/output
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/blas/Matrix.hh"
#include "hpro/matrix/TMatrix.hh"

namespace HLIB
{

//!
//! \{
//! \name Matrix I/O
//! Functions for input/output of matrices.
//!

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMatrixIO
//! \brief    Base class for Matrix IO defining interface
//!
class TMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TMatrixIO ();

    virtual ~TMatrixIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a name
    virtual void  write ( const TMatrix *      A,
                          const std::string &  name ) const;
    
    //! read and return matrix from file \a name
    virtual auto  read  ( const std::string &  name ) const -> std::unique_ptr< TMatrix >;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TAutoMatrixIO
//! \brief    Class for matrix I/O with automatic
//!           file format detection
//!
class TAutoMatrixIO : public TMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TAutoMatrixIO () {}

    virtual ~TAutoMatrixIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a filename
    virtual void
    write ( const TMatrix *      A,
            const std::string &  filename ) const;

    //! write matrix \a A to file \a filename with optional
    //! matrix name \a matname (if file format has corresponding support)
    virtual void
    write ( const TMatrix *      A,
            const std::string &  filename,
            const std::string &  matname ) const;

    //! read and return matrix from file \a filename
    virtual std::unique_ptr< TMatrix >
    read  ( const std::string &  name ) const;

    //! read and return matrix from file \a filename with optional
    //! matrix name \a matname (if file format has corresponding support)
    virtual std::unique_ptr< TMatrix >
    read  ( const std::string &  filename,
            const std::string &  matname ) const;
};

//!
//! \ingroup  IO_Module
//! \brief    Read matrix from file with automatic file format detection.
//!
//!           Read matrix from file \a filename while the format of the file
//!           is detected automatically. If the file format supports storage
//!           of several matrices (or other objects), only the first matrix
//!           will be read.
//!
std::unique_ptr< TMatrix >
read_matrix  ( const std::string &  filename );

//!
//! \ingroup  IO_Module
//! \brief    Read matrix from file with automatic file format detection.
//!
//!           Read matrix from file \a filename while the format of the file
//!           is detected automatically. If the file format supports named
//!           matrices, the matrix with the name \a matname will be read.
//!           Furthermore, if the file format supports storage of several matrices
//!           (or other objects), only the first matrix will be read.
//!
std::unique_ptr< TMatrix >
read_matrix  ( const std::string &  filename,
               const std::string &  matname );

//!
//! \ingroup  IO_Module
//! \brief    Write matrix to file with automatic choice of file format.
//!
//!           Write matrix \a A to file \a filename while the format of the file
//!           is chosen automatically based upon the filename suffix.
//!           - hm                    : HLIBpro format
//!           - mat, m                : Matlab format
//!           - hb, rb, rua, rsa, psa : Harwell/Boeing format
//!           - amg                   : SAMG format
//!           - mtx                   : MatrixMarket format
//!           - hdf, h5               : HDF5 format
//!
void
write_matrix  ( const TMatrix *      A,
                const std::string &  filename );

//!
//! \ingroup  IO_Module
//! \brief    Write matrix to file with automatic choice of file format.
//!
//!           Write matrix \a A to file \a filename while the format of the file
//!           is chosen automatically based upon the filename suffix. If the file
//!           format supports named matrices, the matrix will be stored under the
//!           name \a matname.
//!
void
write_matrix  ( const TMatrix *      A,
                const std::string &  filename,
                const std::string &  matname );

//!
//! \ingroup  IO_Module
//! \brief    Read linear operator from file
//!
//!           Read linear operator from file \a filename. Currently only the
//!           HLIBpro file format supports linear operators.
//!
std::unique_ptr< TLinearOperator >
read_linop  ( const std::string &  filename );

//!
//! \ingroup  IO_Module
//! \brief    Write linear operator to file.
//!
//!           Write linear operator \a A to file \a filename in HLIBpro format.
//!
void
write_linop ( const TLinearOperator *  A,
              const std::string &      filename );

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TOctaveMatrixIO
//! \brief    Class for matrix I/O in octave format
//!
class TOctaveMatrixIO : public TMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TOctaveMatrixIO ();

    virtual ~TOctaveMatrixIO ();
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a name
    virtual void
    write ( const TMatrix *      A,
            const std::string &  name ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TSAMGMatrixIO
//! \brief    Class for matrix I/O in SAMG format
//!
class TSAMGMatrixIO : public TMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TSAMGMatrixIO () {}

    virtual ~TSAMGMatrixIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a name
    virtual void
    write ( const TMatrix *      A,
            const std::string &  name ) const;

    //! read and return matrix from file \a name
    virtual std::unique_ptr< TMatrix >
    read  ( const std::string &  name ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMatlabMatrixIO
//! \brief    Class for matrix I/O in Matlab format
//!
class TMatlabMatrixIO : public TMatrixIO
{
private:
    //! if true (default), any permutation of the matrices will be
    //! applied before writing them to the file
    const bool  _permute;
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct Matlab matrix IO object with internal permutation set to \a perm
    TMatlabMatrixIO ( const bool perm = CFG::IO::permute_save )
            : _permute( perm )
    {}

    virtual ~TMatlabMatrixIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A with name "M" to file \a fname
    virtual void
    write ( const TMatrix *      A,
            const std::string &  fname ) const;

    // write matrix \a A with name \a mname to file \a fname
    virtual void
    write ( const TMatrix *      A,
            const std::string &  fname,
            const std::string &  mname ) const;

    // write linear operator \a A with name \a mname to file \a fname
    virtual void
    write ( const TLinearOperator *  A,
            const std::string &      fname,
            const std::string &      mname ) const;

    // write BLAS matrix \a A with name \a mname to file \a fname
    virtual void
    write ( const BLAS::Matrix< float > &    A,
            const std::string &              fname,
            const std::string &              mname ) const;
    virtual void
    write ( const BLAS::Matrix< double > &   A,
            const std::string &              fname,
            const std::string &              mname ) const;
    virtual void
    write ( const BLAS::Matrix< Complex< float > > &  A,
            const std::string &                       fname,
            const std::string &                       mname ) const;
    virtual void
    write ( const BLAS::Matrix< Complex< double > > &  A,
            const std::string &                        fname,
            const std::string &                        mname ) const;

    //! read and return matrix from file \a fname (assuming only one entry available)
    virtual std::unique_ptr< TMatrix >
    read  ( const std::string &  fname ) const;

    //! read matrix with name \a mname from file \a fname
    virtual std::unique_ptr< TMatrix >
    read  ( const std::string &  fname,
            const std::string &  mname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THLibMatrixIO
//! \brief    Class for matrix I/O in HLIB format
//!
class THLibMatrixIO : public TMatrixIO
{
protected:
    //! flag for using compression
    bool  _compressed;
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct HLIB matrix IO object with compression set to \a compressed
    THLibMatrixIO ( const bool compressed = false )
            : _compressed(compressed)
    {}

    virtual ~THLibMatrixIO () {}

    //! set internal compression usage to \a compressed
    void set_compression ( const bool  compressed ) { _compressed = compressed; }
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a fname
    virtual void
    write ( const TMatrix *      A,
            const std::string &  fname ) const;

    //! read and return matrix from file \a fname
    virtual std::unique_ptr< TMatrix >
    read  ( const std::string &  fname ) const;

    ///////////////////////////////////////////////
    //
    // interface for linear operators
    //

    //! write linear operator \a A to file \a fname
    virtual void
    write_linop ( const TLinearOperator * A,
                  const std::string &     fname ) const;

    //! read and return linear operator from file \a fname
    virtual std::unique_ptr< TLinearOperator >
    read_linop  ( const std::string &     fname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TPLTMGMatrixIO
//! \brief    Class for matrix I/O in PLTMG format
//!
class TPLTMGMatrixIO : public TMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TPLTMGMatrixIO () {}

    virtual ~TPLTMGMatrixIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a fname
    virtual void
    write ( const TMatrix *      A,
            const std::string &  fname ) const;

    //! read and return matrix from file \a fname
    virtual std::unique_ptr< TMatrix >
    read  ( const std::string &  fname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THBMatrixIO
//! \brief    Class for matrix I/O in Harwell-Boeing and
//!           Rutherford-Boeing format
//!
class THBMatrixIO : public TMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    THBMatrixIO () {}

    virtual ~THBMatrixIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a fname
    virtual void
    write ( const TMatrix *      A,
            const std::string &  fname ) const;

    //! read and return matrix from file \a fname
    virtual std::unique_ptr< TMatrix >
    read  ( const std::string &  fname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMMMatrixIO
//! \brief    Class for matrix I/O in MatrixMarket format.
//!
class TMMMatrixIO : public TMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TMMMatrixIO () {}

    virtual ~TMMMatrixIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a fname
    virtual void
    write ( const TMatrix *      A,
            const std::string &  fname ) const;

    //! read and return matrix from file \a fname
    virtual std::unique_ptr< TMatrix >
    read  ( const std::string &  fname ) const;
};


///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THDF5MatrixIO
//! \brief    Class for matrix I/O in HDF5 format.
//!
class THDF5MatrixIO : public TMatrixIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    THDF5MatrixIO () {}

    virtual ~THDF5MatrixIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write matrix \a A to file \a fname
    virtual void
    write ( const TMatrix *      A,
            const std::string &  fname ) const;

    // write matrix \a A with name \a mname to file \a fname
    virtual void
    write ( const TMatrix *      A,
            const std::string &  fname,
            const std::string &  mname ) const;

    // write matrix \a A with name \a mname to file \a fname
    void
    write ( const BLAS::Matrix< real > &     A,
            const std::string &              fname,
            const std::string &              mname ) const;
    void
    write ( const BLAS::Matrix< complex > &  A,
            const std::string &              fname,
            const std::string &              mname ) const;
    
    //! read and return matrix from file \a fname
    virtual std::unique_ptr< TMatrix >
    read  ( const std::string &  fname ) const;
};

//! \}

}// namespace HLIB

#endif  // __HLIB_TMATRIXIO_HH
