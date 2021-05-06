#ifndef __HLIB_TVECTORIO_HH
#define __HLIB_TVECTORIO_HH
//
// Project     : HLib
// File        : TVectorIO.hh
// Description : class for vector input/output
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/blas/Vector.hh"
#include "hpro/vector/TVector.hh"

namespace HLIB
{

//!
//! \{
//! \name Vector I/O
//! Functions for input/output of vectors.
//!

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TVectorIO
//! \brief    Base class for vector IO defining interface.
//!
class TVectorIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TVectorIO ();

    virtual ~TVectorIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a v to file \a fname
    virtual void
    write ( const TVector *      v,
            const std::string &  fname ) const;

    //! read and return vector from file \a fname
    virtual std::unique_ptr< TVector >
    read  ( const std::string &  fname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TAutoVectorIO
//! \brief    Class for vector I/O with automatic
//!           file format detection
//!
class TAutoVectorIO : public TVectorIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TAutoVectorIO () {}

    virtual ~TAutoVectorIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a A to file \a filename
    virtual void
    write ( const TVector *      A,
            const std::string &  filename ) const;

    //! write vector \a A to file \a filename with optional
    //! vector name \a vecname (if file format has corresponding support)
    virtual void
    write ( const TVector *      A,
            const std::string &  filename,
            const std::string &  vecname ) const;

    //! read and return vector from file \a filename
    virtual std::unique_ptr< TVector >
    read  ( const std::string &  filename ) const;

    //! read and return vector from file \a filename with optional
    //! vector name \a vecname (if file format has corresponding support)
    virtual std::unique_ptr< TVector >
    read  ( const std::string &  filename,
            const std::string &  vecname ) const;
};

//!
//! \ingroup  IO_Module
//! \brief    Read vector from file with automatic file format detection.
//!
//!           Read vector from file \a filename while the format of the file
//!           is detected automatically. If the file format supports storage
//!           of several matrices (or other objects), only the first vector
//!           will be read.
//!
std::unique_ptr< TVector >
read_vector  ( const char *  filename );

//!
//! \ingroup  IO_Module
//! \brief    Read vector from file with automatic file format detection.
//!
//!           Read vector from file \a filename while the format of the file
//!           is detected automatically. If the file format supports named
//!           matrices, the vector with the name \a vecname will be read.
//!           Furthermore, if the file format supports storage of several matrices
//!           (or other objects), only the first vector will be read.
//!
std::unique_ptr< TVector >
read_vector  ( const char *  filename,
               const char *  vecname );

//!
//! \ingroup  IO_Module
//! \brief    Write vector to file with automatic choice of file format.
//!
//!           Write vector \a A to file \a filename while the format of the file
//!           is chosen automatically based upon the filename suffix.
//!           - hm                    : HLIBpro format
//!           - mat, m                : Matlab format
//!           - hb, rb, rua, rsa, psa : Harwell/Boeing format
//!           - rhs, sol              : SAMG format
//!           - mtx                   : VectorMarket format
//!
void
write_vector  ( const TVector *  A,
                const char *     filename );

//!
//! \ingroup  IO_Module
//! \brief    Write vector to file with automatic choice of file format.
//!
//!           Write vector \a A to file \a filename while the format of the file
//!           is chosen automatically based upon the filename suffix. If the file
//!           format supports named matrices, the vector will be stored under the
//!           name \a vecname.
//!
void
write_vector  ( const TVector *  A,
                const char *     filename,
                const char *     vecname );

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TSAMGVectorIO
//! \brief    Class for vector I/O in SAMG format
//!
class TSAMGVectorIO : public TVectorIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TSAMGVectorIO () {}

    virtual ~TSAMGVectorIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a v to file \a fname
    virtual void
    write ( const TVector *      v,
            const std::string &  fname ) const;

    //! read and return vector from file \a fname
    virtual std::unique_ptr< TVector >
    read  ( const std::string &  fname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMatlabVectorIO
//! \brief    Class for vector I/O in Matlab format
//!
class TMatlabVectorIO : public TVectorIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TMatlabVectorIO () {}

    virtual ~TMatlabVectorIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a v with name "v" to file \a fname
    virtual void
    write ( const TVector *      v,
            const std::string  & fname ) const;
    
    //! write vector \a v with name \a vname to file \a fname
    virtual void
    write ( const TVector *      v,
            const std::string  & fname,
            const std::string  & vname ) const;

    //! read first vector in file \a fname
    virtual std::unique_ptr< TVector >
    read  ( const std::string & fname ) const;

    //! read and return vector named \a vname from file \a fname
    //! (if vname = "", return first vector)
    virtual std::unique_ptr< TVector >
    read  ( const std::string & fname,
            const std::string & vname ) const;

    // write BLAS vector \a v with name \a mname to file \a fname
    virtual void
    write ( const BLAS::Vector< float > &   v,
            const std::string &             fname,
            const std::string &             mname ) const;
    virtual void
    write ( const BLAS::Vector< double > &  v,
            const std::string &             fname,
            const std::string &             mname ) const;
    virtual void
    write ( const BLAS::Vector< Complex< float > > &   v,
            const std::string &                        fname,
            const std::string &                        mname ) const;
    virtual void
    write ( const BLAS::Vector< Complex< double > > &  v,
            const std::string &                        fname,
            const std::string &                        mname ) const;

};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THLibVectorIO
//! \brief    Class for vector I/O in HLIB format
//!
class THLibVectorIO : public TVectorIO
{
protected:
    // use compression
    const bool  _compressed;
    
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    THLibVectorIO ( const bool compressed = false )
            : _compressed(compressed)
    {}

    virtual ~THLibVectorIO () {}
    
    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a v to file \a fname
    virtual void
    write ( const TVector *      v,
            const std::string &  fname ) const;

    //! read and return vector from file \a fname
    virtual std::unique_ptr< TVector >
    read  ( const std::string &  fname ) const;
};

///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    THBVectorIO
//! \brief    Class for vector I/O in Harwell-Boeing and
//!           Rutherford-Boeing format
//!
class THBVectorIO : public TVectorIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    THBVectorIO () {}

    virtual ~THBVectorIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a v to file \a fname
    virtual void
    write ( const TVector *      v,
            const std::string &  fname ) const;

    //! read and return vector from file \a fname
    virtual std::unique_ptr< TVector >
    read  ( const std::string &  fname ) const;
};


///////////////////////////////////////////////////
//!
//! \ingroup  IO_Module
//! \class    TMMVectorIO
//! \brief    Class for vector I/O in MatrixMarket format
//!
class TMMVectorIO : public TVectorIO
{
public:
    ///////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TMMVectorIO () {}

    virtual ~TMMVectorIO () {}

    ///////////////////////////////////////////////
    //
    // interface for I/O
    //

    //! write vector \a v to file \a fname
    virtual void
    write ( const TVector *      v,
            const std::string &  fname ) const;

    //! read and return vector from file \a fname
    virtual std::unique_ptr< TVector >
    read  ( const std::string &  fname ) const;
};

//! \}

}// namespace

#endif  // __HLIB_TVECTORIO_HH
