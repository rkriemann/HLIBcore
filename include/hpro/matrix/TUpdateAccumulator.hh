#ifndef __HLIB_TUPDATEACCUMULATOR_HH
#define __HLIB_TUPDATEACCUMULATOR_HH
//
// Project     : HLib
// File        : TUpdateAccumulator.hh
// Description : class for handling updates to matrix blocks
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <list>
#include <deque>
#include <memory>

#include "hpro/parallel/TMutex.hh"

namespace HLIB
{

// forward decl.
class TMatrix;
class TBlockMatrix;

//
// base class for direct (e.g. non-recursive) updates
//
class TDirectMatrixUpdate
{
public:
    // dtor
    virtual ~TDirectMatrixUpdate () {}
    
    // apply internal update to matrix \a M with accuracy \a acc
    virtual auto         compute          ( const TTruncAcc &  acc ) -> std::unique_ptr< TMatrix > = 0;

    // return linear operator representing update
    virtual auto         linear_op        () const -> std::unique_ptr< TLinearOperator > = 0;
    
    // return true if result of update is dense
    virtual bool         is_dense_result  () const = 0;
    
    // return descriptive info
    virtual std::string  to_string        () const { return "TDirectMatrixUpdate"; };
};

//
// base class for recursive updates
//
class TRecursiveMatrixUpdate
{
public:
    // dtor
    virtual ~TRecursiveMatrixUpdate () {}
    
    // apply internal update to matrix \a M with accuracy \a acc
    virtual void         apply      ( TMatrix *          M,
                                      const TTruncAcc &  acc ) = 0;

    // return descriptive info
    virtual std::string  to_string  () const { return "TRecursiveMatrixUpdate"; };
};

//!
//! \ingroup Matrix_Module
//! \class   TUpdateAccumulator
//! \brief   Handles updates for a single matrix block by accumulating
//!          direct updates and recursive (pending) updates
//!
class TUpdateAccumulator : public TLockable
{
    //! @cond
    
public:

    // set of direct updates
    using  direct_updates_t    = std::deque< TDirectMatrixUpdate * >;
    
    // set of recursive updates
    using  recursive_updates_t = std::list< TRecursiveMatrixUpdate * >;
    
private:
    // accumulated updates
    std::unique_ptr< TMatrix >  _accumulated;

    // list of pending direct updates
    direct_updates_t            _pending_direct;

    // list of pending recursive updates
    recursive_updates_t         _pending_recursive;

    //! @endcond
    
public:
    //! dtor
    ~TUpdateAccumulator ()
    {
        clear_updates();
    }
    
    //! initialise matrix for accumulated updates
    void       init                   ( const TMatrix *    M );

    //! compute and apply local direct updates;
    //! - use \a dest as hint for destination type, e.g. to choose format of accumulator
    //! - directly update \a dest if \a update_dest == true (accumulator is zero afterwards)
    void       apply_direct           ( const TTruncAcc &  acc,
                                        TMatrix *          dest        = nullptr,
                                        const bool         update_dest = false );

    //! return true if accumulator holds any updates
    bool       has_updates            () const;
    
    //
    // access local updates
    //
    
    //! access accumulated updates
    auto       accumulated_updates    ()       ->       TMatrix * { return _accumulated.get(); }
    auto       accumulated_updates    () const -> const TMatrix * { return _accumulated.get(); }

    //! access set of direct pending updates
    auto       pending_direct         ()       ->       direct_updates_t &    { return _pending_direct; }
    auto       pending_direct         () const -> const direct_updates_t &    { return _pending_direct; }
    
    //! access set of recursive pending updates
    auto       pending_recursive      ()       ->       recursive_updates_t & { return _pending_recursive; }
    auto       pending_recursive      () const -> const recursive_updates_t & { return _pending_recursive; }
    
    //
    // add updates
    //
    
    //! add update matrix
    void       add_update             ( const TMatrix *    M,
                                        const TTruncAcc &  acc,
                                        const TMatrix *    dest = nullptr );
    
    //! add update from parent matrix
    void       add_parent_update      ( const TMatrix *    M,
                                        const TTruncAcc &  acc );
    
    //! add update U to set of recursive pending updates
    void       add_pending_direct     ( TDirectMatrixUpdate *  U )
    {
        if ( U == nullptr )
            return;
        
        TScopedLock  mlock( *this );

        _pending_direct.push_back( U );
    }

    //! add update U to set of recursive pending updates
    void       add_pending_recursive  ( TRecursiveMatrixUpdate *  U )
    {
        if ( U == nullptr )
            return;
        
        TScopedLock  mlock( *this );

        _pending_recursive.push_back( U );
    }
    
    //
    // clear stored updates
    //
    
    //! remove matrix with accumulated updates
    void       clear_accumulated      ()
    {
        TScopedLock  mlock( *this );

        _accumulated.reset( nullptr );
    }
    
    //! remove list of pending updates
    void       clear_pending          ();
    
    //! clear all updates
    void       clear_updates          ()
    {
        clear_accumulated();
        clear_pending();
    }

};

}// namespace HLIB

#endif  // __HLIB_TUPDATEACCUMULATOR_HH
