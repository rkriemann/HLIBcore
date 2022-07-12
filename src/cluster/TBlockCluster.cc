//
// Project     : HLIBpro
// File        : TBlockCluster.cc
// Description : class for a block cluster tree
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2022. All Rights Reserved.
//

#include <iostream>
#include <mutex>

#include "unordered_map.hh"

#include "hpro/base/System.hh"
#include "hpro/parallel/TMutex.hh"

#include "treealg.hh"

#include "hpro/cluster/TBlockCluster.hh"

namespace Hpro
{

using std::list;
using std::vector;

namespace
{

//
// return globally unique ID
//
int
get_id ()
{
    static std::mutex  id_mutex;
    static int         id_counter = 0;
    int                ret_val    = 0;

    {
        std::lock_guard< std::mutex >  lock( id_mutex );

        ret_val = id_counter++;
    }

    return ret_val;
}

}// namespace anonymous

////////////////////////////////////////////////////////
//
// constructor and destructor
//

TBlockCluster::TBlockCluster ( TBlockCluster *  aparent )
        : _id( get_id() )
        , _parent( aparent )
        , _rowcl( nullptr )
        , _colcl( nullptr )
        , _nrows( 0 )
        , _ncols( 0 )
        , _sons( 0 )
        , _adm( false )
        , _procs( PROCSET_INVALID )
{}

TBlockCluster::TBlockCluster ( TBlockCluster *  aparent,
                               TCluster *       arow_ct,
                               TCluster *       acol_ct )
        : _id( get_id() )
        , _parent( aparent )
        , _rowcl( arow_ct )
        , _colcl( acol_ct )
        , _nrows( arow_ct->nsons() )
        , _ncols( acol_ct->nsons() )
        , _sons( _nrows * _ncols )
        , _adm( false )
        , _procs( PROCSET_INVALID )
{}

TBlockCluster::~TBlockCluster ()
{
    for ( uint i = 0; i < nsons(); i++ )
        delete _sons[i];
}

//
// set row cluster
//
void
TBlockCluster::set_rowcl ( TCluster *  cl )
{
    _rowcl = cl;

    if ( _rowcl != nullptr )
        set_layout( _rowcl->nsons(), ncols() );
    else
        set_layout( 0, ncols() );
}

//
// set column cluster
//
void
TBlockCluster::set_colcl ( TCluster *  cl )
{
    _colcl = cl;

    if ( _colcl != nullptr )
        set_layout( nrows(), _colcl->nsons() );
    else
        set_layout( nrows(), 0 );
}

//
// set row and column cluster
void
TBlockCluster::set_clusters ( TCluster *  row_cl,
                              TCluster *  col_cl )
{
    _rowcl = row_cl;
    _colcl = col_cl;

    if ( _rowcl != nullptr )
    {
        if ( _colcl != nullptr )
            set_layout( _rowcl->nsons(), _colcl->nsons() );
        else
            set_layout( _rowcl->nsons(), 0 );
    }// if
    else
    {
        if ( _colcl != nullptr )
            set_layout( 0, _colcl->nsons() );
        else
            set_layout( 0, 0 );
    }// else
}

//
// change block layout
//
void
TBlockCluster::set_layout  ( const uint  anrows,
                             const uint  ancols )
{
    if (( anrows == _nrows ) && ( ancols == _ncols ))
        return;

    _nrows = anrows;
    _ncols = ancols;
    
    _sons.resize( _nrows * _ncols );
}

//
// make node a leaf
//
void
TBlockCluster::make_leaf ()
{
    for ( auto  node : _sons )
    {
        if ( node != nullptr )
            delete node;
    }// for

    set_layout( 0, 0 );
}

////////////////////////////////////////////////
//
// tree interface
//

//
// set number of sons
//

// void
// TBlockCluster::set_nsons ( const uint n )
// {
//     uint  old = nsons();
    
//     _sons.resize( n );

//     for ( uint i = old; i < n; i++ )
//         _sons[i] = nullptr;
// }

// void
// TBlockCluster::adjust_nsons ()
// {
//     if ((_rowcl != nullptr) && (_colcl != nullptr))
//     {
//         uint  tn = _rowcl->nsons();
//         uint  sn = _colcl->nsons();

//         // if one of the clusters has a son, the other cluster
//         // automatically also has a son (itself)
//         if      (( tn != 0 ) && ( sn == 0 )) sn = 1;
//         else if (( tn == 0 ) && ( sn != 0 )) tn = 1;
        
//         set_nsons( tn * sn );
//     }// if
// }

//
// return son wrt to given clusters
//
TBlockCluster *
TBlockCluster::son_cl ( const TCluster *  row_cl,
                        const TCluster *  col_cl )
{
    for ( uint i = 0; i < nsons(); i++ )
    {
        if ( _sons[i] == nullptr )
            continue;

        if (( _sons[i]->rowcl() == row_cl ) && (_sons[i]->colcl() == col_cl ))
            return _sons[i];
    }// for
    
    return nullptr;
}

const TBlockCluster *
TBlockCluster::son_cl ( const TCluster *  row_cl,
                        const TCluster *  col_cl ) const
{
    for ( uint i = 0; i < nsons(); i++ )
    {
        if ( _sons[i] == nullptr )
            continue;

        if (( _sons[i]->rowcl() == row_cl ) && (_sons[i]->colcl() == col_cl ))
            return _sons[i];
    }// for

    return nullptr;
}

//
// access sons block-wise
//

TBlockCluster *
TBlockCluster::son ( const uint  i,
                     const uint  j )
{
    const auto  idx = ( j * nrows()) + i;

    if ( idx >= nsons() )
        HERROR( ERR_ARG, "(TBlockCluster) son", Hpro::to_string( "son_%d,%d not available", i, j ) );
    
    return _sons[ idx ];
}
    
const TBlockCluster *
TBlockCluster::son ( const uint  i,
                     const uint  j ) const
{
    const auto  idx = ( j * nrows()) + i;
    
    if ( idx >= nsons() )
        HERROR( ERR_ARG, "(TBlockCluster) son", Hpro::to_string( "son_%d,%d not available", i, j ) );
    
    return _sons[ idx ];
}

//
// set i'th son
//
void
TBlockCluster::set_son ( const uint     i,
                         TBlockCluster *  son_bct,
                         const bool       del_son )
{
    if ( i >= nsons() )
        HERROR( ERR_ARG, "(TBlockCluster) set_son", "argument index too high" );
    
    if ( del_son && (( _sons[i] != nullptr ) && ( _sons[i] != son_bct )))
        delete _sons[i];

    _sons[i] = son_bct;
}

//
// add a son (gets number i+1)
//
void
TBlockCluster::add_son ( TBlockCluster *  son_bct )
{
    if ( son_bct == nullptr )
        return;

    //
    // check sons and assign first free son
    //

    for ( uint i = 0; i < nsons(); i++ )
    {
        if ( _sons[i] == nullptr )
        {
            _sons[i] = son_bct;
            return;
        }// if
    }// for

    HERROR( ERR_CONSISTENCY, "(TBlockCluster) add_son", "no free slot" );
}

//
// set son block-wise
//
void
TBlockCluster::set_son ( const uint     i,
                         const uint     j,
                         TBlockCluster *  son_bct,
                         const bool       del_son )
{
    const auto  idx = ( j * nrows()) + i;
    
    if ( idx >= nsons() )
        HERROR( ERR_ARG, "(TBlockCluster) set_son", Hpro::to_string( "son_%d,%d not available", i, j ) );
    
    if ( del_son && (( _sons[ idx ] != nullptr ) && ( _sons[ idx ] != son_bct )))
        delete _sons[ idx ];

    _sons[ idx ] = son_bct;
}

//
// return number of nodes in tree
//
uint
TBlockCluster::nnodes () const
{
    return uint(Hpro::Tree::nnodes( this ));
}

//
// return depth of tree
//
uint
TBlockCluster::depth () const
{
    return uint(Hpro::Tree::depth( this ));
}

//
// collect leaves in tree according
//
void
TBlockCluster::collect_leaves ( list< TBlockCluster * > &  leaves,
                                const int                  tdepth,
                                const int                  level ) const
{
    if ( is_leaf() || ((tdepth != -1) && (level == tdepth)))
        leaves.push_back( const_cast< TBlockCluster * >( this ) );
    else
    {
        //
        // collect children
        //
        
        for ( uint i = 0; i < nsons(); i++ )
        {
            if ( _sons[i] != nullptr )
                _sons[i]->collect_leaves( leaves, tdepth, level + 1 );
        }// for
    }// else
}

//
// assign cluster (and sons) to processor
//
void
TBlockCluster::set_procs ( const TProcSet &  ps,
                           const bool        recursive )
{
    _procs = ps;

    if ( recursive )
    {
        for ( uint i = 0; i < nsons(); i++ )
        {
            if ( _sons[i] != nullptr )
                _sons[i]->set_procs( ps, true );
        }// for
    }// if
}

//
// assign processor set according to processor sets of sons
//
void
TBlockCluster::assign_procs ()
{
    if ( is_leaf() )
        return;
    
    //
    // first assign processor on sons and check,
    // if all sons are on the same processor
    //

    if ( procs() != PROCSET_INVALID )
    {
        //
        // broadcast processor set to sons
        //
        
        for ( uint i = 0; i < nsons(); i++ )
        {
            if (( _sons[i] != nullptr ) && ( _sons[i]->procs() == PROCSET_INVALID ))
                _sons[i]->set_procs( procs() );

            _sons[i]->assign_procs();
        }// for
    }// if
    else
    {
        //
        // use joined processor set from sons
        //
        
        TProcSet  ps = PROCSET_INVALID;

        for ( uint i = 0; i < nsons(); i++ )
        {
            if ( _sons[i] == nullptr )
                continue;
        
            _sons[i]->assign_procs();
            
            if ( _sons[i]->procs() == PROCSET_INVALID )
                HERROR( ERR_CONSISTENCY, "(TBlockCluster) assign_procs", "son with invalid processor set" );
                
            if ( ps != PROCSET_INVALID )
            {
                ps = join( ps, _sons[i]->procs() );
            }// if
            else
            {
                ps = _sons[i]->procs();
            }// else
        }// for

        set_procs( ps );
    }// else
}

//
// compute sparsity/sharing constant of tree
//

uint
TBlockCluster::compute_c_sp () const
{
    using  cluster_map_t = std::unordered_map< const TCluster *, uint >;

    //
    // assign IDs in rowcl and colcl
    //

    cluster_map_t   rowcl_no, colcl_no;
    const bool      same_ct = (rowcl() == colcl());
    
    //
    // count block clusters per cluster
    //

    list< const TBlockCluster * > nodes;
    uint                          c_sp = 0;
    
    nodes.push_back( this );
    
    while ( nodes.size() > 0 )
    {
        const TBlockCluster * bc  = behead( nodes );
        const TCluster *      row = bc->rowcl();
        const TCluster *      col = bc->colcl();
        
        rowcl_no[ row ]++;
        c_sp = std::max( c_sp, rowcl_no[ row ] );
        
        if ( ! same_ct )
        {
            colcl_no[ col ]++;
            c_sp = std::max( c_sp, colcl_no[ col ] );
        }// if
        
        if ( ! bc->is_leaf() )
        {
            for ( uint i = 0; i < bc->nsons(); i++ )
            {
                if ( bc->son(i) != nullptr )
                    nodes.push_front( bc->son(i) );
            }// for
        }// else
    }// while
    
    return c_sp;
}

uint
TBlockCluster::compute_c_sh ( const uint nprocs ) const
{
    const auto  n_rowcl = uint( rowcl()->size() );
    const auto  n_colcl = uint( colcl()->size() );

    //
    // first compute c_sh^rowcl
    //
    
    vector< int >  pinfo( n_rowcl * nprocs, 0 );
    uint           c_sh_rowcl = 1;

    list< TBlockCluster * >  leaves;

    collect_leaves( leaves );

    for ( auto  leaf : leaves )
    {
        idx_t   ofs  = leaf->rowcl()->first();
        size_t  size = leaf->rowcl()->size();

        for ( size_t i = ofs; i < ofs+size; i++ )
        {
            HERROR( ERR_NOT_IMPL, "", "" );
            if ( ! leaf->procs().empty() )
//            if ( leaf->procs().master() != NO_PROC )
                pinfo[ i * nprocs + leaf->procs().master() ] = 1;
        }// for
    }// while

    for ( uint i = 0; i < n_rowcl; i++ )
    {
        uint loc_p = 0;

        for ( uint j = 0; j < nprocs; j++ )
            loc_p += pinfo[ i * nprocs + j ];

        c_sh_rowcl = std::max( c_sh_rowcl, loc_p );
    }// for

    //
    // now compute c_sh^colcl
    //

    uint c_sh_colcl = 1;
    
    pinfo.resize( n_colcl * nprocs, 0 );

    for ( uint i = 0; i < n_colcl*nprocs; i++ )
        pinfo[i] = 0;

    for ( auto  leaf : leaves )
    {
        idx_t   ofs  = leaf->colcl()->first();
        size_t  size = leaf->colcl()->size();

        for ( size_t i = ofs; i < ofs+size; i++ )
        {
            HERROR( ERR_NOT_IMPL, "", "" );
            if ( ! leaf->procs().empty() )
//            if ( leaf->procs().master() != NO_PROC )
                pinfo[ i * nprocs + leaf->procs().master() ] = 1;
        }// for
    }// while

    for ( uint i = 0; i < n_colcl; i++ )
    {
        uint loc_p = 0;

        for ( uint j = 0; j < nprocs; j++ )
            loc_p += pinfo[ i * nprocs + j ];

        c_sh_colcl = std::max( c_sh_colcl, loc_p );
    }// for

    return std::max( c_sh_rowcl, c_sh_colcl );
}

//
// return true if given cluster is a subcluster of this
//
bool
TBlockCluster::is_sub_cluster ( const TBlockCluster * c ) const
{
    if ( c->rowcl()->first() < rowcl()->first() ) return false;
    if ( c->rowcl()->last()  > rowcl()->last()  ) return false;
    if ( c->colcl()->first() < colcl()->first() ) return false;
    if ( c->colcl()->last()  > colcl()->last()  ) return false;

    return true;
}
    
//
// stream output
//
void
TBlockCluster::print ( const uint ofs ) const
{
    for ( uint i = 0; i < ofs; i++ )
        std::cout << ' ';

    std::cout << is();

    if ( _adm ) std::cout << ", adm";
    else        std::cout << ", not adm";

    std::cout << _procs.to_string()
              << std::endl;

    for ( uint i = 0; i < nsons(); i++ )
    {
        if ( _sons[i] != nullptr )
            _sons[i]->print( ofs + 4 );
    }// if
}

//
// virtual constructor
//
TBlockCluster *
TBlockCluster::copy  () const
{
    std::unique_ptr< TBlockCluster >  bc( new TBlockCluster( const_cast< TBlockCluster * >( parent() ),
                                                             const_cast< TCluster * >( rowcl() ),
                                                             const_cast< TCluster * >( colcl() ) ) );

    bc->set_adm( is_adm() );
    bc->set_procs( procs() );
    
    bc->set_layout( nrows(), ncols() );

    for ( uint i = 0; i < nrows(); i++ )
    {
        for ( uint j = 0; j < ncols(); j++ )
        {
            if ( son(i,j) != nullptr )
                bc->set_son( i, j, son(i,j)->copy() );
        }// for
    }// for
    
    return bc.release();
}

//
// return size in bytes used by this object
//
size_t
TBlockCluster::byte_size () const
{
    size_t count = 0;

    count += ( sizeof(TBlockCluster*) +
               sizeof(TCluster*) * 2 +
               sizeof(_sons) + sizeof(TBlockCluster*) * _sons.size() +
               sizeof(bool) +
               sizeof(uint) );

    for ( uint  i = 0; i < nsons(); i ++ )
        if (son(i) != nullptr)
            count += son(i)->byte_size();

    return count;
}

//
// The first levels in the block cluster tree with no leaves are eliminated
// such that the tree root will directly have sons on this level.
//
void
flatten_leaf ( TBlockCluster *  root )
{
    if (( root == nullptr ) || ( root->is_leaf() ))
        return;

    //
    // recurse and flatten in sons and collect new son clusters
    //

    std::vector< TBlockCluster * >  bcls( 1 );
    std::vector< TBlockCluster * >  sons;
    uint                            nrows      = 1;
    uint                            ncols      = 1;
    uint                            level      = 0;
    bool                            found_leaf = false;
        
    bcls[0] = root;

    while ( ! found_leaf )
    {
        //
        // first, determine layout of next level
        //

        vector< uint >  son_nrows( nrows, 0 );
        vector< uint >  son_ncols( ncols, 0 );
        
        for ( uint  i = 0; i < nrows; ++i )
        {
            for ( uint  j = 0; j < ncols; ++j )
            {
                auto  bcl = bcls[ j * nrows + i ];

                if ( bcl == nullptr )
                    continue;
                
                if ( son_nrows[i] == 0 )
                    son_nrows[i] = bcl->nrows();
                else if ( son_nrows[i] != bcl->nrows() )
                    HERROR( ERR_CONSISTENCY, "flatten_leaf", "nodes in same row have different son layout" );
                
                if ( son_ncols[j] == 0 )
                    son_ncols[j] = bcl->ncols();
                else if ( son_ncols[j] != bcl->ncols() )
                    HERROR( ERR_CONSISTENCY, "flatten_leaf", "nodes in same row have different son layout" );
            }// for
        }// for

        // sum up nrows/ncols and replace values by offsets
        uint  new_nrows = 0;
        uint  new_ncols = 0;
        
        for ( uint  i = 0; i < nrows; ++i )
        {
            const auto  ofs = new_nrows;
            
            new_nrows    += son_nrows[i];
            son_nrows[i]  = ofs;
        }// for

        for ( uint  i = 0; i < ncols; ++i )
        {
            const auto  ofs = new_ncols;
            
            new_ncols    += son_ncols[i];
            son_ncols[i]  = ofs;
        }// for

        sons.resize( new_nrows * new_ncols );

        for ( uint  i = 0; i < new_nrows; ++i )
            for ( uint  j = 0; j < new_ncols; ++j )
                sons[ j * new_nrows + i ] = nullptr;

        //
        // now copy sons into new layout
        //
        
        for ( uint  i = 0; i < nrows; ++i )
        {
            for ( uint  j = 0; j < ncols; ++j )
            {
                auto  bcl = bcls[ j * nrows + i ];

                if ( bcl == nullptr )
                    continue;
                
                for ( uint  ii = 0; ii < bcl->nrows(); ++ii )
                {
                    for ( uint  jj = 0; jj < bcl->ncols(); ++jj )
                    {
                        auto  son_ij = bcl->son( ii, jj );

                        if ( son_ij == nullptr )
                            continue;
                        
                        if ( son_ij->is_leaf() )
                            found_leaf = true;

                        sons[ ( son_ncols[j] + jj ) * new_nrows + ( son_nrows[i] + ii ) ] = son_ij;

                        if ( level > 0 )
                            bcl->set_son( ii, jj, nullptr, false );
                    }// for
                }// for

                if ( level > 0 )
                    delete bcl;
            }// for
        }// for

        if ( found_leaf )
        {
            //
            // set collected sons as direct sons of root node
            // - but only if hierarchy was actually changed (level > 0)
            //

            if ( level > 0 )
            {
                root->set_layout( new_nrows, new_ncols );

                for ( uint  i = 0; i < new_nrows; ++i )
                {
                    for ( uint  j = 0; j < new_ncols; ++j )
                    {
                        root->set_son( i, j, sons[ j * new_nrows + i ], false );

                        if ( sons[ j * new_nrows + i ] != nullptr )
                            sons[ j * new_nrows + i ]->set_parent( root );
                    }// for
                }// for
            }// if
            
            break;
        }// if
        else
        {
            // proceed to next level
            bcls  = std::move( sons );
            nrows = new_nrows;
            ncols = new_ncols;

            ++level;
        }// else
    }// while
}

}// namespace Hpro
