//
// Project     : HLib
// File        : TBSPCTBuilder.cc
// Description : build clustertrees via BSP
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include <vector>
#include <set>
#include <unordered_map>

#include "list.hh"
#include "treealg.hh"

#include "hpro/base/error.hh"

#include "hpro/cluster/TBSPCTBuilder.hh"

namespace HLIB
{

using std::vector;
using std::deque;
using std::set;
using std::unordered_map;
using std::list;
using std::unique_ptr;
using std::make_unique;

namespace
{

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// local defines
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// nodeset functions
//#define FOREACH( i, set )  for ( uint i = 0; i < set.n_nodes(); i++ )
                                  
// for debugging
#define PRINT if ( false )

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// local constants
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// labels for interface and undefined marks
enum { UNDEF     = 0x0,
       LEFT      = 0x1,
       RIGHT     = 0x2,
       INTERFACE = 0x4 };

// undefined index (for permutation consistency check)
const idx_t   UNDEF_IDX          = -1;

// minimal ratio of size of son clusters
const double  CLUSTER_SIZE_RATIO = 0.2;

// minimal size for parallel calls
const size_t  MIN_PAR_SIZE       = 250000;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//
// local functions
//
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// compute avergage depth of domain sub trees
double    avg_dom_depth ( const TCluster *  node );

}// namespace anonymous

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TGeomCTBuilder (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////
//
// constructor and destructor
//

TGeomCTBuilder::TGeomCTBuilder ( const uint  an_min,
                                 const uint  amin_leaf_lvl )
{
    _n_min         = std::max( an_min, uint(1) );
    _min_leaf_lvl  = amin_leaf_lvl;
    _adjust_bb     = CFG::Cluster::adjust_bbox;
    _sort_wrt_size = CFG::Cluster::sort_wrt_size;
}

//////////////////////////////////////////////
//
// clustertree creation
//

//
// build tree out of point cloud
//
unique_ptr< TClusterTree >
TGeomCTBuilder::build ( const TCoordinate *  coord,
                        const idx_t          index_ofs ) const
{
    if ( coord == nullptr )
        return nullptr;

    //
    // setup array for son assignment and permutation
    //

    const size_t    max_dof = coord->ncoord();
    TPermutation    perm_e2i( max_dof );

    for ( size_t  i = 0; i < max_dof; ++i )
        perm_e2i[i] = UNDEF_IDX;
    
    //
    // build root of clustertree (this node ;)
    // copy indices into local list and divide
    //

    unique_ptr< TGeomCluster > root;
    data_t                     data = { coord, & perm_e2i, _n_min, _min_leaf_lvl, uint(max_dof / 2) };
    TNodeSet                   dofs( max_dof );
    TBBox                      bbox;
    TOptClusterSize            csize;
    
    for ( uint i = 0; i < max_dof; i++ )
        dofs.append( i );

    // comput bounding box of root cluster (only if not automatically adjusted)
    if ( ! _adjust_bb )
        bbox = compute_bb( dofs, data );
    
    // and start partitioning
    root = divide( dofs, 0, bbox, csize, index_ofs, data );
    
    if ( root.get() == nullptr )
        HERROR( ERR_NULL, "(TGeomCTBuilder) build", "root cluster" );
    
    //
    // consistency check: are all DoFs covered by clustertree
    //
    
    bool  all_handled = true;
    
    for ( uint i = 0; i < max_dof; i++ )
        if ( perm_e2i[ i ] == UNDEF_IDX )
        {
            perm_e2i[ i ] = 0;
            all_handled   = false;
            HERROR( ERR_CONSISTENCY, "(TGeomCTBuilder) build",
                    to_string( "DoF %d is not handled", i ) );
        }// if

    if ( ! all_handled )
    {
        root.reset( nullptr );
        return nullptr;
    }// if
    
    //
    // finally put all together in a tree
    //

    auto  permutation_e2i = make_unique< TPermutation >( perm_e2i );
    auto  permutation_i2e = make_unique< TPermutation >( perm_e2i );

    permutation_i2e->invert();
    
    auto  ct = make_unique< TClusterTree >( root.release(), permutation_e2i.release(), permutation_i2e.release() );

    return ct;
}

//
// create a leaf in a clustertree
//
unique_ptr< TGeomCluster >
TGeomCTBuilder::build_leaf ( const TNodeSet &  dofs,
                             const uint        lvl,
                             const idx_t       index_ofs,
                             const TBBox &     bbox,
                             data_t &          data ) const
{
    if ( dofs.nnodes() == 0 )
        HERROR( ERR_ARG, "(TGeomCTBuilder) build_leaf", "empty indexset" );

    if ( lvl < data.min_leaf_lvl )
        HWARNING( to_string( "(TGeomCTBuilder) build_leaf : depth of tree to low, minimal = %d, actual = %d",
                             data.min_leaf_lvl, lvl ) );

    const idx_t  lb = index_ofs;
    idx_t        ub = lb;
    
    for ( auto dof : dofs )
    {
        if ( data.perm != nullptr )
            (*data.perm)[ dof ] = ub++;
    }// for

    TBBox  cl_bbox( bbox );
    
    if ( _adjust_bb )
        cl_bbox = compute_bb( dofs, data );
    
    update_bb( dofs, cl_bbox, data );

    if ( cl_bbox.min().dim() == 0 )
        HERROR( ERR_DIM, "(TGeomCTBuilder) build_leaf", "no bbox" );
    
    if ( size_t( ub - lb ) != dofs.nnodes() )
        HERROR( ERR_CONSISTENCY, "(TGeomCTBuilder) build_leaf", "missing nodes" );
    
    return make_unique< TGeomCluster >( lb, ub-1, cl_bbox );
}

namespace 
{

//
// return bbox of support of given node
//
TBBox
support_size ( const node_t         node,
               const bool           only_idx,
               const TCoordinate *  coord )
{
    TBBox       bbox;
    const uint  dim = coord->dim();

    bbox.min().set_dim( dim );
    bbox.max().set_dim( dim );
    
    // start with coordinate of index
    for ( uint j = 0; j < dim; j++ )
    {
        bbox.min()[j] = bbox.max()[j] = coord->coord( node )[j];

        // and add corresponding bounding box data if available
        if ( ! only_idx && coord->has_bbox() )
        {
            bbox.min()[j] = std::min( bbox.min()[j], coord->bbmin( node )[j] );
            bbox.max()[j] = std::max( bbox.max()[j], coord->bbmax( node )[j] );
        }// if
    }// for

    return bbox;
}

//
// compute bbox of <dofs> in parallel
//
TBBox
par_compute_bb ( const TNodeSet &  dofs,
                 const size_t      lb,
                 const size_t      ub,
                 const TCoordinate *  coord )
{
    if ( ub - lb > 10000 )
    {
        const size_t  mid = (ub + lb) / 2;
        TBBox  bbox1, bbox2;

        bbox1 = par_compute_bb( dofs, lb, mid, coord );
        bbox2 = par_compute_bb( dofs, mid, ub, coord );

        bbox1.join( bbox2 );

        return  bbox1;
    }// if
    else
    {
        TBBox       bbox;
        const uint  dim = coord->dim();
    
        bbox.min().set_dim( dim );
        bbox.max().set_dim( dim );
    
        bbox = support_size( dofs[lb], true, coord );
    
        for ( size_t  i = lb; i < ub; ++i )
            bbox.join( support_size( dofs[i], true, coord ) );

        return  bbox;
    }// else
}

}// namespace anonymous

//
// compute bounding box of a cluster
//
TBBox
TGeomCTBuilder::compute_bb ( const TNodeSet &  dofs,
                             const data_t &    data ) const
{
    return  par_compute_bb( dofs, 0, dofs.nnodes(), data.coord );

    // TBBox       bbox;
    // const uint  dim = data.coord->dim();
    
    // bbox.min().set_dim( dim );
    // bbox.max().set_dim( dim );
    
    // bbox = support_size( dofs[0], true, data.coord );

    // for ( auto  dof : dofs )
    //     bbox.join( support_size( dof, true, data.coord ) );

    // return bbox;
}

//
// update bounding box of cluster
//
void
TGeomCTBuilder::update_bb ( const TNodeSet &  dofs,
                            TBBox &           bbox,
                            const data_t &    data ) const
{
    //
    // update bbox with support of indices
    //
    
    for ( auto  dof : dofs )
        bbox.join( support_size( dof, false, data.coord ) );

    check_bb( bbox, data );
}

//
// check bbox for degenerate axis
//
void
TGeomCTBuilder::check_bb ( TBBox &         bbox,
                           const data_t &  data ) const
{
    //
    // check for degenerate axis and enhance box if necessary,
    // e.g. if axis-length is too small compared with diameter
    //

    const double BBOX_EPS = 1e-5;
    double       diam     = Math::square( bbox.min()[0] - bbox.max()[0] );
    const uint   dim      = data.coord->dim();
    
    for ( uint i = 1; i < dim; i++ )
        diam += Math::square( bbox.min()[i] - bbox.max()[i] );

    diam = Math::sqrt( diam );

    for ( uint i = 0; i < dim; i++ )
    {
        const double diff = bbox.max()[i] - bbox.min()[i];
        
        if ( diff < BBOX_EPS * diam )
        {
            bbox.min()[i] -= 0.5 * ( BBOX_EPS * diam - diff );
            bbox.max()[i] += 0.5 * ( BBOX_EPS * diam - diff );
        }// if
    }// for
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TBSPCTBuilder (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////
//
// constructor and destructor
//

TBSPCTBuilder::TBSPCTBuilder ( const TBSPPartStrat *  part_strat,
                               const uint             anmin,
                               const uint             amin_leaf_lvl )
        : TGeomCTBuilder( anmin, amin_leaf_lvl )
        , _part_strat( part_strat ) 
{
    if ( _part_strat == nullptr )
        HERROR( ERR_ARG, "(TBSPCTBuilder) ctor", "partition strategy is nullptr" );
}

TBSPCTBuilder::~TBSPCTBuilder ()
{
}

//////////////////////////////////////////////
//
// clustertree creation
//

//
// divide given cluster into sons
//
unique_ptr< TGeomCluster >
TBSPCTBuilder::divide ( const TNodeSet &         dofs,
                        const uint               lvl,
                        const TBBox &            bbox,
                        const TOptClusterSize &  csize,
                        const idx_t              index_ofs,
                        data_t &                 data ) const
{
    //
    // check if number of indices in this cluster is too small to divide
    // (and minimal leaf level was reached)
    //

    if (( lvl >= data.min_leaf_lvl ) && ( dofs.nnodes() <= data.nmin ))
        return build_leaf( dofs, lvl, index_ofs, bbox, data );

    if ( lvl > data.max_lvl )
    {
        HWARNING( to_string( "in (TBSPCTBuilder) divide : maximal tree depth reached; depth = %d", lvl ) );
        return build_leaf( dofs, lvl, index_ofs, bbox, data );
    }// if
    
    //
    // if size of cluster leaves is at most optimal size don't split cluster
    //

    if ( ! csize.is_optimal( dofs.nnodes() ) )
    {
        auto  son = divide( dofs, lvl+1, bbox, csize.recurse(), index_ofs, data );

        if ( son.get() == nullptr )
            HERROR( ERR_NULL, "(TBSPCTBuilder) divide", "son cluster" );

        auto  cluster = make_unique< TGeomCluster >( son->first(), son->last(), son->bbox() );
        
        cluster->set_nsons( 1 );
        cluster->set_son( 0, son.release() );

        return cluster;
    }// if
    
    //
    // divide sons
    //

    TBBox            cl_bbox( bbox );
    TNodeSet         son_dofs[2];
    vector< TBBox >  son_bbox( 2 );

    if ( _adjust_bb )
        cl_bbox = compute_bb( dofs, data );

    _part_strat->partition( data.coord, dofs, son_dofs[0], son_dofs[1], cl_bbox, son_bbox, lvl );
    
    HDEBUG( to_string( "(TBSPCTBuilder) divide : left = %d, right = %d",
                       son_dofs[0].nnodes(), son_dofs[1].nnodes() ) );
    
    // if one of the subsets is empty => continue with this cluster
    if      ( son_dofs[0].nnodes() == 0 ) return divide( son_dofs[1], lvl+1, son_bbox[1], csize.recurse(), index_ofs, data );
    else if ( son_dofs[1].nnodes() == 0 ) return divide( son_dofs[0], lvl+1, son_bbox[0], csize.recurse(), index_ofs, data );

    //
    // sort sets according to size:
    //   if size differs much, move smaller set to the right
    //

    if ( _sort_wrt_size )
    {
        if ( double(son_dofs[0].nnodes()) / double(son_dofs[1].nnodes()) < CLUSTER_SIZE_RATIO )
        {
            std::swap( son_dofs[0], son_dofs[1] );
            std::swap( son_bbox[0], son_bbox[1] );
        }// if
    }// if
         
    // consistency check
    if ( son_dofs[0].nnodes() + son_dofs[1].nnodes() != dofs.nnodes() )
        HERROR( ERR_CONSISTENCY, "(TBSPCTBuilder) divide", "lost nodes during divide" );
    
    //
    // recursive call for building clustertrees with sons
    //

    unique_ptr< TGeomCluster >  sons[2];
    const idx_t                 left_ofs  = index_ofs;
    const idx_t                 right_ofs = index_ofs + idx_t( son_dofs[0].nnodes() );
        

    auto  build_soncl =
        [this,lvl,&csize,&data] ( TNodeSet &     sdofs,
                                  const TBBox &  sbbox,
                                  const idx_t    sofs ) -> unique_ptr< TGeomCluster >
        {
            auto  cl = divide( sdofs, lvl+1, sbbox, csize.recurse(), sofs, data );
                
            if ( cl.get() == nullptr )
                HERROR( ERR_NULL, "(TBSPCTBuilder) divide", "son cluster" );
                
            sdofs.remove_all();
            sdofs.resize( 0 );

            return cl;
        };

    sons[0] = build_soncl( son_dofs[0], son_bbox[0], left_ofs  );
    sons[1] = build_soncl( son_dofs[1], son_bbox[1], right_ofs );

    //
    // set bounding box of cluster
    //

    if ( _adjust_bb )
    {
        //
        // as union of bb of sons
        //

        cl_bbox = sons[0]->bbox();
        cl_bbox.join( sons[1]->bbox() );
    }// if
    else
    {
        //
        // as provided by argument and updated with son-bboxes
        //
        
        cl_bbox.join( sons[0]->bbox() );
        cl_bbox.join( sons[1]->bbox() );
    }// else

    check_bb( cl_bbox, data );

    //
    // finally build cluster
    //
    
    auto  cluster = make_unique< TGeomCluster >( std::min( sons[0]->first(), sons[1]->first() ),
                                                 std::max( sons[0]->last(),  sons[1]->last()  ),
                                                 cl_bbox );

    cluster->set_nsons( 2 );
    cluster->set_son( 0, sons[0].release() );
    cluster->set_son( 1, sons[1].release() );
    
    return cluster;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//
// TBSPNDCTBuilder (implementation)
//
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////
//
// constructor and destructor
//

TBSPNDCTBuilder::TBSPNDCTBuilder ( const TSparseMatrix *  S,
                                   const TBSPPartStrat *  part_strat,
                                   const uint             an_min,
                                   const uint             amin_leaf_lvl )
        : TBSPCTBuilder( part_strat, an_min, amin_leaf_lvl )
        , _sparse_mat( S )
        , _sync_interface_depth( CFG::Cluster::sync_interface_depth )
{}

TBSPNDCTBuilder::TBSPNDCTBuilder ( const TBSPPartStrat *  part_strat,
                                   const uint             an_min,
                                   const uint             amin_leaf_lvl )
        : TBSPCTBuilder( part_strat, an_min, amin_leaf_lvl )
        , _sparse_mat( nullptr )
        , _sync_interface_depth( CFG::Cluster::sync_interface_depth )
{}

TBSPNDCTBuilder::~TBSPNDCTBuilder ()
{
}

//////////////////////////////////////////////
//
// clustertree creation
//

//
// build tree out of point cloud with additional information
// about connectivity provided by a sparse matrix
//
unique_ptr< TClusterTree >
TBSPNDCTBuilder::build ( const TCoordinate * coord,
                         const idx_t         index_ofs ) const
{
    return build( coord, _sparse_mat, index_ofs );
}

unique_ptr< TClusterTree >
TBSPNDCTBuilder::build ( const TCoordinate *    coord,
                         const TSparseMatrix *  S,
                         const idx_t            index_ofs ) const
{
    if ( coord == nullptr )
        return nullptr;
    
    if ( S == nullptr )
        HERROR( ERR_ARG, "(TBSPNDCTBuilder) build", "sparse matrix is nullptr" );
    
    //
    // setup array for son assignment and permutation
    //

    const size_t    max_dof = coord->ncoord();
    TPermutation    perm( max_dof );

    for ( size_t  i = 0; i < max_dof; ++i )
        perm[i] = UNDEF_IDX;
    
    //
    // count number of edges per node and total number of edges
    //
    
    const size_t    nnodes = max_dof;
    vector< uint >  nedges( nnodes, 0 );
    size_t          nentries = 0;
    
    for ( node_t  node = 0; node < node_t(nnodes); ++node )
    {
        const idx_t  lb = S->rowptr(node);
        const idx_t  ub = S->rowptr(node+1);
            
        for ( idx_t  j = lb; j < ub; j++ )
        {
            const node_t  neigh = S->colind(j);

            // no single loops
            if ( node == neigh )
                continue;
                
            nedges[node]++;
            nentries++;
        }// for
    }// for
    
    //
    // look for highly connected nodes and remove them
    //
    
    const size_t    avg_degree = nentries / nnodes;
    const size_t    max_degree = 5 * avg_degree;
    list< node_t >  high_degree_nodes;
    vector< bool >  node_del( nnodes, false );

    for ( node_t  node = 0; node < node_t(nnodes); ++node )
    {
        if ( nedges[node] > max_degree )
        {
            HINFO( to_string( "(TBSPNDCTBuilder) build : high degree node : %d (degree = %d vs. %d average)",
                              node, avg_degree ) );
            high_degree_nodes.push_back( node );
            node_del[node] = true;
        }// if
    }// for

    //
    // build root of clustertree (this node ;)
    // copy indices into local list and divide
    //

    unique_ptr< TGeomCluster > root;
    data_t                     data = { coord, & perm, _n_min, _min_leaf_lvl, uint(max_dof / 2) };
    TNodeSet                   dofs( max_dof );
    TBBox                      bbox;
    TOptClusterSize            csize;

    for ( size_t  i = 0; i < max_dof; i++ )
    {
        if ( node_del[i] )
            continue;
        
        dofs.append( node_t(i) );
    }// for

    // comput bounding box of root cluster
    bbox = compute_bb( dofs, data );

    // and start partitioning
    root = divide( dofs, 0, bbox, csize.recurse(), index_ofs, data );
    
    if ( root.get() == nullptr )
        HERROR( ERR_NULL, "(TBSPNDCTBuilder) build", "root cluster" );
    
    //
    // add highly connected nodes as another cluster and extend tree
    //

    if ( ! high_degree_nodes.empty() )
    {
        const idx_t  lb = index_ofs + idx_t( root->size() );
        idx_t        ub = lb;

        dofs.remove_all();
        
        while ( ! high_degree_nodes.empty() )
        {
            const node_t  node = behead( high_degree_nodes );

            dofs.append( node );
            perm[ node ] = ub++;
        }// while
        
        bbox = compute_bb( dofs, data );
        update_bb(  dofs, bbox, data );
        
        auto  high_cl  = make_unique< TGeomCluster >( lb, ub-1 );
        auto  new_root = make_unique< TGeomCluster >( index_ofs, ub-1 );

        high_cl->bbox() = bbox;
        
        new_root->set_nsons( 2 );
        new_root->set_son( 0, root.release() );
        new_root->set_son( 1, high_cl.release() );

        // assign new root
        root = std::move( new_root );
    }// if
    
    //
    // consistency check: are all DoFs covered by clustertree
    //
    
    bool  all_handled = true;
    
    for ( size_t  i = 0; i < max_dof; i++ )
        if ( perm[ i ] == UNDEF_IDX )
        {
            perm[ i ]   = 0;
            all_handled = false;
            HERROR( ERR_CONSISTENCY, "(TBSPNDCTBuilder) build",
                    to_string( "DoF %d is not handled", i ) );
        }// if

    if ( ! all_handled )
    {
        root.reset( nullptr );
        return nullptr;
    }// if
    
    //
    // finally put all together in a tree
    //

    auto  perm_e2i = make_unique< TPermutation >( perm );
    auto  perm_i2e = make_unique< TPermutation >( perm );

    perm_i2e->invert();
    
    return make_unique< TClusterTree >( root.release(), perm_e2i.release(), perm_i2e.release() );
}

//
// divide given cluster into sons
//
unique_ptr< TGeomCluster >
TBSPNDCTBuilder::divide ( const TNodeSet &         dofs,
                          const uint               lvl,
                          const TBBox &            bbox,
                          const TOptClusterSize &  csize,
                          const idx_t              index_ofs,
                          data_t &                 data ) const
{
    //
    // check if number of indices in this cluster is too small to divide
    // (and minimal leaf level was reached)
    //
    
    if (( lvl >= data.min_leaf_lvl ) && ( dofs.nnodes() <= data.nmin ))
        return build_leaf( dofs, lvl, index_ofs, bbox, data );

    if ( lvl > data.max_lvl )
    {
        HWARNING( to_string( "in (TBSPNDPCTBuilder) divide : maximal tree depth reached; depth = %d", lvl ) );
        return build_leaf( dofs, lvl, index_ofs, bbox, data );
    }// if

    //
    // divide sons
    //

    TBBox            cl_bbox( bbox );
    TNodeSet         son_dofs[2];
    vector< TBBox >  son_bbox( 2 );

    if ( _adjust_bb )
        cl_bbox = compute_bb( dofs, data );
        
    _part_strat->partition( data.coord, dofs, son_dofs[0], son_dofs[1], cl_bbox, son_bbox, lvl );
    
    // if one of the subsets is empty => continue with this cluster
    if      ( son_dofs[0].nnodes() == 0 ) return divide( son_dofs[1], lvl+1, son_bbox[1], csize.recurse(), index_ofs, data );
    else if ( son_dofs[1].nnodes() == 0 ) return divide( son_dofs[0], lvl+1, son_bbox[0], csize.recurse(), index_ofs, data );

    //
    // if an interface should be build, remove indices from
    // son-sets until they are decoupled, e.g. no index-pair in
    // different clusters is connected by a common cell
    //

    TNodeSet  if_dofs;

    //
    // mark local nodes and choose set with most dofs first
    //
    
    const uint                    max_son = ( son_dofs[0].nnodes() > son_dofs[1].nnodes() ? 0 : 1 );
    unordered_map< idx_t, bool >  local;

    // mark local nodes
    for ( auto  dof : dofs ) local[ dof ] = true;
    
    if_dofs.resize( son_dofs[max_son].size() );
    
    //
    // detect interface by comparing labels of neighbours;
    // if they differ, node is on interface;
    // do this for both sets but start with larger one
    //

    unordered_map< idx_t, char >  label;

    for ( auto  dof : son_dofs[0] ) label[ dof ] = LEFT;
    for ( auto  dof : son_dofs[1] ) label[ dof ] = RIGHT;
    
    for ( uint l = 0; l < 2; l++ )
    {
        idx_t  idx;
        
        if ( l == 0 ) idx = max_son;
        else          idx = (max_son == 0 ? 1 : 0);
        
        for ( auto  dof : son_dofs[idx] )
        {
            if ( label[ dof ] == INTERFACE )
                continue;
                
            const idx_t  lb = _sparse_mat->rowptr( dof   ); // TODO: replace by options.S
            const idx_t  ub = _sparse_mat->rowptr( dof+1 );
                
            for ( idx_t j = lb; j < ub; j++ )
            {
                const idx_t  neigh = _sparse_mat->colind(j);
                    
                if (( dof != neigh ) &&
                    ( local[ neigh ] ) &&   // neighbour has to be local node
                    // ( local.find( neigh ) != local.end() ) &&   // neighbour has to be local node
                    ( label[ neigh ] != label[ dof ] ) &&       // then if labels differ
                    ( label[ neigh ] != INTERFACE ))            // and neighbour is not on interface
                {
                    //
                    // remove dof from son list and put it into interface
                    //
                    
                    label[ dof ] = INTERFACE;
                    if_dofs.append( dof );
                    break;
                }// if
            }// for
        }// for
    }// for

    //
    // rebuild list of reduced indexsets
    //

    son_dofs[0].remove_all();
    son_dofs[1].remove_all();
    
    for ( auto dof : dofs )
    {
        if      ( label[ dof ] == LEFT  ) son_dofs[0].append( dof );
        else if ( label[ dof ] == RIGHT ) son_dofs[1].append( dof );
    }// for

    // recheck empty sons
    if ( son_dofs[0].nnodes() == 0 ) return build_leaf( dofs, lvl, index_ofs, cl_bbox, data );
    if ( son_dofs[1].nnodes() == 0 ) return build_leaf( dofs, lvl, index_ofs, cl_bbox, data );

    // consistency check
    if ( son_dofs[0].nnodes() + son_dofs[1].nnodes() + if_dofs.nnodes() != dofs.nnodes() )
        HERROR( ERR_CONSISTENCY, "(TBSPNDCTBuilder) divide", "lost nodes during divide" );
    
    //
    // recursive call for building clustertrees with sons
    //

    unique_ptr< TGeomCluster >  sons[2];
    double                      son_lvl[2] = { 0, 0 };
    const idx_t                 left_ofs  = index_ofs;
    const idx_t                 right_ofs = index_ofs + idx_t( son_dofs[0].nnodes() );

    auto  build_soncl =
        [this,lvl,&csize,&data] ( TNodeSet &     sdofs,
                                  const TBBox &  sbbox,
                                  const idx_t    sofs,
                                  double &       slvl ) -> unique_ptr< TGeomCluster >
        {
            if ( sdofs.nnodes() > 0 )
            {
                auto  cl = divide( sdofs, lvl+1, sbbox, csize.recurse(), sofs, data );

                if ( cl.get() == nullptr )
                    HERROR( ERR_NULL, "(TBSPNDCTBuilder) divide", "son cluster" );
        
                slvl = avg_dom_depth( cl.get() );

                sdofs.remove_all();
                sdofs.resize( 0 );

                return cl;
            }// if

            return nullptr;
        };

    sons[0] = build_soncl( son_dofs[0], son_bbox[0], left_ofs,  son_lvl[0] );
    sons[1] = build_soncl( son_dofs[1], son_bbox[1], right_ofs, son_lvl[1] );

    if (( sons[0].get() == nullptr ) || ( sons[1].get() == nullptr ))
        HERROR( ERR_CONSISTENCY, "(TBSPNDCTBuilder) divide", "missing son cluster" );
        
    //
    // build interface cluster
    //
    
    unique_ptr< TGeomCluster >  if_son;
    
    if ( if_dofs.nnodes() > 0 )
    {
        const idx_t   if_ofs    = index_ofs + idx_t(sons[0]->size()) + idx_t(sons[1]->size());
        double        reduction = 0.0;
        TBBox         if_bbox;
        const double  avg_lvl   = std::max( son_lvl[0], son_lvl[1] );

        // compute optimal reduction of clustersize to achieve desired depth
        // reduction = (n_min / |if|)^(1/max_lvl)
        if ( avg_lvl > 0.0 )
            reduction = Math::pow( double(std::max( data.nmin, 1u )) / ( double(if_dofs.nnodes()) ),
                                   1.0 / avg_lvl );

        if_bbox = compute_bb( if_dofs, data );
        
        if_son = divide_if( if_dofs, lvl+1, std::min( uint(lvl + 1 + avg_lvl), data.max_lvl ), if_bbox,
                            TOptClusterSize( if_dofs.nnodes(), reduction ), if_ofs, data );
        
        if ( if_son == nullptr )
            HERROR( ERR_NULL, "(TBSPNDCTBuilder) divide", "interface cluster" );
    }// if

    // mark domain sons as such
    sons[0]->set_domain( true );
    sons[1]->set_domain( true );
    
    //
    // set bounding box of cluster
    //

    if ( _adjust_bb )
    {
        //
        // as union of bb of sons
        //

        cl_bbox = sons[0]->bbox();

        cl_bbox.join( sons[1]->bbox() );

        if ( if_son.get() != nullptr )
            cl_bbox.join( if_son->bbox() );
    }// if
    else
    {
        //
        // as provided by argument and updated with son-bboxes
        //
        
        cl_bbox.join( sons[0]->bbox() );
        cl_bbox.join( sons[1]->bbox() );

        if ( if_son.get() != nullptr )
            cl_bbox.join( if_son->bbox() );
    }// else

    check_bb( cl_bbox, data );
    
    //
    // create indexset of cluster as union of son-iss
    //

    TIndexSet  is;

    if ( if_son == nullptr )
    {
        is.set_first_last( std::min( sons[0]->first(), sons[1]->first() ),
                           std::max( sons[0]->last(),  sons[1]->last()  ) );
    }// if
    else
    {
        is.set_first_last( std::min( std::min( sons[0]->first(), sons[1]->first() ), if_son->first() ),
                           std::max( std::max( sons[0]->last(),  sons[1]->last() ),  if_son->last()  ) );
    }// else
    
    auto  cluster = make_unique< TGeomCluster >( is.first(), is.last(), cl_bbox );

    if ( if_son.get() != nullptr )
        cluster->set_nsons( 3 );
    else
        cluster->set_nsons( 2 );
    
    cluster->set_son( 0, sons[0].release() );
    cluster->set_son( 1, sons[1].release() );

    if ( if_son.get() != nullptr )
        cluster->set_son( 2, if_son.release() );

    return cluster;
}

//
// recursively build cluster for interfaces
//
unique_ptr< TGeomCluster >
TBSPNDCTBuilder::divide_if ( const TNodeSet &         dofs,
                             const uint               lvl,
                             const uint               max_lvl,
                             const TBBox &            bbox,
                             const TOptClusterSize &  csize,
                             const idx_t              index_ofs,
                             data_t &                 data ) const
{
    auto  recurse =
        [this,lvl,max_lvl,index_ofs,&dofs,&bbox,&csize,&data] () -> unique_ptr< TGeomCluster >
        {
            auto  son = divide_if( dofs, lvl+1, max_lvl, bbox, csize.recurse(), index_ofs, data );

            if ( son.get() == nullptr )
                HERROR( ERR_NULL, "(TBSPNDCTBuilder) divide_if", "son cluster" );

            auto  cluster = make_unique< TGeomCluster >( son->first(), son->last(), son->bbox() );

            cluster->set_nsons( 1 );
            cluster->set_son( 0, son.release() );
            
            return cluster;
        };
        
    //
    // check if number of indices in this cluster
    // is too small to divide
    //
    
    if (( dofs.nnodes() <= data.nmin ) || ( lvl >= max_lvl ))
        return build_leaf( dofs, lvl, index_ofs, bbox, data );

    //
    // if size of cluster leaves is at most optimal size don't split cluster
    //

    if ( _sync_interface_depth && ! csize.is_optimal( dofs.nnodes() ) )
        return recurse();
    
    //
    // divide sons
    //

    TBBox           cl_bbox( bbox );
    TNodeSet        son_dofs[2];
    vector< TBBox > son_bbox( 2 );

    if ( _adjust_bb )
        cl_bbox = compute_bb( dofs, data );
        
    _part_strat->partition( data.coord, dofs, son_dofs[0], son_dofs[1], cl_bbox, son_bbox, lvl );
    
    // if one of the subsets is empty => recurse
    if (( son_dofs[0].nnodes() == 0 ) || ( son_dofs[1].nnodes() == 0 ))
        return recurse();
    
    // consistency check
    if ( son_dofs[0].nnodes() + son_dofs[1].nnodes() != dofs.nnodes() )
    {
        _part_strat->partition( data.coord, dofs, son_dofs[0], son_dofs[1], cl_bbox, son_bbox, lvl );
        HERROR( ERR_CONSISTENCY, "(TBSPNDCTBuilder) divide_if", "lost nodes during divide" );
    }// if
    
    //
    // recursive call for building clustertrees with sons
    //

    unique_ptr< TGeomCluster >  sons[2];
    const idx_t                 left_ofs  = index_ofs;
    const idx_t                 right_ofs = index_ofs + idx_t( son_dofs[0].nnodes() );

    auto  build_soncl =
        [this,lvl,max_lvl,&csize,&data] ( TNodeSet &     sdofs,
                                          const TBBox &  sbbox,
                                          const idx_t    sofs ) -> unique_ptr< TGeomCluster >
        {
            if ( sdofs.nnodes() > 0 )
            {
                auto  cl = divide_if( sdofs, lvl+1, max_lvl, sbbox, csize.recurse(), sofs, data );
                
                if ( cl.get() == nullptr )
                    HERROR( ERR_NULL, "(TBSPNDCTBuilder) divide_if", "son cluster is nullptr" );

                sdofs.remove_all();
                sdofs.resize( 0 );

                return cl;
            }// if

            return nullptr;
        };

    sons[0] = build_soncl( son_dofs[0], son_bbox[0], left_ofs  );
    sons[1] = build_soncl( son_dofs[1], son_bbox[1], right_ofs );

    //
    // set bounding box of cluster
    //

    if ( _adjust_bb )
    {
        //
        // as union of bb of sons
        //
        
        cl_bbox = sons[0]->bbox();
        cl_bbox.join( sons[1]->bbox() );
    }// if
    else
    {
        //
        // as provided by argument and updated with son-bboxes
        //
        
        cl_bbox.join( sons[0]->bbox() );
        cl_bbox.join( sons[1]->bbox() );
    }// else

    check_bb( cl_bbox, data );
    
    //
    // create indexset of cluster as union of son-iss
    //
    
    auto  cluster = make_unique< TGeomCluster >( std::min( sons[0]->first(), sons[1]->first() ),
                                                 std::max( sons[0]->last(),  sons[1]->last()  ),
                                                 cl_bbox );

    cluster->set_nsons( 2 );
    cluster->set_son( 0, sons[0].release() );
    cluster->set_son( 1, sons[1].release() );

    return cluster;
}

namespace
{

//
// compute avergage depth of domain sub trees
//
double
avg_dom_depth ( const TCluster *  node )
{
    using  cl_list_t = list< const TCluster * >;

    size_t     depth_sum = 0;  // sum of depth of all sub trees
    size_t     n_leafs   = 0;  // number of leafs
    cl_list_t  nodes;
    size_t     lvl = 0;

    nodes.push_back( node );

    while ( ! nodes.empty() )
    {
        cl_list_t  sons;

        while ( ! nodes.empty() )
        {
            const TCluster *  cl = behead( nodes );

            if ( cl->is_leaf() )
            {
                depth_sum += lvl;
                ++n_leafs;
            }// if
            else
            {
                for ( uint  i = 0; i < cl->nsons(); ++i )
                {
                    const TCluster *  son_i = cl->son( i );

                    if (( son_i == nullptr ) || ! son_i->is_domain() )
                        continue;

                    sons.push_back( son_i );
                }// for
            }// else
        }// while

        nodes = sons;
        ++lvl;
    }// while

    // average size is algebraic mean
    if ( n_leafs > 0 )
        return double( depth_sum ) / double( n_leafs );
    else
        return 0;
}

}// namespace anonymous

}// namespace HLIB
