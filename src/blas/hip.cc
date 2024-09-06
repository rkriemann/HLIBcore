//
// Project     : HLIBpro
// File        : blas/hip.cc
// Description : HIP related functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2024. All Rights Reserved.
//

#include <array>
#include <mutex>

#include <hpro/config.h>

#if HPRO_USE_HIP == 1
#  ifdef __HIP_PLATFORM_NVIDIA__
#    include <cuda_runtime.h>
#  endif

#  include <hipblas/hipblas.h>
#  include <hipsolver/hipsolver.h>
#endif

#include <hpro/blas/Algebra.hh>
#include <hpro/blas/hip.hh>

namespace Hpro { namespace HIP {

////////////////////////////////////////////////////////////////////////////////
//
// management functions
//
////////////////////////////////////////////////////////////////////////////////

#if HPRO_USE_HIP == 1

// wrapper for hip, cuBlas and hipsolver functions
#define HPRO_HIP_CHECK( func, args )                                    \
    {                                                                   \
        auto  result = func args ;                                      \
        if ( result != hipSuccess )                                     \
        {                                                               \
            HERROR( ERR_HIP, #func, Hpro::to_string( "HIP code: %d", result ) ); \
        }                                                               \
    }

#define HPRO_HIPSOLVER_CHECK( func, args )                              \
    {                                                                   \
        auto  result = func args ;                                      \
        if ( result != HIPSOLVER_STATUS_SUCCESS )                       \
        {                                                               \
            HERROR( ERR_HIP, #func, Hpro::to_string( "HIPSOLVER code: %d", result ) ); \
        }                                                               \
    }

//
// mapping of default types to cuBLAS types
//
template < typename T > struct hip_traits
{
    using  hip_type = T;
    using  ptr_type  = T *;
};

template <>
struct hip_traits< float >
{
    using  hip_type = float;
    using  ptr_type = float *;
};

template <>
struct hip_traits< double >
{
    using  hip_type = double;
    using  ptr_type = double *;
};

template <>
struct hip_traits< std::complex< float > >
{
    using  hip_type = cuFloatComplex;
    using  ptr_type = cuFloatComplex *;
};

template <>
struct hip_traits< std::complex< double > >
{
    using  hip_type = cuDoubleComplex;
    using  ptr_type = cuDoubleComplex *;
};

//
// joined handle for cuBLAS and hipsolverDn
//
struct handle_t
{
    hipStream_t        stream;
    hipsolverHandle_t  solver;
};

#endif

namespace
{

#if HPRO_USE_HIP == 1

//
// signals availability of HIP
//
bool  has_hip = false;

//
// concurrent streams
//
constexpr int                        max_handles = 16;
std::array< handle_t, max_handles >  handles;
std::array< bool, max_handles >      handle_used;
std::mutex                           handle_mtx;

#endif

}// namespace anonymous

//
// initialize HIP related data structures
//
void
init ()
{
    #if HPRO_USE_HIP == 1

    //
    // first check number of devices
    //
    
    int  ndev = 0;
    
    if ( hipGetDeviceCount( & ndev ) != hipSuccess )
        return;

    HINFO( to_string( "(HIP) init : found %d HIP device(s)", ndev ) );
    
    if ( ndev == 0 )
        return;

    // get number of compute units
    int  dev_id = 0;
    
    if ( hipGetDeviceProperties( & dev_prop, dev_id ) != hipSuccess )
        return;

    const auto  cu_count = dev_prop.multiProcessorCount;
    
    //
    // now try to set up streams
    // - 4 CUs per stream (try to keep 32 divisible by it)
    //

    uint32_t  cu_num      = 4;
    uint32_t  cu_pattern  = 0b1111; // four CUs
    uint32_t  cu_mask[ 8 ];         // place for 8 * 32 compute units
    uint32_t  cu_idx      = 0;      // start position for active CUs in mask
    bool      found_error = false;
    uint      error_index = 0;
    
    for ( uint  i = 0; i < max_handles; ++i )
    {
        const auto  arr_ofs = cu_idx / 32;
        const auto  bit_ofs = ( cu_idx - arr_ofs * 32 );

        cu_mask[ arr_ofs ] = cu_pattern << bit_ofs;
        cu_idx            += cu_num;
        
        if ( hipExtStreamCreateWithCUMask( & handle.stream, 8, cu_mask ) != hipSuccess )
        {
            HWARNING( "(HIP) init : error in hipExtStreamCreateWithCUMask; disabling HIP" );
            found_error = true;
            error_index = i;
            break;
        }// if

        // after this, assumption is that HIP is present and working!!!
        
        HPRO_HIPSOLVER_CHECK( hipsolverDnCreate,    ( & handles[i].solver ) );
        HPRO_HIPSOLVER_CHECK( hipsolverDnSetStream, (   handles[i].solver, handles[i].stream ) );

        handle_used[i] = false;
    }// for

    if ( found_error )
    {
        // clean up all initialized data
        for ( uint  i = 0; i < error_index; ++i )
        {
            HPRO_HIPSOLVER_CHECK( hipsolverDnDestroy,  ( handles[i].solver ) );
            HPRO_CUBLAS_CHECK(   cublasDestroy,      ( handles[i].blas ) );
            HPRO_HIP_CHECK(     hipStreamDestroy,  ( handles[i].stream ) );
        }// for
    }// if
    else
    {
        has_hip = true;
    }// else

    #endif
}

//
// finish HIP related data structures
//
void
done ()
{
    #if HPRO_USE_HIP == 1

    if ( ! has_hip )
        return;

    has_hip = false;
    
    for ( uint  i = 0; i < max_handles; ++i )
    {
        if ( handle_used[i] )
            HERROR( ERR_HIP, "(HIP) done", "handle still in use" );
        
        HPRO_HIPSOLVER_CHECK( hipsolverDnDestroy,  ( handles[i].solver ) );
        HPRO_CUBLAS_CHECK(   cublasDestroy,      ( handles[i].blas ) );
        HPRO_HIP_CHECK(     hipStreamDestroy,  ( handles[i].stream ) );
    }// for

    #endif
}

//
// get free handle for HIP (stream/cuBlas/hipsolver)
//
hip_handle_t
request_handle ()
{
    #if HPRO_USE_HIP == 1
    
    if ( ! has_hip )
        return INVALID_HANDLE;
    
    auto  lock = std::scoped_lock( handle_mtx );

    for ( uint  i = 0; i < max_handles; ++i )
    {
        if ( ! handle_used[i] )
        {
            handle_used[i] = true;
            return i;
        }// if
    }// for

    #endif
    
    return INVALID_HANDLE;
}

//
// release previously aquired handle
//
#if HPRO_USE_HIP == 1
void
release_handle ( const hip_handle_t  handle )
{
    if ( handle == INVALID_HANDLE )
        return;

    auto  lock = std::scoped_lock( handle_mtx );

    if (( handle < 0 ) || ( handle >= max_handles ))
        HERROR( ERR_HIP, "(HIP) release_handle", "unknown handle given" );

    if ( ! handle_used[ handle ] )
        HERROR( ERR_HIP, "(HIP) release_handle", "given handle not in use" );

    handle_used[ handle ] = false;

}
#else
void release_handle ( const hip_handle_t  ) {}
#endif

////////////////////////////////////////////////////////////////////////////////
//
// device handling
//
////////////////////////////////////////////////////////////////////////////////

#if HPRO_USE_HIP == 1

namespace
{

//
// local breakpoint
//
void
hip_break ()
{
}

//
// allocate device memory
//
template < typename value_t >
value_t *
device_alloc ( const size_t  n )
{
    value_t *  ptr    = nullptr;
    auto       retval = hipMalloc( & ptr, n * sizeof(value_t) );

    if ( retval != hipSuccess )
    {
        hip_break();
        HWARNING( "(HIP) device_alloc : error in \"hipMalloc\" : " + std::string( hipGetErrorString( retval ) ) );
        return nullptr;
    }// if

    return  ptr;
}

//
// deallocate device memory
//
template < typename value_t >
void
device_free ( value_t *  ptr )
{
    const auto  retval = hipFree( ptr );

    if ( retval != hipSuccess )
    {
        hip_break();
        HWARNING( "(HIP) device_alloc : error in \"hipFree\" : " + std::string( hipGetErrorString( retval ) ) );
    }// if
}

//
// copy data to device
//
template < typename value_t >
bool
to_device ( handle_t                                   handle,
            const BLAS::Matrix< value_t > &            M_host,
            typename hip_traits< value_t >::ptr_type  M_dev )
{
    const auto  retval = hipMemcpyAsync( M_dev, M_host.data(), sizeof(value_t) * M_host.nrows() * M_host.ncols(), hipMemcpyHostToDevice, handle.stream );
    
    if ( retval != hipSuccess )
    {
        hip_break();
        HWARNING( to_string( "(HIP) from_device : error in \"cublasSetMatrixAsync\" (%d)", retval ) );
        return false;
    }// if
    else
        return true;
}

//
// copy data from device
//
template < typename value_t >
bool
from_device ( handle_t                                   handle,
              typename hip_traits< value_t >::ptr_type  M_dev,
              BLAS::Matrix< value_t > &                  M_host )
{
    const auto  retval = hipMemcpyAsync( M_host.data(), M_dev, sizeof(value_t) * M_host.nrows() * M_host.ncols(), hipMemcpyDeviceToHost, handle.stream );

    if ( retval != hipSuccess )
    {
        hip_break();
        HWARNING( to_string( "(HIP) from_device : error in \"cublasGetMatrixAsync\" (%d)", retval ) );
        return false;
    }// if
    else
        return true;
}

template < typename value_t >
bool
from_device ( handle_t                                   handle,
              typename hip_traits< value_t >::ptr_type  v_dev,
              BLAS::Vector< value_t > &                  v_host )
{
    const auto  retval = hipMemcpyAsync( v_host.data(), v_dev, sizeof(value_t) * v_host.length(), hipMemcpyDeviceToHost, handle.stream );

    if ( retval != hipSuccess )
    {
        hip_break();
        HWARNING( to_string( "(HIP) from_device : error in \"cublasGetVectorAsync\" (%d)", retval ) );
        return false;
    }// if
    else
        return true;
}

template < typename value_t >
value_t
from_device ( handle_t                                   handle,
              typename hip_traits< value_t >::ptr_type  dev_data )
{
    value_t  data;

    const auto  retval = hipMemcpyAsync( & data, dev_data, sizeof(typename hip_traits< value_t >::hip_type), hipMemcpyDeviceToHost, handle.stream );

    if ( retval != hipSuccess )
    {
        hip_break();
        HWARNING( to_string( "(HIP) from_device : error in \"hipMemcpAsync\" (%d)", retval ) );
    }// if
    
    return data;
}

}// namespace anonymous

#endif

////////////////////////////////////////////////////////////////////////////////
//
// SVD related functions
//
////////////////////////////////////////////////////////////////////////////////

//
// compute SVD decomposition \f$ A = U·S·V^H \f$ of the nxm matrix \a A with
// n×min(n,m) matrix U, min(n,m)×min(n,m) matrix S (diagonal)
// and m×min(n,m) matrix V; \a A will be overwritten with U upon exit
//
#if HPRO_USE_HIP == 1
    
template < typename value_t >
bool
svd  ( const hip_handle_t                       handle_idx,
       BLAS::Matrix< value_t > &                 A,
       BLAS::Vector< real_type_t< value_t > > &  S,
       BLAS::Matrix< value_t > &                 V )
{
    if ( ! has_hip )
        return false;

    if ( handle_idx == INVALID_HANDLE )
        return false;

    // only support for nrows >= ncols
    if ( A.nrows() < A.ncols() )
        return false;

    using  real_t       = real_type_t< value_t >;
    using  hip_value_t = typename hip_traits< value_t >::hip_type;
    using  hip_real_t  = typename hip_traits< real_t >::hip_type;
    
    auto        handle = handles[ handle_idx ];
    const auto  nrows  = A.nrows();
    const auto  ncols  = A.ncols();
    const auto  minnm  = std::min( nrows, ncols );

    //
    // allocate data memory on device
    //

    const auto  type_A = hip_traits< value_t >::blas_type;
    const auto  type_S = hip_traits< real_t >::blas_type;
    
    auto  dev_A    = device_alloc< hip_value_t >( nrows * ncols );
    auto  dev_U    = device_alloc< hip_value_t >( nrows * nrows );
    auto  dev_S    = device_alloc< hip_real_t >( minnm );
    auto  dev_VT   = device_alloc< hip_value_t >( ncols * ncols );
    auto  dev_info = device_alloc< int >( 1 );

    if (! (( dev_A    != nullptr ) &&
           ( dev_U    != nullptr ) &&
           ( dev_S    != nullptr ) &&
           ( dev_VT   != nullptr ) &&
           ( dev_info != nullptr )) )
    {
        hip_break();
        HWARNING( "(HIP) svd : could not allocate memory on device" );
        
        device_free( dev_info );
        device_free( dev_VT );
        device_free( dev_S );
        device_free( dev_U );
        device_free( dev_A );

        return false;
    }// if
        
    to_device( handle, A, dev_A );

    hipsolverDnParams_t  params;

    HPRO_HIPSOLVER_CHECK( hipsolverDnCreateParams, ( & params ) ); // TODO: make it global?
    
    size_t      lwork_dev = 0;
    size_t      lwork_hst = 0;
    const char  jobU      = 'S';
    const char  jobV      = 'S';
    
    // double  pertubation_error;
    
    {
        // auto  retval = hipsolverDnXgesvdp_bufferSize( handle.solver, params,
        //                                              HIPSOLVER_EIG_MODE_VECTOR,
        //                                              1, // economy mode
        //                                              nrows, ncols,
        //                                              type_A, dev_A, nrows,
        //                                              type_S, dev_S,
        //                                              type_A, dev_U, nrows,
        //                                              type_A, dev_VT, minnm,
        //                                              type_A,
        //                                              & lwork_dev, & lwork_hst );
        
        auto  retval = hipsolverDnXgesvd_bufferSize( handle.solver, params,
                                                    jobU, jobV,
                                                    nrows, ncols,
                                                    type_A, dev_A, nrows,
                                                    type_S, dev_S,
                                                    type_A, dev_U, nrows,
                                                    type_A, dev_VT, minnm,
                                                    type_A,
                                                    & lwork_dev, & lwork_hst );

        if ( retval != HIPSOLVER_STATUS_SUCCESS )
        {
            hip_break();
            HWARNING( to_string( "(HIP) svd : error in \"hipsolverDnXgesvd*_bufferSize\" (code: %d)", retval ) );
            
            hipsolverDnDestroyParams( params );
            device_free( dev_info );
            device_free( dev_VT );
            device_free( dev_S );
            device_free( dev_U );
            device_free( dev_A );
            return false;
        }// if
    }

    auto  hst_work = std::vector< value_t >( lwork_hst );
    auto  dev_work = device_alloc< value_t >( lwork_dev );
        
    {
        // auto  retval = hipsolverDnXgesvdp( handle.solver, params,
        //                                   HIPSOLVER_EIG_MODE_VECTOR,
        //                                   1, // economy mode
        //                                   nrows, ncols,
        //                                   type_A, dev_A, nrows,
        //                                   type_S, dev_S,
        //                                   type_A, dev_U, nrows,
        //                                   type_A, dev_VT, minnm,
        //                                   type_A,
        //                                   dev_work, lwork_dev,
        //                                   hst_work.data(), lwork_hst,
        //                                   dev_info, & pertubation_error );

        auto  retval = hipsolverDnXgesvd( handle.solver,
                                         params,
                                         jobU, jobV,
                                         nrows, ncols,
                                         type_A, dev_A, nrows,
                                         type_S, dev_S,
                                         type_A, dev_U, nrows,
                                         type_A, dev_VT, minnm,
                                         type_A,
                                         dev_work, lwork_dev,
                                         hst_work.data(), lwork_hst,
                                         dev_info );

        if ( retval != HIPSOLVER_STATUS_SUCCESS )
        {
            hip_break();
            HWARNING( to_string( "(HIP) svd : error in \"hipsolverDnXgesvd*\" (code: %d)", retval ) );
            
            hipsolverDnDestroyParams( params );
            device_free( dev_info );
            device_free( dev_VT );
            device_free( dev_S );
            device_free( dev_U );
            device_free( dev_A );
            return false;
        }// if
    }

    // if ( hipStreamSynchronize( handle.stream ) != hipSuccess )
    //     HWARNING( "(HIP) svd : error in cuStreamSynchronize %d" );
        
    auto  info = from_device< int >( handle, dev_info );

    if ( info != 0 )
    {
        hip_break();
        if ( info < 0 ) { HWARNING( "(HIP) svd : " + to_string( "error in argument %d", info ) ); }
        else            { HWARNING( "(HIP) svd : " + to_string( "no convergence during SVD: %d", info ) ); }

        device_free( dev_work );
        hipsolverDnDestroyParams( params );
        device_free( dev_info );
        device_free( dev_VT );
        device_free( dev_S );
        device_free( dev_U );
        device_free( dev_A );
        return false;
    }// if
    
    //
    // get data from device
    //

    bool  retval = true;
    auto  VT     = BLAS::Matrix< value_t >( minnm, ncols );
    
    if (( A.nrows() != nrows ) || ( A.ncols() != minnm ))
        A = std::move( BLAS::Matrix< value_t >( nrows, minnm ) );
    
    if ( S.length() != minnm )
        S = std::move( BLAS::Vector< real_t >( minnm ) );

    if ( ! from_device( handle, dev_U,  A  ) ) retval = false;
    if ( ! from_device( handle, dev_S,  S  ) ) retval = false;
    if ( ! from_device( handle, dev_VT, VT ) ) retval = false;

    if (( V.nrows() != ncols ) || ( V.ncols() != minnm ))
        V = std::move( BLAS::Matrix< value_t >( ncols, minnm ) );
    
    BLAS::copy( BLAS::adjoint( VT ), V );
    
    HPRO_HIPSOLVER_CHECK( hipsolverDnDestroyParams, ( params ) ); // TODO: make it global?
    device_free( dev_work );
    device_free( dev_info );
    device_free( dev_VT );
    device_free( dev_S );
    device_free( dev_U );
    device_free( dev_A );

    return retval;
}

#else

template < typename value_t >
bool
svd  ( const hip_handle_t                       ,
       BLAS::Matrix< value_t > &                 ,
       BLAS::Vector< real_type_t< value_t > > &  ,
       BLAS::Matrix< value_t > &                  )
{
    return false;
}

#endif

#define INST_SVD( T )                                           \
    template bool                                               \
    svd< T >  ( const hip_handle_t                 handle_idx, \
                BLAS::Matrix< T > &                 U,          \
                BLAS::Vector< real_type_t< T > > &  S,          \
                BLAS::Matrix< T > &                 V )

INST_SVD( float );
INST_SVD( double );
INST_SVD( std::complex< float > );
INST_SVD( std::complex< double > );

#undef INST_SVD

}}// namespace Hpro::HIP
