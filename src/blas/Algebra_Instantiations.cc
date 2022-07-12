//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// explicit instantiations
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

template void  fill_rand ( Matrix< float > &             M );
template void  fill_rand ( Matrix< double > &            M );
template void  fill_rand ( Matrix< std::complex< float > > &  M );
template void  fill_rand ( Matrix< std::complex< double > > & M );


template float  norm2 ( const Matrix< float > &             M );
template double norm2 ( const Matrix< double > &            M );
template float  norm2 ( const Matrix< std::complex< float > > &  M );
template double norm2 ( const Matrix< std::complex< double > > & M );


template float  cond ( const Matrix< float > &              M );
template double cond ( const Matrix< double > &             M );
template float  cond ( const Matrix< std::complex< float > > &   M );
template double cond ( const Matrix< std::complex< double > > &  M );


template void  invert ( Matrix< float > &             M );
template void  invert ( Matrix< double > &            M );
template void  invert ( Matrix< std::complex< float > > &  M );
template void  invert ( Matrix< std::complex< double > > & M );


template void  invert ( Matrix< float > &             M,
                        const tri_type_t              tri_type,
                        const diag_type_t             diag_type );
template void  invert ( Matrix< double > &            M,
                        const tri_type_t              tri_type,
                        const diag_type_t             diag_type );
template void  invert ( Matrix< std::complex< float > > &  M,
                        const tri_type_t              tri_type,
                        const diag_type_t             diag_type );
template void  invert ( Matrix< std::complex< double > > & M,
                        const tri_type_t              tri_type,
                        const diag_type_t             diag_type );


template void pseudo_invert< float >            ( Matrix< float > &            A,
                                                  const TTruncAcc &            acc );
template void pseudo_invert< double >           ( Matrix< double > &           A,
                                                  const TTruncAcc &            acc );
template void pseudo_invert< std::complex<float> >   ( Matrix< std::complex<float> > &   A,
                                                  const TTruncAcc &            acc );
template void pseudo_invert< std::complex<double> >  ( Matrix< std::complex<double> > &  A,
                                                  const TTruncAcc &            acc );

template void lu< Matrix< float > >                ( Matrix< float > &            A );
template void lu< Matrix< double > >               ( Matrix< double > &           A );
template void lu< Matrix< std::complex<float> > >       ( Matrix< std::complex<float> > &   A );
template void lu< Matrix< std::complex<double> > >      ( Matrix< std::complex<double> > &  A );

template void llt< Matrix< float > >               ( Matrix< float > &            A );
template void llt< Matrix< double > >              ( Matrix< double > &           A );
template void llt< Matrix< std::complex<float> > >      ( Matrix< std::complex<float> > &   A );
template void llt< Matrix< std::complex<double> > >     ( Matrix< std::complex<double> > &  A );

template void llh< Matrix< float > >               ( Matrix< float > &            A );
template void llh< Matrix< double > >              ( Matrix< double > &           A );
template void llh< Matrix< std::complex<float> > >      ( Matrix< std::complex<float> > &   A );
template void llh< Matrix< std::complex<double> > >     ( Matrix< std::complex<double> > &  A );

template void ldlt< Matrix< float > >              ( Matrix< float > &            A );
template void ldlt< Matrix< double > >             ( Matrix< double > &           A );
template void ldlt< Matrix< std::complex<float> > >     ( Matrix< std::complex<float> > &   A );
template void ldlt< Matrix< std::complex<double> > >    ( Matrix< std::complex<double> > &  A );

template void ldlh< Matrix< float > >              ( Matrix< float > &            A );
template void ldlh< Matrix< double > >             ( Matrix< double > &           A );
template void ldlh< Matrix< std::complex<float> > >     ( Matrix< std::complex<float> > &   A );
template void ldlh< Matrix< std::complex<double> > >    ( Matrix< std::complex<double> > &  A );

template void qr< Matrix< float > >                ( Matrix< float > &            A,
                                                     Matrix< float > &            R );
template void qr< Matrix< double > >               ( Matrix< double > &           A,
                                                     Matrix< double > &           R );
template void qr< Matrix< std::complex<float> > >       ( Matrix< std::complex<float> > &   A,
                                                     Matrix< std::complex<float> > &   R );
template void qr< Matrix< std::complex<double> > >      ( Matrix< std::complex<double> > &  A,
                                                     Matrix< std::complex<double> > &  R );

template void tsqr< Matrix< float > >              ( Matrix< float > &            A,
                                                     Matrix< float > &            R,
                                                     const size_t                 ntile );
template void tsqr< Matrix< double > >             ( Matrix< double > &           A,
                                                     Matrix< double > &           R,
                                                     const size_t                 ntile );
template void tsqr< Matrix< std::complex<float> > >     ( Matrix< std::complex<float> > &   A,
                                                     Matrix< std::complex<float> > &   R,
                                                     const size_t                 ntile );
template void tsqr< Matrix< std::complex<double> > >    ( Matrix< std::complex<double> > &  A,
                                                     Matrix< std::complex<double> > &  R,
                                                     const size_t                 ntile );

template void eigen< Matrix< float > >   ( Matrix< float > &            M,
                                           Vector< float > &            eig_val,
                                           Matrix< float > &            eig_vec );
template void eigen< Matrix< double > >  ( Matrix< double > &           M,
                                           Vector< double > &           eig_val,
                                           Matrix< double > &           eig_vec );

template void eigen< Matrix< float > >   ( Matrix< float > &            M,
                                           const Range &                eig_range,
                                           Vector< float > &            eig_val,
                                           Matrix< float > &            eig_vec );
template void eigen< Matrix< double > >  ( Matrix< double > &           M,
                                           const Range &                eig_range,
                                           Vector< double > &           eig_val,
                                           Matrix< double > &           eig_vec );

template void eigen< Vector< float >,
                     Vector< float > >   ( Vector< float > &            diag,
                                           Vector< float > &            subdiag,
                                           Vector< float > &            eig_val,
                                           Matrix< float > &            eig_vec );
template void eigen< Vector< double >,
                     Vector< double > >  ( Vector< double > &           diag,
                                           Vector< double > &           subdiag,
                                           Vector< double > &           eig_val,
                                           Matrix< double > &           eig_vec );

template void svd_double< double >          ( Matrix< double > &           U,
                                              Vector< double > &           S,
                                              Matrix< double > &           V );
template void svd_double< std::complex<double> > ( Matrix< std::complex<double> > &  U,
                                              Vector< double > &           S,
                                              Matrix< std::complex<double> > &  V );


#define INST_SVD( T )                                                   \
    template void                                                       \
    svd< Matrix< T > >  ( Matrix< T > &                                U, \
                          Vector< typename real_type< T >::type_t > &  S, \
                          Matrix< T > &                                V )
INST_SVD( float );
INST_SVD( double );
INST_SVD( std::complex< float > );
INST_SVD( std::complex< double > );
#undef INST_SVD


#define INST_SVD( T )                                                   \
    template void                                                       \
    svd< Matrix< T > >  ( Matrix< T > &                                U, \
                          Vector< typename real_type< T >::type_t > &  S, \
                          const bool left )
INST_SVD( float );
INST_SVD( double );
INST_SVD( std::complex< float > );
INST_SVD( std::complex< double > );
#undef INST_SVD


#define INST_SV( T )                                                    \
    template void                                                       \
    sv< Matrix< T > >  ( Matrix< T > &                                U, \
                         Vector< typename real_type< T >::type_t > &  S )
INST_SV( float );
INST_SV( double );
INST_SV( std::complex< float > );
INST_SV( std::complex< double > );
#undef INST_SV


#define INST_SV( T )                                                    \
    template void                                                       \
    sv< Matrix< T >, Matrix< T > >  ( Matrix< T > &                                A, \
                                      Matrix< T > &                                B, \
                                      Vector< typename real_type< T >::type_t > &  S )
INST_SV( float );
INST_SV( double );
INST_SV( std::complex< float > );
INST_SV( std::complex< double > );
#undef INST_SV

#define INST_APPROX( T )                                            \
    template size_t                                                 \
    approx< Matrix< T > > ( Matrix< T > &      M,                   \
                            const TTruncAcc &  acc,                 \
                            Matrix< T > &      A,                   \
                            Matrix< T > &      B )
INST_APPROX( float );
INST_APPROX( double );
INST_APPROX( std::complex< float > );
INST_APPROX( std::complex< double > );
#undef INST_APPROX


#define INST_APPROX_RRQR( T )                                           \
    template size_t                                                     \
    approx_rrqr< Matrix< T > >  ( Matrix< T > &      M,                 \
                                  const TTruncAcc &  acc,               \
                                  Matrix< T > &      A,                 \
                                  Matrix< T > &      B )
INST_APPROX_RRQR( float );
INST_APPROX_RRQR( double );
INST_APPROX_RRQR( std::complex< float > );
INST_APPROX_RRQR( std::complex< double > );
#undef INST_APPROX_RRQR


#define INST_APPROX_RANDSVD( T )                                        \
    template size_t                                                     \
    approx_randsvd< Matrix< T > >  ( Matrix< T > &      M,              \
                                     const TTruncAcc &  acc,            \
                                     Matrix< T > &      A,              \
                                     Matrix< T > &      B,              \
                                     const uint         power_steps,    \
                                     const uint         oversampling )
INST_APPROX_RANDSVD( float );
INST_APPROX_RANDSVD( double );
INST_APPROX_RANDSVD( std::complex< float > );
INST_APPROX_RANDSVD( std::complex< double > );
#undef INST_APPROX_RANDSVD


#define INST_TRUNCATE( T )                      \
    template size_t                             \
    truncate< T > ( Matrix< T > &      A,       \
                    Matrix< T > &      B,       \
                    const TTruncAcc &  acc )
INST_TRUNCATE( float );
INST_TRUNCATE( double );
INST_TRUNCATE( std::complex< float > );
INST_TRUNCATE( std::complex< double > );
#undef INST_TRUNCATE


#define INST_TRUNCATE_RRQR( T )                     \
    template size_t                                 \
    truncate_rrqr< T > ( Matrix< T > &      A,      \
                         Matrix< T > &      B,      \
                         const TTruncAcc &  acc )
INST_TRUNCATE_RRQR( float );
INST_TRUNCATE_RRQR( double );
INST_TRUNCATE_RRQR( std::complex< float > );
INST_TRUNCATE_RRQR( std::complex< double > );
#undef INST_TRUNCATE_RRQR


#define INST_TRUNCATE_RAND( T )                     \
    template size_t                                 \
    truncate_rand< T > ( Matrix< T > &      A,      \
                         Matrix< T > &      B,      \
                         const TTruncAcc &  acc )
INST_TRUNCATE_RAND( float );
INST_TRUNCATE_RAND( double );
INST_TRUNCATE_RAND( std::complex< float > );
INST_TRUNCATE_RAND( std::complex< double > );
#undef INST_TRUNCATE_RAND


#define INST_FACTORISE_ORTHO( T )                         \
    template void                                         \
    factorise_ortho< Matrix< T > >  ( Matrix< T > &  A,   \
                                      Matrix< T > &  R )
INST_FACTORISE_ORTHO( float );
INST_FACTORISE_ORTHO( double );
INST_FACTORISE_ORTHO( std::complex< float > );
INST_FACTORISE_ORTHO( std::complex< double > );
#undef INST_FACTORISE_ORTHO


#define INST_FACTORISE_ORTHO( T )                       \
    template void                                       \
    factorise_ortho< T >  ( Matrix< T > &      A,       \
                            Matrix< T > &      R,       \
                            const TTruncAcc &  acc )
INST_FACTORISE_ORTHO( float );
INST_FACTORISE_ORTHO( double );
INST_FACTORISE_ORTHO( std::complex< float > );
INST_FACTORISE_ORTHO( std::complex< double > );
#undef INST_FACTORISE_ORTHO


#define INST_TRUNCATE2_SVD( T )                     \
    template                                        \
    std::pair< Matrix< T >, Matrix< T > >           \
    truncate2_svd< T >  ( const Matrix< T > &  A,   \
                          const Matrix< T > &  B,   \
                          const TTruncAcc &    acc )
INST_TRUNCATE2_SVD( float );
INST_TRUNCATE2_SVD( double );
INST_TRUNCATE2_SVD( std::complex< float > );
INST_TRUNCATE2_SVD( std::complex< double > );
#undef INST_TRUNCATE2_SVD


