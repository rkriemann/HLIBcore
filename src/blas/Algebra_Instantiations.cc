//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// explicit instantiations
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

template void  fill_rand ( Matrix< float > &             M );
template void  fill_rand ( Matrix< double > &            M );
template void  fill_rand ( Matrix< Complex< float > > &  M );
template void  fill_rand ( Matrix< Complex< double > > & M );


template float  norm2 ( const Matrix< float > &             M );
template double norm2 ( const Matrix< double > &            M );
template float  norm2 ( const Matrix< Complex< float > > &  M );
template double norm2 ( const Matrix< Complex< double > > & M );


template float  cond ( const Matrix< float > &              M );
template double cond ( const Matrix< double > &             M );
template float  cond ( const Matrix< Complex< float > > &   M );
template double cond ( const Matrix< Complex< double > > &  M );


template void  invert ( Matrix< float > &             M );
template void  invert ( Matrix< double > &            M );
template void  invert ( Matrix< Complex< float > > &  M );
template void  invert ( Matrix< Complex< double > > & M );


template void  invert ( Matrix< float > &             M,
                        const tri_type_t              tri_type,
                        const diag_type_t             diag_type );
template void  invert ( Matrix< double > &            M,
                        const tri_type_t              tri_type,
                        const diag_type_t             diag_type );
template void  invert ( Matrix< Complex< float > > &  M,
                        const tri_type_t              tri_type,
                        const diag_type_t             diag_type );
template void  invert ( Matrix< Complex< double > > & M,
                        const tri_type_t              tri_type,
                        const diag_type_t             diag_type );


template void pseudo_invert< float >            ( Matrix< float > &            A,
                                                  const TTruncAcc &            acc );
template void pseudo_invert< double >           ( Matrix< double > &           A,
                                                  const TTruncAcc &            acc );
template void pseudo_invert< Complex<float> >   ( Matrix< Complex<float> > &   A,
                                                  const TTruncAcc &            acc );
template void pseudo_invert< Complex<double> >  ( Matrix< Complex<double> > &  A,
                                                  const TTruncAcc &            acc );

template void lu< Matrix< float > >                ( Matrix< float > &            A );
template void lu< Matrix< double > >               ( Matrix< double > &           A );
template void lu< Matrix< Complex<float> > >       ( Matrix< Complex<float> > &   A );
template void lu< Matrix< Complex<double> > >      ( Matrix< Complex<double> > &  A );

template void llt< Matrix< float > >               ( Matrix< float > &            A );
template void llt< Matrix< double > >              ( Matrix< double > &           A );
template void llt< Matrix< Complex<float> > >      ( Matrix< Complex<float> > &   A );
template void llt< Matrix< Complex<double> > >     ( Matrix< Complex<double> > &  A );

template void llh< Matrix< float > >               ( Matrix< float > &            A );
template void llh< Matrix< double > >              ( Matrix< double > &           A );
template void llh< Matrix< Complex<float> > >      ( Matrix< Complex<float> > &   A );
template void llh< Matrix< Complex<double> > >     ( Matrix< Complex<double> > &  A );

template void ldlt< Matrix< float > >              ( Matrix< float > &            A );
template void ldlt< Matrix< double > >             ( Matrix< double > &           A );
template void ldlt< Matrix< Complex<float> > >     ( Matrix< Complex<float> > &   A );
template void ldlt< Matrix< Complex<double> > >    ( Matrix< Complex<double> > &  A );

template void ldlh< Matrix< float > >              ( Matrix< float > &            A );
template void ldlh< Matrix< double > >             ( Matrix< double > &           A );
template void ldlh< Matrix< Complex<float> > >     ( Matrix< Complex<float> > &   A );
template void ldlh< Matrix< Complex<double> > >    ( Matrix< Complex<double> > &  A );

template void qr< Matrix< float > >                ( Matrix< float > &            A,
                                                     Matrix< float > &            R );
template void qr< Matrix< double > >               ( Matrix< double > &           A,
                                                     Matrix< double > &           R );
template void qr< Matrix< Complex<float> > >       ( Matrix< Complex<float> > &   A,
                                                     Matrix< Complex<float> > &   R );
template void qr< Matrix< Complex<double> > >      ( Matrix< Complex<double> > &  A,
                                                     Matrix< Complex<double> > &  R );

template void tsqr< Matrix< float > >              ( Matrix< float > &            A,
                                                     Matrix< float > &            R,
                                                     const size_t                 ntile );
template void tsqr< Matrix< double > >             ( Matrix< double > &           A,
                                                     Matrix< double > &           R,
                                                     const size_t                 ntile );
template void tsqr< Matrix< Complex<float> > >     ( Matrix< Complex<float> > &   A,
                                                     Matrix< Complex<float> > &   R,
                                                     const size_t                 ntile );
template void tsqr< Matrix< Complex<double> > >    ( Matrix< Complex<double> > &  A,
                                                     Matrix< Complex<double> > &  R,
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
template void svd_double< Complex<double> > ( Matrix< Complex<double> > &  U,
                                              Vector< double > &           S,
                                              Matrix< Complex<double> > &  V );


#define INST_SVD( T )                                                   \
    template void                                                       \
    svd< Matrix< T > >  ( Matrix< T > &                                U, \
                          Vector< typename real_type< T >::type_t > &  S, \
                          Matrix< T > &                                V )
INST_SVD( float );
INST_SVD( double );
INST_SVD( Complex< float > );
INST_SVD( Complex< double > );
#undef INST_SVD


#define INST_SVD( T )                                                   \
    template void                                                       \
    svd< Matrix< T > >  ( Matrix< T > &                                U, \
                          Vector< typename real_type< T >::type_t > &  S, \
                          const bool left )
INST_SVD( float );
INST_SVD( double );
INST_SVD( Complex< float > );
INST_SVD( Complex< double > );
#undef INST_SVD


#define INST_SV( T )                                                    \
    template void                                                       \
    sv< Matrix< T > >  ( Matrix< T > &                                U, \
                         Vector< typename real_type< T >::type_t > &  S )
INST_SV( float );
INST_SV( double );
INST_SV( Complex< float > );
INST_SV( Complex< double > );
#undef INST_SV


#define INST_SV( T )                                                    \
    template void                                                       \
    sv< Matrix< T >, Matrix< T > >  ( Matrix< T > &                                A, \
                                      Matrix< T > &                                B, \
                                      Vector< typename real_type< T >::type_t > &  S )
INST_SV( float );
INST_SV( double );
INST_SV( Complex< float > );
INST_SV( Complex< double > );
#undef INST_SV

#define INST_APPROX( T )                                            \
    template size_t                                                 \
    approx< Matrix< T > > ( Matrix< T > &      M,                   \
                            const TTruncAcc &  acc,                 \
                            Matrix< T > &      A,                   \
                            Matrix< T > &      B )
INST_APPROX( float );
INST_APPROX( double );
INST_APPROX( Complex< float > );
INST_APPROX( Complex< double > );
#undef INST_APPROX


#define INST_APPROX_RRQR( T )                                           \
    template size_t                                                     \
    approx_rrqr< Matrix< T > >  ( Matrix< T > &      M,                 \
                                  const TTruncAcc &  acc,               \
                                  Matrix< T > &      A,                 \
                                  Matrix< T > &      B )
INST_APPROX_RRQR( float );
INST_APPROX_RRQR( double );
INST_APPROX_RRQR( Complex< float > );
INST_APPROX_RRQR( Complex< double > );
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
INST_APPROX_RANDSVD( Complex< float > );
INST_APPROX_RANDSVD( Complex< double > );
#undef INST_APPROX_RANDSVD


#define INST_TRUNCATE( T )                      \
    template size_t                             \
    truncate< T > ( Matrix< T > &      A,       \
                    Matrix< T > &      B,       \
                    const TTruncAcc &  acc )
INST_TRUNCATE( float );
INST_TRUNCATE( double );
INST_TRUNCATE( Complex< float > );
INST_TRUNCATE( Complex< double > );
#undef INST_TRUNCATE


#define INST_TRUNCATE_RRQR( T )                     \
    template size_t                                 \
    truncate_rrqr< T > ( Matrix< T > &      A,      \
                         Matrix< T > &      B,      \
                         const TTruncAcc &  acc )
INST_TRUNCATE_RRQR( float );
INST_TRUNCATE_RRQR( double );
INST_TRUNCATE_RRQR( Complex< float > );
INST_TRUNCATE_RRQR( Complex< double > );
#undef INST_TRUNCATE_RRQR


#define INST_TRUNCATE_RAND( T )                     \
    template size_t                                 \
    truncate_rand< T > ( Matrix< T > &      A,      \
                         Matrix< T > &      B,      \
                         const TTruncAcc &  acc )
INST_TRUNCATE_RAND( float );
INST_TRUNCATE_RAND( double );
INST_TRUNCATE_RAND( Complex< float > );
INST_TRUNCATE_RAND( Complex< double > );
#undef INST_TRUNCATE_RAND


#define INST_FACTORISE_ORTHO( T )                         \
    template void                                         \
    factorise_ortho< Matrix< T > >  ( Matrix< T > &  A,   \
                                      Matrix< T > &  R )
INST_FACTORISE_ORTHO( float );
INST_FACTORISE_ORTHO( double );
INST_FACTORISE_ORTHO( Complex< float > );
INST_FACTORISE_ORTHO( Complex< double > );
#undef INST_FACTORISE_ORTHO


#define INST_FACTORISE_ORTHO( T )                       \
    template void                                       \
    factorise_ortho< T >  ( Matrix< T > &      A,       \
                            Matrix< T > &      R,       \
                            const TTruncAcc &  acc )
INST_FACTORISE_ORTHO( float );
INST_FACTORISE_ORTHO( double );
INST_FACTORISE_ORTHO( Complex< float > );
INST_FACTORISE_ORTHO( Complex< double > );
#undef INST_FACTORISE_ORTHO


#define INST_TRUNCATE2_SVD( T )                     \
    template                                        \
    std::pair< Matrix< T >, Matrix< T > >           \
    truncate2_svd< T >  ( const Matrix< T > &  A,   \
                          const Matrix< T > &  B,   \
                          const TTruncAcc &    acc )
INST_TRUNCATE2_SVD( float );
INST_TRUNCATE2_SVD( double );
INST_TRUNCATE2_SVD( Complex< float > );
INST_TRUNCATE2_SVD( Complex< double > );
#undef INST_TRUNCATE2_SVD


