//
// Project     : HLib
// File        : TPermutation.cc
// Description : class for a permutation matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2020. All Rights Reserved.
//

#include "hpro/vector/TScalarVector.hh"

#include "hpro/cluster/TPermutation.hh"

namespace HLIB
{

///////////////////////////////////////////
//
// constructor and destructor
//

TPermutation::TPermutation ( const size_t  asize )
{
    resize( asize );
}

TPermutation::TPermutation ( const std::vector< idx_t > &  perm )
{
    *this = perm;
}

TPermutation::TPermutation ( const TPermutation & perm )
        : std::vector< idx_t >( perm )
{}

TPermutation::~TPermutation ()
{
}

///////////////////////////////////////////
//
// misc. methods for permutation handling
//

//
// permute given vectors with source data in \a x and destination \a y
//
void
TPermutation::permute ( const TVector *  x,
                        TVector *        y ) const
{
    if (( x == nullptr ) || ( y == nullptr ))
        HERROR( ERR_ARG, "(TPermutation) permute", "vector is NULL" );
    
    if ( ! ( IS_TYPE( x, TScalarVector ) && IS_TYPE( y, TScalarVector ) ) )
        HERROR( ERR_VEC_TYPE, "(TPermutation) permute", "only supporting TScalarVector" );

    if (( x->size() != size() ) || ( y->size() != size() ))
        HERROR( ERR_VEC_SIZE, "(TPermutation) permute", "" );

    if ( x->is_complex() != y->is_complex() )
        HERROR( ERR_REAL_CMPLX, "(TPermutation) permute", "" );
        
    const TScalarVector *  sx = cptrcast( x, TScalarVector );
    TScalarVector *        sy = ptrcast( y, TScalarVector );

    if ( x->is_complex() )
    {
        for ( idx_t  i = 0; i < idx_t(size()); ++i )
        {
            sy->blas_cvec()( permute( i ) ) = sx->blas_cvec()( i );
        }// for
    }// if
    else
    {
        for ( idx_t  i = 0; i < idx_t(size()); ++i )
        {
            sy->blas_rvec()( permute( i ) ) = sx->blas_rvec()( i );
        }// for
    }// else
}

//
// use inverse permutation to reorder given vector \a x and write result to \a y
//
void
TPermutation::permute_inv ( const TVector *  x,
                            TVector *        y ) const
{
    if (( x == nullptr ) || ( y == nullptr ))
        HERROR( ERR_ARG, "(TPermutation) permute_inv", "vector is NULL" );
    
    if ( ! ( IS_TYPE( x, TScalarVector ) && IS_TYPE( y, TScalarVector ) ) )
        HERROR( ERR_VEC_TYPE, "(TPermutation) permute_inv", "only supporting TScalarVector" );

    if (( x->size() != size() ) || ( y->size() != size() ))
        HERROR( ERR_VEC_SIZE, "(TPermutation) permute_inv", "" );

    if ( x->is_complex() != y->is_complex() )
        HERROR( ERR_REAL_CMPLX, "(TPermutation) permute_inv", "" );
        
    const TScalarVector *  sx = cptrcast( x, TScalarVector );
    TScalarVector *        sy = ptrcast( y, TScalarVector );

    if ( x->is_complex() )
    {
        for ( idx_t  i = 0; i < idx_t(size()); ++i )
        {
            sy->blas_cvec()( i ) = sx->blas_cvec()( permute( i ) );
        }// for
    }// if
    else
    {
        for ( idx_t  i = 0; i < idx_t(size()); ++i )
        {
            sy->blas_rvec()( i ) = sx->blas_rvec()( permute( i ) );
        }// for
    }// else
}

namespace
{

//
// sort entries of vector according to given permutation
//
template <typename T>
void
vector_sort ( BLAS::Vector< T > &  x,
              TPermutation &    perm,
              const idx_t       lb,
              const idx_t       ub )
{
    if ( lb >= ub ) return;

    if ( (ub - lb) < 20 )
    {
        //
        // apply insertion sort for small ranges
        //

        for ( idx_t  i = lb+1; i <= ub; i++ )
        {
            const idx_t  v   = perm[i];
            idx_t        j   = i-1;
            T            tmp = x(i);
            
            while (( j >= 0 ) && ( perm[j] > v ))
            {
                x(j+1)    = x(j);
                perm[j+1] = perm[j];
                j--;
            }// if

            x(j+1)    = tmp;
            perm[j+1] = v;
        }// for
    }// if
    else
    {
        //
        // apply quick sort for standard ranges
        //

        idx_t        i         = lb;
        idx_t        j         = ub;
        const idx_t  mid       = (lb + ub) / 2;
        idx_t        choice[3] = { perm[lb], perm[mid], perm[ub] };
        idx_t        pivot;

        // choose pivot (median-of-three)
        if ( choice[0] > choice[1] ) std::swap( choice[0], choice[1] );
        if ( choice[0] > choice[2] ) std::swap( choice[0], choice[2] );
        if ( choice[1] > choice[2] ) std::swap( choice[1], choice[2] );
        pivot = choice[1];

        // partition
        while ( i < j )
        {
            while ( perm[i] < pivot   ) i++;
            while ( pivot   < perm[j] ) j--;

            if ( i < j )
            {
                std::swap( x(i), x(j) );
                std::swap( perm[i], perm[j] );
            }// if
        }// while

        // recursion
        vector_sort( x, perm, lb, i-1 );
        vector_sort( x, perm, i+1, ub );
    }// else
}

}// namespace anonymous

//
// permute given vector inplace
//
void
TPermutation::permute ( TVector *  x ) const
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TPermutation) permute", "vector is NULL" );
    
    if ( ! IS_TYPE( x, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TPermutation) permute", "only supporting TScalarVector" );

    if ( x->size() != size() )
        HERROR( ERR_VEC_SIZE, "(TPermutation) permute", "" );

    TScalarVector *  sx = ptrcast( x, TScalarVector );

    #if 0
    if ( x->is_complex() )
    {
        TPermutation  perm = *this;
        
        vector_sort( sx->blas_cvec(), perm, 0, idx_t(x->size())-1 );
    }// if
    else
    {
        TPermutation  perm = *this;
        
        vector_sort( sx->blas_rvec(), perm, 0, idx_t(x->size())-1 );
    }// else
    #else
    if ( x->is_complex() )
    {
        TScalarVector  tmp( * sx );

        permute( sx, & tmp );
        BLAS::copy( tmp.blas_cvec(), sx->blas_cvec() );
    }// if
    else
    {
        TScalarVector  tmp( * sx );

        permute( sx, & tmp );
        BLAS::copy( tmp.blas_rvec(), sx->blas_rvec() );
    }// else
    #endif
}
    
//
// apply inverse permutation to given vector \x
//
void
TPermutation::permute_inv ( TVector *  x ) const
{
    if ( x == nullptr )
        HERROR( ERR_ARG, "(TPermutation) permute_inv", "vector is NULL" );
    
    if ( ! IS_TYPE( x, TScalarVector ) )
        HERROR( ERR_VEC_TYPE, "(TPermutation) permute_inv", "only supporting TScalarVector" );

    if ( x->size() != size() )
        HERROR( ERR_VEC_SIZE, "(TPermutation) permute_inv", "" );

    TScalarVector *  sx = ptrcast( x, TScalarVector );

    if ( x->is_complex() )
    {
        TScalarVector  tmp( * sx );

        permute_inv( sx, & tmp );
        BLAS::copy( tmp.blas_cvec(), sx->blas_cvec() );
    }// if
    else
    {
        TScalarVector  tmp( * sx );

        permute_inv( sx, & tmp );
        BLAS::copy( tmp.blas_rvec(), sx->blas_rvec() );
    }// else
}
    
//
// invert permutation
//
void
TPermutation::invert ()
{
    // TODO: check maximal index in permutation and init with corresponding size
	const idx_t   s = idx_t(size());
	TPermutation  P( *this );

	for ( idx_t  i = 0; i < s; i++ )
		(*this)[ P[i] ] = i;
}

//
// return inverse permutation
//
TPermutation
TPermutation::inverse ()
{
	const idx_t   s = idx_t(size());
	TPermutation  P( *this );

	for ( idx_t  i = 0; i < s; i++ )
		P[ (*this)[ i ] ] = i;

    return P;
}

/////////////////////////////////////////////////
//
// serialisation
//

void
TPermutation::read ( TByteStream & s )
{
    size_t  n;
    
    s.get( & n, sizeof(size_t) );

    resize( n );
    s.get( this->data(), n * sizeof(idx_t) );
}

void
TPermutation::write ( TByteStream & s ) const
{
    size_t  n = size();
    
    s.put( & n, sizeof(size_t) );
    s.put( (* const_cast< TPermutation * >( this )).data(), n * sizeof(idx_t) );
}

//
// returns size of object in bytestream
//
size_t
TPermutation::bs_size () const
{
    return sizeof(size_t) + size() * sizeof(idx_t);
}

/////////////////////////////////////////////////
//
// misc.
//

//
// return size in bytes used by this object
//
size_t
TPermutation::byte_size () const
{
    return sizeof( std::vector< idx_t > ) + size() * sizeof(idx_t);
}

}// namespace
