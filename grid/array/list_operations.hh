/**
@file list_operations.hh
@brief defines sorting, unique, reindexing operations for the ManagedArray class.
*/

#ifndef LIST_OPERATIONS_HH
#define LIST_OPERATIONS_HH

namespace conformingsimplexgrid {

template< class Array, class ... Arrays >
void array_conjoin( Array &a0, Array & a1, Arrays & ... ai )        ///glues the Array a_0, a_1, ... , a_n together, such that a0 = [a0; a1; ... ; an].
{
    const int N = Array::DIM;
    const std::size_t pos0 = a0.size()*N;
    const std::size_t pos1 = pos0 + a1.size()*N;
    a0.resize( a0.size() + a1.size() );
    for( std::size_t k = pos0; k < pos1; ++k )
        a0[k] = a1[k-pos0];
    array_conjoin(a0, ai ...);
}

template< class Array >
void array_conjoin( Array & ) {}
    
template< class Array, class BoolVector>
void remove_false( Array &a, BoolVector &idx)   ///removes all rows k for with idx[k] is false.
{
    const int N = Array::DIM; 
    std::size_t ctr = 0;
#if N==1
    for( std::size_t el = 0; el < idx.size(); ++el)
        if( idx[el] )
            a[ctr++] = a[el];
#else
    for( std::size_t el = 0; el < idx.size(); ++el)
        if( idx[el] ){
            for( std::size_t k = 0; k < N; ++k)
                a(ctr,k) = a(el,k);
            ctr++;
        }
#endif
    a.resize(ctr);
}
    
template< class Array, class BoolVector>
void remove_false( const Array &a, BoolVector &idx, Array &b)       ///removes all rows k for with idx[k] is false, output is stored in b and a remains untouched.
{
    const int N = Array::DIM; 
    std::size_t ctr = 0;
    b.resize(idx.size());
#if N==1
    for( std::size_t el = 0; el < idx.size(); ++el)
        if( idx[el] )
            b[ctr++] = a[el];
#else
    for( std::size_t el = 0; el < idx.size(); ++el)
        if( idx[el] ){
            for( std::size_t k = 0; k < N; ++k)
                b(ctr,k) = a(el,k);
            ctr++;
        }
#endif
    b.resize(ctr);
} 
    
    
template< class Array, class IndexVector>
void remove_negative( Array &a, IndexVector &idx)       ///removes all rows k for with idx[k] is negative.
{
    const int N = Array::DIM; 
    std::size_t ctr = 0;
#if N==1
    for( std::size_t el = 0; el < idx.size(); ++el)
        if( !(idx[el] < 0) )
            a[ctr++] = a[el];
#else
    for( std::size_t el = 0; el < idx.size(); ++el)
        if( !(idx[el] < 0) ){
            for( std::size_t k = 0; k < N; ++k)
                a(ctr,k) = a(el,k);
            ctr++;
        }
#endif
    a.resize(ctr);
} 
    
template< class Array, class IndexVector>
void remove_negative( const Array &a, IndexVector &idx, Array &b)       ///removes all rows k for with idx[k] is negative, output is stored in b and a remains untouched.
{
    const int N = Array::DIM;
    b.resize(idx.size());
    std::size_t ctr = 0;
#if N==1
    for( std::size_t el = 0; el < idx.size(); ++el)
        if( !(idx[el] < 0) )
            b[ctr++] = a[el];
#else
    for( std::size_t el = 0; el < idx.size(); ++el)
        if( !(idx[el] < 0) ){
            for( std::size_t k = 0; k < N; ++k)
                b(ctr,k) = a(el,k);
            ctr++;
        }
#endif
    b.resize(ctr);
} 
    
template< class Array, class IndexVector >
void index2newIndex( Array &a, IndexVector &index2newIndex )      ///Sets new Indices for an array of indices A, such that A(i,j) = index2newIndex[A(i,j)].
{
    const int N = Array::DIM;
    for( std::size_t i = 0; i < a.size(); ++i)
        for( std::size_t j = 0; j < N; ++j){
            assert( index2newIndex[ a(i,j) ] >= 0);
            a(i,j) = index2newIndex[ a(i,j) ];
        }
}

template< class Array, class IndexVector >
Array rearange( Array &a, const IndexVector &idx )        ///Copys contents of A via index vector, such that B = A(idx,:). Returns array B.
{
    const int N = Array::DIM;
    Array b( idx.size() );
#if N==1
    for( std::size_t el = 0; el < idx.size(); ++el)
        b[el] = a[idx[el]];
#else
    for( std::size_t i = 0; i < idx.size(); ++i)
        for( std::size_t j = 0; j < N; ++j)
            b(i,j) = a(idx[i],j);
#endif         
    return b;
}

template< class Array, class IndexVector >
Array rearange_inverse( Array &a, const IndexVector &idx )        ///Copys contents of A into an array B, such that B(idx,:) = A. Returns array B.
{
    const int N = Array::DIM;
    Array b( a.size() );
#if N==1
    for( std::size_t el = 0; el < idx.size(); ++el)
        b[idx[el]] = a[el];
#else
    for( std::size_t i = 0; i < a.size(); ++i)
        for( std::size_t j = 0; j < N; ++j)
            b(idx[i],j) = a(i,j);
#endif 
    return b;
}

template< class Array, class IndexVector >
Array rearange_inverse_delete( Array &a, const IndexVector &idx ) ///Copys contents of A into an array B and delets content where idx < 0, such that B(idx(idx >= 0),:) = A(idx >= 0,:). Returns array B.
{
    const int N = Array::DIM;
    Array b( a.size() );
    std::size_t counter = 0;
    for( std::size_t i = 0; i < a.size(); ++i){
        if( idx[i] < 0 ) continue;
        for( std::size_t j = 0; j < N; ++j)
            b(idx[i],j) = a(i,j);
        counter++;
    }
    b.resize( counter );
    return b;
}
    
template< class IndexVector >
IndexVector invert_indexVector( IndexVector &idx )
{
    IndexVector idx_inv( idx.size() );
    for( std::size_t i = 0; i < idx.size(); ++i)
        idx_inv[ idx[i] ]= i;
    return idx_inv;
}
    
template< class IndexVector > 
std::size_t indexVector_enumerate( IndexVector &idx )   ///Populates vector with increasing indices, starting from 0.
{
    std::size_t ctr = 0;
    for( std::size_t i = 0; i < idx.size(); ++i)
        if( idx[i] >= 0 )
            idx[i] = ctr++;
    return ctr;
}

template< class Array >
void sortRows( Array &a ) ///Sort all rows of an array.
{         
    const int N = Array::DIM;
    for( size_t i = 0; i < a.size(); ++i )
        std::sort( a.begin() + i*N, a.begin() + (i+1)*N );
}

template< class Array >
Array sortEachRow( const Array &a ) ///Sort all rows of an array, returns result and leaves content unchanged.
{
    Array b(a);
    for( size_t i = 0; i < a.size(); ++i )
        std::sort( b.begin() + i*b.dim(), b.begin() + (i+1)*b.dim() );
    return b;
}
    
template< class Array, class BoolVector >
void sortRows_2d( Array &a, BoolVector &swap )    ///Sort all rows of an array with 2 columns, the value swap[i] is true iff initially a(i,0) > a(i,1).
{
    swap.resize( a.size() );
    for( size_t i = 0; i < a.size(); ++i )
        if( a( i, 0 ) > a( i, 1) ){
            std::swap( a( i, 0 ), a( i, 1 ) );
            swap[i] = true;
        } else {
            swap[i] = false;
        }
}
 
///@cond
template< class Array >
struct Array_Compare
{
    explicit Array_Compare ( const Array &a )
    : a_( a ) {}

    template< class Integer_Type >
    bool operator() ( Integer_Type j, Integer_Type k )
    { return std::lexicographical_compare( a_.data() + j*a_.dim(), a_.data() + (j+1)*a_.dim(), a_.data() + k*a_.dim(), a_.data() + (k+1)*a_.dim() ); }
    
    template< class Integer_Type >
    bool equal ( Integer_Type j, Integer_Type k )
    { return std::equal( a_.data() + j*a_.dim(), a_.data() + (j+1)*a_.dim(), a_.data() + k*a_.dim() ); }
    
    private:
    const Array &a_;
};
///@endcond
    
template< class Array, class IndexVector >
void sortByRow_void( const Array &a, IndexVector &ia ) ///sort rows of array via lexicographical comparison. The input array a is not changed, the results are given via the index array ia.
{
    Array_Compare<Array> comp( a );
    ia.resize( a.size() );
    for( std::size_t i = 0; i < a.size(); ++i ) 
        ia[i] = i;
    std::sort( ia.begin(), ia.end(), comp );
}
    
template< class Array, class IndexVector >
void sortByRow_void( const Array &a, IndexVector &ia, IndexVector &ib ) ///sort rows of array via lexicographical comparison. The input array a is not changed, the results are given via the index array ia. The reversed index array of ia is returned by ib.
{
    Array_Compare<Array> comp( a );
    ia.resize( a.size() );
    for( std::size_t i = 0; i < a.size(); ++i ) 
        ia[i] = i;
    std::sort( ia.begin(), ia.end(), comp );
    ib.resize( a.size() );
    for( std::size_t i = 0; i < a.size(); ++i )
        ib[ ia[i] ] = i;
}
    
template< class Array, class IndexVector >
Array sortByRow( const Array &a, IndexVector &ia ) ///sort rows of array via lexicographical comparison. The performed sorting is stored in ia.
{
    const int N = Array::DIM;
    Array b(a);
    sortByRow_void( a, ia );
    for( std::size_t i = 0; i < a.size(); ++i )
        for( std::size_t k = 0; k < N; ++k )
            b( i, k ) = a( ia[i], k );
    return b;
}

template< class Array, class IndexVector_1, class IndexVector_2 >
Array sortByRow( const Array &a, IndexVector_1 &ia, IndexVector_2 &ib )   ///sort rows of array via lexicographical comparison. The performed sorting is stored in ia. ib is the reversed index vector of ia.
{
    Array b = sortByRow( a, ia );
    ib.resize( a.size() );
    for( std::size_t i = 0; i < a.size(); ++i )
        ib[ ia[i] ] = i;
    return b;
}

template< class Vector  >
void vector_unique( Vector &v )
{   
    std::sort( v.begin(), v.end() );
    auto last = std::unique( v.begin(), v.end() );
    v.erase( last, v.end() );
}

template< class IndexVector >
void indexVector_unique( IndexVector &v )
{
    size_t m = *std::max_element(v.begin(), v.end());
    std::vector<bool> member(m, false);
    for(size_t k = 0; k < v.size(); ++k)
        member[v[k]] = true;
    std::size_t counter = 0;
    for(size_t k = 0; k < v.size(); ++k)
        if( member[k] ) 
            v[counter++] = v[k];
    v.resize(counter);
    for(size_t k = 0; k < v.size(); ++k)
        printf("%d ", v[k]);
    printf("\n");
}
     
#if 1
template< class Array, class IndexVector_1, class IndexVector_2 >
void unique_rows_void( const Array &a, IndexVector_1 &ia, IndexVector_2 &ib )
{
    std::vector<int> idx_sort( a.size() );
    ia.resize( a.size() );
    ib.resize( a.size() );
        
    sortByRow_void( a, idx_sort );
    Array_Compare<Array> comp( a );
    
    ia[0] = idx_sort[0]; ib[idx_sort[0]] = 0;
    std::size_t count = 0;
    for( std::size_t i = 1; i < a.size(); i++ ){
        if( !comp.equal( idx_sort[i-1], idx_sort[i] ) )
            ia[++count] = idx_sort[i];
        ib[idx_sort[i]] = count;
    }
    ia.resize( count+1 );    
}
#endif
     
#if 0
template< class Array, class IndexVector_1, class IndexVector_2 >
void unique_rows_void( const Array &a, IndexVector_1 &ia, IndexVector_2 &ib )
{
    Array b(a);
    std::vector<int> idx_sort( b.size() );
    ia.resize( b.size() );
    ib.resize( b.size() );
        
    //sortRows( b );
    sortByRow_void( b, idx_sort );
    Array_Compare<Array> comp( b );
    
    ia[0] = idx_sort[0]; ib[idx_sort[0]] = 0;
    std::size_t count = 0;
    for( std::size_t i = 1; i < b.size(); i++ ){
        if( !comp.equal( idx_sort[i-1], idx_sort[i] ) )
            ia[++count] = idx_sort[i];
        ib[idx_sort[i]] = count;
    }
    ia.resize( count+1 );    
}
#endif
    
template< class Array, class IndexVector_1, class IndexVector_2 >
Array unique_rows( const Array &a, IndexVector_1 &ia, IndexVector_2 &ib ) ///The returned array b contains the unique rows of a 
{
    const int N = Array::DIM;
    Array b;
    unique_rows_void( a, ia, ib );
    b.resize( ia.size() );
    for( std::size_t i = 0; i < b.size(); i++)
        for( std::size_t k = 0; k < N; k++)
            b( i, k ) = a( ia[i], k);
    return b;
}

template< class Array, class IndexVector_1, class IndexVector_2 >
void unique_rows_void_nosort( const Array &a, IndexVector_1 &ia, IndexVector_2 &ib )
{
    ManagedArray<int> idx_sort( a.size() );
    ia.resize( a.size() );
    ib.resize( a.size() );
        
    Array_Compare<Array> comp( a );
    for(std::size_t k = 0; k < idx_sort.size(); ++k ) idx_sort[k] = k;
    std::sort( idx_sort.begin(), idx_sort.end(), comp );
    
    ia[0] = idx_sort[0]; ib[idx_sort[0]] = 0;
    std::size_t count = 0;
    for( std::size_t i = 1; i < a.size(); i++ ){
        if( !comp.equal( idx_sort[i-1], idx_sort[i] ) )
            ia[++count] = idx_sort[i];
        ib[idx_sort[i]] = count;
    }

    ia.resize( count+1 ); 
    idx_sort.resize( ia.size() );
    ia = sortByRow( ia, idx_sort );
    idx_sort = invert_indexVector( idx_sort );
    index2newIndex( ib, idx_sort );     
}

template< class Array, class IndexVector_1, class IndexVector_2 >
Array unique_rows_nosort( const Array &a, IndexVector_1 &ia, IndexVector_2 &ib )
{
    const int N = Array::DIM;
    Array b;
    unique_rows_void_nosort( a, ia, ib );        
    b.resize( ia.size() );
    for( std::size_t i = 0; i < b.size(); i++)
        for( std::size_t k = 0; k < N; k++)
            b( i, k ) = a( ia[i], k);            
    return b;
}

template < class Array, class T >
void fill_random( Array &a, int n, T min, T max )
{
    srand (time(NULL));
    const int N = Array::DIM;
    a.resize( n );
    for( std::size_t i = 0; i < a.size(); ++i)
        for( std::size_t j = 0; j < N; ++j)
            a(i,j) = ((T)rand() / RAND_MAX ) * ((max+1) - min) + min;
}

template < class Array >
void fill_random( Array &a, int n, int min, int max )
{
    srand (time(NULL));
    const int N = Array::DIM;
    a.resize( n );
    for( std::size_t i = 0; i < a.size(); ++i)
        for( std::size_t j = 0; j < N; ++j)
            a(i,j) = (int)(rand()) % ((max+1) - min) + min;
}
    
} // namespace conformingsimplexgrid
#endif // LIST_OPERATIONS_HH
 
 
