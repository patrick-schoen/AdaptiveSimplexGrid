#ifndef ARRAY_TEST_HH
#define ARRAY_TEST_HH

#include "../simplexGrid.hh"

namespace conformingsimplexgrid {
    
    template < class Vector >
    void vector_print( Vector &A )
    {
        std::cout << "vector: " << A.size() << std::endl;
        for( std::size_t i = 0; i < A.size(); ++i){
            std::cout << i << " : "; 
                std::cout << A[i] << " ";
            std::cout << std::endl;
        }    
        std::cout << std::endl;
    }
    
    template < class Array >
    void array_print( Array &A )
    {
        std::cout << "array: " << A.size() << std::endl;
        for( std::size_t i = 0; i < A.size(); ++i){
            std::cout << i << " : "; 
            for( std::size_t j = 0; j < A.dim(); ++j){
                std::cout << A(i,j) << " ";
            }
            std::cout << std::endl;
        }    
        std::cout << std::endl;
    }
    
    template< class T, int k, int N >
    void test_sortRow() {
        std::vector< int > ia;
        std::vector< int > ia1;
        std::vector< int > ib;
        ManagedArray< T, k > A( 0 );
        ManagedArray< T, k > B( 0 );
        ManagedArray< T, k > C( 0 );
        ManagedArray< T, k > D( 0 );
        
        Array_Compare< ManagedArray< T, k > > comp_A( A );
        Array_Compare< ManagedArray< T, k > > comp_B( B );
        
        fill_random(A, N, (T)0, (T)5 );
        
        sortRows( A );
        
        for( size_t i = 0; i < A.size(); ++i )
            for( size_t j = 0; j < A.dim() - 1; ++j )
                assert( A(i,j) <= A(i,j+1) );
            
        B = sortByRow( A, ia );  
        B = sortByRow( A, ia );
        ib = invert_indexVector( ia);
        
        //assert( ia1 == ia );
      
        for( size_t i = 0; i < B.size() - 1; ++i )
            assert( comp_B(i,i+1) || comp_B.equal(i,i+1) ); 
            
 //       C = rearange( A, ia );
 //       D = rearange( B, ib );

  //      assert( C == B );
 //       assert( D == A );
        
   //     C = rearange_inverse( A, ib );
   //     D = rearange_inverse( B, ia );  
 //       
   //     assert( C == B );
  //      assert( D == A );
        
   //     ManagedArray<T,2> E( 0 );
    //    ManagedArray<T,2> F( 0 );
    //    std::vector<bool> swap;
    //    fill_random(E, N, (T)0, (T)5 );
        
    //    F = E;

    //    sortRows_2d( E, swap );
        
    //    for( size_t i = 0; i < E.size(); ++i ){
    //        assert( E(i,0) <= E(i,1) );
    //        if( swap[i] == true ){
    //            assert( F(i,0) >= E(i,0) );
    //        }
    //    }
    }    
    
    template< class T, int k, int N >
    void test_unique_row() {
        
        ManagedArray< int > ia;
        ManagedArray< int > ib;
        ManagedArray< T, k > A( 0 );
        ManagedArray< T, k > B( 0 );
        ManagedArray< T, k > C( 0 );
        ManagedArray< T, k > D( 0 );
        
        fill_random(A, N, (T)0, (T)5 );
        
        //sortRows(A);

        B = unique_rows( A, ia, ib );

        
        Array_Compare< ManagedArray< T, k > > comp_A( A );
        Array_Compare< ManagedArray< T, k > > comp_B( B );
  
        //for( size_t i = 0; i < B.size() - 1; ++i )
        //    assert( comp_B(i,i+1) ); 
        
        C = rearange( A, ia );
        D = rearange( B, ib );
        
        sortRows(A);
        sortRows(B);
        sortRows(C);
        sortRows(D);
       
        assert( C == B );
        assert( D == A );
 
    }
    
    
    template< int k, int N >
    bool test_unique_row2() {
        
        ManagedArray< int > ia(N);
        ManagedArray< int > ib(N);
        ManagedArray< int, k > A(N);
        ManagedArray< int, k > B(N);
 //       ManagedArray< T, k > C( 0 );
 //       ManagedArray< T, k > D( 0 );
        
        fill_random(A, N, 0, 3 );
        
        //sortRows(A);

        B = unique_rows_nosort( A, ia, ib );

        
//         Array_Compare< ManagedArray< T, k > > comp_A( A );
//         Array_Compare< ManagedArray< T, k > > comp_B( B );
//   
//         //for( size_t i = 0; i < B.size() - 1; ++i )
//         //    assert( comp_B(i,i+1) ); 
//         
//         C = rearange( A, ia );
//         D = rearange( B, ib );
//         
//         sortRows(A);
//         sortRows(B);
//         sortRows(C);
//         sortRows(D);
//        
//         assert( C == B );
//         assert( D == A );
        
        return true;
 
    }
    
    template< int k, int N >
    void test_index_function() {
        
        ManagedArray< int > ia;
        ManagedArray< int > ib;
        ManagedArray< int, k > A( 0 );
        ManagedArray< int, k > B( 0 );
        
        fill_random(A, N, (int)0, (int)7 );
        
        //B = sortByRow( A, ia, ib );
        B = sortByRow( A, ia );
        ib = invert_indexVector( ia );
        
        ib = invert_indexVector( ia );
        ib = invert_indexVector( ib );

        assert( ib == ia );
        
        B = A;
        index2newIndex( A, ia );
        ia = invert_indexVector( ia );
        index2newIndex( A, ia );
        assert( B == A );
        
        A = rearange( A, ia );
        A = rearange_inverse( A, ia );
        assert( B == A );
    }
    
    template<class T, int k, int N >
    void test_array_numerics() {
        ManagedArray< T, k > A( 0 );
        ManagedArray< T, k > B( 0 );
        ManagedArray< T, k > C( 0 );
        ManagedArray< T, k > D( 0 );

        fill_random(A, N, (int)0, (int)7 );
        fill_random(B, N, (int)0, (int)7 );
        
        C = A;
        C += A + A + A;
        D = 4*A;
        assert( C == D );
        
        D.fill(2);
        C = A;
        C *= D;
        D = 2*A;
        assert( C == D );
        
        A *= B;
        A += B;
        A -= B;
        A = A + B;
        A = A - B;
        A = A * B;
        
        A + 1;
        A - 1;
        A * 1;
        A / 1;
        A ^ 1;        
        
        1 + A + 1;
        1 - A - 1;
        1 * A * 1;
                
        A + 1.;
        A - 1.;
        A * 1.;
        A / 1.;
        A ^ 1.;
        
        B.fill(2);
        1 / B;
        A = A / B;
    }
    
    void test_array_lib()
    {
        #if 0
        std::cout << "testing array lib functionality ...";
    
        unsigned int seed = time(NULL);
        srand(seed);

        static const int N = 111;
        
        test_sortRow<double, 1, N>();
        test_sortRow<double, 5, N>();
        test_sortRow<double, 7, N>();
        test_sortRow<int, 1, N>();
        test_sortRow<int, 4, N>();
        test_sortRow<int, 17, N>();
        
        test_index_function<1, N>();
        test_index_function<2, N>();
        test_index_function<3, N>();
      
        test_unique_row<int, 1, N>();
        test_unique_row<int, 4, N>();
        test_unique_row<int, 18, N>();
        
        test_array_numerics<int, 1, 5>();
        test_array_numerics<double, 1, 6>();
        test_array_numerics<int, 4, 4>();
        test_array_numerics<double, 4, 3>();
        
        std::cout << " done.\n";
        #endif
    }

} // namespace conformingsimplexgrid
#endif // ARRAY_TEST_HH

