#include <stdio.h>
#include <stdlib.h>
#include "../grid/array/array_lib.hh"
using namespace conformingsimplexgrid;
main(int argc, char **argv){ 
    ManagedArray<int,3> A, B;
    fill_random(A,5,0,10);
    fill_random(B,7,0,5);
    array_conjoin(A,B);           //stacks the two array A and B over another
    
    ManagedArray<int> newIndices;
    fill_random(newIndices,A.size(),6,10);
    index2newIndex(A,newIndices); //sets new Indices for A according to newIndices
    A.print();
    
    A = sortEachRow(A);           //sort values in every row of A
    ManagedArray<int> iA,iB;
    A = sortByRow(A,iA,iB);       //sort rows in A via lexicographical comp.
    
    B = unique_rows(A,iA,iB);     //get unique rows in A
    B.print();
    
    ManagedArray<int,3> C, D;
    C = rearange(A,iA);           //regain A,B from index vectors iA,iB
    D = rearange(B,iB);
    if( B==C ) printf("B == C\n");
    if( A==D ) printf("A == D\n");
    
}