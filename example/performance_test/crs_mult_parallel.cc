// compile: g++ -O3 -DNDEBUG -std=c++11 array_test.cc
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <time.h>
#include <assert.h>

#include "mpi.h"

#include "../grid/mesh.hh"

using namespace conformingsimplexgrid;

//Sparse Matrix Class, stored in CRS format. Initialization is done via Triplet (I,J,X) (COO format).
template< class IntVec, class DoubleVec >
class CRSMatrix {

public:    
    CRSMatrix( IntVec &I, IntVec &J, DoubleVec &X, int &N) : N_(N) { 
        convert_ijx_to_crs(I, J, X, C_, R_, S_, N_);
        getDiagonalEntries(I, J, X, diagonalValues, N_);
    }
    
    DoubleVec operator * ( const DoubleVec &vec ) const
    {
        DoubleVec res(vec);
        res.fill(0.);
        for (int i = 0; i < N_; i++)
            for (int k = R_[i]; k < R_[i+1]; k++)
                res[i] += S_[k] * vec[C_[k]];
        return res;
    }

    DoubleVec precond( const DoubleVec &vec ) const
    {
        DoubleVec res(vec.size(),0.);
        for (int i = 0; i < diagonalValues.size(); i++)
            res[i] = vec[i] / diagonalValues[i];
        return res;
    }
    
private:
    
    void convert_ijx_to_crs(IntVec &I, IntVec &J, DoubleVec &X,
                        IntVec &C, IntVec &R, DoubleVec &S, int &N){
        C.resize(N*3); R.resize(N+1); S.resize(N*3); 
        int cumsum = 0, tmp = 0, last = 0;
        for (int k=0; k< (int)I.size(); k++) 
            R[I[k]]++;
        for (int i=0; i<N+1; i++){
            tmp = R[i]; R[i] = cumsum; cumsum += tmp;}
        for (int k=0; k< (int)I.size(); k++){
            C[R[I[k]]] = J[k]; S[R[I[k]]] = X[k]; R[I[k]]++;}
        for (int i=0; i<=N; i++){
            tmp = R[i]; R[i] = last; last = tmp;}
    }
    
    void getDiagonalEntries(IntVec &I, IntVec &J, DoubleVec &X, 
                            DoubleVec &D, int &N){
        D.resize(N); D.fill(0.);
        for (int k = 0; k < (int)I.size(); k++) 
            if( I[k] == J[k] )
                D[I[k]] += X[k];
    }

    int N_;
    IntVec C_, R_;
    DoubleVec S_;
    DoubleVec diagonalValues;
};

//CG Algorithm.
template< class SystemMatrix, class Vector >
size_t cg( SystemMatrix &A, Vector &rhs, Vector &x, double TOL, int iter_max, bool preconditioning ){
    double alpha, sp_r_z_new, sp_r_z_old, d_A_d;
    size_t k;
    Vector d(x.size()), z(x.size()), A_d(x.size()), r(x.size());
          
    r = rhs - A * x;
    if(preconditioning)
        z = A.precond(r);
    else
        z = r;
    d = z;
    sp_r_z_old = scalar_product( r , z );
    
    for( k = 0; k < iter_max; k++){
        A_d = A * d;
        d_A_d = scalar_product( d, A_d );
        alpha = sp_r_z_old / d_A_d;       
        x += alpha * d;
        r -= alpha * A_d;
        if(preconditioning)
            z = A.precond(r);
        else
            z = r;
        sp_r_z_new = scalar_product( r, z ); 
        if( sp_r_z_new < TOL )
            break;
        d = z + (sp_r_z_new/sp_r_z_old)*d;
        sp_r_z_old = sp_r_z_new;
    }
    printf("\nCG ALGORIHTM: sp_r = %lf, iter = %d, preconditioning = %d\n", sp_r_z_new, (int)(k),(int)preconditioning );
    return k;
}

template< class IntVec, class DoubleVec >
void fdm_poisson(IntVec &I, IntVec &J, DoubleVec &X, int N) {
  int ctr = 0;
  I[ctr] = 0; J[ctr] = 0; X[ctr] =  2.; ctr++;
  I[ctr] = 0; J[ctr] = 1; X[ctr] = -1.; ctr++;
  for (int k=1; k<N-1; k++){
    I[ctr] = k; J[ctr] = k-1; X[ctr] = -1.; ctr++;
    I[ctr] = k; J[ctr] = k+0; X[ctr] =  2.; ctr++;
    I[ctr] = k; J[ctr] = k+1; X[ctr] = -1.; ctr++;
  }
  I[ctr] = N-1; J[ctr] = N-2; X[ctr] = -1.; ctr++;
  I[ctr] = N-1; J[ctr] = N-1; X[ctr] =  2.; ctr++;
}

int main(int argc, char* argv[]){ 
  int N = 1e03;
  ManagedArray<double> u(N,0.), Au(N,0.), b(N,1.), X(3*N-2,0.);
  ManagedArray<int> I(3*N-2,0), J(3*N-2,0);
  fdm_poisson(I,J,X,N);
  CRSMatrix<ManagedArray<int>, ManagedArray<double>> A(I,J,X,N);
  cg(A,b,u,1e-3,10000,0);
  Au = A * u;
  double err = 0.;
  for (int j = 0; j < N; j++) err += pow((Au[j] - b[j])/(N+1),2.);
  printf("||Au-b||_L2 = %f\n",pow(err,.5));
}

