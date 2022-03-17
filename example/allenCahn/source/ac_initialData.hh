#ifndef AC_INITIALDATA_HH
#define AC_INITIALDATA_HH

#include<math.h>
#include "../../../grid/mesh.hh"

using namespace conformingsimplexgrid;

template<int DIM>
ManagedArray<double,1> u_0( ManagedArray<double,DIM> &x, double eps)
{
    int nC = x.size();
    ManagedArray<double,1> val(nC);
    for(int i = 0; i < nC; i++){
        for(int k = 0; k < DIM; k++)
            val[i] += x[i*DIM + k]*x[i*DIM + k];
        val[i] = tanh( (sqrt(val[i]) - 1))/(2*sqrt(2) );
    }
    return val;
}

double dist1( double *x ){
    double val;
    double alpha[4] = {-4.1600,  3.1200,  0.,  0.1200};
    val = sqrt(x[1]*x[1] + x[2]*x[2]) - (alpha[0]*(fabs(x[0])*x[0]*x[0]) + alpha[1]*x[0]*x[0] + alpha[3]);
    return val;
}

double dist2( double *x, double* m){
    double radius = .38;
    double val;
    val = sqrt( (x[0]-m[0])*(x[0]-m[0]) + (x[1]-m[1])*(x[1]-m[1]) + (x[2]-m[2])*(x[2]-m[2]) ) - radius;
    return val;
}  

template<int DIM>
ManagedArray<double,1> u_0_dumbbell(ManagedArray<double,DIM> &x, double eps){
    double *z;
    int nC = x.size();
    ManagedArray<double,1> val(nC);
    z = (double*)std::malloc( 3 * sizeof(double));
    double delta = 0.12; double radius = 0.38;
    double m_1[3] = {delta + radius, 0., 0.};
    double m_2[3] = {-delta - radius, 0., 0.};
    for(int i = 0; i < nC; i++){
        val[i] = 0.;
        z[0] = x[i*3]; z[1] = x[i*3+1]; z[2] = x[i*3+2];
        if( x[i*3] < m_1[0]  && x[i*3] > m_2[0]  )
            val[i] = tanh( dist1( z )/eps );
        else if( x[i*3] >= m_1[0] )
            val[i] = tanh( dist2( z, m_1) /eps);
        else if( x[i*3] <= m_2[0] )
            val[i] = tanh( dist2( z, m_2) /eps);
    }
    free(z);
    return val;	
}

double dist_circle( std::array<double,3> &x, std::array<double,3> &m, double &radius){
    return sqrt( (x[0]-m[0])*(x[0]-m[0]) + (x[1]-m[1])*(x[1]-m[1]) + (x[2]-m[2])*(x[2]-m[2]) ) - radius;
}

template<int DIM>
ManagedArray<double,1> u_0_circle(ManagedArray<double,DIM> &x, std::array<double,3> &m, double &radius, double &eps){
    int nC = x.size();
    ManagedArray<double,1> val(nC);
    std::array<double,3> z;
    for(int i = 0; i < nC; i++){
        z[0] = x(i,0); z[1] = x(i,1); z[2] = x(i,2);
        val[i] = tanh( dist_circle(z,m,radius) / eps );
    }
    return val;
}

#endif // AC_INITIALDATA_HH


