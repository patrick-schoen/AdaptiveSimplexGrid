#ifndef CG_POISSON_HH
#define CG_POISSON_HH

#include "../../../grid/mesh.hh"
#include "../../../grid/parallel.hh"

using namespace conformingsimplexgrid;

template<  class SystemMatrix, class Vector >
size_t pcg( SystemMatrix &A, Vector &rhs, Vector &x0, double TOL, Vector &x, bool usePrecond ){
    
    double alpha, sp_r_z_new, sp_r_z_old, d_A_d;
    size_t iter_max = 1000;
    size_t k;
   
    x.resize( x0.size() );
    
    ManagedArray<double,1> d(x.size()), z(x.size());
    ManagedArray<double,1> A_d( x.size() );
    ManagedArray<double,1> r(rhs);
      
    r -= A * x0;
    if(usePrecond)
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
        if(usePrecond)
            z = A.precond(r);
        else
            z = r;
        sp_r_z_new = scalar_product( r, z ); 
        if( sp_r_z_new < TOL )
            break;
        d = z + (sp_r_z_new/sp_r_z_old)*d;
        sp_r_z_old = sp_r_z_new;
    }
    printf("\nPCG ALGORIHTM: sp_r = %lf, iter = %d\n", sp_r_z_new, (int)(k+1) );
    return k+1;
}

template<  class Decomposition, class SystemMatrix, class Vector >
size_t parallel_pcg( Decomposition &decomp, SystemMatrix &A, Vector &rhs, Vector &x0, double TOL, Vector &x, bool usePrecond  ){
    
    double alpha, sp_r_z_new, sp_r_z_old, d_A_d;
    size_t iter_max = 1000;
    size_t k;
   
    x.resize( x0.size() );
    
    ManagedArray<double,1> d( x.size() );
    ManagedArray<double,1> z( x.size() );
    ManagedArray<double,1> A_d( x.size() );
    ManagedArray<double,1> r(rhs);
      
    r -= A * x0;
    decomp.sumNodeVector( r );
    if(usePrecond)
        z = A.precond(r);
    else
        z = r;
    d = z;
    sp_r_z_old = decomp.scalar_product( r , z );
    for( k = 0; k < iter_max; k++){
        A_d = A * d;
        decomp.sumNodeVector( A_d );     
        d_A_d = decomp.scalar_product( d, A_d );
        alpha = sp_r_z_old / d_A_d;     
        x += alpha * d;
        r -= alpha * A_d;
        if(usePrecond)
            z = A.precond(r);
        else
            z = r;
        sp_r_z_new = decomp.scalar_product( r, z ); 
        if( sp_r_z_new < TOL )
            break;
        d = z + (sp_r_z_new/sp_r_z_old)*d;
        sp_r_z_old = sp_r_z_new;
    }
    if(decomp.master() )
        printf("\nPCG ALGORIHTM: sp_r = %lf, iter = %d\n", sp_r_z_new, (int)(k+1) );
    return k+1;
}

#endif
