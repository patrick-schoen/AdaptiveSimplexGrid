#ifndef AC_CG_HH
#define AC_CG_HH

#include "../../../grid/mesh.hh"
#include "../../../grid/parallel.hh"

using namespace conformingsimplexgrid;

template< class SystemMatrix, class Vector >
size_t pcg( SystemMatrix &A, Vector &r, Vector &x, const double TOL ){
    double alpha, sp_r_z_new, sp_r_z_old;
    size_t iter_max = 10000;
    size_t k;

    //ManagedArray<double,1> r( x.size() );
    ManagedArray<double,1> d( x.size() );
    ManagedArray<double,1> z( x.size() );
    ManagedArray<double,1> A_d( x.size() );

    //rhs.get(r);
    r -= A * x;
    z = A.precond(r);
    d = z;
    sp_r_z_old = scalar_product(r,z);

    for( k = 0; k < iter_max; k++){
        A_d = A * d;
        alpha = sp_r_z_old / scalar_product(d,A_d);
        x += alpha * d;
        r -= alpha * A_d;
            z = A.precond(r);
            sp_r_z_new = scalar_product(r,z);
            if( sp_r_z_new < TOL )
                break;
        d = z + (sp_r_z_new/sp_r_z_old)*d;
        sp_r_z_old = sp_r_z_new;
    }

	//printf("\nPCG ALGORIHTM: sp_r = %lf, iter = %d\n", sp_r_z_new, (int)(k+1) );
	return k+1;
}

template<  class Decomposition, class SystemMatrix, class Vector >
size_t parallel_pcg( const Decomposition &decomp, const SystemMatrix &A, const Vector &rhs, const Vector &x0, const double TOL, Vector &x ){
    
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
    
    z = A.precond(r);
    //z = r;
    
    d = z;
    sp_r_z_old = decomp.scalar_product( r , z );
    

    for( k = 0; k < iter_max; k++){
        A_d = A * d;
        
        decomp.sumNodeVector( A_d );
                  
        d_A_d = decomp.scalar_product( d, A_d );
        
        alpha = sp_r_z_old / d_A_d;
                
        x += alpha * d;
        r -= alpha * A_d;
        
        z = A.precond(r);
        //z = r;
        
        sp_r_z_new = decomp.scalar_product( r, z ); 

        if( sp_r_z_new < TOL )
            break;
        d = z + (sp_r_z_new/sp_r_z_old)*d;
        sp_r_z_old = sp_r_z_new;
    }
    //if(decomp.master() )
    //    printf("\nPCG ALGORIHTM: sp_r = %lf, iter = %d \n ", sp_r_z_new, (int)(k+1) );
    return k+1;
}

template<  class Decomposition, class Communicator, class SystemMatrix, class Vector >
size_t parallel_pcg_opt( Decomposition &decomp, Communicator &comm, SystemMatrix &A, Vector &r, Vector &x, const double TOL ){
    
    double alpha, sp_r_z_new, sp_r_z_old, d_A_d;
    size_t iter_max = 1000;
    size_t k;
    
    assert( x.size() );
    
    ManagedArray<double,1> d( x.size() );
    ManagedArray<double,1> z( x.size() );
    ManagedArray<double,1> A_d( x.size() );
 
    A.get_bdy_entries( decomp );
      
    r -= A.matrix_vector_bdy(x);  
    ManagedArray< ManagedArray<double> > buffer(comm.size());
    ManagedArray< ManagedArray<double> > recv_buffer(comm.size());
    ManagedArray< typename Communicator::Request > requests(comm.size());
    decomp.sumNodeVector_send( r, requests, buffer  );
    decomp.sumNodeVector_recv( r, requests, recv_buffer );
    r -= A.matrix_vector_inner(x);
    decomp.sumNodeVector_wait( r, requests, recv_buffer );
    
    //z = r;
    z = A.precond(r);
    d = z;
    sp_r_z_old = decomp.scalar_product( r , z );
    

    for( k = 0; k < iter_max; k++){
        //A_d = A * d;
        A_d = A.matrix_vector_bdy(d); 
        
        ///sumNodeVector
        ManagedArray< typename Communicator::Request > requests(comm.size());
        ManagedArray< ManagedArray<double> > buffer(comm.size());
        ManagedArray< ManagedArray<double> > recv_buffer(comm.size());
        decomp.sumNodeVector_send( A_d, requests, buffer  );
        decomp.sumNodeVector_recv( A_d, requests, recv_buffer );
        A_d += A.matrix_vector_inner(d); 
        decomp.sumNodeVector_wait( A_d, requests, recv_buffer );
        
        //decomp.sumNodeVector( A_d );
                  
        d_A_d = decomp.scalar_product( d, A_d );

        alpha = sp_r_z_old / d_A_d;
                
        x += alpha * d;
        r -= alpha * A_d;
        z = A.precond(r);
        //z = r;
        sp_r_z_new = decomp.scalar_product( r, z ); 
        if( sp_r_z_new < TOL )
            break;
        d = z + (sp_r_z_new/sp_r_z_old)*d;
        sp_r_z_old = sp_r_z_new;
    }
    //if(decomp.master() )
    //    printf("\nPCG ALGORIHTM: sp_r = %lf, iter = %d\n ", sp_r_z_new, (int)(k+1) );
    return k+1;
}


#endif // AC_CG_HH
