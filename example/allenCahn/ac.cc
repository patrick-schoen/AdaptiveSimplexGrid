#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <string>
#include <limits>
#include <time.h>
#include <cmath>
        
#include "../../grid/mesh.hh"
#include "../../grid/parallel.hh"

#include "source/ac_dataContainer.hh"
#include "source/ac_initialData.hh"
#include "source/ac_pdeSystem.hh"
#include "source/ac_cg.hh"
#include "source/doerfler.hh"

#include <math.h>

using namespace conformingsimplexgrid;
using Communicator = MpiCommunicator;
using Grid = Grid_3D;
using Decomposition = Decomposition_3D<Grid,Communicator>;
using GlobalId = GlobalNodeId<Grid, Communicator, Decomposition>;

template<class Matrix, class MassMatrix, class Vector, class Decomposition>
double vector_iteration(Matrix &Y, MassMatrix &M, Vector &w, Decomposition & subDecomposition ){
    double mu = 0.; 
    double mu_old = 0.;
    double diff_mu = 1.0; 
    double eps_stop = .1;
    double TOL_CG = 1e-7;
        
    ManagedArray<double> rhs_val(w.size(),0.);
    ManagedArray<double> v(w.size(),0.);
    ManagedArray<double> w_start(w.size(),0.);
    while (fabs(diff_mu) > eps_stop) {
        rhs_val = M*w;
        subDecomposition.sumNodeVector(rhs_val);
        w_start = w;
        parallel_pcg( subDecomposition, Y, rhs_val, w_start, TOL_CG, w );
        v = M*w;
        subDecomposition.sumNodeVector(v);
        double tmp = std::sqrt(subDecomposition.scalar_product(w,v));
        for(size_t k = 0; k < w.size(); k++) w[k] = w[k]/tmp;
        v = Y*w;
        subDecomposition.sumNodeVector(v);
        mu = subDecomposition.scalar_product(w,v);
        diff_mu = mu - mu_old;
        mu_old = mu;
        if(subDecomposition.master())
            printf("rank = %d, mu = %lf tmp = %lf,diff_mu = %lf\n", subDecomposition.rank(), mu, tmp, diff_mu);
    } 
    return mu;
}


int main(int argc, char **argv)
{
    Communicator comm(&argc, &argv);
    
    Grid subgrid;
    Decomposition subDecomposition( subgrid, comm );
    GlobalId subGlobalId( subgrid, comm, subDecomposition );
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    double time1, start_time;
    
    if(comm.rank() == 0){
        Grid grid;
        importFromFile( grid, "3dmesh.txt" );
        //test_grid(grid);
        
        start_time = clock();
        for(int k = 0; k < 12; ++k){
            ManagedArray<int> marked( grid.nElements(), 1);
            refine(grid,marked);
        }
           
        ManagedArray<int,1> partition;    
        metisDecomposition(grid, partition, (int)comm.size());
        
        Decomposition decomposition( grid, comm );
        GlobalId globalId( grid, comm, decomposition );
        globalId.init_macro( partition );
        decomposition.init_partition( partition );
   
        for( int rank = 0; rank < (int)comm.size(); rank++){
            if( comm.rank() == rank) 
                continue;
            getSubGrid( subgrid, subDecomposition, subGlobalId, rank, grid, decomposition, globalId, partition );
            auto sentMesh = sendMesh( subgrid, rank, 1, comm, subGlobalId, subDecomposition );
        }
        getSubGrid( subgrid, subDecomposition, subGlobalId, comm.rank(), grid, decomposition, globalId, partition );
        getFacesForElements( subgrid.elements, subgrid.faces, subgrid.element2faces, subgrid.level );
        getFacesForElements( subgrid.faces, subgrid.edges, subgrid.face2edges );  
    } else {
        recvMesh( 0, 1, comm, subgrid, subGlobalId, subDecomposition);
    }
    
    subDecomposition.init_communication();
    
    //get geoemtry data
    ManagedArray<double,1> mass;
    ManagedArray<double,NV*DIM> normals; 
    outerNormals( mass, normals, subgrid.elements, subgrid.coordinates );
    
    //set solution vectors
    ManagedArray<double> u(subgrid.nNodes());
    ManagedArray<double> u_old(subgrid.nNodes());
    ManagedArray<double> rhs_val(subgrid.nNodes());
    
    //set pde parameter
    double T = .1;
    double eps = 1./16.;
    double tau = pow(eps,2)/10;
    double TOL_CG = 1e-9;
    double TOL = 1.;
    
    //set indicator variables
    ManagedArray<int> marked( subgrid.nElements(), 1);
    local_norm<Grid, Decomposition> norm_T( subgrid, subDecomposition, normals, mass );
    ManagedArray<double, 1> eta_space( subgrid.elements.size() );
    ManagedArray<double, 1> eta_time( subgrid.elements.size() );
    ManagedArray<double, 1> eta_coarse( subgrid.elements.size() );
    ManagedArray<double, 1> h_T( subgrid.elements.size() );

    //get system matrix and right hand side
    SystemMatrix_AC<NV, DIM, Decomposition > A( subDecomposition, subgrid.elements, normals, mass, u, u_old, eps, tau );
    rhs_AC<NV,DIM> rhs( subgrid.elements, normals, mass, u_old, eps, tau );

    //set up manager to adapt data
    AC_Manager2<Communicator,Grid> ac_manager( comm, subgrid, u, u_old );
    
    //set initial state
    u = u_0_dumbbell<DIM>( subgrid.coordinates, eps);
    u_old = u;
    
    Timer<Decomposition> timer(subDecomposition);
    
    int nRefinements = 3;
    int nIterations_total = 50;
      
    if(argc > 1) nRefinements = atoi(argv[2]);
    if(argc > 2) nIterations_total = atoi(argv[3]);
    
    for( int rr = 0; rr < nRefinements; rr++ ) {
        ManagedArray<int> marked( subgrid.nElements(), 1);
        refine_parallel(subDecomposition,subgrid,marked,ac_manager);
    }
    
    ManagedArray<int> partition;
    parMetisDecomposition(comm, subDecomposition, subgrid, partition);
    
    subGlobalId.communicate();
    redistributeNewPartition( subgrid, comm, subDecomposition, subGlobalId, partition, ac_manager );
   
    //recalculate geoemtry data
    outerNormals( mass, normals, subgrid.elements, subgrid.coordinates );
    
    u = u_0_dumbbell<DIM>( subgrid.coordinates, eps);
    u_old = u;

    int iter;
    std::vector<double> lambda_values;
    double lambda = 0.;
    double lambda_neg = 0.;
    
    
    //ManagedArray<double> f_u; f_u = f(u); 
    //subDecomposition.sumNodeVector(f_u);
#if 0
    double c_shift = 0.;
    //subDecomposition.scalar_max(c_shift);
    ManagedArray<double> w(subgrid.nNodes(),0.);
    ManagedArray<double> w_init(subgrid.nNodes(),0.);
    double min = -.5; double max = .5;
    srand (time(NULL));
    for( std::size_t i = 0; i < w.size(); ++i )
        w[i] = ((double)rand() / RAND_MAX ) * ((max+1) - min) + min;
    subDecomposition.sumNodeVector(w);
    for( std::size_t i = 0; i < w.size(); ++i){
        w[i] = w[i]/(double)subDecomposition.node2ranks[i].size();
    }
    w_init = w;
    
    EigenvalueMatrix_AC<NV, DIM, Decomposition > Y ( subDecomposition, subgrid.elements, normals, mass, u, u, eps, c_shift );
    MassMatrix<NV, DIM, Decomposition > M ( subDecomposition, subgrid.elements, mass );
#endif

    //for( int j = 0; j < (T/tau); j++ ) { 
    timer.begin();
    for( int j = 0; j < nIterations_total; j++ ) {  
       
        rhs.get( rhs_val );
        ManagedArray<double> xstart(u);
        parallel_pcg( subDecomposition, A, rhs_val, xstart, TOL_CG, u );
        
        u_old = u;
        //estimate
        ManagedArray<double> h_T_2 = norm_T.h_T_2();
        ManagedArray<double> f_u = f(u);
        ManagedArray<double> xi = (1/tau)*(u-u_old) + (1/eps/eps)*f_u;
        
        ManagedArray<double> norm_l2  = norm_T.l2( xi );
        eta_space.resize( subgrid.nElements() );
        eta_space.fill(0.);
        for( std::size_t i = 0; i < eta_space.size(); i++)
            eta_space[i] = h_T_2[i]*norm_l2[i];
        
        eta_space += norm_T.h1_jump( u );

        double eta = eta_space.sum();
        subDecomposition.scalar_sum( eta );
        eta = sqrt(eta);
        
        int nE_global;
        double nC_global;
        
        nE_global = (int)subgrid.nElements();
        subDecomposition.scalar_sum( nE_global );

        nC_global = 0.;
        for( int k = 0; k < (int)subgrid.nNodes(); k++)
            nC_global += (1./(double)subDecomposition.node2ranks[k].size());
        subDecomposition.scalar_sum( nC_global );
        
        if( comm.master() )
            printf("eta = %lf || nC = %d || nE = %d \n", eta, (int)nC_global, nE_global);
        
     
        std::string file_name = "vtu/ac";
        file_name.append(std::to_string(comm.rank()));
        file_name.append("_");
        file_name.append(std::to_string(j));
        file_name.append(".vtu");
  
        ManagedArray<int,1> partition(subgrid.nElements(), (int)comm.rank() );
        export2vtk(subgrid,u,partition,partition,file_name.c_str());
#if 0        
        //w = w_init;
        lambda_neg = vector_iteration(Y,M,w,subDecomposition);
        lambda = -lambda_neg + c_shift/eps/eps;
        printf("\n lambda = %lf on rank %d\n, c_shift = %lf",lambda,(int)comm.rank(),c_shift);
        lambda_values.push_back(lambda);
#endif
        
    }
 
#if 0 
    if(comm.master()){
        printf("\n lambda_values = \n",lambda);
        for( int k = 0; k < lambda_values.size(); k++)
            printf("%lf\n",lambda_values[k]);
    }
#endif

    double time_local;
    double time_average;
    double time_max;
    timer.getTime(time_local, time_average, time_max);
    
    if(comm.master())
        printf("adaptive algorithm, ranks = %d, time = %lf\n",(int)comm.size(),time_max/(double)nIterations_total);

    comm.finalize();
    
    return 0;
}

