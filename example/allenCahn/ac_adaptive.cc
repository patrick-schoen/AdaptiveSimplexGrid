#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <string>
#include <limits>
#include <time.h>
        
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

int main(int argc, char **argv)
{
    Communicator comm(&argc, &argv);
    
    int ref = 12;

    Grid subgrid;
    Decomposition subDecomposition( subgrid, comm );
    GlobalId subGlobalId( subgrid, comm, subDecomposition );
    SpaceFillingCurve<Communicator,Decomposition,Grid> curve_( comm, subDecomposition, subgrid );
    
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    double time1, start_time;
    
    if(comm.rank() == 0){

        printf("\nAllen Cahn Example, parallel and adaptive with %d ranks.\n\n", (int)comm.size() );
        
        Grid grid;
        
        importFromFile( grid, "3dmesh_sorted.txt" );
        
        SpaceFillingCurve<Communicator,Decomposition,Grid> curve( comm, subDecomposition, grid );
        
        ManagedArray<int> indices_curve(grid.nElements());
        for(int k = 0; k < grid.nElements(); k++)
            indices_curve[k] = k;
        curve.init(indices_curve);
  
        for(int k = 0; k < ref; ++k){
            ManagedArray<int> marked( grid.nElements(), 1);
            refine(grid,marked,curve);
        }
        curve_.init(curve.indices);
  
        ManagedArray<int,1> partition;    
        //metisDecomposition(grid, partition, (int)comm.size());
        spaceFillingCurvePartition_master(comm, grid, curve, partition); 
         
        Decomposition decomposition( grid, comm );
        GlobalId globalId( grid, comm, decomposition );
        globalId.init_macro( partition );

        decomposition.init_partition( partition );
        
        std::string file_name = "vtu/marco";
        file_name.append(std::to_string(comm.rank())); 
        file_name.append(".vtu");
        ManagedArray<double> v(grid.nNodes());
        for(size_t k = 0; k < grid.nEdges(); k++){
            v[grid.edges(k,0)] += (double)decomposition.edge2ranks[k].size();
            v[grid.edges(k,1)] += (double)decomposition.edge2ranks[k].size();
        }
        for(size_t k = 0; k < grid.nNodes(); k++)
           if( v[k] > 0 ) v[k] = 1.; 
        
        ManagedArray<int,1> marked(grid.nElements(), (int)comm.rank() );
        export2vtk(grid,v,partition,partition,file_name.c_str());
                    
        start_time = clock();
        for( int rank = 0; rank < (int)comm.size(); rank++){
            if( comm.rank() == rank) 
                continue;
            getSubGrid( subgrid, subDecomposition, subGlobalId, rank, grid, decomposition, globalId, partition, curve_ );
            auto sentMesh = sendMesh( subgrid, rank, 1, comm, subGlobalId, subDecomposition, curve_ );
        }
        getSubGrid( subgrid, subDecomposition, subGlobalId, comm.rank(), grid, decomposition, globalId, partition, curve_ );
        curve_.init(curve.indices_buffer);
        getFacesForElements( subgrid.elements, subgrid.faces, subgrid.element2faces, subgrid.level );
        getFacesForElements( subgrid.faces, subgrid.edges, subgrid.face2edges );
        time1 = (float) (clock() - start_time) / CLOCKS_PER_SEC;
  
    } else {
        recvMesh( 0, 1, comm, subgrid, subGlobalId, subDecomposition, curve_);
    }
    
    curve_.conjoin();
        
    ManagedArray<double,1> mass;
    ManagedArray<double,NV*DIM> normals; 
    
    outerNormals( mass, normals, subgrid.elements, subgrid.coordinates );
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
    AC_Manager<Communicator,Grid,SpaceFillingCurve<Communicator,Decomposition,Grid>> ac_manager( comm, subgrid, u, u_old, curve_ );
    
    //set initial state
    u = u_0_dumbbell<DIM>( subgrid.coordinates, eps);
    u_old = u;
    
    double t = 0;
    double time_global;
    
    subDecomposition.init_communication();
    
    Timer<Decomposition> timer(subDecomposition);
    
    int nRefinements = 1;
    int nIterations_total = 50;
    
    double TOL_COARSE = .5;
    double TOL_REFINE = 2.5;
    if(argc > 1) TOL_REFINE = atof(argv[2]);
    if(argc > 2) TOL_COARSE = atof(argv[3]);
    
    if(argc > 3) nRefinements = atoi(argv[4]);
    if(argc > 4) nIterations_total = atoi(argv[5]);
    
    if(comm.rank() == 0){ printf("initial refinement steps: %d\n",ref);}
    
    bool adaptive = true;
    
    if(adaptive) {
    
        while(1) {
            u = u_0_dumbbell<DIM>( subgrid.coordinates, eps);  

            ManagedArray<double> h_T_2 = norm_T.h_T_2();
            ManagedArray<double> f_u = f(u);
            ManagedArray<double> xi = (1/eps/eps)*f_u;
            
            ManagedArray<double> norm_l2  = norm_T.l2( xi );            
            eta_space.resize( subgrid.nElements() );
            eta_space.fill(0.);
            for( std::size_t i = 0; i < eta_space.size(); i++)
                eta_space[i] = h_T_2[i]*norm_l2[i];
            
            eta_space += norm_T.h1_jump( u );

            double eta = eta_space.sum();
            subDecomposition.scalar_sum( eta );
            eta = sqrt(eta);
                        
            if( eta < TOL_REFINE ){
                if( comm.master() )
                    printf("accept initial refinement u_0 \n");
                break;
            }
            
            marked.resize( subgrid.nElements() );
            marked.fill(0);
            double eta_max = eta_space.max();
            subDecomposition.scalar_max( eta_max );
            
            for( int el = 0; el < subgrid.nElements(); el++ )
                if( eta_space[el] > eta_max / 3. )
                    marked[el] = 1;
            
            int nE_global;
            double nC_global;
            
            nE_global = (int)subgrid.nElements();
            subDecomposition.scalar_sum( nE_global );

            nC_global = 0.;
            for( int k = 0; k < (int)subgrid.nNodes(); k++)
                nC_global += (1./(double)subDecomposition.node2ranks[k].size());
            subDecomposition.scalar_sum( nC_global );            
            
            double u_max = u.max();
            double u_min = u.min();
            
            subDecomposition.scalar_max( u_max ); 
            subDecomposition.scalar_min( u_min );
            
            refine_parallel(subDecomposition,subgrid,marked,ac_manager);
            
            ManagedArray<int> partition;
            //parMetisDecomposition(comm, subDecomposition, subgrid, partition);
            spaceFillingCurvePartition(comm, subgrid, curve_, partition);

            //subGlobalId.communicate();
            redistributeNewPartition( subgrid, comm, subDecomposition, subGlobalId, partition, ac_manager);
                 
            outerNormals( mass, normals, subgrid.elements, subgrid.coordinates );
        }
    } else {
        if(comm.master()) printf("\n\n initial uniform refinement\n");
        for( int rr = 0; rr < 10; rr++ ) {
            marked.resize( subgrid.nElements() );
            marked.fill(1);
            refine_parallel(subDecomposition,subgrid,marked,ac_manager);
        }
    }
    
    ManagedArray<int> partition;
    //parMetisDecomposition(comm, subDecomposition, subgrid, partition);
    spaceFillingCurvePartition(comm, subgrid, curve_, partition);

    //subGlobalId.communicate();
    redistributeNewPartition( subgrid, comm, subDecomposition, subGlobalId, partition, ac_manager);
        
    outerNormals( mass, normals, subgrid.elements, subgrid.coordinates );
    
    u = u_0_dumbbell<DIM>( subgrid.coordinates, eps);
    u_old = u;

    int iter;
    
    std::vector<int> numberOfNodes;
    std::vector<int> numberOfElements;
    
    for( int j = 0; j < nIterations_total; j++ ) { 
        
        if(adaptive) {
        outerNormals( mass, normals, subgrid.elements, subgrid.coordinates ); 
        //solve and refine
        iter = 0;
        u_old = u;
        
        while(1) {
            //solve 
            rhs.get( rhs_val );
            
            ManagedArray<double> xstart(u);
            //parallel_pcg( subDecomposition, A, rhs_val, xstart, TOL_CG, u );
            parallel_pcg_opt( subDecomposition, comm, A, rhs_val, u, TOL_CG );

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
            
            iter++;
            if( eta < TOL_REFINE ){
                if( iter > 1 ){
                    if( comm.master() )
                        printf("eta = %lf, iter = %d, accept refinement j = %d \n", eta, iter,j);
                    break;
                }
            }

            //mark
            marked.resize( subgrid.nElements() );
            marked.fill(0);
            double eta_max = eta_space.max();
            subDecomposition.scalar_max( eta_max );  
            for( int el = 0; el < subgrid.nElements(); el++ )
                if( eta_space[el] > eta_max / 2. )
                    marked[el] = 1;

            //refine
            refine_parallel(subDecomposition,subgrid,marked,ac_manager);
            outerNormals( mass, normals, subgrid.elements, subgrid.coordinates ); 
        }

        u_old = u;
        
        int nE_global;
        double nC_global;

        nE_global = (int)subgrid.nElements();
        subDecomposition.scalar_sum( nE_global );

        nC_global = 0.;
        for( int k = 0; k < (int)subgrid.nNodes(); k++)
            nC_global += (1./(double)subDecomposition.node2ranks[k].size());
        subDecomposition.scalar_sum( nC_global );
        
        if(comm.master())
            printf("after refinement : nElements_global = %d, dof_global = %d \n", (int)nE_global, (int)nC_global);
        
        numberOfElements.push_back(nE_global);
        numberOfNodes.push_back((int)nC_global);
        
        std::string file_name = "vtu/ac";
        file_name.append(std::to_string(comm.rank()));
        file_name.append("_");
        file_name.append(std::to_string(j));
        file_name.append(".vtu");
        ManagedArray<int,1> part(subgrid.nElements(), (int)comm.rank() );
        assert(u.size() == subgrid.nNodes());
        export2vtk(subgrid,u,part,part,file_name.c_str());
        t = t + tau;
        if( comm.master() )
            printf("t = %lf\n", t);
        
        //coarse one step
        ManagedArray<double> h_T_2 = norm_T.h_T_2();
        ManagedArray<double> f_u = f(u);
        ManagedArray<double> xi = (1/eps/eps)*f_u;
        ManagedArray<double> norm_l2  = norm_T.l2( xi );
        
        eta_space.resize( subgrid.nElements() );
        eta_space.fill(0.);
        for( std::size_t i = 0; i < eta_space.size(); i++)
            eta_space[i] = h_T_2[i]*norm_l2[i];
        
        eta_space += norm_T.h1_jump( u ); 
        double eta_max = eta_space.max();
        subDecomposition.scalar_max( eta_max ); 
        
        marked.resize( subgrid.nElements() );
        marked.fill(0);
        for( int el = 0; el < subgrid.nElements(); el++)
            marked[el] = (eta_space[el] < .5 * eta_max ) ? 1 : 0;
        
        coarse_parallel(subDecomposition,subgrid,marked,ac_manager);
#if 0
        //coarse until eta < TOL_COARSE
        iter = 0;
        double eta_c = 0.;
        while(eta_c < TOL_COARSE) {
            
            outerNormals( mass, normals, subgrid.elements, subgrid.coordinates ); 
            
            ManagedArray<double> h_T_2 = norm_T.h_T_2();
            ManagedArray<double> f_u = f(u);
            ManagedArray<double> xi = (1/eps/eps)*f_u;
            ManagedArray<double> norm_l2  = norm_T.l2( xi );
            
            eta_space.resize( subgrid.nElements() );
            eta_space.fill(0.);
            for( std::size_t i = 0; i < eta_space.size(); i++)
                eta_space[i] = h_T_2[i]*norm_l2[i];
            
            eta_space += norm_T.h1_jump( u );
            double eta_max = eta_space.max();
            subDecomposition.scalar_max( eta_max );
            
            if( comm.master() )
                printf("coarse, eta_max = %lf, nE = %d\n",eta_max,subgrid.nElements());
            
            if( iter > 0)
                break;
            iter++;
            
            marked.resize( subgrid.nElements() );
            
            double delta_coarse = 0.0;
            double eta_cc = 0.0;
            while( delta_coarse < .5 & (eta_c + eta_cc) < TOL_COARSE ) {
                marked.fill(0);
                delta_coarse += .1;
                for( int el = 0; el < subgrid.nElements(); el++)
                    marked[el] = (eta_space[el] < delta_coarse * eta_max ) ? 1 : 0;
                eta_coarse = norm_T.coarse_norm( subgrid, u, marked );
                eta_cc = eta_coarse.sum();
                subDecomposition.scalar_sum( eta_cc );
                if( comm.master() ) printf("coarse marking, eta_c = %lf\n",eta_cc);
                eta_c += eta_cc;
                if( eta_c >= TOL_COARSE ){
                    for( int el = 0; el < subgrid.nElements(); el++)
                        marked[el] = (eta_space[el] < (delta_coarse - .1) * eta_max ) ? 1 : 0;
                }
            }
            //repartitionForCoarsening( subgrid, comm, subDecomposition, subGlobalId, marked ,ac_manager );
            coarse_parallel(subDecomposition,subgrid,marked,ac_manager);
        }      
#endif
        nE_global = (int)subgrid.nElements();
        subDecomposition.scalar_sum( nE_global );

        nC_global = 0.;
        for( int k = 0; k < (int)subgrid.nNodes(); k++)
        nC_global += (1./(double)subDecomposition.node2ranks[k].size());
        subDecomposition.scalar_sum( nC_global );   

        if( comm.master() )
            printf("after coarsening, nElements = %d, dof = %d, nElements_global = %d, dof_global = %d\n",(int)subgrid.nElements(),(int)subgrid.nNodes(),(int)nE_global,(int)nC_global);
    
        ManagedArray<int> partition;
        //parMetisDecomposition(comm, subDecomposition, subgrid, partition);
        spaceFillingCurvePartition(comm, subgrid, curve_, partition);

        subGlobalId.communicate();

        redistributeNewPartition( subgrid, comm, subDecomposition, subGlobalId, partition, ac_manager);

        u = u_old;
        
        // uniform refinement
        } else {
            rhs.get( rhs_val );
            ManagedArray<double> xstart(u);
            parallel_pcg( subDecomposition, A, rhs_val, xstart, TOL_CG, u );
            //parallel_pcg_opt( subDecomposition, comm, A, rhs_val, u, TOL_CG );
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
            
            double u_max = u.max();
            double u_min = u.min();

            subDecomposition.scalar_max( u_max ); 
            subDecomposition.scalar_min( u_min );

            if( comm.master() )
                printf("u_max = %lf, u_min = %lf, eta = %lf, nE = %d, nC = %d, nE_global = %d, nC_global = %d\n",u_max,u_min,eta,
                    (int)subgrid.nElements(),(int)subgrid.nNodes(),(int)nE_global,(int)nC_global);
            
        }
    }
    
    comm.finalize();
    return 0;
}


