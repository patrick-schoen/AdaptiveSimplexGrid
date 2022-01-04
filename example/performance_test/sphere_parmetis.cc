#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <time.h>
        
#include "../grid/mesh.hh"
#include "../grid/parallel.hh"

#include <math.h>

using namespace conformingsimplexgrid;
using Communicator = MpiCommunicator;
using Grid = Grid_3D;
using Decomposition = Decomposition_3D<Grid,Communicator>;
using GlobalId = GlobalNodeId<Grid, Communicator, Decomposition>;

double dist_circle_0( std::array<double,3> &x, std::array<double,3> &m, double &radius){
    return sqrt( (x[0]-m[0])*(x[0]-m[0]) + (x[1]-m[1])*(x[1]-m[1]) + (x[2]-m[2])*(x[2]-m[2]) ) - radius;
}

template<int DIM>
ManagedArray<double,1> u_0_circle( ManagedArray<double,DIM> &x, std::array<double,3> &m, double &radius)
{
    double eps = 1/16;
    int nC = x.size();
    std::array<double,3> z;
    ManagedArray<double,1> val(nC);
    for(int i = 0; i < nC; i++){
       z[0] = x[i*3]; z[1] = x[i*3+1]; z[2] = x[i*3+2];
       val[i] = tanh( dist_circle_0( z, m, radius ) / eps );
    }
    return val;
}

template< class Grid, class Communicator, class Decomposition, class GlobalId >
void circle_test_local(Grid &subgrid, Communicator &comm, Decomposition &subDecomposition, GlobalId &subGlobalId )
{
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    if(comm.master()) printf("circle_test with %d ranks, parmetis \n",(int)comm.size());
    
    double nE_fine = 5.*1E06;
    double nE_coarse = 5.*1E06 / 8.;
     
    subDecomposition.init_communication();
    
    Timer<Decomposition> timer(subDecomposition);
    
    for( int k = 0; k < subgrid.coordinates.size(); k++)
        for (int d = 0; d < 3; d++)
            subgrid.coordinates(k,d) = (subgrid.coordinates(k,d) + 1.)/2.;
    
    std::array<double,3> Y;
    double tau = .1;
    double t = 0;
    double T = 1;
    
    Y[0] = .5 + (1./3.);
    Y[1] = .5;
    Y[2] = .5;
    
    int nE_global, nC_global;
    
    double radius = .175;

    ManagedArray<double> u(subgrid.nElements(),0.);
    ManagedArray<char> marked( subgrid.nElements(), 0);
    
    ManagedArray<double,1> mass;
    ManagedArray<double,NV*DIM> normals; 
    
    outerNormals( mass, normals, subgrid.elements, subgrid.coordinates );
    
    local_norm<Grid, Decomposition> norm_T( subgrid, subDecomposition, normals, mass );
    ManagedArray<double, 1> eta_jump( subgrid.elements.size() );
    
    int nRanks = (int)comm.size();
    int thisRank = (int)comm.rank();
    
    
    std::vector<double> refine_time, coarse_time, repartition_time, parmetis_time, marking_time, nE_global_vec, nC_global_vec, numberIterations_refine, numberIterations_coarse, nE_partition, repartition_time_data, repartition_time_communication, repartition_coarse;
    
    //for( int j = 0; j < T/tau; j++ )
    for( int j = 0; j < 5; j++ )
    {        
        
        if(comm.master()) printf("iteration %d\n", j);
        
        nC_global = (int)subgrid.nNodes();
        subDecomposition.scalar_sum( nC_global );
        
        nE_global = (int)subgrid.nElements();
        subDecomposition.scalar_sum( nE_global );
        
        outerNormals( mass, normals, subgrid.elements, subgrid.coordinates );
        
        int counter = 0;
        double Time_refine = 0.;
        double Time_marking = 0.;
        double Time_parMetis = 0.;
        double Time_repart = 0.;
        
        double Time_repart_data = 0.;
        double Time_repart_communcation = 0.;
        
        int nIterations_ref = 0;
        
        int nElementsPartition = 0;
        
        while( nE_global < nE_fine ) {
                    
            u = u_0_circle(subgrid.coordinates, Y, radius); 
                  
            eta_jump = norm_T.h1( u );

            double eta = eta_jump.sum();
            subDecomposition.scalar_sum( eta );
            eta = sqrt(eta);
                                    
            marked.resize( subgrid.nElements() );
            marked.fill(0);
            double eta_max = eta_jump.max();
            subDecomposition.scalar_max( eta_max );
            
            for( int el = 0; el < subgrid.nElements(); el++ )
                if( eta_jump[el] > eta_max / 2. )
                    marked[el] = 1;
          
            int nMarked = 0;
            for(int el = 0; el < subgrid.nElements(); el++)
               if (marked[el] > 0 ) nMarked++;
            int nMarked_global = nMarked;
            subDecomposition.scalar_sum( nMarked_global );
            
            double refineTime = 0.;
            double markingTime = 0.;
            timer.begin();
            refine_parallel_markingTime(markingTime,subDecomposition,subgrid,marked);
            timer.getTime_local(refineTime);
            //if(comm.master()) printf("refine Time %lf - %lf nElements_old = %d nC_old = %d \n", refineTime,markingTime, nE_global, nC_global);
            
            Time_refine += refineTime;
            Time_marking += markingTime;
            counter++;
            
            comm.barrier();
            
            nE_global = (int)subgrid.nElements();
            subDecomposition.scalar_sum( nE_global );

            nC_global = (int)subgrid.nNodes();
            subDecomposition.scalar_sum( nC_global );
            
            //if(comm.master()) printf("nElements = %d nC = %d \n", nE_global, nC_global);

            ManagedArray<int> partition;
            double time_mestis = 0.;
            timer.begin();
            parMetisDecomposition(comm, subDecomposition, subgrid, partition);
            timer.getTime_local(time_mestis);
            Time_parMetis += time_mestis;
            //if(comm.master()) printf("time metis %lf\n", time_mestis);
            
            for(std::size_t k = 0; k < partition.size(); k++){
                if( partition[k] == comm.rank() )
                    continue;
                nElementsPartition++;
            }
  
            double time_repart = 0.;
            double time_repart_data = 0.;
            double time_repart_communcation = 0.;
            double time_repart_communcation1 = 0.;
            timer.begin();
            subGlobalId.communicate();
            timer.getTime_local(time_repart_communcation1);
            timer.begin();
            redistributeNewPartition_Timer( subgrid, comm, subDecomposition, subGlobalId, partition, time_repart_data, time_repart_communcation);
            timer.getTime_local(time_repart);
            Time_repart += time_repart;
            Time_repart_data += time_repart_data;
            Time_repart_communcation += time_repart_communcation;
            Time_repart_communcation += time_repart_communcation1;
            //if(comm.master()) printf("time repart %lf\n", time_repart);
            
            outerNormals( mass, normals, subgrid.elements, subgrid.coordinates );
            
            nIterations_ref++;
        }
        
        subDecomposition.scalar_sum( nElementsPartition );
        
        if(comm.master()) printf("nElements = %d nC = %d, nIterations = %d, time refine = %lf (%lf), time metis = %lf, time repart = %lf (data %lf, comm %lf, nElements = %d) \n", nE_global, nC_global, nIterations_ref, Time_refine, Time_marking, Time_parMetis, Time_repart, Time_repart_data, Time_repart_communcation, nElementsPartition);
        
        nE_global_vec.push_back(nE_global);
        nC_global_vec.push_back(nC_global);
        refine_time.push_back(Time_refine);
        repartition_time.push_back(Time_repart);
        repartition_time_data.push_back(Time_repart_data);
        repartition_time_communication.push_back(Time_repart_communcation);
        parmetis_time.push_back(Time_parMetis);
        marking_time.push_back(Time_marking);
        numberIterations_refine.push_back(nIterations_ref);
        
        
        nE_partition.push_back(nElementsPartition);
          
        outerNormals( mass, normals, subgrid.elements, subgrid.coordinates );
             
           
        t = t + tau;
        Y[0] = .5 + (1./3.)*cos(2.*M_PI*t);
        Y[1] = .5 + (1./3.)*sin(2.*M_PI*t);

        counter = 0;
        double Time_coarse = 0.;
        double Time_repartition_coarse = 0.;
        
        nC_global = (int)subgrid.nNodes();
        subDecomposition.scalar_sum( nC_global );

        nE_global = (int)subgrid.nElements();
        subDecomposition.scalar_sum( nE_global );
        
        int nIterations_coarse = 0;
        
        while (nE_global > nE_coarse) {
         
            u = u_0_circle(subgrid.coordinates, Y, radius );         
            eta_jump = norm_T.h1( u );

            double eta = eta_jump.sum();
            subDecomposition.scalar_sum( eta );
            eta = sqrt(eta);
                                    
            marked.resize( subgrid.nElements() );
            marked.fill(0);
            double eta_max = eta_jump.max();
            subDecomposition.scalar_max( eta_max );
            
            for( int el = 0; el < subgrid.nElements(); el++ )
                //if( eta_jump[el] < eta_max / 2. )
                    marked[el] = 1;  
            
            double time_local1;
            timer.begin();
            //repartitionForCoarsening( subgrid, comm, subDecomposition, subGlobalId, marked );
            timer.getTime_local(time_local1);
            Time_repartition_coarse += time_local1,
            
            marked.fill(1);
                            
            double time_local;
            timer.begin();
            coarse_parallel(subDecomposition,subgrid,marked);
            assert(subgrid.level.size() == subgrid.elements.size());
            timer.getTime_local(time_local);
            
            Time_coarse += time_local,
            counter++;
            
            outerNormals( mass, normals, subgrid.elements, subgrid.coordinates );
            
            nIterations_coarse++;
            
            int nC_global_old = nC_global;
            
            nE_global = (int)subgrid.nElements();
            subDecomposition.scalar_sum( nE_global );

            nC_global = (int)subgrid.nNodes();
            subDecomposition.scalar_sum( nC_global ); 
            
            if(nC_global_old == nC_global) break;
        }
        coarse_time.push_back(Time_coarse);
        repartition_coarse.push_back(Time_repartition_coarse);
        
        numberIterations_coarse.push_back(nIterations_coarse);
        
        if(comm.master()) printf("nElements = %d nC = %d, nIterations = %d, time coarse = %lf, time repartition coarse = %lf\n", nE_global, nC_global, nIterations_coarse, Time_coarse, Time_repartition_coarse);
    
    }
   
        
    std::ofstream ofs_refine, ofs_coarse, ofs_repartition, ofs_parmetis,ofs_marking ; 
    std::string file_refine = "out/refine_" + std::to_string(nRanks) + ".txt";
    std::string file_coarse = "out/coarse_" + std::to_string(nRanks) + ".txt";
    std::string file_repartition = "out/repartition_" + std::to_string(nRanks) + ".txt";
    std::string file_parmetis = "out/parmetis_" + std::to_string(nRanks) + ".txt";
    std::string file_marking = "out/marking_" + std::to_string(nRanks) + ".txt";
    
    
    ofs_refine.open(file_refine, std::ofstream::out | std::ofstream::trunc);
    ofs_coarse.open(file_coarse, std::ofstream::out | std::ofstream::trunc);
    ofs_repartition.open(file_repartition, std::ofstream::out | std::ofstream::trunc);
    ofs_parmetis.open(file_parmetis, std::ofstream::out | std::ofstream::trunc);
    ofs_marking.open(file_marking, std::ofstream::out | std::ofstream::trunc);
    
    for(int k = 0; k < refine_time.size(); k++){
        double time_global = refine_time[k];
        double time_max = time_global;
        subDecomposition.scalar_sum( time_global );
        subDecomposition.scalar_max( time_max );
        double time_average = time_global / subDecomposition.size();
        if( comm.master() ){
            ofs_refine << k << "," <<  time_average << "," << time_max << "," << nC_global_vec[k] << "," << nE_global_vec[k] << "," << numberIterations_refine[k] << std::endl; 
            ofs_refine.flush();
        }  
    }
    
    for(int k = 0; k < coarse_time.size(); k++){
        double time_global = coarse_time[k];
        double time_max = time_global;
        subDecomposition.scalar_sum( time_global );
        subDecomposition.scalar_max( time_max );
        double time_average = time_global / subDecomposition.size();
        double time_global1 = repartition_coarse[k];
        subDecomposition.scalar_max( time_global1 );
        if( comm.master() ){
            ofs_coarse << k << "," << time_average << "," << time_max << "," << nC_global_vec[k] << "," << nE_global_vec[k] << "," << numberIterations_coarse[k] << "," << time_global1 << std::endl; 
            ofs_coarse.flush();
        }  
    }
    
    for(int k = 0; k < repartition_time.size(); k++){
        double time_global = repartition_time[k];
        double time_max = time_global;
        subDecomposition.scalar_sum( time_global );
        subDecomposition.scalar_max( time_max );
        double time_global1 = repartition_time_data[k];
        subDecomposition.scalar_max( time_global1);
        double time_global2 = repartition_time_communication[k];
        subDecomposition.scalar_max( time_global2);
        double time_average = time_global / subDecomposition.size();
        if( comm.master() ){
            ofs_repartition << k << "," << time_average << "," << time_max << "," << nC_global_vec[k] << "," << nE_global_vec[k] << "," << nE_partition[k] <<  "," << time_global1 << "," << time_global2 << std::endl;
            ofs_repartition.flush();
        }  
    }
    
    for(int k = 0; k < parmetis_time.size(); k++){
        double time_global = parmetis_time[k];
        double time_max = time_global;
        subDecomposition.scalar_sum( time_global );
        subDecomposition.scalar_max( time_max );
        double time_average = time_global / subDecomposition.size();
        if( comm.master() ){
            ofs_parmetis << k << "," << time_average << "," << time_max << "," << nC_global_vec[k] << "," << nE_global_vec[k] << std::endl; 
            ofs_parmetis.flush();
        }  
    }
    
    for(int k = 0; k < marking_time.size(); k++){
        double time_global = marking_time[k];
        double time_max = time_global;
        subDecomposition.scalar_sum( time_global );
        subDecomposition.scalar_max( time_max );
        double time_average = time_global / subDecomposition.size();
        if( comm.master() ){
            ofs_marking << k << "," << time_average << "," << time_max << "," << nC_global_vec[k] << "," << nE_global_vec[k] << std::endl; 
            ofs_marking.flush();
        }  
    }
    
    ofs_refine.close();
    ofs_coarse.close();
    ofs_repartition.close();
    ofs_parmetis.close();
    ofs_marking.close();
}


int main(int argc, char **argv)
{
    Communicator comm(&argc, &argv);
    
    int ref = 12;
    //if(argc > 0) ref = atoi(argv[1]);
   
    
    Grid subgrid;
    Decomposition subDecomposition( subgrid, comm );
    GlobalId subGlobalId( subgrid, comm, subDecomposition );
    SpaceFillingCurve<Communicator,Decomposition,Grid> curve_( comm, subDecomposition, subgrid );
    
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    double time1, start_time;
  
    
    if(comm.rank() == 0){
        Grid grid;
        SpaceFillingCurve<Communicator,Decomposition,Grid> curve( comm, subDecomposition, grid );
      
        importFromFile( grid, "3dmesh_sorted.txt" );
        
        ManagedArray<int> indices_curve(grid.nElements());
        for(int k = 0; k < grid.nElements(); k++)
            indices_curve[k] = k;
        curve.init(indices_curve);
        
        start_time = clock();
        for(int k = 0; k < ref; ++k){
            ManagedArray<int> marked( grid.nElements(), 1);
            refine(grid,marked,curve);
        };
        curve_.init(curve.indices);
        
        time1 = (float) (clock() - start_time) / CLOCKS_PER_SEC; 
        //printf("time for macro refine: %f seconds.\n", time1);
        
        start_time = clock();
        double time1 = (float) (clock() - start_time) / CLOCKS_PER_SEC; 
           
        //printf("\n parallel grid with %d nodes, %d elements and %d ranks \n", (int)grid.nNodes(),(int)grid.nElements(), (int)comm.size() );
        
        ManagedArray<int,1> partition;    
        start_time = clock();
        metisDecomposition(grid, partition, (int)comm.size());
        //spaceFillingCurvePartition_master(comm, grid, curve, partition); 
        time1 = (float) (clock() - start_time) / CLOCKS_PER_SEC; 
        //printf("time for metisDecomposition: %f seconds.\n", time1);
        
        start_time = clock();
        Decomposition decomposition( grid, comm );
        GlobalId globalId( grid, comm, decomposition );
        globalId.init_macro( partition );
        
        start_time = clock();
        decomposition.init_partition( partition );
        time1 = (float) (clock() - start_time) / CLOCKS_PER_SEC; 
        //printf("time setting up decomposition: %f seconds.\n", time1);  
        
         
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
        //printf("time sending submeshs: %f seconds.\n", time1);
        
    } else {
        recvMesh( 0, 1, comm, subgrid, subGlobalId, subDecomposition, curve_);
    }
    curve_.conjoin();
#if 1
    ManagedArray<int> partition(subgrid.nElements());
    partition.fill(comm.rank());
    std::string file_name = "vtu/marco";
    file_name.append(std::to_string(comm.rank())); 
    file_name.append(".vtu");
    ManagedArray<double> v(subgrid.nNodes());
    export2vtk(subgrid,v,partition,partition,file_name.c_str());
#endif
    circle_test_local(subgrid, comm, subDecomposition, subGlobalId );
    
    comm.finalize();
    return 0;
}


