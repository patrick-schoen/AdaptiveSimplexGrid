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

template< class Grid, class Communicator, class Decomposition, class GlobalId >
void bar_test(Grid &subgrid, Communicator &comm, Decomposition &subDecomposition, GlobalId &subGlobalId )
{
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    if(comm.master()) printf("bar_test with %d ranks \n",(int)comm.size());
    
    double nE_fine = 1E07;
    double nE_coarse = 1E07 / 8.;
     
    subDecomposition.init_communication();
    
    Timer<Decomposition> timer(subDecomposition);
    
    double tau = .1;
    double t = 0.;
    double T = 1.;
    
    double x_j = 0.;
    double delta_t = .25;
    x_j += delta_t;
    
    int nE_global, nC_global;
    
    int nRanks = (int)comm.size();
    int thisRank = (int)comm.rank();
    
    ManagedArray<char> marked( subgrid.nElements(), 0);

    std::vector<double> refine_time, coarse_time, marking_time, nE_global_vec, nC_global_vec, numberIterations_refine, numberIterations_coarse;
    
    for( int j = 0; j < 5; j++ )
    {        

        if(comm.master()) printf("iteration j = %d\n", j);
        
        nC_global = (int)subgrid.nNodes();
        subDecomposition.scalar_sum( nC_global );
        
        nE_global = (int)subgrid.nElements();
        subDecomposition.scalar_sum( nE_global );
       
        int counter = 0;
        double Time_refine = 0.;
        double Time_marking = 0.;        
        int nIterations_ref = 0;
        
        for(int iter_refine = 0; iter_refine < 20; iter_refine++) {                  
            marked.resize( subgrid.nElements() );
            marked.fill(0);
            for( int el = 0; el < subgrid.nElements(); el++){
                ManagedArray<double> z(4,0.);
                z[0] = subgrid.coordinates( subgrid.elements(el,0), 0 );
                z[1] = subgrid.coordinates( subgrid.elements(el,1), 0 );
                z[2] = subgrid.coordinates( subgrid.elements(el,2), 0 );
                z[3] = subgrid.coordinates( subgrid.elements(el,3), 0 );
                if( z.max() >= x_j && z.min() <= x_j )
                    marked[el] = 1;
            }
            
            int nMarked_interfaces = 0;
            for (int el = 0; el < subgrid.nElements(); el++)
                if( marked[el]> 0 )
                    nMarked_interfaces++;
            
            double refineTime = 0.;
            double markingTime = 0.;
            
            nC_global = (int)subgrid.nNodes();
            subDecomposition.scalar_sum( nC_global );
            
            timer.begin();
            refine_parallel_markingTime(markingTime,subDecomposition,subgrid,marked);
            timer.getTime_local(refineTime);
            
            Time_refine += refineTime;
            Time_marking += markingTime;
            counter++;
            nC_global = (int)subgrid.nNodes();
            subDecomposition.scalar_sum( nC_global );
            nIterations_ref++;
        }
         
        nE_global = (int)subgrid.nElements();
        subDecomposition.scalar_sum( nE_global );

        
        
        if(comm.master()) printf("nElements = %d nC = %d, nIterations = %d, time refine = %lf (%lf)\n", nE_global, nC_global, nIterations_ref, Time_refine, Time_marking);
        
//         std::string file_name = "out_bar/bar_";
//         file_name.append(std::to_string(j)); 
//         file_name.append("_"); 
//         file_name.append(std::to_string(comm.rank())); 
//         file_name.append(".vtu");
//         ManagedArray<double> v(subgrid.nNodes(),0.);
//         //for(size_t k = 0; k < grid.nEdges(); k++){
//         //    v[grid.edges(k,0)] += (double)decomposition.edge2ranks[k].size();
//         //    v[grid.edges(k,1)] += (double)decomposition.edge2ranks[k].size();
//         //}
//         //for(size_t k = 0; k < grid.nNodes(); k++)
//         //   if( v[k] > 0 ) v[k] = 1.; 
//         ManagedArray<int> partition( subgrid.nElements(), comm.rank());
//         export2vtk(subgrid,v,partition,partition,file_name.c_str());
//         
//         
        
        nE_global_vec.push_back(nE_global);
        nC_global_vec.push_back(nC_global);
        refine_time.push_back(Time_refine);
        marking_time.push_back(Time_marking);
        numberIterations_refine.push_back(nIterations_ref);
          
        t = t + tau;
        x_j += delta_t;
        
        counter = 0;
        double Time_coarse = 0.;
        
        nC_global = (int)subgrid.nNodes();
        subDecomposition.scalar_sum( nC_global );

        nE_global = (int)subgrid.nElements();
        subDecomposition.scalar_sum( nE_global );
        
        int nIterations_coarse = 0;
        
        while (true) {
         
            marked.resize( subgrid.nElements() );
            marked.fill(1);
            for( int el = 0; el < subgrid.nElements(); el++){
                ManagedArray<double> z(4,0.);
                z[0] = subgrid.coordinates( subgrid.elements(el,0), 0 );
                z[1] = subgrid.coordinates( subgrid.elements(el,1), 0 );
                z[2] = subgrid.coordinates( subgrid.elements(el,2), 0 );
                z[3] = subgrid.coordinates( subgrid.elements(el,3), 0 );
                if( z.max() >= x_j && z.min() <= x_j )
                    marked[el] = 0;
            }
                            
            double time_local1,time_average1,time_max1;
   
            double time_local;
            timer.begin();
            coarse_parallel(subDecomposition,subgrid,marked);
            assert(subgrid.level.size() == subgrid.elements.size());
            timer.getTime_local(time_local);
            
            Time_coarse += time_local,
            counter++;
            
            nIterations_coarse++;
            
            int nC_global_old = nC_global;
            nC_global = (int)subgrid.nNodes();
            subDecomposition.scalar_sum( nC_global ); 
            
            if( nC_global_old == nC_global )
                break;
        }
        coarse_time.push_back(Time_coarse);
        numberIterations_coarse.push_back(nIterations_coarse);
        
        nE_global = (int)subgrid.nElements();
        subDecomposition.scalar_sum( nE_global );

        
        if(comm.master()) printf("nElements = %d nC = %d, nIterations = %d, time coarse = %lf\n", nE_global, nC_global, nIterations_coarse, Time_coarse);
    
    }
   
        
    std::ofstream ofs_refine, ofs_coarse,ofs_marking ; 
    std::string file_refine = "out_bar/refine_" + std::to_string(nRanks) + ".txt";
    std::string file_coarse = "out_bar/coarse_" + std::to_string(nRanks) + ".txt";
    std::string file_marking = "out_bar/marking_" + std::to_string(nRanks) + ".txt";
    
    
    ofs_refine.open(file_refine, std::ofstream::out | std::ofstream::trunc);
    ofs_coarse.open(file_coarse, std::ofstream::out | std::ofstream::trunc);
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
        if( comm.master() ){
            ofs_coarse << k << "," << time_average << "," << time_max << "," << nC_global_vec[k] << "," << nE_global_vec[k] << "," << numberIterations_coarse[k] << std::endl; 
            ofs_coarse.flush();
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
    ofs_marking.close();
}


int main(int argc, char **argv)
{
    Communicator comm(&argc, &argv);
    
    int ref = 12;
 
    Grid subgrid;
    Decomposition subDecomposition( subgrid, comm );
    GlobalId subGlobalId( subgrid, comm, subDecomposition );
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    double time1, start_time;
    
    
    if(comm.rank() == 0){
        Grid grid;
        
        std::string file_bar = "3Dbar_";
        file_bar.append(std::to_string(comm.size())); 
        file_bar.append(".txt");   
        importFromFile( grid, file_bar.data() );

        ManagedArray<int,1> partition( grid.nElements(), 0 );
        
        std::size_t counter = 0;
        std::size_t elementsPerRank = (grid.nElements() / comm.size()) - 1;
        std::size_t partition_counter = 0;
        for( int k = 0; k < grid.nElements(); k++) {
                partition[k] = partition_counter;
                counter++;
                if(counter > elementsPerRank){
                    counter = 0;
                    partition_counter++;
                }
        }
      
        Decomposition decomposition( grid, comm );
        GlobalId globalId( grid, comm, decomposition );
        globalId.init_macro( partition );
        decomposition.init_partition( partition );
          
        std::string file_name = "out_bar/marco";
        file_name.append(std::to_string(comm.rank())); 
        file_name.append(".vtu");
        ManagedArray<double> v(grid.nNodes(),0.);
        //for(size_t k = 0; k < grid.nEdges(); k++){
        //    v[grid.edges(k,0)] += (double)decomposition.edge2ranks[k].size();
        //    v[grid.edges(k,1)] += (double)decomposition.edge2ranks[k].size();
        //}
        //for(size_t k = 0; k < grid.nNodes(); k++)
        //   if( v[k] > 0 ) v[k] = 1.; 
        
        export2vtk(grid,v,partition,partition,file_name.c_str());
         
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
    
    bar_test(subgrid, comm, subDecomposition, subGlobalId );
    
    comm.finalize();
    return 0;
}


