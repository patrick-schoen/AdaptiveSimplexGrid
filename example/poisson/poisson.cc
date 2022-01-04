/**
@file poisson_seriell.cc
* \page Examples
* \section P0 seriell poisson problem
* \subsection P1 Implementation
* \include poisson_seriell.cc
* \subsection P2 CG algorithm
* \include cg_seriell.hh
* \subsection P3 System Matrix
* \include pdeSystem_seriell.hh
*/
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <string>
#include <limits>
#include <time.h>
        
#include "../../grid/mesh.hh"
#include "../../grid/parallel.hh"

#include "source/cg.hh"
#include "source/pdeSystem.hh"

using namespace conformingsimplexgrid;
using Communicator = MpiCommunicator;

void poisson2D(int argc, char **argv);
void poisson3D(int argc, char **argv);
void poisson2D_parallel(int argc, char **argv);
void poisson3D_parallel(int argc, char **argv);
void poisson2D_adaptive(int argc, char **argv);
void poisson3D_hierarchical(int argc, char **argv);

int main(int argc, char **argv)
{
    int test = atoi(argv[1]);
    switch( test ) {
        case 1 : printf("poisson2D\n"); poisson2D(argc, argv); break;
        case 2 : printf("poisson3D\n"); poisson3D(argc, argv); break;
        case 3 : printf("poisson2D_parallel\n"); poisson2D_parallel(argc, argv); break;
        case 4 : printf("poisson3D_parallel\n"); poisson3D_parallel(argc, argv); break;
        case 5 : printf("poisson3D_hierarchical\n"); poisson3D_hierarchical(argc, argv); break;
        case 6 : printf("poisson2D_adaptive\n"); poisson2D_adaptive(argc, argv); break;
    }
    return 0;
}

void poisson2D(int argc, char **argv) /// poisson example
{
    static const int NV = Grid_2D::NV;
    static const int DIM = Grid_2D::dim;

    static const int nRef = 12;
    static bool usePrecond = true;
    static double tol = 10E-06;

    Grid_2D grid;

    importFromFile( grid, "2dmesh_fichera.txt" );

    for(int k = 0; k < nRef; k++){
        ManagedArray<int> marked( grid.nElements(), 1);
        refine2D(grid,marked);
    }

    ManagedArray<double,1> mass;
    ManagedArray<double,NV*DIM> normals;
    outerNormals( mass, normals, grid.elements, grid.coordinates );

    //get system matrix and right hand side
    SystemMatrix<Grid_2D> A( grid, normals, mass );
    Rhs<Grid_2D> rhs( grid,  normals, mass);

    ManagedArray<double,1> u( grid.nNodes(), 0. );
    ManagedArray<double,1> x0( grid.nNodes(), 0. );
    ManagedArray<double,1> rhs_val( grid.nNodes(), 0. );
    rhs.get( rhs_val );

    pcg( A, rhs_val, x0, tol, u, usePrecond);

    std::string file_name = "vtu/poisson_2D.vtu";
    export2vtk(grid,u/u.max(),file_name.c_str());
}

void poisson3D(int argc, char **argv) /// poisson example
{
    static const int NV = Grid_3D::NV;
    static const int DIM = Grid_3D::dim;

    static const int nRef = 12;
    static bool usePrecond = true;
    static double tol = 10E-06;

    Grid_3D grid;

    importFromFile( grid, "3dmesh_fichera.txt" );

    for(int k = 0; k < nRef; k++){
        ManagedArray<int> marked( grid.nElements(), 1);
        refine(grid,marked);
    }

    ManagedArray<double,1> mass;
    ManagedArray<double,NV*DIM> normals;
    outerNormals( mass, normals, grid.elements, grid.coordinates );

    //get system matrix and right hand side
    SystemMatrix<Grid_3D> A( grid, normals, mass );
    Rhs<Grid_3D> rhs( grid, normals, mass);

    ManagedArray<double,1> u( grid.nNodes(), 0. );
    ManagedArray<double,1> x0( grid.nNodes(), 0. );
    ManagedArray<double,1> rhs_val( grid.nNodes(), 0. );
    rhs.get( rhs_val );

    pcg( A, rhs_val, x0, tol, u, usePrecond);

    std::string file_name = "vtu/poisson_3D.vtu";
    export2vtk(grid,u,file_name.c_str());
}


void poisson2D_parallel(int argc, char **argv) {
    static const int NV = Grid_2D::NV;
    static const int DIM = Grid_2D::dim;

    static const int nRef = 12;
    static bool usePrecond = true;
    static double tol = 10E-06;

    Communicator comm(&argc, &argv);
    Grid_2D grid;
    Decomposition_2D<Grid_2D,Communicator> decomposition( grid, comm );

    importFromFile( grid, "2dmesh_fichera.txt" );

    for(int k = 0; k < nRef; k++){
        ManagedArray<int> marked( grid.nElements(), 1);
        refine2D(grid,marked);
    }

    ManagedArray<double,1> mass;
    ManagedArray<double,NV*DIM> normals;
    outerNormals( mass, normals, grid.elements, grid.coordinates );

    ManagedArray<int> partition( grid.nElements(),0);
    metisDecomposition(grid, partition, (int)comm.size());  
    decomposition.init_partition( partition );
    
    {// runs in parallel
        Grid_2D subGrid;
        Decomposition_2D<Grid_2D,Communicator> subDecomposition( subGrid, comm );

        getSubGrid2D( subGrid, subDecomposition, comm.rank(), grid, decomposition, partition );
        getFacesForElements( subGrid.elements, subGrid.edges, subGrid.element2edges );

        subDecomposition.init_communication();

        ManagedArray<double,1> mass;
        ManagedArray<double,NV*DIM> normals; 
        outerNormals( mass, normals, subGrid.elements, subGrid.coordinates );

        //get system matrix and right hand side
        SystemMatrix_parallel<Grid_2D, Decomposition_2D<Grid_2D,Communicator> > A( subGrid, subDecomposition, normals, mass );
        Rhs_parallel<Grid_2D, Decomposition_2D<Grid_2D,Communicator> > rhs( subGrid, subDecomposition, normals, mass);

        ManagedArray<double,1> u( subGrid.nNodes(), 0. );
        ManagedArray<double,1> x0( subGrid.nNodes(), 0. );
        ManagedArray<double,1> rhs_val( subGrid.nNodes(), 0. );
        rhs.get( rhs_val );

        parallel_pcg( subDecomposition, A, rhs_val, x0, tol, u, usePrecond );

        std::string file_name = "vtu/poisson_2D_parallel_";
        file_name.append(std::to_string(comm.rank()));
        file_name.append(".vtu");     
        ManagedArray<int> part( subGrid.nElements(), comm.rank() );
        export2vtk(subGrid,u,part,part,file_name.c_str());
    }
    comm.finalize();
}

void poisson3D_parallel(int argc, char **argv) {
    
    Communicator comm(&argc, &argv);

    Grid_3D subGrid;
    Decomposition_3D<Grid_3D,Communicator> subDecomposition( subGrid, comm );
    GlobalNodeId<Grid_3D, Communicator, Decomposition_3D<Grid_3D,Communicator> > subGlobalId( subGrid, comm, subDecomposition );
    static const int NV = Grid_3D::NV;
    static const int DIM = Grid_3D::dim;
        
    static const int nRef = 12;
    static bool usePrecond = true;
    static double tol = 10E-06;
        
    double time1, start_time;

    if(comm.rank() == 0){
        Grid_3D grid;
        importFromFile( grid, "3dmesh_fichera.txt" );
        
        start_time = clock();
        for(int k = 0; k < nRef; ++k){
            ManagedArray<int> marked( grid.nElements(), 1);
            refine(grid,marked);
        }
            
        ManagedArray<int,1> partition;    
        metisDecomposition(grid, partition, (int)comm.size());
        
        Decomposition_3D<Grid_3D,Communicator> decomposition( grid, comm );
        GlobalNodeId<Grid_3D, Communicator, Decomposition_3D<Grid_3D,Communicator>> globalId( grid, comm, decomposition );
        globalId.init_macro( partition );
        decomposition.init_partition( partition );

        for( int rank = 0; rank < (int)comm.size(); rank++){
            if( comm.rank() == rank) 
                continue;
            getSubGrid( subGrid, subDecomposition, subGlobalId, rank, grid, decomposition, globalId, partition );
            auto sentMesh = sendMesh( subGrid, rank, 1, comm, subGlobalId, subDecomposition );
        }
        getSubGrid( subGrid, subDecomposition, subGlobalId, comm.rank(), grid, decomposition, globalId, partition );
        getFacesForElements( subGrid.elements, subGrid.faces, subGrid.element2faces, subGrid.level );
        getFacesForElements( subGrid.faces, subGrid.edges, subGrid.face2edges );  
    } else {
        recvMesh( 0, 1, comm, subGrid, subGlobalId, subDecomposition);
    }

    subDecomposition.init_communication();

    ManagedArray<double,1> mass;
    ManagedArray<double,NV*DIM> normals; 
    outerNormals( mass, normals, subGrid.elements, subGrid.coordinates );

    //get system matrix and right hand side
    SystemMatrix_parallel<Grid_3D, Decomposition_3D<Grid_3D,Communicator> > A( subGrid, subDecomposition, normals, mass );
    Rhs_parallel<Grid_3D, Decomposition_3D<Grid_3D,Communicator> > rhs( subGrid, subDecomposition, normals, mass);

    ManagedArray<double,1> u( subGrid.nNodes(), 0. );
    ManagedArray<double,1> x0( subGrid.nNodes(), 0. );
    ManagedArray<double,1> rhs_val( subGrid.nNodes(), 0. );
    rhs.get( rhs_val );

    parallel_pcg( subDecomposition, A, rhs_val, x0, tol, u, usePrecond );

    std::string file_name = "vtu/poisson_3D_parallel_";
    file_name.append(std::to_string(comm.rank()));
    file_name.append(".vtu");     
    ManagedArray<int> part( subGrid.nElements(), comm.rank() );
    export2vtk(subGrid,u,part,part,file_name.c_str());

    comm.finalize();
}

template< class Grid, class Normals, class Mass, class Vector>
Vector jump_norm( Grid &grid, Normals &normals, Mass & mass, Vector &u ){          
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;

    Faculty<DIM> faculty;
    double fac = (double)(1./faculty.value);

    ManagedArray<double, 1> res( grid.nElements(), 0. );
    ManagedArray<double, 1 > Du_nE( grid.nFaces(), 0. );

    for( size_t el = 0; el < grid.elements.size(); el ++ )
        for (size_t k = 0; k < NV; k++){
            const int face = grid.element2faces(el,(NV-1)-k);
            double temp = 0;
            for (size_t j = 0; j < NV; j++){
                for( size_t p = 0; p < DIM; p++){
                    temp += u[grid.elements(el,j)]*normals[el*NV*DIM+j*DIM+p]*normals[el*NV*DIM+k*DIM+p];
                }
            }
            Du_nE[face] += (1/mass[el]) * fac * temp;
        }
        
    res.fill(0.);
    for( size_t el = 0; el < grid.elements.size(); el ++ )
        for (size_t k = 0; k < NV; k++){
            const int face = grid.element2faces(el,k);
            res[el] += Du_nE[face]*Du_nE[face];
        }
    return res;
}

void poisson2D_adaptive(int argc, char **argv) /// poisson example
{
    static const int NV = Grid_2D::NV;
    static const int DIM = Grid_2D::dim;

    static const int nRef = 12;
    static bool usePrecond = true;
    static double tol = 10E-06;

    Grid_2D grid;

    importFromFile( grid, "2dmesh_fichera.txt" );

    ManagedArray<int> marked( grid.nElements(), 1);
    ManagedArray<double> u( grid.nNodes(), 0. );

    for(int k = 0; k < 5; k++){
        marked.resize(grid.nElements());
        marked.fill(1);
        refine2D(grid,marked);
    }
        
    for(int k = 0; k < nRef; k++){
        
        ManagedArray<double,1> mass;
        ManagedArray<double,NV*DIM> normals;
        outerNormals( mass, normals, grid.elements, grid.coordinates );
        SystemMatrix<Grid_2D> A( grid, normals, mass );
        Rhs<Grid_2D> rhs( grid,  normals, mass);
        
        u.resize(grid.nNodes()); u.fill(0.);
        ManagedArray<double,1> x0( grid.nNodes(), 0. );
        ManagedArray<double,1> rhs_val( grid.nNodes(), 0. );
        rhs.get( rhs_val );
        
        pcg( A, rhs_val, x0, tol, u, usePrecond);
        
        ManagedArray<double, 1> eta;
        eta = jump_norm(grid,normals,mass,u);
        double eta_max = eta.max();
        
        printf("eta_max = %lf\n", eta_max);
        if( eta_max < 1E-3 || k > nRef-2 ) break;
        
        marked.resize( grid.nElements() );
        marked.fill(0);

        for( int el = 0; el < grid.nElements(); el++ )
            if( eta[el] > eta_max / 2. )
                marked[el] = 1;
            
        refine2D(grid,marked);
    }

    std::string file_name = "vtu/poisson_2D_adaptive.vtu";
    export2vtk(grid,u/u.max(),file_name.c_str());
}

template< class Grid >
class DataManagerPoisson
{
    public:
        
    DataManagerPoisson( Grid &grid_in, ManagedArray<double,1> &u_in ) : grid(grid_in), u(u_in) {}

    static const int NV = Grid::NV;
    static const int DIM = Grid::DIM;

    Grid &grid;
    ManagedArray<double> &u;
        
    template< class Refine_Data >
    void pre_refine( Refine_Data &data )
    {
        int nNewNodes = 0;
        for( size_t k = 0; k < data.edge2newNode.size(); k++ )
            if( data.edge2newNode[k] >= 0 )
                nNewNodes++;
        u.resize( u.size() + nNewNodes );
            
        std::vector<bool> nodeFlag( u.size() );
        std::fill( nodeFlag.begin(), nodeFlag.end(), false );
        for(size_t el = 0; el < data.element2newNode.size(); el++)
            if( data.element2newNode[el] >= 0 ){
                const int y = data.element2newNode[el];
                if( nodeFlag[ y ] == false ){
                    const int e1 = grid.elements(el,0);
                    const int e2 = grid.elements(el,NV-1);
                    u[y] = .5*(u[e1] + u[e2]);
                    nodeFlag[y] = true;
                }
            }
    }

    template< class Refine_Data >
    void refine( Refine_Data & )
    {
    }

    template< class Coarse_Data >
    void coarse( Coarse_Data &data ) {
        u = rearange_inverse_delete( u, data.node2newNodes );
    }

    template< class Vector >
    void newIndices( Vector &, Vector &, Vector &node2newNodes  ) {
        u = rearange( u, node2newNodes );
    }

};

void poisson3D_hierarchical(int argc, char **argv) /// poisson example
{
    static const int NV = Grid_3D::NV;
    static const int DIM = Grid_3D::dim;

    static const int nRef = 12;
    static bool usePrecond = true;
    static double tol = 10E-06;

    Grid_3D grid;

    importFromFile( grid, "3dmesh_fichera.txt" );

    //get system matrix and right hand side
    ManagedArray<double,1> mass;
    ManagedArray<double,NV*DIM> normals;

    ManagedArray<double,1> u( grid.nNodes(), 0. );
    DataManagerPoisson<Grid_3D> dataManager(grid, u);    

    for(int k = 0; k < 4; k++){    
        ManagedArray<int> marked( grid.nElements(), 1);
        refine(grid,marked,dataManager);
    }
        
    for(int k = 0; k < nRef; k++){    
        ManagedArray<int> marked( grid.nElements(), 1);
        refine(grid,marked,dataManager);
        
        outerNormals( mass, normals, grid.elements, grid.coordinates );
        SystemMatrix<Grid_3D> A( grid, normals, mass );
        Rhs<Grid_3D> rhs( grid, normals, mass);
        
        ManagedArray<double,1> x0(u);
        ManagedArray<double,1> rhs_val;
        rhs.get( rhs_val );

        pcg( A, rhs_val, x0, tol, u, usePrecond);
    }

    std::string file_name = "vtu/poisson_3D.vtu";
    export2vtk(grid,u,file_name.c_str());
}


