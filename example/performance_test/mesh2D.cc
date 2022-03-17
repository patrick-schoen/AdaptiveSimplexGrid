#include <stdio.h>
#include <stdlib.h>
 
#include "../grid/mesh.hh"

using namespace conformingsimplexgrid;
using Grid = Grid_2D;

int main(int argc, char **argv)
{ 
    Grid grid;
    importFromFile( grid, "grid_files/2dmesh.txt" );
    
    for( int rr = 0; rr < 12; rr++ ) {
        ManagedArray<int> marked( grid.nElements(), 1);
        marked[0] = 1;
        refine2D(grid,marked);
    }
    
    ManagedArray<double> u(grid.nNodes(),0.);
    
    export2vtk(grid,u,"output/refine.vtu");
    
    for( int rr = 0; rr < 12; rr++ ) {
        ManagedArray<int> marked( grid.nElements(),1);
        coarse2D(grid,marked);
    }
    
    ManagedArray<double> u(grid.nNodes(),0.);
    
    export2vtk(grid,u,"output/coarse.vtu");
    
}