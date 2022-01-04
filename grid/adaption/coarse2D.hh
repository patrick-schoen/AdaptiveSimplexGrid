#ifndef COARSE_PARALLEL_HH
#define COARSE_PARALLEL_HH

#include "../simplexGrid.hh"

namespace conformingsimplexgrid {

    template< class Grid, class VectorBool, class Vector > 
    struct ProlongationData_2D_coarse
    { 
        ProlongationData_2D_coarse( Grid & grid_in, VectorBool & remainingElements_in, Vector &node2newNodes_in )
        : grid(grid_in), remainingElements(remainingElements_in), node2newNodes(node2newNodes_in)
        {}
        
        Grid &grid;
        VectorBool &remainingElements;
        Vector &node2newNodes;
    };
          
    template< class Grid, class ElementVector, class ... Manager >
    void coarse2D( Grid &grid, ElementVector &marked, Manager & ... manager )
    {  
        grid.numNodes = grid.coordinates.size();
        grid.numEdges = grid.edges.size();
        
        ManagedArray<int> markedNodes( grid.nNodes(), 0 );
        findRemovableNodes2D( grid, marked, markedNodes, manager ... );
        marked.fill(0);
        for(size_t el = 0; el < marked.size(); ++el)
            if( markedNodes[ grid.elements(el,1) ] == 1 )
                marked[el] = 1; 
        coarseNodes2D( grid, marked, manager ... );
    }
    
    template< class Grid, class ElementVector, class NodeVector >
    void findRemovableNodes2D( Grid &grid, ElementVector &marked, NodeVector &markedNodes)
    {  
        const int NV = Grid::NV;
        
        for(size_t el = 0; el < marked.size(); ++el)
            if( marked[el] > 0 && grid.level[el] > 0 )
                markedNodes[ grid.elements(el,1) ] = 1;
            
        for(size_t el = 0; el < marked.size(); ++el){
            markedNodes[ grid.elements(el,0) ] = 0;
            for( int k = 2; k < NV; ++k)
                markedNodes[ grid.elements(el,k) ] = 0;
        }
        marked.fill(0);
        for(size_t el = 0; el < marked.size(); ++el)
            if( markedNodes[ grid.elements(el,1) ] == 1 )
                marked[el] = 1;  
    }
    
    template< class Grid, class ElementVector, class Vector >
    void getBrothers2D( Grid &grid, ElementVector &marked, Vector &left, Vector &right )
    {     
        const int NV = Grid::NV;
        
        int nMarked = 0;
        for(size_t el = 0; el < marked.size(); ++el)
            if( marked[el] > 0 )
                nMarked++;
        
        left.resize( nMarked / 2 ); left.fill(0);
        right.resize( nMarked / 2 ); right.fill(0);
        ManagedArray<int> newestEdgesData( grid.nEdges(), -1);
        for(size_t el = 0, counter = 0; el < marked.size(); ++el){
            if( marked[el] > 0 ){
                const int e0 = grid.element2edges(el,2);
                if( newestEdgesData[e0] > -1){
                    left[counter] = newestEdgesData[e0];
                    right[counter] = el;
                    counter++;
                } else {
                    newestEdgesData[e0] = el;
                }
            }
        } 

    }
    
    template< class Grid, class ElementVector, class ... Manager >
    void coarseNodes2D( Grid &grid, ElementVector &marked, Manager & ... manager )
    {        
        const int NV = Grid::NV;
        const int DIM = Grid::DIM;

        ManagedArray<int> left, right;
        getBrothers2D( grid, marked, left, right );
        
        std::vector<bool> remainingElements( grid.nElements(), true );
        ManagedArray<int> node2newNodes( grid.nNodes(), 0 );
        
        for( size_t k = 0; k < right.size(); k++)
            remainingElements[right[k]] = false;

        coarse_elements2D(grid,left,right,remainingElements,node2newNodes);
        coarse_level2D(grid,left,right,remainingElements);
        getElement2Faces( grid.elements, grid.edges, grid.element2edges );
        
        //ProlongationData_coarse<Grid, std::vector<bool>, ManagedArray<int> > pdata( grid, remainingElements,
        //    face2newFaces, edge2newEdges, node2newNodes );
        
        //coarse_Manager( pdata, manager ... );
        
        size_t counter = 0;
        for(size_t k = 0; k < node2newNodes.size(); ++k){
            if( node2newNodes[k] < 0 )
                continue;
            for(size_t d = 0; d < DIM; ++d)
                grid.coordinates(counter,d) = grid.coordinates(k,d);
            counter++;
        }
        grid.coordinates.resize(counter);
    }
      
    template< class Grid, class Vector, class BoolVector, class NodeVector >
    void coarse_elements2D( Grid &grid, Vector &left, Vector &right, BoolVector &remainingElements, NodeVector &nodes2newNodes ){
        const int NV = Grid::NV;
       
        for(size_t k = 0; k < left.size(); ++k){
            grid.elements(left[k],1) = grid.elements(left[k],2);
            grid.elements(left[k],2) = grid.elements(right[k],0);              
        }
        
        int counter = 0;
        for(size_t el = 0; el < grid.nElements(); ++el)
            if( remainingElements[el] ){
                for( size_t k = 0; k < NV; ++k)
                    grid.elements(counter,k) = grid.elements(el,k);
                counter++;
            }
        grid.elements.resize(counter);
        
        nodes2newNodes.fill(-1);
        for(size_t el = 0; el < grid.elements.size(); ++el)
            for( size_t k = 0; k < NV; ++k)
                nodes2newNodes[ grid.elements(el,k) ] = 1;
            
        for(size_t k = 0, counter = 0; k < nodes2newNodes.size(); ++k){
            if( nodes2newNodes[k] < 0 )
                continue;
            nodes2newNodes[k] = counter++;
        }
         
        for(size_t el = 0; el < grid.elements.size(); ++el)
            for( size_t k = 0; k < NV; ++k){
                assert( nodes2newNodes[grid.elements(el,k)] >= 0 );
                grid.elements(el,k) = nodes2newNodes[grid.elements(el,k)];
            }
    }
     
    template< class Grid, class Vector, class BoolVector >
    void coarse_level2D( Grid &grid, Vector &left, Vector &right, BoolVector &remainingElements ){
        for( size_t k = 0; k < left.size(); ++k)
            grid.level[left[k]] = grid.level[left[k]] - 1;
        
        int counter = 0;
        for(size_t k = 0; k < grid.level.size(); ++k)
            if( remainingElements[k] )
                grid.level[counter++] = grid.level[k];
        grid.level.resize(counter);          
    }
    
    
}

#endif //COARSE_HH