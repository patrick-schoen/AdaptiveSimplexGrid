#ifndef REFINE_2D_HH
#define REFINE_2D_HH

#include "../array/managed_array.hh"

namespace conformingsimplexgrid {
    
    template< class Vector > 
    struct ProlongationData_2D_pre_refine
    {
        ProlongationData_2D_pre_refine( Vector & element2newNode_in, Vector & edge2newNode_in )
        : element2newNode(element2newNode_in), edge2newNode(edge2newNode_in)
        {}
        
        Vector &element2newNode;
        Vector &edge2newNode;
    };
        
    class refine_scheme_2D {
    public:
        std::array< int, 3 > idx_left;
        std::array< int, 3 > idx_right;
        refine_scheme_2D() {
            idx_left[0] = 0; idx_left[1] = -1; idx_left[2] = 1;
            idx_right[0] = 2; idx_right[1] = -1; idx_right[2] = 1;
        }
    };
    
    template< class Grid, class ElementVector, class ... Manager >
    void refine2D( Grid &grid, ElementVector &marked, Manager & ... manager )
    {
        const int NV = Grid::NV;
                
        grid.numNodes = grid.coordinates.size();
        grid.numEdges = grid.edges.size();
        
        bool stillElementsToRefine = true;
        
        while( stillElementsToRefine ) {
            ManagedArray<int> markedElements( grid.nElements(), 0 );
            ManagedArray<int> markedEdges( grid.nEdges(), 0 );
            
            for(size_t el = 0; el < marked.size(); el++)
                markedElements[el] = marked[el];
            
            markEntities( grid, markedElements, markedEdges );
            
            int lmin = 0;
            for(size_t el = 0; el < grid.nElements(); el++)
                if( markedElements[el] > 0 ){
                    lmin = grid.level[el];
                    break;
                }
            for(size_t el = 0; el < grid.nElements(); el++)
                if( markedElements[el] > 0 && lmin > grid.level[el] )
                    lmin = grid.level[el];
                
            int nMarked = 0;
            for(size_t el = 0; el < grid.nElements(); el++)
                if( markedElements[el] > 0){
                    if( lmin < grid.level[el] )
                        markedElements[el] = 0;
                    else
                        nMarked++;
                }
            for(size_t el = 0; el < marked.size(); el++)    
               if( markedElements[el] > 0 && lmin == grid.level[el])
                   marked[el] = 0;    
       
            markedEdges.fill( 0 );
            for(size_t el = 0; el < grid.nElements(); el++)
                if( markedElements[el] > 0 )
                    markedEdges[grid.element2edges(el,1)] = 1;
            
            ManagedArray<int> newNodesOnElements( grid.nElements(), -1 );
            ManagedArray<int> newNodesOnEdges( grid.nEdges(), -1 );
            
            getNewNodes2D(grid,markedElements,markedEdges,newNodesOnElements,newNodesOnEdges);

            //ProlongationData_2D_pre_refine<ManagedArray<int> > pdata0(newNodesOnElements,newNodesOnEdges );
            //pre_refine_Manager( pdata0, manager ... );
            
            const size_t nOld = grid.nElements();
            grid.resizeElements( nOld + nMarked );
            
            for( size_t el = 0; el < nOld; el++)
                if( markedElements[el] > 0 )
                    grid.level[el] = lmin + 1;
            for(size_t k = nOld; k < nOld + nMarked; k++)
                grid.level[k] = lmin + 1;
            
            ManagedArray<int> newInnerEdges( grid.nEdges(), -1 );
            ManagedArray<int> left2rightElements( grid.nElements(), -1 );
            ManagedArray<int> left2rightEdges( grid.nEdges(), -1 );
               
            refineElements2D( grid, markedElements, newNodesOnElements, left2rightElements, nOld, lmin );
           
            getElement2Faces( grid.elements, grid.edges, grid.element2edges );
            
            //refineEntity2SubEntities2D( grid, markedElements, markedFaces, markedEdges, lmin, newInnerFaces, newInnerEdges, left2rightFaces, left2rightEdges );
            
            //ProlongationData_refine2D<ManagedArray<int> > pdata( newNodesOnElements, newNodesOnEdges, newInnerEdges, left2rightElements, left2rightEdges );
            
            //refine_Manager( pdata, manager ... );
            
            //refineLeftChild2D( grid, markedElements, lmin );

            stillElementsToRefine = false;
            for(size_t el = 0; el < marked.size(); el++)
                if(marked[el] > 0){
                    stillElementsToRefine = true;
                    break;
                }
        }
        
    }

    template< class Grid, class ElementVector, class EdgeVector >
    void markEntities( Grid &grid, ElementVector &markedElements, EdgeVector &markedEdges )
    {
        const int NV = Grid::NV;

        int nMarkedElements = 0;
        for(size_t el = 0; el < grid.nElements(); el++)
            if( markedElements[el] > 0 )
                nMarkedElements++;
        
        int stillElementsMarked = 1;
        int iterMarkEntities = 0;
        while( stillElementsMarked > 0 ) {
            for(size_t el = 0; el < grid.nElements(); el++)
                if( markedElements[el] > 0 ){
                    for(size_t k = 1; k < NV-1; k++)
                        markedEdges[grid.element2edges(el,k)] = 1;
                }
                
            getMarkedEntitiesFromMarkedEdges(grid,markedElements,markedEdges );

            const int nMarkedElements_old = nMarkedElements;
            nMarkedElements = 0;
            for(size_t el = 0; el < grid.nElements(); el++)
                if( markedElements[el] > 0 )
                    nMarkedElements++;
                
            stillElementsMarked = nMarkedElements - nMarkedElements_old;
        }
        for(size_t el = 0; el < grid.nElements(); el++)
            if( markedElements[el] > 0 )
                for(size_t k = 1; k < NV-1; k++)
                    markedEdges[grid.element2edges(el,k)] = 1;          
    }
    
    template< class Grid, class ElementVector, class EdgeVector >
    void getMarkedEntitiesFromMarkedEdges( Grid &grid, ElementVector &markedElements, EdgeVector &markedEdges ){
        const int NV = Grid::NV;
        
        ManagedArray<char, NV> element2markedEdges(grid.nElements(),0);
        
        for(size_t f = 0; f < grid.element2edges.size(); f++)
            for(size_t k = 0; k < NV; k++)
                element2markedEdges(f,k) = markedEdges[ grid.element2edges(f,k)];
                    
        for(size_t el = 0; el < grid.nElements(); el++)
            for(size_t k = 0; k < NV; k++)
                if(element2markedEdges(el,k) > 0){
                    markedElements[el] = 1;
                    break;
                }
    }
    
    template< class Grid, class ElementVector, class EdgeVector, class Vector >
    void getNewNodes2D( Grid &grid, ElementVector &markedElements, EdgeVector &markedEdges, Vector &newNodesOnElements, Vector &newNodesOnEdges )
    {
        const int DIM = Grid::DIM;
        const int NV = Grid::NV;
        
        int nC = grid.nNodes();
        int nNewNodes = 0;
        for(size_t i = 0; i < grid.nEdges(); i++)
            if( markedEdges[i] > 0 )
                newNodesOnEdges[i] = nC + nNewNodes++;
                            
       for(size_t i = 0; i < grid.nElements(); i++)
            if( markedElements[i] > 0 )
                newNodesOnElements[i] = newNodesOnEdges[grid.element2edges(i,1)];
            
       grid.resizeNodes( grid.numNodes + nNewNodes );
       
       std::vector<bool> nodeFlag( grid.nNodes() );
       std::fill( nodeFlag.begin(), nodeFlag.end(), false );
       for(size_t el = 0; el < markedElements.size(); el++)
           if( markedElements[el] > 0 ){
               const int y = newNodesOnElements[el];
               if( nodeFlag[ y ] == false ){
                    int e1 = grid.elements(el,0);
                    int e2 = grid.elements(el,NV-1);
                    for( int d = 0; d < DIM; ++d)
                        grid.coordinates(y,d) = .5*(grid.coordinates(e1,d) + grid.coordinates(e2,d)); 
                    nodeFlag[y] = true;
               }
           }
    }
       
    template< class Grid, class ElementVector, class NodeVector, class Vector >
    void refineElements2D( Grid &grid, ElementVector &marked, NodeVector &newNodesOnElements, Vector & left2rightElements,
                         int nOld, int lmin )
    {
        const int NV = Grid::NV;
                
        refine_scheme_2D scm;
        
        int counter = nOld;
        for(int el = 0; el < nOld; ++el){
            if(marked[el] > 0){
                grid.elements(counter,0) = grid.elements(el,scm.idx_right[0]);
                for(size_t k = 2; k < NV; ++k )
                    grid.elements(counter,k) = grid.elements(el,scm.idx_right[k]);
                grid.elements(counter,1) = newNodesOnElements[el];
                left2rightElements[el] = counter;
                counter++;
            }
        }
        for(int el = 0; el < nOld; ++el){
            if(marked[el] > 0){
                grid.elements(el,0) = grid.elements(el,scm.idx_left[0]);
                for(size_t k = NV-1; k > 1; --k )
                    grid.elements(el,k) = grid.elements(el,scm.idx_left[k]);
                grid.elements(el,1) = newNodesOnElements[el];
            }
        }
   
    }

} // namespace conformingsimplexgrid

#endif //REFINE_HH