/** 
@file coarse3D.hh
@brief implementation of the coarsening algorithm for a tetrahedral mesh. All marked elements, where the newest node can be removed from the mesh, are coarsed.
*/

#ifndef COARSE_3D_HH
#define COARSE_3D_HH

#include "../simplexGrid.hh"

namespace conformingsimplexgrid {

template< class Grid, class VectorBool, class Vector > 
struct ProlongationData_coarse
{ 
    ProlongationData_coarse( Grid & grid_in, VectorBool & remainingElements_in, Vector &face2newFaces_in,
        Vector &edge2newEdges_in, Vector &node2newNodes_in, Vector &left_in, Vector &right_in )
    : grid(grid_in), remainingElements(remainingElements_in), face2newFaces(face2newFaces_in),
    edge2newEdges(edge2newEdges_in), node2newNodes(node2newNodes_in), left(left_in), right(right_in)
    {}
    
    Grid &grid;
    VectorBool &remainingElements;
    Vector &face2newFaces;
    Vector &edge2newEdges;
    Vector &node2newNodes;
    Vector &left;
    Vector &right;
};
          
    
    template< class Grid, class ElementVector, class ... Manager >
    void coarse3D( Grid &grid, ElementVector &marked, Manager & ... manager )
    {  
        coarse(grid, marked, manager ... );
    }
        
    template< class Grid, class ElementVector, class ... Manager >
    void coarse( Grid &grid, ElementVector &marked, Manager & ... manager )
    {  
        grid.numNodes = grid.coordinates.size();
        grid.numEdges = grid.edges.size();
        grid.numFaces = grid.faces.size();
        
        ManagedArray<int> markedNodes( grid.nNodes(), 0 );
        findRemovableNodes( grid, marked, markedNodes, manager ... );
        marked.fill(0);
        for(size_t el = 0; el < marked.size(); ++el)
            if( markedNodes[ grid.elements(el,1) ] == 1 )
                marked[el] = 1; 
        coarseNodes( grid, marked, manager ... );
    }
    
    template< class Decomposition, class Grid, class ElementVector, class ... Manager >
    void coarse_parallel( Decomposition &decomp, Grid &grid, ElementVector &marked, Manager & ... manager )
    {  
        const int NV = Grid::NV;
        grid.numNodes = grid.coordinates.size();
        grid.numEdges = grid.edges.size();
        grid.numFaces = grid.faces.size();
        
        ManagedArray<int> markedNodes( grid.nNodes(), 0 );
        
        for(size_t el = 0; el < marked.size(); ++el)
            if( marked[el] > 0 && grid.level[el] > 0 )
                markedNodes[ grid.elements(el,1) ] = 1;
                
        for(size_t el = 0; el < grid.nElements(); ++el){
            markedNodes[ grid.elements(el,0) ] = 0;
            for( int k = 2; k < NV; ++k)
                markedNodes[ grid.elements(el,k) ] = 0;
        }
        
        
        //findRemovableNodes( grid, marked, markedNodes );
        //ManagedArray<int> notRemovable( grid.nNodes(), 1 );
        //for(size_t k = 0; k < markedNodes.size(); k++)
        //    if( markedNodes[k] > 0 )
        //        notRemovable[k] = 0; 
                        
        unmarkedNodesWhereBrotherIsOnOtherRank( grid, markedNodes );
        
        //for(size_t k = 0; k < markedNodes.size(); k++)
        //    if( notRemovable[k] > 0 )
        //        markedNodes[k] = 0;
   
        decomp.sumNodeVector( markedNodes );
        for(size_t k = 0; k < markedNodes.size(); k++){
            if( (int)decomp.node2ranks[k].size() == markedNodes[k] )
                markedNodes[k] = 1;
            else
                markedNodes[k] = 0;
        }
        
        marked.fill(0);
        for(size_t el = 0; el < marked.size(); ++el)
            if( markedNodes[ grid.elements(el,1) ] == 1 )
                marked[el] = 1;
            
        coarseNodes( grid, marked, decomp, manager ... );
    }
        
     
    template< class Grid, class ElementVector, class NodeVector >
    void findRemovableNodes( Grid &grid, ElementVector &marked, NodeVector &markedNodes)
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
    
  
    template< class Grid, class Communicator, class Decomposition, class GlobalId, class Vector, class ... Manager >
    void repartitionForCoarsening( Grid &grid, Communicator &comm, Decomposition &decomp, GlobalId &globalId, Vector &markedElements, Manager & ... manager )
    {
        const int NV = Grid::NV;
        assert( grid.elements.size() == markedElements.size());
        assert( grid.elements.size() == grid.level.size() );
        ManagedArray<char> markedNodes( grid.nNodes(), 0);
        for(size_t el = 0; el < markedElements.size(); ++el)
            if( markedElements[el] > 0 && grid.level[el] > 0 )
                markedNodes[ grid.elements(el,1) ] = 1;
        for(size_t el = 0; el < grid.nElements(); ++el){
            markedNodes[ grid.elements(el,0) ] = 0;
            for( int k = 2; k < NV; ++k)
                markedNodes[ grid.elements(el,k) ] = 0;
        }           
        //findRemovableNodes( grid, markedElements, markedNodes);
        decomp.sumNodeVector( markedNodes );
        for(size_t k = 0; k < markedNodes.size(); k++){
            if( (int)decomp.node2ranks[k].size() == markedNodes[k] )
                markedNodes[k] = 1;
            else
                markedNodes[k] = 0;
        }
        ManagedArray<int> partition(grid.nElements());
        getRightChildElements( grid, decomp, markedNodes, partition );
        std::vector<bool> remainingElements(grid.nElements(), false);   
        for( size_t el = 0; el < grid.nElements(); ++el)
            if( partition[el] == comm.rank() )
                remainingElements[el] = true;
        remove_false( markedElements , remainingElements );
        getRightChildElements( grid, decomp, markedNodes, partition );
        
        redistributeNewPartition( grid, comm, decomp, globalId, partition, manager ... );
        
        const int N = markedElements.size();
        markedElements.resize(grid.elements.size());
        for(int k = N; k < grid.elements.size(); k++)
            markedElements[k] = 1;
    }
      
    template< class Grid, class Decomposition, class NodeVector, class ElementVector >
    void getRightChildElements( Grid &grid, Decomposition &decomp, NodeVector &markedNodes, ElementVector &partition )
    {
        partition.fill(decomp.rank());
        
        const int NV = Grid::NV;
        ManagedArray<int> faceData( grid.nFaces(), 0 );
        for(size_t el = 0; el < grid.nElements(); ++el)
            if( markedNodes[ grid.elements(el,1) ] == 1 ){
                const int f0 = grid.element2faces(el,0);
                const int f1 = grid.element2faces(el,NV-1);
                int f;
                if( grid.levelFaces[f1] > grid.levelFaces[f0]){
                    f = f1;
                } else {
                    f = f0;
                }
                faceData[f]++;
            }
       for(size_t el = 0; el < grid.nElements(); ++el)
            if( markedNodes[ grid.elements(el,1) ] == 1 ){
               const int f0 = grid.element2faces(el,0);
               const int f1 = grid.element2faces(el,NV-1);
               int f;
               int k;
               if( grid.levelFaces[f1] > grid.levelFaces[f0]){
                   f = f1;
                   k = 3;
               } else {
                   f = f0;
                   k = 0;
               }
               if( faceData[f] < 2 ){
                   const bool isRightChild = (bool)grid.rightChild(el,k);
                   if( isRightChild ){
                       int otherRank = (decomp.face2ranks[f][0] == decomp.rank()) ? decomp.face2ranks[f][1] : decomp.face2ranks[f][0];  
                       partition[el] = otherRank;
                   }
               }
            }
    }
    
    template< class Grid, class NodeVector >
    void unmarkedNodesWhereBrotherIsOnOtherRank( Grid &grid, NodeVector &markedNodes )
    {
        const int NV = Grid::NV;
        ManagedArray<int> faceData( grid.nFaces(), 0 );
        for(size_t el = 0; el < grid.nElements(); ++el)
            if( markedNodes[ grid.elements(el,1) ] == 1 ){
               const int f0 = grid.element2faces(el,0);
               const int f1 = grid.element2faces(el,NV-1);
               int f;
               if( grid.levelFaces[f1] > grid.levelFaces[f0]){
                   f = f1;
               } else {
                   f = f0;
               }
               faceData[f]++;
            }
       for(size_t el = 0; el < grid.nElements(); ++el)
            if( markedNodes[ grid.elements(el,1) ] == 1 ){
               const int f0 = grid.element2faces(el,0);
               const int f1 = grid.element2faces(el,NV-1);
               int f;
               if( grid.levelFaces[f1] > grid.levelFaces[f0]){
                   f = f1;
               } else {
                   f = f0;
               }
               if( faceData[f] < 2)
                   markedNodes[ grid.elements(el,1) ] = 0;
            }    
    }
    
    template< class Grid, class ElementVector, class Vector >
    void getBrothers( Grid &grid, ElementVector &marked, Vector &left, Vector &right )
    {     
        const int NV = Grid::NV;
        
        int nMarked = 0;
        for(size_t el = 0; el < marked.size(); ++el)
            if( marked[el] > 0 )
                nMarked++;
        
        left.resize( nMarked / 2 ); left.fill(0);
        right.resize( nMarked / 2 ); right.fill(0);
        ManagedArray<int> newestFacesData( grid.nFaces(), -1);
        for(size_t el = 0, counter = 0; el < marked.size(); ++el){
            if( marked[el] > 0 ){
                const int f0 = grid.element2faces(el,0);
                const int f1 = grid.element2faces(el,NV-1);
                int f;
                int k;
                if( grid.levelFaces[f1] > grid.levelFaces[f0]){
                    f = f1;
                    k = 3;
                } else {
                    f = f0;
                    k = 0;
                }
                if( newestFacesData[f] < 0 ){
                    newestFacesData[f] = el;
                } else {
                    const bool isRightChild = (bool)grid.rightChild(el,k);
                    const int t1 = newestFacesData[f];
                    const int t2 = el;
                    if( isRightChild ){
                        left[counter] = t1;
                        right[counter] = t2;
                    } else {
                        left[counter] = t2;
                        right[counter] = t1;
                    }
                    counter++;
                }
            }
        } 
    }
    
    template< class Grid, class ElementVector, class ... Manager >
    void coarseNodes( Grid &grid, ElementVector &marked, Manager & ... manager )
    {        
        const int DIM = Grid::DIM;

        ManagedArray<int> left, right;
        getBrothers( grid, marked, left, right );
        
        std::vector<bool> remainingElements( grid.nElements(), true );
        ManagedArray<int> face2newFaces( grid.nFaces(), 0 );
        ManagedArray<int> edge2newEdges( grid.nEdges(), 0 );
        ManagedArray<int> node2newNodes( grid.nNodes(), 0 );
        
        for( size_t k = 0; k < right.size(); k++)
            remainingElements[right[k]] = false;

        coarse_nodes2newNodes(grid,left,right,remainingElements,node2newNodes);
        coarse_elements(grid,node2newNodes);
        coarse_entity2subEntities(grid,left,right,remainingElements,face2newFaces,edge2newEdges);
        coarse_rightChild(grid,left,right,remainingElements);
        coarse_level(grid,left,right,remainingElements,face2newFaces);
        
        getFacesForElements(grid.elements,grid.faces,grid.element2faces,grid.level); 
        getFacesForElements(grid.faces,grid.edges,grid.face2edges);
        
        ProlongationData_coarse<Grid, std::vector<bool>, ManagedArray<int> > pdata( grid, remainingElements,
            face2newFaces, edge2newEdges, node2newNodes, left, right );
        
        coarse_Manager( pdata, manager ... );
        
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
    void coarse_nodes2newNodes( Grid &grid, Vector &left, Vector &right, BoolVector &remainingElements, NodeVector &nodes2newNodes ){
        const int NV = Grid::NV;
        
        ManagedArray<int> dataOnFace( grid.nFaces(), 0);
        for(size_t k = 0; k < right.size(); ++k){
            dataOnFace[ grid.element2faces(right[k],0) ] = grid.elements(right[k],3);
            dataOnFace[ grid.element2faces(right[k],3) ] = grid.elements(right[k],0);
        }
        for(size_t k = 0; k < left.size(); ++k){
            grid.elements(left[k],1) = grid.elements(left[k],2);
            grid.elements(left[k],2) = grid.elements(left[k],3);
            grid.elements(left[k],3) = dataOnFace[ grid.element2faces(left[k],3)];              
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
    }
        
    template< class Grid, class NodeVector >
    void coarse_elements( Grid &grid, NodeVector &nodes2newNodes ){
        const int NV = Grid::NV;   
        
        for(size_t el = 0; el < grid.elements.size(); ++el)
            for( size_t k = 0; k < NV; ++k){
                assert( nodes2newNodes[grid.elements(el,k)] >= 0 );
                grid.elements(el,k) = nodes2newNodes[grid.elements(el,k)];
            }
    }
    
    template< class Grid, class Vector, class BoolVector >
    void coarse_rightChild( Grid &grid, Vector &left, Vector &right, BoolVector &remainingElements ){
        const int NV = Grid::NV;
        
        for(size_t k = 0; k < left.size(); ++k){
            const int tmp0 = grid.rightChild(left[k],0);
            const int tmp1 = grid.rightChild(left[k],1);
            const int tmp2 = grid.rightChild(left[k],2);
            grid.rightChild(left[k],0) =  tmp2;
            grid.rightChild(left[k],1) =  tmp0;
            grid.rightChild(left[k],2) =  tmp1;
            grid.rightChild(left[k],3) =  grid.rightChild(right[k],2);
        }
        
        int counter = 0;
        for(size_t el = 0; el < grid.rightChild.size(); ++el)
            if( remainingElements[el] ){
                for( size_t k = 0; k < NV; ++k)
                    grid.rightChild(counter,k) = grid.rightChild(el,k);
                counter++;
            }
        grid.rightChild.resize(counter);
        
    }
    
    template< class Grid, class Vector, class BoolVector, class FaceVector >
    void coarse_level( Grid &grid, Vector &left, Vector &, BoolVector &remainingElements, FaceVector &face2newFaces ){
        for( size_t k = 0; k < left.size(); ++k)
            grid.level[left[k]] = grid.level[left[k]] - 1;
        
        int counter = 0;
        for(size_t k = 0; k < grid.level.size(); ++k)
            if( remainingElements[k] )
                grid.level[counter++] = grid.level[k];
        grid.level.resize(counter);
                
        counter = 0;
        for(size_t k = 0; k < grid.levelFaces.size(); ++k){
            if( face2newFaces[k] < 0 )
                continue;
            grid.levelFaces[counter++] = grid.levelFaces[k];
        }
        grid.levelFaces.resize(counter);
    }
    
    template< class Grid, class Vector, class BoolVector, class FaceVector, class EdgeVector >
    void coarse_entity2subEntities( Grid &grid, Vector &left, Vector &right, 
                                    BoolVector &remainingElements, FaceVector &face2newFaces, EdgeVector &edge2newEdges  ){
        const int NV = Grid::NV;
        int counter = 0;
        assert(NV == 4);
        
        face2newFaces.fill(0);

        //mark inner faces with -2
        for(size_t k = 0; k < left.size(); ++k){
            face2newFaces[grid.element2faces(left[k],3)] = -2; 
        }
                
        //mark right children faces with -1 
        for(size_t k = 0; k < right.size(); ++k){
            const int type = (grid.level[right[k]] % (NV-1));
            if( type == 0 ) {
                face2newFaces[grid.element2faces(right[k],0)] = -1; 
                face2newFaces[grid.element2faces(right[k],1)] = -1; 
            } else {
                face2newFaces[grid.element2faces(right[k],1)] = -1; 
                face2newFaces[grid.element2faces(right[k],3)] = -1; 
            }
        }
        
        //mark left children faces with +1
        for(size_t k = 0; k < left.size(); ++k){
            face2newFaces[ grid.element2faces(left[k],0) ] = 1;
            face2newFaces[ grid.element2faces(left[k],1) ] = 1;
        }

        for(size_t k = 0; k < right.size(); ++k){
            const int tmp0 = grid.element2faces(left[k],0);
            const int tmp1 = grid.element2faces(left[k],1);
            grid.element2faces(left[k],0) = grid.element2faces(left[k],2);
            grid.element2faces(left[k],1) = tmp0;
            grid.element2faces(left[k],2) = tmp1;
            grid.element2faces(left[k],3) = grid.element2faces(right[k],2);  
        }
        
        ManagedArray<int> dataOnEdge( grid.nEdges(), 0);
        for(size_t k = 0; k < face2newFaces.size(); ++k){
            if( face2newFaces[k] == -1 ){  
                dataOnEdge[ grid.face2edges(k,0) ] = grid.face2edges(k,1);
                dataOnEdge[ grid.face2edges(k,2) ] = grid.face2edges(k,1);       
            }
        }
        
        for(size_t k = 0; k < face2newFaces.size(); ++k){
            if( face2newFaces[k] == 1 ){
                const int tmp0 = grid.face2edges(k,0);
                grid.face2edges(k,0) = grid.face2edges(k,1);
                grid.face2edges(k,1) = tmp0;
                grid.face2edges(k,2) = dataOnEdge[ grid.face2edges(k,2) ];    
            }
        }
        
        counter = 0;
        for(size_t el = 0; el < grid.element2faces.size(); ++el)
            if( remainingElements[el] ){
                for( size_t k = 0; k < NV; ++k)
                    grid.element2faces(counter,k) = grid.element2faces(el,k);
                counter++;
            }
        grid.element2faces.resize(counter);
                
        face2newFaces.fill(-1);
        for(size_t k = 0; k < grid.element2faces.size(); ++k)
            for( size_t j = 0; j < NV; ++j)
                face2newFaces[grid.element2faces(k,j)] = 1;

        for(size_t k = 0, counter = 0; k < face2newFaces.size(); ++k){
            if( face2newFaces[k] < 0 )
                continue;
            face2newFaces[k] = counter++;
        }
        
        for(size_t el = 0, counter = 0; el < grid.element2faces.size(); ++el){
            for( size_t k = 0; k < NV; ++k){
                grid.element2faces(counter,k) = face2newFaces[ grid.element2faces(el,k) ];
            }
            counter++;
        }
            
        counter = 0;
        for(size_t el = 0; el < grid.face2edges.size(); ++el){
            if( face2newFaces[el] < 0 )
                continue;
            for( size_t k = 0; k < NV - 1; ++k)
                grid.face2edges(counter,k) = grid.face2edges(el,k);
            counter++;
        }
        grid.face2edges.resize(counter);
        
        edge2newEdges.fill(-1);
        for(size_t k = 0; k < grid.face2edges.size(); ++k)
            for( size_t j = 0; j < NV-1; ++j)
                edge2newEdges[grid.face2edges(k,j)] = 1;
  
        for(size_t k = 0, counter = 0; k < edge2newEdges.size(); ++k){
            if( edge2newEdges[k] < 0 )
                continue;
            edge2newEdges[k] = counter++;
        }

        for(size_t el = 0; el < grid.face2edges.size(); ++el){
            for( size_t k = 0; k < NV - 1; ++k)
                grid.face2edges(el,k) = edge2newEdges[ grid.face2edges(el,k) ];
        }
    }
}

#endif //COARSE_HH