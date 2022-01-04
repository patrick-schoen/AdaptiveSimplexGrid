#ifndef GLOBALNODEID_HH
#define GLOBALNODEID_HH

namespace conformingsimplexgrid {
    
template< class Grid, class Communicator, class Decomposition >
class GlobalNodeId
{
public:       
    
    typedef  GlobalNodeId<Grid, Communicator, Decomposition> ThisType; 
    const int NV = Grid::verticesPerElement;
    
    ManagedArray<int,2> globalNodeIndex;
    
    GlobalNodeId( Grid &grid, Communicator &comm, Decomposition &decomp)
    : grid_(grid), comm_(comm), decomp_(decomp)
    {}
    
    void communicate()
    {
        globalNodeIndex.resize( grid_.nNodes() );
        globalNodeIndex.fill(0);
        ManagedArray<int> nodes( grid_.nNodes() );
        
        for(int k = 0; k < (int)nodes.size(); k++)
            if( decomp_.node2ranks[k].size() && decomp_.node2ranks[k][0] < (int)comm_.rank() ) {
                globalNodeIndex(k,0) = decomp_.node2ranks[k][0];
                nodes[k] = 0;
            } else {
                globalNodeIndex(k,0) = (int)comm_.rank();
                nodes[k] = k;
            }        
  
        decomp_.sumNodeVector(nodes);
        
        for(int k = 0; k < (int)nodes.size(); k++)
            globalNodeIndex(k,1) = nodes[k];
    }
    
    //call on macor grid to initiate global Node-Index information
    template< class ElementVector >
    void init_macro( ElementVector &partition )
    {
        globalNodeIndex.resize(grid_.nNodes());
        globalNodeIndex.fill(-1);
        for( int k = 0; k < (int)grid_.nNodes(); k++)
            globalNodeIndex(k,1) = k;
        
        for( int el = 0; el < (int)grid_.nElements(); ++el){
            const int p = partition[el];
            for(int k = 0; k < NV; ++k)
                globalNodeIndex(grid_.elements(el,k),0) = p;
        }
        
        for( int el = 0; el < (int)grid_.nElements(); ++el){
            const int p = partition[el];
            for(int k = 0; k < NV; ++k )
                if(globalNodeIndex(grid_.elements(el,k),0) > p)
                    globalNodeIndex(grid_.elements(el,k),0) = p;
        }
    }
    
    void pack( void *buffer, std::size_t bufferSize, int &position )
    {            
        int nNodes = (int)globalNodeIndex.size();
        comm_.pack( &nNodes, 1, buffer, bufferSize, position);
        comm_.pack( globalNodeIndex.data(), 2*nNodes, buffer, bufferSize, position);
    }
    
    void unpack( void *buffer, std::size_t bufferSize, int &position )
    {
        int nNodes; 
        comm_.unpack( &nNodes, 1, buffer, bufferSize, position);
        globalNodeIndex.resize(nNodes);
        comm_.unpack( globalNodeIndex.data(), 2*nNodes, buffer, bufferSize, position);
    }

    std::size_t bufferSize()
    {
        return comm_.template bufferSize<int>(1+2*globalNodeIndex.size());
    }
    
    template< class NodeVector >
    void removeNodes( NodeVector & nodeVector, ThisType & subGlobalId) const
    {
        std::size_t ctr = 0;
        subGlobalId.globalNodeIndex.resize(globalNodeIndex.size()); 
        for( std::size_t k = 0; k < globalNodeIndex.size(); ++k){
            if( nodeVector[k] == false )
                continue;
            subGlobalId.globalNodeIndex(ctr,0) = globalNodeIndex(k,0);
            subGlobalId.globalNodeIndex(ctr,1) = globalNodeIndex(k,1);
            ctr++;
        }
        subGlobalId.globalNodeIndex.resize(ctr);
    }
    
    template< class NodeVector >
    void removeNodes( NodeVector & nodeVector )
    {
        std::size_t ctr = 0;
        for( std::size_t k = 0; k < globalNodeIndex.size(); ++k){
            if( nodeVector[k] == false )
                continue;
            globalNodeIndex(ctr,0) = globalNodeIndex(k,0);
            globalNodeIndex(ctr,1) = globalNodeIndex(k,1);
            ctr++;
        }
        globalNodeIndex.resize(ctr);
    }
    
    
private:
Grid &grid_;
Communicator &comm_;
Decomposition &decomp_;
};

template < class Grid, class Decomposition, class ... Manager  >
void sortByGlobalCoordinates( Grid &grid, Decomposition &decomp, Manager & ... manager )
{
    static const int dim = Grid::dimensionworld;
    ManagedArray<double,dim> mpF(grid.nFaces());
    ManagedArray<double,dim> mpE(grid.nEdges());
    
    for(int el = 0; el < (int)grid.nFaces(); el++)
        for( int k = 0; k < dim; k++)
            mpF(el,k) = (grid.coordinates(grid.faces(el,0),k)+grid.coordinates(grid.faces(el,1),k)+grid.coordinates(grid.faces(el,2),k) ) / 3.;
    
    for(int el = 0; el < (int)grid.nEdges(); el++)
        for( int k = 0; k < dim; k++)
            mpE(el,k) = (grid.coordinates(grid.edges(el,0),k)+grid.coordinates(grid.edges(el,1),k))  / 2.;    
    
    ManagedArray<int> idx, faces2newFaces, edges2newEdges, nodes2newNodes;
    sortByRow_void( mpF, idx, faces2newFaces );
    sortByRow_void( mpE, idx, edges2newEdges );
    sortByRow_void( grid.coordinates, idx, nodes2newNodes );
    
    grid.coordinates = rearange_inverse( grid.coordinates, nodes2newNodes );
    index2newIndex( grid.elements, nodes2newNodes );
    index2newIndex( grid.element2faces, faces2newFaces );
    grid.face2edges = rearange_inverse( grid.face2edges, faces2newFaces );
    getFacesForElements( grid.elements, grid.faces, grid.element2faces, grid.level);
    
    index2newIndex( grid.face2edges, edges2newEdges );
    getFacesForElements( grid.faces , grid.edges, grid.face2edges );
        
    grid.levelFaces = rearange_inverse( grid.levelFaces, faces2newFaces );
    decomp.node2ranks = rearange_inverse( decomp.node2ranks, nodes2newNodes );
    decomp.edge2ranks = rearange_inverse( decomp.edge2ranks, edges2newEdges );
    decomp.face2ranks = rearange_inverse( decomp.face2ranks, faces2newFaces );
    
    rearange_Manager( faces2newFaces, edges2newEdges, nodes2newNodes, manager ... );
    
    decomp.init_communication();
    
}

template < class Grid, class GlobalId, class Decomposition, class ... Manager >
void sortByGlobalIndices( Grid &grid, GlobalId &globalId, Decomposition &decomp, Manager & ... manager )
{
    ManagedArray<int> faces2newFaces, edges2newEdges, nodes2newNodes;
   
    faces2newFaces.resize(grid.nFaces());
    edges2newEdges.resize(grid.nEdges());
    nodes2newNodes.resize(grid.nNodes());
    
    ManagedArray<int> idx;
    globalId.globalNodeIndex = sortByRow( globalId.globalNodeIndex, idx, nodes2newNodes );
    index2newIndex( grid.elements, nodes2newNodes );
    index2newIndex( grid.faces, nodes2newNodes );
    index2newIndex( grid.edges, nodes2newNodes );
    
    sortRows( grid.faces );
    sortByRow_void( grid.faces, idx, faces2newFaces );
    index2newIndex( grid.element2faces, faces2newFaces );
    grid.face2edges = rearange_inverse( grid.face2edges, faces2newFaces );
    getFacesForElements( grid.elements, grid.faces, grid.element2faces, grid.level);
    
    sortRows( grid.edges );
    sortByRow_void( grid.edges, idx, edges2newEdges );
    index2newIndex( grid.face2edges, edges2newEdges );
    getFacesForElements( grid.faces , grid.edges, grid.face2edges );
 
    grid.coordinates = rearange_inverse( grid.coordinates, nodes2newNodes );
    
    grid.levelFaces = rearange_inverse( grid.levelFaces, faces2newFaces );
    decomp.node2ranks = rearange_inverse( decomp.node2ranks, nodes2newNodes );
    decomp.edge2ranks = rearange_inverse( decomp.edge2ranks, edges2newEdges );
    decomp.face2ranks = rearange_inverse( decomp.face2ranks, faces2newFaces );
    
    newIndices_Manager( faces2newFaces, edges2newEdges, nodes2newNodes, manager ... );
   
    decomp.init_communication(); 
}
} // namespace conformingsimplexgrid
#endif // GLOBALNODEID_HH
