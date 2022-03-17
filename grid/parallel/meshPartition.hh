#ifndef MESH_PARTITION_HH
#define MESH_PARTITION_HH

#include "../simplexGrid.hh"
#include "../parallel.hh"


namespace conformingsimplexgrid{

    class DefaultContainer;     
    template< class Communicator > struct SentMesh;

    template < class Grid, class Communicator, class Decomposition, class GlobalId, class ElementVector >
    void redistributeNewPartition( Grid & grid, Communicator & comm, Decomposition & decomp, GlobalId & globalId, ElementVector &partition )
    {
        DefaultContainer manager;
        redistributeNewPartition( grid, comm, decomp, globalId, partition, manager );
    }
    
    template < class Grid, class Communicator, class Decomposition, class GlobalId, class ElementVector, class Manager >
    void redistributeNewPartition( Grid & grid, Communicator & comm, Decomposition & decomp, GlobalId & globalId, ElementVector &partition, Manager & manager )
    {       
        globalId.communicate();
        decomp.init_partition( partition );
        communicateDecompositionData( decomp, comm );
                    
        ManagedArray<int> nElementsPerRank_send( comm.size(), 0 );
        for( size_t k = 0; k < partition.size(); k++ )
            nElementsPerRank_send[ partition[k] ]++;
        
        ManagedArray< char > SendMessageToRank(comm.size(),0);  
        ManagedArray< char > RecvMessageFromRank(comm.size(),0); 
        std::vector< int > processes;
        
        for( int p = 0; p < (int)comm.size(); p++ )
            if( nElementsPerRank_send[p] > 0 && p != (int)comm.rank() ){
                processes.push_back(p);
                SendMessageToRank[p] = 1;
            }
        
        comm.all2all( SendMessageToRank.data(), RecvMessageFromRank.data(), 1 );
 
        std::size_t recvMessages = 0;
        for( int p = 0; p < (int)comm.size(); p++ )
            if( RecvMessageFromRank[p] > 0 )
                recvMessages++;
        
        std::vector< SentMesh<Communicator> > SentMeshContainer;
        for( auto p : processes ) {
            Grid subgrid;
            Decomposition subDecomposition( subgrid, comm );
            GlobalId subGlobalId( subgrid, comm, subDecomposition );
            getSubGrid( subgrid, subDecomposition, subGlobalId, p, grid, decomp, globalId, partition, manager );
            SentMeshContainer.push_back( sendMesh( subgrid, p, 1, comm, subGlobalId, subDecomposition, manager ) );
        }
 
        removeEntitiesfromGrid( grid, decomp, globalId, partition, comm.rank(), manager );

        while( recvMessages > 0 ){
            Grid subgrid;
            Decomposition subDecomposition( subgrid, comm );
            GlobalId subGlobalId( subgrid, comm, subDecomposition );
            if( recvMesh( comm, subgrid, subGlobalId, subDecomposition, manager ) ) {
                conjoinTwoGrids( grid, decomp, globalId, subgrid, subDecomposition, subGlobalId, manager );
                recvMessages--;
            }
        }
        if(grid.nElements())
            removeMultipleIndices( grid, decomp, globalId, manager );
        decomp.init_communication();
        decomp.barrier();
    }
    
    template < class Grid, class Communicator, class Decomposition, class GlobalId, class ElementVector >
    void redistributeNewPartition_Timer( Grid & grid, Communicator & comm, Decomposition & decomp, GlobalId & globalId, ElementVector &partition, double &time_data, double &time_communication )
    {
        DefaultContainer manager;
        redistributeNewPartition_Timer( grid, comm, decomp, globalId, partition, manager, time_data, time_communication );
    }
    
    template < class Grid, class Communicator, class Decomposition, class GlobalId, class ElementVector, class Manager >
    void redistributeNewPartition_Timer( Grid & grid, Communicator & comm, Decomposition & decomp, GlobalId & globalId, ElementVector &partition, Manager & manager, double &time_data, double &time_communication )
    {      
       
        time_data = 0.;
        time_communication = 0.;
        double start_time, stop_time;
        
        start_time = clock();
        globalId.communicate(); 
        stop_time = (float)(clock() - start_time) / CLOCKS_PER_SEC;
        time_communication += stop_time; 
            
        start_time = clock();
        decomp.init_partition( partition ); 
        stop_time = (float)(clock() - start_time) / CLOCKS_PER_SEC;
        time_data += stop_time;
        
        start_time = clock();
        communicateDecompositionData( decomp, comm ); 
       
        ManagedArray<int> nElementsPerRank_send( comm.size(), 0 );
        for( size_t k = 0; k < partition.size(); k++ )
            nElementsPerRank_send[ partition[k] ]++;
        
        ManagedArray< char > SendMessageToRank(comm.size(),0);  
        ManagedArray< char > RecvMessageFromRank(comm.size(),0); 
        std::vector< int > processes;
        
        for( int p = 0; p < (int)comm.size(); p++ )
            if( nElementsPerRank_send[p] > 0 && p != (int)comm.rank() ){
                processes.push_back(p);
                SendMessageToRank[p] = 1;
            }
        
        comm.all2all( SendMessageToRank.data(), RecvMessageFromRank.data(), 1 );
        stop_time = (float)(clock() - start_time) / CLOCKS_PER_SEC;
        time_communication += stop_time;      
        
        std::size_t recvMessages = 0;
        for( int p = 0; p < (int)comm.size(); p++ )
            if( RecvMessageFromRank[p] > 0 )
                recvMessages++;
        
        std::vector< SentMesh<Communicator> > SentMeshContainer;
        for( auto p : processes ) {
            Grid subgrid;
            Decomposition subDecomposition( subgrid, comm );
            GlobalId subGlobalId( subgrid, comm, subDecomposition );
            
            start_time = clock();
            getSubGrid( subgrid, subDecomposition, subGlobalId, p, grid, decomp, globalId, partition, manager );
            stop_time = (float)(clock() - start_time) / CLOCKS_PER_SEC;
            time_data += stop_time;
            
            start_time = clock();
            SentMeshContainer.push_back( sendMesh( subgrid, p, 1, comm, subGlobalId, subDecomposition, manager ) );
            stop_time = (float)(clock() - start_time) / CLOCKS_PER_SEC;
            time_communication += stop_time;   
        }
        
        start_time = clock();
        removeEntitiesfromGrid( grid, decomp, globalId, partition, comm.rank(), manager );
        stop_time = (float)(clock() - start_time) / CLOCKS_PER_SEC;
        time_data += stop_time;
        
        start_time = clock();
        while( recvMessages > 0 ){
            Grid subgrid;
            Decomposition subDecomposition( subgrid, comm );
            GlobalId subGlobalId( subgrid, comm, subDecomposition );
            start_time = clock();
            bool flag = recvMesh( comm, subgrid, subGlobalId, subDecomposition, manager );
            stop_time = (float)(clock() - start_time) / CLOCKS_PER_SEC;
            time_communication += stop_time;   
            if( flag ) {
                start_time = clock();
                conjoinTwoGrids( grid, decomp, globalId, subgrid, subDecomposition, subGlobalId, manager );
                stop_time = (float)(clock() - start_time) / CLOCKS_PER_SEC;
                time_data += stop_time;
                recvMessages--;
            }
        }
        
        start_time = clock();
        removeMultipleIndices( grid, decomp, globalId, manager );
        decomp.init_communication();
        stop_time = (float)(clock() - start_time) / CLOCKS_PER_SEC;
        time_data += stop_time;
        
        decomp.barrier();
    }
    
    
    template < class Subgrid, class SubDecomposition, class SubGlobalId, class Grid, class Decomposition, class GlobalId, class ElementVector  >
    void getSubGrid( Subgrid & subgrid, SubDecomposition & subDecomposition, SubGlobalId & subGlobalId, int rank, 
                     const Grid &grid, const Decomposition &decomposition, const GlobalId & globalId, const ElementVector &partition )
    {
        DefaultContainer manager;
        getSubGrid(subgrid,subDecomposition,subGlobalId,rank,grid,decomposition,globalId,partition,manager);
    }
    
    
    template < class Subgrid, class SubDecomposition, class SubGlobalId, class Grid, class Decomposition, class GlobalId, class ElementVector, class Manager  >
    void getSubGrid( Subgrid & subgrid, SubDecomposition & subDecomposition, SubGlobalId & subGlobalId, int rank, 
                     const Grid &grid, const Decomposition &decomposition, const GlobalId & globalId, const ElementVector &partition,
                     Manager & manager )
    {
        const int NV = Grid::verticesPerElement;
        std::vector<bool> remainingElements(grid.nElements(), false);

        for( size_t el = 0; el < grid.nElements(); ++el)
            if( partition[el] == rank )
                remainingElements[el] = true;
            
        remove_false( grid.elements, remainingElements, subgrid.elements);        
        remove_false( grid.rightChild, remainingElements, subgrid.rightChild);   
        remove_false( grid.element2faces, remainingElements, subgrid.element2faces);   
        remove_false( grid.level, remainingElements, subgrid.level);
        
#if 0
        int ctr_size = 0;
        for( size_t el = 0; el < grid.nElements(); ++el )
            if( partition[el] == rank )
                ctr_size++;
  
        subgrid.resizeElements( ctr_size );
        size_t ctr_el = 0;
        for( size_t el = 0; el < grid.nElements(); ++el)
            if( partition[el] == rank ){
                for( size_t j = 0; j < NV; ++j)
                    subgrid.elements(ctr_el,j) = grid.elements(el,j);
                for( size_t j = 0; j < NV; ++j)
                    subgrid.rightChild(ctr_el,j) = grid.rightChild(el,j);
                for( size_t j = 0; j < NV; ++j)
                    subgrid.element2faces(ctr_el,j) = grid.element2faces(el,j);
                subgrid.level[ctr_el] = grid.level[el];
                ctr_el++;
            }
#endif

        ManagedArray<int> nodes2newNodes( grid.nNodes(), -1 );
        for( size_t el = 0; el < subgrid.nElements(); ++el)
            for( size_t j = 0; j < NV; ++j)
                nodes2newNodes[ subgrid.elements(el,j) ] = 1;
  
        int nNewNodes = indexVector_enumerate( nodes2newNodes );
        index2newIndex( subgrid.elements, nodes2newNodes );
        
        remove_negative( grid.coordinates, nodes2newNodes, subgrid.coordinates );

        ManagedArray<int> faces2newFaces( grid.nFaces(), -1 );
        for( size_t el = 0; el < grid.nElements(); ++el)
            if( partition[el] == rank )
                for( size_t k = 0; k < NV; ++k )
                    faces2newFaces[ grid.element2faces(el,k) ] = 1;

        int nNewFaces = indexVector_enumerate( faces2newFaces );
            
        ManagedArray<int> edges2newEdges( grid.nEdges(), -1 );
        for( size_t el = 0; el < grid.nFaces(); ++el)
            if( faces2newFaces[el] >= 0 )
                for( size_t k = 0; k < NV-1; ++k )
                    edges2newEdges[ grid.face2edges(el,k) ] = 1;  

        int nNewEdges = indexVector_enumerate( edges2newEdges );
        
        subgrid.resizeFaces(nNewFaces);
 
        for( size_t el = 0, ctr = 0; el < grid.nFaces(); ++el)
            if( faces2newFaces[el] >= 0 ){
                for( size_t k = 0; k < NV-1; ++k )
                    subgrid.face2edges(ctr,k) = grid.face2edges(el,k);
                subgrid.levelFaces[ctr] = grid.levelFaces[el];
                ctr++;
            }
            
        index2newIndex( subgrid.element2faces, faces2newFaces );
        index2newIndex( subgrid.face2edges, edges2newEdges );
        
        remove_negative( globalId.globalNodeIndex, nodes2newNodes, subGlobalId.globalNodeIndex );
        remove_negative( decomposition.face2ranks, faces2newFaces, subDecomposition.face2ranks );
        remove_negative( decomposition.edge2ranks, edges2newEdges, subDecomposition.edge2ranks );
        remove_negative( decomposition.node2ranks, nodes2newNodes, subDecomposition.node2ranks );
        
        manager.preparePacking( remainingElements, faces2newFaces, edges2newEdges, nodes2newNodes );
        
        assert( subgrid.element2faces.min() == 0 );
        assert( subgrid.face2edges.min() == 0 );
    }

    template < class Grid, class Decomposition, class GlobalId, class ElementVector >
    void removeEntitiesfromGrid( Grid & grid, Decomposition & decomposition, GlobalId & globalId, ElementVector & partition, int rank )
    {
        DefaultContainer manager;
        removeEntitiesfromGrid(grid,decomposition,globalId,partition,rank,manager);
    }
    
    template < class Grid, class Decomposition, class GlobalId, class ElementVector, class Manager >
    void removeEntitiesfromGrid( Grid & grid, Decomposition & decomposition, GlobalId & globalId, ElementVector & partition, int rank, Manager & manager )
    {
        const int NV = Grid::verticesPerElement;
        
        std::vector<bool> remainingElements(grid.nElements(), false);   
        ManagedArray<int> faces2newFaces( grid.nFaces(), -1 );
        ManagedArray<int> edges2newEdges( grid.nEdges(), -1 ); 
        ManagedArray<int> nodes2newNodes( grid.nNodes(), -1 );
        
        for( size_t el = 0; el < grid.nElements(); ++el)
            if( partition[el] == rank )
                remainingElements[el] = true;
            
        remove_false( grid.elements , remainingElements );        
        remove_false( grid.rightChild , remainingElements );   
        remove_false( grid.element2faces , remainingElements );   
        remove_false( grid.level , remainingElements );

        for( size_t el = 0; el < grid.element2faces.size(); ++el)
            for( size_t k = 0; k < NV; ++k )
                faces2newFaces[ grid.element2faces(el,k) ] = 1;
        
        remove_negative( grid.face2edges, faces2newFaces );    
            
        for( size_t el = 0; el < grid.face2edges.size(); ++el)
            for( size_t k = 0; k < NV-1; ++k )
                edges2newEdges[ grid.face2edges(el,k) ] = 1;
               
        for( size_t el = 0; el < grid.nElements(); ++el)
            for( size_t j = 0; j < NV; ++j)
                nodes2newNodes[ grid.elements(el,j) ] = 1;      
  
        indexVector_enumerate( faces2newFaces );        
        indexVector_enumerate( edges2newEdges );
        indexVector_enumerate( nodes2newNodes );

        remove_negative( grid.coordinates, nodes2newNodes );
        
        index2newIndex( grid.element2faces, faces2newFaces );
        index2newIndex( grid.face2edges, edges2newEdges );
        index2newIndex( grid.elements, nodes2newNodes );
                
        remove_negative( grid.levelFaces, faces2newFaces );
        remove_negative( globalId.globalNodeIndex, nodes2newNodes );
        remove_negative( decomposition.face2ranks, faces2newFaces );
        remove_negative( decomposition.edge2ranks, edges2newEdges );
        remove_negative( decomposition.node2ranks, nodes2newNodes );
        
        ManagedArray<int> left,right;
        ProlongationData_coarse<Grid,std::vector<bool>,ManagedArray<int>> pData( grid, remainingElements, faces2newFaces, edges2newEdges, nodes2newNodes, left, right );
        coarse_Manager( pData, manager );
               
        if( grid.elements.size() > 0 ) {
            getFacesForElements( grid.elements, grid.faces, grid.element2faces, grid.level );
            getFacesForElements( grid.faces, grid.edges, grid.face2edges );
        } else {
            grid.resizeElements( 0 ); 
            grid.resizeFaces( 0 );
            grid.resizeEdges( 0 );
        }
    }

    template < class Ranks,  class Communicator, class Grid, class Decomposition, class GlobalId>
    void sendSubDomains(Ranks & ranks, Communicator & comm, std::vector<Grid> & subGrids,
                        std::vector<Decomposition> & subDecomposition, std::vector<GlobalId> & subGlobalId )
    {
        for( int k = 0; k < ranks.size(); k++ ) {
            const int p = ranks[k];
            auto sentMesh = sendMesh( subGrids[k], p, 1, comm, subDecomposition[k], subGlobalId[k] );
        }
    }
    
    template < class Ranks,  class Communicator, class Grid, class Decomposition, class GlobalId>
    void recvSubDomains(Ranks & ranks, Communicator & comm, std::vector<Grid> & subGrids,
                        std::vector<Decomposition> & subDecomposition, std::vector<GlobalId> & subGlobalId )
    {
        for( int k = 0; k < ranks.size(); k++ ) {
            const int p = ranks[k];
            subGrids[k] = recvMesh( p, 1, comm, subDecomposition[k], subGlobalId[k]  );
        }
    }
    
    template < class Grid, class Decomposition, class GlobalId, class Manager >
    void conjoinTwoGrids( Grid &grid1, Decomposition &decomposition1, GlobalId & globalID1,
                          Grid &grid2, Decomposition &decomposition2, GlobalId & globalID2, Manager & manager )
    {
        grid2.elements += grid1.nNodes();
        grid2.faces += grid1.nNodes();
        grid2.edges += grid1.nNodes();
        grid2.element2faces += grid1.nFaces();
        grid2.face2edges += grid1.nEdges();
        
        
        grid1.elements.conjoin( grid2.elements );        
        grid1.faces.conjoin( grid2.faces );
        grid1.edges.conjoin( grid2.edges );
        grid1.element2faces.conjoin( grid2.element2faces );
        grid1.face2edges.conjoin( grid2.face2edges );        
        grid1.level.conjoin( grid2.level );
        grid1.rightChild.conjoin( grid2.rightChild );
        grid1.levelFaces.conjoin( grid2.levelFaces );
        
        grid1.coordinates.conjoin(grid2.coordinates);

        decomposition1.node2ranks.conjoin( decomposition2.node2ranks );
        decomposition1.edge2ranks.conjoin( decomposition2.edge2ranks );
        decomposition1.face2ranks.conjoin( decomposition2.face2ranks );
        
        globalID1.globalNodeIndex.conjoin( globalID2.globalNodeIndex );
        
        manager.conjoin();
    }
    
    template < class Grid, class Decomposition, class GlobalId, class Manager >
    void removeMultipleIndices( Grid &grid, Decomposition &decomposition, GlobalId & globalID, Manager & manager )
    {
        ManagedArray<int> uniqueNodes, nodes2newNodes, uniqueFaces, uniqueEdges, edges2newEdges, faces2newFaces;
        
        globalID.globalNodeIndex = unique_rows( globalID.globalNodeIndex, uniqueNodes, nodes2newNodes );

        index2newIndex( grid.elements, nodes2newNodes );
        index2newIndex( grid.faces, nodes2newNodes );
        index2newIndex( grid.edges, nodes2newNodes );

        sortRows( grid.faces );
        unique_rows_void( grid.faces, uniqueFaces, faces2newFaces );
        index2newIndex( grid.element2faces, faces2newFaces );
        grid.face2edges = rearange( grid.face2edges, uniqueFaces );
        getFacesForElements( grid.elements, grid.faces, grid.element2faces, grid.level ); 
            
        sortRows( grid.edges );
        unique_rows_void( grid.edges, uniqueEdges, edges2newEdges );
        index2newIndex( grid.face2edges, edges2newEdges );
        getFacesForElements( grid.faces, grid.edges, grid.face2edges );
              
        grid.coordinates = rearange( grid.coordinates, uniqueNodes ); 

        decomposition.node2ranks = rearange( decomposition.node2ranks, uniqueNodes );
        decomposition.edge2ranks = rearange( decomposition.edge2ranks, uniqueEdges );
        decomposition.face2ranks = rearange( decomposition.face2ranks, uniqueFaces );
        
        grid.levelFaces = rearange( grid.levelFaces, uniqueFaces );
      
        manager.newIndices( uniqueFaces, uniqueEdges, uniqueNodes );
    }
    
    
#if 0
    template < class Grid, class Decomposition, class GlobalId, class Manager >
    void mergeTwoGrids( Grid &grid1, Decomposition &decomposition1, GlobalId & globalID1,
                        Grid &grid2, Decomposition &decomposition2, GlobalId & globalID2, Manager & manager )
    {
        ManagedArray<int> nodes2newNodes1(grid1.nNodes(),-1);
        ManagedArray<int> nodes2newNodes2(grid2.nNodes(),-1);
          
        ManagedArray<int,2> globalIndices;
        
        globalID1.globalNodeIndex.conjoin(globalID2.globalNodeIndex);
        
        ManagedArray<int> uniqueNodes, nodes2newNodes;
          
        globalID1.globalNodeIndex = unique_rows_nosort( globalID1.globalNodeIndex, uniqueNodes, nodes2newNodes );
        
        assert( globalID1.globalNodeIndex.size() == uniqueNodes.size() );
    
        for( size_t k = 0; k < nodes2newNodes1.size(); k++)
            nodes2newNodes1[k] = nodes2newNodes[k];
        for( size_t k = 0; k < nodes2newNodes2.size(); k++)
            nodes2newNodes2[k] = nodes2newNodes[k + grid1.nNodes()];
        
        index2newIndex( grid1.elements, nodes2newNodes1 );
        index2newIndex( grid1.faces, nodes2newNodes1 );
        index2newIndex( grid1.edges, nodes2newNodes1 );
        
        index2newIndex( grid2.elements, nodes2newNodes2 );
        index2newIndex( grid2.faces, nodes2newNodes2 );
        index2newIndex( grid2.edges, nodes2newNodes2 );  
        
        grid1.elements.conjoin( grid2.elements );
        size_t nFaces1 = grid1.nFaces();
        size_t nEdges1 = grid1.nEdges();
        
        grid1.faces.conjoin( grid2.faces );
        grid1.edges.conjoin( grid2.edges );
        
        grid2.element2faces += nFaces1;
        grid2.face2edges += nEdges1;
        
        grid1.element2faces.conjoin( grid2.element2faces );
        grid1.face2edges.conjoin( grid2.face2edges );        
        
        ManagedArray<int> uniqueFaces, uniqueEdges, edges2newEdges, faces2newFaces;
        
        sortRows(grid1.faces);
        unique_rows_void( grid1.faces, uniqueFaces, faces2newFaces );
        sortRows(grid1.edges);
        unique_rows_void( grid1.edges, uniqueEdges, edges2newEdges );
 
        index2newIndex( grid1.element2faces, faces2newFaces );
        index2newIndex( grid1.face2edges, edges2newEdges );
        
        grid1.face2edges = rearange( grid1.face2edges, uniqueFaces );
        
        grid1.level.conjoin( grid2.level );
        grid1.rightChild.conjoin( grid2.rightChild );
        
        grid1.coordinates.conjoin(grid2.coordinates);
        grid1.coordinates = rearange( grid1.coordinates, uniqueNodes );

        decomposition1.node2ranks.conjoin( decomposition2.node2ranks );
        decomposition1.node2ranks = rearange( decomposition1.node2ranks, uniqueNodes );
        decomposition1.edge2ranks.conjoin( decomposition2.edge2ranks );
        decomposition1.edge2ranks = rearange( decomposition1.edge2ranks, uniqueEdges );
        decomposition1.face2ranks.conjoin( decomposition2.face2ranks );
        decomposition1.face2ranks = rearange( decomposition1.face2ranks, uniqueFaces );
        
        grid1.levelFaces.conjoin( grid2.levelFaces );
        grid1.levelFaces = rearange( grid1.levelFaces, uniqueFaces );
  
        getFacesForElements( grid1.elements, grid1.faces, grid1.element2faces, grid1.level );
        getFacesForElements( grid1.faces, grid1.edges, grid1.face2edges );
        
        assert( grid1.face2edges.max() + 1 == grid1.edges.size() );
        assert( uniqueEdges.size() == grid1.edges.size() ); 
        
        manager.merge( faces2newFaces, edges2newEdges, nodes2newNodes, uniqueFaces, uniqueEdges, uniqueNodes );
    }
#endif
#if 0 
    template < class Grid, class Decomposition, class GlobalId, class ... Manager >
    void mergeSubgrids( int myRank, std::vector<int> & subgridRanks, std::vector<Grid> & subgrids, std::vector<Decomposition> & decompositions, 
                        std::vector<GlobalId> & globalIDs, std::vector<Manager> & ... manager )
    {
        std::vector< ManagedArray<int> > nodes2newNodes( globalIDs.size() );
        size_t counter = 0;
        for( size_t k = 0; k < globalIDs.size(); k++){
            counter += globalIDs[k].globalNodeIndex.size();
            nodes2newNodes[k].resize( globalIDs[k].globalNodeIndex.size() );
        }
        ManagedArray<int,2> combinedGlobalIds(counter);
        for( size_t k = 0, ctr = 0; k < globalIDs.size(); k++)
            for(size_t j = 0; j < globalIDs[k].globalNodeIndex.size(); j++){
                combinedGlobalIds(ctr,0) = globalIDs[k].globalNodeIndex(j,0);
                combinedGlobalIds(ctr,1) = globalIDs[k].globalNodeIndex(j,1);
            }
        ManagedArray<int> index(counter);    
        ManagedArray<int> combinedNodes2newNodes(counter);
        unique_rows_void( combinedGlobalIds, index, combinedNodes2newNodes);
        
        for( size_t k = 0, ctr = 0; k < globalIDs.size(); k++)
            for(size_t j = 0; j < globalIDs[k].globalNodeIndex.size(); j++)
                nodes2newNodes[k][j] = combinedNodes2newNodes[ctr];
    }
#endif

#if 0  
    template < class ElementVector, class Mesh, class ParallelData, class GlobalIndex >
    void cutSubMeshes( int p, int numElements_p, ElementVector &partition, Mesh &mesh, Mesh &mesh_p,
                       ParallelData &pm, ParallelData &pm_p, GlobalIndex &globalIndices, GlobalIndex &globalIndices_p )
    {
        const int NV = Mesh::verticesPerElement;
        const int dimworld = Mesh::dimension;
            
        ManagedArray<int> nodes2newNodes( mesh.c4n.size(), -1 );
        ManagedArray<int> edges2newEdges( mesh.edgeData.nEdges, -1 );
        ManagedArray<int> faces2newFaces( mesh.neighbors.nFaces, -1 );

        for( int el = 0; el < (int)mesh.n4e.size(); el++){
            if( partition[el] == p ){
                for ( int k = 0; k < NV; k++)
                    nodes2newNodes[ mesh.n4e(el,k) ] = 1;
                for ( int k = 0; k < mesh.edgeData.element2edges.dimension(); k++)
                    edges2newEdges[ mesh.edgeData.element2edges(el,k) ] = 1;
                for ( int k = 0; k < NV; k++)
                    faces2newFaces[ mesh.neighbors.faceIndex(el,k) ] = 1;
            }
        }
        
        int nNodes = 0;
        for( int k = 0; k < (int)nodes2newNodes.size(); k++)
            if( nodes2newNodes[k] > -1 )
                nodes2newNodes[k] = nNodes++;
        
        int nEdges = 0;
        for( int k = 0; k < (int)edges2newEdges.size(); k++)
            if( edges2newEdges[k] > -1 )
                edges2newEdges[k] = nEdges++;
            
        mesh_p.resizeMesh( numElements_p );
        mesh_p.edgeData.element2edges.resize( numElements_p );
        mesh_p.edgeData.nEdges = nEdges;
        mesh_p.resizeNodes( nNodes );
        
        int nElements = 0;
        
        for( int el = 0; el < (int)mesh.n4e.size(); el++)
            if( partition[el] == p ){
                for ( int k = 0; k < NV; k++)
                    mesh_p.n4e(nElements,k) = nodes2newNodes[ mesh.n4e(el,k)];
                mesh_p.level[nElements] = mesh.level[el];
                for( int k = 0; k < (int)mesh.edgeData.element2edges.dimension(); ++k )
                    mesh_p.edgeData.element2edges(nElements,k) = edges2newEdges[mesh.edgeData.element2edges(el,k)];
                nElements++;
            }
        
        for( int k = 0; k < (int)nodes2newNodes.size(); k++)
            if( nodes2newNodes[k] > -1 ) 
                for( int d = 0; d < dimworld; d++ )
                    mesh_p.c4n( nodes2newNodes[k], d ) = mesh.c4n( k, d);
                
        pm_p.vertex2ranks.resize( nNodes );
        for( int k = 0; k < (int)nodes2newNodes.size(); ++k)
            if( nodes2newNodes[k] > -1 )
                pm_p.vertex2ranks[ nodes2newNodes[k] ] = pm.vertex2ranks[ k ];
        
        pm_p.edge2ranks.resize( nEdges );
        for( int k = 0; k < (int)edges2newEdges.size(); ++k)
            if( edges2newEdges[k] > -1 )
                pm_p.edge2ranks[ edges2newEdges[k] ] = pm.edge2ranks[ k ];
                
        mesh.neighbors.getNeighborForOneSubmesh( mesh_p, p, numElements_p, partition, faces2newFaces );
        
        globalIndices_p.globalNodeIndex.resize(nNodes);
        for( int k = 0; k < (int)nodes2newNodes.size(); k++)
            if( nodes2newNodes[k] > -1 ){
                globalIndices_p.globalNodeIndex(nodes2newNodes[k],0) = globalIndices.globalNodeIndex(k,0);
                globalIndices_p.globalNodeIndex(nodes2newNodes[k],1) = globalIndices.globalNodeIndex(k,1);
            }   
    }
#endif  
    
#if 0
    template < class Mesh, class Communicator, class ParallelData >
    class TemporalGlobalIndices
    {
    public:
        const int NV = Mesh::verticesPerElement;
        
        ManagedArray<int,2> globalNodeIndex;
        ManagedArray<int,2> globalEdgeIndex;
        ManagedArray<int,2> globalFaceIndex;
        
        ManagedArray<int,2> globalNodeIndex_temp;
        ManagedArray<int,2> globalEdgeIndex_temp;
        ManagedArray<int,2> globalFaceIndex_temp;
        
        TemporalGlobalIndices( Mesh &mesh, Communicator &comm, ParallelData &pm)
        : mesh_(mesh), comm_(comm), pm_(pm)
        {
            globalNodeIndex.resize( mesh.c4n.size() );
            globalEdgeIndex.resize( mesh_.edgeData.nEdges );
            globalFaceIndex.resize( mesh_.neighbors.nFaces );
        }
        
        void communicate()
        {   
            //nodes
            ManagedArray<int> nodes( mesh_.c4n.size(), -1 );
            globalNodeIndex.resize( mesh_.c4n.size() );
            
            for( int k = 0; k < (int)nodes.size(); ++k)
                globalNodeIndex(k,0) = (int)comm_.rank();
                        
            for( int p = 0; p < (int)comm_.size(); ++p )
                for( int k = 0; k < (int)pm_.rank2bdyNodes[p].size(); ++k ){
                    if( pm_.vertex2ranks[pm_.rank2bdyNodes[p][k]][0] > (int)comm_.rank() ){
                        nodes[ pm_.rank2bdyNodes[p][k] ] = pm_.rank2bdyNodes[p][k];
                    } else {
                        nodes[ pm_.rank2bdyNodes[p][k] ] = 0;
                        globalNodeIndex( pm_.rank2bdyNodes[p][k], 0 ) = pm_.vertex2ranks[pm_.rank2bdyNodes[p][k]][0];
                    }
                }
            for( int k = 0; k < (int)nodes.size(); k++)
                if( nodes[k] < 0 )
                    nodes[k] = k;
            
            pm_.sumVertexVector( nodes );
                
            for( int k = 0; k < (int)nodes.size(); k++)
                globalNodeIndex(k,1) = nodes[k];
            
            //edges 
            ManagedArray<int> edges( mesh_.edgeData.nEdges, -1 );
            globalEdgeIndex.resize( edges.size() );
            
            for( int k = 0; k < (int)edges.size(); ++k)
                globalEdgeIndex(k,0) = (int)comm_.rank();
                        
            for( int p = 0; p < (int)comm_.size(); ++p )
                for( int k = 0; k < (int)pm_.rank2bdyEdges[p].size(); ++k ){
                    if( pm_.edge2ranks[pm_.rank2bdyEdges[p][k]][0] > (int)comm_.rank() ){
                        edges[ pm_.rank2bdyEdges[p][k] ] = pm_.rank2bdyEdges[p][k];
                    } else {
                        edges[ pm_.rank2bdyEdges[p][k] ] = 0;
                        globalEdgeIndex( pm_.rank2bdyEdges[p][k], 0 ) = pm_.edge2ranks[pm_.rank2bdyEdges[p][k]][0];
                    }
                }
                
            for( int k = 0; k < (int)edges.size(); k++)
                if( edges[k] < 0 )
                    edges[k] = k;
            
            pm_.sumEdgeVector( edges );
                
            for( int k = 0; k < (int)edges.size(); k++)
                globalEdgeIndex(k,1) = edges[k];
            
            //faces
            ManagedArray<int> faces( mesh_.neighbors.nFaces );
            globalFaceIndex.resize( mesh_.neighbors.nFaces );
            for( int k = 0; k < (int)globalFaceIndex.size(); ++k ){
                globalFaceIndex(k,0) = (int)comm_.rank();
                faces[k] = k;
            }
            
            for( int el = 0; el < (int)mesh_.n4e.size(); ++el )
                for( int k = 0; k < NV; ++k )
                    if( mesh_.neighbors.isParallel(el,k) && mesh_.neighbors.rank(el,k) < (int)comm_.rank()) {
                        globalFaceIndex( mesh_.neighbors.faceIndex(el,k), 0 ) = mesh_.neighbors.rank(el,k);
                        faces[ mesh_.neighbors.faceIndex(el,k) ] = 0;
                    }
                    
            pm_.sumFaceVector( faces );

            for( int k = 0; k < (int)globalFaceIndex.size(); ++k )
                globalFaceIndex(k,1) = faces[k];

        } 
        
        template<class NodeVector, class EdgeVector, class FaceVector >
        void storeTemporaryActiveEntities( NodeVector &nodes2newNodes, EdgeVector &edges2newEdges, FaceVector &faces2newFaces )
        {
            assert( globalNodeIndex.size() == nodes2newNodes.size() );
            assert( globalEdgeIndex.size() == edges2newEdges.size() );
            assert( globalFaceIndex.size() == faces2newFaces.size() );
            
            globalNodeIndex_temp.resize( globalNodeIndex.size() );
            globalEdgeIndex_temp.resize( globalEdgeIndex.size() );
            globalFaceIndex_temp.resize( globalFaceIndex.size() );

            int count_nodes = 0;
            for( int k = 0; k < (int)nodes2newNodes.size(); ++k )
                if( nodes2newNodes[k] > -1 ){
                    globalNodeIndex_temp( nodes2newNodes[k], 0 )  = globalNodeIndex(k,0);
                    globalNodeIndex_temp( nodes2newNodes[k], 1 )  = globalNodeIndex(k,1);
                    count_nodes++;
                }
            globalNodeIndex_temp.resize( count_nodes );
            
            int count_edges = 0;
            for( int k = 0; k < (int)edges2newEdges.size(); ++k )
                if( edges2newEdges[k] > -1 ){
                    globalEdgeIndex_temp( edges2newEdges[k], 0 )  = globalEdgeIndex(k,0);
                    globalEdgeIndex_temp( edges2newEdges[k], 1 )  = globalEdgeIndex(k,1);
                    count_edges++;
                }
            globalEdgeIndex_temp.resize( count_edges );
            
            int count_faces = 0;
            for( int k = 0; k < (int)faces2newFaces.size(); ++k )
                if( faces2newFaces[k] > -1 ){
                    globalFaceIndex_temp( faces2newFaces[k], 0 )  = globalFaceIndex(k,0);
                    globalFaceIndex_temp( faces2newFaces[k], 1 )  = globalFaceIndex(k,1);
                    count_faces++;
                }
            globalFaceIndex_temp.resize( count_faces );
            
        }
                
        void pack( void *buffer, std::size_t bufferSize, int &position ) 
        {            
            int nNodes = (int)globalNodeIndex_temp.size();
            int nEdges = (int)globalEdgeIndex_temp.size();
            int nFaces = (int)globalFaceIndex_temp.size();
            
            comm_.pack( &nNodes, 1, buffer, bufferSize, position);
            comm_.pack( &nEdges, 1, buffer, bufferSize, position);
            comm_.pack( &nFaces, 1, buffer, bufferSize, position);
            comm_.pack( globalNodeIndex_temp.data()(), 2*nNodes, buffer, bufferSize, position);
            comm_.pack( globalEdgeIndex_temp.data()(), 2*nEdges, buffer, bufferSize, position);
            comm_.pack( globalFaceIndex_temp.data()(), 2*nFaces, buffer, bufferSize, position);
        }
        void unpack( void *buffer, std::size_t bufferSize, int &position ) 
        {
            int nNodes;
            int nEdges;
            int nFaces;
            
            comm_.unpack( &nNodes, 1, buffer, bufferSize, position);
            comm_.unpack( &nEdges, 1, buffer, bufferSize, position);
            comm_.unpack( &nFaces, 1, buffer, bufferSize, position);
            
            globalNodeIndex.resize(nNodes);
            globalEdgeIndex.resize(nEdges);
            globalFaceIndex.resize(nFaces);
            
            comm_.unpack( globalNodeIndex.data()(), 2*nNodes, buffer, bufferSize, position);
            comm_.unpack( globalEdgeIndex.data()(), 2*nEdges, buffer, bufferSize, position);
            comm_.unpack( globalFaceIndex.data()(), 2*nFaces, buffer, bufferSize, position);
        }

        std::size_t bufferSize()
        {
            return comm_.template bufferSize<int>( 3 + 2*globalNodeIndex_temp.size() + 2*globalEdgeIndex_temp.size() + 2*globalFaceIndex_temp.size()); 
        }
        
    private:
        Mesh &mesh_;
        Communicator &comm_;
        ParallelData &pm_;
        
    };
#endif


#if 0  
    template < class Mesh, class Communicator, class ParallelData, class ElementVector, class ... Manager >
    void redistributeParallelMesh ( Mesh &mesh, Communicator &comm, ParallelData &pm, ElementVector &partition, Manager & ... manager )
    {
        Global_Id<Mesh, Communicator, ParallelData> globalIndices(mesh, comm, pm);
        globalIndices.communicate();
        
        updateEntity2Ranks( mesh, comm, pm, partition );
     
        sendEntity2Ranks3( mesh, comm, pm, partition );

        //ManagedArray<int> nodes2newNodes, edges2newEdges, faces2newFaces;
        //sortByGlobalIndices( mesh, pm, globalIndices, nodes2newNodes, edges2newEdges, faces2newFaces );              
        
        ManagedArray< int > numElements;
        ManagedArray< int > numElementsOtherRank;
        
        communicateNumbersOfElementsToSend( mesh, comm, partition, numElements, numElementsOtherRank );
                
        sendSubMeshes( mesh, comm, pm, globalIndices, partition, numElements, manager ... );

        deleteInactiveElements( mesh, comm, pm, partition, globalIndices, manager ... );

        mergeWithRecvMeshes( mesh, comm, pm, globalIndices, partition, numElementsOtherRank, manager ... );
        
        pm.init_bdyData();
        
        printf(" merging on %d finished ... \n", (int)comm.rank() );
            
        assert( mesh.edgeData.check( mesh ) );
        assert( mesh.neighbors.check( ) );
        assert( pm.vertex2ranks.size() == mesh.c4n.size() );
        assert( (int)pm.edge2ranks.size() == mesh.edgeData.nEdges );
        assert( pm.testVertexComm() );
    }
#endif 

#if 0
    template< class Mesh, class Communicator, class ElementVector, class RankVector >
    void communicateNumbersOfElementsToSend( Mesh & mesh, Communicator & comm, ElementVector & partition, RankVector & numElements, 
                                             RankVector & numElementsOtherRank )
    {
        numElements.resize( comm.size());
        numElements.fill(0);
        numElementsOtherRank.resize( comm.size() );
        numElementsOtherRank.fill(0);
       
        for( int el = 0; el < (int)mesh.n4e.size(); el++){
            const int p = partition[el];
            assert( p > -1 );
            assert( p < comm.size() );
            numElements[p]++;
        }
        
        std::vector< typename Communicator::Request > requests( comm.size() );
        std::vector< int > flag( comm.size(), 0 );
        
        for( int p = 0; p < comm.size(); p++ ){
            if( p != comm.rank() )
                comm.Isend( &numElements[p], 1, p, 1, &requests[p] );
        }
        
        for( int p = 0; p < comm.size(); p++ )
            if( p != comm.rank() )
                comm.Irecv( &numElementsOtherRank[p], 1, p, 1, &requests[p] );
        
        int counter = 1;
        while( counter > 0){
            counter = 0;
            for( int p = 0; p < comm.size(); p++ )
                if( p != comm.rank() && flag[p] == 0 ){
                    comm.test( &requests[p], &flag[p] );
                    if( flag[p] == 0 )
                        counter++;
                }
        }
    }
#endif

#if 0
    template < class Mesh, class Communicator, class ParallelData, class GlobalIndex, class ElementVector, class RankVector, class ... Manager >
    void sendSubMeshes ( Mesh &mesh, Communicator &comm, ParallelData &pm, GlobalIndex &globalIndices, ElementVector &partition, 
                         RankVector & numElements, Manager & ... manager )
    {
        for( int p = 0; p < comm.size(); p++ ) {
            if( p == comm.rank() )
                continue;
            if( numElements[p] == 0 )
                continue;
            Mesh mesh_p;
            ParallelMeshData< Mesh, Communicator > pm_p( mesh_p, comm);
            Global_Id< Mesh, Communicator, ParallelData > globalIndices_p( mesh_p, comm, pm_p );
            cutSubMeshes(p, numElements[p], partition, mesh, mesh_p, pm, pm_p, globalIndices, globalIndices_p );
            auto sentMesh = sendMesh( mesh_p, p, 1, comm, pm_p, globalIndices_p, manager ... );
        }
    }
    


    template < class Mesh, class Communicator, class ParallelData, class GlobalIndices, class ElementVector, class RankVector, class ... Manager >
    void mergeWithRecvMeshes( Mesh &mesh, Communicator &comm, ParallelData &pm, GlobalIndices &globalIndices, 
                              ElementVector &partition, RankVector &numElementsOtherRank, Manager & ... manager )
    {
        const int NV = Mesh::verticesPerElement;
        const int dimworld = Mesh::dimension;
                
        for( int p = 0; p < comm.size(); p++ ) {
            if( p == comm.rank() )
                continue;
            if( numElementsOtherRank[p] == 0 )
                continue;

            Mesh mesh_p;
            ParallelData pm_p( mesh_p, comm );
            Global_Id<Mesh, Communicator, ParallelData> globalIndices_p(mesh_p, comm, pm_p);
            
            mesh_p = recvMesh<Mesh>( p, 1, comm, pm_p, globalIndices_p, manager ... ); 
            
            ManagedArray<int> nodes2newNodes_1, nodes2newNodes_2;
            ManagedArray<int> edges2newEdges_1, edges2newEdges_2;
            ManagedArray<int> faces2newFaces_1, faces2newFaces_2;
            
            int nNodes, nEdges, nFaces;
            
            globalIndices.getUniqueEntitieIndices( mesh, mesh_p, globalIndices_p, nodes2newNodes_1, nodes2newNodes_2,
                                                   edges2newEdges_1, edges2newEdges_2, faces2newFaces_1, faces2newFaces_2,
                                                   nNodes, nEdges, nFaces );

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            //conjoin coordinate data 
            
            ManagedArray<double, dimworld> c4n( nNodes );
            for( int k = 0; k < (int)nodes2newNodes_1.size(); ++k )
                for( int d = 0; d < dimworld; ++d )
                    c4n( nodes2newNodes_1[k], d ) = mesh.c4n(k,d);
            
            for( int k = 0; k < (int)nodes2newNodes_2.size(); ++k )
                for( int d = 0; d < dimworld; ++d )
                    c4n( nodes2newNodes_2[k], d ) = mesh_p.c4n(k,d);
                
            mesh.c4n = c4n;  
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            //conjoin entity2ranks parallel data
            ManagedArray< ManagedArray<int> > vertex2ranks( nNodes );
            for( int k = 0; k < (int)nodes2newNodes_1.size(); ++k)
                vertex2ranks[ nodes2newNodes_1[k] ] = pm.vertex2ranks[k];
            for( int k = 0; k < (int)nodes2newNodes_2.size(); ++k)
                vertex2ranks[ nodes2newNodes_2[k] ].conjoin( pm_p.vertex2ranks[k] );

            for(int k = 0; k < (int)vertex2ranks.size(); ++k)
                if( vertex2ranks[k].size() )
                    vertex2ranks[k] = vertex2ranks[k].unique();
                
            pm.vertex2ranks = vertex2ranks;    
                
            ManagedArray< ManagedArray<int> > edge2ranks( nEdges );
            for( int k = 0; k < (int)edges2newEdges_1.size(); ++k)
                edge2ranks[ edges2newEdges_1[k] ] = pm.edge2ranks[k];
            for( int k = 0; k < (int)edges2newEdges_2.size(); ++k)
                edge2ranks[ edges2newEdges_2[k] ].conjoin( pm_p.edge2ranks[k] );
                
            for(int k = 0; k < (int)edge2ranks.size(); ++k)
                if( edge2ranks[k].size() )
                    edge2ranks[k] = edge2ranks[k].unique();
                               
            pm.edge2ranks = edge2ranks;           
                
            for( int el = 0; el < (int)mesh_p.n4e.size(); ++el )
                for( int k = 0; k < NV; ++k){
                    if( mesh_p.neighbors.isBoundary(el,k) )
                        continue;
                    mesh_p.neighbors.traits_(el,k).setIndex( mesh_p.neighbors.index(el,k) + (int)mesh.n4e.size() );
                }
            
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // conjoin mesh data
            mesh.n4e.conjoin( mesh_p.n4e );
            mesh.level.conjoin( mesh_p.level );
            mesh.neighbors.traits_.conjoin( mesh_p.neighbors.traits_ );
            mesh.edgeData.element2edges.conjoin( mesh_p.edgeData.element2edges );

            mesh.marked.resize( mesh_p.n4e.size() );
            
            mesh.neighbors.nFaces = nFaces;
            mesh.edgeData.nEdges = nEdges;
            
            //glue new neigbors together
            ManagedArray<int,2> faces2Elements( nFaces, -1 );
            for( int el = 0; el < (int)mesh.n4e.size(); ++el )
                for( int k = 0; k < NV; ++k)
                    if( mesh.neighbors.isBoundary(el,k) && mesh.neighbors.isParallel(el,k) ){
                        const int f = mesh.neighbors.faceIndex(el,k);
                        if( faces2Elements( f, 0 ) == -1 ) {
                            faces2Elements( f, 0 ) = el;
                            faces2Elements( f, 1 ) = k;
                        } else {
                            const int N = faces2Elements(f,0);
                            const int N_k = faces2Elements(f,1);
                            mesh.neighbors.traits_(N,N_k).setIndex( el );
                            mesh.neighbors.traits_(N,N_k).setBoundary(false);
                            mesh.neighbors.traits_(N,N_k).setParallel(false);
                            mesh.neighbors.traits_(el,k).setIndex( N );
                            mesh.neighbors.traits_(el,k).setBoundary(false);
                            mesh.neighbors.traits_(el,k).setParallel(false); 
                        }
                    }

           //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           
        } // for( int p = 0; p < comm.size(); p++ )
    }

    template < class Mesh, class Communicator, class ParallelData, class ElementVector, class Global_Id, class ... Manager >
    void deleteInactiveElements( Mesh &mesh, Communicator &comm, ParallelData &pm, ElementVector &partition, Global_Id &globalIndices, Manager & ... manager )
    {
        const int NV = Mesh::verticesPerElement;
        const int dimworld = Mesh::dimension;
        
        ManagedArray<int> elements2newElements( mesh.n4e.size(), -1 );
        ManagedArray<int> faces2newFaces( mesh.neighbors.nFaces, -1 );
        ManagedArray<int> edges2newEdges( mesh.edgeData.nEdges, -1 );
        ManagedArray<int> nodes2newNodes( mesh.c4n.size(), -1 );
        
        //update rank information
        for( int el = 0; el < (int)partition.size(); ++el ){
            if( partition[el] != comm.rank() )
                continue;
            for( int k = 0; k < NV; ++k )
                if( !mesh.neighbors.isBoundary(el,k) ){
                    const int N = mesh.neigh(el,k);
                    if( !(N < (int)partition.size() ) )
                        continue;
                    if(  partition[N] == comm.rank() )
                        continue;
                    mesh.neighbors.traits_(el,k).setIndex( partition[N] );
                    mesh.neighbors.traits_(el,k).setBoundary(true);
                    mesh.neighbors.traits_(el,k).setParallel(true);
                }
        }
        
        int nElements = 0;
        for( int el = 0; el < (int)mesh.n4e.size(); ++el ){
            if( el < (int)partition.size() && partition[el] != (int)comm.rank() )
                continue;
            
            for( int k = 0; k < NV; ++k )
                mesh.n4e( nElements, k ) = mesh.n4e( el, k );
            mesh.level[ nElements ] = mesh.level[el];
            for( int k = 0; k < NV; ++k )
                mesh.neighbors.traits_(nElements,k) = mesh.neighbors.traits_(el,k);
            for( int k = 0; k < mesh.edgeData.element2edges.dimension(); ++k )
                mesh.edgeData.element2edges(nElements,k) = mesh.edgeData.element2edges(el,k);
            
            elements2newElements[el] = nElements++;
            for( int k = 0; k < NV; ++k )
                faces2newFaces[ mesh.neighbors.faceIndex(el,k) ] = 1;
            for( int k = 0; k < mesh.edgeData.element2edges.dimension(); ++k )
                edges2newEdges[ mesh.edgeData.element2edges( el, k ) ] = 1;
            for( int k = 0; k < NV; ++k )
                nodes2newNodes[ mesh.n4e(el,k) ] = 1;
        }
        mesh.n4e.resize( nElements );
        mesh.level.resize( nElements );
        mesh.neighbors.resize( nElements);
        mesh.edgeData.element2edges.resize( nElements );
                
        int nFaces = 0; int nEdges = 0; int nNodes = 0;
        for( int k = 0; k < (int)faces2newFaces.size(); ++k )
            if( faces2newFaces[k] > -1 )
                faces2newFaces[k] = nFaces++;
        for( int k = 0; k < (int)edges2newEdges.size(); ++k )
            if( edges2newEdges[k] > -1 )
                edges2newEdges[k] = nEdges++;
        for( int k = 0; k < (int)nodes2newNodes.size(); ++k )
            if( nodes2newNodes[k] > -1 )
                nodes2newNodes[k] = nNodes++;

                        
        for( int j = 0; j < (int)mesh.c4n.size(); ++j )
            if( nodes2newNodes[j] > -1 ){
                pm.vertex2ranks[ nodes2newNodes[j] ] = pm.vertex2ranks[ j ];
                for( int d = 0; d < dimworld; ++d )
                    mesh.c4n( nodes2newNodes[j], d ) = mesh.c4n( j, d );
            }
        mesh.c4n.resize( nNodes );
        pm.vertex2ranks.resize( nNodes );
        
        for( int j = 0; j < (int)pm.edge2ranks.size(); ++j )
            if( edges2newEdges[j] > -1 )
                pm.edge2ranks[ edges2newEdges[j] ] = pm.edge2ranks[ j ];
        pm.edge2ranks.resize( nEdges );
        
        for( int el = 0; el < (int)mesh.n4e.size(); ++el ){
            for( int k = 0; k < NV; ++k )
                mesh.n4e( el, k ) = nodes2newNodes[ mesh.n4e( el, k ) ];
            for( int k = 0; k < NV; ++k )
                if( !mesh.neighbors.isBoundary(el,k) )
                    mesh.neighbors.traits_(el,k).setIndex( elements2newElements[ mesh.neighbors.traits_(el,k).index() ] );
            for( int k = 0; k < NV; ++k )
                mesh.neighbors.traits_(el,k).setFaceNumber( faces2newFaces[ mesh.neighbors.traits_(el,k).faceNumber() ] );
            for( int k = 0; k < mesh.edgeData.element2edges.dimension(); ++k )
                mesh.edgeData.element2edges(el,k) = edges2newEdges[ mesh.edgeData.element2edges(el,k) ];
        }
        
        for( int k = 0; k < (int)nodes2newNodes.size(); k++)
            if( nodes2newNodes[k] > -1 ){
                globalIndices.globalNodeIndex(nodes2newNodes[k],0) = globalIndices.globalNodeIndex(k,0);
                globalIndices.globalNodeIndex(nodes2newNodes[k],1) = globalIndices.globalNodeIndex(k,1);
            } 
        globalIndices.globalNodeIndex.resize(nNodes);
        
        mesh.neighbors.nFaces = nFaces;
        mesh.edgeData.nEdges = nEdges;
    }
  #endif   
#if 0   
    //TODO this methode might vanish --> sendNewEntity2RankInformation
    template < class Mesh, class ElementVector, class ParallelData, class Communicator >
    void sendEntity2Ranks( Mesh &mesh, Communicator &comm, ParallelData &pm, ElementVector &partition )
    {
        const int NV = Mesh::verticesPerElement;
        
        ManagedArray< ManagedArray<int> > buffer_vertices( comm.size() );
        ManagedArray< ManagedArray<int> > buffer_edges( comm.size() );

        ManagedArray< ManagedArray<int> > buffer_vertex_size( comm.size() );
        ManagedArray< ManagedArray<int> > buffer_edge_size( comm.size() );
        
        ManagedArray<int> buffer_size( comm.size(), 0 );
        
        for( int p = 0; p < (int)comm.size(); p++) {
            buffer_vertex_size[p].resize( pm.rank2bdyNodes[p].size() );
            for( int k = 0; k < (int)pm.rank2bdyNodes[p].size(); k++){
               buffer_vertex_size[p][k] = (int)pm.vertex2ranks[ pm.rank2bdyNodes[p][k] ].size();
               buffer_size[p] += (int)pm.vertex2ranks[ pm.rank2bdyNodes[p][k] ].size();
            }
        }
        
        for( int p = 0; p < (int)comm.size(); p++)
            buffer_vertices[p].resize( buffer_size[p] );
        
        buffer_size.fill( 0 );
        for( int p = 0; p < (int)comm.size(); p++) {
            buffer_edge_size[p].resize( pm.rank2bdyEdges[p].size() );
            for( int k = 0; k < (int)pm.rank2bdyEdges[p].size(); k++){
               buffer_edge_size[p][k] = (int)pm.edge2ranks[ pm.rank2bdyEdges[p][k] ].size();
               buffer_size[p] += (int)pm.edge2ranks[ pm.rank2bdyEdges[p][k] ].size();
            }
        }
        
        for( int p = 0; p < (int)comm.size(); p++)
            buffer_edges[p].resize( buffer_size[p] );
        
        for( int p = 0; p < (int)comm.size(); p++){
            int j = 0;
            for( int k = 0; k < (int)pm.rank2bdyNodes[p].size(); k++ )
                for( int r = 0; r < buffer_vertex_size[p][k]; r++ ) 
                    buffer_vertices[p][j++] = pm.vertex2ranks[ pm.rank2bdyNodes[p][k] ][r];
            j = 0;
            for( int k = 0; k < (int)pm.rank2bdyEdges[p].size(); k++ )
                for( int r = 0; r < buffer_edge_size[p][k]; r++ ) 
                    buffer_edges[p][j++] = pm.edge2ranks[ pm.rank2bdyEdges[p][k] ][r];
        }
        
        //faces //TODO
        ManagedArray<int> Face2InnerRank( mesh.neighbors.nFaces );
        for( int el = 0; el < (int)mesh.n4e.size(); el++ )
            for( int k = 0; k < NV; k++ )
                Face2InnerRank[ mesh.neighbors.faceIndex(el,k) ] = partition[el];
            
        ManagedArray< ManagedArray<int> > buffer_faces( comm.size() );
        for( int p = 0; p < (int)comm.size(); p++ ) {
            if( pm.rank2bdyFaces[p].size() > 0 ){
                buffer_faces[p].resize( pm.rank2bdyFaces[p].size() );
                for( int k = 0; k < (int)pm.rank2bdyFaces[p].size(); k++ )
                    buffer_faces[p][k] = Face2InnerRank[ pm.rank2bdyFaces[p][k] ];
            }
        }

        //communication
        std::vector< typename Communicator::Request > requests( comm.size() );
        std::vector< typename Communicator::Status > status( comm.size() );
        
        std::vector<bool> message( comm.size() );
        for(int p = 0; p < comm.size(); p++ ) message[p] = true;
        
        for(int p = 0; p < (int)comm.size(); p++ ){
            if( buffer_vertices[p].size() < 1 ){
            template<class Mesh, class ParallelData, class GlobalIndices, class NodeVector, class EdgeVector, class FaceVector >
    void sortByGlobalIndices( Mesh &mesh, ParallelData &pm, GlobalIndices &gI,
                              NodeVector &nodes2newNodes, EdgeVector &edges2newEdges, FaceVector &faces2newFaces )
    {                
        //const int dimworld = Mesh::dimension;
        const int NV = Mesh::verticesPerElement;
        
        nodes2newNodes.resize( mesh.c4n.size() );
        edges2newEdges.resize( mesh.edgeData.nEdges );
        faces2newFaces.resize( mesh.neighbors.nFaces );
        
        
        ManagedArray<int> idxA;
        
        gI.globalNodeIndex = gI.globalNodeIndex.sortByRow( idxA, nodes2newNodes );
        
        mesh.n4e.index2newIndex( nodes2newNodes );
        mesh.c4n.rearange( nodes2newNodes );
   /*     
        */
        //edges
        ManagedArray<int,2> edges( mesh.edgeData.nEdges );
        edges.fill(-1);
        for(int el = 0; el < (int)mesh.edgeData.element2edges.size(); el++)
            for(int k = 0, j = 0; k < (int)mesh.edgeData.element2edges.dimension(); k++) {
                const int e = mesh.edgeData.element2edges(el,k);
                if( edges(e,0) < 0 ) {
                    const int x = mesh.n4e(el, mesh.edgeData.referenceEdges[j++]);
                    const int y = mesh.n4e(el, mesh.edgeData.referenceEdges[j++]);
                    edges(e,0) = x;
                    edges(e,1) = y;
                } else {
                    j += 2;
                }
            }
        edges.sortRows();
        edges.sortByRow( idxA, edges2newEdges );
        
        mesh.edgeData.element2edges.index2newIndex( edges2newEdges );
        
        //faces
        ManagedArray<int,NV-1> faces( mesh.neighbors.nFaces, -1 );
        for(int el = 0; el < (int)mesh.n4e.size(); el++)
            for( int k = 0; k < NV; k++){
                const int f = mesh.neighbors.faceIndex(el,k);
                if( faces(f,0) < 0 )
                    for( int j = 0; j < NV - 1; j++ )
                        faces(f,j) = mesh.n4e(el, j < k ? j : j + 1 );
            }
       
        faces.sortRows();
        faces.sortByRow( idxA, faces2newFaces );  
            
        for( int el = 0; el < (int)mesh.n4e.size(); ++el )
            for( int k = 0; k < NV; ++k)
                mesh.neighbors.traits_(el,k).setFaceNumber( faces2newFaces[ mesh.neighbors.faceIndex(el,k) ] ); 
            
        pm.rearangeEntities( nodes2newNodes, edges2newEdges );
        
        //rI.rearangeEntities( nodes2newNodes, edges2newEdges );
    }        message[p] = false;
                continue;
            }

            int position = 0;
            int bufferSize = comm.template bufferSize<int>( 2 + buffer_vertex_size[p].size() + buffer_edge_size[p].size() 
                                                              + buffer_vertices[p].size() + buffer_edges[p].size()
                                                              + buffer_faces[p].size() );
            void * Buffer = std::malloc( bufferSize );
            
            int size_vertices = (int)buffer_vertices[p].size();
            int size_edges = (int)buffer_edges[p].size();
            
            comm.pack( &size_vertices, 1, Buffer, bufferSize, position ); 
            comm.pack( &size_edges, 1, Buffer, bufferSize, position ); 
            comm.pack( buffer_vertex_size[p].data()(), buffer_vertex_size[p].size(), Buffer, bufferSize, position);
            comm.pack( buffer_edge_size[p].data()(), buffer_edge_size[p].size(), Buffer, bufferSize, position);
            comm.pack( buffer_vertices[p].data()(), buffer_vertices[p].size(), Buffer, bufferSize, position);
            comm.pack( buffer_edges[p].data()(), buffer_edges[p].size(), Buffer, bufferSize, position);
            comm.pack( buffer_faces[p].data()(), buffer_faces[p].size(), Buffer, bufferSize, position);
            
            comm.Isend_packed( Buffer, bufferSize, p, 1, &requests[p] );
        }
        
        for(int p = 0; p < (int)comm.size(); p++){
            if( !message[p] )
                continue;
            int position = 0;
            comm.probe( p, 1, &status[p]);
            int bufferSize = (int)comm.getCount_packed( &status[p] );
        
            void * Buffer = std::malloc( bufferSize );
            
            comm.recv_packed( Buffer, bufferSize, p, 1 ); 
            
            int size_vertices;
            int size_edges;
            
            comm.unpack( &size_vertices, 1, Buffer, bufferSize, position );
            comm.unpack( &size_edges, 1, Buffer, bufferSize, position );
            comm.unpack( buffer_vertex_size[p].data()(), pm.rank2bdyNodes[p].size(), Buffer, bufferSize, position);
            comm.unpack( buffer_edge_size[p].data()(), pm.rank2bdyEdges[p].size(), Buffer, bufferSize, position);
            buffer_vertices[p].resize( size_vertices );
            buffer_edges[p].resize( size_edges);
            comm.unpack( buffer_vertices[p].data()(), size_vertices, Buffer, bufferSize, position);
            comm.unpack( buffer_edges[p].data()(), size_edges, Buffer, bufferSize, position);
            comm.unpack( buffer_faces[p].data()(), pm.rank2bdyFaces[p].size(), Buffer, bufferSize, position);
        }
        
        for( int p = 0; p < (int)comm.size(); p++ ){
            if( pm.rank2bdyNodes[p].size() > 0 ){
                for( int k = 0; k < (int)pm.rank2bdyNodes[p].size(); k++ )
                    for( int r = 0; r < buffer_vertex_size[p][k]; r++ ) {
                        if( buffer_vertices[p][r] != (int)comm.rank() )
                            pm.vertex2ranks[ pm.rank2bdyNodes[p][k] ].push_back( buffer_vertices[p][r] );
                    }
            }
            if( pm.rank2bdyEdges[p].size() > 0 ){
                for( int k = 0; k < (int)pm.rank2bdyEdges[p].size(); k++ )
                    for( int r = 0; r < buffer_edge_size[p][k]; r++ ) {
                        if( buffer_edges[p][r] != (int)comm.rank() )
                            pm.edge2ranks[ pm.rank2bdyEdges[p][k] ].push_back( buffer_edges[p][r] );
                    }
            }
        }
        
        for( int i = 0; i < (int)pm.vertex2ranks.size(); i++ )
            if( pm.vertex2ranks[i].size() )
                pm.vertex2ranks[i] = pm.vertex2ranks[i].unique();
            
        for( int i = 0; i < (int)pm.edge2ranks.size(); i++ )
            if( pm.edge2ranks[i].size() )
                pm.edge2ranks[i] = pm.edge2ranks[i].unique();
        
        //unpack new faces ranks //TODO
        ManagedArray<int> newFaceRanks( mesh.neighbors.nFaces, -1 );
        
        for( int p = 0; p < (int)comm.size(); p++ )
             for( int k = 0; k < (int)pm.rank2bdyFaces[p].size(); k++ )
                 newFaceRanks[ pm.rank2bdyFaces[p][k] ] = buffer_faces[p][k];    

        for( int el = 0; el < (int)mesh.n4e.size(); el++ )
             for( int k = 0; k < NV; k++ )
                 if( mesh.neighbors.isParallel(el,k) && newFaceRanks[ mesh.neighbors.faceIndex(el,k) ] > 0)
                     mesh.neighbors.traits_(el,k).setIndex( newFaceRanks[ mesh.neighbors.faceIndex(el,k) ] );

    }
#endif

} // namespace conformingsimplexgrid

#endif //MESH_PARTITION_HH
