#ifndef REFINE_PARALLEL_HH
#define REFINE_PARALLEL_HH

#include "../simplexGrid.hh"
#include "../mesh.hh"
#include "../parallel.hh"

namespace conformingsimplexgrid {
        
    template< class Vector > struct ProlongationData_pre_refine;
    template< class Vector > struct ProlongationData_refine;
    
    template< class Decomposition, class Grid, class ElementVector, class ... Manager >
    void refine3D_parallel( Decomposition &decomp, Grid &grid, ElementVector &marked, Manager & ... manager )
    {
        refine_parallel( decomp, grid, marked, manager ... );
    }
    
    template< class Decomposition, class Grid, class ElementVector, class ... Manager >
    void refine_parallel( Decomposition &decomp, Grid &grid, ElementVector &marked, Manager & ... manager )
    {
        
        grid.numNodes = grid.coordinates.size();
        grid.numEdges = grid.edges.size();
        grid.numFaces = grid.faces.size();
        
        const int NV = Grid::NV;
        
        ManagedArray<char> markedElements( grid.nElements(), 0 );
        ManagedArray<char> markedFaces( grid.nFaces(), 0 );
        ManagedArray<char> markedFaces_initial( grid.nFaces(), 0 );
        ManagedArray<char> markedEdges_initial( grid.nEdges() + grid.nFaces(), 0 );
        
        markEntities_parallel( decomp, grid, marked, markedFaces_initial , markedEdges_initial ); 
      
        int stillElementsToRefine = 1;
                
        while( stillElementsToRefine > 0 ) {
            ManagedArray<char> markedElements( grid.nElements(), 0 );
            ManagedArray<char> markedFaces( grid.nFaces(), 0 );
            ManagedArray<char> markedEdges( grid.nEdges(), 0 );

            for(size_t k = 0; k < markedEdges.size() & k < markedEdges_initial.size() ; k++)
                markedEdges[k] = markedEdges_initial[k];

            getMarkedEntitiesFromMarkedEdges(grid,markedElements,markedFaces,markedEdges ); 
            
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
            
            if( nMarked > 0 ){
                markedFaces.fill(0);
                markedEdges.fill(0);
                for(size_t el = 0; el < grid.nElements(); el++)
                    if( markedElements[el] > 0 )
                        for(size_t k = 1; k < NV-1; k++)
                            markedFaces[grid.element2faces(el,k)] = 1;
                
                size_t nMarkedFaces = 0;
                for(size_t el = 0; el < markedFaces.size(); el++)
                    if( markedFaces[el] > 0 )
                        nMarkedFaces++;
                    
                for(size_t f = 0; f < grid.nFaces(); f++)
                    if( markedFaces[f] > 0 )
                        for(size_t k = 1; k < NV-2; k++)
                            markedEdges[grid.face2edges(f,k)] = 1;
                        
                for(size_t k = 0; k < markedEdges.size(); k++)
                    if( markedEdges[k] > 0 )
                        markedEdges_initial[k] = 0;
                        
                ManagedArray<int> newNodesOnElements( grid.nElements(), -1 );
                ManagedArray<int> newNodesOnFaces( grid.nFaces(), -1 );
                ManagedArray<int> newNodesOnEdges( grid.nEdges(), -1 );
                            
                getNewNodes(grid,markedElements,markedFaces,markedEdges,newNodesOnElements,newNodesOnFaces,newNodesOnEdges);
                
                ProlongationData_pre_refine< ManagedArray<int> > pdata0( newNodesOnElements, newNodesOnFaces, newNodesOnEdges );
                pre_refine_Manager( pdata0, manager ... );
                
                const size_t nOld = grid.nElements();
                grid.resizeElements( nOld + nMarked );
                
                for( size_t el = 0; el < nOld; el++)
                    if( markedElements[el] > 0 )
                        grid.level[el] = lmin + 1;
                for(size_t k = nOld; k < nOld + nMarked; k++)
                    grid.level[k] = lmin + 1;
                
                const size_t nlevelFaces = grid.levelFaces.size();
    
                grid.levelFaces.resize( nlevelFaces + nMarkedFaces + nMarked );
                for( size_t el = 0, counter = nlevelFaces; el < nlevelFaces; el++)
                    if( markedFaces[el] > 0 )
                        grid.levelFaces[counter++] = grid.levelFaces[el];
                for( size_t el = nlevelFaces + nMarkedFaces; el < nlevelFaces + nMarkedFaces + nMarked; el++)
                    grid.levelFaces[el] = lmin+1; 
                
                ManagedArray<int> newInnerFaces( grid.nElements(), -1 );
                ManagedArray<int> newInnerEdges( grid.nFaces(), -1 );
                ManagedArray<int> left2rightElements( grid.nElements(), -1 );
                ManagedArray<int> left2rightFaces( grid.nFaces(), -1 );
                ManagedArray<int> left2rightEdges( grid.nEdges(), -1 );
                
                refineElements( grid, markedElements, newNodesOnElements, left2rightElements, nOld, lmin );
             
                refineEntity2SubEntities3( grid, markedElements, markedFaces, markedEdges, lmin, newInnerFaces,
                                           newInnerEdges, left2rightFaces, left2rightEdges );

                refineLeftChild( grid, markedElements, lmin );
                
                ProlongationData_refine<ManagedArray<int> > pdata( newNodesOnElements, newNodesOnFaces, newNodesOnEdges,
                    newInnerFaces, newInnerEdges, left2rightElements, left2rightFaces, left2rightEdges );
         
                decomp.refine( pdata );
                refine_Manager( pdata, manager ... );   
                                
                //mark edge of inner faces of marked Elements (x_0, x_1, x_2, x_3) of type 2, where the edge (x_1, x_2) is marked initially
                for( int el = 0; el < markedFaces_initial.size(); el++)
                    if( markedFaces_initial[el] > 0 && newInnerEdges[el] > 0){
                        markedFaces_initial[el] = 0;
                        markedEdges_initial[newInnerEdges[el]] = 1; 
                    }
            }

            stillElementsToRefine = 0;
            for(size_t el = 0; el < markedEdges_initial.size(); el++)
                if(markedEdges_initial[el] > 0){
                    stillElementsToRefine = 1;
                    break;
                }
        }
        getFacesForElements( grid.elements, grid.faces, grid.element2faces, grid.level );
        getFacesForElements( grid.faces, grid.edges, grid.face2edges );
    }
    
    template< class Decomposition, class Grid, class ElementVector, class ... Manager >
    void refine_parallel_markingTime( double &markingTime, Decomposition &decomp, Grid &grid, ElementVector &marked, Manager & ... manager )
    {
        
        grid.numNodes = grid.coordinates.size();
        grid.numEdges = grid.edges.size();
        grid.numFaces = grid.faces.size();
        
        const int NV = Grid::NV;
        
        ManagedArray<char> markedElements( grid.nElements(), 0 );
        ManagedArray<char> markedFaces( grid.nFaces(), 0 );
        ManagedArray<char> markedFaces_initial( grid.nFaces(), 0 );
        ManagedArray<char> markedEdges_initial( grid.nEdges() + grid.nFaces(), 0 );
        
        /////////////////
 
        double start_time = clock();
        markEntities_parallel( decomp, grid, marked, markedFaces_initial , markedEdges_initial ); 
        markingTime = (float)(clock() - start_time) / CLOCKS_PER_SEC;
        //////////////////
        
        int stillElementsToRefine = 1;
                
        while( stillElementsToRefine > 0 ) {
            ManagedArray<char> markedElements( grid.nElements(), 0 );
            ManagedArray<char> markedFaces( grid.nFaces(), 0 );
            ManagedArray<char> markedEdges( grid.nEdges(), 0 );

            for(size_t k = 0; k < markedEdges.size() & k < markedEdges_initial.size() ; k++)
                markedEdges[k] = markedEdges_initial[k];

            getMarkedEntitiesFromMarkedEdges(grid,markedElements,markedFaces,markedEdges ); 
            
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
            
            if( nMarked > 0 ){
                
                
                markedFaces.fill(0);
                markedEdges.fill(0);
                for(size_t el = 0; el < grid.nElements(); el++)
                    if( markedElements[el] > 0 )
                        for(size_t k = 1; k < NV-1; k++)
                            markedFaces[grid.element2faces(el,k)] = 1;
                
                size_t nMarkedFaces = 0;
                for(size_t el = 0; el < markedFaces.size(); el++)
                    if( markedFaces[el] > 0 )
                        nMarkedFaces++;
                    
                for(size_t f = 0; f < grid.nFaces(); f++)
                    if( markedFaces[f] > 0 )
                        for(size_t k = 1; k < NV-2; k++)
                            markedEdges[grid.face2edges(f,k)] = 1;
                        
                for(size_t k = 0; k < markedEdges.size(); k++)
                    if( markedEdges[k] > 0 )
                        markedEdges_initial[k] = 0;
                        
                ManagedArray<int> newNodesOnElements( grid.nElements(), -1 );
                ManagedArray<int> newNodesOnFaces( grid.nFaces(), -1 );
                ManagedArray<int> newNodesOnEdges( grid.nEdges(), -1 );
                            
                getNewNodes(grid,markedElements,markedFaces,markedEdges,newNodesOnElements,newNodesOnFaces,newNodesOnEdges);
                
                ProlongationData_pre_refine< ManagedArray<int> > pdata0( newNodesOnElements, newNodesOnFaces, newNodesOnEdges );
                pre_refine_Manager( pdata0, manager ... );
                
                const size_t nOld = grid.nElements();
                grid.resizeElements( nOld + nMarked );
                
                for( size_t el = 0; el < nOld; el++)
                    if( markedElements[el] > 0 )
                        grid.level[el] = lmin + 1;
                for(size_t k = nOld; k < nOld + nMarked; k++)
                    grid.level[k] = lmin + 1;
                
                const size_t nlevelFaces = grid.levelFaces.size();
    
                grid.levelFaces.resize( nlevelFaces + nMarkedFaces + nMarked );
                for( size_t el = 0, counter = nlevelFaces; el < nlevelFaces; el++)
                    if( markedFaces[el] > 0 )
                        grid.levelFaces[counter++] = grid.levelFaces[el];
                for( size_t el = nlevelFaces + nMarkedFaces; el < nlevelFaces + nMarkedFaces + nMarked; el++)
                    grid.levelFaces[el] = lmin+1; 
                
                ManagedArray<int> newInnerFaces( grid.nElements(), -1 );
                ManagedArray<int> newInnerEdges( grid.nFaces(), -1 );
                ManagedArray<int> left2rightElements( grid.nElements(), -1 );
                ManagedArray<int> left2rightFaces( grid.nFaces(), -1 );
                ManagedArray<int> left2rightEdges( grid.nEdges(), -1 );
                
                
                refineElements( grid, markedElements, newNodesOnElements, left2rightElements, nOld, lmin );
             
              
                refineEntity2SubEntities3( grid, markedElements, markedFaces, markedEdges, lmin, newInnerFaces,
                                           newInnerEdges, left2rightFaces, left2rightEdges );
  
                refineLeftChild( grid, markedElements, lmin );
           
                
                
                ProlongationData_refine<ManagedArray<int> > pdata( newNodesOnElements, newNodesOnFaces, newNodesOnEdges,
                    newInnerFaces, newInnerEdges, left2rightElements, left2rightFaces, left2rightEdges );
                
                //double start_time = clock();
                decomp.refine( pdata );
               // double time = (float)(clock() - start_time) / CLOCKS_PER_SEC;
              //  printf("rank = %d, refineDecomposition = %lf\n",decomp.rank(),time); 

                refine_Manager( pdata, manager ... );                 
                
                
                                
                //mark edge of inner faces of marked Elements (x_0, x_1, x_2, x_3) of type 2, where the edge (x_1, x_2) is marked initially
                for( int el = 0; el < markedFaces_initial.size(); el++)
                    if( markedFaces_initial[el] > 0 && newInnerEdges[el] > 0){
                        markedFaces_initial[el] = 0;
                        markedEdges_initial[newInnerEdges[el]] = 1; 
                    }
            }

            stillElementsToRefine = 0;
            for(size_t el = 0; el < markedEdges_initial.size(); el++)
                if(markedEdges_initial[el] > 0){
                    stillElementsToRefine = 1;
                    break;
                }
        }
        getFacesForElements( grid.elements, grid.faces, grid.element2faces, grid.level );
        getFacesForElements( grid.faces, grid.edges, grid.face2edges );
    }
    
    
    template< class Decomposition, class Grid, class ElementVector, class Vector >
    void markEntities_parallel( Decomposition &decomp, Grid &grid, ElementVector &markedElements, Vector &markedFaces_initial , Vector &markedEdges_initial )
    {
        const int NV = Grid::NV;
        
        int nMarkedEdges = 0;
        int nMarkedEdges_old = -1;
        ManagedArray<char> markedFaces( grid.nFaces(), 0 );
        
        const int preMarking = 5;
        
        for( int n = 0; n < preMarking; n++ ) {
            markEntities( grid, markedElements, markedFaces, markedEdges_initial );
            decomp.sumEdgeVector( markedEdges_initial );
            for( int k = 0; k < markedEdges_initial.size(); k++)
                markedEdges_initial[k] > 0 ? markedEdges_initial[k] = 1 : markedEdges_initial[k] = 0;
            getMarkedEntitiesFromMarkedEdges(grid,markedElements,markedFaces,markedEdges_initial );
        }
        for(size_t k = 0; k < markedEdges_initial.size(); ++k)
            if(markedEdges_initial[k] > 0)
                nMarkedEdges++;
            
        int iter = 0; 
        int stillEdgesToMark = 1;
        while( stillEdgesToMark > 0 ) {
            markEntities( grid, markedElements, markedFaces, markedEdges_initial );
            decomp.sumEdgeVector( markedEdges_initial );
            for( int k = 0; k < markedEdges_initial.size(); k++)
                markedEdges_initial[k] > 0 ? markedEdges_initial[k] = 1 : markedEdges_initial[k] = 0;
            getMarkedEntitiesFromMarkedEdges(grid,markedElements,markedFaces,markedEdges_initial );
            nMarkedEdges_old = nMarkedEdges;
            nMarkedEdges = 0;
            for(size_t k = 0; k < markedEdges_initial.size(); ++k)
                if(markedEdges_initial[k] > 0)
                    nMarkedEdges++;
            stillEdgesToMark = 1;
            if( nMarkedEdges == nMarkedEdges_old ) stillEdgesToMark = 0;
            decomp.scalar_sum( stillEdgesToMark );
            if( stillEdgesToMark == 0) break;
            printf("reiteration needed: iter = %d, stillEdgesToMark = %d, rank = %d\n",++iter,stillEdgesToMark,decomp.rank());
        }  
        
        markEntities( grid, markedElements, markedFaces, markedEdges_initial );
                        
        //mark 1 face of marked Elements (x_0, x_1, x_2, x_3) of type 2, where the edge (x_1, x_2) is marked
        for( int el = 0; el < markedElements.size(); el++)
            if( markedElements[el] > 0 && (grid.level[el] % (NV-1)) == 2 && markedEdges_initial[grid.face2edges(grid.element2faces(el,0),2)] > 0 )
                markedFaces_initial[grid.element2faces(el,2)] = 1;    

        decomp.sumFaceVector( markedFaces_initial );    
        for( int k = 0; k < markedFaces_initial.size(); k++)
                markedFaces_initial[k] > 0 ? markedFaces_initial[k] = 1 : markedFaces_initial[k] = 0;
    }
    
    template< class Grid, class Decomposition, class Vector >
    bool testMarkedEntities(Grid &grid, Decomposition &decomp, Vector &markedFaces, Vector &markedEdges) {
        bool returnValue = true;
        ManagedArray<int> markedFacesOnRank1(decomp.size(), 0);
        ManagedArray<int> markedEdgesOnRank1(decomp.size(), 0);
        ManagedArray<int> markedFacesOnRank2(decomp.size(), 0);
        ManagedArray<int> markedEdgesOnRank2(decomp.size(), 0);

        for( int k = 0; k < decomp.edge2ranks.size(); k++)
            if( markedEdges[k] > 0 )
                for( int r = 0; r < decomp.edge2ranks[k].size(); r++){
                    const int p = decomp.edge2ranks[k][r];
                    markedEdgesOnRank1[p]++;
                }
        for( int k = 0; k < decomp.face2ranks.size(); k++)
            if( markedFaces[k] > 0 )
                for( int r = 0; r < decomp.face2ranks[k].size(); r++){
                    const int p = decomp.face2ranks[k][r];
                    markedFacesOnRank1[p]++;
                }
        
        decomp.all2all( markedFacesOnRank1.data(), markedFacesOnRank2.data(), 1 );
        decomp.all2all( markedEdgesOnRank1.data(), markedEdgesOnRank2.data(), 1 );
        
        for( int p = 0; p < decomp.size(); p++ ){
            if( p == decomp.rank() )
                continue;
            if( !(markedFacesOnRank1[p] == markedFacesOnRank2[p])){
                printf("Faces, %d on rank %d, %d -- %d \n", p, decomp.rank(), markedFacesOnRank1[p], markedFacesOnRank2[p]);
                returnValue = false;
            }
        }
        
        for( int p = 0; p < decomp.size(); p++ ){
            if( p == decomp.rank() )
                continue;
            if( !(markedEdgesOnRank1[p] == markedEdgesOnRank2[p])){
                printf("Edges, %d on rank %d, %d -- %d \n", p, decomp.rank(), markedEdgesOnRank1[p], markedEdgesOnRank2[p]);
                returnValue = false;
            }
        }
        
        return returnValue;
        
    }

} // namespace conformingsimplexgrid

#endif //REFINE_PARALLEL_HH