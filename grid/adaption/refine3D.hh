/** 
@file refine3D.hh
@brief implementation of the recursive refinement algorithm for a tetrahedral mesh. All marked elements will be refined at least one time. To ensure conformity further elements are refined during the recursive process.
*/

#ifndef REFINE_HH
#define REFINE_HH

#include "../array/managed_array.hh"
#include "../auxiliary/element2subentity.hh"

namespace conformingsimplexgrid {
    
    template< class Vector > 
    struct ProlongationData_pre_refine
    {
        ProlongationData_pre_refine( Vector & element2newNode_in, Vector & face2newNode_in, Vector & edge2newNode_in )
        : element2newNode(element2newNode_in), face2newNode(face2newNode_in), edge2newNode(edge2newNode_in)
        {}
        
        Vector &element2newNode;
        Vector &face2newNode;
        Vector &edge2newNode;
    };
    
    template< class Vector > 
    struct ProlongationData_refine
    {
        ProlongationData_refine( Vector & element2newNode_in, Vector & face2newNode_in, Vector & edge2newNode_in,
            Vector & element2innerFace_in, Vector & face2innerEdge_in,
            Vector & left2rightElements_in, Vector & left2rightFaces_in, Vector & left2rightEdges_in )
        : element2newNode(element2newNode_in), face2newNode(face2newNode_in), edge2newNode(edge2newNode_in),
        element2innerFace(element2innerFace_in), face2innerEdge(face2innerEdge_in), left2rightElements(left2rightElements_in),
        left2rightFaces(left2rightFaces_in), left2rightEdges(left2rightEdges_in)
        {}
        
        Vector &element2newNode;
        Vector &face2newNode;
        Vector &edge2newNode;
        Vector &element2innerFace;
        Vector &face2innerEdge;
        Vector &left2rightElements;
        Vector &left2rightFaces;
        Vector &left2rightEdges;
    };
    
    template< std::size_t NV > class ReferenceSimplex;
    
    template< size_t NV >
    class refine_scheme {
    public:
        std::array< int, NV > idx_left;
        std::array< int, NV > idx_right;
        int type;
        explicit refine_scheme ( int level ) {
            type = (level % (NV-1));
            
            idx_left[1] = -1;
            idx_left[0] = 0;
            for(size_t k = 2; k < NV; ++k)
                idx_left[k] = k - 1;
            
            idx_right[0] = type + 1; idx_right[1] = -1;
            for(size_t k = 1, ctr = 2; k < NV; ++k )
                if( size_t(type + 1) != k) idx_right[ctr++] = k;
        }
    };
    
    template< class Grid, class ElementVector, class ... Manager >
    void refine3D( Grid &grid, ElementVector &marked, Manager & ... manager ) \
    /// Call this function the refine all elements k where marked[k] > 0. The recursive refinement algorithm will refine further elements to ensure conformity.
    {
        refine( grid, marked, manager ... );
    }
    
    
    template< class Grid, class ElementVector, class ... Manager >
    void refine( Grid &grid, ElementVector &marked, Manager & ... manager )
    {
        const int NV = Grid::NV;
                
        grid.numNodes = grid.coordinates.size();
        grid.numEdges = grid.edges.size();
        grid.numFaces = grid.faces.size();
        
        bool stillElementsToRefine = true;
        
        while( stillElementsToRefine ) {
            ManagedArray<int> markedElements( grid.nElements(), 0 );
            ManagedArray<int> markedFaces( grid.nFaces(), 0 );
            ManagedArray<int> markedEdges( grid.nEdges(), 0 );
            
            for(size_t el = 0; el < marked.size(); el++)
                markedElements[el] = marked[el];
            
            markEntities( grid, markedElements, markedFaces, markedEdges );
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
       
            markedFaces.fill( 0 );
            markedEdges.fill( 0 );
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
                    
            ManagedArray<int> newNodesOnElements( grid.nElements(), -1 );
            ManagedArray<int> newNodesOnFaces( grid.nFaces(), -1 );
            ManagedArray<int> newNodesOnEdges( grid.nEdges(), -1 );
            
            getNewNodes(grid,markedElements,markedFaces,markedEdges,newNodesOnElements,newNodesOnFaces,newNodesOnEdges);

            ProlongationData_pre_refine<ManagedArray<int> > pdata0( newNodesOnElements, newNodesOnFaces, newNodesOnEdges );
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
           
            refineEntity2SubEntities3( grid, markedElements, markedFaces, markedEdges, lmin, newInnerFaces, newInnerEdges, left2rightFaces, left2rightEdges );
            
            ProlongationData_refine<ManagedArray<int> > pdata( newNodesOnElements, newNodesOnFaces, newNodesOnEdges,
                newInnerFaces, newInnerEdges, left2rightElements, left2rightFaces, left2rightEdges );
            refine_Manager( pdata, manager ... );
            
            refineLeftChild( grid, markedElements, lmin );

            stillElementsToRefine = false;
            for(size_t el = 0; el < marked.size(); el++)
                if(marked[el] > 0){
                    stillElementsToRefine = true;
                    break;
                }
        }
        getFacesForElements( grid.elements, grid.faces, grid.element2faces, grid.level );
        getFacesForElements( grid.faces, grid.edges, grid.face2edges );
    }

    template< class Grid, class ElementVector, class FaceVector, class EdgeVector >
    void markEntities( Grid &grid, ElementVector &markedElements, FaceVector & markedFaces, EdgeVector &markedEdges )
    ///Function called by the recursive refinement routine to mark elements, faces, edges for refinement in order to conformingly refine elements indicated by markedElements.
    {
        const int NV = Grid::NV;

        int nMarkedElements = 0;
        for(size_t el = 0; el < grid.nElements(); el++)
            if( markedElements[el] > 0 )
                nMarkedElements++;
        
        int stillElementsMarked = 1;
        while( stillElementsMarked > 0 ) {
            for(size_t el = 0; el < grid.nElements(); el++)
                if( markedElements[el] > 0 ){
                    for(size_t k = 1; k < NV-1; k++)
                        markedFaces[grid.element2faces(el,k)] = 1;
                }
            for(size_t f = 0; f < grid.nFaces(); f++)
                if( markedFaces[f] > 0 )
                    for(size_t k = 1; k < NV-2; k++)
                        markedEdges[grid.face2edges(f,k)] = 1;
            
            getMarkedEntitiesFromMarkedEdges(grid,markedElements,markedFaces,markedEdges );

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
                    markedFaces[grid.element2faces(el,k)] = 1;
        for(size_t f = 0; f < grid.nFaces(); f++)
            if( markedFaces[f] > 0 )
                for(size_t k = 1; k < NV-2; k++)
                    markedEdges[grid.face2edges(f,k)] = 1;       
    }
    
    template< class Grid, class ElementVector, class FaceVector, class EdgeVector >
    void getMarkedEntitiesFromMarkedEdges( Grid &grid, ElementVector &markedElements, FaceVector & markedFaces, EdgeVector &markedEdges )
    ///Function called by the recursive refinement routine to get generate the markedFaces and markedEdges information from marked Elements.
    {
        const int NV = Grid::NV;
        
        ManagedArray<char, NV> element2markedFaces(grid.nElements(),0);
        ManagedArray<char, NV-1> faces2markedEdges(grid.face2edges.size(),0);
        
        for(size_t f = 0; f < grid.face2edges.size(); f++)
            for(size_t k = 0; k < NV-1; k++)
                faces2markedEdges(f,k) = markedEdges[ grid.face2edges(f,k)];
                
        for(size_t f = 0; f < grid.face2edges.size(); f++)
            for(size_t k = 0; k < NV-1; k++)
                if(faces2markedEdges(f,k) > 0){
                    markedFaces[f] = 1;
                    break;
                }

        for(size_t el = 0; el < grid.nElements(); el++)
            for(size_t k = 0; k < NV; k++)
                element2markedFaces(el,k) = markedFaces[ grid.element2faces(el,k) ];
                    
        for(size_t el = 0; el < grid.nElements(); el++)
            for(size_t k = 0; k < NV; k++)
                if(element2markedFaces(el,k) > 0){
                    markedElements[el] = 1;
                    break;
                }
    }
    
    template< class Grid, class ElementVector, class FaceVector, class EdgeVector, class Vector >
    void getNewNodes( Grid &grid, ElementVector &markedElements, FaceVector &markedFaces, EdgeVector &markedEdges,
                      Vector &newNodesOnElements, Vector &newNodesOnFaces, Vector &newNodesOnEdges  )
    ///Function called by the recursive refinement routine to generate new Node indices and coordinates.
    {
        const int DIM = Grid::DIM;
        const int NV = Grid::NV;
        
        int nC = grid.nNodes();
        int nNewNodes = 0;
        for(size_t i = 0; i < grid.nEdges(); i++)
            if( markedEdges[i] > 0 )
                newNodesOnEdges[i] = nC + nNewNodes++;
                            
        for(size_t i = 0; i < grid.nFaces(); i++)
            if( markedFaces[i] > 0 )
                newNodesOnFaces[i] = newNodesOnEdges[grid.face2edges(i,1)];
                 
       for(size_t i = 0; i < grid.nElements(); i++)
            if( markedElements[i] > 0 )
                newNodesOnElements[i] = newNodesOnFaces[grid.element2faces(i,1)];
            
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
    void refineElements( Grid &grid, ElementVector &marked, NodeVector &newNodesOnElements, Vector & left2rightElements,
                         int nOld, int lmin )
    ///Function called by the recursive refinement routine to refine markedElements via the refinement scheme.
    {
        const int NV = Grid::NV;
                
        refine_scheme<NV> scm( lmin );
        
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
    
    template< class Grid, class ElementVector, class FaceVector, class EdgeVector, class Vector >
    void refineEntity2SubEntities3( Grid &grid, ElementVector &markedElements, FaceVector &markedFaces, EdgeVector &markedEdges, int lmin, 
                                    Vector & newInnerFaces, Vector & newInnerEdges, Vector &right_faces, Vector &right_edges)
    ///Function called by the recursive refinement routine to refine element2faces and faces2edges structure.
    {
        
        const int NV = Grid::NV;
        assert(NV == 4);
        refine_scheme<NV> scm( lmin );
        int type = scm.type;
        
        int nMarkedEdges = 0;
        for( size_t k = 0; k < markedEdges.size(); ++k)
            if( markedEdges[k] > 0 )
                right_edges[k] = grid.nEdges() + nMarkedEdges++;

        int nMarkedFaces = 0;
        for( size_t k = 0; k < markedFaces.size(); ++k)
            if( markedFaces[k] > 0 )
                right_faces[k] = grid.nFaces() + nMarkedFaces++;
         
        int nMarkedElements = 0;   
        for( size_t k = 0; k < markedElements.size(); ++k)
            if( markedElements[k] > 0 )
                nMarkedElements++;
        
        ManagedArray<int,3> innerFaces2edges( nMarkedElements, 0 );
        
        int nElementsOld = grid.nElements() - nMarkedElements;
        int nFacesOld = grid.nFaces();
        int nEdgesOld = grid.nEdges();
        
        int nFacesNew = grid.nFaces() + nMarkedFaces + nMarkedElements;
        
        int count_faces = nFacesOld + nMarkedFaces;
        for( size_t el = 0; el < markedElements.size(); el++ )
            if( markedElements[el] > 0 )
                newInnerFaces[el] = count_faces++; 
    
        grid.resizeFaces( count_faces );    
            
        int count_edges = nEdgesOld + nMarkedEdges;
        for( size_t el = 0; el < markedFaces.size(); el++ )
            if( markedFaces[el] > 0 )
                newInnerEdges[el] = count_edges++;
        
        grid.resizeEdges( count_edges );
                
        //refine faces2edges information for right children simplices
        int counter = nFacesOld;
        switch (type) {
            case 0:
                for(size_t el = 0; el < markedFaces.size(); ++el){
                    if(markedFaces[el] > 0){ 
                        grid.face2edges(counter,0) = newInnerEdges[el];
                        grid.face2edges(counter,1) = grid.face2edges(el,2);
                        grid.face2edges(counter,2) = right_edges[grid.face2edges(el,1)];
                    counter++;
                    }
                }
                for(size_t el = 0; el < markedElements.size(); ++el){
                    if(markedElements[el] > 0){ 
                        grid.face2edges(counter,0) = newInnerEdges[ grid.element2faces(el,1) ];
                        grid.face2edges(counter,1) = grid.face2edges( grid.element2faces(el,0) ,2);
                        grid.face2edges(counter,2) = newInnerEdges[ grid.element2faces(el,2) ];
                    counter++;
                    }
                }
            break;
            
            case 1:
            {
                ManagedArray<int> idx_12( nFacesNew, 0 );
                for(size_t el = 0; el < markedElements.size(); ++el){
                    if(markedElements[el] > 0){
                        idx_12[grid.element2faces(el,1)] = 1;
                        idx_12[grid.element2faces(el,2)] = 2;
                    }
                }
                
                for(size_t el = 0; el < markedFaces.size(); ++el){
                    if(markedFaces[el] > 0){
                        assert( idx_12[el] > 0 );
                        if( idx_12[el] == 1 ){
                            grid.face2edges(counter,0) = right_edges[grid.face2edges(el,1)];
                            grid.face2edges(counter,1) = grid.face2edges(el,2);
                            grid.face2edges(counter,2) = newInnerEdges[el];
                        } else {
                            grid.face2edges(counter,0) = newInnerEdges[el];
                            grid.face2edges(counter,1) = grid.face2edges(el,2);
                            grid.face2edges(counter,2) = right_edges[grid.face2edges(el,1)];
                        }
                        counter++;
                    }
                }
                
                for(size_t el = 0; el < markedElements.size(); ++el){
                    if(markedElements[el] > 0){ 
                        grid.face2edges(counter,0) = newInnerEdges[ grid.element2faces(el,2) ];
                        grid.face2edges(counter,1) = grid.face2edges( grid.element2faces(el,0) ,2);
                        grid.face2edges(counter,2) = newInnerEdges[ grid.element2faces(el,1) ];
                    counter++;
                    }
                }
            }
            break;
            
            case 2:
                for(size_t el = 0; el < markedFaces.size(); ++el){
                    if(markedFaces[el] > 0){
                        grid.face2edges(counter,0) = right_edges[grid.face2edges(el,1)];
                        grid.face2edges(counter,1) = grid.face2edges(el,2);
                        grid.face2edges(counter,2) = newInnerEdges[el];
                    counter++;
                    }
                }
                
                for(size_t el = 0; el < markedElements.size(); ++el){
                    if(markedElements[el] > 0){ 
                        grid.face2edges(counter,0) = newInnerEdges[ grid.element2faces(el,1) ];
                        grid.face2edges(counter,1) = newInnerEdges[ grid.element2faces(el,2) ];
                        grid.face2edges(counter,2) = grid.face2edges( grid.element2faces(el,0) ,2);
                    counter++;
                    }
                }
            break;
        }
        
        //refine faces2edges information for left children simplices
        for(size_t el = 0; el < markedFaces.size(); ++el){
            if(markedFaces[el] > 0){
                const int tmp = grid.face2edges(el,0);
                grid.face2edges(el,0) = grid.face2edges(el,1);
                grid.face2edges(el,1) = tmp;
                grid.face2edges(el,2) = newInnerEdges[el];
            }
        }
        
        counter = nElementsOld;
        switch (type) {
            case 0:
                for(size_t el = 0; el < markedElements.size(); ++el){
                    if(markedElements[el] > 0){ 
                        grid.element2faces(counter,0) = newInnerFaces[el];
                        grid.element2faces(counter,1) = right_faces[grid.element2faces(el,1)];
                        grid.element2faces(counter,2) = grid.element2faces(el,3);
                        grid.element2faces(counter,3) = right_faces[grid.element2faces(el,2)];
                    counter++;
                    }
                }
                
            break;
            
            case 1:
                for(size_t el = 0; el < markedElements.size(); ++el){
                    if(markedElements[el] > 0){ 
                        grid.element2faces(counter,0) = newInnerFaces[el];
                        grid.element2faces(counter,1) = right_faces[grid.element2faces(el,2)];
                        grid.element2faces(counter,2) = grid.element2faces(el,3);
                        grid.element2faces(counter,3) = right_faces[grid.element2faces(el,1)];
                    counter++;
                    }
                }
                
            break;
                
            case 2:
                for(size_t el = 0; el < markedElements.size(); ++el){
                    if(markedElements[el] > 0){ 
                        grid.element2faces(counter,0) = right_faces[grid.element2faces(el,1)];
                        grid.element2faces(counter,1) = right_faces[grid.element2faces(el,2)];
                        grid.element2faces(counter,2) = grid.element2faces(el,3);
                        grid.element2faces(counter,3) = newInnerFaces[el];
                    counter++;
                    }
                }
            break;    
        }
        
        for(size_t el = 0; el < markedElements.size(); ++el){
            if(markedElements[el] > 0){
                const int tmp = grid.element2faces(el,0);
                grid.element2faces(el,0) = grid.element2faces(el,1);
                grid.element2faces(el,1) = grid.element2faces(el,2);
                grid.element2faces(el,2) = tmp;
                grid.element2faces(el,3) = newInnerFaces[el];
            }
        }
          
    }
    
    template< class Grid, class ElementVector >
    void refineLeftChild( Grid &grid, ElementVector &markedElements, int lmin )
    {
        const int NV = Grid::NV;
        assert(NV == 4);
        refine_scheme<NV> scm( lmin );
        int type = scm.type;
                
        switch(type){
            case 0:
               for(size_t el = 0, counter = markedElements.size(); el < markedElements.size(); ++el){
                    if(markedElements[el] > 0){
                        grid.rightChild(counter,0) = 1;
                        grid.rightChild(counter,1) = grid.rightChild(el,1);
                        grid.rightChild(counter,2) = grid.rightChild(el,3);
                        grid.rightChild(counter,3) = grid.rightChild(el,2);
                    counter++;
                    }
                }
            break;
            
            case 1:
               for(size_t el = 0, counter = markedElements.size(); el < markedElements.size(); ++el){
                    if(markedElements[el] > 0){
                        grid.rightChild(counter,0) = 1;
                        grid.rightChild(counter,1) = grid.rightChild(el,2);
                        grid.rightChild(counter,2) = grid.rightChild(el,3);
                        grid.rightChild(counter,3) = grid.rightChild(el,1);
                    counter++;
                    }
                }
            break;
            
            case 2:
                for(size_t el = 0, counter = markedElements.size(); el < markedElements.size(); ++el){
                    if(markedElements[el] > 0){
                        grid.rightChild(counter,0) = grid.rightChild(el,1);
                        grid.rightChild(counter,1) = grid.rightChild(el,2);
                        grid.rightChild(counter,2) = grid.rightChild(el,3);
                        grid.rightChild(counter,3) = 1;
                    counter++;
                    }
                }
            break;
        }
        for(size_t el = 0; el < markedElements.size(); ++el){
            if(markedElements[el] > 0){ 
                const int tmp0 = grid.rightChild(el,0);
                grid.rightChild(el,0) = grid.rightChild(el,1);
                grid.rightChild(el,1) = grid.rightChild(el,2);
                grid.rightChild(el,2) = tmp0;
                grid.rightChild(el,3) = 0;
            }
        }
    }

} // namespace conformingsimplexgrid

#endif //REFINE_HH