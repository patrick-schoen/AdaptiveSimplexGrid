/** 
@file mesh.hh
@brief this header file contains all required headers to perform seriell mesh refinement. Additionally it defines a 2D and a 3D grid class, containing index vectors for elements, faces and edges, that define the adaptive mesh structure.
*/

#ifndef MESH_HH
#define MESH_HH

#include "simplexGrid.hh"

#include <vector>
#include <array>
#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <math.h>
#include <stdio.h>

#include "array/array_lib.hh"
#include "auxiliary/auxiliary.hh"
#include "adaption/adaption_seriell.hh"


namespace conformingsimplexgrid {
    
    class Grid_3D /// mesh structure of a 3-dimensional grid of tetrahedron.
    {        
    public:
        static const int dimensionworld = 3;
        static const int verticesPerElement = 4;
        static const int NV = verticesPerElement;
        static const int dim = dimensionworld;
        static const int DIM = dimensionworld;

        ManagedArray<int,NV> elements; ///<each row contains the nodes indices of one element.
        ManagedArray<int,NV> element2faces; ///<each row contains the faces indices of one element.
        ManagedArray<int,NV-1> face2edges; ///<each row contains the edge indices of one face.
        
        ManagedArray<unsigned char,1> level; ///<level[k] gives the level of the k-th element.
        ManagedArray<unsigned char,NV> rightChild; ///<rightChild(el,k) is 1 iff the the element el was a rightChild upon creating the face k.
        ManagedArray<unsigned char,1> levelFaces; ///<levelFaces[k] gives the level of the k-th face.
 
        ManagedArray<int,NV-1> faces; ///<each row contains the nodes indices of one face.
        ManagedArray<int,NV-2> edges; ///<each row contains the nodes indices of one edge.
        
        std::size_t numFaces, numEdges, numNodes;
        
        ManagedArray<double,DIM> coordinates; ///<each row contains contains the coordinates of one node.
        
        size_t nNodes() const  { return coordinates.size(); } 
        ///<returns the number of nodes of the mesh.
       
        size_t nEdges() const  { return edges.size(); } 
        ///<returns the edges of nodes of the mesh.
       
        size_t nFaces() const  { return faces.size(); }
        ///<returns the faces of nodes of the mesh.
        
        size_t nElements() const { return elements.size(); } 
        ///<returns the elements of nodes of the mesh.
        
        int& element2edges( std::size_t col, std::size_t row ) {return face2edges[element2faces[col,(int)row/4],row % 3];}
        const int& element2edges( std::size_t col, std::size_t row ) const {return face2edges[element2faces[col,(int)row/4],row % 3];}    
        
        void initialTriangulation() ///initiates the mesh arrays, when only elements and coordinates are given. Thus it initializes level, rightChild, faces, edges, element2faces, faces2edges, levelFaces.
        {
            level.resize( elements.size() );
            rightChild.resize( elements.size() );
            level.fill(0);
            std::fill(rightChild.begin(),rightChild.end(),false); 
            getElement2Faces ( elements , faces, element2faces );
            resizeFaces( faces.size() );
            getElement2Faces ( faces , edges, face2edges );
            resizeEdges( edges.size() );    
            levelFaces.fill(0);
            numNodes = coordinates.size();
        }
        
        template< class Communicator >
        void pack( Communicator &comm, void *buffer, std::size_t bufferSize, int &position ) {    
            if( nElements() > 0 ) {
                std::array<int, 4> mesh_sizes;
                
                mesh_sizes[0] = (int)elements.size();
                mesh_sizes[1] = (int)face2edges.size();
                mesh_sizes[2] = (int)edges.size();
                mesh_sizes[3] = (int)coordinates.size();
                
                comm.pack( mesh_sizes.data(), mesh_sizes.size(), buffer, bufferSize, position );
                comm.pack( elements.data(), elements.size() * NV, buffer, bufferSize, position );
                comm.pack( level.data(), level.size(), buffer, bufferSize, position );
                comm.pack( levelFaces.data(), levelFaces.size(), buffer, bufferSize, position );
                comm.pack( rightChild.data(), rightChild.size() * NV, buffer, bufferSize, position );
                comm.pack( element2faces.data(), element2faces.size() * NV, buffer, bufferSize, position );
                comm.pack( face2edges.data(), face2edges.size() * (NV-1), buffer, bufferSize, position );
                comm.pack( coordinates.data(), coordinates.size() * dim, buffer, bufferSize, position );
            } else {
                int zero_value = 0;
                comm.pack( &zero_value, 1, buffer, bufferSize, position); 
            }
        }
        
        template< class Communicator >
        void unpack( Communicator &comm, void *buffer, std::size_t bufferSize, int &position ) {
            std::array<int, 4> mesh_sizes = {{0,0,0,0}};
            comm.unpack( mesh_sizes.data(), mesh_sizes.size(), buffer, bufferSize, position );
            if( mesh_sizes[0] > 0 ) {
                resizeElements(mesh_sizes[0]);
                resizeFaces(mesh_sizes[1]);
                resizeEdges(mesh_sizes[2]);
                resizeNodes(mesh_sizes[3]);

                comm.unpack( elements.data(), elements.size() * NV, buffer, bufferSize, position );
                comm.unpack( level.data(), level.size(), buffer, bufferSize, position );
                comm.unpack( levelFaces.data(), levelFaces.size(), buffer, bufferSize, position );
                comm.unpack( rightChild.data(), rightChild.size() * NV, buffer, bufferSize, position );
                comm.unpack( element2faces.data(), element2faces.size() * NV, buffer, bufferSize, position );
                comm.unpack( face2edges.data(), face2edges.size() * (NV-1), buffer, bufferSize, position );
                comm.unpack( coordinates.data(), coordinates.size() * dim, buffer, bufferSize, position );
            }   
        }
        
        template< class Communicator >
        std::size_t bufferSize( Communicator &comm ) {
            std::size_t bufferSize = 0;
            if( nElements() > 0 ) {
                bufferSize += comm.template bufferSize<int>( 4 );
                bufferSize += comm.template bufferSize<int>( nElements() * NV );
                bufferSize += comm.template bufferSize<unsigned char>( nElements() );
                bufferSize += comm.template bufferSize<unsigned char>( nFaces() );
                bufferSize += comm.template bufferSize<unsigned char>( nElements() * NV );
                bufferSize += comm.template bufferSize<int>( nElements() * NV );
                bufferSize += comm.template bufferSize<int>( nFaces() * (NV-1) );
                bufferSize += comm.template bufferSize<double>( nNodes() * dim );
            } else {
                bufferSize += comm.template bufferSize<int>( 1 );
            }
            return bufferSize;
        }
        
        void resizeElements( int newSize ) ///sets new mesh element size. Resizes elements, element2faces, level, rightChild.
        { 
            elements.resize( newSize );
            element2faces.resize( newSize );
            level.resize( newSize );
            rightChild.resize( newSize );
        }
        
        void resizeFaces( int newSize ) ///sets new number of faces. Resizes faces, faces2edges, levelFaces.
        { 
            faces.resize( newSize );
            face2edges.resize( newSize );
            levelFaces.resize( newSize );
            numFaces = newSize;
        }
        
        void resizeEdges( int newSize ) ///sets new number of edges. Resizes edges.
        {
            edges.resize( newSize );
            numEdges = newSize;
        }
        
        void resizeNodes( int newSize ) ///sets new number of nodes. Resizes coordinates array.
        {
            coordinates.resize( newSize ); 
            numNodes = newSize;
        }

        template< class NodeVector >
        void newNodeIndices( NodeVector & nodes2newNodes ) ///coordinates are rearanged via the given new nodes indices.
        {
            coordinates = rearange_inverse_delete( coordinates, nodes2newNodes );
        }
        template< class EdgeVector >
        void newEdgeIndices( EdgeVector & edge2newEdges ) {} 
        template< class FaceVector >
        void newFaceIndices( FaceVector & face2newFaces ) {} 
                
        void print() {
            coordinates.print("coordinates");
            elements.print("elements");
            faces.print("faces");
            edges.print("edges");
            element2faces.print("element2faces");
            face2edges.print("face2edges");
        }
        
        
    };
    
    class Grid_2D /// mesh structure of a 2-dimensional grid of triangles. See the documentation of Grid_3D.
    {        
    public:
        static const int dimensionworld = 2;
        static const int verticesPerElement = 3;
        static const int NV = verticesPerElement;
        static const int dim = dimensionworld;
        static const int DIM = dimensionworld;

        ManagedArray<int,NV> elements;
        ManagedArray<int,NV> element2edges;
        
        int& element2faces( std::size_t col, std::size_t row ) {return element2edges(col,row);}
        const int& element2faces( std::size_t col, std::size_t row ) const {return element2edges(col,row);}
        
        int face2edges( std::size_t col, std::size_t row ) {return 0;}
        const int face2edges( std::size_t col, std::size_t row ) const {return 0;}

        ManagedArray<unsigned char,1> level;
        ManagedArray<unsigned char,NV> rightChild;
 
        ManagedArray<int,NV-1> edges;
        
        std::size_t numEdges, numNodes;
        
        ManagedArray<double,DIM> coordinates;
        
        size_t nNodes() const { return coordinates.size(); }
        
        size_t nEdges() const { return edges.size(); }
        
        size_t nFaces() const { return edges.size(); }   
        
        size_t nElements() const { return elements.size(); }
        
        void initialTriangulation() {
            level.fill(0);
            std::fill(rightChild.begin(),rightChild.end(),false); 
            getElement2Faces ( elements , edges, element2edges );
            resizeEdges( edges.size() );    
            numNodes = coordinates.size();
        }
        
        void resizeElements( int newSize ){ 
            elements.resize( newSize );
            element2edges.resize( newSize );
            level.resize( newSize );
            rightChild.resize( newSize );
        }
        
        void resizeEdges( int newSize ){ 
            edges.resize( newSize );
            numEdges = newSize;
        }
        
        void resizeNodes( int newSize ){ 
            coordinates.resize( newSize ); 
            numNodes = newSize;
        }

        template< class NodeVector >
        void newNodeIndices( NodeVector & nodes2newNodes ) {
            coordinates = rearange_inverse_delete( coordinates, nodes2newNodes );
        }
        template< class EdgeVector >
        void newEdgeIndices( EdgeVector & edge2newEdges ) {} 
                
        void print() {
            coordinates.print("coordinates");
            elements.print("elements");
            edges.print("edges");
            element2edges.print("element2edges");
        }
        
    };
    
    template< class Grid, class Decomposition >
    void test_grid( Grid & grid, Decomposition &decomp ){
        test_grid(grid);
        assert( decomp.face2ranks.size() == grid.nFaces() );
        assert( decomp.edge2ranks.size() == grid.nEdges() );
        assert( decomp.node2ranks.size() == grid.nNodes() );
    }
    
    template< class Grid >
    void test_grid( Grid & grid ){
        assert( (int)grid.nFaces() == (int)grid.element2faces.max() + 1);
        assert( (int)grid.nEdges() == (int)grid.face2edges.max() + 1);
        assert( (int)grid.nElements() == (int)grid.level.size() );
        
        assert( (int)grid.elements.max() == (int)grid.edges.max() );
        assert( (int)grid.faces.max() == (int)grid.edges.max() );
        
        assert( (int)grid.element2faces.min() == 0 );
        assert( (int)grid.face2edges.min() == 0 );
        assert( (int)grid.elements.min() == 0 );
        assert( (int)grid.faces.min() == 0 );
        assert( (int)grid.edges.min() == 0 );
        
        ManagedArray<int> FaceCounter( grid.nFaces(), 0 );
        for( size_t el = 0; el < grid.element2faces.size(); ++el)
            for( size_t k = 0; k < 4; ++k )
                FaceCounter[  grid.element2faces(el,k) ]++;
        for( size_t i = 0; i < FaceCounter.size(); ++i){
            assert(FaceCounter[i] > 0 );
            assert(FaceCounter[i] < 3 );
        }
        
        ManagedArray<int> EdgeCounter( grid.nEdges(), 0 );
        for( size_t el = 0; el < grid.face2edges.size(); ++el)
            for( size_t k = 0; k < 3; ++k )
                EdgeCounter[ grid.face2edges(el,k) ]++;
        for( size_t i = 0; i < EdgeCounter.size(); ++i)
            assert(EdgeCounter[i] > 0 );
        
        for( size_t el = 0; el < grid.nElements(); ++el ){
            std::vector<int> tmp1( 12 );
            std::vector<int> tmp2( 4 );
            size_t ctr = 0;
            for( size_t k = 0; k < 4; ++k )
                for( size_t r = 0; r < 3; ++r)
                    tmp1[ctr++] = grid.faces(grid.element2faces(el,k),r);
            for( size_t k = 0; k < 4; ++k )
                tmp2[k] = grid.elements(el,k);
            vector_unique(tmp1);
            vector_unique(tmp2);
            assert( tmp1.size() == 4 );
            assert( tmp2.size() == 4 );
            assert( tmp1 == tmp2 );
        }
        
        for( size_t el = 0; el < grid.nFaces(); ++el ){
            std::vector<int> tmp1( 6 );
            std::vector<int> tmp2( 3 );
            size_t ctr = 0;
            for( size_t k = 0; k < 3; ++k )
                for( size_t r = 0; r < 2; ++r)
                    tmp1[ctr++] = grid.edges(grid.face2edges(el,k),r);
            for( size_t k = 0; k < 3; ++k )
                tmp2[k] = grid.faces(el,k);
            vector_unique(tmp1);
            vector_unique(tmp2);
            
            assert( tmp1.size() == 3 );
            assert( tmp2.size() == 3 );
            assert( tmp1 == tmp2 );
        }
            
        
    }
} // namespace conformingsimplexgrid
#endif // MESH_HH
