/**
@file element2subentity.hh
@brief Functions to generated mesh entity structures, like element2faces and faces.
*/
#ifndef _CONFORMINGSIMPLEXGRID_ELEMENT2SUBENTITY_HH
#define _CONFORMINGSIMPLEXGRID_ELEMENT2SUBENTITY_HH

//#include "reference_element.hh"

namespace conformingsimplexgrid {

template< std::size_t NV >
class ReferenceSimplex 
{
public:
    ReferenceSimplex() 
    {
        for( std::size_t i = 0; i < NV; ++i )
            for( std::size_t j = 0; j < NV-1; ++j )
                referenceSimplex[ i*(NV-1) + j ] = j + ( (j + i) >= NV-1 );
    }
    
    constexpr std::size_t size() const { return NV*(NV-1); }
    
    constexpr int operator[] (const int k) const { return referenceSimplex[k]; }

private:
    std::array< int, NV*(NV-1) > referenceSimplex;
};
    
template< class Grid >
void getEntity2Subentity( Grid & grid )
{
    getElement2Faces ( grid.elements , grid.faces, grid.element2faces );
    grid.resizeFaces( grid.faces.size() );
    getElement2Faces ( grid.faces , grid.edges, grid.face2edges );
    grid.resizeEdges( grid.edges.size() );
}

template< class Elements, class Faces, class Element2Faces  >
void getFacesForElements( Elements & elements, Faces & faces, Element2Faces & element2faces )
{
    static const std::size_t NV = Elements::DIM;
    faces.resize(element2faces.max() + 1 );        
    assert( (int)faces.size() == (int)element2faces.max() + 1 );
    assert( elements.size() == element2faces.size() );
    ReferenceSimplex<NV> referenceFaces; 
    for( std::size_t el = 0; el < element2faces.size(); ++el )
        for( std::size_t k = 0, J = 0; k < NV; ++k )
            for( std::size_t j = 0; j < NV - 1; ++j )
                faces(element2faces(el,k),j) = elements(el,referenceFaces[J++]);
}

template< class Elements, class Faces, class Element2Faces, class Level  >
void getFacesForElements( Elements & elements, Faces & faces, Element2Faces & element2faces, Level & level )
{
    static const std::size_t NV = Elements::DIM;
    faces.resize(element2faces.max() + 1 );        
    assert( (int)faces.size() == (int)element2faces.max() + 1 );
    assert( elements.size() == element2faces.size() );
    ReferenceSimplex<NV> referenceFaces;
    const int refGamma[9] = {1,2,3,2,1,3,3,1,2};
    for( std::size_t el = 0; el < element2faces.size(); ++el ){
        const int type_idx = (level[el] % (NV-1)-1)*3 + 3;
        for( std::size_t k = 0, J = 0; k < NV-1; ++k )
            for( std::size_t j = 0; j < NV - 1; ++j )
                faces(element2faces(el,k),j) = elements(el,referenceFaces[J++]);
        for( std::size_t j = 0, J = type_idx; j < NV - 1; ++j )
                faces(element2faces(el,NV-1),j) = elements(el,refGamma[J++]);    
    }
}
    
template< class Elements, class Faces, class Element2Faces >
void getElement2Faces ( Elements & elements, Faces & faces, Element2Faces & element2faces )
/// generated element2faces and faces array from mesh elements. Here elements can also be faces. In this case the function returns the edges and faces2edges structure.
{
    static const std::size_t NV = Elements::DIM;
    ReferenceSimplex<NV> referenceFaces;
    
    faces.resize( elements.size() * NV );

    for( std::size_t el = 0, I = 0; el < elements.size(); ++el )
        for( std::size_t J = 0; J < NV*(NV-1); I++ ) 
            for( std::size_t j = 0; j < NV - 1; ++j )
                faces(I,j) = elements(el, referenceFaces[J++] );
    
    ManagedArray<int> ia;
    ManagedArray<int,NV-1> faces_sorted = sortEachRow( faces );
    unique_rows_void( faces_sorted, ia, element2faces.array_ );
    faces = rearange( faces, ia);
}
     
template< class Grid, class Neighbors >
void getNeighbors ( Grid &grid, Neighbors &neighbors ) 
/// get neighboring elements.
{
    static const std::size_t NV = Grid::NV;
    ManagedArray<int,2> neighbors_on_face( grid.nFaces() );
    neighbors_on_face.fill(-1);
    neighbors.resize( grid.nElements() );
    for(size_t el = 0; el < grid.nElements(); ++el)
        for(size_t k = 0; k < NV; ++k){
            if(neighbors_on_face(grid.element2faces(el,k),0) < 0 ){
                neighbors_on_face(grid.element2faces(el,k),0) = el;
                neighbors_on_face(grid.element2faces(el,k),1) = k;
                neighbors(el,k) = -1;
            } else {
                const int N = neighbors_on_face(grid.element2faces(el,k),0);
                const int kN = neighbors_on_face(grid.element2faces(el,k),1);
                neighbors(el,k) = N;
                neighbors(N,kN) = el;
            }
        }
}
    
template<class Communicator, class Decomposition, class Vector, class Grid, class Neighbors >
void getGlobalNeighbors( Communicator &comm, Decomposition &decomposition, Vector &vtxdist, Grid &grid, Neighbors &neighbors )
/// get neighboring elements in a parallel distributed mesh.
{
    const int p = (int)comm.rank();
    static const std::size_t NV = Grid::NV;
    
    getNeighbors( grid, neighbors );

    for(size_t el = 0; el < grid.nElements(); ++el)
        for(size_t k = 0; k < NV; ++k)
            if( neighbors(el,k) >= 0 )
                neighbors(el,k) += vtxdist[p];
        
    ManagedArray<int,1> neighbors_on_bdyFace_1( grid.nFaces(), -1 );
    ManagedArray<int,1> neighbors_on_bdyFace_2( grid.nFaces(), 0 );
    
    for(size_t k = 0; k < decomposition.face2ranks.size(); ++k)
        if( decomposition.isBoundaryFace(k) ) {
            neighbors_on_bdyFace_1[k] = 1;
        }
            
    for(size_t el = 0; el < grid.nElements(); ++el)
        for(size_t k = 0; k < NV; ++k)
            if( neighbors_on_bdyFace_1[grid.element2faces(el,k)] >= 0 ){
                neighbors_on_bdyFace_1[grid.element2faces(el,k)] = el + vtxdist[p];
                neighbors_on_bdyFace_2[grid.element2faces(el,k)] = el + vtxdist[p];
            }
            
    decomposition.sumFaceVector( neighbors_on_bdyFace_1 );
    neighbors_on_bdyFace_1 -= neighbors_on_bdyFace_2;
    
    for(size_t el = 0; el < grid.nElements(); ++el)
        for(size_t k = 0; k < NV; ++k)
            if( neighbors_on_bdyFace_1[grid.element2faces(el,k)] >= 0 ){
                assert( neighbors(el,k) < 0 );
                neighbors(el,k) = neighbors_on_bdyFace_1[grid.element2faces(el,k)];
            }
}
    
} // namespace conformingsimplexgrid

#endif // #ifndef _CONFORMINGSIMPLEXGRID_ELEMENT2SUBENTITY_HH
