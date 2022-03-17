#ifndef AC_DATA_CONTAINER_HH
#define AC_DATA_CONTAINER_HH

#include "../../../grid/mesh.hh"
#include "../../../grid/parallel.hh"

using namespace conformingsimplexgrid;

template< class Communicator, class Grid, class Curve >
class AC_Manager 
{
public:
    
    explicit AC_Manager( Communicator &comm_in, Grid &grid_in, ManagedArray<double,1> &u_in, ManagedArray<double,1> &u_old_in, Curve &curve_in)
    : comm(comm_in), grid(grid_in), u(u_in), u_old(u_old_in), curve(curve_in)
    {}
    typedef AC_Manager<Communicator, Grid, Curve> ThisType;
    
    static const int NV = Grid::NV;
    static const int DIM = Grid::DIM;
    
    Communicator &comm;
    Grid &grid;
    ManagedArray<double> &u;
    ManagedArray<double> &u_old;
    Curve &curve;
        
    template< class Refine_Data >
    void pre_refine( Refine_Data &data )
    {
        int nNewNodes = 0;
        for( size_t k = 0; k < data.edge2newNode.size(); k++ )
            if( data.edge2newNode[k] >= 0 )
                nNewNodes++;
        u.resize( u.size() + nNewNodes );
        u_old.resize( u_old.size() + nNewNodes );
            
        std::vector<bool> nodeFlag( u.size() );
        std::fill( nodeFlag.begin(), nodeFlag.end(), false );
        for(size_t el = 0; el < data.element2newNode.size(); el++)
            if( data.element2newNode[el] >= 0 ){
                const int y = data.element2newNode[el];
                if( nodeFlag[ y ] == false ){
                    const int e1 = grid.elements(el,0);
                    const int e2 = grid.elements(el,NV-1);
                    u[y] = .5*(u[e1] + u[e2]);
                    u_old[y] = .5*(u_old[e1] + u_old[e2]);
                    nodeFlag[y] = true;
                }
            }
        curve.pre_refine(data);
    }
    
    template< class Refine_Data >
    void refine( Refine_Data & )
    {
    }
    
    template< class Coarse_Data >
    void coarse( Coarse_Data &data ) {
        u = rearange_inverse_delete( u, data.node2newNodes );
        u_old = rearange_inverse_delete( u_old, data.node2newNodes );
        curve.coarse(data);
    }
    
    template< class Vector >
    void newIndices( Vector &, Vector &, Vector &node2newNodes  ) {
        u = rearange( u, node2newNodes );
        u_old = rearange( u_old, node2newNodes );        
    }
        
    template<class BoolVector, class Vector> 
    void preparePacking(BoolVector &remainingElements, Vector &faces2newFaces, Vector &edges2newEdges, Vector &nodes2subNodes)
    {
        u_buffer.resize(nodes2subNodes.size());
        int counter = 0;
        for( size_t k = 0; k < nodes2subNodes.size(); k++ )
            if( nodes2subNodes[k] >= 0 ) {
                u_buffer[nodes2subNodes[k]] = u[k];
                counter++;
            }
        u_buffer.resize(counter);
        curve.preparePacking(remainingElements, faces2newFaces, edges2newEdges, nodes2subNodes);
    }
    
    size_t size() {return u.size() + curve.size();}
    
    void pack( void *buffer, std::size_t bufferSize, int &position ) {
        int message_size = (int)u_buffer.size();
        comm.pack( &message_size, 1, buffer, bufferSize, position );
        comm.pack( u_buffer.data(), u_buffer.size(), buffer, bufferSize, position );
        curve.pack(buffer,bufferSize,position);
    }
    
    void unpack( void *buffer, std::size_t bufferSize, int &position ) {
        int message_size;
        comm.unpack( &message_size, 1, buffer, bufferSize, position );
        u_buffer.resize( message_size );
        comm.unpack( u_buffer.data(), u_buffer.size(), buffer, bufferSize, position );
        curve.unpack(buffer,bufferSize,position);
    }
        
    std::size_t bufferSize() {
        return comm.template bufferSize<int>(1)
             + comm.template bufferSize<double>(u_buffer.size())
             + curve.bufferSize();
    }

    template< class Vector > 
    void merge( Vector &, Vector &, Vector &, Vector &, Vector &, Vector &uniqueNodes ){
        u.conjoin(u_buffer);
        u = array_copyIndices( u, uniqueNodes );
        u_old = u;
        
        
    }

    void conjoin(){
        u.conjoin(u_buffer);
        u_old = u;
        curve.conjoin();
    }
    
private:
    ManagedArray<double> u_buffer;
};

template< class Communicator, class Grid >
class AC_Manager2
{
public:
    
    explicit AC_Manager2( Communicator &comm_in, Grid &grid_in, ManagedArray<double,1> &u_in, ManagedArray<double,1> &u_old_in)
    : comm(comm_in), grid(grid_in), u(u_in), u_old(u_old_in)
    {}
    typedef AC_Manager2<Communicator, Grid> ThisType;
    
    static const int NV = Grid::NV;
    static const int DIM = Grid::DIM;
    
    Communicator &comm;
    Grid &grid;
    ManagedArray<double> &u;
    ManagedArray<double> &u_old;
        
    template< class Refine_Data >
    void pre_refine( Refine_Data &data )
    {
        int nNewNodes = 0;
        for( size_t k = 0; k < data.edge2newNode.size(); k++ )
            if( data.edge2newNode[k] >= 0 )
                nNewNodes++;
        u.resize( u.size() + nNewNodes );
        u_old.resize( u_old.size() + nNewNodes );
            
        std::vector<bool> nodeFlag( u.size() );
        std::fill( nodeFlag.begin(), nodeFlag.end(), false );
        for(size_t el = 0; el < data.element2newNode.size(); el++)
            if( data.element2newNode[el] >= 0 ){
                const int y = data.element2newNode[el];
                if( nodeFlag[ y ] == false ){
                    const int e1 = grid.elements(el,0);
                    const int e2 = grid.elements(el,NV-1);
                    u[y] = .5*(u[e1] + u[e2]);
                    u_old[y] = .5*(u_old[e1] + u_old[e2]);
                    nodeFlag[y] = true;
                }
            }
    }
    
    template< class Refine_Data >
    void refine( Refine_Data & )
    {
    }
    
    template< class Coarse_Data >
    void coarse( Coarse_Data &data ) {
        u = rearange_inverse_delete( u, data.node2newNodes );
        u_old = rearange_inverse_delete( u_old, data.node2newNodes );
    }
    
    template< class Vector >
    void newIndices( Vector &, Vector &, Vector &node2newNodes  ) {
        u = rearange( u, node2newNodes );
        u_old = rearange( u_old, node2newNodes );        
    }
        
    template<class BoolVector, class Vector> 
    void preparePacking(BoolVector &remainingElements, Vector &faces2newFaces, Vector &edges2newEdges, Vector &nodes2subNodes)
    {
        u_buffer.resize(nodes2subNodes.size());
        int counter = 0;
        for( size_t k = 0; k < nodes2subNodes.size(); k++ )
            if( nodes2subNodes[k] >= 0 ) {
                u_buffer[nodes2subNodes[k]] = u[k];
                counter++;
            }
        u_buffer.resize(counter);
    }
    
    size_t size() {return u.size();}
    
    void pack( void *buffer, std::size_t bufferSize, int &position ) {
        int message_size = (int)u_buffer.size();
        comm.pack( &message_size, 1, buffer, bufferSize, position );
        comm.pack( u_buffer.data(), u_buffer.size(), buffer, bufferSize, position );
    }
    
    void unpack( void *buffer, std::size_t bufferSize, int &position ) {
        int message_size;
        comm.unpack( &message_size, 1, buffer, bufferSize, position );
        u_buffer.resize( message_size );
        comm.unpack( u_buffer.data(), u_buffer.size(), buffer, bufferSize, position );
    }
        
    std::size_t bufferSize() {
        return comm.template bufferSize<int>(1)
             + comm.template bufferSize<double>(u_buffer.size());
    }

    template< class Vector > 
    void merge( Vector &, Vector &, Vector &, Vector &, Vector &, Vector &uniqueNodes ){
        u.conjoin(u_buffer);
        u = array_copyIndices( u, uniqueNodes );
        u_old = u;  
    }

    void conjoin(){
        u.conjoin(u_buffer);
        u_old = u;
    }
    
private:
    ManagedArray<double> u_buffer;
};

#endif // AC_DATA_CONTAINER_HH
