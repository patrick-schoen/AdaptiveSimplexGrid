#ifndef CS_DECOMPOSITION_2D_HH
#define CS_DECOMPOSITION_2D_HH

#include "mpi.h"

#include "../simplexGrid.hh"
#include "../parallel.hh"

namespace conformingsimplexgrid {

template< class Grid, class Communicator>
    class Decomposition_2D {
        typedef  Decomposition_2D<Grid,Communicator> ThisType; 
        
    public:
        static const int DIM = Grid::dimensionworld;
        static const int NV = Grid::verticesPerElement;
        
        ManagedArray< ManagedArray<char> > node2ranks;
        ManagedArray< ManagedArray<char> > edge2ranks;
        
        ManagedArray< ManagedArray<int> > rank2bdyNodes; 
        ManagedArray< ManagedArray<int> > rank2bdyEdges;
        
        Decomposition_2D( Grid &grid, Communicator &comm )
        : grid_(grid), comm_(comm)
        {
            rank2bdyNodes.resize( comm_.size() );
            rank2bdyEdges.resize( comm_.size() );
            communication_is_valid = false;
        }
                
        bool master() const {
            return comm_.master();
        }
        
        int rank() const {
            return (int)comm_.rank();
        }
        
        int size() const {
            return (int)comm_.size();
        }
        
        bool isBoundaryNode( int k ) const {
            return node2ranks[k].size() > 1;
        }
        
        bool isBoundaryEdge( int k ) const {
            return edge2ranks[k].size() > 1;
        }
        
        
        template< class T>
        int all2all( T *sendbuf, T *recvbuf, int size ) const {
            return comm_.all2all( sendbuf, recvbuf, size );
        }
        
        template< class T >
        void sumNodeVector( ManagedArray<T> &vec ) { 
            assert( vec.size() == node2ranks.size() );
            assert( node2ranks.size() == grid_.nNodes() );
            sumParallelDataVector( vec, rank2bdyNodes ); 
        }
                
        template< class T >
        void sumEdgeVector( ManagedArray<T> &vec ) { 
            sumParallelDataVector( vec, rank2bdyEdges ); 
        }
        
        template< class T >
        void sumParallelDataVector( ManagedArray<T> &vec, ManagedArray< ManagedArray<int> > &rank2bdydata ) 
        {
            init_communication();
            
            assert( communication_is_valid == true);
            assert(rank2bdydata[comm_.rank()].size() == 0);
            
            ManagedArray<int> processes(comm_.size(),0);
            ManagedArray<int> processes_k(comm_.size(),0);
            size_t counter = 0;
            for(int p = 0; p < (int)comm_.size(); p++)
                if( rank2bdydata[p].size() > 0 ){
                    processes[counter] = p;
                    processes_k[p] = counter;
                    counter++;
                }
            processes.resize(counter);
            

            ManagedArray< typename Communicator::Request > requests(processes.size());
            ManagedArray< ManagedArray<T> > buffer(processes.size());
            ManagedArray< typename Communicator::Request > recv_requests(processes.size());
            ManagedArray< ManagedArray<T> > recv_buffer(processes.size());
            
            for(int k = 0; k < (int)processes.size(); k++){
                const int p = processes[k];
                buffer[k].resize(rank2bdydata[p].size());
                recv_buffer[k].resize(rank2bdydata[p].size());   
            }
            
            for(int k = 0; k < (int)processes.size(); k++){
                const int p = processes[k];
                comm_.Irecv(recv_buffer[k].data(),recv_buffer[k].size(),p,1,&recv_requests[k]);
            }
            
            for(int k = 0; k < (int)processes.size(); k++){
                const int p = processes[k];
                for( int j = 0; j < (int)rank2bdydata[p].size(); j++){
                    vec[rank2bdydata[p][j]];
                    buffer[k][j] = vec[rank2bdydata[p][j]];
                }
            }
            
            for(int k = 0; k < (int)processes.size(); k++){
                const int p = processes[k];
                comm_.Isend(buffer[k].data(),buffer[k].size(),p,1,&requests[k]);
            }
            
            std::vector<bool> completed_messages( processes.size(), false );
            int nRecv = 0;
            while( nRecv < (int)processes.size() ) {               
                for(int k = 0; k < (int)processes.size(); k++){
                    const int p = processes[k];
                    if( !completed_messages[k] ){
                        int flag = 0;
                        comm_.test( &recv_requests[k], &flag);
                        if( flag ){
                            for(int j = 0; j < (int)rank2bdydata[p].size(); j++)
                                vec[rank2bdydata[p][j]] += recv_buffer[k][j];
                            nRecv++;
                            completed_messages[k] = true;
                        }
                    }
                }
            }   

            for(int k = 0; k < (int)processes.size(); k++){
                comm_.wait(&requests[k]);
                comm_.wait(&recv_requests[k]);
            }

        }

        template< class T >
        void scalar_min( T & x ) const
        {
            T temp;
            comm_.reduce_min( &x, &temp, 1);
            x = temp;
        }
        
        template< class T >
        void scalar_max( T & x ) const
        {
            T temp;
            comm_.reduce_max( &x, &temp, 1);
            x = temp;
        }
        
 
        template< class T >
        void scalar_sum( T & x ) const
        {
            T temp;
            comm_.reduce_sum( &x, &temp, 1);
            x = temp;
        }
        
        template< class T >
        void scalar_reduce_sum( T *send, T *recv ) const
        {
            comm_.reduce_sum( send, recv, 1);
        }
        
        double scalar_product( ManagedArray< double > &v, ManagedArray< double > &w ) const
        {
            assert( v.size() == w.size() );
            double value = 0;
            double value_local = 0.; 
            for( int i = 0; i < (int)v.size(); i++){
                if( node2ranks[i].size() ){
                    value_local += (1./(double)node2ranks[i].size()) * v[i] * w[i];
                }
                else
                    value_local += v[i] * w[i];
            }
            scalar_reduce_sum( &value_local, &value);
            return value;
        }

        template< class Partition>
        void init_partition( Partition &partition )
        { 
            init_partition( partition, node2ranks, edge2ranks );  
        }
        
        template< class Partition, class NodeArray, class EdgeArray >
        void init_partition( Partition &partition, NodeArray &node2ranks, EdgeArray &edge2ranks )
        {            
            node2ranks.resize( grid_.nNodes() );
            edge2ranks.resize( grid_.nEdges() );
            
            for( int k = 0; k < (int)node2ranks.size(); ++k )
                node2ranks[k].clear();
            for( int k = 0; k < (int)edge2ranks.size(); ++k )
                edge2ranks[k].clear();

            for( int el = 0; el < (int)grid_.nElements(); ++el ){
                const int p = partition[el];
                for( int k = 0; k < NV; ++k ){
                    node2ranks[ grid_.elements(el,k) ].push_back( p );
                    edge2ranks[ grid_.element2edges(el,k) ].push_back( p );
                }
            }
            
            for( int k = 0; k < (int)node2ranks.size(); ++k )
                vector_unique(node2ranks[k].array_);
            
            for( int k = 0; k < (int)edge2ranks.size(); ++k )
                vector_unique(edge2ranks[k].array_);
                   
            communication_is_valid = false;
        }
                
        void init_communication()
        {
            for( int k = 0; k < (int)rank2bdyEdges.size(); ++k )
                rank2bdyEdges[k].clear();
            for( int k = 0; k < (int)rank2bdyNodes.size(); ++k )
                rank2bdyNodes[k].clear();
                                         
            for( int k = 0; k < (int)edge2ranks.size(); ++k )
                for( int j = 0; j < (int)edge2ranks[k].size(); ++j){
                    if( edge2ranks[k][j] == comm_.rank() )
                        continue;
                    rank2bdyEdges[edge2ranks[k][j]].push_back( k ); 
                }
                
            for( int k = 0; k < (int)node2ranks.size(); ++k )
                for( int j = 0; j < (int)node2ranks[k].size(); ++j){
                    if( node2ranks[k][j] == comm_.rank() )
                        continue;
                    rank2bdyNodes[node2ranks[k][j]].push_back( k ); 
                }
            communication_is_valid = true;
        }
           
        template< class NodeVector >
        void removeNodes( NodeVector & remainingNodes )
        {
           std::size_t counter = 0;
           for(std::size_t k = 0; k < remainingNodes.size(); k++){
               if( remainingNodes[k] < 0 )
                   continue;
              node2ranks[counter++] = node2ranks[k];
           }
           node2ranks.resize(counter);
        }
        
        template< class EdgeVector >
        void removeEdges( EdgeVector & remainingEdges )
        {
           std::size_t counter = 0;          
           for(std::size_t k = 0; k < remainingEdges.size(); k++){
               if( remainingEdges[k] < 0 )
                   continue;
               edge2ranks[counter++] = edge2ranks[k];
           }
           edge2ranks.resize(counter);
        }
                
        template< class NodeVector >
        void removeNodes( NodeVector & remainingNodes, ThisType & subDecomposition) const
        {
           subDecomposition.node2ranks.resize( node2ranks.size() ); 
           std::size_t counter = 0;
           for( std::size_t k = 0; k < remainingNodes.size(); k++){
               if( remainingNodes[k] < 0 )
                   continue;
               subDecomposition.node2ranks[counter++] = node2ranks[k];
           }
           subDecomposition.node2ranks.resize(counter);
        }
        
        template< class EdgeVector >
        void removeEdges( EdgeVector & remainingEdges, ThisType & subDecomposition) const
        {
           subDecomposition.edge2ranks.resize( remainingEdges.size() ); 
           std::size_t counter = 0;          
           for( std::size_t k = 0; k < remainingEdges.size(); k++){
               if( remainingEdges[k] < 0 )
                   continue;
               subDecomposition.edge2ranks[counter++] = edge2ranks[k];
           }
           subDecomposition.edge2ranks.resize(counter);
        }
 
    private:
        Grid &grid_;
        Communicator &comm_;
        bool communication_is_valid;
    };
      
    template < class Subgrid, class SubDecomposition, class Grid, class Decomposition, class ElementVector >
    void getSubGrid2D( Subgrid & subgrid, SubDecomposition & subDecomposition, int rank, 
                    const Grid &grid, const Decomposition &decomposition, const ElementVector &partition )
    {
        const int NV = Grid::verticesPerElement;
        std::vector<bool> remainingElements(grid.nElements(), false);

        for( size_t el = 0; el < grid.nElements(); ++el)
            if( partition[el] == rank )
                remainingElements[el] = true;
            
        remove_false( grid.elements, remainingElements, subgrid.elements);        
        remove_false( grid.element2edges, remainingElements, subgrid.element2edges);   
        
        ManagedArray<int> nodes2newNodes( grid.nNodes(), -1 );
        for( size_t el = 0; el < subgrid.nElements(); ++el)
            for( size_t j = 0; j < NV; ++j)
                nodes2newNodes[ subgrid.elements(el,j) ] = 1;

        indexVector_enumerate( nodes2newNodes );
        index2newIndex( subgrid.elements, nodes2newNodes );
        
        remove_negative( grid.coordinates, nodes2newNodes, subgrid.coordinates );

        ManagedArray<int> edges2newEdges( grid.nEdges(), -1 );
        for( size_t el = 0; el < grid.nElements(); ++el)
            if( partition[el] == rank )
                for( size_t k = 0; k < NV; ++k )
                    edges2newEdges[ grid.element2edges(el,k) ] = 1;

        indexVector_enumerate( edges2newEdges );        
        
            
        index2newIndex( subgrid.element2edges, edges2newEdges );
        
        remove_negative( decomposition.edge2ranks, edges2newEdges, subDecomposition.edge2ranks );
        remove_negative( decomposition.node2ranks, nodes2newNodes, subDecomposition.node2ranks );
    }
    
}
 
#endif // 
    