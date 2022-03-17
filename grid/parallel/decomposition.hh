#ifndef CS_DECOMPOSITION_HH
#define CS_DECOMPOSITION_HH

#include "mpi.h"
#include<bitset>

#include "../simplexGrid.hh"
#include "../parallel.hh"

namespace conformingsimplexgrid {
    
    template< class Grid, class Communicator>
    class Decomposition_3D {
        typedef  Decomposition_3D<Grid,Communicator> ThisType; 
        
    public:
        static const int DIM = Grid::dimensionworld;
        static const int NV = Grid::verticesPerElement;
        
        ManagedArray< ManagedArray<char> > node2ranks;
        ManagedArray< ManagedArray<char> > edge2ranks;
        ManagedArray< ManagedArray<char> > face2ranks;
        
        ManagedArray< ManagedArray<int> > rank2bdyNodes; 
        ManagedArray< ManagedArray<int> > rank2bdyEdges;
        ManagedArray< ManagedArray<int> > rank2bdyFaces;
        
        Decomposition_3D( Grid &grid, Communicator &comm )
        : grid_(grid), comm_(comm)
        {
            rank2bdyNodes.resize( comm_.size() );
            rank2bdyEdges.resize( comm_.size() );
            rank2bdyFaces.resize( comm_.size() );
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
        
        bool isBoundaryFace( int k ) const {
            return face2ranks[k].size() > 1;
        }
        
        template< class T>
        int all2all( T *sendbuf, T *recvbuf, int size ) const {
            return comm_.all2all( sendbuf, recvbuf, size );
        }
        
        template< class T >
        void sumNodeVector( ManagedArray<T> &vec ) const { 
            assert( vec.size() == node2ranks.size() );
            assert( node2ranks.size() == grid_.nNodes() );
            sumParallelDataVector( vec, rank2bdyNodes ); 
        }
        
        template< class T, class Request >
        void sumNodeVector_send( ManagedArray<T> &vec, ManagedArray< Request > &requests, ManagedArray< ManagedArray<T> > & buffer ) { 
            assert( vec.size() == node2ranks.size() );
            assert( node2ranks.size() == grid_.nNodes() );
            sumParallelDataVector_send( vec, rank2bdyNodes, requests, buffer ); 
        }
        
        template< class T, class Request >
        void sumNodeVector_recv( ManagedArray<T> &vec, ManagedArray< Request > &requests, ManagedArray< ManagedArray<T> > & buffer ) { 
            assert( vec.size() == node2ranks.size() );
            assert( node2ranks.size() == grid_.nNodes() );
            sumParallelDataVector_recv( vec, rank2bdyNodes, requests, buffer); 
        }
        
        template< class T, class Request >
        void sumNodeVector_wait( ManagedArray<T> &vec, ManagedArray< Request > &requests, ManagedArray< ManagedArray<T> > & buffer ) { 
            assert( vec.size() == node2ranks.size() );
            assert( node2ranks.size() == grid_.nNodes() );  
            sumParallelDataVector_wait( vec, rank2bdyNodes, requests, buffer); 
        }
        
        template< class T >
        void sumEdgeVector( ManagedArray<T> &vec ) const { 
            sumParallelDataVector( vec, rank2bdyEdges ); 
        }
        
        template< class T >
        void sumFaceVector( ManagedArray<T> &vec ) const { 
            sumParallelDataVector( vec, rank2bdyFaces );
        }
       
        template< class T >
        void sumParallelDataVector( ManagedArray<T> &vec, const ManagedArray< ManagedArray<int> > &rank2bdydata ) const
        {            
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
            while( nRecv < (int)processes.size() ){               
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
        
        template< class T, class Request>
        void sumParallelDataVector_send( ManagedArray<T> &vec, ManagedArray< ManagedArray<int> > &rank2bdydata, 
                                         ManagedArray< Request > &requests, ManagedArray< ManagedArray<T> > & buffer ) 
        {
            if(communication_is_valid == false)
                init_communication();
            assert(rank2bdydata[comm_.rank()].size() == 0);

            for(int p = 0; p < (int)comm_.size(); p++)
                buffer[p].resize(rank2bdydata[p].size());
            
            for(int p = 0; p < (int)comm_.size(); p++){
                if( rank2bdyNodes[p].size() < 1 ){
                    continue;
                }
                for( int k = 0; k < (int)rank2bdydata[p].size(); k++)
                    buffer[p][k] = vec[rank2bdydata[p][k]];
                comm_.Isend(buffer[p].data(),buffer[p].size(),p,1,&requests[p]);
            }
        }
        
        template< class T, class Request >
        void sumParallelDataVector_recv( ManagedArray<T> &, ManagedArray< ManagedArray<int> > &rank2bdydata, 
                                         ManagedArray< Request > &requests, ManagedArray< ManagedArray<T> > & recv_buffer )
        {
            if(communication_is_valid == false)
                init_communication();
            assert(rank2bdydata[comm_.rank()].size() == 0);
                        
            for(int p = 0; p < (int)comm_.size(); p++)
                recv_buffer[p].resize(rank2bdydata[p].size());   

            for(int p = 0; p < (int)comm_.size(); p++){
                if( rank2bdyNodes[p].size() < 1 )
                    continue;
                comm_.Irecv(recv_buffer[p].data(),recv_buffer[p].size(),p,1,&requests[p]);
            }
        }
            
        template< class T, class Request >
        void sumParallelDataVector_wait( ManagedArray<T> &vec, ManagedArray< ManagedArray<int> > &rank2bdydata, 
                                         ManagedArray< Request > &requests, ManagedArray< ManagedArray<T> > & recv_buffer )
        {
            for(int p = 0; p < (int)comm_.size(); p++){
                if( rank2bdyNodes[p].size() < 1 )
                    continue;
                comm_.wait(&requests[p]);
            }
            
            for(int p = 0; p < (int)comm_.size(); p++){
                if( rank2bdyNodes[p].size() < 1 )
                    continue;
                for(int k = 0; k < (int)rank2bdydata[p].size(); k++)
                    vec[rank2bdydata[p][k]] += recv_buffer[p][k];
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
        
        template<class ProlongationData>
        void pre_refine( ProlongationData &data ){}
        
        template<class ProlongationData>
        void refine( ProlongationData &data ){
            
            node2ranks.resize( grid_.nNodes() );
            edge2ranks.resize( grid_.nEdges() );
            face2ranks.resize( grid_.nFaces() );
            
            for(size_t k = 0; k < data.edge2newNode.size(); ++k){
                if( data.edge2newNode[k] < 0 )
                    continue;
                node2ranks[ data.edge2newNode[k] ] = edge2ranks[k];
            }
                       
            for(size_t k = 0; k < data.left2rightEdges.size(); ++k){
                if( data.left2rightEdges[k] < 0 )
                    continue;
                edge2ranks[ data.left2rightEdges[k] ] = edge2ranks[k];
            }
           
            for(size_t k = 0; k < data.left2rightFaces.size(); ++k){
                if( data.left2rightFaces[k] < 0 )
                    continue;
                face2ranks[ data.left2rightFaces[k] ] = face2ranks[k];
            }
            
            for(size_t k = 0; k < data.face2innerEdge.size(); ++k){
                if( data.face2innerEdge[k] < 0 )
                    continue;
                edge2ranks[ data.face2innerEdge[k] ] = face2ranks[k];
            }
            
            
            for(size_t k = 0; k < data.element2innerFace.size(); ++k){
                if( data.element2innerFace[k] < 0 )
                    continue;
                face2ranks[ data.element2innerFace[k] ].resize(1);
                face2ranks[ data.element2innerFace[k] ][0] = (int)rank();
            }
                  
            init_communication();
          
        }
        
        template<class ProlongationData>
        void coarse( ProlongationData &data) {
            
            removeNodes( data.node2newNodes );
            removeEdges( data.edge2newEdges );
            removeFaces( data.face2newFaces );
            
            init_communication();
        }
                
        void pack( void *buffer, std::size_t bufferSize, int &position ) 
        {
            assert(node2ranks.size());
            assert(edge2ranks.size());
            assert(face2ranks.size());

            pack_entity( comm_, buffer, bufferSize, position, node2ranks );
            pack_entity( comm_, buffer, bufferSize, position, edge2ranks );
            pack_entity( comm_, buffer, bufferSize, position, face2ranks );
        }
        
        void unpack( void *buffer, std::size_t bufferSize, int &position ) 
        {
            unpack_entity( comm_, buffer, bufferSize, position, node2ranks );
            unpack_entity( comm_, buffer, bufferSize, position, edge2ranks );
            unpack_entity( comm_, buffer, bufferSize, position, face2ranks );
        }
        

        std::size_t bufferSize()
        {   
            std::size_t bufferSize = 0;
            bufferSize += bufferSize_entity( comm_, node2ranks );
            bufferSize += bufferSize_entity( comm_, edge2ranks );
            bufferSize += bufferSize_entity( comm_, face2ranks );
            return bufferSize;
        }
        
        template< class Partition>
        void init_partition( Partition &partition )
        { 
            init_partition( partition, node2ranks, edge2ranks, face2ranks );  
        }
        
        template< class Partition, class NodeArray, class EdgeArray, class FaceArray >
        void init_partition( Partition &partition, NodeArray &node2ranks, EdgeArray &edge2ranks, FaceArray &face2ranks )
        {            
            node2ranks.resize( grid_.nNodes() );
            edge2ranks.resize( grid_.nEdges() );
            face2ranks.resize( grid_.nFaces() );
            
            for( int k = 0; k < (int)node2ranks.size(); ++k ) node2ranks[k].clear();
            for( int k = 0; k < (int)edge2ranks.size(); ++k ) edge2ranks[k].clear();
            for( int k = 0; k < (int)face2ranks.size(); ++k ) face2ranks[k].clear();

            for( int el = 0; el < (int)grid_.nElements(); ++el ){
                const int p = partition[el];
                for( int k = 0; k < NV; ++k ){
                    node2ranks[ grid_.elements(el,k) ].push_back( p );
                    face2ranks[ grid_.element2faces(el,k) ].push_back( p );
                }
            }
            
            for( int k = 0; k < (int)node2ranks.size(); ++k ) vector_unique(node2ranks[k].array_);
            for( int k = 0; k < (int)face2ranks.size(); ++k ) vector_unique(face2ranks[k].array_);
            
            for( int el = 0; el < (int)grid_.nFaces(); ++el )
                for( size_t j = 0; j < face2ranks[el].size(); ++j ){
                    const int p = face2ranks[el][j];
                    for( int k = 0; k < NV - 1; ++k )
                        edge2ranks[grid_.face2edges(el,k)].push_back( p );
                }
     
            for( int k = 0; k < (int)edge2ranks.size(); ++k ) vector_unique(edge2ranks[k].array_);
       
            communication_is_valid = false;
        }
          
        void print ()
        {     
            
            printf("node2ranks \n");       
            for( size_t k = 0; k < node2ranks.size(); ++k){
                printf("%d :",(int)k);
                for( size_t j = 0; j < node2ranks[k].size(); ++j)
                    printf(" %d",(int)node2ranks[k][j]);
                printf("\n");
            }

            printf("edge2ranks \n");       
            for( size_t k = 0; k < edge2ranks.size(); ++k){
                printf("%d :",(int)k);
                for( size_t j = 0; j < edge2ranks[k].size(); ++j)
                    printf(" %d",(int)edge2ranks[k][j]);
                printf("\n");
            }

            printf("face2ranks \n");       
            for( size_t k = 0; k < face2ranks.size(); ++k){
                printf("%d :",(int)k);
                for( size_t j = 0; j < face2ranks[k].size(); ++j)
                    printf(" %d",(int)face2ranks[k][j]);
                printf("\n");
            }

        }
        
        void init_communication()
        {
            for( int k = 0; k < (int)rank2bdyFaces.size(); ++k )
                rank2bdyFaces[k].clear();
            for( int k = 0; k < (int)rank2bdyEdges.size(); ++k )
                rank2bdyEdges[k].clear();
            for( int k = 0; k < (int)rank2bdyNodes.size(); ++k )
                rank2bdyNodes[k].clear();
                        
            for( int k = 0; k < (int)face2ranks.size(); ++k )
                for( int j = 0; j < (int)face2ranks[k].size(); ++j){
                    if( face2ranks[k][j] == comm_.rank() )
                        continue;
                    rank2bdyFaces[face2ranks[k][j]].push_back( k );
                }
                 
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
                
        template< class FaceVector  >
        void removeFaces( FaceVector & remainingFaces )
        {      
           ManagedArray< ManagedArray<int> > face2ranks_temp( face2ranks.size() );
           for(std::size_t k = 0; k < face2ranks.size(); k++)
               face2ranks_temp[k] = face2ranks[k];         
           std::size_t counter = 0;  
           for(std::size_t k = 0; k < remainingFaces.size(); k++){
               if( remainingFaces[k] < 0 )
                   continue;
               face2ranks[counter++] = face2ranks_temp[k];
           } 
           face2ranks.resize(counter);
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
        
        template< class FaceVector  >
        void removeFaces( FaceVector & remainingFaces, ThisType & subDecomposition ) const
        {
           assert(remainingFaces.size());
           subDecomposition.face2ranks.resize( remainingFaces.size() );          
           std::size_t counter = 0;  
           for( std::size_t k = 0; k < remainingFaces.size(); k++){
               if( remainingFaces[k] < 0 )
                   continue;
               subDecomposition.face2ranks[counter++] = face2ranks[k];
           } 
           subDecomposition.face2ranks.resize(counter);
        } 
        
        bool testNodeCommunication()
        { 
            bool returnBool = true;
            
            init_communication();
            
            
            for( int k = 0; k < face2ranks.size(); k++){
                assert(face2ranks[k].size() > 0 );
                assert(face2ranks[k].size() < 3 );
            }
            
            for( int k = 0; k < edge2ranks.size(); k++){
                assert(edge2ranks[k].size() > 0 );
                assert(edge2ranks[k].size() <= size() );
            }
            
            for( int k = 0; k < node2ranks.size(); k++){
                assert(node2ranks[k].size() > 0 );
                assert(node2ranks[k].size() <= size() );
            }
            
            for( int k = 0; k < face2ranks.size(); k++)
                for( int r = 0; r < face2ranks[k].size() - 1; r++)
                    assert(face2ranks[k][r] < face2ranks[k][r+1] );
            
            for( int k = 0; k < edge2ranks.size(); k++)
                for( int r = 0; r < edge2ranks[k].size() - 1; r++)
                    assert(edge2ranks[k][r] < edge2ranks[k][r+1] ); 
                
            for( int k = 0; k < node2ranks.size(); k++)
                for( int r = 0; r < node2ranks[k].size() - 1; r++)
                    assert(node2ranks[k][r] < node2ranks[k][r+1] ); 

                

            ManagedArray<int> FacesOnRank1(size(), 0);
            ManagedArray<int> EdgesOnRank1(size(), 0);
            ManagedArray<int> NodesOnRank1(size(), 0);
            ManagedArray<int> FacesOnRank2(size(), 0);
            ManagedArray<int> EdgesOnRank2(size(), 0);
            ManagedArray<int> NodesOnRank2(size(), 0);

            for( int k = 0; k < node2ranks.size(); k++)
                for( int r = 0; r < node2ranks[k].size(); r++){
                    const int p = node2ranks[k][r];
                    NodesOnRank1[p]++;
                }
            for( int k = 0; k < edge2ranks.size(); k++)
                for( int r = 0; r < edge2ranks[k].size(); r++){
                    const int p = edge2ranks[k][r];
                    EdgesOnRank1[p]++;
                }
            for( int k = 0; k < face2ranks.size(); k++)
                for( int r = 0; r < face2ranks[k].size(); r++){
                    const int p = face2ranks[k][r];
                    FacesOnRank1[p]++;
                }

            all2all( NodesOnRank1.data(), NodesOnRank2.data(), 1 );
            all2all( EdgesOnRank1.data(), EdgesOnRank2.data(), 1 );
            all2all( FacesOnRank1.data(), FacesOnRank2.data(), 1 );
    
            for( int p = 0; p < size(); p++ ){
                if( p == rank() )
                    continue;
                if( !(NodesOnRank1[p] == NodesOnRank2[p])){
                    
                    printf("Nodes, other rank = %d, this rank = %d, %d -- %d \n", p, rank(), NodesOnRank1[p], NodesOnRank2[p]);
                    for( int k = 0; k < node2ranks.size(); k++)
                        for( int r = 0; r < node2ranks[k].size(); r++){
                            const int q = node2ranks[k][r];
                            if( q == p)
                                printf("node = %d, %lf %lf %lf \n", k, grid_.coordinates(k,0),grid_.coordinates(k,1), grid_.coordinates(k,2) );
                        }
                    returnBool = false;
                }
            }
            for( int p = 0; p < size(); p++ ){
                if( p == rank() )
                    continue;
                if( !(EdgesOnRank1[p] == EdgesOnRank2[p])){
                    printf("Edges, other rank = %d, this rank = %d, %d -- %d \n", p, rank(), EdgesOnRank1[p], EdgesOnRank2[p]);
                    returnBool = false;
                }
            }
            
            for( int p = 0; p < size(); p++ ){
                if( p == rank() )
                    continue;
                if( !(FacesOnRank1[p] == FacesOnRank2[p])){
                    
                    
                    printf("Faces, other rank = %d, this rank = %d, %d -- %d \n", p, rank(), FacesOnRank1[p], FacesOnRank2[p]);
                    returnBool = false;
                }
            }
            
            if( returnBool == false )
                printf("failed ahead of comm test, level max = %d\n", grid_.level.max());
            assert(returnBool);
      
            ManagedArray<int> NodeNumberRanks(node2ranks.size(), 0);
            ManagedArray<int> EdgeNumberRanks(edge2ranks.size(), 0);
            ManagedArray<int> FaceNumberRanks(face2ranks.size(), 0);
            
            for( int k = 0; k < node2ranks.size(); k++)
                NodeNumberRanks[k] = (int)node2ranks[k].size();
            ManagedArray<int> NodeNumberRanks2( NodeNumberRanks );
            sumNodeVector(NodeNumberRanks);
            for( int k = 0; k < node2ranks.size(); k++)
                NodeNumberRanks[k] = NodeNumberRanks[k] / (int)node2ranks[k].size();
            
            for( int k = 0; k < edge2ranks.size(); k++)
                EdgeNumberRanks[k] = (int)edge2ranks[k].size();
            ManagedArray<int> EdgeNumberRanks2( EdgeNumberRanks ); 
            sumEdgeVector(EdgeNumberRanks);
            
            for( int k = 0; k < edge2ranks.size(); k++)
                EdgeNumberRanks[k] = EdgeNumberRanks[k] / (int)edge2ranks[k].size();
            
            for( int k = 0; k < face2ranks.size(); k++)
                FaceNumberRanks[k] = (int)face2ranks[k].size();
            ManagedArray<int> FaceNumberRanks2( FaceNumberRanks ); 
            sumFaceVector(FaceNumberRanks);
            for( int k = 0; k < face2ranks.size(); k++)
                FaceNumberRanks[k] = FaceNumberRanks[k] / (int)face2ranks[k].size();
            
            assert(NodeNumberRanks == NodeNumberRanks2);
            assert(EdgeNumberRanks == EdgeNumberRanks2);
            assert(FaceNumberRanks == FaceNumberRanks2);
            
  
            
            printf( "test Nodes Communication on rank %d \n", (int)comm_.rank());
            

            ManagedArray<double> pcoord1(grid_.coordinates.size() );
            for( size_t k = 0; k < node2ranks.size(); ++k){
                pcoord1[k] = grid_.coordinates(k,0);
                if( node2ranks[k].size() > 0 ){
                    pcoord1[k] =  pcoord1[k] / (double)node2ranks[k].size();
                }
            }
            
            sumNodeVector( pcoord1 );
            
            for( int k = 0; k < (int)grid_.coordinates.size(); k++)
                if( pcoord1[k] - grid_.coordinates(k,0)  > 1E-12 ){
                    printf( "Vertex Communication failed %d,  rank = %d\n",(int)k, rank() );
                    node2ranks[k].print();
                    printf( "%lf %lf %lf\n", grid_.coordinates(k,0), grid_.coordinates(k,1) , grid_.coordinates(k,2));
                    returnBool = false;
                }
                
            ManagedArray<double> pcoord2(grid_.coordinates.size() );
            for( size_t k = 0; k < node2ranks.size(); ++k){
                pcoord2[k] = grid_.coordinates(k,1);
                if( node2ranks[k].size() > 0 )
                    pcoord2[k] =  pcoord2[k] / (double)node2ranks[k].size();
            }
            
            sumNodeVector( pcoord2 );
            
            for( int k = 0; k < (int)grid_.coordinates.size(); k++)
                if( pcoord2[k] - grid_.coordinates(k,1)  > 1E-12 ){                     
                    printf( "Vertex Communication failed %d,  rank = %d\n",(int)k, rank() );
                    node2ranks[k].print();
                    printf( "%lf %lf %lf\n", grid_.coordinates(k,0), grid_.coordinates(k,1) , grid_.coordinates(k,2));
                    returnBool = false;
                }
                
            ManagedArray<double> pcoord3(grid_.coordinates.size() );
            for( size_t k = 0; k < node2ranks.size(); ++k){
                pcoord3[k] = grid_.coordinates(k,2);
                if( node2ranks[k].size() > 0 )
                    pcoord3[k] =  pcoord3[k] / (double)node2ranks[k].size();
            }
            
            sumNodeVector( pcoord3 );
            
            for( int k = 0; k < (int)grid_.coordinates.size(); k++)
                if( pcoord3[k] - grid_.coordinates(k,2)  > 1E-12 ){         
                    printf( "Vertex Communication failed %d,  rank = %d\n",(int)k, rank() );
                    node2ranks[k].print();
                    returnBool = false;
                    printf( "%lf %lf %lf\n", grid_.coordinates(k,0), grid_.coordinates(k,1) , grid_.coordinates(k,2));
                }
     
            printf( "test Edge Communication on rank %d \n", (int)comm_.rank());
 
            ManagedArray<double> mpE1( grid_.nEdges(), 0. );
            ManagedArray<double> mpE2( grid_.nEdges(), 0. );
            ManagedArray<double> mpE3( grid_.nEdges(), 0. );
            for( size_t k = 0; k < grid_.nEdges(); ++k){
                const int e0 = grid_.edges(k,0);
                const int e1 = grid_.edges(k,1);
                mpE1[k] = .5*( grid_.coordinates(e0,0) + grid_.coordinates(e1,0) );
                mpE2[k] = .5*( grid_.coordinates(e0,1) + grid_.coordinates(e1,1) );
                mpE3[k] = .5*( grid_.coordinates(e0,2) + grid_.coordinates(e1,2) );
                if( edge2ranks[k].size() > 1 ){
                    mpE1[k] =  mpE1[k] / (double)edge2ranks[k].size();
                    mpE2[k] =  mpE2[k] / (double)edge2ranks[k].size();
                    mpE3[k] =  mpE3[k] / (double)edge2ranks[k].size();
                }
            }
            sumEdgeVector( mpE1 );
            sumEdgeVector( mpE2 );
            sumEdgeVector( mpE3 );
            
            ManagedArray<double> mpE12( grid_.nEdges(), 0. );
            ManagedArray<double> mpE22( grid_.nEdges(), 0. );
            ManagedArray<double> mpE32( grid_.nEdges(), 0. );
            for( size_t k = 0; k < grid_.nEdges(); ++k){
                const int e0 = grid_.edges(k,0);
                const int e1 = grid_.edges(k,1);
                mpE12[k] = .5*( grid_.coordinates(e0,0) + grid_.coordinates(e1,0) );
                mpE22[k] = .5*( grid_.coordinates(e0,1) + grid_.coordinates(e1,1) );
                mpE32[k] = .5*( grid_.coordinates(e0,2) + grid_.coordinates(e1,2) );
            }
            for( size_t k = 0; k < grid_.nEdges(); ++k){
                if( mpE1[k] - mpE12[k]  > 1E-12 ){
                     printf( "Edge Communication failed %d\n",(int)k);
                     returnBool = false;
                }
                if( mpE2[k] - mpE22[k]  > 1E-12 ){
                     printf( "Edge Communication failed %d\n",(int)k);
                     returnBool = false;
                }
                if( mpE3[k] - mpE32[k]  > 1E-12 ){
                     printf( "Edge Communication failed %d\n",(int)k);
                     returnBool = false;
                }
            }
            
            printf( "test face Communication on rank %d \n", (int)comm_.rank()); 
              
            ManagedArray<double> mpF1( grid_.nFaces(), 0. );
            ManagedArray<double> mpF2( grid_.nFaces(), 0. );
            ManagedArray<double> mpF3( grid_.nFaces(), 0. );
            for( size_t k = 0; k < grid_.nFaces(); ++k){
                const int f0 = grid_.faces(k,0);
                const int f1 = grid_.faces(k,1);
                const int f2 = grid_.faces(k,2);
                mpF1[k] = (1./3.)*(grid_.coordinates(f0,0) + grid_.coordinates(f1,0) + grid_.coordinates(f2,0));
                mpF2[k] = (1./3.)*(grid_.coordinates(f0,1) + grid_.coordinates(f1,1) + grid_.coordinates(f2,1));
                mpF3[k] = (1./3.)*(grid_.coordinates(f0,2) + grid_.coordinates(f1,2) + grid_.coordinates(f2,2));
                if( face2ranks[k].size() > 1 ){
                    mpF1[k] =  mpF1[k] / (double)face2ranks[k].size();
                    mpF2[k] =  mpF2[k] / (double)face2ranks[k].size();
                    mpF3[k] =  mpF3[k] / (double)face2ranks[k].size();
                }
            }
            
            sumFaceVector( mpF1 );
            sumFaceVector( mpF2 );
            sumFaceVector( mpF3 );
              
            ManagedArray<double> mpF12( grid_.nFaces(), 0. );
            ManagedArray<double> mpF22( grid_.nFaces(), 0. );
            ManagedArray<double> mpF32( grid_.nFaces(), 0. );  
               for( size_t k = 0; k < grid_.nFaces(); ++k){
                const int f0 = grid_.faces(k,0);
                const int f1 = grid_.faces(k,1);
                const int f2 = grid_.faces(k,2);
                mpF12[k] = (1./3.)*(grid_.coordinates(f0,0) + grid_.coordinates(f1,0) + grid_.coordinates(f2,0));
                mpF22[k] = (1./3.)*(grid_.coordinates(f0,1) + grid_.coordinates(f1,1) + grid_.coordinates(f2,1));
                mpF32[k] = (1./3.)*(grid_.coordinates(f0,2) + grid_.coordinates(f1,2) + grid_.coordinates(f2,2));
            }
            
            for( size_t k = 0; k < grid_.nFaces(); ++k){
                if( mpF1[k] - mpF12[k]  > 1E-12 ){
                     printf( "Face Communication failed %d\n",(int)k);
                     returnBool = false;
                }
                if( mpF2[k] - mpF22[k]  > 1E-12 ){
                     printf( "Face Communication failed %d\n",(int)k);
                     returnBool = false;
                }
                if( mpF3[k] - mpF32[k]  > 1E-12 ){
                     printf( "Face Communication failed %d\n",(int)k);
                     returnBool = false;
                }
            }
            
            if(returnBool)
                printf( "... done testing communication %d \n", (int)comm_.rank());   
            else
                printf("error in node communication %d \n",(int)comm_.rank());
            
            assert( returnBool );
            
            return returnBool;
        }
        
        void barrier() { comm_.barrier(); }
 
    private:
        Grid &grid_;
        Communicator &comm_;
        bool communication_is_valid;
    };

    template<class Communicator, class Entity2rank>
    void pack_entity( Communicator & comm_, void *buffer, std::size_t bufferSize, int &position, Entity2rank & entity2ranks )
    {
        int nEntities = entity2ranks.size();
        std::vector<char> entity2ranks_size( entity2ranks.size() );
        int entity_buffer_size = 0;
        for( size_t i = 0; i < entity2ranks.size(); i++ ){
            entity2ranks_size[i] = (char)entity2ranks[i].size();
            entity_buffer_size += entity2ranks[i].size();
        }
        std::vector<char> entity2ranks_buffer( entity_buffer_size );
        std::size_t ctr_buffer = 0;
        for( size_t i = 0; i < entity2ranks.size(); i++ )
            for( size_t k = 0; k < entity2ranks[i].size(); k++ )
                entity2ranks_buffer[ctr_buffer++] = entity2ranks[i][k];
            
        comm_.pack( &nEntities, 1, buffer, bufferSize, position );
        comm_.pack( &entity_buffer_size, 1, buffer, bufferSize, position );
        comm_.pack( entity2ranks_size.data(), nEntities, buffer, bufferSize, position );  
        comm_.pack( entity2ranks_buffer.data(), entity_buffer_size, buffer, bufferSize, position ); 
    }

    template<class Communicator, class Entity2rank>
    void unpack_entity( Communicator & comm_,  void *buffer, std::size_t bufferSize, int &position, Entity2rank & entity2ranks )
    {
        int nEntities;
        int entity_buffer_size;
        comm_.unpack( &nEntities, 1, buffer, bufferSize, position);
        comm_.unpack( &entity_buffer_size, 1, buffer, bufferSize, position);

        std::vector<char> entity2ranks_size(nEntities);
        std::vector<char> entity2ranks_buffer(entity_buffer_size);
        comm_.unpack( entity2ranks_size.data(), nEntities, buffer, bufferSize, position);  
        comm_.unpack( entity2ranks_buffer.data(), entity_buffer_size, buffer, bufferSize, position); 
        
        entity2ranks.resize(nEntities);
        std::size_t ctr_buffer = 0;
        for(int k = 0; k < (int)entity2ranks.size(); ++k){
            entity2ranks[k].resize(entity2ranks_size[k]);
            for(int j = 0; j < (int)entity2ranks_size[k]; ++j)
                entity2ranks[k][j] = entity2ranks_buffer[ctr_buffer++];
        }
    }

    template<class Communicator, class Entity2rank>
    std::size_t bufferSize_entity( Communicator & comm_, Entity2rank & entity2ranks )
    {   
        std::size_t bufferSize = comm_.template bufferSize<int>( 2  );
        bufferSize += comm_.template bufferSize<char>( entity2ranks.size());
        for(size_t i = 0; i < entity2ranks.size(); i++ ){
                bufferSize += comm_.template bufferSize<char>(entity2ranks[i].size());
        }
        return bufferSize;
    }

    template< class Decomposition, class Communicator>
    void communicateDecompositionData( Decomposition & decomposition, Communicator &comm ) {
        communicateDecompositionData_opt( decomposition, comm, decomposition.node2ranks, decomposition.edge2ranks, decomposition.face2ranks);
    }
    
    template< class Decomposition, class Communicator, class NodeArray, class EdgeArray, class FaceArray >
    void communicateDecompositionData( Decomposition & decomposition, Communicator &comm, 
            NodeArray &node2ranks, EdgeArray &edge2ranks, FaceArray &face2ranks ) {
        
        std::vector< ManagedArray<int> > buffer_send( comm.size() );
        std::vector< typename Communicator::Request > requests_send(comm.size());

        for( int p = 0; p < comm.size(); p++)
            if( decomposition.rank2bdyNodes[p].size() > 0 ) {
                std::size_t bufferSize = 0;
                
                bufferSize += decomposition.rank2bdyNodes[p].size();
                for( size_t k = 0; k < decomposition.rank2bdyNodes[p].size(); k++)
                    bufferSize += node2ranks[ decomposition.rank2bdyNodes[p][k] ].size();
                
                bufferSize += decomposition.rank2bdyEdges[p].size();
                for( size_t k = 0; k < decomposition.rank2bdyEdges[p].size(); k++)
                    bufferSize += edge2ranks[ decomposition.rank2bdyEdges[p][k] ].size();
                
                bufferSize += decomposition.rank2bdyFaces[p].size();
                for( size_t k = 0; k < decomposition.rank2bdyFaces[p].size(); k++)
                    bufferSize += face2ranks[ decomposition.rank2bdyFaces[p][k] ].size();
                
                buffer_send[p].resize(bufferSize, 0);

                size_t counter = 0;
                for( size_t k = 0; k < decomposition.rank2bdyNodes[p].size(); k++)
                    buffer_send[p][counter++] = (char)node2ranks[ decomposition.rank2bdyNodes[p][k] ].size();
                for( size_t k = 0; k < decomposition.rank2bdyNodes[p].size(); k++)
                    for( size_t j = 0; j < node2ranks[ decomposition.rank2bdyNodes[p][k] ].size(); j++)
                        buffer_send[p][counter++] = node2ranks[ decomposition.rank2bdyNodes[p][k] ][j];
                                  
                for( size_t k = 0; k < decomposition.rank2bdyEdges[p].size(); k++)
                    buffer_send[p][counter++] = (char)edge2ranks[ decomposition.rank2bdyEdges[p][k] ].size();
                for( size_t k = 0; k < decomposition.rank2bdyEdges[p].size(); k++)
                    for( size_t j = 0; j < edge2ranks[ decomposition.rank2bdyEdges[p][k] ].size(); j++)
                        buffer_send[p][counter++] = edge2ranks[ decomposition.rank2bdyEdges[p][k] ][j];
            
                for( size_t k = 0; k < decomposition.rank2bdyFaces[p].size(); k++)
                    buffer_send[p][counter++] = (char)face2ranks[ decomposition.rank2bdyFaces[p][k] ].size();    
                for( size_t k = 0; k < decomposition.rank2bdyFaces[p].size(); k++)
                    for( size_t j = 0; j < face2ranks[ decomposition.rank2bdyFaces[p][k] ].size(); j++)
                        buffer_send[p][counter++] = face2ranks[ decomposition.rank2bdyFaces[p][k] ][j];
                 
                comm.Isend( buffer_send[p].data(), buffer_send[p].size(), p, 1, &requests_send[p] );
            }
        
        for( size_t p = 0; p < (size_t)comm.size(); p++)
            if( decomposition.rank2bdyNodes[p].size() > 0 ) {
                typename Communicator::Status status;
                comm.probe( p, 1, &status );
                std::size_t nRecv = comm.template getCount<int>( &status );
                ManagedArray<int> buffer_recv( nRecv );
                comm.recv( buffer_recv.data(), buffer_recv.size(), p, 1);
        
                size_t counter = 0;
                ManagedArray<char> sizes_nodes(decomposition.rank2bdyNodes[p].size(),0);
                for( size_t k = 0; k < decomposition.rank2bdyNodes[p].size(); k++){
                    sizes_nodes[k] = buffer_recv[counter++];
                    node2ranks[k].reserve( node2ranks[k].size() + sizes_nodes[k] ); 
                }
                for( size_t k = 0; k < decomposition.rank2bdyNodes[p].size(); k++)
                    for( int j = 0; j < sizes_nodes[k]; j++)
                        node2ranks[ decomposition.rank2bdyNodes[p][k] ].push_back( buffer_recv[counter++] );
                        
                ManagedArray<char> sizes_edges(decomposition.rank2bdyEdges[p].size(),0);
                for( size_t k = 0; k < decomposition.rank2bdyEdges[p].size(); k++){
                    sizes_edges[k] = buffer_recv[counter++];
                    edge2ranks[k].reserve( edge2ranks[k].size() + sizes_edges[k] ); 
                }
                for( size_t k = 0; k < decomposition.rank2bdyEdges[p].size(); k++)
                    for( int j = 0; j < sizes_edges[k]; j++)
                        edge2ranks[ decomposition.rank2bdyEdges[p][k] ].push_back( buffer_recv[counter++] );
                    
                ManagedArray<char> sizes_faces(decomposition.rank2bdyFaces[p].size(),0);
                for( size_t k = 0; k < decomposition.rank2bdyFaces[p].size(); k++){
                    sizes_faces[k] = buffer_recv[counter++];
                    face2ranks[k].reserve( face2ranks[k].size() + sizes_faces[k] ); 
                }
                for( size_t k = 0; k < decomposition.rank2bdyFaces[p].size(); k++)
                    for( int j = 0; j < sizes_faces[k]; j++)
                        face2ranks[ decomposition.rank2bdyFaces[p][k] ].push_back( buffer_recv[counter++] );    
   
                assert(counter == buffer_recv.size());               
            }
         
        for(int p = 0; p < (int)comm.size(); p++)
            if( decomposition.rank2bdyNodes[p].size() > 0 )
                comm.wait(&requests_send[p]);

        for( int k = 0; k < (int)node2ranks.size(); ++k )
            vector_unique(node2ranks[k].array_);
        for( int k = 0; k < (int)edge2ranks.size(); ++k )
            vector_unique(edge2ranks[k].array_);
        for( int k = 0; k < (int)face2ranks.size(); ++k )
            vector_unique(face2ranks[k].array_);            
    }
#if 0
    template< class Decomposition, class Communicator, class NodeArray, class EdgeArray, class FaceArray>
    class DecompositionContainer
    {
    public:
        DecompositionContainer(Decomposition & decomposition, Communicator &comm, 
            NodeArray &node2ranks, EdgeArray &edge2ranks, FaceArray &face2ranks ) :
            decomposition_(decomposition), comm_(comm), node2ranks_(node2ranks), edge2ranks_(edge2ranks), face2ranks_(face2ranks)
            {}
        
        void pack( void *buffer, std::size_t bufferSize, int &position ) {
            //int message_size = (int)u_buffer.size();
            //comm_.pack( &message_size, 1, buffer, bufferSize, position );
            //comm_.pack( u_buffer.data(), u_buffer.size(), buffer, bufferSize, position );
        }
    
        void unpack( void *buffer, std::size_t bufferSize, int &position ) {
            //int message_size;
            //comm.unpack( &message_size, 1, buffer, bufferSize, position );
            //u_buffer.resize( message_size );
            //comm_.unpack( u_buffer.data(), u_buffer.size(), buffer, bufferSize, position );
        }
        
        std::size_t bufferSize() {
            //return comm_.template bufferSize<int>(1)
            //    + comm_.template bufferSize<double>(u_buffer.size());     
        }
    private:
        Decomposition & decomposition_;
        Communicator &comm_;
        NodeArray &node2ranks_;
        EdgeArray &edge2ranks_;
        FaceArray &face2ranks_;
        
    }
#endif  
    template< class Decomposition, class Communicator, class NodeArray, class EdgeArray, class FaceArray >
    void communicateDecompositionData_opt( Decomposition & decomposition, Communicator &comm, 
            NodeArray &node2ranks, EdgeArray &edge2ranks, FaceArray &face2ranks ) {
        
        std::vector< ManagedArray<int> > buffer_send( comm.size() );
        std::vector< typename Communicator::Request > requests_send(comm.size());
        
        std::size_t nMessages = 0;
        
        for( int p = 0; p < comm.size(); p++)
            if( decomposition.rank2bdyNodes[p].size() > 0 ) {
                std::size_t bufferSize = 0;
                
                bufferSize += decomposition.rank2bdyNodes[p].size();
                for( size_t k = 0; k < decomposition.rank2bdyNodes[p].size(); k++)
                    bufferSize += node2ranks[ decomposition.rank2bdyNodes[p][k] ].size();
                
                bufferSize += decomposition.rank2bdyEdges[p].size();
                for( size_t k = 0; k < decomposition.rank2bdyEdges[p].size(); k++)
                    bufferSize += edge2ranks[ decomposition.rank2bdyEdges[p][k] ].size();
                
                bufferSize += decomposition.rank2bdyFaces[p].size();
                for( size_t k = 0; k < decomposition.rank2bdyFaces[p].size(); k++)
                    bufferSize += face2ranks[ decomposition.rank2bdyFaces[p][k] ].size();
                
                buffer_send[p].resize(bufferSize, 0);

                size_t counter = 0;
                for( size_t k = 0; k < decomposition.rank2bdyNodes[p].size(); k++)
                    buffer_send[p][counter++] = (char)node2ranks[ decomposition.rank2bdyNodes[p][k] ].size();
                for( size_t k = 0; k < decomposition.rank2bdyNodes[p].size(); k++)
                    for( size_t j = 0; j < node2ranks[ decomposition.rank2bdyNodes[p][k] ].size(); j++)
                        buffer_send[p][counter++] = node2ranks[ decomposition.rank2bdyNodes[p][k] ][j];
                                  
                for( size_t k = 0; k < decomposition.rank2bdyEdges[p].size(); k++)
                    buffer_send[p][counter++] = (char)edge2ranks[ decomposition.rank2bdyEdges[p][k] ].size();
                for( size_t k = 0; k < decomposition.rank2bdyEdges[p].size(); k++)
                    for( size_t j = 0; j < edge2ranks[ decomposition.rank2bdyEdges[p][k] ].size(); j++)
                        buffer_send[p][counter++] = edge2ranks[ decomposition.rank2bdyEdges[p][k] ][j];
            
                for( size_t k = 0; k < decomposition.rank2bdyFaces[p].size(); k++)
                    buffer_send[p][counter++] = (char)face2ranks[ decomposition.rank2bdyFaces[p][k] ].size();    
                for( size_t k = 0; k < decomposition.rank2bdyFaces[p].size(); k++)
                    for( size_t j = 0; j < face2ranks[ decomposition.rank2bdyFaces[p][k] ].size(); j++)
                        buffer_send[p][counter++] = face2ranks[ decomposition.rank2bdyFaces[p][k] ][j];
                 
                comm.Isend( buffer_send[p].data(), buffer_send[p].size(), p, 1, &requests_send[p] );
                
                nMessages++;
            }
        
        
        while(nMessages > 0){
            std::size_t nRecv = 0;
            int p;
            int tag = 1;
            if( comm.Iprobe_any_int( p, tag, nRecv) ) {
                ManagedArray<int> buffer_recv( nRecv );
                comm.recv( buffer_recv.data(), buffer_recv.size(), p, 1);
        
                size_t counter = 0;
                ManagedArray<char> sizes_nodes(decomposition.rank2bdyNodes[p].size(),0);
                for( size_t k = 0; k < decomposition.rank2bdyNodes[p].size(); k++){
                    sizes_nodes[k] = buffer_recv[counter++];
                    node2ranks[k].reserve( node2ranks[k].size() + sizes_nodes[k] ); 
                }
                for( size_t k = 0; k < decomposition.rank2bdyNodes[p].size(); k++)
                    for( int j = 0; j < sizes_nodes[k]; j++)
                        node2ranks[ decomposition.rank2bdyNodes[p][k] ].push_back( buffer_recv[counter++] );
                        
                ManagedArray<char> sizes_edges(decomposition.rank2bdyEdges[p].size(),0);
                for( size_t k = 0; k < decomposition.rank2bdyEdges[p].size(); k++){
                    sizes_edges[k] = buffer_recv[counter++];
                    edge2ranks[k].reserve( edge2ranks[k].size() + sizes_edges[k] ); 
                }
                for( size_t k = 0; k < decomposition.rank2bdyEdges[p].size(); k++)
                    for( int j = 0; j < sizes_edges[k]; j++)
                        edge2ranks[ decomposition.rank2bdyEdges[p][k] ].push_back( buffer_recv[counter++] );
                    
                ManagedArray<char> sizes_faces(decomposition.rank2bdyFaces[p].size(),0);
                for( size_t k = 0; k < decomposition.rank2bdyFaces[p].size(); k++){
                    sizes_faces[k] = buffer_recv[counter++];
                    face2ranks[k].reserve( face2ranks[k].size() + sizes_faces[k] ); 
                }
                for( size_t k = 0; k < decomposition.rank2bdyFaces[p].size(); k++)
                    for( int j = 0; j < sizes_faces[k]; j++)
                        face2ranks[ decomposition.rank2bdyFaces[p][k] ].push_back( buffer_recv[counter++] );    
                
                assert(counter == buffer_recv.size());
                nMessages--;
            }
        }
         
        for(int p = 0; p < (int)comm.size(); p++)
            if( decomposition.rank2bdyNodes[p].size() > 0 )
                comm.wait(&requests_send[p]);

        for( int k = 0; k < (int)node2ranks.size(); ++k )
            vector_unique(node2ranks[k].array_);
        for( int k = 0; k < (int)edge2ranks.size(); ++k )
            vector_unique(edge2ranks[k].array_);
        for( int k = 0; k < (int)face2ranks.size(); ++k )
            vector_unique(face2ranks[k].array_);            
    }
    
    
    template< class Grid, class Communicator>
    class Entity2ranks {
        typedef  Entity2ranks<Grid,Communicator> ThisType;
    public:
        
        Entity2ranks( Grid &grid, Communicator &comm ) : grid_(grid), comm_(comm) {}
        
        void resize(std::size_t size_) { data.resize(size_); }
        
        size_t nRanks( std::size_t entity ){ return data[entity].count(); }
        
        void set( std::size_t entity, std::size_t rank ) {
            data[entity].set(rank);  
        }
        
        void reset( std::size_t entity, std::size_t rank ) {
            data[entity].reset(rank);  
        }
        
        bool test(std::size_t entity, std::size_t rank) {
            return data[entity].test(rank);
        }
        
        void pack( void *buffer, std::size_t bufferSize, int &position ) {}
        
        void unpack( void *buffer, std::size_t bufferSize, int &position ) {}
        
        std::size_t bufferSize() { return 0; }
        
    private:
        std::vector<char> neighboring_ranks_;
        std::vector< std::bitset<8> > data;
        Grid &grid_;
        Communicator &comm_;
    };
    
} // namespace conformingsimplexgrid
#endif // CS_DECOMPOSITION_HH

