#ifndef CS_MPICOMMUNICATION_HH
#define CS_MPICOMMUNICATION_HH

#include "mpi.h"

#include "../mesh.hh"
#include "../parallel.hh"


namespace conformingsimplexgrid {
 
   template<class T>
   struct mpi_traits;
   
   template<>
   struct mpi_traits<double>
   { static MPI_Datatype type() { return MPI_DOUBLE; } };
   
   template<>
   struct mpi_traits<int>
   { static MPI_Datatype type() { return MPI_INT; } };
   
   template<>
   struct mpi_traits<unsigned int>
   { static MPI_Datatype type() { return MPI_UNSIGNED; } };
   
   template<>
   struct mpi_traits<char>
   { static MPI_Datatype type() { return MPI_CHAR; } };
   
   template<>
   struct mpi_traits<unsigned char>
   { static MPI_Datatype type() { return MPI_UNSIGNED_CHAR; } };
///@endcond   
   //template<>
   //struct mpi_traits<bool>
   //{ static MPI_Datatype type() { return MPI_BOOL; } };
    
   class MpiCommunicator {
   
   public:
    typedef MPI_Request Request;
    typedef MPI_Status Status;
       
    MpiCommunicator (int *argc, char ***argv) {  
        MPI_Init( argc, argv );
        MPI_Comm_rank ( MPI_COMM_WORLD, &rank_ );
        MPI_Comm_size ( MPI_COMM_WORLD, &numP_ );
    }
    
    void finalize () {
        MPI_Finalize();
    }
    
    int rank() {
        return rank_;
    }
    
    bool master() {
        return rank_ == 0;
    }
    
    int size()  {
        return numP_;
    }
    
    template <class T>
    int send( T *buffer, int size, int dest, int tag ) { 
        return MPI_Send( buffer, size, mpi_traits<T>::type(), dest, tag, MPI_COMM_WORLD );
    }
    
    template <class T>
    int recv( T *buffer, int size, int source, int tag ) {
        return MPI_Recv( buffer, size, mpi_traits<T>::type(), source, tag, MPI_COMM_WORLD, 0 );
    }
    
    template< class T >
    int Isend( T *buffer, int size, int dest, int tag, Request *request) {
        return MPI_Isend( buffer, size, mpi_traits<T>::type(), dest, tag, MPI_COMM_WORLD, request );
    }
    
    template< class T >
    int Isend_Byte( T *buffer, int size, int dest, int tag, Request *request) {
        return MPI_Isend( buffer, size, MPI_BYTE, dest, tag, MPI_COMM_WORLD, request );
    }
    
    template< class T >
    int Irecv( T *buffer, int size, int source, int tag, Request *request) {
        return MPI_Irecv( buffer, size, mpi_traits<T>::type(), source, tag, MPI_COMM_WORLD, request );
    }
    
    template< class T >
    int Irecv_Byte( T *buffer, int size, int source, int tag, Request *request) {
        return MPI_Irecv( buffer, size, MPI_BYTE, source, tag, MPI_COMM_WORLD, request );
    }
    
    int Isend_packed( void *buffer, size_t bufferSize, int dest, int tag, Request *request ){
        return MPI_Isend( buffer, bufferSize, MPI_PACKED, dest, tag, MPI_COMM_WORLD, request );
    }
    
    int recv_packed( void *buffer, size_t bufferSize, int source, int tag ){
        return MPI_Recv( buffer, bufferSize, MPI_PACKED, source, tag, MPI_COMM_WORLD, 0 );
    }
    
    int Irecv_packed( void *buffer, size_t bufferSize, int source, int tag, Request *request ){
        return MPI_Irecv( buffer, bufferSize, MPI_PACKED, source, tag, MPI_COMM_WORLD, request );
    }
    
    int probe( int source, int tag, Status *status ){
        return MPI_Probe( source, tag, MPI_COMM_WORLD, status );   
    }
    
    bool Iprobe( int &source, int &tag, std::size_t &bufferSize ){
        Status status;   
        int flag = 0;
        MPI_Iprobe( source, tag, MPI_COMM_WORLD, &flag, &status);
        if( !flag )
            return false;
        bufferSize = getCount_packed( &status );
        return true;
    }
    
    bool Iprobe_any( int &source, int &tag, std::size_t &bufferSize ){
        Status status;
        int flag = 0;
        MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        if( !flag )
            return false;
        source = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        bufferSize = getCount_packed( &status );
        return true;
    }
    
    bool Iprobe_any_int( int &source, int &tag, std::size_t &bufferSize ){
        Status status;
        int flag = 0;
        MPI_Iprobe( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
        if( !flag )
            return false;
        source = status.MPI_SOURCE;
        tag = status.MPI_TAG;
        bufferSize = getCount<int>( &status );
        return true;
    }

    template< class T>
    std::size_t getCount( Status *status ) {
        int count;
        MPI_Get_count( status, mpi_traits<T>::type(), &count );
        return count;
    }
    
    template< class T>
    std::size_t getCount_Byte( Status *status ) {
        int count;
        MPI_Get_count( status, MPI_BYTE, &count );
        return count;
    }
    
    std::size_t getCount_packed( Status *status ) {
        int count;
        MPI_Get_count( status, MPI_PACKED, &count );
        return count;
    }
    
    template< class T >
    int pack( T *data, std::size_t dataSize, void *buffer, std::size_t bufferSize, int &position ) 
    {
        //printf("rank %d --> packing %d , position %d \n", (int)rank(), (int)dataSize, position);
        return MPI_Pack( data, dataSize, mpi_traits<T>::type(), buffer, bufferSize, &position, MPI_COMM_WORLD );
    }
    
    template< class T >
    int unpack( T *data, std::size_t dataSize, void *buffer, std::size_t bufferSize, int &position ) 
    {
        //printf("rank %d --> unpacking %d , position %d \n", (int)rank(), (int)dataSize, position);
        return MPI_Unpack( buffer, bufferSize, &position, data, dataSize, mpi_traits<T>::type(), MPI_COMM_WORLD );
    }
    
    template< class T >
    std::size_t bufferSize( std::size_t dataSize ){
        int size(0);
        MPI_Pack_size( dataSize, mpi_traits<T>::type(), MPI_COMM_WORLD, &size);  
        return size;
    }
    
    int wait( Request *request ) {
     	return MPI_Wait( request, MPI_STATUS_IGNORE );
    }
    
    int test( Request *request, int *flag ) {
    	return MPI_Test( request, flag, MPI_STATUS_IGNORE );
    }
    
    template< class T>
    int reduce_sum( T *sendbuf, T *recvbuf, int size ) {
    	return MPI_Allreduce( sendbuf, recvbuf, size, mpi_traits<T>::type(), MPI_SUM, MPI_COMM_WORLD);
    }
    
    template< class T>
    int reduce_min( T *sendbuf, T *recvbuf, int size ) {
        return MPI_Allreduce( sendbuf, recvbuf, size, mpi_traits<T>::type(), MPI_MIN, MPI_COMM_WORLD);
    }
    
    template< class T>
    int reduce_max( T *sendbuf, T *recvbuf, int size ) {
        return MPI_Allreduce( sendbuf, recvbuf, size, mpi_traits<T>::type(), MPI_MAX, MPI_COMM_WORLD);
    }
    
    template< class T>
    int all2all( T *sendbuf, T *recvbuf, int size ) {
        return MPI_Alltoall( sendbuf, size, mpi_traits<T>::type(), recvbuf, size, mpi_traits<T>::type(), MPI_COMM_WORLD );
    }
    
    void barrier() const { MPI_Barrier( MPI_COMM_WORLD ); }
    
   private:
    int rank_;
    int numP_;
   };
   
   //###################################################################################
#if 0
   template< class Mesh, class Communicator >
   class faceCommunication{
   
   public:
    explicit faceCommunication( Mesh &mesh, ManagedArray<int> &bdyFace2Rank, Communicator &mpiComm )
    : mesh_(mesh), bdyFace2Rank_(bdyFace2Rank), mpiComm_(mpiComm)
    {  
        const int NV = Mesh::verticesPerElement;
        
        ManagedArray<int> rankMap( mpiComm.size() );
        
        
        int count = 1;
        rankMap[ bdyFace2Rank[0] ] = 0;
        neighborRanks.push_back( bdyFace2Rank[0]);
        for(int i = 1; i < bdyFace2Rank_.size(); i++){
            if( bdyFace2Rank[i-1] != bdyFace2Rank[i] ) {
                ManagedArray<int> nodes_here( count * (NV-1) );
                commNodes.push_back(nodes_here);
                rankMap[ bdyFace2Rank[i] ] = commNodes.size();
                neighborRanks.push_back( bdyFace2Rank[i] );
                count = 1;
            } else {
                count++;
            }
        } 

        ManagedArray<int> nodes_here( count*(NV-1) );
        commNodes.push_back(nodes_here);
        
        ManagedArray<int> counter( commNodes.size(), 0);
         for(int el = 0; el < mesh_.n4e.size(); el++) {
             for(int k = 0; k < NV; k++ ) {
                 if( !mesh_.neighbors.isParallel(el,k) )
                     continue;
                           
                 const int bdyNum = mesh_.neigh(el,k);
                 const int rank = bdyFace2Rank[ bdyNum ];
                 const int rank_index = rankMap[rank];
                 
                 for(int j = 0; j < NV-1; j++)
                    commNodes[rank_index][counter[rank_index]*(NV-1) + j] = mesh_.n4e(el, (j < k) ? j : j + 1);
                 
                 counter[rank_index]++;
             }
         }

         
         for(int i = 0; i < commNodes.size(); i++)
             commNodes[i] = commNodes[i].unique_unsorted();
                  
            
    }

   private:
    Mesh &mesh_;
    ManagedArray<int> &bdyFace2Rank_;
    Communicator &mpiComm_;
    std::vector<int> neighborRanks;
    std::vector< ManagedArray<int> > commNodes;
   };   
   
   template< class Mesh, class Communicator >
   class coarseSpace{
       explicit coarseSpace( Mesh &mesh, ManagedArray<int> &bdyFace2Rank, Communicator &mpiComm )
        : mesh_(mesh), bdyFace2Rank_(bdyFace2Rank), mpiComm_(mpiComm)
        {}
       
   
   private:
        Mesh &mesh_;
        ManagedArray<int> &bdyFace2Rank_;
        Communicator &mpiComm_;
   };
#endif
//###################################################################################
} // namespace conformingsimplexgrid

#endif // CS_MPICOMMUNICATION_HH
