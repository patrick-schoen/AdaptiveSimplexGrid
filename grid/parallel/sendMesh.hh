#ifndef SEND_MESH_HH
#define SEND_MESH_HH

#include<functional>

#include "mpi.h"

#include "../simplexGrid.hh"
#include "../parallel.hh"


namespace conformingsimplexgrid {
       
    template< class Communicator >
    struct SentMesh 
    {        
        SentMesh( Communicator &comm, int rank, int tag, void* buffer, std::size_t bufferSize ) 
        : comm_(comm), buffer_(buffer)
        {
            comm_.get().Isend_packed( buffer_, bufferSize, rank, tag, &request_ ); 
        }
        
        ~SentMesh() { wait(); }
        
        SentMesh( const SentMesh & ) = delete;
        
        SentMesh( SentMesh && other)
        : comm_(other.comm_), buffer_(other.buffer_), request_( other.request_ )
        {
            other.buffer_ = nullptr;   
        }
        
        SentMesh & operator = (const SentMesh & ) = delete;
        
//         SentMesh & operator = ( SentMesh && other )
//         {
//             wait();
//             comm_ = other.comm_;
//             buffer_ = other.buffer_;
//             request_ = other.request_;
//             other.buffer_ = nullptr;
//         }

        bool pending () const { return bool( buffer_ ); }
        
        void wait ()
        {
            if( pending() ) {
                comm_.get().wait( &request_);   
                std::free(buffer_);
                buffer_ = nullptr;
            }
        }
        
    private: 
        std::reference_wrapper< Communicator > comm_;
        void * buffer_;
        typename Communicator::Request request_;
        
    };
   
    template < class Mesh, class Communicator, class ... Manager >
    SentMesh< Communicator > sendMesh( Mesh &mesh, int rank, int tag, Communicator &comm, Manager & ... manager) 
    {
        std::size_t bufferSize = 0;
                    
        bufferSize += mesh.bufferSize( comm );
        bufferSize += bufferSize_Manager( manager ... );
        
        void * buffer = std::malloc( bufferSize );
        int position = 0;
        
        mesh.pack( comm, buffer, bufferSize, position );
        pack_Manager( buffer, bufferSize, position, manager ... );
    
        return SentMesh< Communicator >( comm, rank, tag, buffer, position);      
    }
    
    template < class Mesh, class Communicator, class ... Manager >
    int recvMesh( Communicator &comm, Mesh & mesh, Manager & ... manager)
    {
        std::size_t bufferSize;
        int rank;
        int tag;
        
        if( !comm.Iprobe_any( rank, tag, bufferSize ) )
            return 0;

        void* buffer = std::malloc(bufferSize);
        comm.recv_packed( buffer, bufferSize, rank, tag  );
        
        int position = 0;
        
        mesh.unpack( comm, buffer, bufferSize, position );
        unpack_Manager( buffer, bufferSize, position, manager ... );
        
        std::free(buffer);
        
        getFacesForElements( mesh.elements, mesh.faces, mesh.element2faces, mesh.level );
        getFacesForElements( mesh.faces, mesh.edges, mesh.face2edges );
        
        return 1;
    }
    
    template < class Mesh, class Communicator, class ... Manager >
    void recvMesh( int rank, int tag, Communicator &comm, Mesh & mesh, Manager & ... manager)
    {
        typename Communicator::Status status;
        static const int NV = Mesh::verticesPerElement;
        static const int dim = Mesh::dimensionworld;

        comm.probe( rank, tag, &status);
        
        std::size_t bufferSize = comm.getCount_packed( &status );
        
        void* buffer = std::malloc(bufferSize);
        comm.recv_packed( buffer, bufferSize, rank, tag  );
        
        int position = 0;
        
        mesh.unpack( comm, buffer, bufferSize, position );
        unpack_Manager( buffer, bufferSize, position, manager ... );
        
        std::free(buffer);
        
        getFacesForElements( mesh.elements, mesh.faces, mesh.element2faces, mesh.level );
        getFacesForElements( mesh.faces, mesh.edges, mesh.face2edges );
    }

} // namespace conformingsimplexgrid

#endif // SEND_MESH_HH