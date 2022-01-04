#ifndef MESSAGE_BUFFER_HH
#define MESSAGE_BUFFER_HH

#include<functional>
#include "mpi.h"

namespace conformingsimplexgrid{
    
    
    //template< Communicator &comm, Decomposition &decomp >
    //void getMessagesViaAllToAll() {}
    
    
    template< class Communicator, class ... Manager >
    struct SentBuffer 
    {        
        SentBuffer( Communicator &comm, int dest, int tag, Manager & ... manager ) 
        : comm_(comm)
        {
            std::size_t bufferSize = 0;        
            bufferSize += bufferSize_Manager( manager ... );
            buffer_ = std::malloc( bufferSize );
            std::size_t position = 0;
            pack_Manager( buffer_, bufferSize, position, manager ... );
            comm_.get().Isend_packed( buffer_, bufferSize, dest, tag, &request_ ); 
        }       

        ~SentBuffer() { wait(); }
        
        SentBuffer( const SentBuffer & ) = delete;
        
        SentMesh( SentMesh && other)
        : comm_(other.comm_), buffer_(other.buffer_),request_( other.request_ )
        {
            other.buffer_ = nullptr;   
        }
        
        SentBuffer & operator = (const SentBuffer & ) = delete;

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
    
    template< class Communicator >
    struct RecvBuffer 
    {        
        RecvBuffer( Communicator &comm ) 
        : comm_(comm)
        {
            haveMessageSize = false;
        }

        bool Iprobe( int source, int tag ) {
            if(haveMessageSize == false) {
                bool flag = comm.Iprobe( source, tag, bufferSize_ );
                if( flag ) buffer_ = std::malloc(bufferSize_);
                haveMessageSize = true;
                return flag;
            }
        }
   
        void Irecv( int source, int tag ) {
            assert( haveMessageSize == true );
            comm_.get().Irecv_packed( buffer_, bufferSize_, source, tag, &request_ );
        }
        
        template< class ... Manager>
        void unpack( Manager & ... manager )
        {
            std::size_t position = 0;
            unpack_Manager( buffer_, bufferSize_, position, manager ... );
            std::free(buffer);
        }
        
        ~RecvBuffer() { wait(); }
        
        RecvBuffer( const RecvBuffer & ) = delete;
        
        RecvMesh( SentMesh && other)
        : comm_(other.comm_), buffer_(other.buffer_),request_( other.request_ )
        {
            other.buffer_ = nullptr;   
        }
        
        RecvBuffer & operator = (const RecvBuffer & ) = delete;

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
        std::bool haveMessageSize;
        void * buffer_;
        std::size_t bufferSize_;
        typename Communicator::Request request_;
    };
    
#if 0
    template< class Communicator >
    struct RecvBuffer 
    {        
        RecvBuffer( Communicator &comm, int source, int tag, void* buffer, std::size_t bufferSize ) 
        : comm_(comm), buffer_(buffer)
        {
            comm_.get().Irecv_packed( buffer_, bufferSize, source, tag, &request_ ); 
        }
        
        ~RecvBuffer() { wait(); }
        
        RecvBuffer( const RecvBuffer & ) = delete;
        
        RecvBuffer( RecvBuffer && other)
        : comm_(other.comm_), buffer_(other.buffer_), request_( other.request_ )
        {
            other.buffer_ = nullptr;   
        }
        
        RecvBuffer & operator = ( const RecvBuffer & ) = delete;
      
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
    
    
#endif   
    
    
    
    
    
} //namespace conformingsimplexgrid

#endif //MESSAGE_BUFFER_HH