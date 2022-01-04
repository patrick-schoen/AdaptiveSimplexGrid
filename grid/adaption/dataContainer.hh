#ifndef DATA_CONTAINER_HH
#define DATA_CONTAINER_HH

#include "../simplexGrid.hh"
#include "../parallel.hh"

// Standard data manager. Provides the minimal interface requiered by mesh's adaption routines

namespace conformingsimplexgrid {
      
    class DefaultContainer
    {
    public:
        template< class Pre_Refine_Data >
        void pre_refine( Pre_Refine_Data & ) {}
        
        template< class Refine_Data >
        void refine( Refine_Data & ) {}
    
        template< class Coarse_Data >
        void coarse( Coarse_Data & ) {}
        
        template< class Vector >
        void newIndices( Vector &, Vector &, Vector & ) {}
           
        size_t size() { return 0; }
        
        template< class BoolVector, class ... Vector > 
        void preparePacking( BoolVector &, Vector & ... ) {}
                
        void pack( void*, std::size_t, int& ) {}
        
        void unpack( void*, std::size_t, int& ) {}
        
        template< class ... Vector >
        void merge( Vector & ...  ) {}
        
        void conjoin() {}
            
        std::size_t bufferSize() { return 0; }  
    };
    
    template< class PrologationData, class FirstManager, class ... Manager >
    void pre_refine_Manager( PrologationData &data, FirstManager & first, Manager & ... manager ) {
        first.pre_refine( data );
        pre_refine_Manager( data, manager ... );
    };

    template< class PrologationData >
    void pre_refine_Manager( PrologationData & ) {};
    
    template< class PrologationData, class FirstManager, class ... Manager >
    void refine_Manager( PrologationData &data, FirstManager & first, Manager & ... manager ) {
        first.refine( data );
        refine_Manager( data, manager ... );
    };

    template< class PrologationData >
    void refine_Manager( PrologationData & ) {};

    template< class PrologationData, class FirstManager, class ... Manager >
    void coarse_Manager( PrologationData &data, FirstManager & first, Manager & ... manager ) {
        first.coarse( data );
        coarse_Manager( data, manager ... );
    };

    template< class PrologationData >
    void coarse_Manager( PrologationData & ) {};
    
    template< class Vector,class FirstManager, class ... Manager >
    void newIndices_Manager( Vector &faces2newFaces, Vector &edges2newEdges, Vector &nodes2newNodes, FirstManager & first, Manager & ... manager  ) {
        first.newIndices( faces2newFaces, edges2newEdges, nodes2newNodes );
        newIndices_Manager( faces2newFaces, edges2newEdges, nodes2newNodes, manager ... );
    }
    
    template< class Vector >
    void newIndices_Manager( Vector &, Vector &, Vector & ) {};
          
    void pack_Manager(void *, std::size_t, int &) {}
    
    template< class First, class ... Manager >
    void pack_Manager( void *buffer, std::size_t bufferSize, int &position, First & first, Manager & ... manager )
    {
        first.pack( buffer, bufferSize, position);
        pack_Manager( buffer, bufferSize, position, manager ... );
    }
  
    void unpack_Manager(void *, std::size_t, int &){}
    
    template< class First, class ... Manager >
    void unpack_Manager( void *buffer, std::size_t bufferSize, int &position, First & first, Manager & ... manager )
    {
        first.unpack( buffer, bufferSize, position);
        unpack_Manager( buffer, bufferSize, position, manager ... );
    }
  
    std::size_t bufferSize_Manager() { return 0; }
    
    template< class First, class ... Manager >
    std::size_t bufferSize_Manager( First & first, Manager & ... manager )
    {
        return first.bufferSize() + bufferSize_Manager( manager ... );
    }
    

} // namespace conformingsimplexgrid
    
#endif // DATA_CONTAINER_HH
