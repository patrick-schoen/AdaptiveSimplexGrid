#ifndef SPACE_FILLING_CURVE_HH
#define SPACE_FILLING_CURVE_HH

#include "../simplexGrid.hh"
#include "../parallel.hh"

namespace conformingsimplexgrid{
    
template< class Communicator, class Decomposition, class Grid >
class SpaceFillingCurve {
public:
        
    SpaceFillingCurve( Communicator &comm_in, Decomposition &decomp_in, Grid &grid_in )
    : comm(comm_in), decomp(decomp_in), grid(grid_in)
    {}
    
    template< class ElementVector >
    void init( ElementVector &indices_in){
        indices.resize(indices_in.size());
        for( size_t k = 0; k < indices_in.size(); k++)
            indices[k] = (double)indices_in[k];
    }
    
    template< class Refine_Data >
    void pre_refine( Refine_Data &data )
    {
        std::size_t nMarked = 0;
        for( size_t k = 0; k < data.element2newNode.size(); k++)
            if( data.element2newNode[k] > 0 )
                nMarked++;
        size_t counter = indices.size();    
        indices.resize( indices.size() + nMarked );
        for( size_t k = 0; k < data.element2newNode.size(); k++)
            if( data.element2newNode[k] > 0 ){
                indices[counter++] = indices[k] + pow( 2., -(grid.level[k]+1) );
            }
    }
    
    template< class Refine_Data >
    void refine( Refine_Data & )
    {
    }
    
    template< class Coarse_Data >
    void coarse( Coarse_Data &data ) 
    {
        for(size_t k = 0; k < data.left.size(); k++)
            if(indices[data.left[k]] > indices[data.right[k]]){
                indices[data.left[k]] = indices[data.right[k]];
                printf("hier\n");
            }
        remove_false(indices,data.remainingElements);
        
    }
    
    template< class Vector >
    void newIndices( Vector &, Vector &, Vector & ) {}

    template<class BoolVector, class Vector> 
    void preparePacking(BoolVector &remainingElements, Vector &faces2newFaces, Vector &edges2newEdges, Vector &nodes2newNodes)
    {
        std::size_t nElements = 0;
        for( size_t k = 0; k < remainingElements.size(); k++)
            if( remainingElements[k] )
                nElements++;
        indices_buffer.resize(remainingElements.size());
        std::size_t counter = 0;
        for( size_t k = 0; k < remainingElements.size(); k++)
            if( remainingElements[k] )
                indices_buffer[counter++] = indices[k];
        indices_buffer.resize(counter);
    }
    
    size_t size() {return indices.size();}
    
    void pack( void *buffer, std::size_t bufferSize, int &position ) {
        int message_size = (int)indices_buffer.size();
        comm.pack( &message_size, 1, buffer, bufferSize, position );
        comm.pack( indices_buffer.data(), indices_buffer.size(), buffer, bufferSize, position );
    }
    
    void unpack( void *buffer, std::size_t bufferSize, int &position ) {
        int message_size;
        comm.unpack( &message_size, 1, buffer, bufferSize, position );
        indices_buffer.resize( message_size );
        comm.unpack( indices_buffer.data(), indices_buffer.size(), buffer, bufferSize, position );
    }
        
    std::size_t bufferSize() {
        return comm.template bufferSize<int>(1)
             + comm.template bufferSize<double>(indices_buffer.size());     
    }

    void conjoin(){
        indices.conjoin(indices_buffer);
        indices_buffer.clear();
    }
    
    void sort( ManagedArray<int> &idx ) {
        sortByRow_void(indices,idx);
    }
    
    void test(){
        double min_value = indices.min();
        double max_value = indices.max();
        ManagedArray<double> min_values(comm.size());
        ManagedArray<double> max_values(comm.size());
        ManagedArray<double> min_values_(comm.size());
        ManagedArray<double> max_values_(comm.size());
        min_values.fill(min_value);
        max_values.fill(max_value);
        int flag1, flag2;
        flag1 = comm.all2all(min_values.data(),min_values_.data(),1);
        flag2 = comm.all2all(max_values.data(),max_values_.data(),1);
        printf("rank = %d min %lf, max %lf\n",comm.rank(),min_value,max_value);
        if(comm.master()){
            printf("min %lf, max %lf\n",min_value,max_value);
            for( int k = 0; k < comm.size(); k++)
                printf("| %lf | %lf |\n",min_values_[k],max_values_[k]);
            printf("\n");
        }
            
    }
    
//private:
    Communicator &comm;
    Decomposition &decomp;
    Grid &grid;
    ManagedArray<double> indices;
    ManagedArray<double> indices_buffer;
};
    
template< class Communicator, class Grid, class SpaceFillingCurve, class Partition >
void spaceFillingCurvePartition_master(Communicator &comm, Grid &grid, SpaceFillingCurve &curve, Partition &partition) 
{
    bool communicate = false;
    spaceFillingCurvePartition(comm, grid, curve, partition, communicate);
}

template< class Communicator, class Grid, class SpaceFillingCurve, class Partition >
void spaceFillingCurvePartition(Communicator &comm, Grid &grid, SpaceFillingCurve &curve, Partition &partition) 
{
    bool communicate = true;
    spaceFillingCurvePartition(comm, grid, curve, partition, communicate);
}

template< class Communicator, class Grid, class SpaceFillingCurve, class Partition >
void spaceFillingCurvePartition(Communicator &comm, Grid &grid, SpaceFillingCurve &curve, Partition &partition, bool communicate) 
{
    const int NumP = (int)comm.size();
    partition.resize(grid.nElements());
    partition.fill(comm.rank());
    ManagedArray<int> idx(grid.elements.size());
    curve.sort(idx);
    
    ManagedArray< int > elementsOnRank_(comm.size(),0);
    ManagedArray< int > elementsOnRank(comm.size(),0);
    elementsOnRank_.fill(grid.nElements());
    
    if(communicate)
        comm.all2all( elementsOnRank_.data(), elementsOnRank.data(), 1 );
    else
        elementsOnRank[comm.rank()] = grid.nElements();

    std::size_t nElementsGlobal = 0;
    for( std::size_t k = 0; k < elementsOnRank.size(); k++)
        nElementsGlobal += elementsOnRank[k];
    
    std::size_t nElementsPerRank = std::size_t((double)nElementsGlobal / (double)comm.size());
       
    ManagedArray<int> newElementsOnRank(comm.size(),nElementsPerRank);
    std::size_t totalElements = 0;
    for( std::size_t k = 0; k < newElementsOnRank.size(); k++)
        totalElements += newElementsOnRank[k];
    
    std::size_t iter1 = 0;
    while(totalElements > nElementsGlobal){
        newElementsOnRank[iter1++] -= 1;
        totalElements--;
        if( iter1 >= comm.size() ) iter1 = 0;
    }
    std::size_t iter2 = 0;
    while(totalElements < nElementsGlobal){
        newElementsOnRank[iter2++] += 1;
        totalElements++;
        if( iter2 >= comm.size() ) iter2 = 0;
    }
    
    std::size_t counter = 0;
    std::size_t rank_counter = 0;
    std::size_t counter1 = 0;
    std::size_t rank_counter1 = 0;
    std::size_t counter2 = 0;
    for( std::size_t el = 0; el < nElementsGlobal; el++ ){
        if( counter >= elementsOnRank[rank_counter] ){
            counter = 0;
            rank_counter++;
        }
        if( counter1 >= newElementsOnRank[rank_counter1] ){
            counter1 = 0;
            rank_counter1++;
        }
        if( rank_counter == comm.rank() && counter2 < partition.size() ){
            partition[idx[counter2]] = rank_counter1;
            counter2++;
        }
        counter++;
        counter1++;
    }
    
};

} // namespace conformingsimplexgrid

#endif // SPACE_FILLING_CURVE_HH
