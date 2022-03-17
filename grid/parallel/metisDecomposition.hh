#ifndef METIS_DECOMPOSITION_HH
#define METIS_DECOMPOSITION_HH

#include <metis.h>
#include <parmetis.h>

#include "../simplexGrid.hh"
#include "../parallel.hh"

namespace conformingsimplexgrid{
 
template< class Grid, class Partition >
void metisDecomposition (Grid &grid, Partition &partition, int NumP) 
{
        ManagedArray<int> vtxdist( 2 );
        vtxdist[0] = 0;
        vtxdist[1] = grid.nElements();
    
        if( NumP <= 1){
            partition.resize( grid.nElements() );
            partition.fill(0);
            return;
        }
        const int NV = Grid::NV;
        ManagedArray<int,NV> neighbors;
        getNeighbors ( grid, neighbors );
        
        partition.resize( grid.nElements() );
        
        ManagedArray< idx_t > adjncy( grid.nElements() * NV);
        ManagedArray< idx_t > xadj( grid.nElements() + 1);

        std::vector< int > buffer( NV );
    
        int count_xadj = 0;
        int count = 0;
        xadj[0] = 0;                                                //Assumption: jedes Element besitzt mindestens einen Nachbarn
        for( int el = 0; el < (int)grid.nElements(); el++){
            count = 0;
            std::fill( buffer.begin(), buffer.end(), 0 );
            for( int k = 0; k < NV; k++){
                if( neighbors(el,k) < 0 )
                    continue;
                buffer[count] = neighbors(el,k);
                count++;
            }
            std::sort( buffer.begin(), buffer.begin() + count );
            std::copy( buffer.begin(), buffer.begin() + count, adjncy.array_.begin() + count_xadj );
            count_xadj += count;
            xadj[el+1] = count_xadj;
        }
        adjncy.resize(count_xadj);
        idx_t ncon = 1;
        idx_t nparts = NumP;
        idx_t wgtflag = 0;
        idx_t numflag = 0;
        ManagedArray< real_t > tpwgts( ncon*nparts, 1. / (double)nparts ); 
        ManagedArray< real_t > ubvec( ncon, 1.05 );
        ManagedArray< idx_t > options(3, 0);
        
        idx_t nvtxs = grid.nElements();
        idx_t objval = 0;
        idx_t edgecut;
        
        //METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(), 0, 0, 0, &nparts, 0, 0, 0, &objval, partition.data() );
        MPI_Comm communicator = MPI_COMM_SELF;
        
        ParMETIS_V3_PartKway( vtxdist.data(), xadj.data(), adjncy.data(), 
                              0, 0, &wgtflag, &numflag, &ncon, &nparts, tpwgts.data(), ubvec.data(), options.data(),
                              &edgecut, partition.data(), &communicator);

};

template< class Communicator, class Decomposition, class Grid, class Partition >
void parMetisDecomposition (Communicator &comm, Decomposition &decomposition, Grid &grid, Partition &partition) 
{
        const int NumP = (int)comm.size();
        ManagedArray<int> vtxdist( comm.size() + 1);
        ManagedArray<int> nElementsPerRank_send( comm.size(), grid.nElements() );
        ManagedArray<int> nElementsPerRank_recv( comm.size(), 0 );
        
        comm.all2all( nElementsPerRank_send.data(), nElementsPerRank_recv.data(), 1 );
        vtxdist[0] = 0;
        for(size_t p = 0, ctr = 0; p < comm.size(); p++)
            vtxdist[p+1] = ( ctr += nElementsPerRank_recv[p] );
        
    
        if( NumP <= 1){
            partition.resize( grid.nElements() );
            partition.fill(0);
            return;
        }
        const int NV = Grid::NV;
        ManagedArray<int,NV> neighbors;       
        getGlobalNeighbors( comm, decomposition, vtxdist, grid,neighbors );
        
        partition.resize( grid.nElements() );
        
        ManagedArray< idx_t > adjncy( grid.nElements() * NV);
        ManagedArray< idx_t > xadj( grid.nElements() + 1);

        std::vector< int > buffer( NV );
    
        int count_xadj = 0;
        int count = 0;
        xadj[0] = 0;                                                //Assumption: jedes Element besitzt mindestens einen Nachbarn
        for( int el = 0; el < (int)grid.nElements(); el++){
            count = 0;
            std::fill( buffer.begin(), buffer.end(), 0 );
            for( int k = 0; k < NV; k++){
                if( neighbors(el,k) < 0 )
                    continue;
                buffer[count] = neighbors(el,k);
                count++;
            }
            std::sort( buffer.begin(), buffer.begin() + count );
            std::copy( buffer.begin(), buffer.begin() + count, adjncy.array_.begin() + count_xadj );
            count_xadj += count;
            xadj[el+1] = count_xadj;
        }
        adjncy.resize(count_xadj);

         
        idx_t ncon = 1;
        idx_t nparts = NumP;
        
        idx_t nvtxs = grid.nElements();
        idx_t objval = 0;
        
        idx_t wgtflag = 0;
        idx_t numflag = 0;
        ManagedArray< real_t > tpwgts( ncon*nparts, 1. / (double)nparts ); 
        ManagedArray< real_t > ubvec( ncon, 1.01 );
        real_t itr = 100000.;
        ManagedArray< idx_t > options(4, 0);
        
        
        idx_t edgecut;
        MPI_Comm communicator = MPI_COMM_WORLD;
        
        ParMETIS_V3_AdaptiveRepart( vtxdist.data(), xadj.data(), adjncy.data(), 
                                    0, 0, 0, &wgtflag, &numflag, &ncon, &nparts, tpwgts.data(), ubvec.data(), &itr, options.data(),
                                    &edgecut, partition.data(), &communicator);

};

} // namespace conformingsimplexgrid

#endif // METIS_DECOMPOSITION_HH
