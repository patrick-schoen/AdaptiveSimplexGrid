
#ifndef POISSON_PDE_SYSTEM_HH
#define POISSON_PDE_SYSTEM_HH

#include "../../../grid/mesh.hh"
#include "../../../grid/parallel.hh"

#include <math.h>

using namespace conformingsimplexgrid;

template< class Grid >
class DirichletNodes
{
public:
    explicit DirichletNodes( Grid &grid )
    : grid_(grid)
    {}
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    void get( ManagedArray<int, 1> &val ){      
        ManagedArray<int> diriNodes(grid_.nNodes(),0);
        ManagedArray<int> diriEdges(grid_.nEdges(),0);
        
        if( NV == 3 ) {
            for( int el = 0; el < (int)grid_.elements.size(); el++)
                for( int k = 0; k < 3; k++)
                    diriEdges[grid_.element2edges(el,k)]++;
        } else {
            ManagedArray<int> diriFaces(grid_.nFaces(),0);
            for( int el = 0; el < (int)grid_.elements.size(); el++)
                for( int k = 0; k < 4; k++)
                    diriFaces[grid_.element2faces(el,k)]++;
            
            for( int el = 0; el < (int)diriFaces.size(); el++)
                if( diriFaces[el] == 1 )
                    for( int k = 0; k < 3; k++)
                        diriEdges[grid_.face2edges(el,k)] = 1;
        }
        
        for( int el = 0; el < (int)diriEdges.size(); el++)
            if( diriEdges[el] == 1 ){
                diriNodes[ grid_.edges(el,0) ] = 1;
                diriNodes[ grid_.edges(el,1) ] = 1;
            }
            
        int counter = 0;
        val.resize(diriNodes.size());
        for( int k = 0; k < (int)diriNodes.size(); k++)
            if(diriNodes[k] > 0 )
                val[counter++] = k;
        
        val.resize(counter);
    }
 
private:
    Grid &grid_;
};

template<class Grid>
class SystemMatrix
{
public:
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    explicit SystemMatrix( Grid &grid, ManagedArray<double, NV*DIM> &normals, ManagedArray<double, 1> &mass)
    : grid_(grid), normals_( normals ), mass_( mass) 
    {         
        DirichletNodes<Grid> DN( grid_ );
        DN.get(diriNodes);
    }
        
    //Matrix-Vector Multiplication of System Matrix
    ManagedArray<double> operator * ( const ManagedArray<double, 1> &vec ) const
    {   
        ManagedArray<double, 1> res( vec.size() , 0.);
        for( size_t el = 0; el < grid_.elements.size(); el++ ){
            for( size_t i = 0; i < NV; i++ ){
                size_t I = grid_.elements(el,i);
                for( size_t j = 0; j < NV; j++ ){
                    size_t J = grid_.elements(el,j);
                    double S = 0.;
                    for(size_t p = 0; p < DIM; p++)
                        S += (1./(mass_[el]*DIM*DIM))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+j*DIM+p];
                    res[I] += S * vec[J];
                }
            }
        }
        
        for(size_t k = 0; k < diriNodes.size(); k++)
            res[diriNodes[k]] = vec[diriNodes[k]];

        return res;
    }
    
    // Jacobian Preconditioner
    ManagedArray<double, 1> precond( const ManagedArray<double, 1> &vec ) const
    {        
        size_t I;
        double S;
        ManagedArray<double, 1> res( vec.size() );
        res.fill(0.);
        for( size_t el = 0; el < grid_.elements.size(); el++){
            for( size_t i = 0; i < NV; i++){
                I = grid_.elements(el,i);
                S = 0.;
                for(size_t p = 0; p < DIM; p++)
                    S += (1./(mass_[el]*DIM*DIM*NV))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+i*DIM+p];
                res[I] += S; 
            }
        }
        return vec / res;
    }

private:
    Grid &grid_;
    const ManagedArray<double, NV*DIM> &normals_;
    const ManagedArray<double, 1> &mass_;
    ManagedArray<int> diriNodes;
};


template<class Grid>
class Rhs
{
public:
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    explicit Rhs( Grid &grid, ManagedArray<double, NV*DIM> &normals, ManagedArray<double, 1> &mass )
    : grid_(grid), normals_( normals ), mass_( mass)
    {}
    
    void get( ManagedArray<double, 1> &val ) {
        int I;
        val.resize(grid_.coordinates.size());
        val.fill(0.);
        
        ManagedArray<double,DIM> mp_T(grid_.elements.size(), 0.);
        for( int el = 0; el < (int)grid_.elements.size(); el++)
            for( int k = 0; k < NV; k++)
                for( int d = 0; d < DIM; d++)
                    mp_T(el,d) += (1./double(NV))*( grid_.coordinates(grid_.elements(el,k),d) );
                    
        for(size_t j = 0; j < grid_.elements.size(); j++){
            for(size_t k = 0; k < NV; k++){   
                I = grid_.elements[j*NV + k];
                double tmp[NV-1];
                for( size_t r = 0; r < NV-1; r++)
                    tmp[r] = mp_T(j,r);
                val[I] += (1./double(NV)) * mass_[j] * f(tmp);
            }
        }
        
        DirichletNodes<Grid> DN( grid_ );
        ManagedArray<int> diriNodes;
        DN.get(diriNodes);
        
        for( size_t k = 0; k < diriNodes.size(); k++){
            double tmp[NV-1];
            for( size_t r = 0; r < NV-1; r++)
                tmp[r] = grid_.coordinates(diriNodes[k],r);
            val[diriNodes[k]] = u_D(tmp);
        }
    }
    
    double f( double x[] ) {
        //return (8.*M_PI*M_PI) * sin(2.*M_PI*x[0])*sin(2.*M_PI*x[1]);
        return 1.;
    }
    
    double u_D( double x[] ){
        return 0.;
    }
    

private:
    Grid &grid_;
    const ManagedArray<double, NV*DIM> &normals_;
    const ManagedArray<double, 1> &mass_;
};


template< class Grid, class Decomposition >
class DirichletNodes_parallel
{
public:
    explicit DirichletNodes_parallel( Grid &grid, Decomposition &decomp)
    : grid_(grid), decomp_(decomp)
    {}
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    void get( ManagedArray<int, 1> &val ){      
        ManagedArray<int> diriNodes(grid_.nNodes(),0);
        ManagedArray<int> diriEdges(grid_.nEdges(),0);
        
        for( int el = 0; el < (int)grid_.elements.size(); el++){
            diriEdges[grid_.element2edges(el,0)]++;
            diriEdges[grid_.element2edges(el,1)]++;
            diriEdges[grid_.element2edges(el,2)]++;
        }
            
        for( int p = 0; p < (int)decomp_.size(); p++){
            if ( p == decomp_.rank() )
                continue;
            for( int k = 0; k < (int)decomp_.rank2bdyEdges[p].size(); k++){
                diriEdges[decomp_.rank2bdyEdges[p][k]]++;
            }
        }
        
        for( int el = 0; el < (int)diriEdges.size(); el++)
            if( diriEdges[el] == 1 ){
                diriNodes[ grid_.edges(el,0) ] = 1;
                diriNodes[ grid_.edges(el,1) ] = 1;
            }
            
        int counter = 0;
        val.resize(diriNodes.size());
        for( int k = 0; k < (int)diriNodes.size(); k++)
            if(diriNodes[k] > 0 )
                val[counter++] = k;
        
        val.resize(counter);
    }
 
private:
    Grid &grid_;
    Decomposition &decomp_;
};



template<class Grid, class Decomposition>
class SystemMatrix_parallel
{
public:
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    explicit SystemMatrix_parallel( Grid &grid, Decomposition &decomp, ManagedArray<double, NV*DIM> &normals, ManagedArray<double, 1> &mass)
    : grid_(grid), decomp_(decomp), normals_( normals ), mass_( mass) 
    {         
        DirichletNodes_parallel<Grid, Decomposition > DN( grid_, decomp_ );
        DN.get(diriNodes);
    }
        
    //Matrix-Vector Multiplication of System Matrix
    ManagedArray<double> operator * ( const ManagedArray<double, 1> &vec ) const
    {   
        ManagedArray<double, 1> res( vec.size() , 0.);
        for( size_t el = 0; el < grid_.elements.size(); el++ ){
            for( size_t i = 0; i < NV; i++ ){
                size_t I = grid_.elements(el,i);
                for( size_t j = 0; j < NV; j++ ){
                    size_t J = grid_.elements(el,j);
                    double S = 0.;
                    for(size_t p = 0; p < DIM; p++)
                        S += (1./(mass_[el]*DIM*DIM))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+j*DIM+p];
                    res[I] += S * vec[J];
                }
            }
        }
        
        for(size_t k = 0; k < diriNodes.size(); k++)
            res[diriNodes[k]] = vec[diriNodes[k]];

        return res;
    }
    
    // Jacobian Preconditioner
    ManagedArray<double, 1> precond( const ManagedArray<double, 1> &vec ) const
    {        
        size_t I;
        double S;
        ManagedArray<double, 1> res( vec.size() );
        res.fill(0.);
        for( size_t el = 0; el < grid_.elements.size(); el++){
            for( size_t i = 0; i < NV; i++){
                I = grid_.elements(el,i);
                S = 0.;
                for(size_t p = 0; p < DIM; p++)
                    S += (1./(mass_[el]*DIM*DIM*NV))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+i*DIM+p];
                res[I] += S; 
            }
        }
        decomp_.sumNodeVector( res );
        return vec / res;
    }

private:
    Grid &grid_;
    Decomposition &decomp_;
    const ManagedArray<double, NV*DIM> &normals_;
    const ManagedArray<double, 1> &mass_;
    ManagedArray<int> diriNodes;
};



template<class Grid, class Decomposition>
class Rhs_parallel
{
public:
    static const int NV = Grid::NV;
    static const int DIM = Grid::dim;
    
    explicit Rhs_parallel( Grid &grid, Decomposition &decomp, ManagedArray<double, NV*DIM> &normals, ManagedArray<double, 1> &mass )
    : grid_(grid), decomp_(decomp), normals_( normals ), mass_( mass)
    {}
    
    void get( ManagedArray<double, 1> &val ) {
        int I;
        val.resize(grid_.coordinates.size());
        val.fill(0.);
        
        ManagedArray<double,2> mp_T(grid_.elements.size(), 0.);
        for( int el = 0; el < (int)grid_.elements.size(); el++)
            for( int k = 0; k < NV; k++)
                for( int d = 0; d < DIM; d++)
                    mp_T(el,d) += (1./double(NV))*( grid_.coordinates(grid_.elements(el,k),d) );
   
        for(size_t j = 0; j < grid_.elements.size(); j++){
            for(size_t k = 0; k < NV; k++){   
                I = grid_.elements[j*NV + k];
                double tmp[NV-1];
                for( size_t r = 0; r < NV-1; r++)
                    tmp[r] = mp_T(j,r);
                val[I] += (1./double(NV)) * mass_[j] * f(tmp);
            }
        }
        
        DirichletNodes_parallel<Grid, Decomposition > DN( grid_, decomp_ );
        ManagedArray<int> diriNodes;
        DN.get(diriNodes);
        
        for( size_t k = 0; k < diriNodes.size(); k++){
            double tmp[NV-1];
            for( size_t r = 0; r < NV-1; r++)
                tmp[r] = grid_.coordinates(diriNodes[k],r);
            val[diriNodes[k]] = u_D(tmp);
        }
        
        decomp_.sumNodeVector(val);
    }
    
    double f( double x[] ) {
        return 1.;
    }
    
    double u_D( double x[] ){
        return 0.;
    }
    

private:
    Grid &grid_;
    Decomposition &decomp_;
    const ManagedArray<double, NV*DIM> &normals_;
    const ManagedArray<double, 1> &mass_;
};

#endif // POISSON_PDE_SYSTEM_HH
