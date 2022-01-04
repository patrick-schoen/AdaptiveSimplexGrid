#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "algos.hh"

// scaled outer nomal vectors, consider unoriented(!) n-Simplex T = (x_0, ... x_n)
// mass = | (1/n!)*det( x_1 - x_0, ... , x_n - x_0 ) | = | (1/n!)*<x_1 - x_0, n_0>|
// sign = signum( <x_1 - x_0, n_0> )
// n_0_i = sign * det( e_i, x_2 - x_1, ..., x_n - x_1 )
// n_k_i = sign * det( e_i, x_1 - x_0, ..., x_{k-1} - x_0, x_{k+1} - x_0, ..., x_n - x_0), for k = 1,..n
// it holds that |n_k| = (n-1)! * |E_k|, for k = 0,...,n and E_k beegin the Hypersurface E_k = {x_0, ... , x_{k-1}, x_{k+1}, ... , x_n}
namespace conformingsimplexgrid { 

template< class Mass, class Normals, class Elements, class Coord>
void outerNormals( Mass &mass, Normals &normals, Elements &elements, Coord &coordinates) {
    static const int DIM = Coord::DIM;
    static const int NV = Elements::DIM;

    mass.resize( elements.size());
    normals.resize( elements.size());
    Faculty<DIM> faculty;
    double fac = (double)(1./faculty.value);    
    if( NV == 3 && DIM == 2) {
        double x_0_0, x_0_1, x_1_0, x_1_1, x_2_0, x_2_1;
        for( size_t j = 0; j < elements.size(); j++){
            x_0_0 = coordinates( elements(j,0), 0); x_0_1 = coordinates( elements(j,0), 1);
            x_1_0 = coordinates( elements(j,1), 0); x_1_1 = coordinates( elements(j,1), 1);
            x_2_0 = coordinates( elements(j,2), 0); x_2_1 = coordinates( elements(j,2), 1);
            //n_0_i = det(e_i,x_2 - x_1)
            normals(j,0) = (x_2_1 - x_1_1);
            normals(j,1) = -(x_2_0 - x_1_0);
            //n_1_i = det(e_i,x_2 - x_0)
            normals(j,2) = -(x_2_1 - x_0_1);
            normals(j,3) = (x_2_0 - x_0_0);
            //n_2_i = det(e_i,x_1 - x_0)
            normals(j,4) = (x_1_1 - x_0_1);
            normals(j,5) = -(x_1_0 - x_0_0);
            //mass_signed = <x_0 - x_1, n_0>
            mass[j] = (1./2.)*((x_1_0 - x_0_0)*normals(j,0)+(x_1_1 - x_0_1)*normals(j,1));
            if(mass[j] < 0){
                mass[j] *= -1.;
                for(size_t i = 0; i < NV*DIM; i++) normals(j,i) *= -1.;
            }
        }
    } else if( NV == 4 && DIM == 3) {
        double x_0_0, x_0_1, x_0_2, x_1_0, x_1_1, x_1_2, x_2_0, x_2_1, x_2_2, x_3_0, x_3_1, x_3_2;
        for( size_t j = 0; j < elements.size(); j++){
            x_0_0 = coordinates( elements(j,0), 0); x_0_1 = coordinates( elements(j,0), 1); x_0_2 = coordinates( elements(j,0), 2);
            x_1_0 = coordinates( elements(j,1), 0); x_1_1 = coordinates( elements(j,1), 1); x_1_2 = coordinates( elements(j,1), 2);
            x_2_0 = coordinates( elements(j,2), 0); x_2_1 = coordinates( elements(j,2), 1); x_2_2 = coordinates( elements(j,2), 2);
            x_3_0 = coordinates( elements(j,3), 0); x_3_1 = coordinates( elements(j,3), 1); x_3_2 = coordinates( elements(j,3), 2);
            //n_0_i = det(e_i,x_2 - x_1,x_3 - x_1)
            normals(j,0) = -(x_2_1-x_1_1)*(x_3_2-x_1_2)+(x_2_2-x_1_2)*(x_3_1-x_1_1);
            normals(j,1) = -(x_2_2-x_1_2)*(x_3_0-x_1_0)+(x_2_0-x_1_0)*(x_3_2-x_1_2);		
            normals(j,2) = -(x_2_0-x_1_0)*(x_3_1-x_1_1)+(x_2_1-x_1_1)*(x_3_0-x_1_0);
            //n_1_i = det(e_i,x_2 - x_0,x_3 - x_0)
            normals(j,3) = (x_2_1-x_0_1)*(x_3_2-x_0_2)-(x_2_2-x_0_2)*(x_3_1-x_0_1);
            normals(j,4) = (x_2_2-x_0_2)*(x_3_0-x_0_0)-(x_2_0-x_0_0)*(x_3_2-x_0_2);		
            normals(j,5) = (x_2_0-x_0_0)*(x_3_1-x_0_1)-(x_2_1-x_0_1)*(x_3_0-x_0_0);
            //n_2_i = det(e_i,x_1 - x_0,x_3 - x_0)
            normals(j,6) = (x_1_2-x_0_2)*(x_3_1-x_0_1)-(x_1_1-x_0_1)*(x_3_2-x_0_2);
            normals(j,7) = (x_1_0-x_0_0)*(x_3_2-x_0_2)-(x_1_2-x_0_2)*(x_3_0-x_0_0);		
            normals(j,8) = (x_1_1-x_0_1)*(x_3_0-x_0_0)-(x_1_0-x_0_0)*(x_3_1-x_0_1);
            //n_3_i = det(e_i,x_1 - x_0,x_2 - x_0)
            normals(j,9) = (x_1_1-x_0_1)*(x_2_2-x_0_2)-(x_1_2-x_0_2)*(x_2_1-x_0_1);
            normals(j,10) = (x_1_2-x_0_2)*(x_2_0-x_0_0)-(x_1_0-x_0_0)*(x_2_2-x_0_2);		
            normals(j,11) = (x_1_0-x_0_0)*(x_2_1-x_0_1)-(x_1_1-x_0_1)*(x_2_0-x_0_0);
            // mass_signed = (1/n!)*<x_1 - x_0, n_0>
            mass[j]= -fac*((x_1_0-x_0_0)*normals(j,0)+(x_1_1-x_0_1)*normals(j,1)+(x_1_2-x_0_2)*normals(j,2));	//TODO: check the signs
            if(mass[j] < 0){
                mass[j] *= -1.;
                for(size_t i = 0; i < NV*DIM; i++) normals(j,i) *= -1.;
            }
        }
    } else {
            assert( "OuterNormals not implemented for Number of Vertices or Dimension");
    }
}

template< class Grid, class Decomposition >
class local_norm
/// This class allows the evaluation of local norms, like ||v||_L^2{T).
{
    static const size_t NV = Grid::NV;
    static const size_t DIM = Grid::dim;
    
public:
    explicit local_norm ( Grid &mesh, Decomposition &decomp, ManagedArray<double, NV*DIM> &normals, ManagedArray<double, 1> &mass)
    : mesh_( mesh ), decomp_(decomp), normals_( normals ), mass_( mass)
    {
        for( size_t j = 0; j < NV; j++)
            for( size_t k = 0; k < NV; k++){
                if( k == j){
                    Mloc[j*NV + k] = 2. / ((NV + 1.) * NV);
                } else {
                    Mloc[j*NV + k] = 1. / ((NV + 1.) * NV);
                }
        }

        for( size_t j = 0; j < NV; j++)
            for( size_t k = 0; k < NV; k++){
                if( k == j){
                    Mloc_lumped[j*NV + k] = 1. / NV;		
                } else {
                    Mloc_lumped[j*NV + k] = 0.;		
                }
        }	

        Faculty<DIM> faculty;
        fac = (double)(1./faculty.value);

    }
    
    // || xi ||_T^2
    ManagedArray<double, 1> l2( const ManagedArray<double, 1> &u ) 
    ///returns a vector of local L^2 norms, ||v||_L^2{T).
    {
        ManagedArray<double> res( mesh_.elements.size(), 0. );
        for( size_t el = 0; el < mesh_.elements.size(); el++ )
            for( size_t i = 0; i < NV; i++ ){
                const size_t I = mesh_.elements(el,i);
                for( size_t j = 0; j < NV; j++ ){
                    const size_t J = mesh_.elements(el,j);
                    const double M = mass_[el] * Mloc[i*NV + j];
                    res[el] += u[I] * u[J] * M; 
                }
            }
        return res;
    }

    // || D(xi) ||_T^2			
    ManagedArray<double, 1> h1( const ManagedArray<double, 1> &u )
     ///returns a vector of local H^1 norms, ||Dv||_L^2{T).
    {
        size_t I; size_t J;
        double S;
        ManagedArray<double, 1> res( mesh_.elements.size() );
        res.fill(0.);
        for( size_t el = 0; el < mesh_.elements.size(); el++){
            for( size_t i = 0; i < NV; i++){
                I = mesh_.elements(el,i);
                    for( size_t j = 0; j < NV; j++){
                        J = mesh_.elements(el,j);					
                        S = 0.;
                        for(size_t p = 0; p < DIM; p++) 
                            S += (1./(mass_[el]*DIM*DIM*NV))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+j*DIM+p];
                            res[el] += u[I] * u[J] * S; 
                    }
            }
        }
        return res;
    }

    // sum_{E \in T} h_E * ||[D(u)*n_E]||_L^2(E) = sum_{E \in T} |[D(u)*n_E*h_E]|^2	
    ManagedArray<double, 1> h1_jump( const ManagedArray<double, 1> &u )
    ///returns a vector of local jumps, scaled by h_E, h_E*||[Dv*nE]||_L^2{T).
    {
            
        ManagedArray<double, 1> res( mesh_.nElements(), 0. );
        ManagedArray<double, 1 > Du_nE( mesh_.nFaces(), 0. );

        for( size_t el = 0; el < mesh_.elements.size(); el ++ )
            for (size_t k = 0; k < NV; k++){
                const int face = mesh_.element2faces(el,(NV-1)-k);
                double temp = 0;
                for (size_t j = 0; j < NV; j++){
                    for( size_t p = 0; p < DIM; p++){
                        temp += u[mesh_.elements(el,j)]*normals_[el*NV*DIM+j*DIM+p]*normals_[el*NV*DIM+k*DIM+p];
                    }
                }
                Du_nE[face] += (1/mass_[el]) * fac * temp;
            }
        
        decomp_.sumFaceVector( Du_nE );
            
        res.fill(0.);
        for( size_t el = 0; el < mesh_.elements.size(); el ++ )
            for (size_t k = 0; k < NV; k++){
                const int face = mesh_.element2faces(el,k);
                res[el] += Du_nE[face]*Du_nE[face];
            }
        return res;
    }
    
    template< class uVector, class Vector >
    ManagedArray<double,1> coarse_norm( Grid &grid, uVector &u, Vector &marked ){
        const int NV = Grid::NV;
        
        ManagedArray<int> left, right;
        ManagedArray<int> markedNodes(grid.nNodes(), 0);
        
        findRemovableNodes( grid, marked, markedNodes);
        
        getBrothers( grid, marked, left, right );
        
        ManagedArray<double,1> res(grid.nElements(), 0.);
        if( left.size() == 0)
            return res;
        
        int type = (grid.level[left[0]] % (NV-1));
        for(int k = 0; k < (int)left.size(); k++){
            int z0 = grid.elements(left[k],0);
            int z1 = grid.elements(right[k],NV-1);
            if( type == 0 )
                int z1 = grid.elements(right[k],0);
            int y = grid.elements(left[k],1);
            res[left[k]] = mass_[left[k]] * Mloc[0] * (u[y] - .5*(u[z0] + u[z1]))*(u[y] - .5*(u[z0] + u[z1]));
            res[right[k]] = mass_[right[k]] * Mloc[0] * (u[y] - .5*(u[z0] + u[z1]))*(u[y] - .5*(u[z0] + u[z1]));
        }
        return res;
    }

    ManagedArray<double, 1> h_T_2 ( ){		
        ManagedArray<double, 1> res( mesh_.elements.size() );
        res.fill(0);			
        for( size_t el = 0; el < mesh_.elements.size(); el++){
            for( size_t k = 0; k < DIM; k++ ){
                res[el] += (mesh_.coordinates(mesh_.elements(el,0),k)-mesh_.coordinates(mesh_.elements(el,NV-1),k))
                         * (mesh_.coordinates(mesh_.elements(el,0),k)-mesh_.coordinates(mesh_.elements(el,NV-1),k));	
            }
        }		
        return res;
    }

    int getFaceNumberInNeighbor( int nb, int element ){
    for( size_t k = 0; k < NV; k++ )
        if( mesh_.neigh(nb,k) == element )
        return k;
    printf("Error: Element not found in Neighbor\n");
    abort();
    return 0;            
}

private:
    Grid &mesh_;
    Decomposition &decomp_;
    const ManagedArray<double, NV*DIM> &normals_;
    const ManagedArray<double, 1> &mass_;

    double fac;
    double Mloc[(DIM+1)*(DIM+1)];
    double Mloc_lumped[(DIM+1)*(DIM+1)];
};

} // namespace conformingsimplexgrid

#endif /*GEOMETRY_H_*/
