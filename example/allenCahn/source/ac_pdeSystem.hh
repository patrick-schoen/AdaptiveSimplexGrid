#ifndef AC_PDE_SYSTEM_HH
#define AC_PDE_SYSTEM_HH

using namespace conformingsimplexgrid;

#include "../../../grid/mesh.hh"
#include "../../../grid/parallel.hh"

template< class Vector >
Vector f( Vector &u )
{
    ManagedArray<double,1> f_u( u.size() );
    for( int k = 0; k < u.size(); ++k)
        f_u[k] = 2.*((u[k]*u[k]*u[k]) - u[k]);
    return f_u;
}

template<int NV, int DIM, class Decomposition>
class SystemMatrix_AC
{
public:
    explicit SystemMatrix_AC( Decomposition &decomp, ManagedArray<int, NV> &n4e, ManagedArray<double, NV*DIM> &normals, 
                              ManagedArray<double, 1> &mass, ManagedArray<double, 1> &u, ManagedArray<double, 1> &u_old, 
                              double& eps, double &tau )
    : decomp_(decomp), n4e_( n4e ), normals_( normals ), mass_( mass), u_( u ), u_old_(u_old), eps_( eps), tau_( tau ) 
    {
        for( size_t j = 0; j < DIM + 1; j++)
            for( size_t k = 0; k < DIM + 1; k++){
                if( k == j){
                    Mloc[j*(DIM+1) + k] = 2. / ((DIM + 2.) * (DIM + 1.));		
                } else {
                    Mloc[j*(DIM+1) + k] = 1. / ((DIM + 2.) * (DIM + 1.));		
                }
        }

        for( int j = 0; j < DIM + 1; j++)
            for( int k = 0; k < DIM + 1; k++){
                if( k == j){
                    Mloc_lumped[j*(DIM+1) + k] = 1. / (DIM + 1.);		
                } else {
                    Mloc_lumped[j*(DIM+1) + k] = 0.;		
                }
        }
    }
        
    //Matrix-Vector Multiplication of System Matrix
    ManagedArray<double> operator * ( const ManagedArray<double, 1> &vec ) const
    {
        assert( vec.size() == u_old_.size() );
        size_t I; size_t J;
        double S, M, M_lumped, u_i, df, A_ij_T;
        ManagedArray<double, 1> res( vec.size() );
        res.fill(0.);
        for( size_t el = 0; el < n4e_.size(); el++ ){
            for( size_t i = 0; i < NV; i++ ){
                I = n4e_(el,i);
                u_i = u_old_[I];
                df = 2.*(3.*(u_i*u_i) - 1.);
                for( size_t j = 0; j < NV; j++ ){
                    J = n4e_(el,j);
                    S = 0.;
                    for(size_t p = 0; p < DIM; p++)
                        S += (1./(mass_[el]*DIM*DIM*NV))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+j*DIM+p];
                    M = mass_[el] * Mloc[i*NV+j];
                    M_lumped = mass_[el] * Mloc_lumped[i*NV+j];
                    A_ij_T = (1./tau_)*M + S + (1./eps_/eps_)*M_lumped*df;
                    res[I] += A_ij_T * vec[J]; 
                }
            }
        }
        return res;
    }
    
    //precompute bdy entries for parallel algorithms
    void get_bdy_entries( Decomposition &decomp )
    {
        isBoundaryElement.resize(n4e_.size());
        std::fill(isBoundaryElement.begin(), isBoundaryElement.end(), false);
        for( size_t el = 0; el < n4e_.size(); el++)
            for( size_t k = 0; k < NV; k++)
                if( decomp.isBoundaryNode( n4e_(el,k) ) )
                    isBoundaryElement[el] = true;
    }
        
    ManagedArray<double, 1> matrix_vector_bdy( const ManagedArray<double, 1> &vec)
    {
        assert( vec.size() == u_old_.size() );
        size_t I; size_t J;
        double S, M, M_lumped, u_i, df, A_ij_T;
        ManagedArray<double, 1> res( vec.size() );
        res.fill(0.);
        for( size_t el = 0; el < n4e_.size(); el++ ){
            if( !isBoundaryElement[el] )
                continue;
            for( size_t i = 0; i < NV; i++ ){
                I = n4e_(el,i); u_i = u_old_[I]; df = 2.*(3.*(u_i*u_i) - 1.);
                for( size_t j = 0; j < NV; j++ ){
                    J = n4e_(el,j);
                    S = 0.;
                    for(size_t p = 0; p < DIM; p++)
                        S += (1./(mass_[el]*DIM*DIM*NV))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+j*DIM+p];
                    M = mass_[el] * Mloc[i*(DIM+1)+j];
                    M_lumped = mass_[el] * Mloc_lumped[i*(DIM+1)+j];
                    A_ij_T = (1./tau_)*M + S + (1./eps_/eps_)*M_lumped*df;
                    res[I] += A_ij_T * vec[J]; 
                }
            }
        }
        return res; 
    }
    
    ManagedArray<double, 1> matrix_vector_inner( const ManagedArray<double, 1> &vec )
    {
        assert( vec.size() == u_old_.size() );
        size_t I; size_t J;
        double S, M, M_lumped, u_i, df, A_ij_T;
        ManagedArray<double, 1> res( vec.size() );
        res.fill(0.);
        for( size_t el = 0; el < n4e_.size(); el++){
            for( size_t i = 0; i < NV; i++){
                if( isBoundaryElement[el] )
                    continue;             
                I = n4e_(el,i); u_i = u_old_[I]; df = 2.*(3.*(u_i*u_i) - 1.);
                for( size_t j = 0; j < NV; j++){
                    J = n4e_(el,j);
                    S = 0.;
                    for(size_t p = 0; p < DIM; p++)
                        S += (1./(mass_[el]*DIM*DIM*NV))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+j*DIM+p];
                    M = mass_[el] * Mloc[i*(DIM+1)+j];
                    M_lumped = mass_[el] * Mloc_lumped[i*(DIM+1)+j];
                    A_ij_T = (1./tau_)*M + S + (1./eps_/eps_)*M_lumped*df;
                    res[I] += A_ij_T * vec[J];
                }
            }
        }
        return res; 
    }
    
    // Preconditiner
    ManagedArray<double, 1> precond( const ManagedArray<double, 1> &vec ) const
    {        
        size_t I;
        double S, M, M_lumped, u_i, df, A_ij_T;
        ManagedArray<double, 1> res( vec.size() );
        res.fill(0.);
        for( size_t el = 0; el < n4e_.size(); el++){
            for( size_t i = 0; i < NV; i++){
                I = n4e_(el,i); 			
                u_i = u_old_[I]; df = 2.*(3.*(u_i*u_i) - 1.);
                S = 0.;
                for(size_t p = 0; p < DIM; p++)
                    S += (1./(mass_[el]*DIM*DIM*NV))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+i*DIM+p];
                M = mass_[el] * Mloc[i*(DIM+1)+i];    
                M_lumped = mass_[el] * Mloc_lumped[i*(DIM+1)+i];   
                A_ij_T = (1./tau_)*M + S + (1./eps_/eps_)*M_lumped*df;
                res[I] += A_ij_T; 
            }
        }
        decomp_.sumNodeVector( res );
        return vec / res;
    }

    const ManagedArray<double, 1> uData () const { return u_; }
    
    ManagedArray<double, 1> uData () { return u_; }

private:
    Decomposition &decomp_;
    
    const ManagedArray<int, NV> &n4e_;
    const ManagedArray<double, NV*DIM> &normals_;
    const ManagedArray<double, 1> &mass_;
    const ManagedArray<double, 1> &u_;
    const ManagedArray<double, 1> &u_old_;
    
    const double &eps_;
    const double &tau_;

    std::vector<bool> isBoundaryElement;
    
    double Mloc[(DIM+1)*(DIM+1)];
    double Mloc_lumped[(DIM+1)*(DIM+1)];
};

template<int NV, int DIM, class Decomposition>
class EigenvalueMatrix_AC
{
public:
    explicit EigenvalueMatrix_AC( Decomposition &decomp, ManagedArray<int, NV> &n4e, ManagedArray<double, NV*DIM> &normals, 
                              ManagedArray<double, 1> &mass, ManagedArray<double, 1> &u, ManagedArray<double, 1> &u_old, 
                              double& eps, double &c_shift )
    : decomp_(decomp), n4e_( n4e ), normals_( normals ), mass_( mass), u_( u ), u_old_(u_old), eps_( eps), c_shift_( c_shift ) 
    {
        for( size_t j = 0; j < DIM + 1; j++)
            for( size_t k = 0; k < DIM + 1; k++){
                if( k == j){
                    Mloc[j*(DIM+1) + k] = 2. / ((DIM + 2.) * (DIM + 1.));               
                } else {
                    Mloc[j*(DIM+1) + k] = 1. / ((DIM + 2.) * (DIM + 1.));               
                }
        }

        for( int j = 0; j < DIM + 1; j++)
            for( int k = 0; k < DIM + 1; k++){
                if( k == j){
                    Mloc_lumped[j*(DIM+1) + k] = 1. / (DIM + 1.);               
                } else {
                    Mloc_lumped[j*(DIM+1) + k] = 0.;            
                }
        }
    }
        
    //Matrix-Vector Multiplication of System Matrix
    ManagedArray<double> operator * ( const ManagedArray<double, 1> &vec ) const
    {
        size_t I; size_t J;
        double S, M, M_lumped, u_i, df, A_ij_T;
        ManagedArray<double, 1> res( vec.size() );
        res.fill(0.);
        for( size_t el = 0; el < n4e_.size(); el++ ){
            for( size_t i = 0; i < NV; i++ ){
                I = n4e_(el,i);
                u_i = u_old_[I];
                df = 2.*(3.*(u_i*u_i) - 1.);
                for( size_t j = 0; j < NV; j++ ){
                    J = n4e_(el,j);
                    S = 0.;
                    for(size_t p = 0; p < DIM; p++)
                        S += (1./(mass_[el]*DIM*DIM*NV))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+j*DIM+p];
                    M = mass_[el] * Mloc[i*NV+j];
                    M_lumped = mass_[el] * Mloc_lumped[i*NV+j];
                    A_ij_T = S + (1./eps_/eps_)*M_lumped*(df) + (1./eps_/eps_)*M*(c_shift_);
                    res[I] += A_ij_T * vec[J]; 
                }
            }
        }
        return res;
    }
    
    //precompute bdy entries for parallel algorithms
    void get_bdy_entries( Decomposition &decomp )
    {
        isBoundaryElement.resize(n4e_.size());
        std::fill(isBoundaryElement.begin(), isBoundaryElement.end(), false);
        for( size_t el = 0; el < n4e_.size(); el++)
            for( size_t k = 0; k < NV; k++)
                if( decomp.isBoundaryNode( n4e_(el,k) ) )
                    isBoundaryElement[el] = true;
    }
        
    ManagedArray<double, 1> matrix_vector_bdy( const ManagedArray<double, 1> &vec)
    {
        size_t I; size_t J;
        double S, M, M_lumped, u_i, df, A_ij_T;
        ManagedArray<double, 1> res( vec.size() );
        res.fill(0.);
        for( size_t el = 0; el < n4e_.size(); el++ ){
            if( !isBoundaryElement[el] )
                continue;
            for( size_t i = 0; i < NV; i++ ){
                I = n4e_(el,i); u_i = u_old_[I]; df = 2.*(3.*(u_i*u_i) - 1.);
                for( size_t j = 0; j < NV; j++ ){
                    J = n4e_(el,j);
                    S = 0.;
                    for(size_t p = 0; p < DIM; p++)
                        S += (1./(mass_[el]*DIM*DIM*NV))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+j*DIM+p];
                    M = mass_[el] * Mloc[i*(DIM+1)+j];
                    M_lumped = mass_[el] * Mloc_lumped[i*(DIM+1)+j];
                    A_ij_T = S + (1./eps_/eps_)*M_lumped*(df) + (1./eps_/eps_)*M*(c_shift_);
                    res[I] += A_ij_T * vec[J]; 
                }
            }
        }
        return res; 
    }
    
    ManagedArray<double, 1> matrix_vector_inner( const ManagedArray<double, 1> &vec )
    {
        assert( vec.size() == u_old_.size() );
        size_t I; size_t J;
        double S, M, M_lumped, u_i, df, A_ij_T;
        ManagedArray<double, 1> res( vec.size() );
        res.fill(0.);
        for( size_t el = 0; el < n4e_.size(); el++){
            for( size_t i = 0; i < NV; i++){
                if( isBoundaryElement[el] )
                    continue;             
                I = n4e_(el,i); u_i = u_old_[I]; df = 2.*(3.*(u_i*u_i) - 1.);
                for( size_t j = 0; j < NV; j++){
                    J = n4e_(el,j);
                    S = 0.;
                    for(size_t p = 0; p < DIM; p++)
                        S += (1./(mass_[el]*DIM*DIM*NV))*normals_[el*NV*DIM+i*DIM+p]*normals_[el*NV*DIM+j*DIM+p];
                    M = mass_[el] * Mloc[i*(DIM+1)+j];
                    M_lumped = mass_[el] * Mloc_lumped[i*(DIM+1)+j];
                    A_ij_T = S + (1./eps_/eps_)*M_lumped*(df) + (1./eps_/eps_)*M*(c_shift_);
                    res[I] += A_ij_T * vec[J];
                }
            }
        }
        return res; 
    }
    
    // Preconditiner
    ManagedArray<double, 1> precond( const ManagedArray<double, 1> &vec ) const
    {        
        return vec;
    }

    const ManagedArray<double, 1> uData () const { return u_; }
    
    ManagedArray<double, 1> uData () { return u_; }

private:
    Decomposition &decomp_;
    
    const ManagedArray<int, NV> &n4e_;
    const ManagedArray<double, NV*DIM> &normals_;
    const ManagedArray<double, 1> &mass_;
    const ManagedArray<double, 1> &u_;
    const ManagedArray<double, 1> &u_old_;
    
    const double &eps_;
    const double &c_shift_;

    std::vector<bool> isBoundaryElement;
    
    double Mloc[(DIM+1)*(DIM+1)];
    double Mloc_lumped[(DIM+1)*(DIM+1)];
};

template<int NV, int DIM, class Decomposition>
class MassMatrix
{
public:
    explicit MassMatrix( Decomposition &decomp, ManagedArray<int, NV> &n4e, ManagedArray<double, 1> &mass )
    : decomp_(decomp), n4e_( n4e ),  mass_( mass)
    {
        for( size_t j = 0; j < DIM + 1; j++)
            for( size_t k = 0; k < DIM + 1; k++){
                if( k == j){
                    Mloc[j*(DIM+1) + k] = 2. / ((DIM + 2.) * (DIM + 1.));               
                } else {
                    Mloc[j*(DIM+1) + k] = 1. / ((DIM + 2.) * (DIM + 1.));               
                }
        }

        for( int j = 0; j < DIM + 1; j++)
            for( int k = 0; k < DIM + 1; k++){
                if( k == j){
                    Mloc_lumped[j*(DIM+1) + k] = 1. / (DIM + 1.);               
                } else {
                    Mloc_lumped[j*(DIM+1) + k] = 0.;            
                }
        }
    }
        
    //Matrix-Vector Multiplication of System Matrix
    ManagedArray<double> operator * ( const ManagedArray<double, 1> &vec ) const
    {
        size_t I; size_t J;
        double M, M_lumped, A_ij_T;
        ManagedArray<double, 1> res( vec.size() );
        res.fill(0.);
        for( size_t el = 0; el < n4e_.size(); el++ ){
            for( size_t i = 0; i < NV; i++ ){
                I = n4e_(el,i);
                for( size_t j = 0; j < NV; j++ ){
                    J = n4e_(el,j);
                    M = mass_[el] * Mloc[i*NV+j];
                    M_lumped = mass_[el] * Mloc_lumped[i*NV+j];
                    A_ij_T = M;
                    res[I] += A_ij_T * vec[J]; 
                }
            }
        }
        return res;
    }
private:
    Decomposition &decomp_;
    
    const ManagedArray<int, NV> &n4e_;
    const ManagedArray<double, 1> &mass_;
    
    double Mloc[(DIM+1)*(DIM+1)];
    double Mloc_lumped[(DIM+1)*(DIM+1)];
};



template<int NV, int DIM>
class rhs_AC
{
public:
    explicit rhs_AC( ManagedArray<int, NV> &n4e, ManagedArray<double, NV*DIM> &normals, ManagedArray<double, 1> &mass, ManagedArray<double,1> &u, double& eps, double &tau )
    : n4e_( n4e ), normals_( normals ), mass_( mass), u_( u ), eps_(eps), tau_(tau)
    {
        for( size_t j = 0; j < DIM + 1; j++)
            for( size_t k = 0; k < DIM + 1; k++){
                if( k == j){
                    Mloc[j*(DIM+1) + k] = 2. / ((DIM + 2.) * (DIM + 1.));		
                } else {
                    Mloc[j*(DIM+1) + k] = 1. / ((DIM + 2.) * (DIM + 1.));		
                }
        }

        for( int j = 0; j < DIM + 1; j++)
            for( int k = 0; k < DIM + 1; k++){
                if( k == j){
                    Mloc_lumped[j*(DIM+1) + k] = 1. / (DIM + 1.);		
                } else {
                    Mloc_lumped[j*(DIM+1) + k] = 0.;		
                }
        }
    }
    
    void get( ManagedArray<double, 1> &val ) {
        int I,J;
        double M, M_lumped, f, df;
        val.resize(u_.size());
        val.fill(0.);
        for(size_t j = 0; j < n4e_.size(); j++){
            for(size_t k = 0; k < NV; k++){   
                I = n4e_[j*NV + k];  
                f = 2.*(u_[I]*u_[I]*u_[I] - u_[I]); df = 2.*(3.*(u_[I]*u_[I]) - 1.);	
                for(size_t m = 0; m < NV; m++){
                    J = n4e_[j*NV + m];
                    M = mass_[j] * Mloc[k*(DIM+1) + m];    
                    M_lumped = mass_[j] * Mloc_lumped[k*(DIM+1) + m];
                    val[J] += (1./tau_)*M*u_[I] +(1./eps_/eps_)*M_lumped*(df*u_[I] - f);
                }
            }
        }
    }

private:
    const ManagedArray<int, NV> &n4e_;
    const ManagedArray<double, NV*DIM> &normals_;
    const ManagedArray<double, 1> &mass_;
    const ManagedArray<double, 1> &u_;
    const double &eps_;
    const double &tau_;
    
    double Mloc[(DIM+1)*(DIM+1)];
    double Mloc_lumped[(DIM+1)*(DIM+1)];
};

#endif // AC_PDE_SYSTEM_HH
