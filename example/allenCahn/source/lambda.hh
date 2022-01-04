#ifndef LAMBDA_HH
#define LAMBDA_HH

#include<math.h>

namespace conformingsimplexgrid {
    template< int NV, int DIM>
    double calcEigenValue( MeshType &mesh, ManagedArray<double, NV*DIM> &normals, ManagedArray<double, 1> &mass, ManagedArray<double,1> &u, ManagedArray<double, 1> &u_old, double eps, double tau, ManagedArray<double, 1> &w)
    {
        MassMatrix_AC<NV, DIM > M( mesh.n4e, normals, mass, u, u_old, eps, tau );
        Y_Matrix_AC<NV, DIM > Y( mesh.n4e, normals, mass, u, u_old, eps, tau );

        double mu = 0.;
        double mu_old = 0.;
        double diff_mu = 1.;
        double eps_stop = .1;
        double TOL_CG = 1e-06;
        
        double alpha;
        
        ManagedArray<double,1> rhs( w.size(), 0.);
        ManagedArray<double,1> temp( w.size(), 0.);
        
        temp = M * w;
        alpha = w * temp;
        for(int i = 0; i < w.size(); i++)
            w[i] = w[i] / sqrt( alpha );

        while ( sqrt( diff_mu * diff_mu) > eps_stop ) {
            
            rhs = M * w;
            pcg(Y, rhs, w, TOL_CG );
            
            temp = M * w;
            alpha = w * temp;
            for(int i = 0; i < w.size(); i++)
                w[i] = w[i] / sqrt( alpha );
            
            temp = Y * w;
            mu = w * temp;
            diff_mu = mu - mu_old;
            mu_old = mu;
            printf(" mu = %lf, diff_mu = %lf\n", mu, diff_mu);
        }
        
        double c_shift = Y.get_c_shift();
        
        return -mu + c_shift;
    }   
    
    void createEigValFile( const char *file ){
        FILE *fid;
        fid = fopen(file,"w");
        fclose(fid);
    }
    
    void printEigVal( const char *file, double t,  double lambda ){
        FILE *fid;
        fid = fopen(file,"a");  
        fprintf(fid,"%3.9f %3.9f",t, lambda);
        fprintf(fid,"\n");
        fclose(fid);
    }
    
    template< int NV, int DIM, class Decomposition>
    double calcEigenValue_parallel( MeshType &mesh, Decomposition &meshDecomp, ManagedArray<double, NV*DIM> &normals, ManagedArray<double, 1> &mass, ManagedArray<double,1> &u, ManagedArray<double, 1> &u_old, double eps, double tau, ManagedArray<double, 1> &w)
    {
        MassMatrix_AC<NV, DIM > M( mesh.n4e, normals, mass, u, u_old, eps, tau );
        Y_Matrix_AC<NV, DIM > Y( mesh.n4e, normals, mass, u, u_old, eps, tau );

        double mu = 0.;
        double mu_old = 0.;
        double diff_mu = 1.;
        double eps_stop = .5;
        double TOL_CG = 1e-06;
        
        double alpha;
        
        ManagedArray<double,1> rhs( w.size(), 0.);
        ManagedArray<double,1> temp( w.size(), 0.);
        
        temp = M * w;
        //meshDecomp.sumVertexVector(temp);
        //alpha = meshDecomp.scalar_product( w, temp);
        alpha = w * temp;
        for(int i = 0; i < w.size(); i++)
            w[i] = w[i] / sqrt( alpha );
            

        while ( sqrt( diff_mu * diff_mu) > eps_stop ) {
            
            rhs = M * w;
            meshDecomp.sumVertexVector(rhs);
            parallel_pcg(meshDecomp, Y, rhs, w, TOL_CG );
            
            temp = M * w;
            meshDecomp.sumVertexVector(temp);
            alpha = meshDecomp.scalar_product( w, temp );
            for(int i = 0; i < w.size(); i++)
                w[i] = w[i] / sqrt( alpha );
            
            temp = Y * w;
            meshDecomp.sumVertexVector(temp);
            mu = meshDecomp.scalar_product( w, temp );
            diff_mu = mu - mu_old;
            mu_old = mu;
        }
        
        double c_shift = Y.get_c_shift();
        
        return -mu + c_shift;
    }       
    
} // LAMBDA_HH conformingsimplexgrid
#endif // LAMBDA_HH

