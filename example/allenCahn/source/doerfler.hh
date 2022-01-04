#ifndef DOERFLER_HH
#define DOERFLER_HH

#include<math.h>

namespace conformingsimplexgrid {
  
template< class Vector >
void doerfler_marking( Vector &marked, ManagedArray<double,1> &eta_space, double eta_max)
{
        double delta = .05;
        double theta = .25;
        double theta2 = 1. - delta;
        double eta = eta_space.sum();
        double sigma = 0.;

        while( sigma < theta*eta )
        {
                for( int i = 0; i < (int)marked.size(); i++ ) 
                    marked[i] = ( eta_space[i] > theta2*eta_max ? 1 : 0 );
                sigma = eta_space.sum( marked );
                theta2 -= delta;
        }
        printf("doerfler_marking with theta2 = %lf, sigma = %lf, eta = %lf\n", sqrt(theta2), sqrt(sigma), sqrt(eta));

}

template< class MeshType >
void doerfler_marking_coarse( MeshType &mesh, ManagedArray<double,1> &eta_sc)
{
        double delta = .05;
        double theta = .25;
        double theta2 = 1. - delta;
        double eta = eta_sc.sum();
        double eta_max = eta_sc.max();
        double sigma = 0.;

        while( sigma > theta*eta )
        {               
                for( int i = 0; i < (int)mesh.marked.size(); i++ ) 
                    mesh.marked[i] = ( eta_sc[i] < theta2*eta_max ? 1 : 0 );
                sigma = eta_sc.sum( mesh.marked );
                theta2 += delta;
        }
        printf("doerfler_marking_coarse with theta2 = %lf, sigma = %lf, eta = %lf\n", sqrt(theta2), sqrt(sigma), sqrt(eta));

}

    
} // DOERFLER_HH conformingsimplexgrid
#endif // DOERFLER_HH

