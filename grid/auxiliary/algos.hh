#ifndef ALGOS_HH
#define ALGOS_HH

namespace conformingsimplexgrid {

    template< size_t N >
    struct Faculty {
      static const size_t value = N * Faculty< N-1 >::value;
    };

    template<>
    struct Faculty< 0 > {
      static const size_t value = 1;
    };

    template< size_t N, size_t K >
    struct BinomialCoefficient {
      static const size_t value = Faculty< N >::value / (Faculty< K >::value * Faculty< N-K >::value);
    };

    // Determinante
    template< class Array >
    double det(Array & matrix ) {
        double val;
        switch(matrix.dim()){
            case 1:
                val = matrix[0];
            case 2:
                val = (matrix[0*3+0]*matrix[1*3+1])-(matrix[0*3+1]*matrix[1*3+0]);
            case 3:
                val = (matrix[0*3+0]*matrix[1*3+1]*matrix[2*3+2])+(matrix[0*3+1]*matrix[1*3+2]*matrix[2*3+0])+(matrix[0*3+2]*matrix[1*3+0]*matrix[2*3+1])
                -(matrix[0*3+1]*matrix[1*3+0]*matrix[2*3+2])-(matrix[0*3+0]*matrix[1*3+2]*matrix[2*3+1])-(matrix[0*3+2]*matrix[1*3+1]*matrix[2*3+0]);
            default:
                assert( "det not implemented for N > 3" );
        }
    return val;			
    };
    
} //namespace conformingsimplexgrid
#endif // ALGOS_HH
