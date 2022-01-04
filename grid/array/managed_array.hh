/**
@file managed_array.hh
@brief Implementation of the ManagedArray class.
*/
#ifndef CS_MANAGED_ARRAY_HH
#define CS_MANAGED_ARRAY_HH

///@cond
#include <iomanip>
///@endcond

namespace conformingsimplexgrid {
      
    ///@cond
    // forward declaration
    template< class T, int n > class ManagedArray;
    ///@endcond
    
    template< class T, int n = 1 >
    class ManagedArray /// ManagedArray class defines an array with values of type T, where the number of rows is dynamic and the number of columns is given by the static template argument n.
    {
        typedef ManagedArray < T, n > ThisType;

    public:
        typedef T value_type;
        
        static const int DIM = n;
        
        std::vector< T > array_;
        
        typedef typename std::vector<T>::iterator Iterator;
        typedef typename std::vector<T>::const_iterator ConstIterator;
         
        explicit ManagedArray( std::size_t size, const T &value = T() )
        : array_( size * n, value ) {}  
        
        ManagedArray ( std::vector<T> && values) : array_( std::forward< std::vector<T> >(values) ) {}
        
        ManagedArray ( const std::vector<T> & values) : array_(values) {}
        
        template< class R >
        ManagedArray ( const std::vector<R> & values ) : array_( values.size() )
        {
            static_assert( std::is_convertible<R,T>::value, "unknown konversion (ManagedArray)" );
            for(size_t i = 0; i < size(); ++i)
                array_[i] = T( values[i] );
        }

        ManagedArray () {}                                                      // Default - Constructor 

        ManagedArray( const ThisType &Other ) = default;
        
        ManagedArray & operator= ( const ThisType & ) = default; 
        
        ManagedArray( ThisType &&Other ) = default;
        
        ManagedArray & operator= ( ThisType && ) = default;
        
        template<class R>
        ManagedArray(const ManagedArray< R, n> &Other)                          // Copy - Constructor
        : array_( Other.size()*n ) 
        {
            std::copy( Other.begin(), Other.end(), array_.begin() );
        }
            
        template<class R>
        ThisType& operator= (const ManagedArray< R, n> &Other)                  // Copy - Constructor
        {
            array_.resize( Other.size()*n );
            std::copy( Other.begin(), Other.end(), array_.begin() );
            return(*this);
        }
            
        T& operator() ( std::size_t col, std::size_t row )                      /// Access values by column and row index.
        {
            assert( size() > 0 );
            assert( row < n );
            assert( (col*n + row) < n*size() );
            return array_[ col*n + row ];
        }

        const T& operator() ( std::size_t col, std::size_t row ) const
        {
            assert( size() > 0 );
            assert( row < n );
            assert( (col*n + row) < n*size() );
            return array_[ col*n + row ];
        }
        
        ManagedArray<T,n> operator() ( const std::vector<bool> &idx, bool pred = true ) 
        {
            ManagedArray<T,n> res( size() );
            assert( idx.size() == size() );
            int count = 0;
            for( int i = 0; i < size(); ++i )
                if( idx[i] == pred )
                    std::copy( array_.begin() + i*n, array_.begin() + (i+1)*n, res.begin() + n*count++);
            res.resize( count );
            return res;
        }
     
        T& operator[] ( const std::size_t index )                               /// Access values in the underlying std::vector directly.
        {
            assert( index < n*size() );
            return array_[ index ];
        }

        const T& operator[] ( const std::size_t index ) const    
        {
            assert( size() > 0 );
            assert( index < n*size() );
            return array_[ index ];
        }

    //    const T* rawData () const { return array_.data(); }                     //old function, will vanish in future.
    //    T *rawData () { return array_.data(); }
        
        const T* data () const { return array_.data(); } 
        T* data () /// returns pointer to the underlying std::vector array.
        { return array_.data(); }
         
        static std::size_t dim() /// returns static number of columns..
        { return std::size_t(n); } 

        ConstIterator begin () const { return array_.begin(); }
        Iterator begin () { return array_.begin(); }

        ConstIterator end () const { return array_.end(); }
        Iterator end () { return array_.end(); }

        void resize ( std::size_t newSize, const T &value = T() ) /// resize array by newSize many rows, fill new values by value.
        {
            array_.resize( n*newSize, value );
        }
        
        void push_back( const T &value ) ///adds an element to the end of underlying std::vector, i.e. ex. std::vector::push_back.
        {
            array_.push_back( value ); 
        }
        
        void reserve( std::size_t res ) ///reserves storage for the underlying sdt::vector, i.e. ex. std::vector::reserve.
        {
            array_.reserve( res );
        }
        
        void swap ( ThisType &other ) ///swaps the contents of the array the with ManagedArray other. 
        {
            array_.swap(other);
        }

        bool empty () const ///checks whether the container is empty.
        { return ( array_.empty() ); } 

        void clear () ///clears the contents.
        { array_.clear(); } 

        std::size_t size () const ///returns the number rows.
        { return size_t( array_.size() / n ); } 

        void fill ( T val ) ///sets all values of the array to val.
        {
            std::fill(begin(), end(), val ); 
        }

        void conjoin( ManagedArray<T,n> &other ) ///glues an array of same column size to the array A, such that A = [A; other].
        {
            int old_size = size();
            (*this).resize( size() + other.size() );
            for( int i = old_size; i < (int)size(); i++)
                for( int k = 0; k < n; k++)
                    (*this)(i,k) = other(i - old_size,k);
        }
      
        void print() const
        {
            std::cout << "printing ManagedArray: size = " << size() << " n = " << n << std::endl;
            for( size_t i = 0; i < size(); i++){
                std::cout << std::setw(5) << i << ": [";
                for( size_t k = 0; k < n; k++)
                    std::cout << " " << std::setw(5) << array_[i*n +k];
                std::cout << " ] " << std::endl;
            }
            std::cout << std::endl;
        }
        
        void print(const char * str) const ///prints contents of the arry to the screen for debugging, a string can be passed to distinguish between different outputs. 
        {
            std::cout << str << std::endl;
            (*this).print();
        }
        
        void print( size_t j ) const
        {
          for( size_t k = 0; k < n; k++)
             std::cout << " " << array_[j*n + k];
          std::cout << std::endl;
        }
                
        void print_as_int( ) const
        {
            std::cout << "printing ManagedArray: size = " << size() << " n = " << n << std::endl;
            for( size_t i = 0; i < size(); i++){
                std::cout << std::setw(5) << i << ": [";
                for( size_t k = 0; k < n; k++)
                    std::cout << " " << std::setw(5) << (int)array_[i*n +k];
                std::cout << " ] " << std::endl;
            }
            std::cout << std::endl;
        }

        T max() const ///returns maximal value of array A, ex. std::max_element
        {
            assert(size());
            return *std::max_element( begin(), end());
        }
        
        T min() const ///returns minimal value of array A, ex. std::max_element
        {
            assert(size());
            return *std::min_element( begin(), end());
        }
        
        T sum() const ///returns the sum of all values in A
        {
            T val = 0.;
            for(size_t k = 0; k < n*size(); k++)
                val += array_[k];
            return val;
        }
        
        T sum( const ManagedArray<int> &idx ) ///returns the sum of all values k in the underlying vector, where idx[k] > 0/
        {
            T val = 0.;
            for(size_t k = 0; k < n*size(); k++)
                if(idx[k] > 0 )
                    val += array_[k];
            return val;
        }
   
        bool operator== ( const ThisType &other ) const
        {
            assert( size() == other.size() );
            if( size() != other.size() ) return false;
            for( size_t k = 0; k < n*size(); ++k){
                if( array_[k] != other.array_[k] ) return false;
            }
            return true;
        }

        template <class Other>
        const ManagedArray& operator += (const ManagedArray<Other,n> &v)
        {
            assert( size() == v.size() );
            for( size_t i = 0; i < n*size(); i++)
                array_[i] += v[i];
            return *this;
        }

        template <class POD>
        const ManagedArray& operator += (const POD &value)
        {
            for( size_t i = 0; i < n*size(); i++)
                array_[i] += value;
            return *this;
        }

        template <class Other>
        const ManagedArray& operator -= (const ManagedArray<Other,n> &v)
        {
            assert( size() == v.size() );
            for( size_t i = 0; i < n*size(); i++)
                array_[i] -= v[i];
            return *this;
        }

        template <class POD>
        const ManagedArray& operator -= (const POD &value)
        {
            for( size_t i = 0; i < n*size(); i++)
                array_[i] -= value;
            return *this;
        }

        template <class Other>
        const ManagedArray& operator *= (const ManagedArray<Other,n> &v)
        {
            assert( size() == v.size() );
            for( size_t i = 0; i < n*size(); i++)
                array_[i] *= v[i];
            return *this;
        }

        template <class POD>
        const ManagedArray& operator *= (const POD &value)
        {
            for( size_t i = 0; i < n*size(); i++)
                array_[i] *= value;
            return *this;
        }

        template <class Other>
        const ManagedArray& operator /= (const ManagedArray<Other,n> &v)
        {
            assert( size() == v.size() );
            for( size_t i = 0; i < n*size(); i++)
                array_[i] /= v[i];
            return *this;
        }

        template <class POD>
        const ManagedArray& operator /= (const POD &value)
        {
            for( size_t i = 0; i < n*size(); i++)
                array_[i] /= value;
            return *this;
        }

    };
    
} // namespace conformingsimplexgrid
#endif // DATAMANAGEMENT_HH
