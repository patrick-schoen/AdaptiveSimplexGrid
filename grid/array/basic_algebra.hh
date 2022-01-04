/** 
*
@file basic_algebra.hh
@brief defines algebraic operations for the ManagedArray class.
*/

#ifndef BASIC_ALGEBRA_HH
#define BASIC_ALGEBRA_HH
 
namespace conformingsimplexgrid {
    
    template <class T, class S, int n > T scalar_product(const ManagedArray< T, n> &a1, const ManagedArray< S, n> &a2)
    {
        assert( a1.size() == a2.size() );
        T res = 0;
        for( size_t i = 0; i < a1.size()*n; i++)
            res += a1[i] * a2[i];
        return res;    
    }
    template <class T, class S, int n >  ManagedArray< T, n> operator * (const ManagedArray< T, n> &a1, const ManagedArray< S, n> &a2)
    {
        assert( a1.size() == a2.size() );
        ManagedArray< T, n> res(a1);
        return res *= a2;
    }
    template <class T, int n, class POD > ManagedArray< T, n > operator * (const ManagedArray< T, n> &a1, const POD &value)
    {
        ManagedArray<T, n> res(a1);
        return res *= value;    
    }
    template <class T, int n, class POD > ManagedArray< T, n > operator * (const POD &value, const ManagedArray< T, n> &a1)
    {
        ManagedArray<T, n> res(a1);
        return res *= value;    
    }
    template <class T, int n, class S > ManagedArray< T, n > operator / (const ManagedArray< T, n> &a1, const ManagedArray< S, n> &a2)
    {
        assert( a1.size() == a2.size() );
        ManagedArray<T, n> res(a1);
        return res /= a2;    
    }
    template <class T, int n, class POD > ManagedArray< T, n > operator / (const ManagedArray< T, n> &a1, const POD &value)
    {
        ManagedArray<T, n> res(a1);
        return res /= value;    
    }
    template <class T, int n, class POD > ManagedArray< T, n > operator / (const POD &value, const ManagedArray< T, n> &a1)
    {
        ManagedArray<T, n> res(a1);
        for( size_t i = 0; i < a1.size()*n; i++)
            res[i] = value / a1[i];
        return res;    
    }
    template <class T, int n, class S > ManagedArray< T, n > operator + (const ManagedArray< T, n> &a1, const ManagedArray< S, n> &a2)
    {
        assert( a1.size() == a2.size() );
        ManagedArray<T, n> res(a1);
        return res += a2;    
    }
    template <class T, int n, class POD > ManagedArray< T, n > operator + (const ManagedArray< T, n> &a1, const POD &value)
    {
        ManagedArray<T, n> res(a1);
        return res += value;    
    }
    template <class T, int n, class POD > ManagedArray< T, n > operator + (const POD &value, const ManagedArray< T, n> &a1)
    {
        ManagedArray<T, n> res(a1);
        return res += value;    
    }
    template <class T, int n, class S > ManagedArray< T, n > operator - (const ManagedArray< T, n> &a1, const ManagedArray< S, n> &a2)
    {
        assert( a1.size() == a2.size() );
        ManagedArray<T, n> res(a1);
        return res -= a2;    
    }
    template <class T, int n, class POD > ManagedArray< T, n > operator - (const ManagedArray< T, n> &a1, const POD &value)
    {
        ManagedArray<T, n> res(a1);
        return res -= value;
    }
    template <class T, int n, class POD > ManagedArray< T, n > operator - (const POD &value, const ManagedArray< T, n> &a1)
    {
        ManagedArray<T, n> res(a1);
        res *= -1;
        return res += value;    
    }
    template <class T, int n, class POD > std::vector<bool> operator > (const ManagedArray< T, n> &a1, const POD &value)
    {
        std::vector<bool> res( a1.size()*n );
        for( size_t k = 0; k < n*a1.size(); k++)
            res[k] = ( a1[k] > value );
        return res;
    }
    template <class T, int n, class POD > std::vector<bool> operator > (const POD &value, const ManagedArray< T, n> &a1)
    {
        std::vector<bool> res( a1.size()*n );
        for( size_t k = 0; k < n*a1.size(); k++)
            res[k] = ( value > a1[k] );
        return res;
    }
    template <class T, int n, class POD > std::vector<bool> operator < (const ManagedArray< T, n> &a1, const POD &value)
    {
        std::vector<bool> res( a1.size()*n );
        for( size_t k = 0; k < n*a1.size(); k++)
            res[k] = ( a1[k] < value );
        return res;
    }
    template <class T, int n, class POD > std::vector<bool> operator < (const POD &value, const ManagedArray< T, n> &a1)
    {
        std::vector<bool> res( a1.size()*n );
        for( size_t k = 0; k < n*a1.size(); k++)
            res[k] = ( value < a1[k] );
        return res;
    }
    template <class T, class S, int n > ManagedArray< T, n > operator ^ (const ManagedArray< T, n > &a1, const S power )
    {
        ManagedArray<T, n> res( a1.size() );
        for( size_t k = 0; k < a1.size()*n; k++ )
            res[k] = pow( a1[k], power );
        return res;    
    }
    
} // namespace conformingsimplexgrid
#endif // BASIC_ALGEBRA_HH
 
 
