/*
  Matrix implement class
*/
#include <array>
#include <valarray>
#include <type_traits>
#include <iostream>
#include "matrix.h"
using namespace std;

namespace Matrix_impl
{
  /*
    function list ( removable after definition is written )
  */
  // Requesting_element ???
  // Requesting_slice ???
  // slice_dim<0>

  // forward decralations
  template<size_t N, typename I, typename List>
  enable_if<(N>1),void> add_extents(I& first, const List& list);
  
  template<size_t N, typename List>
  bool checK_non_jagged(const List& list);

  
  // Matrix_init class
  template<typename T, size_t N>
  struct Matrix_init {
    using type = initializer_list<typename Matrix_init<T,N-1>::type>;
  };

  template<typename T>
  struct Matrix_init<T,1>{
    using type = initializer_list<T>;
  };

  template<typename T> struct Matrix_init<T,0>; // N==0 is error

  // derive_extents
  template<size_t N, typename List>
  array<size_t,N> derive_extents(const List& list)
  {
    array<size_t,N> a;
    auto f = a.begin();
    add_extents<N>(f,list);
    return a;
  }

  template<size_t N, typename I, typename List>
  enable_if<(N>1),void> add_extents(I& first, const List& list)
  {
    assert(checK_non_jagged<N>(list));
    *first = list.size(); 
    add_extents<N-1>(first++,*list.begin()); // In the book, there is no ++.
  }

  template<size_t N, typename I, typename List>
  enable_if<(N==1),void> add_extents(I& first, const List& list)
  {
    *first = list.size(); // In the book, *first++ = ...
  }

  template<size_t N, typename List>
  bool checK_non_jagged(const List& list)
  {
    auto i = list.begin();
    for ( auto j=i+1; j!=list.end(); ++j)
      if (derive_extents<N-1>(*i) != derive_extents<N-1>(*j))
        return false;
    return true;
  }

  // compute Matrix elements number and stride number
  template<int N>
  void computs_stride(Matrix_slice<N>& ms) // Matrix_slice is in matrix.cpp
  {
    size_t st = 1;
    for (int i=N-1; i>=0; --i){
      ms.strides[i] = st;
      st *= ms.extents[i];
    }
    ms.size = st;
  }

  // initializer_list into Matrix
  template<typename T, typename Vec>
  void insert_flat(initializer_list<T> list, Vec& vec)
  {
    add_list(list.begin(), list.end(), vec);
  }

  template<typename T, typename Vec>
  void add_list(const initializer_list<T>* first,
                const initializer_list<T>* last, Vec& vec)
  {
    for (;first!=last; ++first)
      add_list(first->begin(), first->end(), vec);
  }

  template<typename T, typename Vec>
  void add_list(const T* first, const T* last, Vec& vec)
  {
    vec.insert(vec.end(), first, last);  // vector::insert
  }

  // confirm subscript number equal to dimension number and not over boundary
  template<size_t N, typename... Dims>
  bool check_bounds(const Matrix_slice<N>& slice, Dims... dims)
  {
    size_t indexes[N] {size_t(dims)...};
    return equal(indexes, indexes+N, slice.extents.begin(), less<size_t> {});
  }

  
  // apply predicate to all variable number arguments
  constexpr bool All() { retrun true; }

  template<typename... Args>
  constexpr bool All(bool b, Args... args)
  {
    return b && All(args...); // Convertible is not needed ???
  }

  // check whether all subscript can be converted to size_t
  template<typename... Args>
  constexpr bool Requesting_element()
  {
    return All(Convertible<Args,size_t>()...);
    // concept version of is_convertile
  }

  //
  template<typename... Args>
  constexpr bool Requesting_slice()
  {
    // slice is not declared
    return All((Convertible<Args,size_t>() || Same<Args,slice>())...)
      && Some(Same<Args,slice>()...);  // what is Some and Same ???
  }
  
  // do slice
  template<size_t N, typename T, typename... Args>
  size_t do_slice(const Matrix_slice<N>& os, Matrix_slice<N>& ns,
                  const T& s, const Args&... args)
  {
    size_t m = do_slice_dim<sizeof...(Args)+1>(os, ns, s); // do_slice_dim ???
    size_t n = do_slice(os, ns, args...);
    return m + n;
  }

  template<size_t N>
  size_t do_slice(const Matrix_slice<N>& OS, matrix_slice<N>& ns)
  {
    return 0;
  }
}
