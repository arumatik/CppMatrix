/*
  Matrix library for C++.
  C++ 11 or heigher version is needed.
  
  Reference
    Birane Stroustrup(2014), The C++ Programming Language
*/

#include "Matrix.h"

/*
  Matrix slice implement
*/
// Matrix_slice constructors
template<size_t N>
Matrix_slice<N>::Matrix_slice(size_t offset, initializer_list<size_t> exts)
  : start(offset), extents(exts)
{
  size = 1;
  for ( size_t x : exts )
    size *= x;
}

template<size_t N>
Matrix_slice<N>::Matrix_slice(size_t offset, initializer_list<size_t> exts,
                              initializer_list<size_t> strs)
  : start(offset), extents(exts), strides(strs)
{
  size = 1;
  for ( size_t x : exts )
    size *= x;
}

// optimization is needed
template<size_t N> template<typename... Dims>
size_t Matrix_slice<N>::operator()(Dims... dims) const
{
  static_assert(sizeof...(Dims)==N,
                "Matrix_slice<N>:::operator(): diemension mismatch");
  
  size_t args[N] { size_t(dims)... };
  return start+inner_product(args,args+N,strides.begin(),size_t{0});
  // inner_product is in standard library <numeric>
}



/*
  namespace "Matrix_impl" definition
*/

namespace Matrix_impl
{

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
  Enable_if<(N>1),void> add_extents(I& first, const List& list)
  {
    assert(checK_non_jagged<N>(list));
    *first = list.size();
    add_extents<N-1>(++first,*list.begin()); // In the book, there is no ++.
  }

  template<size_t N, typename I, typename List>
  Enable_if<(N==1),void> add_extents(I& first, const List& list)
  {
    *first = list.size(); // In the book, *first++ = ...
  }

  template<size_t N, typename List>
  bool checK_non_jagged(const List& list)
  {
    auto i = list.begin();
    for ( auto j=i+1; j!=list.end(); ++j){
      if (derive_extents<N-1>(*i) != derive_extents<N-1>(*j))
        return false;}
    return true;
  }

  // compute Matrix elements number and stride number
  template<int N>
  void compute_strides(Matrix_slice<N>& ms) // Matrix_slice is in matrix.cpp
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

  template<typename T, typename Vec>
  void insert_flat(initializer_list<T> list, Vec& vec)
  {
    add_list(list.begin(), list.end(), vec);
  }

  // confirm subscript number equal to dimension number and not over boundary
  template<size_t N, typename... Dims>
  bool check_bounds(const Matrix_slice<N>& slice, Dims... dims)
  {
    size_t indexes[N] {size_t(dims)...};
    return equal(indexes, indexes+N, slice.extents.begin(), less<size_t> {});
  }

  
  // apply predicate to all variable number arguments
  template<typename... Args>
  constexpr bool All(bool b, Args... args)
  {
    return b && All(args...); 
  }

  //
  template<typename... Args>
  constexpr bool Some(bool b, Args... args)
  {
    return b || Some(args...);
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
    return All((Convertible<Args,size_t>() || Same<Args,slice>())...)
      && Some(Same<Args,slice>()...);  
  }
  
  // do slice
  template<size_t N, typename T, typename... Args>
  size_t do_slice(const Matrix_slice<N>& os, Matrix_slice<N>& ns,
                  const T& s, const Args&... args)
  {
    size_t m = do_slice_dim<sizeof...(Args)+1>(os, ns, s); // do_slice_dim ???
    size_t n = do_slice(os, ns, args...);
    return m * n;
  }

  template<size_t N>
  size_t do_slice(const Matrix_slice<N>& OS, Matrix_slice<N>& ns)
  {
    return 1;
  }

  // My original
  template<size_t N, typename T, size_t M>
  size_t do_slice_dim(const Matrix_slice<N>& os, Matrix_slice<N>& ns,
                      const T& s)
  {
    ns.size *= s.length;
    ns.extents[N-M] = s.length;
    ns.strides[N-M] = os.strides[N-M];
    return s.start;
  }
  
  // My original
  template<size_t N, size_t M>
  void slice_dim(size_t n, Matrix_slice<N>& desc, Matrix_slice<N-1>& row)
  {
    row.size = 1;
    size_t j {0};
    for ( size_t i=0; i!=N; ++i )
      if ( i != M ){
        row.size *= desc.extents[i];
        row.extents[j] = desc.extents[i];
        j++;
      }
    if ( M == 0 )
      row.strides = 1;
    else if ( M == 1 )
      row.strides = row.size / row.extents[0];
  }

}

/*
  Implement of Matrix_ref class
*/



/*
  Matrix class's functions definitions
*/
// Constructor from extent
template<typename T, size_t N> template<typename... Exts>
Matrix<T,N>::Matrix(Exts... exts)
  :desc{exts...}, elems(desc.size) {}

// Constructor from initialization_list
template<typename T, size_t N>
Matrix<T,N>::Matrix(Matrix_initializer<T,N> init)
{
  desc.extents = Matrix_impl::derive_extents<N>(init);
  Matrix_impl::compute_strides<N>(desc);
  elems.reserve(desc.size);
  Matrix_impl::insert_flat(init,elems);
  assert(elems.size() == desc.size);
}

// Assignment operation for initalization_list ???
/*
template<typename T, size_t N, typename U>
Matrix<T,N>& Matrix<T,N>::operator=(Matrix_initializer<T,N> init)
{
  desc.extents = Matrix_impl::derive_extents(init);

  Matrix_impl::compute_strides(desc);
  elems.reserve(desc.size);
  Matrix_impl:: insert_flat(init,elems);
  assert(elems.size() == desc.size);
}
*/

// Constructor from Matrix_ref
template<typename T, size_t N> template<typename U>
Matrix<T,N>::Matrix(const Matrix_ref<U,N>& x)
  :desc{x.desc}, elems(x.begin(), x.end())
{
  static_assert(Convertible<U,T>(),
                "Matrix constructor: incompatible element types");
}

// Assignment operation for Matrix_ref
template<typename T, size_t N> template<typename U>
Matrix<T,N>& Matrix<T,N>::operator=(const Matrix_ref<U,N>& x)
{
  static_assert(Convertible<U,T>(), "Matrix =: incompatible element types");

  desc = x.desc;
  elems.assign(x.begin(), x.end());
  return *this;
}


/*  Arithmetic operations  */
// apply all elements to the function
template<typename T, size_t N> template<typename F>
Matrix<T,N>& Matrix<T,N>::apply(F f)
{
  for (auto& x : elems) f(x);
  return *this;
}

// apply *this and m elements to f(x,mx)
template<typename T, size_t N> template<typename M, typename F>
Enable_if<Matrix_type<M>(),Matrix<T,N>&> Matrix<T,N>::apply(const M& m, F f)
{
  assert(same_extents(desc,m.descriptor()));
  for(auto i=this.begin(), j=m.begin(); i!=this.end(); ++i, ++j)
    f(*i,*j);
  return *this;
}

// scalar arithemtics
// +=
template<typename T, size_t N>
Matrix<T,N>& Matrix<T,N>::operator+=(const T& val)
{
  return apply([&](T& a) { a+=val; });
}
// -=
template<typename T, size_t N>
Matrix<T,N>& Matrix<T,N>::operator-=(const T& val)
{
  return apply([&](T& a) { a-=val; });
}
// *=
template<typename T, size_t N>
Matrix<T,N>& Matrix<T,N>::operator*=(const T& val)
{
  return apply([&](T& a) { a*=val; });
}
// /=
template<typename T, size_t N>
Matrix<T,N>& Matrix<T,N>::operator/=(const T& val)
{
  return apply([&](T& a) { a/=val; });
}
// %=
template<typename T, size_t N>
Matrix<T,N>& Matrix<T,N>::operator%=(const T& val)
{
  return apply([&](T& a) { a%=val; });
}
// = ???

// simple assignment operator (scalar)
// +
template<typename T, size_t N>
Matrix<T,N> operator+(const Matrix<T,N>& m, const T& val)
{
  Matrix<T,N> res = m;
  res += val;
  return res;
}
// -
template<typename T, size_t N>
Matrix<T,N> operator-(const Matrix<T,N>& m, const T& val)
{
  Matrix<T,N> res = m;
  res -= val;
  return res;
}
// *
template<typename T, size_t N>
Matrix<T,N> operator*(const Matrix<T,N>& m, const T& val)
{
  Matrix<T,N> res = m;
  res *= val;
  return res;
}
// /
template<typename T, size_t N>
Matrix<T,N> operator/(const Matrix<T,N>& m, const T& val)
{
  Matrix<T,N> res = m;
  res /= val;
  return res;
}


// Matrix arithmetics
// +=
template<typename T, size_t N> template<typename M>
Enable_if<Matrix_type<M>(),Matrix<T,N>&> Matrix<T,N>::operator+=(const M& m)
{
  static_assert(m.order==N,"+=: mismatched Matrix dimensions");
  assert(same_extents(desc,m.descriptor()));

  return apply(m, [](T& a, const Value_type<M>&b) { a += b; }); // value_type?
}
// -=

// simple assignment operator (matrix)
// +
template<typename T, size_t N>
Matrix<T,N> operator+(const Matrix<T,N>& a, const Matrix<T,N>& b)
{
  Matrix<T,N> res = a;
  res += b;
  return res;
}
/*
template<typename T, typename T2, size_t N,
typename RT = Matrix<Common_type<Value_type<T>,Value_type<T2>>, N>>
Matrix<RT,N> operator+(const Matrix<T,N>& a, const Matrix<T2,N>& b)
{
  Matrix<RT,N> res = a;
  res += b;
  return res;
}
 */
// -
template<typename T, size_t N>
Matrix<T,N> operator-(const Matrix<T,N>& a, const Matrix<T,N>& b)
{
  Matrix<T,N> res = a;
  res -= b;
  return res;
}
//  m x 1 * n x 1   // if you want to be more fast, remove overhead of res
template<typename T>
Matrix<T,2> operator*(const Matrix<T,1>& u, const Matrix<T,1>& v)
{
  const size_t n = u.extent(0);
  const size_t m = v.extent(0);
  Matrix<T,2> res(n,m);
  for (size_t i=0; i!=n; ++i)
    for (size_t j=0; j!=m; ++j)
      res(i,j) = u[i]*v[j];
  return res;
}
// n x m * m x 11
template<typename T>
Matrix<T,1> operator*(const Matrix<T,2>& m, const Matrix<T,1>& v)
{
  assert(m.extent(1)==v.extent(0));
  const size_t nr = m.extent(0);
  const size_t nc = m.extent(1);
  Matrix<T,1> res(nr);
  for (size_t i=0; i!=nr; ++i)
    for (size_t j=0; j!=nc; ++j)
      res(i) += m(i,j)*v(j);
  return res;
}
// n x m * m x p  // by many means, it is possible to be more compatible
template<typename T>
Matrix<T,2> operator*(const Matrix<T,2>& m1, const Matrix<T,2>& m2)
{
  const size_t nr = m1.extent(0);
  const size_t nc = m1.extent(1);
  assert(nc==m2.extent(0));
  const size_t p = m2.extent(1);
  Matrix<T,2> res(nr,p);
  for (size_t i=0; i!=nr; ++i)
    for (size_t j=0; j!=p; ++j)
      for (size_t k=0; k!=nc; ++k)
        res(i,j) += m1(i,k) * m2(k,j);
  return res;
}

// Arithmetic operations for Matrix_ref
// += ???
// -= ???
// *= ???
// /= ???
// %= ???
// = ???
// +
template<typename T, size_t N>
Matrix_ref<T,N> operator+(const Matrix_ref<T,N>& x, const T& n)
{
  Matrix_ref<T,N> res = x;
  res += n;
  return res;
}
// -


/*
  Judge whether type is Matrix.
 */

template<typename T, size_t N, typename M>
constexpr bool Matrix_type()
{
  return Same<M,Matrix<T,N>> || Same<M,Matrix_ref<T,N>>;
}


/*
  Access to Matrix
*/
// access to row
template<typename T, size_t N>
Matrix_ref<T,N-1> Matrix<T,N>::row(size_t n)
{
  assert(n<elems.size());  // no definition of rows() in above.
  Matrix_slice<N-1> row;
  Matrix_impl::slice_dim<0>(n, desc, row);
  return {row,data()}; // Matrix_ref's arguments are Matrix_slice and pointer
}

/*
template<typename T>
T& Matrix<T,1>::row(size_t i)
{
  return elems[i];
}
*/

//access to column
template<typename T, size_t N>
Matrix_ref<T,N-1> Matrix<T,N>::col(size_t n)
{
  assert(n<elems[0]); // cols() ???
  Matrix_slice<N-1> col;
  Matrix_impl::slice_dim<1>(n, desc, col);
  return {col,data()};
}

/*
template<typename T, size_t N>
T& Matrix<T,1>::col(size_t i)
{
  return elems[i];  // same as row, should be corrected later.
}
*/

// subscript operations
template<typename T, size_t N> template<typename... Args>
Enable_if<Matrix_impl::Requesting_element<Args...>(), T&>
Matrix<T,N>::operator()(Args... args)
{
  assert(Matrix_impl::check_bounds(desc, args...));
  return *(data() + desc(args...));
}

// slice subscript operations
template<typename T, size_t N> template<typename... Args>
Enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<T,N>>
Matrix<T,N>::operator()(const Args&... args)
{
  Matrix_slice<N> d;
  d.start = Matrix_impl::do_slice(desc,d,args...);
  return {d, data()};
}


int main(void){

  Matrix<int,2> a {
    {1,2,3},
    {4,5,6}
  };

  Matrix_ref<int,1> b = a(1,slice(1,2));
  /*
  for ( size_t i : a(1,slice(1,2)) )
    cout << i << "\n";
  */
}

