/*
  Header file of C++ Matrix modules

  reference
    Biarane Stroustup(2014), The C++ Programming Language
*/
#ifndef MATRIX_MATRIX_H_
#define MATRIX_MATRIX_H_

#include <string>
#include <algorithm>
#include <numeric>
#include <array>
#include <type_traits>
#include <vector>
#include <iostream>
#include <cassert>
#include "estd.h"
using namespace estd;


// slice struct  // standard library also have slice
struct slice{
  slice() :start(-1), length(-1), stride(1) { }
  explicit slice(size_t s) :start(s), length(-1), stride(1) { }
  slice(size_t s, size_t l, size_t n=1) :start(s), length(l), stride(n) { }

  size_t operator()(size_t i) const { return start+i*stride; }

  static slice all;

  size_t start;
  size_t length;
  size_t stride;
};


// matrix_type judgement
template<typename T, size_t N, typename M>
constexpr bool Matrix_type();


// Matrix slice struct
template<size_t N>
struct Matrix_slice {
  Matrix_slice() = default;

  Matrix_slice(size_t offset, initializer_list<size_t> exts);
  Matrix_slice(size_t offset, initializer_list<size_t> exts,
               initializer_list<size_t> strs);

  template<typename... Dims> Matrix_slice(Dims... dims); // ???

  template<typename... Dims>
           //typename = Enable_if<All(Convertible<Dims,size_t>()...)>>
  size_t operator()(Dims... dims) const;

  size_t size;
  size_t start {0};
  array<size_t,N> extents;
  array<size_t,N> strides;
};


namespace Matrix_impl
{
  /*
    Matrix implementation namespace.
  */

  // derive extents (dimension's elements number)
  template<size_t N, typename List>
  array<size_t,N> derive_extents(const List& list);

  template<size_t N, typename I, typename List>
  Enable_if<(N>1),void> add_extents(I& first, const List& list);
  template<size_t N, typename I, typename List>
  Enable_if<(N==1),void> add_extents(I& first, const List& list);
  
  template<size_t N, typename List>
  bool checK_non_jagged(const List& list);

  // compute Matrix elements number and stride number
  template<int N>
  void compute_strides(Matrix_slice<N>& ms);

  // initializer_list into Matrix
  template<typename T, typename Vec>
  void add_list(const initializer_list<T>* first,
                const initializer_list<T>* last, Vec& vec);
  template<typename T, typename Vec>
  void add_list(const T* first, const T* last, Vec& vec);

  template<typename T, typename Vec>
  void insert_flat(initializer_list<T> list, Vec& vec);

  // confirm subscript number equal to dimension number and not over boundary
  template<size_t N, typename... Dims>
  bool check_bounds(const Matrix_slice<N>& slice, Dims... dims);

  // apply predicate to all variable number arguments
  constexpr bool All() { return true; };

  template<typename... Args>
  constexpr bool All(bool b, Args... args);

  //
  constexpr bool Some() { return false; }

  template<typename... Args>
  constexpr bool Some(bool b, Args... args);
  
  // check whether all subscript can be converted to size_t
  template<typename... Args>
  constexpr bool Requesting_element();

  // check whether slice can be used
  template<typename... Args>
  constexpr bool Requesting_slice();

  
  // slice
  // do slice
  template<size_t N, typename T, typename... Args>
  size_t do_slice(const Matrix_slice<N>& os, Matrix_slice<N>& ns,
                  const T& s, const Args&... args);

  template<size_t N>
  size_t do_slice(const Matrix_slice<N>& OS, Matrix_slice<N>& ns);
  // My original
  template<size_t N, size_t M>
  void slice_dim(size_t n, Matrix_slice<N>& desc, Matrix_slice<N-1>& row);
  
  template<size_t N, size_t M, typename T, typename... Args>
  size_t do_slice_dim(const Matrix_slice<N>& os, Matrix_slice<N>& ns,
                      const T& s);

  // Matrix_init struct declarations
  template<typename T, size_t N>
  struct Matrix_init;

  template<typename T>
  struct Matrix_init<T,1>;

  template<typename T> struct Matrix_init<T,0>; // N==0 is error
}


/*
  Matrix_initializer
*/
template<typename T, size_t N>
using Matrix_initializer = typename Matrix_impl::Matrix_init<T,N>::type;


/*
  Decralation of Matrix_base class
  (Public class of Matrix and Matrix_ref class)
*/


/*
  Decralation of Matrix_ref class
*/
template<typename T, size_t N>
class Matrix_ref {
public:
  Matrix_ref(const Matrix_slice<N>& s, T* p) :desc{s}, ptr{p} {}

  Matrix_ref() = default;
  Matrix_ref(Matrix_ref&&) = default;
  Matrix_ref& operator=(Matrix_ref&&) = default;
  Matrix_ref(const Matrix_ref&) = default;
  Matrix_ref& operator=(const Matrix_ref&) = default;
  ~Matrix_ref() = default;

  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  iterator begin() { ptr->elems.begin();}
  iterator end() { ptr->elems.end();}

  //template<typename U> Matrix& operator=(initializer_list<U>) = delete;

  size_t extent(size_t n) const { return desc.extents[n]; }
  size_t size() const{ return ptr->elems.size(); }
  const Matrix_slice<N>& descriptor() const { return desc; }

  // data access
  T* data() { return ptr->elems.data(); }
  const T* data() const { return ptr->elems.data(); }

  // Subscript and slicing operations
  // subscript access
  template<typename... Args> 
  Enable_if<Matrix_impl::Requesting_element<Args...>(), T&>
  operator()(Args... args);
  
  template<typename... Args>
  Enable_if<Matrix_impl::Requesting_element<Args...>(), const T&>
  operator()(Args... args) const;
  
  // slicing access
  template<typename... Args>
  Enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<T,N>>
  operator()(const Args&... args);
  
  template<typename... Args>
  Enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<const T,N>>
  operator()(const Args&... args) const;

  // row and column access
  Matrix_ref<T,N-1> operator[](size_t i){ return row(i); }
  Matrix_ref<const T,N-1> operator[](size_t i) const { return row(i); }

  Matrix_ref<T,N-1> row(size_t n);
  Matrix_ref<const T,N-1> row(size_t n) const;

  Matrix_ref<T,N-1> col(size_t n);
  Matrix_ref<const T,N-1> col(size_t n) const;


  // Arithmetic operations
  /*
  template<typename F> Matrix& apply(F f); // apply all x to f(x)

  template<typename M, typename F> // apply *this and m elements to f(x,mx)
  Enable_if<Matrix_type<M>(),Matrix<T,N>&> apply(const M& m, F f);

  Matrix& operator=(const T& value);

  Matrix& operator+=(const T& value);
  Matrix& operator-=(const T& value);
  Matrix& operator*=(const T& value);
  Matrix& operator/=(const T& value);
  Matrix& operator%=(const T& value);

  template<typename M>
  Enable_if<Matrix_type<M>(),Matrix&> operator+=(const M& x);
  template<typename M>
  Enable_if<Matrix_type<M>(),Matrix&> operator-=(const M& x);
  */
  
  
private:
  Matrix_slice<N> desc;  
  T* ptr;
};


/*
  Decralation of Matrix class
*/
template<typename T, size_t N>
class Matrix {
public:
  // Constructors and operator overloading
  static constexpr size_t oreder = N;
  using value_type = T;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  Matrix() = default;
  Matrix(Matrix&&) = default;
  Matrix& operator=(Matrix&&) = default;
  Matrix(const Matrix&) = default;
  Matrix& operator=(const Matrix&) = default;
  ~Matrix() = default;

  // construct Matrix class by Matrix_ref class
  template<typename U> Matrix(const Matrix_ref<U,N>&);
  template<typename U> Matrix& operator=(const Matrix_ref<U,N>&);

  // construct Matrix class by defining Matrix_slice class
  template<typename... Exts> explicit Matrix(Exts... exts);

  // construct Matrix class by  defineing all elements
  Matrix(Matrix_initializer<T,N>);
  Matrix& operator=(Matrix_initializer<T,N>);

  template<typename U> Matrix(initializer_list<U>) = delete;
  template<typename U> Matrix& operator=(initializer_list<U>) = delete;

  size_t extent(size_t n) const { return desc.extents[n]; }
  size_t size() const{ return elems.size(); }
  const Matrix_slice<N>& descriptor() const { return desc; }

  iterator begin() { elems.begin();}
  iterator end() { elems.end();}
  
  // data access
  T* data() { return elems.data(); }
  const T* data() const { return elems.data(); }

  // Subscript and slicing operations
  // subscript access
  template<typename... Args> 
  Enable_if<Matrix_impl::Requesting_element<Args...>(), T&>
  operator()(Args... args);

  template<typename... Args>
  Enable_if<Matrix_impl::Requesting_element<Args...>(), const T&>
  operator()(Args... args) const;

  // slicing access
  template<typename... Args>
  Enable_if<Matrix_impl::Requesting_slice<Args...>(), Matrix_ref<T,N>>
  operator()(const Args&... args);

  template<typename... Args>
  Enable_if<Matrix_impl::Requesting_slice<Args...>(),
                      Matrix_ref<const T,N>>
  operator()(const Args&... args) const;

  
  // row and column access
  // C-style subscript operation []
  Matrix_ref<T,N-1> operator[](size_t i){ return row(i); }
  Matrix_ref<const T,N-1> operator[](size_t i) const { return row(i); }

  Matrix_ref<T,N-1> row(size_t n);
  Matrix_ref<const T,N-1> row(size_t n) const;

  Matrix_ref<T,N-1> col(size_t n);
  Matrix_ref<const T,N-1> col(size_t n) const;

  // Arithmetic operations
  template<typename F> Matrix& apply(F f); // apply all x to f(x)

  template<typename M, typename F> // apply *this and m elements to f(x,mx)
  Enable_if<Matrix_type<M>(),Matrix<T,N>&> apply(const M& m, F f);

  Matrix& operator=(const T& value);

  Matrix& operator+=(const T& value);
  Matrix& operator-=(const T& value);
  Matrix& operator*=(const T& value);
  Matrix& operator/=(const T& value);
  Matrix& operator%=(const T& value);

  template<typename M>
  Enable_if<Matrix_type<M>(),Matrix&> operator+=(const M& x);
  template<typename M>
  Enable_if<Matrix_type<M>(),Matrix&> operator-=(const M& x);

private:
  Matrix_slice<N> desc;
  vector<T> elems;
  
};


// Matrix class definition when N == 0
template<typename T>
class Matrix<T,0>{
public:
  static constexpr size_t order = 0;
  using value_type = T;

  Matrix(const T& x) : elem(x) { }
  Matrix& operator=(const T& value) { elem = value; return *this; }

  // no access to row and column
  T& row(size_t i) = delete;
  T& col(size_t i) = delete;

  T& operator()() { return elem; }
  const T& operator()() const { return elem; }
  
  operator T&() { return elem; }   // ???
  operator const T&() { return elem; }
private:
  T elem;
};


// simple assignment operator for Matrix (scalar)
// +
template<typename T, size_t N>
Matrix<T,N> operator+(const Matrix<T,N>& m, const T& val);
// -
template<typename T, size_t N>
Matrix<T,N> operator-(const Matrix<T,N>& m, const T& val);
// *
template<typename T, size_t N>
Matrix<T,N> operator*(const Matrix<T,N>& m, const T& val);
// /
template<typename T, size_t N>
Matrix<T,N> operator/(const Matrix<T,N>& m, const T& val);


// simple assignment operator for Matrix_ref (scalar)
// +
template<typename T, size_t N>
Matrix_ref<T,N> operator+(const Matrix_ref<T,N>& m, const T& val);
// -
template<typename T, size_t N>
Matrix_ref<T,N> operator-(const Matrix_ref<T,N>& m, const T& val);
// *
template<typename T, size_t N>
Matrix_ref<T,N> operator*(const Matrix_ref<T,N>& m, const T& val);
// /
template<typename T, size_t N>
Matrix_ref<T,N> operator/(const Matrix_ref<T,N>& m, const T& val);


// simple assignment operator for Matrix (matrix)
// +
template<typename T, size_t N>
Matrix<T,N> operator+(const Matrix<T,N>& a, const Matrix<T,N>& b);
// -
template<typename T, size_t N>
Matrix<T,N> operator-(const Matrix<T,N>& a, const Matrix<T,N>& b);
//  m x 1 * n x 1 
template<typename T>
Matrix<T,2> operator*(const Matrix<T,1>& u, const Matrix<T,1>& v);
// n x m * m x 1
template<typename T>
Matrix<T,1> operator*(const Matrix<T,2>& m, const Matrix<T,1>& v);
// n x m * m x p
template<typename T>
Matrix<T,2> operator*(const Matrix<T,2>& m1, const Matrix<T,2>& m2);


// simple assigtnment operator for Matrix_ref (matrix)
// +
template<typename T, size_t N>
Matrix_ref<T,N> operator+(const Matrix_ref<T,N>& a, const Matrix_ref<T,N>& b);
// -
template<typename T, size_t N>
Matrix_ref<T,N> operator-(const Matrix_ref<T,N>& a, const Matrix_ref<T,N>& b);
//  m x 1 * n x 1 
template<typename T>
Matrix_ref<T,2> operator*(const Matrix_ref<T,1>& u, const Matrix_ref<T,1>& v);
// n x m * m x 1
template<typename T>
Matrix_ref<T,1> operator*(const Matrix_ref<T,2>& m, const Matrix_ref<T,1>& v);
// n x m * m x p
template<typename T>
Matrix_ref<T,2> operator*(const Matrix_ref<T,2>& m1, const Matrix_ref<T,2>& m2);

#endif // MATRIX_MATRIX_H_
