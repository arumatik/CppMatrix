#include <iostream>
#include <vector>
#include <array>
#include "estd.h"
using namespace estd;


class test{
  using iterator = std::vector<int>::iterator;
private:
  vector<int> a;
public:
  test(initializer_list<int> x) : a{x} {}
  test() = default;
  test(test&&) = default;
  test& operator=(test&&) = default;
  test(const test&) = default;
  test& operator=(const test&) = default;
  ~test() = default;
  

  iterator begin() { return a.begin(); }
  iterator end() { return a.end(); }

};

class test2{
private:
  test x;
  vector<float> y;
public:
  test2(test a, initializer_list<float> b) : x{a}, y{b} {}
};

template<typename X> 
Enable_if<Convertible<X,int>(),int> testfunc(){
  return 10;
}

template<typename X>
Enable_if<Convertible<X,test2>(),int> testfunc(){
  return 20;
}

int main(void){
  bool x {is_convertible<string,test2>::value};
  bool z {Convertible<long,float>()};
  bool y {is_convertible<int,long>::value};
  test a {1,2,3,4};
  const int b {10};
  int xx = testfunc<int>();
  int yy = testfunc<test2>();
  cout << xx << "\n";
  cout << yy << "\n";
}
