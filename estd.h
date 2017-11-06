/*
  Extended standard library

  reference
    Birane Stroustrup(2014), The C++ Programming Language
*/
#include <typeinfo>

namespace estd{
  using namespace std;

  
  /*
    Type Functions' extended definitions
  */
  template<typename X, typename Y>
  constexpr bool Convertible(){
    return is_convertible<X,Y>::value;
  }

  template<typename X, typename Y>
  constexpr bool Same(){
    return is_same<X,Y>::value;
  }

  template<bool X, typename Y = void>
  using Enable_if = typename enable_if<X,Y>::type;

  
  /*
    STL Container Member Types' definitions
  */
  template<typename M>
  using Value_type = typename M::value_type;

}
