

#include<lemon/tolerance.h>
#include<lemon/core.h>
namespace lemon {

  float Tolerance<float>::def_epsilon = static_cast<float>(1e-4);
  double Tolerance<double>::def_epsilon = 1e-10;
  long double Tolerance<long double>::def_epsilon = 1e-14;

#ifndef LEMON_ONLY_TEMPLATES
  const Invalid INVALID = Invalid();
#endif

} //namespace lemon
