
#ifndef NUMERIC_TRAITS_DOUBLE_H
#define NUMERIC_TRAITS_DOUBLE_H

#include<cmath>
#include "Numeric_traits_types/Protected_number_type.h" 

namespace FunctionalMeasures 
{
  struct Numeric_traits_double
  {
    typedef double                                                      Number_type;
    typedef Numeric_traits_double                                        Self;

    class Is_exact
    {
      public:
 
        bool operator()(void)
        { return false;}
    };

    class To_double
    {
      public:

        double operator()(double x)
        { return x; }
    };

    class Power
    {
      public:

        double operator()(double x, int k)
        { return std::pow(x,k); }
    };

    class Ceiling
    {
      public:

        double operator()(double x)
        { return ceil(x); }
    };

    class Square_root
    {
      public:

        double operator()(double x)
        { return std::sqrt(x); }
    };

    class Absolute_value
    {
      public:
 
        double operator()(double x)
        { return std::abs(x);}
    };

    class Is_smaller
    {
      public:

        bool operator()(const double &a, const double &b) const
        { return a<b; }
    };

    class Cosine
    {
      public:
 
        double operator()(double x)
        { return std::cos(x);}
    };

    class Sine
    {
      public:
 
        double operator()(double x)
        { return std::sin(x);}
    };


    /////////////////////////////////////////////////////////////////////
    // Type used for preserving accuracy in polynomial multiplications //
    /////////////////////////////////////////////////////////////////////

    typedef typename FunctionalMeasures::Protected_number_type<Self>  Protected_number_type;

  }; // struct Numeric_traits_double

} //namespace FunctionalMeasures

#endif // NUMERIC_TRAITS_DOUBLE_H
