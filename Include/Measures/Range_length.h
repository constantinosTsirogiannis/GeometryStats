#ifndef RANGE_LENGTH_H
#define RANGE_LENGTH_H

#include<thread>
#include<cmath>

namespace FunctionalMeasures
{

template<class KERNEL_TYPE>
class Range_length : public KERNEL_TYPE::Measure_base_Poisson_binomial
{
 public:

  typedef KERNEL_TYPE                                      Kernel;
  typedef typename Kernel::Measure_base_Poisson_binomial  Base;
  typedef typename Kernel::Distribution_type              Distribution_type;
  typedef typename Kernel::Number_type                    Number_type;
  typedef typename Kernel::Numeric_traits                 Numeric_traits;
  typedef typename Numeric_traits::Protected_number_type  Protected_number_type;
  typedef typename Numeric_traits::Is_smaller             Is_smaller;
  typedef typename Kernel::Polynomial                     Polynomial;
  typedef typename Kernel::Polynomial_multiplication      Polynomial_multiplication;
  typedef typename Kernel::Exception_functor              Exception_functor;
  typedef typename Kernel::Exception_type                 Exception_type;

  Range_length()
  { initialise_measure();}

 private:

  // Computes all together the probability values f(x) = \binom{x}{s-1}/\binom{n}{s}
  void compute_all_hypergeometric_probabilities_s_minus_one( int sample_size);

  Number_type hypergeom_minus_one( int x );

  // Computes all together the probability values f(x) = \binom{x}{s-2}/\binom{n}{s}
  void compute_all_hypergeometric_probabilities_s_minus_two( int sample_size);

  Number_type hypergeom_minus_two( int x );

 public:

  void initialise_measure()
  {
    _last_sample_size_exp = -1;
    _last_sample_size_var = -1;
    _exp = Number_type(-1.0);
    _var = Number_type(-1.0);
  }

  void set_probability_distribution( Distribution_type distrib )
  { 
    initialise_measure();
    Base::set_probability_distribution(distrib); 
  }

  template<class RangeIterator>
  void load_values( RangeIterator r_begin, RangeIterator r_end);

  template<class RangeIterator>
  void load_values( RangeIterator r_begin_a, RangeIterator r_end_a,
                     RangeIterator r_begin_b, RangeIterator r_end_b);


  void execute_precomputations()
  {
    if(Base::probability_distribution() == Kernel::UNIFORM_FIXED_SIZE)
      return _execute_precomputations_uniform_fixed_size();

    if(Base::probability_distribution() == Kernel::POISSON_BINOMIAL)
      return _execute_precomputations_poisson_binomial();
  }

  void _execute_precomputations_uniform_fixed_size();


  void _execute_precomputations_poisson_binomial();

  template<class RangeIterator>
  Number_type operator()(RangeIterator r_begin, RangeIterator r_end);

  Number_type compute_expectation(int sample_size)
  {
    if(Base::probability_distribution() == Kernel::UNIFORM_FIXED_SIZE)
      return _compute_expectation_uniform_fixed_size(sample_size);

    if(Base::probability_distribution() == Kernel::POISSON_BINOMIAL)
      return _compute_expectation_poisson_binomial();

    return Number_type(-1.0);
  }


  Number_type compute_expectation()
  {
    assert(Base::probability_distribution() == Kernel::POISSON_BINOMIAL);
    return _compute_expectation_poisson_binomial();
  }


  Number_type compute_variance(int sample_size)
  {
    if(Base::probability_distribution() == Kernel::UNIFORM_FIXED_SIZE)
      return _compute_variance_uniform_fixed_size(sample_size);

    if(Base::probability_distribution() == Kernel::POISSON_BINOMIAL)
      return _compute_variance_poisson_binomial();

    return Number_type(-1.0);
  }

  Number_type compute_variance()
  {
    assert(Base::probability_distribution() == Kernel::POISSON_BINOMIAL);
    return _compute_variance_poisson_binomial();
  }

  Number_type element(int i)
  { return _vals[this->mapped_indices[i]];}

  int number_of_elements()
  { return _vals.size();}

 private:

  Number_type _compute_expectation_uniform_fixed_size(int sample_size);


  Number_type _compute_variance_uniform_fixed_size(int sample_size);


  Number_type _compute_expectation_poisson_binomial();


  Number_type _compute_variance_poisson_binomial();

 private:

  std::vector<Number_type> _vals;
  std::vector<Number_type>  _hypergeom_minus_one, _hypergeom_minus_two;

  int _last_sample_size_exp, _last_sample_size_var;
  Number_type _exp, _var;

}; // class Range_length

} // namespace FunctionalMeasures

#include "Range_length_impl.h"

#endif //RANGE_LENGTH_H
