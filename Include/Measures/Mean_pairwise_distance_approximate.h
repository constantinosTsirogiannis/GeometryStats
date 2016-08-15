#ifndef MEAN_PAIRWISE_DISTANCE_APPROXIMATE_H
#define MEAN_PAIRWISE_DISTANCE_APPROXIMATE_H

#include<thread>
#include<cmath>

namespace FunctionalMeasures
{

template<class KERNEL_TYPE>
class Mean_pairwise_distance_approximate: public KERNEL_TYPE::Measure_base
{
 public:

  typedef KERNEL_TYPE                                                 Kernel;
  typedef typename Kernel::Number_type                               Number_type; 
  typedef typename Kernel::Point                                     Point;
  typedef typename Kernel::Compute_Euclidean_distance                Compute_Euclidean_distance;
  typedef typename Kernel::Numeric_traits                            Numeric_traits;
  typedef typename Numeric_traits::To_double                         To_double;
  typedef typename Numeric_traits::Square_root                       Square_root;
  typedef typename Numeric_traits::Protected_number_type             Protected_number_type;
  typedef typename Kernel::Mean_pairwise_distance_exact              Mean_pairwise_distance_exact;
  typedef typename Kernel::Well_separated_pair_decomposition         Well_separated_pair_decomposition;
  typedef typename Well_separated_pair_decomposition::Pointset_pair  Pointset_pair;
  typedef typename Kernel::Exception_functor                         Exception_functor;
  typedef typename Kernel::Exception_type                            Exception_type;

  Mean_pairwise_distance_approximate():_dim(-1){}
 
 public:

  template<class RangeIterator>
  void load_point_set( RangeIterator r_begin, RangeIterator r_end) 
  {
    _wspd.clear();    
    _wspd.load_point_set(r_begin, r_end);
    _n = _wspd.split_tree().number_of_sorted_points();
  }

  void set_epsilon(Number_type &epsilon, bool method_b = false)
  {
    _epsilon = epsilon;
    _total_dist_sum = Number_type(0.0);
    _total_dist_sum_sq = Number_type(0.0);
    _cost_interval = std::make_pair(Number_type(0.0),Number_type(0.0));
    _sorted_points.clear();
    _decomp.clear();

    if(method_b == true)
    {
      _wspd.extract_decomposition_epsilon(epsilon, std::back_inserter(_sorted_points), 
                                          std::back_inserter(_decomp));

      _total_dist_sum = _wspd.total_decomposition_cost();
      _cost_interval = _wspd.cost_interval();
    }
    else
    {
      Number_type separation_factor = Number_type(2.0)/_epsilon;

      _wspd.extract_decomposition(separation_factor, std::back_inserter(_sorted_points), std::back_inserter(_decomp));

      Compute_Euclidean_distance compute_distance;

      for(int i=0; i<_decomp.size(); i++)
      {
        int size_a = _decomp[i].set_a.second + 1 - _decomp[i].set_a.first,
            size_b = _decomp[i].set_b.second + 1 - _decomp[i].set_b.first;

        Number_type dist = compute_distance( _sorted_points[_decomp[i].set_a.first],
                                             _sorted_points[_decomp[i].set_b.first]);

        _total_dist_sum += dist*Number_type(size_a)*Number_type(size_b);      
        _total_dist_sum_sq += dist*dist*Number_type(size_a)*Number_type(size_b);      

      } // for(int i=0; i<_decomp.size(); i++)

   } // else of if(method_b == true)

  } // void precompute_decomposition(Number_type &epsilon)

  Number_type compute_expectation(int sample_size)
  { 
    if(sample_size < 2 || sample_size > _wspd.split_tree().number_of_sorted_points())
      return Number_type(0.0);

    return _total_dist_sum*Number_type(2.0)/(Number_type(_n)*Number_type(_n-1)); 
  }

  Number_type compute_expectation_lower_bound(int sample_size)
  { 
    if(sample_size < 2 || sample_size > _wspd.split_tree().number_of_sorted_points())
      return Number_type(0.0);

    return _cost_interval.first*Number_type(2.0)/(Number_type(_n)*Number_type(_n-1)); 
  }

  Number_type compute_expectation_upper_bound(int sample_size)
  { 
    if(sample_size < 2 || sample_size > _wspd.split_tree().number_of_sorted_points())
      return Number_type(0.0);

    return _cost_interval.second*Number_type(2.0)/(Number_type(_n)*Number_type(_n-1)); 
  }

  Number_type compute_variance(int sample_size)
  {
    if(sample_size < 2 || sample_size > _wspd.split_tree().number_of_sorted_points())
      return Number_type(0.0);

    int ps = _wspd.split_tree().number_of_sorted_points();

    Number_type fact = Number_type(4.0)/(Number_type(sample_size)*Number_type(sample_size)
                                         *Number_type(sample_size-1)*Number_type(sample_size-1)),
                coefficient_2 = Number_type(sample_size)*Number_type(sample_size-1)/
                              (Number_type(ps)*Number_type(ps-1)),
                coefficient_3 = Number_type(sample_size)*Number_type(sample_size-1)*Number_type(sample_size-2)/
                              (Number_type(ps)*Number_type(ps-1)*Number_type(ps-2)),
                coefficient_4 = Number_type(sample_size)*Number_type(sample_size-1)
                                *Number_type(sample_size-2)*Number_type(sample_size-3)/
                              (Number_type(ps)*Number_type(ps-1)
                              *Number_type(ps-2)*Number_type(ps-3));
    
    Number_type val(0.0);

    val += (_total_dist_sum*_total_dist_sum)*coefficient_4;
    val += _total_dist_sum_sq*(coefficient_2-coefficient_4);
    val += _wspd.extract_stored_sums()*(coefficient_3-coefficient_4);

    Number_type expec = compute_expectation(sample_size);

    return (fact*val) - (expec*expec);

  } // compute_variance(int sample_size)

  std::pair<Number_type,Number_type> cost_interval()
  { return _cost_interval;}

  int size_of_decomposition()
  { return _decomp.size();}

 private:

  Well_separated_pair_decomposition _wspd;
  std::vector<Pointset_pair> _decomp;
  std::vector<Point> _sorted_points;
  std::pair<Number_type,Number_type> _cost_interval;
  Number_type _epsilon,_total_dist_sum, _total_dist_sum_sq;
  int _dim,_n; 

}; // class Mean_pairwise_distance_exact

} // namespace FunctionalMeasures

//#include "Mean_pairwise_distance_approximate_impl.h"

#endif //MEAN_PAIRWISE_DISTANCE_APPROXIMATE_H
