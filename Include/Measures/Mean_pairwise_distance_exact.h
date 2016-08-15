#ifndef MEAN_PAIRWISE_DISTANCE_EXACT_H
#define MEAN_PAIRWISE_DISTANCE_EXACT_H

#include<thread>
#include<cmath>

namespace FunctionalMeasures
{

template<class KERNEL_TYPE>
class Mean_pairwise_distance_exact: public KERNEL_TYPE::Measure_base
{
 public:

  typedef KERNEL_TYPE                                   Kernel;
  typedef typename Kernel::Number_type                 Number_type; 
  typedef typename Kernel::Point                       Point;
  typedef typename Kernel::Compute_Euclidean_distance  Compute_Euclidean_distance;
  typedef typename Kernel::Numeric_traits              Numeric_traits;
  typedef typename Numeric_traits::To_double           To_double;
  typedef typename Numeric_traits::Square_root         Square_root;
  typedef typename Kernel::Exception_functor           Exception_functor;
  typedef typename Kernel::Exception_type              Exception_type;

  Mean_pairwise_distance_exact():_dim(-1), _total_dist_sum(0.0){}
 
 private:


  struct Pairwise_distance_functor_complete
  {
    Pairwise_distance_functor_complete(){}

    Pairwise_distance_functor_complete
    ( int start, int length, int interval, std::vector<Point> *points, 
      std::vector<Number_type> *points_dist_sum, 
      std::vector<std::pair<Number_type,Number_type> > *vals, int index, bool forward): 
      _start(start), _length(length), _interval(interval), _ppoints(points),  
      _ppoints_dist_sum(points_dist_sum), _vals(vals), _index(index), _forward(forward) {}
 
    void operator()();
    
   private:

    bool _forward;
    int _start, _length, _interval, _index;
    std::vector<Point> *_ppoints; 
    std::vector<Number_type> *_ppoints_dist_sum; 
    std::vector<std::pair<Number_type,Number_type> > *_vals;

  }; // Pairwise_distance_functor_complete


  struct Pairwise_distance_functor_complete_2
  {
    Pairwise_distance_functor_complete_2(){}

    Pairwise_distance_functor_complete_2
    ( int first, int last, std::vector<Point> *points, 
      std::vector<Number_type> *points_dist_sum, 
      std::vector<std::pair<Number_type,Number_type> > *vals, int index, bool forward): 
      _first(first), _last(last), _ppoints(points),  
      _ppoints_dist_sum(points_dist_sum), _vals(vals), _index(index), _forward(forward) {}
 
    void operator()();
    
   private:

    bool _forward;
    int _first, _last, _index;
    std::vector<Point> *_ppoints; 
    std::vector<Number_type> *_ppoints_dist_sum; 
    std::vector<std::pair<Number_type,Number_type> > *_vals;

  }; // Pairwise_distance_functor_complete_2


  struct Pairwise_distance_functor_subset_complete
  {
    Pairwise_distance_functor_subset_complete(){}

    Pairwise_distance_functor_subset_complete
    ( int start, int length, int interval, std::vector<Point> *points, 
      std::vector<Number_type > *vals, int index): 
      _start(start),  _length(length), _interval(interval), 
      _ppoints(points), _vals(vals), _index(index) {}
 
    void operator()();

   private:

    int _start, _length, _interval, _index;
    std::vector<Point> *_ppoints; 
    std::vector<Number_type> *_vals;

  }; // Pairwise_distance_functor_subset_complete


  struct Pairwise_distance_functor_subset_bipartite
  {
    Pairwise_distance_functor_subset_bipartite(){}

    Pairwise_distance_functor_subset_bipartite
    ( int rbegin_a, int rend_a, int rbegin_b, int rend_b,
      std::vector<Point> *points, std::vector<Number_type> *vals, int index): 
      r_begin_a(rbegin_a),  r_end_a(rend_a), 
      r_begin_b(rbegin_b),  r_end_b(rend_b), _ppoints(points),  
      _vals(vals), _index(index) {}
 
    void operator()();

   private:

    int r_begin_a, r_end_a, r_begin_b, r_end_b, _index;
    std::vector<Point> *_ppoints; 
    std::vector<Number_type> *_vals;

  }; // Pairwise_distance_functor_subset_bipartite

  template<class RangeIterator>
  Number_type _operator_simple(RangeIterator r_begin, RangeIterator r_end);

  template<class RangeIterator>
  Number_type _operator_parallel(RangeIterator r_begin, RangeIterator r_end, int number_of_processors);

 public:

  template<class RangeIterator>
  void load_point_set( RangeIterator r_begin, RangeIterator r_end,
                        bool use_parallelisation=false, 
                        int number_of_processors = -1);

  void precompute_distances_simple();

  void precompute_distances_parallel(int number_of_processors);

  template<class RangeIterator>
  Number_type operator()(RangeIterator r_begin, RangeIterator r_end, 
                          bool use_parallelisation = false,
                          int number_of_processors = -1)
  {
    if(r_begin == r_end)
      return Number_type(0.0);

    if(use_parallelisation == false || number_of_processors == 1)
      return _operator_simple(r_begin, r_end);

    return _operator_parallel(r_begin, r_end, number_of_processors);
  }

  Number_type compute_expectation(int sample_size);

  Number_type compute_variance(int sample_size);

  Point element(int i)
  { return _points[i];}

  int number_of_elements()
  { return _points.size();}

 private:

  std::vector<Point> _points, _temp_points;
  std::vector<Number_type> _point_dist_sum;
  Number_type _total_dist_sum, _total_dist_sum_sq, _point_dist_sum_sq;
  int _dim; 

}; // class Mean_pairwise_distance_exact

} // namespace FunctionalMeasures

#include "Mean_pairwise_distance_exact_impl.h"

#endif //MEAN_PAIRWISE_DISTANCE_EXACT_H
