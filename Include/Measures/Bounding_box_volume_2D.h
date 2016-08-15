#ifndef BOUNDING_BOX_VOLUME_2D_H
#define BOUNDING_BOX_VOLUME_2D_H

#include<vector>

namespace FunctionalMeasures
{

template<class KERNEL_TYPE>
class Bounding_box_volume_2D : public KERNEL_TYPE::Measure_base_Poisson_binomial
{

 public:

  typedef KERNEL_TYPE                                         Kernel;
  typedef typename Kernel::Measure_base_Poisson_binomial     Base;
  typedef typename Kernel::Distribution_type                 Distribution_type;
  typedef typename Kernel::Number_type                       Number_type;
  typedef typename Kernel::Numeric_traits                    Numeric_traits;
  typedef typename Numeric_traits::Protected_number_type     Protected_number_type;
  typedef typename Kernel::Point                             Point;
  typedef typename Numeric_traits::Is_smaller                Is_smaller;
  typedef typename Kernel::Exception_functor                 Exception_functor;
  typedef typename Kernel::Exception_type                    Exception_type;
  typedef typename Kernel::Product_tree                      Product_tree;
  typedef typename Kernel::Successor_product_structure       Successor_product_structure;
  typedef typename Kernel::Successor_product_structure_fast  Successor_product_structure_fast;

  struct Extended_point
  {
    Point pt;
    int index;
    Protected_number_type prod;

  }; // Extended_point

  struct Has_smaller_x
  {
    bool operator()(const Extended_point &p, const Extended_point &q ) const
    {
      if(p.pt[0] < q.pt[0])
        return true;

      return false;
    }

  }; // Has_smaller_x

  struct Has_smaller_y
  {
    bool operator()(const Extended_point &p, const Extended_point &q ) const
    {
      if(p.pt[1] < q.pt[1])
        return true;

      return false;
    }

  }; // Has_smaller_y
  
  Bounding_box_volume_2D():_dim(-1)
  { initialise_measure();}

  void initialise_measure()
  {
    _exp = Number_type(-1.0);
    _var = Number_type(-1.0);
  }

  Point element(int i)
  { return _points[i];}

  int number_of_elements()
  { return _points.size();}

  template<class RangeIterator>
  void load_point_set( RangeIterator r_begin, RangeIterator r_end)
  {
     this->prb.clear();
     _pixs.clear(); 
     _piys.clear();
     _map_index_to_rank_y.clear(); 
     _map_rank_x_to_index.clear();

     _points.clear();
     _points.insert(_points.begin(),r_begin,r_end);

     if(_points.size() > 0)
      _dim = _points[0].dim();
     else
      _dim = -1;

     _exp= Number_type(-1.0);
     _var= Number_type(-1.0);

  } // load_point_set(...)

  template<class RangeIteratorPoints, class RangeIteratorProbabilities>
  void load_point_set( RangeIteratorPoints pts_begin, RangeIteratorPoints pts_end, 
                        RangeIteratorProbabilities prbs_begin, RangeIteratorProbabilities prbs_end)
  {
    assert(Base::probability_distribution() == Kernel::POISSON_BINOMIAL || 
           Base::probability_distribution() == Kernel::POISSON_BINOMIAL_FIXED_SIZE);

    this->prb.clear();
    _pixs.clear(); 
    _piys.clear();
    _map_index_to_rank_y.clear(); 
    _map_rank_x_to_index.clear();

    _points.clear();
    _points.insert(_points.begin(),pts_begin,pts_end);
    this->prb.insert(this->prb.begin(),prbs_begin,prbs_end);

    if(_points.size() > 0)
     _dim = _points[0].dim();
    else
     _dim = -1;

    _exp= Number_type(-1.0);
    _var= Number_type(-1.0);

  } // load_point_set(...)

  template<class RangeIterator>
  Number_type operator()(RangeIterator r_begin, RangeIterator r_end);

  Number_type compute_expectation()
  {
    assert(Base::probability_distribution() == Kernel::POISSON_BINOMIAL);

    if(_exp>= Number_type(0.0))
      return _exp;
   
    _execute_precomputations_poisson_binomial();
    return _exp;
  }

  Number_type compute_variance()
  {
    assert(Base::probability_distribution() == Kernel::POISSON_BINOMIAL);
    return Number_type(-1.0);
  }

  void set_probability_distribution( Distribution_type distrib )
  { 
    initialise_measure();
    Base::set_probability_distribution(distrib); 
  }

  void execute_precomputations()
  {
    if(Base::probability_distribution() == Kernel::UNIFORM_FIXED_SIZE)
      return;

    if(Base::probability_distribution() == Kernel::POISSON_BINOMIAL)
      return _execute_precomputations_poisson_binomial();
  }

 private:

  void _execute_precomputations_poisson_binomial()
  { _exp = _compute_all_max_products(); }

  Number_type _compute_all_max_products();

  Protected_number_type _compute_max_products_two_distinct_points();

  Protected_number_type _compute_max_products_single_point();

  Protected_number_type _compute_max_products_single_point_slow();

  Protected_number_type _compute_max_products_single_point_fast();

 private:

  std::vector<Point>  _points;
  std::vector<Protected_number_type>  _pixs, _piys;
  std::vector<int> _map_index_to_rank_y, _map_rank_x_to_index;
  int _dim;
  Number_type _exp,_var;

}; // class Bounding_box_volume_2D

} // namespace FunctionalMeasures

#include "Bounding_box_volume_2D_impl.h"

#endif // BOUNDING_BOX_VOLUME_2D_H
