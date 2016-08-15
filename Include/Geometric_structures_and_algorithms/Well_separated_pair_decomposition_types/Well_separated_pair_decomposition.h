#ifndef WELL_SEPARATED_DECOMPOSITION_H
#define WELL_SEPARATED_DECOMPOSITION_H

namespace FunctionalMeasures
{

template <class KERNEL_TYPE>
class Well_separated_pair_decomposition
{
 public:

  typedef KERNEL_TYPE                          Kernel;
  typedef typename Kernel::Number_type        Number_type;
  typedef typename Kernel::Point              Point;
  typedef typename Kernel::Split_tree         Split_tree;
  typedef typename Split_tree::Pointset_pair  Pointset_pair;

  Well_separated_pair_decomposition(){}
  
  template <class RangeIterator>
  void load_point_set( RangeIterator rbegin, RangeIterator rend)
  { _tree.load_point_set(rbegin, rend);}

  template <class OutputIterator_a, class OutputIterator_b> 
  void extract_decomposition( Number_type &sf, OutputIterator_a ot_a, OutputIterator_b ot_b)
  { _tree.extract_decomposition(sf,ot_a,ot_b); }

  template <class OutputIterator_a, class OutputIterator_b> 
  void extract_decomposition_epsilon( Number_type &epsilon, OutputIterator_a ot_a, OutputIterator_b ot_b)
  { _tree.extract_decomposition_epsilon(epsilon,ot_a,ot_b); }

  bool is_valid_decomposition(Number_type &sf, std::vector<Point> &points, 
                               std::vector<Pointset_pair> &pprs)
  { return _tree.is_valid_decomposition(sf, points, pprs);}


  Number_type total_decomposition_cost()
  { return _tree.total_decomposition_cost();}

  std::pair<Number_type, Number_type> cost_interval()
  { return _tree.cost_interval();}

  Number_type extract_stored_sums()
  { return _tree.extract_stored_sums();}

  Split_tree& split_tree()
  { return _tree;}

  void clear()
  { _tree.clear();}

 private: 

  Split_tree _tree;

}; //class Well_separated_pair_decomposition 

} //namespace FunctionalMeasures

#endif // WELL_SEPARATED_DECOMPOSITION_H
