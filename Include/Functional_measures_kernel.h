#ifndef FUNCTIONAL_MEASURES_KERNEL_H
#define FUNCTIONAL_MEASURES_KERNEL_H

#include "Exception_related_types.h"
#include "Geometric_types/Point_d.h"
#include "Geometric_types/Bbox_d.h"
#include "Measures/Data_structures_for_bbox/Successor_product_structure.h"
#include "Measures/Data_structures_for_bbox/Successor_product_structure_fast.h"
#include "Measures/Data_structures_for_bbox/Product_tree.h"
#include "Geometric_structures_and_algorithms/Well_separated_pair_decomposition_types/Split_tree.h"
#include "Geometric_structures_and_algorithms/Well_separated_pair_decomposition_types/Well_separated_pair_decomposition.h"
#include "Geometric_structures_and_algorithms/Well_separated_pair_decomposition_types/Well_separated_pair_decomposition_slow.h"

#include "Geometric_functions/Compute_Euclidean_distance.h"
#include "Polynomial_related_types/Polynomial_rep.h"
#include "Polynomial_related_types/FFT.h"
#include "Construct_point_set.h"
#include "Calculate_moments_with_sampling.h"
#include "Measures/Measure_base/Measure_base.h"
#include "Measures/Bounding_box_volume_2D.h"
#include "Measures/Mean_pairwise_distance_exact.h"
#include "Measures/Mean_pairwise_distance_approximate.h"
//#include "Measures/Convex_hull_volume.h"
//#include "Measures/Convex_hull_volume_impl.h"
#include "Measures/Range_length.h"
#include "Numeric_traits_double.h"

template< class NTS = typename FunctionalMeasures::Numeric_traits_double>
struct Functional_measures_kernel
{
  typedef NTS                                              Numeric_traits;
  typedef Functional_measures_kernel<Numeric_traits>       Self;
  typedef typename Numeric_traits::Number_type            Number_type;
  typedef typename Numeric_traits::Protected_number_type  Protected_number_type;
  
  enum Distribution_type {UNIFORM_FIXED_SIZE, POISSON_BINOMIAL, POISSON_BINOMIAL_FIXED_SIZE, SINGLE_CONVEX_HULL};
  
  ////////////////////////////
  // Basic geometric types. //
  ////////////////////////////
  
  typedef typename FunctionalMeasures::Point_d<Self>    Point;
  typedef typename FunctionalMeasures::Bbox_d<Self>     Bbox_d;
    
  ///////////////////////////////////////////////
  // Basic Geometric Functions and Predicates  //
  ///////////////////////////////////////////////
  
  typedef typename FunctionalMeasures::Compute_Euclidean_distance<Self>    Compute_Euclidean_distance;


  /////////////////////////////////////////
  // Geometric Structures and Algorithms //
  /////////////////////////////////////////

  typedef typename FunctionalMeasures::Successor_product_structure<Self>      Successor_product_structure;
  typedef typename FunctionalMeasures::Successor_product_structure_fast<Self> Successor_product_structure_fast;
  typedef typename FunctionalMeasures::Product_tree<Self>                     Product_tree;

  typedef typename FunctionalMeasures::Split_tree<Self>  Split_tree;
  typedef typename Split_tree::Pointset_pair             Pointset_pair;

  typedef typename FunctionalMeasures::Well_separated_pair_decomposition<Self> 
                                                           Well_separated_pair_decomposition;

  typedef typename FunctionalMeasures::Well_separated_pair_decomposition_slow<Self> 
                                                           Well_separated_pair_decomposition_slow;

  //////////////////////////////////////
  // Classes that handle polynomials  //
  //////////////////////////////////////

  //typedef typename FunctionalMeasures::FFT<Self>                           FFT;
  //typedef typename FunctionalMeasures::Polynomial_rep<Self>                Polynomial;
  //typedef typename FunctionalMeasures::Polynomial_multiplication<Self>     Polynomial_multiplication;

  ////////////////////////////
  // Supplementary Classes  //
  ////////////////////////////

  typedef typename FunctionalMeasures::Construct_point_set<Self>              Construct_point_set;
  typedef typename FunctionalMeasures::Calculate_moments_with_sampling<Self>  Calculate_moments_with_sampling;

  ///////////////////////////
  // Measure base Classes  //
  ///////////////////////////

  typedef typename FunctionalMeasures::Measure_base<Self>                   Measure_base;
  typedef typename FunctionalMeasures::Measure_base_Poisson_binomial<Self>  Measure_base_Poisson_binomial;

  //////////////////////////////////////////////////////////////////////////////////////
  // Classes that Compute the Values and Moments of Functional Biodiversity Measures  //
  //////////////////////////////////////////////////////////////////////////////////////
  
  typedef typename FunctionalMeasures::Mean_pairwise_distance_exact<Self>        Mean_pairwise_distance_exact;
  typedef Mean_pairwise_distance_exact                                            Mean_pairwise_distance;

  typedef typename FunctionalMeasures::Mean_pairwise_distance_approximate<Self>  
                                                                          Mean_pairwise_distance_approximate;

  typedef typename FunctionalMeasures::Range_length<Self>                        Range_length;
  typedef typename FunctionalMeasures::Bounding_box_volume_2D<Self>              Bounding_box_volume_2D;
  //typedef typename FunctionalMeasures::Convex_hull_volume<Self>                  Convex_hull_volume;

  ///////////////////////////////////////////////////
  // Classes related to polynomial multiplication  //
  ///////////////////////////////////////////////////

  typedef typename FunctionalMeasures::FFT<Self>                        FFT;
  typedef typename FunctionalMeasures::Polynomial_rep<Self>             Polynomial;
  typedef typename FunctionalMeasures::Polynomial_multiplication<Self>  Polynomial_multiplication;

  //////////////////////////////////////////
  // Classes used for handling exceptions //
  //////////////////////////////////////////

  typedef typename ExceptionRelatedTypes::Exception_functor  Exception_functor;
  typedef typename ExceptionRelatedTypes::Exception_type     Exception_type;  

}; // class Functional_measures_kernel

#endif //FUNCTIONAL_MEASURES_KERNEL_H
