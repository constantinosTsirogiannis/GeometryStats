#ifndef MEASURE_BASE_H
#define MEASURE_BASE_H

#include<vector>
#include<cmath>

namespace FunctionalMeasures {

template< class KernelType >
struct Measure_base
{
  typedef KernelType                           Kernel;
  typedef typename Kernel::Number_type        Number_type;
  typedef typename Kernel::Point              Point;
  typedef typename Kernel::Distribution_type  Distribution_type;
  typedef typename Kernel::Exception_type     Exception_type;
  typedef typename Kernel::Exception_functor  Exception_functor; 

 public: 
 
  Measure_base():_distribution(Kernel::UNIFORM_FIXED_SIZE){}

  void set_probability_distribution( Distribution_type distrib )
  { _distribution = distrib; }

  Distribution_type probability_distribution()
  { return _distribution; }

 protected:

  std::vector<int>        mapped_indices; // Stores the original order (ranking)
                                          // of the measure's elements. Needed for the
                                          // measures where input elements have to be rearranged,
                                          // in order to compute the moments efficiently.  
 
 protected:

  // Function that reads a list of integers from a file,
  // considered to be samples sizes that are used as input
  // for computing the expectation and deviation of a measure
  // on a given set of points. 
  template < class OutputIterator >
  void _read_sample_sizes_from_file(char *filename, std::vector<Point> &points, OutputIterator ot);
  
  // Input:  A point set and a range of iterators that indicate a list of species names (in std::string format).
  // Output: The value of the current measure for this set of species.
  template <class RangeIterator, class Measure>    
  Number_type _list_query(std::vector<Point> &points, RangeIterator rbegin, RangeIterator rend, Measure &msr);
  
  // Input: A vector of points, an array with the species names for these points and a matrix such that: 
  // each column corresponds to one of these species and each row indicates a sample of these 
  // species for which we want to compute the distance measure. A certain species is considered 
  // as part of the i-th sample if in the i-th row 
  // of the matrix there is a '1' at the column that corresponds to this species (otherwise a '0').

  // The fifth argument is a boolean value that indicates if we just want to compute the value of 
  // the measure for the indicated sample of species, or if we want to "standardise" this value;
  // that is to subtract from this value the mean of the measure for all samples of the same size
  // and then divide by the standard deviation.

  // The sixth argument is an output iterator of the type std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the (standardised) value of the measure for the sample that is described in the i-th row of the input matrix.

  template < class Measure, class OutputIterator>
  int _matrix_query( std::vector<Point> &points, std::vector<std::string> &names, 
                     std::vector< std::vector<bool> > &matrix, 
                     Measure &msr, bool standardised, OutputIterator ot );
  
  // Input: A csv file which stores a matrix where each column corresponds to a species of the point set
  // and each row indicates a sample of these species for which we want to compute the
  // distance measure. A certain species is considered as part of the i-th sample if in the i-th row 
  // of the matrix there is a '1' at the column that corresponds to this species (otherwise there is a '0').

  // The fourth argument is a boolean value that indicates if we just want to compute the value of 
  // the measure for the indicated sample of species, or if we want to "standardise" this value;
  // that is to subtract from this value the mean of the measure for all samples of the same size
  // and then divide by the standard deviation.

  // The fifth argument is an output iterator of the type std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the (standardised) value of the measure for the sample that is described in the i-th row of the input matrix.
  template <class Measure, class OutputIterator>
  int _csv_matrix_query( std::vector<Point> &points, char *filename, Measure &msr, bool standardised, OutputIterator ot );

 private:
    
  Distribution_type        _distribution;

}; // Measure_base

template< class KernelType >
struct Measure_base_Poisson_binomial: public KernelType::Measure_base
{
  typedef KernelType                     Kernel;
  typedef typename Kernel::Number_type  Number_type;

  std::vector<Number_type> prb;

  /////////////////////////////////////////////////////////////////////////////////////////
  // Next function loads the set of probability values to a vector inside this object.   //
  // RangeIterator should be of type std::vector<Number_type>::iterator                  //
  ///////////////////////////////////////////////////////////////////////////////////////// 

  template< class RangeIterator>
  void insert_probability_values(RangeIterator r_begin, RangeIterator r_end)
  { 
    prb.clear();
    prb.insert(prb.begin(), r_begin, r_end);
  }

  //////////////////////////////////////////////////////////////////
  // Fuction for accessing the probability value of each point.   //
  //////////////////////////////////////////////////////////////////

  Number_type element_probability(int i)
  { 
    if(this->prb.size()<i)
    {
      std::cout<<"Warning: no assigned probability to point "<<i<<std::endl;
      return Number_type(-1.0);
    }

    if(this->mapped_indices.size()==0)
      return prb[i];

    return prb[this->mapped_indices[i]];
  }

}; // struct Measure_base_Poisson_binomial

} // namespace FunctionalMeasures

#include"Measure_base_impl.h"

#endif //MEASURE_BASE_H
