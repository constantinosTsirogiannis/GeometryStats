#ifndef CONSTRUCT_POINT_SET_H
#define CONSTRUCT_POINT_SET_H

//#include<Rcpp.h>
#include<vector>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<set>
#include<map>

namespace FunctionalMeasures{

template<class KERNEL_TYPE>
class Construct_point_set
{
 public:

  typedef KERNEL_TYPE                         Kernel;
  typedef typename Kernel::Number_type        Number_type;
  typedef typename Kernel::Point              Point;
  typedef typename Kernel::Exception_functor  Exception_functor;
  typedef typename Kernel::Exception_type     Exception_type;

 public:

  Construct_point_set(){}
  
  template<class OutputIterator>
  void operator()( char *filename, OutputIterator ot, int d=-1);

  template<class OutputIterator>
  void operator()( std::vector<std::string> &names, std::vector<std::vector<double> > &vals,  
                     OutputIterator ot, int d=-1);

  //template<class OutputIterator>
  //void operator()( std::vector<std::string> &names, NumericMatrix &vals,  
  //                   OutputIterator ot, int d=-1);

  void write_point_set_to_file(std::vector<Point> &points, char *filename);

  template<class OutputIterator>
  void construct_random_point_set(int size, int dim, OutputIterator ot, int seed=0)
  {
    srand(seed);

    for(int i=0; i<size; i++)
    {
      Point p;

      std::ostringstream os;

      os << "p_" << i;

      p.taxon = os.str();

      for(int j=0; j<dim; j++)
        p.push_back(Number_type(rand()%size)+Number_type(double(rand())/double(RAND_MAX))); 

      *ot++ = p;

    } // for(int i=0; i<size; i++)

  } // construct_random_point_set(...)


  template<class OutputIterator>
  void construct_random_values(int size, OutputIterator ot, int seed=0)
  {
    srand(seed);

    for(int i=0; i<size; i++)
      *ot++ = Number_type(rand()%size)/Number_type(size);

  } // construct_random_values(...)

  template<class OutputIterator>
  void construct_random_probability_values(int size, OutputIterator ot, int seed=0)
  {
    srand(seed);

    for(int i=0; i<size; i++)
    {
      Number_type res(0.0);

      while(res == Number_type(0.0) || res == Number_type(1.0))
        res = Number_type(rand()%size)/Number_type(size);

      *ot++ = res;
    }


  } // construct_random_probability_values(...)

  template<class OutputIterator>
  void construct_fixed_probability_value(int size, OutputIterator ot, Number_type prob)
  {
    for(int i=0; i<size; i++)
      *ot++ = prob;

  } // construct_random_probability_values(...)
  
}; // class Construct_point_set

} // namespace FunctionalMeasures

#include "Construct_point_set_impl.h"

#endif //CONSTRUCT_POINT_SET_H
