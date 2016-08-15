#ifndef CALCULATE_MOMENTS_WITH_SAMPLING_H
#define CALCULATE_MOMENTS_WITH_SAMPLING_H

#include<vector>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include<set>
#include<map>

namespace FunctionalMeasures{

template<class KERNEL_TYPE>
class Calculate_moments_with_sampling
{
 public:

  typedef KERNEL_TYPE                          Kernel;
  typedef typename Kernel::Number_type        Number_type;
  typedef typename Kernel::Point              Point;
  typedef typename Kernel::Exception_functor  Exception_functor;
  typedef typename Kernel::Exception_type     Exception_type;

 public:

  Calculate_moments_with_sampling(){}
  
  template<class RangeIterator>
  void load_point_set(RangeIterator r_begin, RangeIterator r_end)
  { 
    _points.clear();
    _points.insert(_points.begin(), r_begin, r_end);
  }

  int number_of_points()
  { return _points.size();}

  void set_random_seed(int seed)
  { srand(seed);}

  template<class OutputIterator>
  void select_random_sample(int sample_size, int n, OutputIterator ot)
  {
    std::vector<bool> vb;

    vb.assign(n,false);    

    int count=0;
   
    while(count < sample_size)
    {
      int index = rand()%n;

      if(vb[index] == false)
      {
        vb[index] = true;
        *ot++ = index;
        count++;
      }

    } // while(count < sample_size)

  } // select_random_sample(int sample_size, int n, OutputIterator ot)


  template<class OutputIterator>
  void select_random_sample(std::vector<Number_type> prbs, OutputIterator ot)
  {   
    for(int i=0; i < prbs.size(); i++)
    {
      Number_type toss = Number_type(rand())/Number_type(RAND_MAX);

      if(toss <= prbs[i])
        *ot++ = i;

    } // while(count < sample_size)

  } // select_random_sample( std::vector<Number_type> prbs, int n)

  /////////////////////////////////////////////////////////////////////
  // Operator that selects samples according to the fixed-size model //
  /////////////////////////////////////////////////////////////////////

  template <class MeasureType>
  std::pair<Number_type, Number_type>
  operator()(int sample_size, int number_of_samples, MeasureType &msr)
  {
    int n = number_of_points();
    Number_type sum(0.0), sum_sq(0.0);

    for(int i=0; i<number_of_samples; i++)
    {
      std::vector<int> vec;
   
      select_random_sample(sample_size, n, std::back_inserter(vec));
      Number_type val = msr(vec.begin(), vec.end());

      sum += val;
      sum_sq += val*val;
    }   

    sum = sum/Number_type(number_of_samples);  
    sum_sq = (sum_sq/Number_type(number_of_samples)) - (sum*sum);
 
    return std::make_pair(sum,sum_sq);

  } // operator()(int sample_size, int number_of_samples, Measure_type &msr)

  ////////////////////////////////////////////////////////////////////
  // Operator that selects samples according to the Bernoulli model //
  ////////////////////////////////////////////////////////////////////

  template <class MeasureType>
  std::pair<Number_type, Number_type>
  operator()(int number_of_samples, MeasureType &msr)
  {
    Number_type sum(0.0), sum_sq(0.0);

    std::vector<Number_type> prbs;

    for(int i=0; i<msr.number_of_elements(); i++)
      prbs.push_back(msr.prb[i]);

    for(int i=0; i<number_of_samples; i++)
    {
      std::vector<int> vec;
   
      select_random_sample(prbs, std::back_inserter(vec));
      Number_type val = msr(vec.begin(), vec.end());

      sum += val;
      sum_sq += val*val;
    }   

    sum = sum/Number_type(number_of_samples);  
    sum_sq = (sum_sq/Number_type(number_of_samples)) - (sum*sum);
 
    return std::make_pair(sum,sum_sq);

  } // operator()(int number_of_samples, Measure_type &msr)

 private:

  std::vector<Point> _points;

}; // class Calculate_moments_with_sampling

} // namespace FunctionalMeasures

#endif //CALCULATE_MOMENTS_WITH_SAMPLING_H
