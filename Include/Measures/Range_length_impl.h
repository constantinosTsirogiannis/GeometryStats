#ifndef RANGE_LENGTH_IMPL_H
#define RANGE_LENGTH_IMPL_H

#include<thread>
#include<cmath>

namespace FunctionalMeasures
{
  template<class KERNEL_TYPE>
  void Range_length<KERNEL_TYPE>::compute_all_hypergeometric_probabilities_s_minus_one( int sample_size)
  {
    _last_sample_size_exp = sample_size;
    int n = _vals.size();

    if( !_hypergeom_minus_one.empty() )
      _hypergeom_minus_one.clear();

    std::vector<Number_type> tempgeom;
    tempgeom.push_back(Number_type(sample_size)/Number_type(n-sample_size+1));

    for( int i= n-1; i>=_last_sample_size_exp-1; i-- )
    {
      Number_type x(i+1);
      tempgeom.push_back( tempgeom.back()/ Number_type(x/(x+Number_type(1.0)-Number_type(_last_sample_size_exp)))  );
    }

    for( int i= tempgeom.size()-1; i>=0; i-- )
      _hypergeom_minus_one.push_back(tempgeom[i]);

  } // compute_all_hypergeometric_probabilities_s_minus_one(...)

  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Number_type 
  Range_length<KERNEL_TYPE>::hypergeom_minus_one( int x )
  {
     if( x < _last_sample_size_exp-1 || x > _vals.size() )
       return Number_type(0.0);

     if( x == _vals.size() )
       return Number_type(_last_sample_size_exp)/Number_type(_vals.size()-_last_sample_size_exp+1);

     return _hypergeom_minus_one[x-_last_sample_size_exp+1];
  }


  // Computes all together the probability values f(x) = \binom{x}{s-2}/\binom{n}{s}

  template<class KERNEL_TYPE>
  void Range_length<KERNEL_TYPE>::
  compute_all_hypergeometric_probabilities_s_minus_two( int sample_size)
  {
    _last_sample_size_var = sample_size;
    int n = _vals.size();

    if( !_hypergeom_minus_two.empty() )
      _hypergeom_minus_two.clear();

    std::vector<Number_type> tempgeom;
    tempgeom.push_back(Number_type(sample_size)/Number_type(n-sample_size+1));

    tempgeom.back() = tempgeom.back()*Number_type(sample_size-1)/Number_type(n-sample_size+2);

    for( int i= n-1; i>=_last_sample_size_var-2; i-- )
    {
      Number_type x(i+1);
      tempgeom.push_back( tempgeom.back()/ Number_type(x/(x+Number_type(2.0)-Number_type(_last_sample_size_var)))  );
    }

    for( int i= tempgeom.size()-1; i>=0; i-- )
      _hypergeom_minus_two.push_back( tempgeom[i] );

  } // compute_all_hypergeometric_probabilities_s_minus_two(...)

  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Number_type 
  Range_length<KERNEL_TYPE>::hypergeom_minus_two( int x )
  {
     if( x < _last_sample_size_var-2 || x > _vals.size() )
       return Number_type(0.0);

     return _hypergeom_minus_two[x-_last_sample_size_var+2];
  }

  template<class KERNEL_TYPE>
  template<class RangeIterator>
  void Range_length<KERNEL_TYPE>::load_values( RangeIterator r_begin, RangeIterator r_end)
  {
    assert(Base::probability_distribution() == Kernel::UNIFORM_FIXED_SIZE);

    _vals.clear();
    _vals.insert(_vals.begin(), r_begin, r_end);

    this->mapped_indices.clear();

    for(int i=0; i<_vals.size(); i++)
      this->mapped_indices.push_back(i);

    initialise_measure();

    execute_precomputations();

  } // load_values( RangeIterator r_begin, RangeIterator r_end)

  template<class KERNEL_TYPE>
  template<class RangeIterator>
  void Range_length<KERNEL_TYPE>::load_values( RangeIterator r_begin_a, RangeIterator r_end_a,
                                                RangeIterator r_begin_b, RangeIterator r_end_b)
  {
    assert(Base::probability_distribution() == Kernel::POISSON_BINOMIAL || 
           Base::probability_distribution() == Kernel::POISSON_BINOMIAL_FIXED_SIZE);

    _vals.clear();
    _vals.insert(_vals.begin(), r_begin_a, r_end_a);
  
    insert_probability_values(r_begin_b, r_end_b);

    this->mapped_indices.clear();

    for(int i=0; i<_vals.size(); i++)
      this->mapped_indices.push_back(i);

    initialise_measure();

    execute_precomputations();

  } // load_values( RangeIterator r_begin, RangeIterator r_end)


  template<class KERNEL_TYPE>
  void Range_length<KERNEL_TYPE>::_execute_precomputations_uniform_fixed_size()
  { 
    std::set<Number_type, Is_smaller> srt_vals;
    std::map<Number_type, int, Is_smaller> indices;

    assert(this->prb.size() == 0);
 
    for(int i=0; i<_vals.size(); i++)
    {
      srt_vals.insert(_vals[i]); 
      indices[_vals[i]] = i;
    }
 
    typename std::set<Number_type, Is_smaller>::iterator it;

    int count=0;

    for(it=srt_vals.begin(); it != srt_vals.end(); it++)
    {
      _vals[count] = *it;
      this->mapped_indices[indices[*it]] = count;
      count++;
    }

  } // void _execute_precomputations_uniform_fixed_size()


  template<class KERNEL_TYPE>
  void Range_length<KERNEL_TYPE>::_execute_precomputations_poisson_binomial()
  { 
    std::map<Number_type, Number_type, Is_smaller> paired_vals;
    std::map<Number_type, int, Is_smaller> indices;

    assert(_vals.size() == this->prb.size());
 
    for(int i=0; i<_vals.size(); i++)
    {
      paired_vals[_vals[i]] = this->prb[i]; 
      indices[_vals[i]] = i;
    }
 
    typename std::map<Number_type, Number_type, Is_smaller>::iterator it;

    int count=0;

    for(it=paired_vals.begin(); it != paired_vals.end(); it++)
    {
      _vals[count] = it->first;
      this->prb[count] = it->second;
      this->mapped_indices[indices[it->first]] = count;
      count++;
    }

  } // void _execute_precomputations_poisson_binomial()

  template<class KERNEL_TYPE>
  template<class RangeIterator>
  typename KERNEL_TYPE::Number_type 
  Range_length<KERNEL_TYPE>::operator()(RangeIterator r_begin, RangeIterator r_end)
  {
    if(r_end-r_begin < 2)
      return Number_type(0.0);

    Number_type min = _vals[this->mapped_indices[*r_begin]],
                max = _vals[this->mapped_indices[*r_begin]];
 
    for(RangeIterator rit = r_begin; rit != r_end; rit++)
    {
      if(_vals[this->mapped_indices[*rit]] < min)
        min = _vals[this->mapped_indices[*rit]]; 

      if(_vals[this->mapped_indices[*rit]] > max)
        max = _vals[this->mapped_indices[*rit]];
    }

    return max-min;

  } // Number_type operator()(RangeIterator r_begin, RangeIterator r_end)


  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Number_type 
  Range_length<KERNEL_TYPE>::_compute_expectation_uniform_fixed_size(int sample_size)
  {
    if(sample_size <= 1 || sample_size > _vals.size())
      return Number_type(0.0);

    if(_last_sample_size_var == sample_size)
      return _exp;

    _last_sample_size_exp = sample_size;

    compute_all_hypergeometric_probabilities_s_minus_one(sample_size);

    Number_type total(0.0);
    int size = _vals.size();

    for(int i=0; i<_vals.size(); i++)
      total += _vals[i]*(hypergeom_minus_one(i) - hypergeom_minus_one(size-i-1)); 

    _exp = total;

    return _exp;

  } // _compute_expectation_uniform_fixed_size(int sample_size)


  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Number_type 
  Range_length<KERNEL_TYPE>::_compute_variance_uniform_fixed_size(int sample_size)
  {
    if(sample_size <= 1 || sample_size > _vals.size()-1)  // > minus one the size because then there is only one subset.
      return Number_type(0.0);

    if(_last_sample_size_exp != sample_size)
      _compute_expectation_uniform_fixed_size(sample_size);

    if(_last_sample_size_var == sample_size)
      return _var;

    _last_sample_size_var = sample_size;
    compute_all_hypergeometric_probabilities_s_minus_two(sample_size);

    Number_type total(0.0);
    int size = _vals.size();

    for(int i=0; i<_vals.size(); i++)
      total += _vals[i]*_vals[i]*(hypergeom_minus_one(i) + hypergeom_minus_one(size-i-1)); 

    Polynomial mins,maxes,res_p;

    mins.assign(_vals.size(), Protected_number_type(Number_type(0.0)));
    maxes.assign(_vals.size(), Protected_number_type(Number_type(0.0)));

    for(int i=0; i<_vals.size(); i++)
    {
      mins[i] = Protected_number_type(_vals[i]);
      maxes[_vals.size()-1-i] = Protected_number_type(_vals[i]);
    }

    Polynomial_multiplication mult;
 
    mult(mins,maxes,res_p);

    for(int i=0; i<_vals.size()-1; i++)
      total = total - (Number_type(2.0)*(res_p[i].to_number_type())*hypergeom_minus_two(size-i-2)); 

    _var = total - (_exp*_exp);

    return _var;

  } // _compute_variance_uniform_fixed_size(int sample_size)


  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Number_type 
  Range_length<KERNEL_TYPE>::_compute_expectation_poisson_binomial()
  {
    Protected_number_type total(Number_type(0.0)),
                          nprbs_before(Number_type(1.0)),
                          nprbs_after(Number_type(1.0));

    if(_exp != Number_type(-1.0))
      return _exp;

    for(int i=0; i<_vals.size(); i++)
      nprbs_after = nprbs_after*(Protected_number_type(1.0)-Protected_number_type(this->prb[i]));
  
    for(int i=0; i<_vals.size(); i++)
    {
      nprbs_after = nprbs_after/(Protected_number_type(1.0)-Protected_number_type(this->prb[i]));

      total += (Protected_number_type(_vals[i]*this->prb[i])*nprbs_after) - 
               (Protected_number_type(_vals[i]*this->prb[i])*nprbs_before); 

      nprbs_before = nprbs_before*(Protected_number_type(1.0)-Protected_number_type(this->prb[i]));
    }

    _exp = total.to_number_type();

    return _exp;

  } // _compute_expectation_poisson_binomial()


  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Number_type 
  Range_length<KERNEL_TYPE>::_compute_variance_poisson_binomial()
  {
    if(_var != Number_type(-1.0))
      return _var; 

    if(_vals.size() < 2)
      return Number_type(0.0);

    Protected_number_type total(Number_type(0.0)),
                          nprbs_before(Number_type(1.0)),
                          nprbs_after(Number_type(1.0)), 
                          none_prob, temp(1.0);

    std::vector<Protected_number_type> sum_before;

    for(int i=0; i<_vals.size(); i++)
    {
      if(i>0)
        nprbs_after = nprbs_after*(Protected_number_type(1.0)-Protected_number_type(this->prb[i-1]));

      sum_before.push_back(nprbs_after);
    }  

    //for(int i=0; i<_vals.size(); i++)
    //  std::cout << "{A} sum_before [" << i << "] :" << sum_before[i].to_number_type() << std::endl; 

    sum_before[0] = sum_before[0]*Protected_number_type(_vals[0]*this->prb[0]);

    for(int i=1; i<_vals.size(); i++)
      sum_before[i] = sum_before[i-1]+(sum_before[i]*Protected_number_type(_vals[i]*this->prb[i]));

    //for(int i=0; i<_vals.size(); i++)
    //  std::cout << "{B} sum_before [" << i << "] :" << sum_before[i].to_number_type() << std::endl;

    nprbs_after = nprbs_after*(Protected_number_type(1.0)-Protected_number_type(this->prb.back()));

    none_prob = nprbs_after;

    // The first elemment can only be a max if it is also a min, resulting in a zero-length
    // interval. Therefore, we just ingore this case in the following loop.

    //nprbs_after = nprbs_after/(Protected_number_type(1.0)-Protected_number_type(prb[0]));
    //nprbs_before = nprbs_before*(Protected_number_type(1.0)-Protected_number_type(prb[0]));

    for(int i=0; i<_vals.size(); i++)
    {
      nprbs_after = nprbs_after/(Protected_number_type(1.0)-Protected_number_type(this->prb[i]));

      //std::cout << " nprbs_after[" << i << "]: " << (nprbs_after.to_number_type()) 
      //          << " , nprbs_before[" << i << "] : " << (nprbs_before.to_number_type()) << std::endl;

      Protected_number_type ttt(Protected_number_type(_vals[i]*_vals[i]*this->prb[i])*(nprbs_after + nprbs_before));

      total = total + ttt;
      //std::cout << " Total squares["<< i <<"]: " << (total.to_number_type()) 
      //          << " , ttt: " << (ttt.to_number_type()) << std::endl; 

      if(i>0)       
        total = total + Protected_number_type(Number_type(-2.0*_vals[i]*this->prb[i]))*nprbs_after*sum_before[i-1];  

      // Account for the case where the same value is both the min and max of the range
      // (The range length is zero in these cases, yet the calculations that we have done
      //  so far have included in the variance a non-zero term which needs to be canceled out)
      
      total = total - (Protected_number_type(Number_type(2.0)*_vals[i]*_vals[i]*this->prb[i])*
                       none_prob/(Protected_number_type(1.0)-Protected_number_type(this->prb[i])));

      nprbs_before = nprbs_before*(Protected_number_type(1.0)-Protected_number_type(this->prb[i]));

    } // for(int i=0; i<_vals.size(); i++)

    if(_exp == Number_type(-1.0))
      _compute_expectation_poisson_binomial();

    _var = total.to_number_type();
    _var = _var-(_exp*_exp);

    //std::cout << " _exp: " << _exp << std::endl;

    return _var;

  } // _compute_variance_poisson_binomial()

} // namespace FunctionalMeasures

#endif //RANGE_LENGTH_IMPL_H
