#ifndef MEAN_PAIRWISE_DISTANCE_EXACT_IMPL_H
#define MEAN_PAIRWISE_DISTANCE_EXACT_IMPL_H

  template<class KERNEL_TYPE>
  void FunctionalMeasures::Mean_pairwise_distance_exact<KERNEL_TYPE>::
        Pairwise_distance_functor_complete::operator()()
  {
    Compute_Euclidean_distance compute_distance;

    //int count = 0;

    if(_forward == true)
      for(int i = 0; i < _length; i++)
      {
        int k=(i*_interval)+_start; 

        Point p = (*_ppoints)[k];
 
        for(int j = 0; j < k; j++)
        {
          //count++;

          Number_type res = compute_distance(p,(*_ppoints)[j]);

          (*_vals)[_index].first += res;
          (*_vals)[_index].second += res*res;
          (*_ppoints_dist_sum)[k] += res;
        
        } // for(int j = 0; j < k; j++)

      } // for(int i = 0; i < _length; i++)
    else // of if(_forward == true)
    {
      int psize = _ppoints->size()-1;

      for(int i = 0; i < _length; i++)
      {
        int k=(i*_interval)+_start; 

        Point p = (*_ppoints)[k];
 
        for(int j = psize; j > k; j--) 
        {
          //count++;
          (*_ppoints_dist_sum)[k] += compute_distance(p,(*_ppoints)[j]);
        }        

      } // for(int j = psize; j > k; j--) 

    } // else of if(_forward == true)

  } // Pairwise_distance_functor_complete::operator()(...)


  template<class KERNEL_TYPE>
  void FunctionalMeasures::Mean_pairwise_distance_exact<KERNEL_TYPE>::
        Pairwise_distance_functor_complete_2::operator()()
  {
    Compute_Euclidean_distance compute_distance;

    //int count = 0;

    //std::cout << " First: " << _first << " , Last: " << _last << std::endl; 

    if(_forward == true)
      for(int i = _first; i <= _last; i++)
      {
        Point p = (*_ppoints)[i];
 
        for(int j = 0; j < i; j++)
        {
          //count++;

          Number_type res = compute_distance(p,(*_ppoints)[j]);

          (*_vals)[_index].first += res;
          (*_vals)[_index].second += res*res;
          (*_ppoints_dist_sum)[i] += res;
        
        } // for(int j = 0; j < i; j++)

      } // for(int i = 0; i < _length; i++)
    else // of if(_forward == true)
    {
      //std::cout << " Entered! " << std::endl;  

      int psize = _ppoints->size()-1;

      for(int i = _first; i <= _last; i++)
      {
        //std::cout << " i: " <<  i << std::endl;

        Point p = (*_ppoints)[i];
 
        for(int j = psize; j > i; j--) 
        {
          //count++;
          (*_ppoints_dist_sum)[i] += compute_distance(p,(*_ppoints)[j]);
        }        

      } // for(int j = psize; j > k; j--) 

    } // else of if(_forward == true)

  } // Pairwise_distance_functor_complete_2::operator()(...)

  template<class KERNEL_TYPE>
  void FunctionalMeasures::Mean_pairwise_distance_exact<KERNEL_TYPE>::
  Pairwise_distance_functor_subset_complete::operator()()
  {
    Compute_Euclidean_distance compute_distance;

    for(int i = 0; i < _length; i++)
    {
      int k=(i*_interval)+_start; 

      Point p = (*_ppoints)[k];

      for(int j = 0; j < k; j++)
        (*_vals)[_index] += compute_distance(p,(*_ppoints)[j]);
    }

  } // Pairwise_distance_functor_subset_complete::operator()(...)


  template<class KERNEL_TYPE>
  template<class RangeIterator>
  typename KERNEL_TYPE::Number_type 
  FunctionalMeasures::Mean_pairwise_distance_exact<KERNEL_TYPE>::
  _operator_simple(RangeIterator r_begin, RangeIterator r_end)
  {
    Number_type val(0.0);
    int count=0;

    Compute_Euclidean_distance compute_distance;

    for(RangeIterator it_1 = r_begin; it_1 != r_end; it_1++)
    {
      for(RangeIterator it_2 = r_begin; it_2 != it_1; it_2++)
        val += compute_distance(_points[*it_1],_points[*it_2]);

      count++;
    }

    return val*Number_type(2)/(Number_type(count)*Number_type(count-1));

  } // _operator_simple(RangeIterator r_begin, RangeIterator r_end)


  template<class KERNEL_TYPE>
  template<class RangeIterator>
  typename KERNEL_TYPE::Number_type
  FunctionalMeasures::Mean_pairwise_distance_exact<KERNEL_TYPE>::
  _operator_parallel(RangeIterator r_begin, RangeIterator r_end, int number_of_processors)
  {
    if(number_of_processors <= 0)
      number_of_processors = std::thread::hardware_concurrency();

    std::vector<std::thread> threads;

    int count=0;

    _temp_points.clear();

    for(RangeIterator it = r_begin; it != r_end; it++)
    {
      _temp_points.push_back(_points[*it]);
      count++;
    } 

    Number_type frac = Number_type(2.0)/(Number_type(count)*Number_type(count-1));

    std::vector<Number_type > vals; 
    std::vector<std::pair<int, int> > ranges;  
    std::pair<Number_type, Number_type> spair;

    vals.assign(number_of_processors, Number_type(0.0)); 

    int quota = _temp_points.size()/(number_of_processors),
        extra = _temp_points.size()%(number_of_processors);

    for(int i=0; i<number_of_processors; i++)
    {
      int current_quota = quota;      

      if(i<extra)
        current_quota++;     

      ranges.push_back(std::make_pair(i,current_quota));
    }

    for(int i=0; i<ranges.size(); i+=1 )
    {
      Pairwise_distance_functor_subset_complete pdf(ranges[i].first, ranges[i].second, 
                                                    number_of_processors, &_temp_points, &vals, i);
      threads.push_back(std::thread(pdf));
    } 

    for(int i=0; i<number_of_processors; i++ )
      threads[i].join();

    threads.clear();

    Number_type total_sum(0.0);

    for(int i=0; i<number_of_processors; i++ )
      total_sum += vals[i];

    return frac*total_sum;

  } // operator_parallel(int number_of_processors)


  template<class KERNEL_TYPE>
  template<class RangeIterator>
  void FunctionalMeasures::Mean_pairwise_distance_exact<KERNEL_TYPE>::
  load_point_set( RangeIterator r_begin, RangeIterator r_end,
                  bool use_parallelisation, 
                  int number_of_processors)
  {
    _points.clear();
    _points.insert(_points.begin(), r_begin, r_end);

    _point_dist_sum.clear();
    _point_dist_sum.assign(_points.size(),Number_type(0.0));

    _total_dist_sum = Number_type(0.0);
    _total_dist_sum_sq = Number_type(0.0);
    _point_dist_sum_sq = Number_type(0.0);

    if(_points.size() > 0)
      _dim = _points[0].size();
    else
      _dim = -1;

    if(use_parallelisation == false || number_of_processors == 1)
      precompute_distances_simple();
    else
      precompute_distances_parallel(number_of_processors);

  } // load_point_set( RangeIterator r_begin, RangeIterator r_end) 

  template<class KERNEL_TYPE>
  void FunctionalMeasures::Mean_pairwise_distance_exact<KERNEL_TYPE>::
  precompute_distances_simple()
  {
    Compute_Euclidean_distance compute_distance;

    for(int i=0; i<_points.size(); i++)
      for(int j=0; j<i; j++)
      {
        Number_type res = compute_distance(_points[i],_points[j]);

        _total_dist_sum += res;
        _total_dist_sum_sq += res*res;
        _point_dist_sum[i] += res;
        _point_dist_sum[j] += res;
      }  

    _point_dist_sum_sq = Number_type(0.0);

    for(int i=0; i<_points.size(); i++)
      _point_dist_sum_sq += _point_dist_sum[i] * _point_dist_sum[i];

  } // precompute_distances_simple() 

  template<class KERNEL_TYPE>
  void FunctionalMeasures::Mean_pairwise_distance_exact<KERNEL_TYPE>::
  precompute_distances_parallel(int number_of_processors)
  {
    if(number_of_processors <= 0)
      number_of_processors = std::thread::hardware_concurrency();

    std::vector<std::thread> threads;

    std::vector<std::pair<Number_type, Number_type> > vals; 
    std::vector<std::pair<int, int> > ranges;  
    std::pair<Number_type, Number_type> spair;

    
    spair.first = Number_type(0.0);
    spair.second = Number_type(0.0);

    vals.assign(number_of_processors, spair); 

    int quota = _points.size()/(number_of_processors),
        extra = _points.size()%(number_of_processors);

    for(int i=0; i<number_of_processors; i++)
    {
      int current_quota = quota;      

      if(i<extra)
        current_quota++;     

      ranges.push_back(std::make_pair(i,current_quota));
    }

    for(int i=0; i<ranges.size(); i++ )
    {
      Pairwise_distance_functor_complete pdf(ranges[i].first, ranges[i].second, number_of_processors,
                                             &_points, &_point_dist_sum, &vals, i, true);
      threads.push_back(std::thread(pdf));
    } 

    for(int i=0; i<number_of_processors; i++ )
      threads[i].join();

    threads.clear();

    for(int i=0; i<ranges.size(); i++ )
    {
      Pairwise_distance_functor_complete pdf(ranges[i].first, ranges[i].second, number_of_processors,
                                             &_points, &_point_dist_sum, &vals, i, false);
      threads.push_back(std::thread(pdf));
    } 

    for(int i=0; i<number_of_processors; i++ )
      threads[i].join();

    threads.clear();


    /*////////////////////////////////////////////////////////////////////
    int prev=0;

    Number_type all =  Number_type(4*_points.size())*Number_type(_points.size()-1)/Number_type(number_of_processors),
                one =  Number_type(1.0);
    Square_root sq_root;
    To_double   t_double;

    for(int i=0; i<number_of_processors; i++)
    {
      Number_type rest;


      if(prev ==0) 
        rest = 0;
      else
        rest = Number_type(4*prev)*Number_type(prev-1);

      Number_type sol = (sq_root(one+all+rest)+one)/Number_type(2.0);

      int k(std::floor(t_double(sol)));

      if(i == 0)    
        ranges.push_back(std::make_pair(0,k));
      else if(i<number_of_processors-1)
        ranges.push_back(std::make_pair(prev+1,k));
      else
        ranges.push_back(std::make_pair(prev+1,_points.size()-1)); 

      //std::cout << " Range first: " << ranges.back().first << "  , last: " << ranges.back().second << std::endl;

      prev = k;
    }


    for(int i=0; i<ranges.size(); i++ )
    {
      Pairwise_distance_functor_complete_2 pdf(ranges[i].first, ranges[i].second,
                                               &_points, &_point_dist_sum, &vals, i, true);
      threads.push_back(std::thread(pdf));
    } 

    for(int i=0; i<number_of_processors; i++ )
      threads[i].join();

    threads.clear();

    for(int i=0; i<number_of_processors; i++)
    {
      ranges[i].first = _points.size()-1-ranges[i].first;
      ranges[i].second = _points.size()-1-ranges[i].second;

      int temp = ranges[i].second;
      ranges[i].second = ranges[i].first;
      ranges[i].first = temp;

      //std::cout << " Range first: " << ranges[i].first << "  , last: " << ranges[i].second << std::endl;
    }

    for(int i=0; i<ranges.size(); i++ )
    {
      Pairwise_distance_functor_complete_2 pdf(ranges[i].first, ranges[i].second,
                                               &_points, &_point_dist_sum, &vals, i, false);
      threads.push_back(std::thread(pdf));
    } 

    for(int i=0; i<number_of_processors; i++ )
      threads[i].join();

    threads.clear();    
    *////////////////////////////////////////////////////////////////////



    for(int i=0; i<number_of_processors; i++ )
    {
      _total_dist_sum += vals[i].first;
      _total_dist_sum_sq += vals[i].second;
    }  

    _point_dist_sum_sq = Number_type(0.0);

    for(int i=0; i<_points.size(); i++)
      _point_dist_sum_sq += _point_dist_sum[i]*_point_dist_sum[i];
    
  } // precompute_distances_parallel(int number_of_processors)

  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Number_type 
  FunctionalMeasures::Mean_pairwise_distance_exact<KERNEL_TYPE>::
  compute_expectation(int sample_size)
  {
    if(sample_size < 2 || sample_size > _points.size())
      return Number_type(0.0);

    return _total_dist_sum*Number_type(2.0)/(Number_type(_points.size())*Number_type(_points.size()-1));
  }

  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Number_type 
  FunctionalMeasures::Mean_pairwise_distance_exact<KERNEL_TYPE>::
  compute_variance(int sample_size)
  {
    if(sample_size < 2 || sample_size > _points.size())
      return Number_type(0.0);

    Number_type fact = Number_type(4.0)/(Number_type(sample_size)*Number_type(sample_size)
                                         *Number_type(sample_size-1)*Number_type(sample_size-1)),
                coefficient_2 = Number_type(sample_size)*Number_type(sample_size-1)/
                              (Number_type(_points.size())*Number_type(_points.size()-1)),
                coefficient_3 = Number_type(sample_size)*Number_type(sample_size-1)*Number_type(sample_size-2)/
                              (Number_type(_points.size())*Number_type(_points.size()-1)*Number_type(_points.size()-2)),
                coefficient_4 = Number_type(sample_size)*Number_type(sample_size-1)
                                *Number_type(sample_size-2)*Number_type(sample_size-3)/
                              (Number_type(_points.size())*Number_type(_points.size()-1)
                              *Number_type(_points.size()-2)*Number_type(_points.size()-3));
    
    Number_type val(0.0);

    val += (_total_dist_sum*_total_dist_sum*coefficient_4);
    val += (_total_dist_sum_sq*(coefficient_2-(Number_type(2.0)*coefficient_3)+coefficient_4)); // plus coeff_4 because
    val += (_point_dist_sum_sq*(coefficient_3-coefficient_4)); // (in this line) we count the total_dist_sum_sq values. twice.  

    Number_type expec = compute_expectation(sample_size);

    return (fact*val) - (expec*expec);
  }

#endif //MEAN_PAIRWISE_DISTANCE_EXACT_IMPL_H
