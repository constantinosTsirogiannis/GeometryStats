#ifndef BOUNDING_BOX_VOLUME_2D_IMPL_H
#define BOUNDING_BOX_VOLUME_2D_IMPL_H

#include<vector>
#include<algorithm>

  template<class KERNEL_TYPE>
  template<class RangeIterator>
  typename KERNEL_TYPE::Number_type 
  FunctionalMeasures::Bounding_box_volume_2D<KERNEL_TYPE>::
  operator()(RangeIterator r_begin, RangeIterator r_end)
  {
    if(r_end-r_begin<2 || _dim < 1)
      return Number_type(0.0);

    Number_type vol(1.0);

    for(int d=0; d<2; d++)
    {
      Number_type min = _points[*r_begin][d],
                  max = _points[*r_begin][d];
 
      for(RangeIterator rit = r_begin; rit != r_end; rit++)
      {
        if(_points[*rit][d] < min)
          min = _points[*rit][d]; 

        if(_points[*rit][d] > max)
          max = _points[*rit][d];
      }

      vol = vol*(max-min);

    } // for(int d=0; d<_dim; d++)

    return vol;

  } // Number_type operator()(RangeIterator r_begin, RangeIterator r_end)


  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Number_type 
  FunctionalMeasures::Bounding_box_volume_2D<KERNEL_TYPE>::_compute_all_max_products()
  {
    Protected_number_type res(0.0);

    if(_points.size() < 2)
      return Number_type(0.0);

    Protected_number_type tmp_res_1, tmp_res_2;

    res += tmp_res_1 = _compute_max_products_single_point();  
    res += tmp_res_2 = _compute_max_products_two_distinct_points();

    //std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    //std::cout << " Res 1: " << tmp_res_1.n() << " , " << tmp_res_1.exp() << std::endl;
    //std::cout << " Res 2: " << tmp_res_2.n() << " , " << tmp_res_2.exp() << std::endl;
    //std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

    for(int i=0; i<_points.size(); i++)
      _points[i][0] = Number_type(-1.0)*_points[i][0];

    res += tmp_res_1 = _compute_max_products_single_point(); 
    res += tmp_res_2 = _compute_max_products_two_distinct_points();

    //std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    //std::cout << " Res 3: " << tmp_res_1.n() << " , " << tmp_res_1.exp() << std::endl;
    //std::cout << " Res 4: " << tmp_res_2.n() << " , " << tmp_res_2.exp() << std::endl;
    //std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

    for(int i=0; i<_points.size(); i++)
      _points[i][1] = Number_type(-1.0)*_points[i][1];   

    res += tmp_res_1 = _compute_max_products_single_point(); 
    res += tmp_res_2 = _compute_max_products_two_distinct_points();

    //std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    //std::cout << " Res 5: " << tmp_res_1.n() << " , " << tmp_res_1.exp() << std::endl;
    //std::cout << " Res 6: " << tmp_res_2.n() << " , " << tmp_res_2.exp() << std::endl;
    //std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

    for(int i=0; i<_points.size(); i++)
      _points[i][0] = Number_type(-1.0)*_points[i][0];

    res += tmp_res_1 = _compute_max_products_single_point(); 
    res += tmp_res_2 = _compute_max_products_two_distinct_points();

    //std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;
    //std::cout << " Res 7: " << tmp_res_1.n() << " , " << tmp_res_1.exp() << std::endl;
    //std::cout << " Res 8: " << tmp_res_2.n() << " , " << tmp_res_2.exp() << std::endl;
    //std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;

    for(int i=0; i<_points.size(); i++)
      _points[i][1] = Number_type(-1.0)*_points[i][1]; 

    //std::cout << " Res: " << res.n() << " , " << res.exp() << std::endl;

    _exp = res.to_number_type();

    return _exp;

  } // _compute_all_max_products()

  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Numeric_traits::Protected_number_type 
  FunctionalMeasures::Bounding_box_volume_2D<KERNEL_TYPE>::
  _compute_max_products_two_distinct_points()
  {
 //std::cout << " $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ " << std::endl;

    std::vector<Number_type> sorted_y_coords, sorted_y_probs;
    std::vector<Protected_number_type> sorted_y_prods;

    sorted_y_coords.assign(_points.size(),Number_type(0.0));
    sorted_y_prods.assign(_points.size(),Protected_number_type(0.0));
    sorted_y_probs.assign(_points.size(),Number_type(0.0));

    for(int i=0; i<_points.size(); i++)
    {
      sorted_y_coords[_map_index_to_rank_y[i]]=_points[i][1];
      sorted_y_probs[_map_index_to_rank_y[i]]=this->prb[i];
      sorted_y_prods[_map_index_to_rank_y[i]]=_piys[_map_index_to_rank_y[i]]*
                                              Protected_number_type(this->prb[i]*_points[i][1]);
    }
    
    Product_tree tree;

    tree.construct_tree(sorted_y_coords,sorted_y_probs,sorted_y_prods);

    Protected_number_type res(0.0);

    for(int i=0; i<_points.size(); i++)
    {
      int index = _map_rank_x_to_index[i];
      int rank = _map_index_to_rank_y[index];
      Number_type y = _points[index][1];
      Protected_number_type val = _pixs[i]*Protected_number_type(this->prb[index]*_points[index][0]);

      //std::cout << " #################### " << std::endl;

      //std::cout << " Rank: " << i << std::endl;
      //std::cout << " Index: " << index << std::endl;
      //std::cout << " x-coord: " << _points[index][0] << std::endl;
      //std::cout << " Prb: " << this->prb[index] << std::endl;            
      //std::cout << " Pixs: " << _pixs[i].n() << " , " << _pixs[i].exp() << std::endl;      

      //std::cout << " Val: " << val.n() << " , " << val.exp() << std::endl;

      Protected_number_type query_res = tree.product_elements_larger_than_y(y);
      tree.add_mark(rank);

      //std::cout << " Query res: " << query_res.n() << " , " << query_res.exp() << std::endl;

      Protected_number_type bv = val*query_res;

      //std::cout << " bv: " << bv.n() << " , " << bv.exp() << std::endl;

      res = res+bv;

      //std::cout << " Current res: " << res.n() << " , " << res.exp() << std::endl;

      //std::cout << " #################### " << std::endl;

    } // for(int i=0; i<_points.size(); i++)

    return res;

  } // _compute_max_products_two_distinct_points()

  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Numeric_traits::Protected_number_type 
  FunctionalMeasures::Bounding_box_volume_2D<KERNEL_TYPE>::_compute_max_products_single_point()
  {
    //return _compute_max_products_single_point_slow();

    return _compute_max_products_single_point_fast();
  }

  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Numeric_traits::Protected_number_type 
  FunctionalMeasures::Bounding_box_volume_2D<KERNEL_TYPE>::_compute_max_products_single_point_fast()
  {
    std::vector<Extended_point> points_x, points_y;
    std::vector< Protected_number_type> ref_vec;    

    _pixs.clear();
    _piys.clear();
    _map_index_to_rank_y.clear();
    _map_rank_x_to_index.clear();

    _pixs.assign(_points.size(), Protected_number_type(0.0));
    _piys.assign(_points.size(), Protected_number_type(0.0));
    _map_index_to_rank_y.assign(_points.size(),0);
    _map_rank_x_to_index.assign(_points.size(),0);

    for(int i=0; i<_points.size(); i++)
    {
      Extended_point ep;

      ep.pt = _points[i];
      ep.index = i;
      ep.prod = (Protected_number_type(1.0));

      points_x.push_back(ep);
      points_y.push_back(ep);

      ref_vec.push_back(Protected_number_type(this->prb[i]*_points[i][0]*_points[i][1]));

    } // for(int i=0; i<_points.size(); i++)

    sort(points_x.begin(), points_x.end(), Has_smaller_x());    
    sort(points_y.begin(), points_y.end(), Has_smaller_y());    

    points_x.back().prod = Protected_number_type(1.0);
    points_y.back().prod = Protected_number_type(1.0);
    _pixs.back() = Protected_number_type(1.0);
    _piys.back() = Protected_number_type(1.0);

    _map_index_to_rank_y[points_y.back().index]=points_y.size()-1;
    _map_rank_x_to_index[points_x.size()-1]=points_x.back().index;

    for(int i=_points.size()-2; i>=0; i--)
    {
      points_x[i].prod = points_x[i+1].prod*Protected_number_type(1.0-this->prb[points_x[i+1].index]);
      points_y[i].prod = points_y[i+1].prod*Protected_number_type(1.0-this->prb[points_y[i+1].index]);
  
      ref_vec[points_x[i].index] = ref_vec[points_x[i].index]*points_x[i].prod;
      ref_vec[points_y[i].index] = ref_vec[points_y[i].index]*points_y[i].prod;


      _pixs[i] = points_x[i].prod;
      _piys[i] = points_y[i].prod;
      _map_index_to_rank_y[points_y[i].index]=i;
      _map_rank_x_to_index[i]=points_x[i].index;
    }  

    Successor_product_structure_fast sps;

    std::vector<Number_type> coords, probs;

    for(int i=0; i<points_y.size(); i++)
    {
      coords.push_back(points_y[i].pt[1]);
      probs.push_back(this->prb[points_y[i].index]);
    }  

    sps.construct_structure(coords,probs);

    for(int i=_points.size()-1; i>=0; i--)
    {
      ref_vec[points_x[i].index] = ref_vec[points_x[i].index]*
                                   sps.product_elements_larger_than_y(points_x[i].pt[1]);

      sps.add_mark(_map_index_to_rank_y[points_x[i].index]);
    }

    Protected_number_type res(0.0);

    for(int i=_points.size()-1; i>=0; i--)
      res += ref_vec[i];

    return res;

  } // _compute_max_products_single_point_fast()

  template<class KERNEL_TYPE>
  typename KERNEL_TYPE::Numeric_traits::Protected_number_type 
  FunctionalMeasures::Bounding_box_volume_2D<KERNEL_TYPE>::_compute_max_products_single_point_slow()
  {
    std::vector<Extended_point> points_x, points_y;
    std::vector< Protected_number_type> ref_vec;    

    _pixs.clear();
    _piys.clear();
    _map_index_to_rank_y.clear();
    _map_rank_x_to_index.clear();

    _pixs.assign(_points.size(), Protected_number_type(0.0));
    _piys.assign(_points.size(), Protected_number_type(0.0));
    _map_index_to_rank_y.assign(_points.size(),0);
    _map_rank_x_to_index.assign(_points.size(),0);

    for(int i=0; i<_points.size(); i++)
    {
      Extended_point ep;

      ep.pt = _points[i];
      ep.index = i;
      ep.prod = (Protected_number_type(1.0));

      points_x.push_back(ep);
      points_y.push_back(ep);

      ref_vec.push_back(Protected_number_type(this->prb[i]*_points[i][0]*_points[i][1]));

    } // for(int i=0; i<_points.size(); i++)

    sort(points_x.begin(), points_x.end(), Has_smaller_x());    
    sort(points_y.begin(), points_y.end(), Has_smaller_y());    

    points_x.back().prod = Protected_number_type(1.0);
    points_y.back().prod = Protected_number_type(1.0);
    _pixs.back() = Protected_number_type(1.0);
    _piys.back() = Protected_number_type(1.0);

    //std::cout << " Points_x back index: " << points_x.back().index << std::endl;
    //std::cout << " Points_x back coord: " << _points[points_x.back().index][0] << std::endl;

    _map_index_to_rank_y[points_y.back().index]=points_y.size()-1;
    _map_rank_x_to_index[points_x.size()-1]=points_x.back().index;

    for(int i=_points.size()-2; i>=0; i--)
    {
      points_x[i].prod = points_x[i+1].prod*Protected_number_type(1.0-this->prb[points_x[i+1].index]);
      points_y[i].prod = points_y[i+1].prod*Protected_number_type(1.0-this->prb[points_y[i+1].index]);
  
      ref_vec[points_x[i].index] = ref_vec[points_x[i].index]*points_x[i].prod;
      ref_vec[points_y[i].index] = ref_vec[points_y[i].index]*points_y[i].prod;


      _pixs[i] = points_x[i].prod;
      _piys[i] = points_y[i].prod;
      _map_index_to_rank_y[points_y[i].index]=i;
      _map_rank_x_to_index[i]=points_x[i].index;

      //std::cout << " Pix #" << i << " :" << _pixs[i].n() << " , " << _pixs[i].exp() << std::endl; 
    }  

    Successor_product_structure sps;

    for(int i=_points.size()-1; i>=0; i--)
    {
      //std::cout << " Element with x-coord: " << _points[points_x[i].index] 
      //          << " had before-query value: " << ref_vec[points_x[i].index].n() << " , " 
      //                                         << ref_vec[points_x[i].index].exp() << std::endl;

      ref_vec[points_x[i].index] = ref_vec[points_x[i].index]*
                                   sps.return_successor_product(points_x[i].pt[1]).first;

      //std::cout << " Get after-query value: " << ref_vec[points_x[i].index].n() << " , " 
      //                                        << ref_vec[points_x[i].index].exp() << std::endl;

      sps.insert(points_x[i].pt[1], this->prb[points_x[i].index]);
    }

    Protected_number_type res(0.0);

    for(int i=_points.size()-1; i>=0; i--)
      res += ref_vec[i];

    //for(int i=0; i<_points.size(); i++)
    //{
      //std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ "<< std::endl;
      //std::cout << " x-coord: " << _points[i][0] << std::endl;
      //std::cout << " Single point maximum: " << ref_vec[i].n() << " , " << ref_vec[i].exp() << std::endl; 
    //}


    return res;

  } // _compute_max_products_single_point_slow()

#endif // BOUNDING_BOX_VOLUME_2D_IMPL_H
