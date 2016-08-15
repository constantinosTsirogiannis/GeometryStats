#ifndef SPLIT_TREE_H
#define SPLIT_TREE_H

#include<vector>
#include<map>
#include<cassert>

namespace FunctionalMeasures
{

template< class KERNEL_TYPE>
class Split_tree
{
 public:

  typedef KERNEL_TYPE                                   Kernel;
  typedef typename Kernel::Number_type                 Number_type;
  typedef typename Kernel::Numeric_traits              Numeric_traits;
  typedef typename Numeric_traits::Is_smaller          Is_smaller;
  typedef typename Kernel::Point                       Point;  
  typedef typename Kernel::Bbox_d                      Bbox_d;
  typedef typename Kernel::Compute_Euclidean_distance  Compute_Euclidean_distance;

 public:

  struct Pointset_pair
  { std::pair<int, int> set_a, set_b; }; 

 private:

  struct Node_list_element
  {
    Number_type coord;
    int p_id, node_index, lfrom, lto;

  }; // struct Node_list_element

  struct Node_lists
  {
    std::vector<std::vector<Node_list_element> > coord_lists;
    std::vector<std::vector<int> > p_indices;

    Node_lists():_n(-1){}

    void initialize_front_back(int d)
    { 
      _front.clear();
      _back.clear();
   
      _front.assign(d,-1);
      _back.assign(d,-1);
    }

    int number_of_active_elements() const
    { return _n;}

    int number_of_points() const
    { return coord_lists[0].size();}

    int dim() const
    { return coord_lists.size();}

    int front(int d) const
    { return _front[d];}

    int back(int d) const
    { return _back[d];}

    int set_front(int d, int f)
    { _front[d] = f;}

    int set_back(int d, int b)
    { _back[d] = b;}

    int set_number_of_active_elements(int n)
    { _n = n;}

    Number_type next_coord(int d, int i)
    { return coord_lists[d][coord_lists[d][i].lto].coord;}

    Number_type previous_coord(int d, int i)
    { return coord_lists[d][coord_lists[d][i].lfrom].coord;}

    Number_type next_element(int d, int i)
    { return coord_lists[d][i].lto;}

    Number_type previous_element(int d, int i)
    { return coord_lists[d][i].lfrom;}

    bool is_empty()
    { return _n == 0;}

    bool is_back(int dim, int index) const
    { return coord_lists[dim][index].lto == -1; }

    bool is_front(int dim, int index) const
    { return coord_lists[dim][index].lfrom == -1; }

    Number_type coord(int dim, int index) const
    { return coord_lists[dim][index].coord;}

    void delete_element(int p_id, int node_index)
    { 
      for(int i=0; i<this->dim(); i++)
      {
        coord_lists[i][p_indices[p_id][i]].node_index = node_index;
        Node_list_element nle = coord_lists[i][p_indices[p_id][i]];
                          
        if(nle.lfrom >= 0)
          coord_lists[i][nle.lfrom].lto = nle.lto;
        else
          _front[i] = nle.lto;

        if(nle.lto >= 0)
          coord_lists[i][nle.lto].lfrom = nle.lfrom;
        else
          _back[i] = nle.lfrom; 

      } // for(int i=0; i<this->dim; i++)

      _n--;

    } // delete_element(int p_id)

    void delete_element(int d, int i, int node_index)
    { delete_element(coord_lists[d][i].p_id, node_index);}

   private:

    int _n;
    std::vector<int> _front, _back;

  }; // struct Node_lists

 public:

  struct Node_type
  {
    Bbox_d bbox;
    Point p;
    int parent, right_child, left_child, p_start, p_end;
    Number_type sum_d_times_s, sum_dsquare_times_s;
   
    Node_type(): parent(-1), right_child(-1), left_child(-1), 
                 p_start(-1), p_end(-1), sum_d_times_s(0.0), 
                 sum_dsquare_times_s(0.0){}

    bool is_leaf() const
    { return (left_child==-1 && right_child==-1);}
 
  }; // struct Node_type

 public:

  Split_tree():_d(-1), _total_sums(0.0), _total_dist_sum(0.0),
               _interval_min(0.0), _interval_max(0.0){}

  Node_type node(int i)
  { return _nodes[i];}

  Point sorted_point(int i)
  { return _sorted_points[i];}

  int number_of_nodes()
  { return _nodes.size();}

  int number_of_sorted_points()
  { return _sorted_points.size();}

  template <class RangeIterator>
  void load_point_set( RangeIterator rbegin, RangeIterator rend)
  {
    _nodes.clear();
    _points.clear();
    _points.insert(_points.begin(), rbegin, rend);
    
    if(_points.size() > 0)
    {
      _d = _points[0].size();
      _construct_split_tree();
    }

  } // load_point_set( std::vector<Point> points)

  int dimension() const
  { return _d;}

  int root_index() const
  { return 0;}

  void _construct_split_tree()
  {
    Node_lists nl;

    _execute_precomputations(nl);
    _nodes.push_back(Node_type());
    _expand_split_tree(nl,_nodes.size()-1);
    _register_node_point_sets();
  }

  void construct_bounding_box(Node_lists &nl, Bbox_d &bbox)
  {
    bbox.set_dimension(nl.dim());  
  
    for(int i=0; i<nl.dim(); i++)
    {
       bbox.min(i) = nl.coord_lists[i][nl.front(i)].coord;
       bbox.max(i) = nl.coord_lists[i][nl.back(i)].coord;
    } 

  } // construct_bounding_box(Node_lists &nl, Bbox_d &bbox)

  std::pair<int,int> _in_order_traversal(int index)
  {
    if(_nodes[index].is_leaf())
    {
      _sorted_points.push_back(_nodes[index].p);

      int p_index = _sorted_points.size()-1;

      _nodes[index].p_start = p_index;
      _nodes[index].p_end = p_index;
      return std::make_pair(p_index,p_index);
    }

    std::pair<int,int> pp_l, pp_r;

    pp_l = _in_order_traversal(_nodes[index].left_child);
    pp_r = _in_order_traversal(_nodes[index].right_child);

    assert(pp_l.second == pp_r.first-1);

    _nodes[index].p_start = pp_l.first;
    _nodes[index].p_end = pp_r.second;

    return std::make_pair(pp_l.first,pp_r.second);

  } // std::pair<int,int> _in_order_traversal(int index)

  void _register_node_point_sets()
  { 
    if(_nodes.size() > 0)
    {
      std::pair<int,int> pp = _in_order_traversal(root_index());
      _nodes[root_index()].p_start = pp.first;
      _nodes[root_index()].p_end = pp.second;
    }
  }

  void _expand_split_tree( Node_lists &nl, int node_index)
  {
    if(nl.number_of_points() == 1)
    {
      Point p(nl.coord_lists.size());
      Bbox_d bbox;

      for(int i=0; i<nl.coord_lists.size(); i++)
        p(i) = nl.coord_lists[i][0].coord;

      assert(p.dim()>0);

      bbox.set_dimension(p.dim());

      for(int i=0; i<p.dim(); i++)
      {
        bbox.min(i) = p(i);
        bbox.max(i) = p(i);
      }

      _nodes[node_index].p = p;
      _nodes[node_index].bbox = bbox;

      return;
 
    } // if(nl.number_of_points() == 1)

    std::vector<int> node_indices; 

    int input_size = nl.number_of_active_elements(),
        current_size = nl.number_of_active_elements(),
        current_node_index = node_index;

    int  node_index_l = _nodes.size(),
         node_index_r = _nodes.size()+1;

    _nodes.push_back(Node_type());
    _nodes.push_back(Node_type());

    _nodes[current_node_index].left_child = node_index_l;
    _nodes[current_node_index].right_child = node_index_r;

    _nodes[node_index_l].parent = current_node_index;
    _nodes[node_index_r].parent = current_node_index;

    int node_count = 0;

    do // while(current_size > input_size/2)
    {
      Bbox_d bbox; 
      construct_bounding_box(nl, bbox);
      _nodes[current_node_index].bbox = bbox;

      std::pair<int, Number_type> dmax_pair = bbox.maximum_dimension();

      int dmax = dmax_pair.first;
      Number_type split_value = dmax_pair.second;

      int j = nl.front(dmax),
          k = nl.back(dmax),
          element_count = 1;

      // By definition neither j can be the last element, nor k can be the first.
      while( nl.next_coord(dmax,j) < split_value && nl.previous_coord(dmax,k) > split_value )
      {
        element_count++;
        j = nl.next_element(dmax,j);
        k = nl.previous_element(dmax,k);
      }

      if(nl.next_coord(dmax,j) >= split_value)
      {
        // Index j indicates the last element of the smallest set.
        // Delete this set by assigning it to a new node. 

        node_indices.push_back(node_index_l);
        current_node_index = node_index_r;

        if(nl.number_of_active_elements()-element_count>input_size/2)
        {
          node_index_l = _nodes.size();
          node_index_r = _nodes.size()+1;     

          _nodes.push_back(Node_type());
          _nodes.push_back(Node_type());

          _nodes[current_node_index].left_child = node_index_l;
          _nodes[current_node_index].right_child = node_index_r;
          _nodes[node_index_l].parent = current_node_index;
          _nodes[node_index_r].parent = current_node_index;
        }

        int next_j = nl.next_element(dmax,j);

        while(nl.front(dmax) != next_j)
        {
          int jj = nl.front(dmax);
          Node_list_element nle = nl.coord_lists[dmax][jj];
          nl.delete_element(dmax,jj,node_count);

        } // while(nl.front(dmax) != next_j)
        
      } // if(nl.is_back(dmax,j) || nl.next_coord(dmax,j) >= split_value)
      else
      {
        // Index k indicates the first element of the smallest set.
        // Delete this set by loading it to a new node. 

        node_indices.push_back(node_index_r);
        current_node_index = node_index_l;

        if(nl.number_of_active_elements()-element_count>input_size/2)
        {
          node_index_l = _nodes.size();
          node_index_r = _nodes.size()+1;     

          _nodes.push_back(Node_type());
          _nodes.push_back(Node_type());

          _nodes[current_node_index].left_child = node_index_l;
          _nodes[current_node_index].right_child = node_index_r;
          _nodes[node_index_l].parent = current_node_index;
          _nodes[node_index_r].parent = current_node_index;
        }

        int previous_k = nl.previous_element(dmax,k);

        while(nl.back(dmax) != previous_k)
        {
          int kk = nl.back(dmax);
          Node_list_element nle = nl.coord_lists[dmax][kk];
          nl.delete_element(dmax,kk,node_count);

        } // while(nl.back(dmax) != previous_k)

      } // else of if(nl.is_back(dmax,j) || nl.next_coord(dmax,j) >= split_value)

      node_count++;   

      current_size = nl.number_of_active_elements();

    }while(current_size > input_size/2);

    // Delete all the elements of the remaining node.
   
    node_indices.push_back(current_node_index);

    while(nl.is_empty() == false)
    {
      int kk = nl.back(0);
      Node_list_element nle = nl.coord_lists[0][kk];
      nl.delete_element(0,kk,node_count);

    } // while(nl.back(dmax) != k)

    // Distribute the elements of the node lists in the current node 
    // to the corresponding lists of the new nodes and recurse.
    
    std::vector<Node_lists> vec_nl;

    vec_nl.assign(node_indices.size(), Node_lists());

    for(int i=0; i < vec_nl.size(); i++)
    {
      vec_nl[i].coord_lists.assign(_d,std::vector<Node_list_element>());    
      vec_nl[i].initialize_front_back(_d);
    }    

    flush_to_new_node_lists(nl,vec_nl);

    for(int i=0; i < vec_nl.size(); i++)
      _expand_split_tree(vec_nl[i], node_indices[i]);

  } // _expand_split_tree( Node_lists &nl, int node_index)

  void flush_to_new_node_lists( Node_lists &nl,
                                 std::vector<Node_lists> &vec_nl )
  { 
    // First compute for each branch node the number 
    // of points that were assigned to this node

    std::vector<int> point_sizes;

    point_sizes.assign(vec_nl.size(),0);

    for(int i=0; i<nl.coord_lists[0].size(); i++)
    {
      Node_list_element nle = nl.coord_lists[0][i];
      
      assert(nle.node_index != -1);
      point_sizes[nle.node_index]++;
    } 
 
    for(int i=0; i<vec_nl.size(); i++)
      vec_nl[i].p_indices.assign(point_sizes[i],std::vector<int>());

    for(int i=0; i<vec_nl.size(); i++)
      for(int j=0; j<vec_nl[i].p_indices.size(); j++)
        vec_nl[i].p_indices[j].assign(_d,-1);

    // Push the elements from the old node lists
    // to the vector with the new node lists.

    for(int i=0; i<_d; i++)
    {
      for(int j=0; j<nl.coord_lists[i].size(); j++)
      {
        Node_list_element ole = nl.coord_lists[i][j];
        Node_list_element nle;

        nle.coord = ole.coord;
        nle.p_id = ole.p_id;
        nle.node_index = -1;

        if(vec_nl[ole.node_index].coord_lists[i].size() == 0)
          nle.lfrom = -1;
        else
          nle.lfrom = vec_nl[ole.node_index].coord_lists[i].size()-1;

        if(vec_nl[ole.node_index].coord_lists[i].size() == point_sizes[ole.node_index]-1)
          nle.lto = -1;
        else
          nle.lto = vec_nl[ole.node_index].coord_lists[i].size()+1;
          
        vec_nl[ole.node_index].coord_lists[i].push_back(nle);

      } // for(int j=0; j<nl.coord_lists[i].size(); j++)

    } // for(int i=0; i<_d; i++)

    for(int i=0; i<vec_nl.size(); i++)
    {
      vec_nl[i].set_number_of_active_elements(vec_nl[i].coord_lists[0].size());
 
      for(int j=0; j<_d; j++)
      {
        vec_nl[i].set_front(j,0);
        vec_nl[i].set_back(j,vec_nl[i].number_of_points()-1);
      }

    } // for(int i=0; i<vec_nl.size(); i++)

    std::vector<int> p_id_map;

    p_id_map.assign(nl.coord_lists[0].size(),-1); 

    // Last step: for every node list, fix the p_ids
    // so that their range extends from 0 to the node's 
    // total number of assigned points minus one.

    for(int i=0; i < vec_nl.size(); i++)
      for(int j=0; j<vec_nl[i].coord_lists[0].size(); j++)
        p_id_map[vec_nl[i].coord_lists[0][j].p_id] = j;

 
    for(int i=0; i < vec_nl.size(); i++)
      for(int j=0; j<vec_nl[i].coord_lists.size(); j++)
        for(int k=0; k<vec_nl[i].coord_lists[j].size(); k++)
        {
          int old_id = vec_nl[i].coord_lists[j][k].p_id;
          vec_nl[i].coord_lists[j][k].p_id = p_id_map[old_id];
          vec_nl[i].p_indices[p_id_map[old_id]][j] = k;
        }  

  } // flush_to_new_node_lists(...)

  void _execute_precomputations(Node_lists &nl)
  {
    nl.coord_lists.assign(_d,std::vector<Node_list_element>());
    nl.p_indices.assign(_points.size(), std::vector<int>());
    nl.initialize_front_back(_d);

    for(int i=0; i<nl.p_indices.size(); i++)
      nl.p_indices[i].assign(_d,-1);

    for(int i=0; i<_d; i++)
    {
      std::map<Number_type, std::vector<int>, Is_smaller> cmap;
      
      for(int j=0; j<_points.size(); j++)
      {
        if(cmap.find(_points[j][i]) == cmap.end())
          cmap[_points[j][i]] = std::vector<int>();

        cmap[_points[j][i]].push_back(j);  
      }

      typename std::map<Number_type, std::vector<int>, Is_smaller>::iterator it;

      for(it = cmap.begin(); it != cmap.end(); it++)
      {
        std::vector<int> vec = it->second;

        for(int j=0; j<vec.size(); j++)
        {
          Node_list_element nle;

          nle.p_id = vec[j];
          nle.coord = it->first;
          nle.node_index = -1;

          if(nl.coord_lists[i].size() == 0)
            nle.lfrom = -1;
          else
            nle.lfrom = nl.coord_lists[i].size()-1;

          if(nl.coord_lists[i].size() == _points.size()-1)
            nle.lto = -1;
          else
            nle.lto = nl.coord_lists[i].size()+1;
 
          nl.coord_lists[i].push_back(nle);
          nl.p_indices[vec[j]][i] = nl.coord_lists[i].size()-1;
        } 

      } // for(it = cmap.begin(); it != cmap.end(); it++)

    } // for(int i=0; i<_d; i++)
 
    nl.set_number_of_active_elements(nl.coord_lists[0].size());

    for(int i=0; i<_d; i++)
    {
      nl.set_front(i,0);
      nl.set_back(i,nl.number_of_points()-1);
    }

  } //  execute_precomputations(Node_lists &nl)

  int total_number_of_checks()
  { return _total_number_of_checks;}

  template <class OutputIterator> 
  void _extract_decomposition( Number_type &sf, int index, OutputIterator ot)
  {   
    if(_nodes[index].is_leaf() == false)
    {
      _extract_decomposition(sf, _nodes[index].left_child, _nodes[index].right_child, ot);
      _extract_decomposition(sf, _nodes[index].left_child, ot);
      _extract_decomposition(sf, _nodes[index].right_child, ot);
    }

  } // extract_decomposition( Number_type &sf, int index, OutputIterator ot)

  template <class OutputIterator_a, class OutputIterator_b> 
  void extract_decomposition( Number_type &sf, OutputIterator_a ot_a, OutputIterator_b ot_b)
  {  
    for(int i=0; i<_sorted_points.size(); i++)
      *ot_a++ = _sorted_points[i];
 
    int index=0;

    _total_number_of_checks = 0;
    _decomposition_size = 0;

    _extract_decomposition(sf, index, ot_b);

    //std::cout  << _d << " , " << _decomposition_size << " , "<<  _total_number_of_checks << " , " 
    //                                         << Number_type(2*_total_number_of_checks)/
    //                                            Number_type(_points.size()*(_points.size()-1)) <<  std::endl;

    flush_sums_to_leaves();

  } // extract_decomposition( Number_type &sf, OutputIterator_a ot_a, OutputIterator_b ot_b)
  

  template <class OutputIterator> 
  void _extract_decomposition_epsilon( Number_type &epsilon, int index, OutputIterator ot)
  {   
    if(_nodes[index].is_leaf() == false)
    {
      _extract_decomposition_epsilon(epsilon, _nodes[index].left_child, _nodes[index].right_child, ot);
      _extract_decomposition_epsilon(epsilon, _nodes[index].left_child, ot);
      _extract_decomposition_epsilon(epsilon, _nodes[index].right_child, ot);
    }

  } // extract_decomposition_epsilon( Number_type &epsilon, int index, OutputIterator ot)

  template <class OutputIterator_a, class OutputIterator_b> 
  void extract_decomposition_epsilon( Number_type &epsilon, OutputIterator_a ot_a, OutputIterator_b ot_b)
  {    
    for(int i=0; i<_sorted_points.size(); i++)
      *ot_a++ = _sorted_points[i];
 
    int index=0;

    _total_dist_sum = Number_type(0.0);
    _interval_min = Number_type(0.0);
    _interval_max = Number_type(0.0);

    _extract_decomposition_epsilon(epsilon, index, ot_b);

    flush_sums_to_leaves();
  }

  Number_type total_decomposition_cost()
  { return _total_dist_sum;}

  std::pair<Number_type, Number_type> cost_interval()
  { return std::make_pair(_interval_min,_interval_max);}

  void clear()
  {
    _nodes.clear();
    _points.clear();
    _sorted_points.clear();
    _d = -1;
    _total_sums = Number_type(0.0);
    _total_dist_sum = Number_type(0.0);
    _interval_min = Number_type(0.0);
    _interval_max = Number_type(0.0);
  }
  
  bool is_well_separated_pair(Number_type &sf, int index_a, int index_b)
  {
    if(_nodes[index_a].is_leaf() && _nodes[index_b].is_leaf())
      return true;

    Number_type  diameter_a = _nodes[index_a].bbox.diameter(),
                 diameter_b = _nodes[index_b].bbox.diameter(),
                 max_diameter;

    if(diameter_a > diameter_b)
      max_diameter = diameter_a;
    else
      max_diameter = diameter_b;

    Compute_Euclidean_distance compute_distance;

    Number_type dist = compute_distance(_nodes[index_a].bbox, _nodes[index_b].bbox);
    
    if(dist >= sf * max_diameter)
      return true;
    
    return false;

  } // bool is_well_separated_pair(Number_type &sf, int index_a, int index_b)

  template <class OutputIterator> 
  void _extract_decomposition( Number_type &sf, int index_a, int index_b, OutputIterator ot)
  {
     _total_number_of_checks++;

     if(is_well_separated_pair(sf,index_a,index_b)==true)
     {
       _decomposition_size++;

       Pointset_pair ppr;

       ppr.set_a.first = _nodes[index_a].p_start;
       ppr.set_a.second = _nodes[index_a].p_end;

       ppr.set_b.first = _nodes[index_b].p_start;
       ppr.set_b.second = _nodes[index_b].p_end;

       int size_a = ppr.set_a.second - ppr.set_a.first + 1,
           size_b = ppr.set_b.second - ppr.set_b.first + 1;

       Number_type dist = Compute_Euclidean_distance()(_sorted_points[ppr.set_a.first], 
                                                       _sorted_points[ppr.set_b.first]);

       _nodes[index_a].sum_d_times_s += dist*Number_type(size_b);
       _nodes[index_a].sum_dsquare_times_s += dist*dist*Number_type(size_b);
       _nodes[index_b].sum_d_times_s += dist*Number_type(size_a);
       _nodes[index_b].sum_dsquare_times_s += dist*dist*Number_type(size_a);

       *ot++ = ppr;
     }
     else
     {
       Number_type res_a = _nodes[index_a].bbox.maximum_dimension_length(),
                   res_b = _nodes[index_b].bbox.maximum_dimension_length();

       if(res_a<= res_b)
       {
         _extract_decomposition(sf, index_a, _nodes[index_b].left_child,  ot);
         _extract_decomposition(sf, index_a, _nodes[index_b].right_child, ot);
       }
       else
       {
         _extract_decomposition(sf, _nodes[index_a].left_child, index_b, ot);
         _extract_decomposition(sf, _nodes[index_a].right_child, index_b, ot);
       }

     } // else of if(is_well_separated(sf,index_a,index_b)==true) 

  } // _extract_decomposition(...)


  bool is_well_separated_pair_epsilon(Number_type &epsilon, int index_a, int index_b, 
                                       Number_type &dist, Number_type &interval_min, Number_type &interval_max)
  {
    Compute_Euclidean_distance compute_distance;

    if(_nodes[index_a].is_leaf() && _nodes[index_b].is_leaf())
    {
      dist = compute_distance(_nodes[index_a].p, _nodes[index_b].p);
      interval_min = interval_max = dist;   

      return true;
    }

    Number_type min_dist = compute_distance(_nodes[index_a].bbox, _nodes[index_b].bbox),
                max_dist = compute_distance.maximum_distance(_nodes[index_a].bbox, _nodes[index_b].bbox);
    
    dist = (max_dist+min_dist)/Number_type(2.0);

    interval_min = min_dist;
    interval_max = max_dist;

    if( (max_dist - min_dist)/(Number_type(2.0)*min_dist) <= epsilon)
      return true;
    
    return false;

  } // bool is_well_separated_pair_epsilon(Number_type &epsilon, int index_a, int index_b)


  template <class OutputIterator> 
  void _extract_decomposition_epsilon( Number_type &epsilon, int index_a, int index_b, OutputIterator ot)
  {
     Number_type dist, interval_min, interval_max;

     if(is_well_separated_pair_epsilon(epsilon,index_a,index_b,dist, 
                                       interval_min, interval_max)==true)
     {
       Pointset_pair ppr;

       ppr.set_a.first = _nodes[index_a].p_start;
       ppr.set_a.second = _nodes[index_a].p_end;

       ppr.set_b.first = _nodes[index_b].p_start;
       ppr.set_b.second = _nodes[index_b].p_end;

       int size_a = ppr.set_a.second - ppr.set_a.first + 1,
           size_b = ppr.set_b.second - ppr.set_b.first + 1;

       _total_dist_sum += dist*Number_type(size_a)*Number_type(size_b);
       _interval_min += interval_min*Number_type(size_a)*Number_type(size_b);
       _interval_max += interval_max*Number_type(size_a)*Number_type(size_b);

       _nodes[index_a].sum_d_times_s += dist*Number_type(size_b);
       _nodes[index_a].sum_dsquare_times_s += dist*dist*Number_type(size_b);
       _nodes[index_b].sum_d_times_s += dist*Number_type(size_a);
       _nodes[index_b].sum_dsquare_times_s += dist*dist*Number_type(size_a);

       *ot++ = ppr;
     }
     else
     {
       Number_type res_a = _nodes[index_a].bbox.maximum_dimension_length(),
                   res_b = _nodes[index_b].bbox.maximum_dimension_length();

       if(res_a<= res_b)
       {
         _extract_decomposition_epsilon(epsilon, index_a, _nodes[index_b].left_child,  ot);
         _extract_decomposition_epsilon(epsilon, index_a, _nodes[index_b].right_child, ot);
       }
       else
       {
         _extract_decomposition_epsilon(epsilon, _nodes[index_a].left_child, index_b, ot);
         _extract_decomposition_epsilon(epsilon, _nodes[index_a].right_child, index_b, ot);
       }

     } // else of if(is_well_separated(epsilon,index_a,index_b)==true) 

  } // _extract_decomposition_epsilon(...)


  void _flush_sums_to_leaves( int index, Number_type sum_d, Number_type sum_dsquare)
  {
    _nodes[index].sum_d_times_s += sum_d;
    _nodes[index].sum_dsquare_times_s += sum_dsquare;

    if(_nodes[index].is_leaf() == false)
    {
      _flush_sums_to_leaves(_nodes[index].left_child, _nodes[index].sum_d_times_s, _nodes[index].sum_dsquare_times_s);
      _flush_sums_to_leaves(_nodes[index].right_child, _nodes[index].sum_d_times_s, _nodes[index].sum_dsquare_times_s);
    }
 
  } // _flush_sums_to_leaves( int index, ...)

  void flush_sums_to_leaves()
  { 
    if(_nodes[root_index()].is_leaf() == false)
    {
      _flush_sums_to_leaves(_nodes[root_index()].left_child, Number_type(0.0), Number_type(0.0));
      _flush_sums_to_leaves(_nodes[root_index()].right_child, Number_type(0.0), Number_type(0.0));
    }

    _total_sums = Number_type(0.0);

    for(int i=0; i<_nodes.size(); i++)
      if(_nodes[i].is_leaf() == true)
        _total_sums += (_nodes[i].sum_d_times_s*_nodes[i].sum_d_times_s) - _nodes[i].sum_dsquare_times_s;

  } // flush_sums_to_leaves()


  Number_type extract_stored_sums()
  { return _total_sums;}

  bool is_valid_decomposition( Number_type &sf, std::vector<Point> &points, 
                               std::vector<Pointset_pair> &pprs)
  {
    Compute_Euclidean_distance compute_distance;

    for(int i=0; i<pprs.size(); i++)
    {
      Number_type max_dist(0.0);

      for(int j=pprs[i].set_a.first; j<=pprs[i].set_a.second; j++) 
        for(int k=pprs[i].set_a.first; k<j; k++)
        {
          Number_type dist = compute_distance(points[k], points[j]);

          if(max_dist < dist)
            max_dist = dist; 
        }

      for(int j=pprs[i].set_b.first; j<=pprs[i].set_b.second; j++)
        for(int k=pprs[i].set_b.first; k<j; k++)
        {
          Number_type dist = compute_distance(points[k], points[j]);

          if(max_dist < dist)
            max_dist = dist; 
        }

      Number_type min_dist_ab = compute_distance(points[pprs[i].set_a.first], 
                                                 points[pprs[i].set_b.first]);
 
      for(int j=pprs[i].set_a.first; j<=pprs[i].set_a.second; j++) 
        for(int k=pprs[i].set_b.first; k<=pprs[i].set_b.second; k++) 
        {
          Number_type dist = compute_distance(points[k], points[j]);

          if(min_dist_ab > dist)
            min_dist_ab = dist; 
        }

      if(min_dist_ab < sf*max_dist )
        return false;

    } // for(int i=0; i<pprs.size(); i++)

    return true;   

  } // is_valid_decomposition(...)

  template< class OutputIterator>
  void get_all_points_in_subtree(int i , OutputIterator ot)
  {
    if(this->node(i).left_child == -1 && this->node(i).right_child == -1 )
      *ot++ = this->node(i).p;
    else
    {
      get_all_points_in_subtree(this->node(i).left_child , ot);
      get_all_points_in_subtree(this->node(i).right_child , ot);
    }
  }

  void print_node(int i)
  {
    std::cout << " /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/" << std::endl;
    std::cout << " Node index: " << i << std::endl;
    std::cout << " Left child: " << _nodes[i].left_child << std::endl;
    std::cout << " Right child: " << _nodes[i].right_child << std::endl;
    std::cout << " Parent: " << _nodes[i].parent << std::endl;
    std::cout << " P_start: " << _nodes[i].p_start << std::endl;
    std::cout << " P_end: " << _nodes[i].p_end << std::endl;
    std::cout << " Point: " ; 

    for(int j=0; j<_nodes[i].p.dim(); j++)
      std::cout << _nodes[i].p(j) << " ";

    std::cout << std::endl;
    std::cout << " /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/" << std::endl;
  }

  void print_tree()
  {
    std::cout << " ################# " << std::endl;

    for(int i=0; i<number_of_nodes(); i++)
      print_node(i);
  }

 private:

  std::vector<Point> _points, _sorted_points;
  int _d, _total_number_of_checks, _decomposition_size;
  std::vector<Node_type> _nodes;
  Number_type _total_sums, _total_dist_sum, _interval_min, _interval_max;

}; // class Split_tree

} // namespace FunctionalMeasures


#endif // SPLIT_TREE_H
