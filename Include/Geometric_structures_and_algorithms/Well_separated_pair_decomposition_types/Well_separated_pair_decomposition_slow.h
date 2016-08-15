#ifndef WELL_SEPARATED_PAIR_DECOMPOSITION_SLOW_H
#define WELL_SEPARATED_PAIR_DECOMPOSITION_SLOW_H

#include<vector>
#include<set>

namespace FunctionalMeasures
{

  template<class KERNEL_TYPE>
  class Split_tree_slow
  {
   
   public:

    typedef KERNEL_TYPE                                          Kernel;
    typedef typename Kernel::Number_type                        Number_type;
    typedef typename Kernel::Point                              Point;
    typedef typename Kernel::Well_separated_pair_decomposition  Well_separated_pair_decomposition;
    typedef typename Kernel::Pointset_pair                      Pointset_pair;
    typedef typename Kernel::Numeric_traits                     Numeric_traits;
    typedef typename Numeric_traits::Square_root                Square_root;

    struct Node_type
    {
      Node_type():left_child(-1),right_child(-1){}  

      int left_child,right_child, p_start, p_end;
      Point p;
    };

   public:

    Split_tree_slow(){}

    Node_type node(int i)
    { return nodes[i];}

    Point sorted_point(int i)
    { return sorted_points[i];}

    int number_of_nodes()
    { return nodes.size();}

    int number_of_sorted_points()
    { return sorted_points.size();}

    int root_index()
    { return 0;}

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

    } // get_all_points_in_subtree(int i , OutputIterator ot)

    template<class RangeIterator>
    void load_point_set(RangeIterator r_begin, RangeIterator r_end)
    {
       points.clear();
       points.insert(points.begin(), r_begin, r_end);

       nodes.push_back(Node_type());

       if(points.size() > 0)
       {
         expand_node(points, 0);
         perform_in_order_traversal(0);
       }

    } // load_point_set(RangeIterator r_begin, RangeIterator r_end)


    std::pair<int,int> perform_in_order_traversal(int index)
    {
      if(nodes[index].left_child == -1 && nodes[index].right_child == -1)
        return std::make_pair(nodes[index].p_start, nodes[index].p_end); 

      std::pair<int,int>  ppr_l =  perform_in_order_traversal(nodes[index].left_child),
                           ppr_r =  perform_in_order_traversal(nodes[index].right_child);

      nodes[index].p_start = ppr_l.first;
      nodes[index].p_end = ppr_r.second;

      return std::make_pair(ppr_l.first, ppr_r.second);

    } // perform_in_order_traversal(int index)
   

    void expand_node(std::vector<Point> &pts, int index)    
    {
      if(pts.size() == 1)
      {
        nodes[index].p = pts[0];
        sorted_points.push_back(pts[0]);
        nodes[index].p_start = sorted_points.size()-1;
        nodes[index].p_end = sorted_points.size()-1;
        return;

      } // if(nps == 1)

      int dmax = find_maximum_dimension(pts);
      Number_type split_value = find_split_value(pts);

      std::vector<Point> pts_left, pts_right;

      int count_smaller=0,
          count_larger=0;

      // Find out whether the number of points with
      // strictly smaller coordinate than the split value
      // are less than/equal to the number of points
      // with coordinate larger than the split value.
      // This is important to make sure that this
      // test implementation produces the same results
      // with the theoretically optimal implementation. 

      for(int i=0; i<pts.size(); i++)
        if(pts[i][dmax] < split_value)
          count_smaller++;
        else if(pts[i][dmax] > split_value)
          count_larger++;

      if(count_smaller <= count_larger)
      {
        for(int i=0; i<pts.size(); i++)
          if(pts[i][dmax] < split_value)
            pts_left.push_back(pts[i]);
          else
            pts_right.push_back(pts[i]);

      } // if(count_smaller <= count_larger)
      else
      {
        for(int i=0; i<pts.size(); i++)
          if(pts[i][dmax] <= split_value) // '<=' instead of '<'
            pts_left.push_back(pts[i]);
          else
            pts_right.push_back(pts[i]);

      } // else of if(count_smaller <= count_larger)

      nodes.push_back(Node_type());  
      nodes.push_back(Node_type());

      nodes[index].left_child = nodes.size()-2;
      nodes[index].right_child = nodes.size()-1;

      expand_node(pts_left, nodes[index].left_child);
      expand_node(pts_right, nodes[index].right_child);

    } // expand_node(std::vector<Point> pts, int index)

    template<class OutputIterator_a, class OutputIterator_b>
    void extract_decomposition(Number_type &sf, OutputIterator_a ot_a, OutputIterator_b ot_b)
    {
      if(points.size()>0)
      {
        extract_decomposition(sf, 0, ot_b);

        for(int i=0; i<sorted_points.size(); i++)
          *ot_a++ = sorted_points[i];
      }

    } // extract_decomposition(Number_type &sf, OutputIterator_a, OutputIterator_b)


    template<class OutputIterator>
    void extract_decomposition(Number_type &sf, int index, OutputIterator ot)
    {
      if(nodes[index].left_child!= -1 && nodes[index].right_child != -1 )
      {
        extract_decomposition(sf, nodes[index].left_child, nodes[index].right_child, ot);
        extract_decomposition(sf, nodes[index].left_child, ot);
        extract_decomposition(sf, nodes[index].right_child, ot);
      }  
    }

    template<class OutputIterator>
    void extract_decomposition(Number_type &sf, int n_a, int n_b, OutputIterator ot)
    {
      std::vector<Point> vec_a, vec_b;

      for(int i=nodes[n_a].p_start; i<=nodes[n_a].p_end; i++)
        vec_a.push_back(sorted_points[i]);

      for(int i=nodes[n_b].p_start; i<=nodes[n_b].p_end; i++)
        vec_b.push_back(sorted_points[i]);

      if(is_well_separated_pair(sf, vec_a, vec_b))
      {
        Pointset_pair pts_pr;

        pts_pr.set_a.first = nodes[n_a].p_start;
        pts_pr.set_a.second = nodes[n_a].p_end;
        pts_pr.set_b.first = nodes[n_b].p_start;
        pts_pr.set_b.second = nodes[n_b].p_end;

        *ot++ = pts_pr;
        
        return;

      } // if(is_well_separated_pair(sf, n_a,n_b))

      Number_type len_a = find_length_of_maximum_dimension(vec_a),
                  len_b = find_length_of_maximum_dimension(vec_b);

      if(len_b >= len_a)
      {
        extract_decomposition(sf, n_a, nodes[n_b].left_child, ot);
        extract_decomposition(sf, n_a, nodes[n_b].right_child, ot);
      }
      else
      {
        extract_decomposition(sf, nodes[n_a].left_child, n_b, ot);
        extract_decomposition(sf, nodes[n_a].right_child, n_b, ot);
      }

    } // extract_decomposition(Number_type &sf, int n_a, int n_b, OutputIterator ot)

    bool is_well_separated_pair( Number_type &sf, 
                                 std::vector<Point> &pts_a, 
                                 std::vector<Point> &pts_b )
    {
      if(pts_a.size() == 1 && pts_b.size() == 1 )
        return true;

      std::vector<std::pair<Number_type, Number_type> > bb_a, bb_b;
      
      int dim = pts_a[0].size();   

      bb_a.assign(dim, std::make_pair(Number_type(0.0),Number_type(0.0)) );
      bb_b.assign(dim, std::make_pair(Number_type(0.0),Number_type(0.0)) );

      for(int i=0; i<dim; i++)
      {
        bb_a[i].first = pts_a[0][i];
        bb_a[i].second = pts_a[0][i];

        for(int j=1; j<pts_a.size(); j++)
        {
          if(pts_a[j][i] < bb_a[i].first) 
            bb_a[i].first = pts_a[j][i];      

          if(pts_a[j][i] > bb_a[i].second) 
            bb_a[i].second = pts_a[j][i];  
        }

      } // for(int i=0; i<dim; i++)  


      for(int i=0; i<dim; i++)
      {
        bb_b[i].first = pts_b[0][i];
        bb_b[i].second = pts_b[0][i];

        for(int j=1; j<pts_b.size(); j++)
        {
          if(pts_b[j][i] < bb_b[i].first) 
            bb_b[i].first = pts_b[j][i];      

          if(pts_b[j][i] > bb_b[i].second) 
            bb_b[i].second = pts_b[j][i];  
        }

      } // for(int i=0; i<dim; i++)   

      Number_type diameter_a(0.0), 
                  diameter_b(0.0), 
                  diameter; 

      for(int i=0; i<dim; i++)
      {
        diameter_a += (bb_a[i].second - bb_a[i].first)*(bb_a[i].second - bb_a[i].first);      
        diameter_b += (bb_b[i].second - bb_b[i].first)*(bb_b[i].second - bb_b[i].first);     
 
      } // for(int i=0; i<dim; i++)  
     
      if(diameter_a > diameter_b)
        diameter = diameter_a;
      else
        diameter = diameter_b;
    
      Number_type min_dist(0.0);

      for(int i=0; i<dim; i++)
      {
        if(bb_a[i].second < bb_b[i].first)     
          min_dist += (bb_b[i].first - bb_a[i].second)*(bb_b[i].first - bb_a[i].second);

        if(bb_b[i].second < bb_a[i].first)     
          min_dist += (bb_a[i].first - bb_b[i].second)*(bb_a[i].first - bb_b[i].second);     
 
      } // for(int i=0; i<dim; i++)

      min_dist = Square_root()(min_dist);
      diameter = Square_root()(diameter);

      if(min_dist >= sf*diameter)
        return true;

      return false;

    } // is_well_separated_pair( Number_type &sf, ... )

    int find_maximum_dimension(std::vector<Point> &pts)
    {
      int dim = pts[0].dim();

      Number_type max_length(0.0);
      int max_dim = 0;

      for(int i=0; i<dim; i++)
      {
        Number_type max_v = pts[0][i], 
                    min_v = pts[0][i];

        for(int j=0; j<pts.size(); j++)
        {
          if(pts[j][i] > max_v)
            max_v = pts[j][i]; 

          if(pts[j][i] < min_v)
            min_v = pts[j][i]; 
        } // for(int j=0; j<nps.size(); j++)

        if(max_v-min_v > max_length)
        {
          max_length = max_v - min_v;
          max_dim = i;
        }

      } // for(int i=0; i<dim; i++)

      return max_dim;

    } // find_maximum_dimension(std::vector<Point> &pts)

    Number_type find_length_of_maximum_dimension(std::vector<Point> &pts)
    {
      int dim = pts[0].dim();

      Number_type max_length(0.0);

      for(int i=0; i<dim; i++)
      {
        Number_type max_v = pts[0][i], 
                    min_v = pts[0][i];

        for(int j=0; j<pts.size(); j++)
        {
          if(pts[j][i] > max_v)
            max_v = pts[j][i]; 

          if(pts[j][i] < min_v)
            min_v = pts[j][i]; 
        } // for(int j=0; j<nps.size(); j++)

        if(max_v-min_v > max_length)
          max_length = max_v - min_v;

      } // for(int i=0; i<dim; i++)

      return max_length;

    } // find_maximum_dimension(std::vector<Point> &pts)


    Number_type find_split_value(std::vector<Point> &pts)
    {
      int dim = pts[0].dim();

      Number_type max_length(0.0), global_min_v, global_max_v;
      int max_dim = 0;

      for(int i=0; i<dim; i++)
      {
        Number_type max_v = pts[0][i], 
                    min_v = pts[0][i];

        for(int j=0; j<pts.size(); j++)
        {
          if(pts[j][i] > max_v)
            max_v = pts[j][i]; 

          if(pts[j][i] < min_v)
            min_v = pts[j][i]; 
        } // for(int j=0; j<nps.size(); j++)

        if(max_v-min_v > max_length)
        {
          max_length = max_v - min_v;
          max_dim = i;
          global_max_v = max_v;
          global_min_v = min_v;
        }

      } // for(int i=0; i<dim; i++)

      return Number_type((global_min_v+global_max_v)/Number_type(2.0));

    } // find_split_value(std::vector<Point> &pts)


    void print_node(int i)
    {
      std::cout << " /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/" << std::endl;
      std::cout << " Node index: " << i << std::endl;
      std::cout << " Left child: " << nodes[i].left_child << std::endl;
      std::cout << " Right child: " << nodes[i].right_child << std::endl;
      std::cout << " Parent: " << nodes[i].parent << std::endl;
      std::cout << " P_start: " << nodes[i].p_start << std::endl;
      std::cout << " P_end: " << nodes[i].p_end << std::endl;
      std::cout << " Point: " ; 

      for(int j=0; j<nodes[i].p.dim(); j++)
        std::cout << nodes[i].p(j) << " ";

      std::cout << std::endl;
      std::cout << " /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/" << std::endl;

    } // void print_node(int i)


    void print_tree()
    {
      std::cout << " ################# " << std::endl;

      for(int i=0; i<number_of_nodes(); i++)
        print_node(i);
    }

    std::vector<Point> points, sorted_points;
    std::vector<Node_type>  nodes; 
    
  }; // class Split_tree_slow



  template<class KERNEL_TYPE>
  class Well_separated_pair_decomposition_slow
  {
   
   public:

    typedef KERNEL_TYPE                      Kernel;
    typedef typename Kernel::Number_type    Number_type;
    typedef typename Kernel::Point          Point;
    typedef typename Kernel::Split_tree     Split_tree_fast;
    typedef typename Kernel::Pointset_pair  Pointset_pair;
    typedef Split_tree_slow<Kernel>          Split_tree_sl;

   public:

    Well_separated_pair_decomposition_slow(){}

    template<class OutputIterator_a, class OutputIterator_b>
    void extract_decomposition(Number_type &sf, OutputIterator_a ot_a, OutputIterator_b ot_b)
    { tree.extract_decomposition(sf,ot_a,ot_b);}

    template<class RangeIterator>
    void load_point_set(RangeIterator r_begin, RangeIterator r_end)
    { tree.load_point_set(r_begin, r_end);}


    bool is_identical_decomposition_tree( Split_tree_fast &tf)
    { return is_identical_decomposition_tree(tf, tf.root_index(), tree.root_index() );}

    bool is_identical_decomposition_tree( Split_tree_fast &tf, int ind_f, int ind_t )
    {  
      std::set<Point> pf, pt;

      for(int tt=tree.node(ind_t).p_start; tt<=tree.node(ind_t).p_end; tt++)
        pt.insert(tree.sorted_point(tt));
         
      for(int tt=tf.node(ind_f).p_start; tt<=tf.node(ind_f).p_end; tt++)
        pf.insert(tf.sorted_point(tt));

      typename std::set<Point>::iterator it_tr=pt.begin(),
                                          it_tf=pf.begin();

      for( ; it_tr !=pt.end(); it_tr++, it_tf++)
        if( *it_tr != *it_tf)
        {
          std::cout<< " There was a discrepancy on the point set appearing in a node's subtree." << std::endl;
          std::cout<< " The node index of the optimal tree where the difference was observed is: " << ind_f << std::endl;
          std::cout<< " The node index of the suboptimal tree where the difference was observed is: "<< ind_t << std::endl;
          std::cout<< std::endl << " The point set in the subtree of the optimal tree is: " << std::endl;

          for( it_tf=pf.begin(); it_tf !=pf.end(); it_tf++)
            std::cout << (*it_tf) << std::endl;

          std::cout<< std::endl << " The point set in the subtree of the suboptimal tree is: " << std::endl;

          for( it_tr=pt.begin(); it_tr !=pt.end(); it_tr++)
            std::cout << (*it_tr) << std::endl;

          return false;

        } // if( *it_tr != *it_tf)


      int ct=0, cf=0;

      if(tree.node(ind_t).left_child == -1)
        ct++;

      if(tree.node(ind_t).right_child == -1)
        ct++;

      if(tf.node(ind_f).left_child == -1)
        cf++;

      if(tf.node(ind_f).right_child == -1)
        cf++;

       if( ct != cf)
       {
         std::cout<< " There was a discrepancy on a node's children indices." << std::endl;
         std::cout<< " The node index of the optimal tree where the difference was observed is: " << ind_f << std::endl;
         std::cout<< " The node index of the suboptimal tree where the difference was observed is: "<< ind_t << std::endl;
         std::cout<< " The number of children for the optimal tree is: " << cf << std::endl;
         std::cout<< " The number of children the suboptimal tree is: " << ct << std::endl;

         return false;
       }

       bool res_l,res_r; 

       if(tf.node(ind_f).left_child != -1)
         res_l = is_identical_decomposition_tree(tf, tf.node(ind_f).left_child, tree.node(ind_t).left_child );
             
       if(tf.node(ind_f).right_child != -1)
         res_r = is_identical_decomposition_tree(tf, tf.node(ind_f).right_child, tree.node(ind_t).right_child );

       if(tf.node(ind_f).left_child == -1 && tf.node(ind_f).right_child == -1)
         return true;
       else
         return (res_l && res_r);

    } // is_identical_decomposition_tree(...)

    bool are_identical_decompositions( std::vector<Point> &pts_a, std::vector<Pointset_pair> &prs_a, 
                                        std::vector<Point> &pts_b, std::vector<Pointset_pair> &prs_b)
    {
      //std::cout << " A1" << std::endl;

      if(pts_a.size() != pts_b.size())
        return false;

      //std::cout << " A2" << std::endl;

      for(int i=0; i<pts_a.size(); i++)
        if(pts_a[i].dim() != pts_b[i].dim())
          return false;

      //std::cout << " A3" << std::endl;

      std::set<Point> sp_a, sp_b;
          
      for(int i=0; i<pts_a.size(); i++)
        sp_a.insert(pts_a[i]);

      for(int i=0; i<pts_b.size(); i++)
        sp_b.insert(pts_b[i]);

      if(sp_a.size() != sp_b.size())
        return false;
      
      //std::cout << " A4" << std::endl;

      int dim=pts_a[0].dim();

      typename std::set<Point>::iterator it_a = sp_a.begin(), it_b;

      for(it_b = sp_b.begin(); it_b != sp_b.end(); it_a++, it_b++ )
        for(int j=0; j<dim; j++)
          if(*it_b != *it_a)
            return false;

      // std::cout << " A5" << std::endl;

      /*
      std::cout << " Pairs of decomposition A: " << std::endl;   

      for(int i=0; i<prs_a.size(); i++)
        std::cout << prs_a[i].set_a.first << " - " << prs_a[i].set_a.second << "||" 
                  << prs_a[i].set_b.first << " - " << prs_a[i].set_b.second << std::endl;

      std::cout << " Pairs of decomposition B: " << std::endl;   

      for(int i=0; i<prs_b.size(); i++)
        std::cout << prs_b[i].set_a.first << " - " << prs_b[i].set_a.second << "||" 
                  << prs_b[i].set_b.first << " - " << prs_b[i].set_b.second << std::endl;
      */

      if(prs_a.size() != prs_b.size())
        return false;

      //std::cout << " A6" << std::endl;

      for(int i=0; i<prs_a.size(); i++)
        if( prs_a[i].set_a.second - prs_a[i].set_a.first != prs_b[i].set_a.second - prs_b[i].set_a.first )
            return false; 

      //std::cout << " A7" << std::endl;

      for(int i=0; i<prs_a.size(); i++)
        if( prs_a[i].set_b.second - prs_a[i].set_b.first != prs_b[i].set_b.second - prs_b[i].set_b.first )
            return false; 

      //std::cout << " A8" << std::endl;

      for(int i=0; i<prs_a.size(); i++)
      {
        int j, k=prs_b[i].set_a.first;

        for(int j=prs_a[i].set_a.first; j<=prs_a[i].set_a.second; j++, k++)
          if( pts_a[j] != pts_b[k] )
            return false; 
      }
 
      //std::cout << " A9" << std::endl;

      for(int i=0; i<prs_a.size(); i++)
      {
        int j, k=prs_b[i].set_b.first;

        for(int j=prs_a[i].set_b.first; j<=prs_a[i].set_b.second; j++, k++)
          if( pts_a[j] != pts_b[k] )
            return false; 
      }

      //std::cout << " A10" << std::endl;

      return true;

    } // are_identical_decompositions(...)

   public:
    
    Split_tree_sl tree;    
 
  }; // class Well_separated_decomposition_slow

} // namespace FunctionalMeasures

#endif // WELL_SEPARATED_PAIR_DECOMPOSITION_SLOW_H
