#ifndef PRODUCT_TREE_H
#define PRODUCT_TREE_H

namespace FunctionalMeasures
{

template<class KERNEL_TYPE>
class Product_tree
{
  
 public:

  typedef KERNEL_TYPE                                   Kernel;
  typedef typename Kernel::Number_type                 Number_type; 
  typedef typename Kernel::Protected_number_type       Protected_number_type; 
  typedef typename Kernel::Exception_functor           Exception_functor;
  typedef typename Kernel::Exception_type              Exception_type;

  struct Container_entry
  {
    bool mark;
    Protected_number_type sprod, iprod;
    Number_type coord;
    int lc,rc,pr;

    Container_entry()
    { lc=rc=pr=-1;}

  }; // struct Container_entry

  typedef std::vector<Container_entry>  Container_type;
  
  Product_tree(){}

  int number_of_leaves()
  { return _container.size();}

  void construct_tree( std::vector<Number_type> &coords, 
                       std::vector<Number_type> &probs, 
                       std::vector<Protected_number_type> &cqs)
  {
    _container.clear();
    
    std::vector<int> top_level_nodes;

    ///////////////////////////
    // Construct first level //
    ///////////////////////////
 
    for(int i=0; i<coords.size(); i++)
    {
      Container_entry en;
   
      en.coord = coords[i];
      en.iprod = Protected_number_type(Number_type(1.0)/Number_type(1.0-probs[i]));
      en.sprod = cqs[i];
      en.mark = true;
   
      _container.push_back(en);
      top_level_nodes.push_back(_container.size()-1);

    } // for(int i=0; i<coords.size(); i++)

    ///////////////////////////////////////
    // Construct rest of the tree levels //
    ///////////////////////////////////////

    while(top_level_nodes.size() > 1)
    {
      int new_level_size = top_level_nodes.size()/2;      
      int i;

      std::vector<int> new_top_level_nodes;

      for(i=0; i<new_level_size; i++)
      {
        Container_entry en;

        en.lc = top_level_nodes[2*i];
        en.rc = top_level_nodes[(2*i)+1];
        Container_entry enL = _container[en.lc];
        Container_entry enR = _container[en.rc];   
        en.coord = enR.coord;
        en.iprod = enL.iprod*enR.iprod;
        en.sprod = Protected_number_type(0.0);
        en.mark = true;
 
        _container.push_back(en);

        _container[top_level_nodes[2*i]].pr = _container.size()-1;
        _container[top_level_nodes[(2*i)+1]].pr = _container.size()-1;

        new_top_level_nodes.push_back(_container.size()-1);

      } // for(i=0; i<new_level_size(); i++)
       
      if(top_level_nodes.size()%2==1)
        new_top_level_nodes.push_back(top_level_nodes.back());

      top_level_nodes = new_top_level_nodes;

    } // while(top_level_nodes.size() > 1)

  } // construct_tree(...)

  Container_entry node(int i)
  { return _container[i];}

  bool is_leaf(int i)
  { return _container[i].lc==-1;}

  bool is_root(int i)
  { return _container[i].pr==-1;}

  int root_index()
  { return _container.size()-1;}

  void add_mark(int i)
  {
    assert(is_leaf(i));

    int index = i;
    Container_entry en = _container[index];
    _container[index].mark = false;

    //std::cout << " IPROD of root before: " << _container.back().iprod.n() 
    //                              << " , " << _container.back().iprod.exp() << std::endl;

    if(en.pr != -1)
    {
      Container_entry fthr = _container[en.pr];

      if(fthr.lc == index)
      {
        Container_entry enR = _container[fthr.rc];
 
        if( is_leaf(fthr.rc) && enR.mark == true) // The right sibling is always a leaf
        {
          _container[en.pr].iprod = enR.iprod;
          _container[en.pr].sprod = en.sprod*enR.iprod; 
        }
        else
        {
          _container[en.pr].iprod = Protected_number_type(1.0);
          _container[en.pr].sprod = en.sprod+enR.sprod;
        } 

      } // if(fthr.lc == index)
      else
      {
        Container_entry enL = _container[fthr.lc];

        if( is_leaf(fthr.lc) && enL.mark == true)
        {
          _container[en.pr].iprod = enL.iprod;
          _container[en.pr].sprod = en.sprod;
        }
        else if(is_leaf(fthr.lc))
        {
          _container[en.pr].iprod = Protected_number_type(1.0);
          _container[en.pr].sprod = enL.sprod+en.sprod;
        }       
        else
        {
          _container[en.pr].iprod = enL.iprod;
          _container[en.pr].sprod = enL.sprod+en.sprod;
        } 

      } // else of if(fthr.lc == index)

      en = _container[en.pr];

    } // if(en.pr != -1)

    // From nowon, en is definitely not a leaf node, 
    // but it's sibling could be, and we have to 
    // account for the case that it is a marked sibling.
    // Given the way that we constructed the tree, 
    // the only case that this can happen is if en is 
    // the left of the two siblings.

    while(en.pr != -1)
    {
      index = en.pr;
      en = _container[index];
      Container_entry enL = _container[en.lc];   
      Container_entry enR = _container[en.rc]; 

      if( is_leaf(en.rc) && enR.mark == true)
      {
        _container[index].iprod = enL.iprod*enR.iprod;
        _container[index].sprod = enL.sprod*enR.iprod;
      }
      else if( is_leaf(en.rc) && enR.mark == false)
      {
        _container[index].iprod = enL.iprod;
        _container[index].sprod = enL.sprod+enR.sprod;
      }  
      else
      {
        _container[index].iprod = enL.iprod*enR.iprod;
        _container[index].sprod = (enL.sprod*enR.iprod)+enR.sprod;      
      }

    } // while(en.pr != -1)

    //std::cout << " IPROD of root after: " << _container.back().iprod.n() 
    //                              << " , " << _container.back().iprod.exp() << std::endl;

  } // add_mark(int i)

  Protected_number_type
  product_elements_larger_than_y(Number_type &y)
  { return _recursive_product_larger_than_y(root_index(),y).second; }

 private:

  std::pair<Protected_number_type, Protected_number_type>
  _recursive_product_larger_than_y(int index, Number_type &y)
  {
     assert(index>=0 && index < _container.size());

     //std::cout << "&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
     //std::cout << " Node index: " << index << std::endl;

     Container_entry en = _container[index];

     //std::cout << " Node coord: " << en.coord << std::endl;
     //std::cout << " Node mark: " << en.mark << std::endl;
     //std::cout << " Node iprod: " << en.iprod.n() << " , " << en.iprod.exp() << std::endl;
     //std::cout << " Node sprod: " << en.sprod.n() << " , " << en.sprod.exp() << std::endl;
     //std::cout << " Query: " << y << std::endl;

     if(en.coord <= y)
     {
       //std::cout << "  case en.coord <= y " << std::endl;

       return std::make_pair(Protected_number_type(1.0),Protected_number_type(0.0));
     }

     if(is_leaf(index) && en.mark == true)
     {
       //std::cout << "  case is_leaf && mark == true " << std::endl;

       return std::make_pair(en.iprod,Protected_number_type(0.0));
     }


     if(is_leaf(index))
     {
       //std::cout << "  case is leaf " << std::endl;

       return std::make_pair(Protected_number_type(1.0),en.sprod);
     }

     std::pair<Protected_number_type, Protected_number_type> resL, resR;

     //std::cout << " Recurse to left child " << std::endl;

     resL = _recursive_product_larger_than_y(en.lc,y);

     if(_container[en.lc].coord >= y && !is_leaf(en.rc))
     {
       //std::cout << " Do not recurse to right child " << std::endl;
       resR = std::make_pair(_container[en.rc].iprod, _container[en.rc].sprod);
     }
     else
     {
       //std::cout << " Recurse to right child " << std::endl;
       resR = _recursive_product_larger_than_y(en.rc,y);
     }

     return std::make_pair(resL.first*resR.first, (resL.second*resR.first)+resR.second);

  } // _recursive_product_larger_than_y(...)

 private:

  Container_type  _container;

}; // Product_tree

} //namespace FunctionalMeasures

#endif //PRODUCT_TREE_H
