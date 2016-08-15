#ifndef SUCCESSOR_PRODUCT_STRUCTURE_FAST
#define SUCCESSOR_PRODUCT_STRUCTURE_FAST

#include<vector>

namespace FunctionalMeasures
{

template<class KERNEL_TYPE>
class Successor_product_structure_fast
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
    Protected_number_type iprod;
    Number_type coord;
    int lc,rc,pr;

    Container_entry()
    { lc=rc=pr=-1;}

  }; // struct Container_entry

  typedef std::vector<Container_entry>  Container_type;

  Successor_product_structure_fast(){}

  void construct_structure( std::vector<Number_type> &coords, 
                             std::vector<Number_type> &probs)
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
      en.mark = false;
   
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
        en.iprod = Protected_number_type(1.0);
        en.mark = false;
 
        _container.push_back(en);

        _container[top_level_nodes[2*i]].pr = _container.size()-1;
        _container[top_level_nodes[(2*i)+1]].pr = _container.size()-1;

        new_top_level_nodes.push_back(_container.size()-1);

      } // for(i=0; i<new_level_size(); i++)
       
      if(top_level_nodes.size()%2==1)
        new_top_level_nodes.push_back(top_level_nodes.back());

      top_level_nodes = new_top_level_nodes;

    } // while(top_level_nodes.size() > 1)

  } // construct_structure(...)

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
    _container[index].mark = true;

    if(en.pr != -1)
    {
      Container_entry fthr = _container[en.pr];

      if(fthr.lc == index)
      {
        Container_entry enR = _container[fthr.rc];
 
        if( is_leaf(fthr.rc) && enR.mark == false) // The right sibling is always a leaf
          _container[en.pr].iprod = en.iprod;
        else
          _container[en.pr].iprod = en.iprod*enR.iprod;

      } // if(fthr.lc == index)
      else
      {
        Container_entry enL = _container[fthr.lc];

        if( is_leaf(fthr.lc) && enL.mark == false)
          _container[en.pr].iprod = en.iprod;
        else
          _container[en.pr].iprod = en.iprod*enL.iprod;

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

      if( is_leaf(en.rc) && enR.mark == false)
        _container[index].iprod = enL.iprod;
      else
        _container[index].iprod = enL.iprod*enR.iprod;

    } // while(en.pr != -1)

  } // add_mark(int i)

  Protected_number_type
  product_elements_larger_than_y(Number_type &y)
  { return _recursive_product_larger_than_y(root_index(),y); }

 private:

  Protected_number_type
  _recursive_product_larger_than_y(int index, Number_type &y)
  {
     assert(index>=0 && index < _container.size());

     Container_entry en = _container[index];

     if(en.coord <= y)
       return Protected_number_type(1.0);

     if(is_leaf(index) && en.mark == false)
       return Protected_number_type(1.0);

     if(is_leaf(index))
       return en.iprod;

     Protected_number_type resL, resR;

     resL = _recursive_product_larger_than_y(en.lc,y);

     if(_container[en.lc].coord >= y && !is_leaf(en.rc))
       resR = _container[en.rc].iprod;
     else
       resR = _recursive_product_larger_than_y(en.rc,y);

     return resL*resR;

  } // _recursive_product_larger_than_y(...)

  Container_type _container;
  
}; //class Successor_product_structure_fast


} // namespace FunctionalMeasures

#endif //SUCCESSOR_PRODUCT_STRUCTURE_FAST
