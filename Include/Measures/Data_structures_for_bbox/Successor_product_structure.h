#ifndef SUCCESSOR_PRODUCT_STRUCTURE
#define SUCCESSOR_PRODUCT_STRUCTURE

#include<vector>

namespace FunctionalMeasures
{

template<class KERNEL_TYPE>
class Successor_product_structure
{
  
 public:

  typedef KERNEL_TYPE                                   Kernel;
  typedef typename Kernel::Number_type                 Number_type; 
  typedef typename Kernel::Protected_number_type       Protected_number_type; 
  typedef typename Kernel::Exception_functor           Exception_functor;
  typedef typename Kernel::Exception_type              Exception_type;

  struct Container_entry
  {
    Number_type coord;
    Protected_number_type reciprocal, prod;

    Container_entry(){}

    Container_entry( Number_type &n, Protected_number_type &r, Protected_number_type &p ):
    coord(n),reciprocal(r), prod(p){};

  }; // struct Container_entry

  typedef std::vector< Container_entry >      Container_type;
  typedef typename Container_type::iterator  Container_iterator;

  Successor_product_structure(){}

  void insert(Number_type &n, Number_type &p)
  {
    Protected_number_type reciprocal(Number_type(1.0)/Number_type(1.0-p));

    if(_container.size() == 0)
    {
      Protected_number_type pt(1.0);
      _container.push_back(Container_entry(n,reciprocal,pt));
      return;
    }

    Container_iterator it;

    for(it=_container.begin(); it!=_container.end(); it++)
      if(n<it->coord)  
        break;

    if(it==_container.end())
    {
      Protected_number_type pt(1.0);
      Container_entry new_en(n,reciprocal, pt);
      it = _container.insert(it, new_en);
    }
    else
    {
      Protected_number_type pt = it->prod*it->reciprocal;
      Container_entry new_en(n,reciprocal, pt);
      it = _container.insert(it, new_en);
    }
      
    if(it!=_container.begin())
      do
      {
        it--;
        it->prod = it->prod*reciprocal;
        
      }while(it !=_container.begin());  

  } // void insert(Number_type &n, Number_type &p)


  std::pair<Protected_number_type, Container_iterator> 
  return_successor_product(Number_type &n)
  {
    Container_iterator it = _container.begin();

    for(; it != _container.end(); it++)
      if(it->coord > n)
        break;

    if(it == _container.end())
      return std::make_pair(Protected_number_type(Number_type(1.0)),it);
 
    //std::cout << " Coontainer status: " << std::endl;

    //for(Container_iterator bit = _container.begin(); bit != _container.end() ; bit++)
    //{
    //std::cout << " bit->coord: " << bit->coord << std::endl; 
    //std::cout << " bit->prod: " << bit->prod.n() << " , " << bit->prod.exp() << std::endl; 
    //std::cout << " bit->reciprocal: " << bit->reciprocal.n() << " , " << bit->reciprocal.exp() << std::endl; 
    //}
 
    //std::cout << " Container size: " << _container.size() << std::endl; 
    //std::cout << " Minimum largest element rank: " << (it-_container.begin()) << std::endl; 
    //std::cout << " it->prod: " << it->prod.n() << " , " << it->prod.exp() << std::endl; 
    //std::cout << " it->reciprocal: " << it->reciprocal.n() << " , " << it->reciprocal.exp() << std::endl; 

    return std::make_pair(it->prod*it->reciprocal,it);

  } // find_successor(Number_type &n)

 private:

  Container_type _container;
  
  // In case you want to start from the last searched item:

  // int _last_search_rank;
  // Number_type _last_search_coord;

}; //class Successor_product_structure


} // namespace FunctionalMeasures

#endif //SUCCESSOR_PRODUCT_STRUCTURE
