#ifndef BBOX_D_H
#define BBOX_D_H

namespace FunctionalMeasures
{

template<class KERNEL_TYPE>
class Bbox_d
{
 public:

  typedef KERNEL_TYPE                            Kernel;
  typedef typename Kernel::Number_type          Number_type;
  typedef typename Kernel::Numeric_traits       Numeric_traits;
  typedef typename Numeric_traits::Square_root  Square_root;

 public: 
 
  Bbox_d(){}

  Number_type min( int i) const
  { return extrema[i].first; }

  Number_type max( int i) const
  { return extrema[i].second; }

  Number_type& min( int i)
  { return extrema[i].first; }

  Number_type& max( int i)
  { return extrema[i].second; }

  int dim() const
  { return extrema.size();}

  Number_type diameter()
  {
    Number_type val(0.0);

    for(int i=0; i<dim(); i++)
    {
      Number_type diff =  extrema[i].second - extrema[i].first;
      val = val + (diff*diff); 
    }

    return Square_root()(val);  
  }

  void set_dimension(int i)
  {  
    extrema.clear();
    extrema.assign(i,std::make_pair<Number_type,Number_type>(0.0,0.0));
  }

  std::pair<int, Number_type> maximum_dimension()
  {
    if(this->dim()<1)
      return std::make_pair(0,Number_type(0.0));

    Number_type cur_max(-1.0);
    int ind_max=0;

    for(int i=0; i<this->dim(); i++)
    {
      Number_type length = extrema[i].second - extrema[i].first;

      if(length > cur_max)
      {
        cur_max = length;
        ind_max = i;
      }

    } // for(int i=0; i<this->dim(); i++)

    Number_type split_value = extrema[ind_max].first + (cur_max/Number_type(2.0));

    return std::make_pair(ind_max,split_value);

  } // int maximum_dimension()

  Number_type maximum_dimension_length()
  {
    if(this->dim()<1)
      return Number_type(0.0);

    Number_type cur_max(-1.0);

    for(int i=0; i<this->dim(); i++)
    {
      Number_type length = extrema[i].second - extrema[i].first;

      if(length > cur_max)
        cur_max = length;

    } // for(int i=0; i<this->dim(); i++)

    return cur_max;

  } // int maximum_dimension_length()

 private:

  std::vector< std::pair<Number_type, Number_type> > extrema; 

}; // class Bbox_d

} // namespace FunctionalMeasures

#endif // BBOX_D_H
