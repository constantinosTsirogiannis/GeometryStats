#ifndef COMPUTE_EUCLIDEAN_DISTANCE_H
#define COMPUTE_EUCLIDEAN_DISTANCE_H

#include<vector>
#include<cassert>

namespace FunctionalMeasures
{

template <class KERNEL_TYPE>
class Compute_Euclidean_distance
{
 public:

  typedef KERNEL_TYPE                            Kernel;
  typedef typename Kernel::Number_type          Number_type;
  typedef typename Kernel::Point                Point;
  typedef typename Kernel::Bbox_d               Bbox_d;
  typedef typename Kernel::Numeric_traits       Numeric_traits;
  typedef typename Numeric_traits::Square_root  Square_root;
 
  Number_type operator()(const Point &p1, const Point &p2)
  {
    assert(p1.dim() == p2.dim());

    Number_type val(0.0);
 
    for(int i=0; i<p1.dim(); i++)
    {
      Number_type diff = p1[i] - p2[i];

      val = val + (diff*diff);
    }

    return Square_root()(val);  

  } // operator()(const Point &p1, const Point &p2)

  Number_type operator()(const Bbox_d &b1, const Bbox_d &b2)
  {
    assert(b1.dim() == b2.dim());

    Number_type val(0.0);
 
    for(int i=0; i<b1.dim(); i++)
    {
      if(b1.min(i) > b2.max(i))
      {
        Number_type diff = b1.min(i) - b2.max(i);
        val = val + (diff*diff);
      }
      else if(b1.max(i) < b2.min(i))
      {
        Number_type diff = b2.min(i) - b1.max(i);
        val = val + (diff*diff);
      }

    } // for(int i=0; i<b1.dim(); i++)

    return Square_root()(val);  

  } // operator()(const Bbox_d &b1, const Bbox_d &b2)

  Number_type maximum_distance(const Bbox_d &b1, const Bbox_d &b2)
  {
    assert(b1.dim() == b2.dim());

    Number_type val(0.0);
 
    for(int i=0; i<b1.dim(); i++)
    {
      if(b1.min(i) >= b2.max(i))
      {
        Number_type diff = b1.max(i) - b2.min(i);
        val = val + (diff*diff);
      }
      else if(b1.max(i) <= b2.min(i))
      {
        Number_type diff = b2.max(i) - b1.min(i);
        val = val + (diff*diff);
      }
      else 
      {
        Number_type diff_1 = b1.max(i) - b2.min(i),
                    diff_2 = b2.max(i) - b1.min(i);

        diff_1 = diff_1*diff_1;
        diff_2 = diff_2*diff_2;

        if(diff_1> diff_2)
          val = val + diff_1;
        else
          val = val + diff_2;
      }

    } // for(int i=0; i<b1.dim(); i++)

    return Square_root()(val);  

  } // maximum_distance(const Bbox_d &b1, const Bbox_d &b2)


  Number_type squared_distance(const Point &p1, const Point &p2)
  {
    assert(p1.dim() == p2.dim());

    Number_type val(0.0);
 
    for(int i=0; i<p1.dim(); i++)
    {
      Number_type diff = p1[i] - p2[i];

      val += diff*diff;
    }

    return val;  

  } // operator()(const Point &p1, const Point &p2)

}; // class Compute_Euclidean_distance

} // namespace FunctionalMeasures

#endif // COMPUTE_EUCLIDEAN_DISTANCE_H
