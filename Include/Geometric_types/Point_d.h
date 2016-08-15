#ifndef POINT_D_H
#define POINT_D_H

#include<vector>
#include<cassert>
#include<string>
#include<iostream>

namespace FunctionalMeasures
{

template<class KERNEL_TYPE>
class Point_d: public std::vector<typename KERNEL_TYPE::Number_type> 
{
 public:

  typedef KERNEL_TYPE                         Kernel;
  typedef typename KERNEL_TYPE::Number_type  Number_type;
  typedef std::vector<Number_type>            Base;
  typedef Point_d<Kernel>                     Self;

 public:
  
  std::string taxon;

 public:

  Point_d() : _index(0) {} // _is_reflected(false)

  Point_d(int d) : _index(0)  // ,_is_reflected(false) 
  { this->assign(d,Number_type(0.0)); }

  template<class RangeIterator>
  Point_d(RangeIterator r_begin, RangeIterator r_end) : _index(0) //,_is_reflected(false) 
  { 
    this->clear();
    this->insert(this->begin(), r_begin, r_end);
  }

  template<class RangeIterator>
  void set_coordinates(RangeIterator r_begin, RangeIterator r_end)
  { 
    this->clear();
    this->insert(this->begin(), r_begin, r_end);

  } // set_coordinates(RangeIterator r_begin, RangeIterator r_end)

  Number_type& operator()(int i)
  { return (*this)[i]; }

  Number_type operator()(int i) const
  { return (*this)[i]; }

  int dim(void) const
  { return this->size();}

  //methods for volume algo
  
  size_t index()const
  { return _index; }

  size_t set_index(size_t newi)
  {
    size_t oldi=_index;
    _index=newi;
    return oldi;
  }

 // bool get_is_reflected()const{
 //   return _is_reflected;
 // }
  
 // bool set_is_reflected(bool new_is){
 //   bool old_is = _is_reflected;
 //   _is_reflected = new_is;
 //   return old_is;
 // }
  
  void operator-=(const Point_d &q) 
  {
    for(int i=0; i<q.dim(); i++)
    {
      (*this)[i] = (*this)[i] - q(i);
    }
  }
  
  void operator+=(const Point_d &q) 
  {
    for(int i=0; i<q.dim(); i++)
    {
      (*this)[i] = (*this)[i] + q(i);
    }
  }
  
  void operator*=(const int k) 
  {
    for(int i=0; i<this->dim(); i++)
    {
      (*this)[i] = (*this)[i] * Number_type(k);
    }
  }
  
  //inner product
  Number_type operator*(const Self &p) 
  {
    Number_type result=0;

    for(int i=0; i<this->dim(); i++)
    {
      result += p[i] * (*this)[i];
    }
    return result;
  }
  
  Self operator/(const int k) 
  {
    Self result(*this);

    for(int i=0; i<this->dim(); i++)
    {
      result[i] = (*this)[i] / Number_type(k);
    }
    return result;
  }
  
  //comparisons
  
  bool operator<(const Point_d &p) const
  {
    assert(this->dim() == p.dim());

    for(int i=0; i<p.dim(); i++)
    {
      if((*this)[i] < p(i))
        return true;

      if((*this)[i] > p(i))
        return false;
    }

    return false;

  } // operator<(const Point_d<KERNEL_TYPE> &p) const


  bool operator!=(const Point_d &p) const
  {
    assert(this->dim() == p.dim());

    for(int i=0; i<p.dim(); i++)
    {
      if((*this)[i] != p(i))
        return true;
    }

    return false;

  } // operator<(const Point_d<KERNEL_TYPE> &p) const

private:

  size_t _index;
  //bool _is_reflected;
  
}; // class Point_d

  template<class KERNEL_TYPE>
  std::ostream& operator<<(std::ostream &os, const Point_d<KERNEL_TYPE> &p)
  {
     os << p.taxon << " ";

    for(int i=0; i<p.dim(); i++)
      os << p[i] << " ";

    return os;
  }

} // namespace FunctionalMeasures

#endif //POINT_D_H
