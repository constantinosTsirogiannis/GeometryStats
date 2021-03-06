#ifndef PROTECTED_NUMBER_TYPE_H
#define PROTECTED_NUMBER_TYPE_H

namespace FunctionalMeasures {

template <class NTS>
class Protected_number_type
{
 public:

  typedef NTS                                    Numeric_traits;
  typedef Protected_number_type<Numeric_traits>  Self;          
  typedef typename Numeric_traits::Number_type  Number_type;
  typedef typename Numeric_traits::Is_exact     Is_exact;
  typedef typename Numeric_traits::To_double    To_double;
  typedef typename Numeric_traits::Power        Power;

 private:

  void _set_vol()
  { 
    _vol = Number_type(1.0);
    _gap = Number_type(10.0); 
    _verbose = false;
  }

 
 public:

  Protected_number_type()
  {
    _exp=0;
    _set_vol(); 
  }

  Protected_number_type(Number_type n)
  {
    _n = n;
    _exp=0;
    _set_vol(); 
    make_canonical();
  }  

  Protected_number_type(Number_type n, int exp, bool verb=false)
  {
    _n = n;
    _exp = exp;
    _set_vol();
    _verbose = verb; 
    make_canonical();
  }  

  void make_verbose()
  { _verbose = true;}

  void set_volume_constant(Number_type vol)
  {_vol = vol;}


/*
  void make_canonical()
  {
    if(_n == Number_type(0.0))
      return;

    while(std::abs(To_double()(_n))<Number_type(1.0))
    {
      if(_verbose)
        std::cout << " *, _n: " << _n << " , _exp: " << _exp << std::endl;

      _n = _n*Number_type(10.0);
      _exp = _exp - 1;
    }

    while(std::abs(To_double()(_n))>Number_type(10.0))
    {
      if(_verbose)
        std::cout << " /, _n: " << _n << " , _exp: " << _exp << std::endl;

      _n = _n/Number_type(10.0);
      _exp = _exp + 1;
    }    

  } // void make_canonical()
*/

  void make_canonical()
  {
    if(_n == Number_type(0.0))
      return;

    double appx = std::abs(To_double()(_n));

    if(appx<_vol)
    {
      Number_type k  = _vol/appx;
      int kexp = int(std::ceil(log10(k)));

      _exp -= kexp;
      _n = _n*std::pow(10.0, kexp);
    }
    else if(appx>_gap*_vol)
    {
      Number_type k  = appx/(_gap*_vol);
      int kexp = int(std::ceil(log10(k)));

      _exp += kexp;
      _n = _n/std::pow(10.0, kexp);
    }

  } // void make_canonical()



  Number_type n() const
  { return _n;}

  int exp() const
  { return _exp;}

  void set_n(Number_type n)
  { _n=n;}

  void set_exp(int exp)
  { _exp=exp;}

  Number_type to_number_type()
  { return _n*Power()(Number_type(10.0),_exp);} 

  Self operator+(const Self& pn) const
  {
    Self pn_l,pn_s;

    if(_n == 0)
      return pn;

    if(pn.n() == 0)
      return *this;

    if(pn.exp() > _exp)
    {
      pn_l = pn;
      pn_s = *this; 
    }
    else
    {
      pn_l = *this;
      pn_s = pn; 
    }

    Number_type nexp = Power()(Number_type(10.0), pn_l.exp()-pn_s.exp());

    Number_type n_s = pn_s.n()/nexp;

    return Self( pn_l.n() + n_s, pn_l.exp()); 

  } // Self operator+(const Self& pn) const

  Self operator+=(const Self& pn)
  { 
    Self pn_tmp = *this + pn; 

    _n = pn_tmp.n();
    _exp = pn_tmp.exp();

    return *this;
  }

  Self operator-(const Self& pn) const
  {
    Self pn_l,pn_s;

    if(_n == 0)
    {
      Self res(-pn.n(), pn.exp());

      return res;
    }

    if(pn.n() == 0)
      return *this;

    if(pn.exp() > _exp)
    {
      pn_l = pn;
      pn_s = *this; 
    }
    else
    {
      pn_l = *this;
      pn_s = pn; 
    }

    Number_type nexp = Power()(Number_type(10.0), pn_l.exp()-pn_s.exp());
    Number_type n_s = pn_s.n()/nexp;

    if(pn.exp() <= _exp)
      return Self( pn_l.n() - n_s, pn_l.exp()); 
    else
      return Self( -pn_l.n() + n_s, pn_l.exp()); 

  } // Self operator-(const Self& pn) const

  Self operator-=(const Self& pn)
  { 
    Self pn_tmp = *this - pn; 

    _n = pn_tmp.n();
    _exp = pn_tmp.exp();

    return *this;
  }

  Self operator*(const Self& pn) const
  {
    Self pn2( _n*pn.n(), _exp + pn.exp() );
    pn2.make_canonical();
  
    return pn2;

  } //   Self operator*(const Self& pn) const

  Self operator/(const Self& pn) const
  {
    Self pn2( _n/pn.n(), _exp - pn.exp());
    pn2.make_canonical();
  
    return pn2;
  }

  Self operator=(const Self pn)
  {
    _n = pn.n();
    _exp = pn.exp();
    _set_vol(); 
  
    return *this;
  } 

  Self operator=(const Number_type pn)
  {
    _n = n;
    _exp = 0;
    _set_vol(); 

    make_canonical();
  
    return *this;
  } 

  bool operator==(Self& pn) const
  { return (_n == pn.n() && _exp == pn.exp()); }

  bool operator==(const Self& pn) const
  { return (_n == pn.n() && _exp == pn.exp()); } 

  bool operator>(const Self& pn) const
  { 
    if(_exp > pn.exp())
      return true; 

    if(_exp < pn.exp())
      return false;

    if(_n > pn.n())
      return true;

    return false;
  } 

  bool operator<(const Self& pn) const
  { 
    if(_exp < pn.exp())
      return true; 

    if(_exp > pn.exp())
      return false;

    if(_n < pn.n())
      return true;

    return false;
  }  

 private:

  Number_type _n,_vol, _gap;
  int        _exp;
  bool _verbose;  

}; // class Protected_number_type

} // namespace FunctionalMeasures

#endif // PROTECTED_NUMBER_TYPE_H
