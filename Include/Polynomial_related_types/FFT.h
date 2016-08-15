#ifndef FFT_H
#define FFT_H

#include<vector>
#include<complex>

namespace FunctionalMeasures
{

template< class KernelType >
struct FFT
{
 public:

  typedef KernelType                                       Kernel;
  typedef FFT<Kernel>                                      Self;
  typedef typename Kernel::Number_type                    Number_type;
  typedef typename Kernel::Numeric_traits                 Numeric_traits;
  typedef typename Numeric_traits::Cosine                 Cosine;
  typedef typename Numeric_traits::Sine                   Sine;
  typedef typename Numeric_traits::To_double              To_double;
  typedef typename Numeric_traits::Protected_number_type  Protected_number_type; 

  typedef std::complex<Protected_number_type>              Complex;

  typedef typename Kernel::Exception_type          Exception_type;
  typedef typename Kernel::Exception_functor       Exception_functor;   


 public:

  FFT(){}

  template <class OutputIterator>
  void bit_reverse(int n, OutputIterator ot)
  {
    if(n<0)
      return;

    int m=1, current=0;
    std::vector<bool> binary_rep;
    std::vector<int> components;
    std::vector<int> reverse;

    while(m<n)
      m = 2*m; 

    if(m>=n)
      m = m/2;

    while(m>1)
    {
      binary_rep.push_back(false);
      components.push_back(m);
      m=m/2;  
    }

    binary_rep.push_back(false);
    components.push_back(m);

    *ot++ = 0;

    for(int i=1; i<=n; i++)
    {
      int index=0;      

      while(binary_rep[index]==true && index < binary_rep.size()-1)
      {
        binary_rep[index]=false;
        current = current - components[index];
        index++;
      }

      binary_rep[index]=true;
      current = current + components[index];

      *ot++ = current;

    } // for(int i=1; i<=n; i++)

  } // bit_reverse(int n, OutputIterator ot)

  // Each of the next three functions computes the trigonometric 
  // values for angles 2*pi*j/m where m,j are positive integers, 
  // m is a power of two, m>4, and j/m <=1/8 .
  // 
  // The first function is the trivial one that calls
  // directly the built-in functions of sine and cosine.
  //
  // The second function computes first the sine and cosine of
  // 2*pi/m and then computes the rest of the values using 
  // the addition formulas for sines/cosines.
  //
  // The third function does the computation by evaluating first
  // the sine and cosine at all divisors of m, using the halving
  // formulas of these functions. Then, the trigonometric
  // values for the rest of the angles are calculated by applying the
  // addition/subtraction formulas on the formerly computed values.

  template<class OutputIterator>
  void compute_trigonometric_values_trivial(int m, bool is_reverse, OutputIterator ot)
  {
      Cosine cosine;
      Sine sine;

      int mhalf = m/2;   
      Number_type pi = double(4.0)*std::atan(double(1.0)); 

      *ot++ = std::make_pair(Number_type(1.0), Number_type(0.0));
 
      for(int i=1; i< m/8; i++)
      {
        Number_type cosm = cosine(Number_type(i)*pi/Number_type(mhalf)), 
                    sinm = sine(Number_type(i)*pi/Number_type(mhalf));

        if(is_reverse == false)
          *ot++ = std::make_pair(cosm, sinm);
        else
          *ot++ = std::make_pair(cosm, -sinm);
      }

      if(is_reverse == false)
        *ot++ = std::make_pair( Number_type(sqrt(double(2.0))/double(2.0)), 
                                Number_type(sqrt(double(2.0))/double(2.0)));
      else
        *ot++ = std::make_pair( Number_type(sqrt(double(2.0))/double(2.0)), 
                                Number_type(-sqrt(double(2.0))/double(2.0)));

  } // compute_trigonometric_values_trivial(int m, bool is_reverse, OutputIterator ot)


  template<class OutputIterator>
  void compute_trigonometric_values_incremental(int m, bool is_reverse, OutputIterator ot)
  {
    Cosine cosine;
    Sine sine;

    int mhalf = m/2;   
    Number_type pi = double(4.0)*std::atan(double(1.0)); 
    Number_type cosm = cosine(pi/Number_type(mhalf)), 
                sinm = sine(pi/Number_type(mhalf));

    Number_type cos_curr = cosm, 
                sin_curr = sinm;

    *ot++ = std::make_pair(Number_type(1.0), Number_type(0.0));

    for(int i=1; i<m/8; i++)
    { 
      if(is_reverse == false)
        *ot++ = std::make_pair(cos_curr, sin_curr);
      else
        *ot++ = std::make_pair(cos_curr, -sin_curr);  
  
      sin_curr = (sin_curr*cosm) + (cos_curr*sinm); 
      cos_curr = sqrt(Number_type(1.0) - (sin_curr*sin_curr) );
    }

    if(is_reverse == false)
      *ot++ = std::make_pair(Number_type(sqrt(double(2.0))/Number_type(2.0)),
                             Number_type(sqrt(double(2.0))/Number_type(2.0)));
    else
      *ot++ = std::make_pair(Number_type(sqrt(double(2.0))/Number_type(2.0)),
                             Number_type(-sqrt(double(2.0))/Number_type(2.0)));    

  } // compute_trigonometric_values_incremental(...)


  template<class OutputIterator>
  void compute_trigonometric_values_halving(int m, bool is_reverse, OutputIterator ot)
  {
    Cosine cosine;
    Sine sine;
    To_double to_double;

    std::vector< std::pair<Number_type, Number_type> > res;

    res.assign((m/8)+1, std::make_pair(Number_type(0.0),Number_type(0.0)));

    res[0] = std::make_pair(Number_type(1.0), Number_type(0.0));

    int k=m/8;

    Number_type cs = Number_type(sqrt(double(2.0))/Number_type(2.0)),
                sn = Number_type(sqrt(double(2.0))/Number_type(2.0));

    while(k>=1)
    {
      if(is_reverse == false)
        res[k] = std::make_pair(cs, sn);
      else
        res[k] = std::make_pair(cs, -sn);      

      if(k==1)
        k=0;
      else
      {
        k=k/2;

        Number_type old_cs = cs; 
        cs = Number_type(sqrt(double(2.0)+(double(2.0)*to_double(old_cs)))/Number_type(2.0)); 
        sn = Number_type(sqrt(double(2.0)-(double(2.0)*to_double(old_cs)))/Number_type(2.0)); 
      }

    } // while(k>=1)


    ////////////////////////////////////////////////////////////////

    int curr_pow = 1;

    for(int i=1; i<m/8; i++)
    {
      if(i>curr_pow)
        curr_pow = curr_pow*2;

      if(i != curr_pow)
      { 
        int rem = curr_pow - i; 

        Number_type sin_curr = (res[curr_pow].second*res[rem].first) - (res[curr_pow].first*res[rem].second), 
                    cos_curr = (res[curr_pow].first*res[rem].first) + (res[curr_pow].second*res[rem].second);

        if(is_reverse == false)
          res[i] = std::make_pair(cos_curr, sin_curr);
        else
          res[i] = std::make_pair(cos_curr, -sin_curr);  

      } // if(res[i].first == Number_type(0.0))

    } // for(int i=1; i<m/8; i++)

    for(int i=0; i<res.size(); i++)
      *ot++ = res[i];

  } // compute_trigonometric_values_halving(...)

  void precompute_low_precision_trigonometric_values(int i)
  {
    _trig.clear();

    if(i==0)
      return;

    _trig.clear();

    int m=1;

    while(m<i)
      m= 2*m;

    m = 4*m;

    Cosine cosine;
    Sine sine;

    int mhalf = m/2;   
    Number_type pi = double(4.0)*std::atan(double(1.0)); 

    _trig.push_back(std::make_pair(Number_type(1.0), Number_type(0.0)));

    for(int i=1; i<m/8; i++)
    { 
      Number_type cosm = cosine(Number_type(i)*pi/Number_type(mhalf)), 
                  sinm = sine(Number_type(i)*pi/Number_type(mhalf));

      _trig.push_back(std::make_pair(cosm, sinm));
    }

    _trig.push_back( std::make_pair(Number_type(sqrt(double(2.0))/Number_type(2.0)),
                                    Number_type(sqrt(double(2.0))/Number_type(2.0))));

  } // compute_low_precision_trigonometric_values(int m)

  /*
  void precompute_high_precision_trigonometric_values(int i)
  {
    std::cout << " Starting precomputations " << std::endl;

    if(i==0)
      return;

    _trig.clear();

    int m=1;

    while(m<i)
      m= 2*m;

    m = 4*m;
   
    namespace mp = boost::multiprecision;     // Reduce the typing a bit later...
    typedef mp::number<mp::mpfr_float_backend<6000> >  Mpfloat;    

    int mhalf = m/2;

    Mpfloat f= double(4.0), 
            o = double(1.0),
            t = 2.0,
            mh = mhalf;

    Mpfloat pi = f*atan(o);

    _trig.push_back(std::make_pair(Number_type(1.0), Number_type(0.0))); 

    for(int i=1; i<m/8; i++)
    { 
      Mpfloat ii=i/mh; 

      Mpfloat cosm = cos(ii*pi), 
              sinm = sin(ii*pi);

      _trig.push_back(std::make_pair(Number_type(cosm.convert_to<double>()), 
                                     Number_type(sinm.convert_to<double>()))); 
    }

    Mpfloat fth = sqrt(t)/t;

    _trig.push_back(std::make_pair(Number_type(fth.convert_to<double>()), 
                                   Number_type(fth.convert_to<double>())));

    std::cout << " Done with precomputations " << std::endl;

  } // compute_high_precision_trigonometric_values(int m)
  */

  template<class OutputIterator>
  void extract_precomputed_trigonometric_values(int m, bool is_reverse, OutputIterator ot)
  {
    if(_trig.size() == 0)
    {
      compute_trigonometric_values_trivial(m,is_reverse,ot);
      return;
    }

    int fc = 8*(_trig.size()-1)/m;

    *ot++ = _trig[0];

    if(fc !=0)
      for(int i=fc; i< _trig.size()-1; i+=fc)
      {
        if(is_reverse == false)
          *ot++ = _trig[i];
        else
          *ot++ = std::make_pair(_trig[i].first, -_trig[i].second);
      }
    
    if(is_reverse == false)
      *ot++ = _trig[_trig.size()-1];
    else
      *ot++ = std::make_pair(_trig[_trig.size()-1].first, -_trig[_trig.size()-1].second);
  }


  template<class OutputIterator>
  void compute_complex_roots_of_unity(int m, bool is_reverse, OutputIterator ot)
  {
    if(m==2)
    {
      *ot++ = Complex(Number_type(1.0), Number_type(0.0));
      *ot++ = Complex(Number_type(-1.0), Number_type(0.0));
      return;
    }
   
    if(m==4)
    {
      *ot++ = Complex(Number_type(1.0), Number_type(0.0));

      if(is_reverse == false)
        *ot++ = Complex(Number_type(0.0), Number_type(1.0));
      else
        *ot++ = Complex(Number_type(0.0), Number_type(-1.0));

      *ot++ = Complex(Number_type(-1.0), Number_type(0.0));
      return;
    }

////////////////////////////
////////////////////////////////
///////////////////////////////////////////
/////////////////////////////
////////////////////

    std::vector< std::pair<Number_type, Number_type> > res, res2;

    extract_precomputed_trigonometric_values(m,is_reverse,std::back_inserter(res));
    //compute_trigonometric_values_trivial(m, is_reverse, std::back_inserter(res));

    if(is_reverse == false)
      for(int i=1; i<m/8; i++)
        res.push_back( std::make_pair(res[(m/8)-i].second,res[(m/8)-i].first) );
    else  
      for(int i=1; i<m/8; i++)
        res.push_back( std::make_pair(-res[(m/8)-i].second,
                                      -res[(m/8)-i].first) );

    if(is_reverse == false)
      res.push_back(std::make_pair(Number_type(0.0), Number_type(1.0)));
    else
      res.push_back(std::make_pair(Number_type(0.0), Number_type(-1.0)));


    for(int i=1; i<m/4; i++)
      res.push_back( std::make_pair(-res[(m/4)-i].first,res[(m/4)-i].second) );


    res.push_back(std::make_pair(Number_type(-1.0), Number_type(0.0)));

    for(int i=0; i<res.size(); i++)  
      *ot++ = Complex(res[i].first, res[i].second);


  } // compute_complex_roots_of_unity(int m, OutputIterator ot)

  void transform(std::vector<Complex> &a,  std::vector<Complex> &res)
  { generic_transform(a, res, false); }

  void reverse_transform(std::vector<Complex> &a,  std::vector<Complex> &res)  
  { generic_transform(a, res, true); }


  void generic_transform(std::vector<Complex> &a, std::vector<Complex> &res, bool is_reverse)
  {
    if(a.size() == 0)
      return;

    Protected_number_type sn;

    if(is_reverse == true)
      sn = Protected_number_type(Number_type(-1.0));
    else
      sn = Protected_number_type(Number_type(1.0));

    // Add padding and arrange coefficients according to bit-reverse order 

    int m=1;

    while(m<a.size())
     m = m*2;

    //for(int i=a.size()+1; i<=m; i++)
    //  a.push_back(Complex(0,0));

    for(int i=0; i<m-a.size(); i++)
      a.push_back(a[a.size()-i-2]);

    std::vector<int> br_order;

    bit_reverse(a.size(), std::back_inserter(br_order));

    res.assign(a.size(),Complex(Number_type(0.0),Number_type(0.0)));

    Cosine cosine;
    Sine sine;

    Number_type pi = double(4.0)*std::atan(double(1.0)); 

    for(int i=0; i<a.size(); i++)
      res[br_order[i]] = a[i];    

    for(m=2; m<=res.size(); m=m*2)
    {
      //Complex omega_m( cosine(Number_type(2.0)*pi/Number_type(m)),
      //                 sine(Number_type(2.0)*pi/Number_type(m))   );

      std::vector<Complex> omegas;
      compute_complex_roots_of_unity(m,is_reverse,std::back_inserter(omegas));

      for(int k = 0; k<=res.size()-1; k+=m)
      {
        Complex omega(Number_type(1.0), Number_type(0.0));
        
        for(int j=0; j<=(m/2)-1; j++)
        {
          Complex t = omegas[j]*res[k+j+(m/2)],
                  u = res[k+j];
         
          res[k+j] = u+t;
          res[k+j+(m/2)] = u-t;

          Protected_number_type cs(cosine(pi*Number_type(Number_type((2*j)+2)/Number_type(m)))),
                                sne(sn*Protected_number_type(sine(pi*Number_type(Number_type((2*j)+2)/Number_type(m)))));

          omega = Complex( cs, sne);

        } // for(int j=0; j<=(m/2)-1; j++)       

      } // for(int k = 0; k<=n-1; k+=m)

    } // for(m=1; m<=res.size(); m=m*2)

    if(is_reverse)
      for(int i=0; i<res.size(); i++)
        res[i] = Complex(res[i].real()/Protected_number_type(Number_type(res.size())),
                         res[i].imag()/Protected_number_type(Number_type(res.size())));

  } // FFT_generic(std::vector<Complex> &a, bool is_reverse)

 private:

  std::vector<std::pair<Number_type, Number_type> > _trig;

}; // class FFT

} // namespace FunctionalMeasures

#endif // FFT_H

