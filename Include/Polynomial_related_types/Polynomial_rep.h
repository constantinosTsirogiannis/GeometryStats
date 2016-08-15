#ifndef POLYNOMIAL_REP_H
#define POLYNOMIAL_REP_H

#include<vector>
#include<numeric>
#include<cassert>
#include<thread>

namespace FunctionalMeasures
{

template< class KernelType >
struct Polynomial_rep: public std::vector<typename KernelType::Protected_number_type>
{
 public:

  typedef KernelType                                  Kernel;
  typedef Polynomial_rep<Kernel>                      Self;
  typedef typename Kernel::Number_type               Number_type;
  typedef typename Kernel::Protected_number_type     Protected_number_type;
  typedef std::vector<Protected_number_type>          Container_type;
  typedef typename Container_type::iterator          Polynomial_iterator;
  typedef typename Container_type::reverse_iterator  Polynomial_reverse_iterator;
  typedef typename Kernel::Numeric_traits            Numeric_traits;

  Polynomial_rep():std::vector<Protected_number_type>(){}

}; // struct Polynomial_rep

template< class KernelType >
class Polynomial_multiplication
{
 public:

  typedef KernelType                                Kernel;
  typedef Polynomial_multiplication<Kernel>         Self;
  typedef typename Kernel::Polynomial              Polynomial;
  typedef typename Kernel::FFT                     FFT;
  typedef typename Kernel::Number_type             Number_type;
  typedef typename Kernel::Protected_number_type   Protected_number_type;
  typedef typename Kernel::Numeric_traits          Numeric_traits;

  typedef typename FFT::Complex                    Complex;

 private:

  struct Inner_product_functor_2
  {
    typedef typename Polynomial::Polynomial_iterator  Iterator;
    typedef typename Polynomial::Polynomial_reverse_iterator  Reverse_iterator;    
 
    Inner_product_functor_2(){}

    Inner_product_functor_2(Iterator a_first, Iterator a_last, Reverse_iterator b_r_first, 
                            Protected_number_type *c): _a_first(a_first), _a_last(a_last), _b_r_first(b_r_first), _c(c){}

   public:           
           
    void operator()(void) 
    {
      Protected_number_type cf(0.0); 
      *_c = std::inner_product(_a_first,_a_last,_b_r_first, cf); 
    }

   private:

    Iterator _a_first, _a_last; 
    Protected_number_type *_c; 
    Reverse_iterator _b_r_first; 

  }; // struct Inner_product_functor_2


  struct Inner_product_functor_3
  {     
    Inner_product_functor_3(){}

    Inner_product_functor_3(int i_from, int max, int step, Protected_number_type *a, 
                           Protected_number_type *b, Protected_number_type *c, int size_a, int size_b): 
                            _i_from(i_from), _max(max), _step(step), _a(a), _b(b), _c(c),
                            _size_a(size_a), _size_b(size_b){}

   public:           
           
    void operator()(void) 
    {
      //int count = 0;

      for(int i=_i_from; i<=_max; i+=_step)
      {
         Protected_number_type cf(Number_type(0.0));

         for(int j=0; j<=i && j<_size_a; j++)
           if(i-j<_size_b)
             cf = cf + (_a[j]*_b[i-j]);

         _c[i] = cf;

        //count++;

      } // for(int i=0; i<=threshold; i++)

      //std::cout << " Count: " << count << std::endl;

    }

   private:

    int _i_from, _step, _max, _size_a, _size_b; 
    Protected_number_type *_a, *_b, *_c; 

  }; // struct Inner_product_functor_3

 public:

  Polynomial_multiplication(bool use_naive=true):_use_naive(use_naive),_cutoff(-1),
                                                  _max_threads(-1), _is_parallel(false),
                                                  _count_in(0),_count_out(0),_count_total(0){}

  void set_cutoff_value( int cutoff)
  { _cutoff = cutoff; }

  void set_max_threads(int max_threads)
  { _max_threads = max_threads;}

  void set_is_parallel(int is_parallel)
  { _is_parallel = is_parallel;}

  void operator()(Polynomial &a, Polynomial &b, Polynomial &c)
  { 
    if(_use_naive == true) // || a.size() < 3 || b.size() < 3 )
    {
      naive_operator(a,b,c);

      if(a[0].n() == Number_type(0.0) || b[0].n() == Number_type(0.0))
        c[0] = Protected_number_type(Number_type(0.0));
    }
    else
      FFT_based_operator(a,b,c);

  } //   void operator()(Polynomial &a, Polynomial &b, Polynomial &c)

  void do_precomputations(int i)
  { 
    if(_use_naive == false)
      fft.precompute_low_precision_trigonometric_values(i);
  }

  void FFT_based_operator_inverse_padding(Polynomial &a, Polynomial &b, Polynomial &c)
  {
    std::vector<Complex> aa, bb, cc, y_a, y_b, y_c;

    int max_size = std::max(a.size(), b.size());
   
    int m=1;

    while(m<max_size)
      m = 2*m;

    m = 2*m;

    for(int i=0; i<m-a.size(); i++) 
      aa.push_back(Complex(Number_type(0.0),Number_type(0.0))); 

    for(int i=0; i<m-b.size(); i++) 
      bb.push_back(Complex(Number_type(0.0),Number_type(0.0)));

    for(int i=0; i<a.size(); i++) 
      aa.push_back(Complex(a[i],Number_type(0.0))); 

    for(int i=0; i<b.size(); i++) 
      bb.push_back(Complex(b[i],Number_type(0.0)));

    fft.transform(aa,y_a); 
    fft.transform(bb,y_b);  

    for(int i=0; i<y_a.size(); i++)
      y_c.push_back(y_a[i]*y_b[i]);

    fft.reverse_transform(y_c,cc);   
   
    for(int i=0; i<a.size()+b.size()-1; i++)
    {
      int index = cc.size()+i-a.size()-b.size();
      c.push_back(cc[index].real());
    }

  } // FFT_based_operator_inverse_padding(Polynomial &a, Polynomial &b, Polynomial &c)


  void FFT_based_operator(Polynomial &a, Polynomial &b, Polynomial &c)
  {
    std::vector<Complex> aa, bb, cc, y_a, y_b, y_c;

    for(int i=0; i<a.size(); i++) 
      aa.push_back(Complex(a[i],Number_type(0.0))); 

    for(int i=0; i<b.size(); i++) 
      bb.push_back(Complex(b[i],Number_type(0.0)));

    int max_size = std::max(aa.size(), bb.size());
   
    int m=1;

    while(m<max_size)
      m = 2*m;

    m = 2*m;

    if(a.size() < m/2)
      for(int i=aa.size()+1; i<=m/2; i++)
        aa.push_back(Complex(Number_type(0.0),Number_type(0.0)));

    if(b.size() < m/2)
      for(int i=bb.size()+1; i<=m/2; i++)
        bb.push_back(Complex(Number_type(0.0),Number_type(0.0)));

    for(int i=aa.size()+1; i<=m; i++)
      aa.push_back(Complex(Number_type(0.0),Number_type(0.0))); 

    for(int i=bb.size()+1; i<=m; i++)
      bb.push_back(Complex(Number_type(0.0),Number_type(0.0))); 

    fft.transform(aa,y_a); 
    fft.transform(bb,y_b);  

    for(int i=0; i<y_a.size(); i++)
      y_c.push_back(y_a[i]*y_b[i]);

    fft.reverse_transform(y_c,cc);   
   
    for(int i=0; i<a.size()+b.size()-1; i++)
      c.push_back(cc[i].real());

  } // FFT_based_operator(Polynomial &a, Polynomial &b, Polynomial &c)

  void naive_operator(Polynomial &a, Polynomial &b, Polynomial &c)
  {

    if(a.size()>b.size())
      naive_operator_1(b,a,c);
    else
      naive_operator_1(a,b,c);
  }

  void naive_operator_2(Polynomial &a, Polynomial &b, Polynomial &c)
  {
    int threshold = a.size()+b.size()-2;

    if(_cutoff > -1 && _cutoff < threshold)
      threshold = _cutoff;
      
    typename Polynomial::Polynomial_iterator a_first = a.begin(), a_last = ++a.begin();
    typename Polynomial::Polynomial_reverse_iterator b_r_first, b_r_last, tmp = b.rbegin();

    b_r_last = b.rbegin();

    while(tmp != b.rend())
    {
      tmp++;

      if(tmp != b.rend())
        b_r_last = tmp;
    }

    b_r_first = b_r_last;
    b_r_last = b.rend();

    int count_a=1, 
         count_b=1;

    Protected_number_type cf(Number_type(0.0)),af;
  
    af= std::inner_product(a_first,a_last,b_r_first,cf);
        

    c.assign(threshold+1, cf);
    c[0] = af;

    if(_is_parallel== false || _max_threads<=0 || a.size()*b.size() > 10000)    
      for(int i=1; i<=threshold; i++)
      {
        if(a_last != a.end())
          a_last++;

        if(b_r_first !=  b.rbegin())
          b_r_first--;
        else
          a_first++;


        if(a_first != a_last && b_r_first != b_r_last)
        {
          af = std::inner_product(a_first,a_last,b_r_first,cf);
          c[i] = af;
        }
        else 
          break;

      } // for(int i=0; i<a.size()+b.size()-1; i++)
    else // of if(_max_threads<=0) 
    {
      Protected_number_type cc[c.size()];

      cc[0] = c[0];

      //int cores = std::thread::hardware_concurrency();

      //typename Polynomial::Polynomial_iterator cit = c.begin();
 
      for(int i=1; i<=threshold; )
      {
        int j_threshold = threshold-i;

        if(_max_threads-1 < j_threshold)
          j_threshold = _max_threads-1;

        std::vector<std::thread> threads;

        for(int j=0; j<=j_threshold; j++)
        {
          //cit++;

          if(a_last != a.end())
            a_last++;

          if(b_r_first !=  b.rbegin())
            b_r_first--;
          else
            a_first++;

          Inner_product_functor_2 ipf(a_first,a_last,b_r_first,&cc[i]);
          //Inner_product_functor_2 ipf;

          i++;

          if(a_first != a_last && b_r_first != b_r_last)
            threads.push_back(std::thread(ipf));
          else
          {
            j=j_threshold+1; 
            i=threshold+1;
          }

        } // for(int j=i; j<=threshold; j++)

        for(int j=0; j<threads.size(); j++)
          threads[j].join();

      } // for(int i=1; i<=threshold; i++)


      //std::cout << " Finished "  << a.size() << std::endl;

      for(int i=1; i<=threshold; i++)
        c[i] = cc[i];
 
    } // else of if(_max_threads<=0) 
  


  } // naive_operator_2(Polynomial &a, Polynomial &b, Polynomial &c)

  int count_total()
  { return _count_total;}

  int count_in()
  { return _count_in;}

  int count_out()
  { return _count_out;}

  void naive_operator_1(Polynomial &a, Polynomial &b, Polynomial &c)
  {
    int threshold = a.size()+b.size()-2;

    if(_cutoff > -1 && _cutoff < threshold)
      threshold = _cutoff;

    int mult_complexity = std::min(a.size(),b.size());      

    _count_total++;

    if(mult_complexity < 500)
      _count_out++;
    else
      _count_in++;

    if(_is_parallel== false || _max_threads<=0 || mult_complexity < 500) 
    {
      for(int i=0; i<=threshold; i++)
      {
         Protected_number_type cf(Number_type(0.0));

         for(int j=0; j<=i && j<a.size(); j++)
           if(i-j<b.size())
             cf = cf + (a[j]*b[i-j]);


         c.push_back(cf);

      } // for(int i=0; i<=threshold; i++)
    }
    else // of if(_is_parallel == false ... )
    {
      Protected_number_type *aa = new Protected_number_type[a.size()], 
                            *bb = new Protected_number_type[b.size()], 
                            *cc = new Protected_number_type[threshold+1];
     
      for(int hh=0; hh<a.size(); hh++)
        aa[hh] = a[hh];

      for(int hh=0; hh<b.size(); hh++)
        bb[hh] = b[hh];


      int full_reps = (threshold+1)/_max_threads,
           extra_rep = (threshold+1)%_max_threads;

      int i=0;

      std::vector<std::thread> threads;

      for(int jj=0; jj<_max_threads; jj++ )
      {
        Inner_product_functor_3 ipf3(jj,threshold,_max_threads,aa, bb, cc, int(a.size()), int(b.size()));
        threads.push_back(std::thread(ipf3));

      }  // for(int jj=0; jj<_max_threads; jj++ )

      for(int jj=0; jj<threads.size(); jj++ )
        threads[jj].join();

      for(int hh=0; hh<=threshold; hh++)
        c.push_back(cc[hh]);

      delete[] aa;  
      delete[] bb;  
      delete[] cc;  

    } // // of if(_is_parallel == false ... )

  } // naive_operator_1(Polynomial &a, Polynomial &b, Polynomial &c)


  void naive_operator_0(Polynomial &a, Polynomial &b, Polynomial &c)
  {
    int threshold = a.size()+b.size()-2;

    if(_cutoff > -1 && _cutoff < threshold)
      threshold = _cutoff;
      

    for(int i=0; i<=threshold; i++)
    {
       Protected_number_type cf(Number_type(0.0));

       for(int j=0; j<=i; j++)
       {
         if(j<a.size() && i-j<b.size())
           cf = cf + (a[j]*b[i-j]);
       }  

       c.push_back(cf);

    } // for(int i=0; i<a.size()+b.size()-1; i++)

  } // naive_operator_0(Polynomial &a, Polynomial &b, Polynomial &c)

  Protected_number_type 
  compute_coefficient(Polynomial &a, Polynomial &b, int i)
  {
    assert(i < a.size()+b.size()-1);

    Protected_number_type cf(Number_type(0.0));    

    for(int j=0; j<=i; j++)
    {
      if(j<a.size() && i-j<b.size())
        cf = cf + (a[j]*b[i-j]);
    }

    return cf;

  } // compute_coefficient(Polynomial &a, Polynomial &b, int i)


  bool test_FFT(int max_d=20, int repeat = 1000, int seed=4)
  {
    srand(seed);
 
    for(int k=0; k<repeat; k++)
    {
      Polynomial a,b, c_FFT, c_slow;

      int size_a = rand()%(max_d+1),
          size_b = rand()%(max_d+1);

      for(int i=0; i<size_a+1; i++)
        a.push_back(Number_type(rand()%(max_d+1))/Number_type(1.0));

      for(int i=0; i<size_b+1; i++)
        b.push_back(Number_type(rand()%(max_d+1))/Number_type(1.0));

      naive_operator(a,b,c_slow);
      FFT_based_operator_inverse_padding(a,b,c_FFT);


      if(c_FFT.size() != c_slow.size())
      {
        std::cout << " The output polynomials of the FFT and the slow operator" << std::endl
                  << " have a different number of coefficients." << std::endl; 
        std::cout << " The output of the FFT has " << c_FFT.size() << " coefficients. " << std::endl;
        std::cout << " The output of the slow operator has " << c_slow.size() << " coefficients. " << std::endl;
        return false;
      }

      for(int i=0; i<c_FFT.size(); i++)
        if( std::abs(c_FFT[i].n() - c_slow[i].n()) > 0.0001 || std::abs(c_FFT[i].exp() - c_slow[i].exp()) != 0  )
        {
          std::cout << " The output polynomials of the FFT and the slow operator" << std::endl
                    << " differ in at least one coefficient." << std::endl; 
          std::cout << " This is the coefficient of the " << i << "-th power of x" << std::endl; 
          std::cout << " The value of this coefficient computed with the FFT is " << c_FFT[i].to_double() << std::endl;
          std::cout << " The value of this coefficient computed with the slow operator is " 
                    << c_slow[i].to_double() << std::endl;

          std::cout << " The two input polynomials are: " << std::endl;

          std::cout << " Polynomial A: ";

          for(int j=0; j<a.size(); j++ )
            std::cout << a[j].to_double() << " || " ;
   
          std::cout << std::endl << std::endl;

          std::cout << " Polynomial B: ";

          for(int j=0; j<b.size(); j++ )
            std::cout << b[j].to_double() << " || " ;
   
          std::cout << std::endl << std::endl;

          std::cout << " The two output polynomials are: " << std::endl;

          std::cout << " FFT: ";

          for(int j=0; j<c_FFT.size(); j++ )
            std::cout << c_FFT[j].to_double() << " || " ;
   
          std::cout << std::endl << std::endl;

          std::cout << " Slow operator: ";

          for(int j=0; j<c_slow.size(); j++ )
            std::cout << c_slow[j].to_double() << " || " ;
   
          std::cout << std::endl << std::endl;

          return false;
        }   

    } // for(int k=0; k<repeat; k++)
    
    
    return true;

  } // test_FFT(int max_d=20, int repeat = 1000)
  
 private:

  FFT fft;
  bool _use_naive, _is_parallel;
  int _cutoff, _max_threads, _count_in, _count_out, _count_total;

}; // class Polynomial_multiplication

} // namespace FunctionalMeasures

#endif // POLYNOMIAL_REP
