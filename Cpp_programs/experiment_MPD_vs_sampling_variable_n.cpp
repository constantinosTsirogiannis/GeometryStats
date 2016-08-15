#include<fstream>
#include<ctime>
#include<sys/time.h>
#include<Functional_measures_kernel.h>


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~ List of all types ~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//


typedef FunctionalMeasures::Numeric_traits_double        Numeric_traits;
typedef Functional_measures_kernel<Numeric_traits>       Kernel;
typedef Kernel::Number_type                              Number_type;
typedef Kernel::Point                                    Point;
typedef Kernel::Construct_point_set                      Construct_point_set;
typedef Kernel::Calculate_moments_with_sampling          Calculate_moments_with_sampling;
typedef Kernel::Distribution_type                        Distribution_type;

int diff_ms(timeval t1, timeval t2);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~ List of measures ~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

typedef Kernel::Mean_pairwise_distance                   Mean_pairwise_distance;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~ Essential Variables ~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

int minimum_n_artificial_data = 500,
    maximum_n_artificial_data = 10000,
    minimum_n_real_data = 465,
    maximum_n_real_data = 7965,
    step = 500;

int dimension = 3;
int number_of_samples = 1000;
int random_seed = 4;

// Sample size s is calculated as half of the given total size n.  

int calculate_s( int n)
{ return n/2; }

bool conduct_experiments_on_artificial_data = false,
      conduct_experiments_on_real_data = true; 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~ Input and Output Filenames ~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

char input_real_data[] = "../Data_sets/PCOA_7965_species_coords_dim_20.csv",
      output_precision_artificial_data[] = "Output/MPD_artificial_data_precision_variable_n.txt",
      output_time_artificial_data[] = "Output/MPD_artificial_data_time_variable_n.txt",
      output_precision_real_data[] = "Output/MPD_real_data_precision_variable_n.txt",
      output_time_real_data[] = "Output/MPD_real_data_time_variable_n.txt";


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

int main(int argc, char *argv[])
{

  std::vector<Point> real_points, artificial_points; 

  double seconds_approx, seconds_exact, 
          expectation_approx, expectation_exact, 
          variance_approx, variance_exact;

  timeval time_start, time_end;

  Construct_point_set point_set_constructor;

  point_set_constructor(input_real_data, std::back_inserter(real_points), dimension);

  point_set_constructor.construct_random_point_set(maximum_n_artificial_data, dimension,
                                                   std::back_inserter(artificial_points), random_seed);


  if(conduct_experiments_on_artificial_data == true)
  {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //~~~~~~~~~~ Experiments on artificial data ~~~~~~~~~~~//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    std::ofstream out_time_art(output_time_artificial_data),
                  out_prec_art(output_precision_artificial_data);

    std::cout << std::endl;
    std::cout << " Commencing experiments on artificial data " << std::endl << std::endl;

    for(int n=minimum_n_artificial_data; n<=maximum_n_artificial_data; n+=step )
    {
      std::cout << " [ Executing computations with n=" << n << " ]" << std::endl;

      std::vector<Point>  points;

      int s = calculate_s(n);

      for(int k=0; k<n; k++)
        points.push_back(artificial_points[k]);

      Mean_pairwise_distance  mpd_exact;
      Calculate_moments_with_sampling approx_calc;

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      gettimeofday(&time_start, NULL);

      mpd_exact.load_point_set(points.begin(), points.end());
      expectation_exact = mpd_exact.compute_expectation(s);
      variance_exact = mpd_exact.compute_variance(s);

      // Compute moments for all sample sizes, because we can.
      for(int j=0; j<=n; j++)
      {
        mpd_exact.compute_expectation(j);
        mpd_exact.compute_variance(j);
      }

      gettimeofday(&time_end, NULL);
      seconds_exact = diff_ms(time_end,time_start)/double(1000.0); 

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

      gettimeofday(&time_start, NULL);

      approx_calc.load_point_set(points.begin(), points.end());
      std::pair<double, double> pr = approx_calc(s,number_of_samples,mpd_exact);
      expectation_approx = pr.first;
      variance_approx = pr.second;

      gettimeofday(&time_end, NULL);
      seconds_approx = diff_ms(time_end,time_start)/double(1000.0); 

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
      out_time_art << n << " " << seconds_exact << " " << seconds_approx << std::endl;
      out_prec_art << n << " " << ((std::abs(expectation_exact-expectation_approx)/expectation_exact)*Number_type(100.0))
                        << " " << ((std::abs( std::sqrt(variance_exact)-std::sqrt(variance_approx))/
                                    std::sqrt(variance_exact))*Number_type(100.0)) << std::endl;    
 
    } // for(int n=minimum_n_artificial_data; n<=maximum_n_artificial_data; n+=step )

    out_time_art.close();
    out_prec_art.close();

  } // if(conduct_experiments_on_artificial_data == true)

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  if(conduct_experiments_on_real_data == true)
  {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //~~~~~~~~~~~~~ Experiments on real data ~~~~~~~~~~~~~~//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

    std::ofstream out_time_real(output_time_real_data),
                  out_prec_real(output_precision_real_data);

    std::cout << std::endl;
    std::cout << " Commencing experiments on real data " << std::endl << std::endl;

    for(int n=minimum_n_real_data; n<=maximum_n_real_data; n+=step )
    {
      std::cout << " [ Executing computations with n=" << n << " ]" << std::endl;

      std::vector<Point>  points;

      int s = calculate_s(n);

      for(int k=0; k<n; k++)
        points.push_back(real_points[k]);

      Mean_pairwise_distance  mpd_exact;
      Calculate_moments_with_sampling approx_calc;

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
      gettimeofday(&time_start, NULL);

      mpd_exact.load_point_set(points.begin(), points.end());
      expectation_exact = mpd_exact.compute_expectation(s);
      variance_exact = mpd_exact.compute_variance(s);

      // Compute moments for all sample sizes, because we can.
      for(int j=0; j<=n; j++)
      {
        mpd_exact.compute_expectation(j);
        mpd_exact.compute_variance(j);
      }

      gettimeofday(&time_end, NULL);
      seconds_exact = diff_ms(time_end,time_start)/double(1000.0); 

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

      gettimeofday(&time_start, NULL);

      approx_calc.load_point_set(points.begin(), points.end());
      std::pair<double, double> pr = approx_calc(s,number_of_samples,mpd_exact);
      expectation_approx = pr.first;
      variance_approx = pr.second;

      gettimeofday(&time_end, NULL);
      seconds_approx = diff_ms(time_end,time_start)/double(1000.0); 

      //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
      out_time_real << n << " " << seconds_exact << " " << seconds_approx << std::endl;
      out_prec_real << n << " " << ((std::abs(expectation_exact-expectation_approx)/expectation_exact)*Number_type(100.0))
                        << " " << ((std::abs( std::sqrt(variance_exact)-std::sqrt(variance_approx))/
                                    std::sqrt(variance_exact))*Number_type(100.0)) << std::endl;    
 
    } // for(int n=minimum_n_real_data; n<=maximum_n_real_data; n+=step )

    out_time_real.close();
    out_prec_real.close();

  } // if(conduct_experiments_on_real_data == true)

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

  return 0;

} // main(...)


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~ Function calculating time difference in seconds ~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

int diff_ms(timeval t1, timeval t2)
{
    return (((t1.tv_sec - t2.tv_sec) * 1000000) + 
            (t1.tv_usec - t2.tv_usec))/1000;
}
