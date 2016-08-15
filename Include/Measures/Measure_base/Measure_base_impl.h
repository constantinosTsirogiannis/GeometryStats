
#ifndef MEASURE_BASE_IMPL_H
#define MEASURE_BASE_IMPL_H

#include<string>
#include<vector>
#include<iostream>
#include<fstream>
#include<cstdlib>

  template < class KernelType >
  template < class OutputIterator >
  void FunctionalMeasures::Measure_base<KernelType>::
  _read_sample_sizes_from_file(char *filename, std::vector<Point> &points, OutputIterator ot)
  {
    std::ifstream in(filename);
    std::vector<int> sample_sizes;

    // Reading first file with queries
    if( !( in.is_open() && in.good() ) )
    {
      std::string exception_msg(" There was a problem with opening the file with the species samples.\n");
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    std::string line;
    std::getline(in,line);
 
    int prev_index=-1, current_index =0;

    while( current_index< line.size()-1 )
    {
      do
      {
        current_index++;
      }
      while(line[current_index]!=',' && current_index < line.size() );

      if(current_index -prev_index < 2)
      {
        std::string exception_msg(" There is a mistake in the syntax of the sample sizes file.\n");
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      } 

      std::string substring = line.substr(prev_index+1,current_index -prev_index-1);

      for(int i=0; i<substring.size(); i++)
        if(isdigit(substring[i]) == false)
        {
          std::string exception_msg;
          exception_msg += " There is an error in the syntax of the sample sizes file.\n";     
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

      int size = atoi(substring.c_str());

      if(size < 0 || size > points.size() )
      {
        std::string exception_msg(" One of the sample sizes in the file is out of range.\n");
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      } 

      sample_sizes.push_back(size);

      prev_index = current_index;

    } // while( current_index< line.size()-1 )

    for(int i=0; i<sample_sizes.size(); i++)
      *ot++ = sample_sizes[i];

  } // _read_sample_sizes_from_file(...)


  // Input:  A point set and a range of iterators that indicate a list of species names 
  //         (in std::string format) in the point set.
  // Output: The value of the current measure for this set of species.

  template < class KernelType >
  template < class RangeIterator, class Measure >    
  typename KernelType::Number_type 
  FunctionalMeasures::Measure_base<KernelType>::
  _list_query(std::vector<Point> &points, RangeIterator rbegin, RangeIterator rend, Measure &msr)
  {
    RangeIterator it;
    std::string str;
    std::set<std::string> in_names;
    std::vector<int> point_indices;

    for( it = rbegin; it != rend; it++ )
      in_names.insert(*it);

    // Find the indices of all the points that correspond to the query species.
    for( int i=0; i<points.size(); i++ )
      if( in_names.find(points[i].taxon) != in_names.end() )
           point_indices.push_back(i);

    return msr(point_indices.begin(), point_indices.end());
	
  } // list_query(std::vector<Point> &points, RangeIterator rbegin, RangeIterator rend)
    
  // Input: A csv file which stores a matrix where each column corresponds to a species of the point set
  // and each row indicates a sample of these species for we want to compute the
  // distance measure. A certain species is considered as part of the i-th sample if in the i-th row 
  // of the matrix there is a '1' at the column that corresponds to this species (otherwise a '0').

  // The second argument is an output iterator of the type std::back_insert_iterator<std::vector< Number_type > >.
  // Output: A vector of numbers (passed in the form of an outputiterator), where the i-th number
  // is the value of the measure for the sample that is described in the i-th row of the input matrix.
  
  template < class KernelType >  
  template < class Measure, class OutputIterator>
  int FunctionalMeasures::Measure_base<KernelType>::
  _csv_matrix_query( std::vector<Point> &points, char *filename, Measure &msr, bool standardised, OutputIterator ot )
  {
    std::ifstream in(filename);

    if( !( in.is_open() && in.good() ) )
    {
      std::string exception_msg(" There was a problem with opening the file with the matrix.\n");
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    // Read the first row, the one that contains the species names

    std::vector<std::string> names;
    std::vector<int>   column_to_node_vec;
    std::string line;
    std::vector<Number_type> mean_values, deviation_values;

    mean_values.assign(points.size()+1, Number_type(-1.0));
    deviation_values.assign(points.size()+1, Number_type(-1.0));

    char a;

    std::getline(in,line);

    int c=0;

    while(c < line.size() )
    {
      std::string str;
      a = line[c];

      while(a != ',' && a != ' ' && a != '\r' && c<line.size() )
      {
        str.push_back(a);
        c++;
        a = line[c];
      }

      if( str.size() > 0 )
        names.push_back(str);

      c++;
    }

    if(names.size() < points.size())
    {
      std::string warning(" Warning: the input matrix has fewer columns than the number of species in the point set.");
      Exception_functor().issue_warning(warning);
    }

    std::vector<bool> checked_names;
    checked_names.assign(points.size(),false);
    std::map<std::string, int> point_taxa;

    for(int i=0; i<points.size(); i++)
      point_taxa[points[i].taxon] = i;

    for( int i=0; i<names.size(); i++ )
    {
      typename std::map<std::string, int>::iterator pt_it = point_taxa.find(names[i]);

      if( pt_it == point_taxa.end() )
      {
        std::string exception_msg;
        exception_msg += " One of the species names in the input matrix was not found in the point set (";
        exception_msg += names[i];
        exception_msg += ") \n";
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }
      else
      { 
        if(checked_names[(*pt_it).second] == true)
        {
          std::string exception_msg;
          exception_msg += " Two or more columns of the input matrix share the same species name (";
          exception_msg += (*pt_it).first;
          exception_msg += ") \n";  
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        } 
        else
          checked_names[(*pt_it).second] = true;

      } // else of if( pt_it == point_taxa.end() ) 

      column_to_node_vec.push_back((*pt_it).second);

    } // for( int i=0; i<names.size(); i++ )
     
    // Read the rest of the matrix, executing a query per line.

    int number_of_queries=0;

    while( in.good() )
    {
      line.clear();

      int count=0, start=0;
      std::vector<int> query_points;

      std::getline(in,line);

      // Exclude the first word of the line if first character is not zero or one.

      if( line[start] != '0' && line[start] != '1' )
        while( start < line.size() && line[start] != ',' )
          start++;

      for( int i=start; i<line.size(); i++)
      {
        a = line[i];

        if( a == '1' )
        {
          if(count >= names.size())
          {
            std::string exception_msg;
            exception_msg += " The matrix file has wrong syntax.\n";
            Exception_type excp;
            excp.get_error_message(exception_msg);
            Exception_functor excf;
            excf(excp);
          }

          query_points.push_back( column_to_node_vec[count] );

          count++;
        }
        else if ( a == '0' )
        {
          if(count >= names.size())
          {
            std::string exception_msg;
            exception_msg += " The matrix file has wrong syntax.\n";
            Exception_type excp;
            excp.get_error_message(exception_msg);
            Exception_functor excf;
            excf(excp);
          }

          count++;
        }
        else if( a != ',' && a !=' ' && a !='\r' && a !='\n')
        {
          std::string exception_msg;
          exception_msg += " The matrix file has wrong syntax.\n";
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

      } // for( int i=start; i<line.size(); i++)

      if( in.good() )
      {
        number_of_queries++;       

        if( query_points.size() < 1 )
          *ot++ = 0.0;
        else
        {
          Number_type single_sample_result = msr(query_points.begin(), query_points.end());

          if(standardised == false)  
            *ot++ = single_sample_result;
          else
          {
            if(mean_values[query_points.size()] == Number_type(-1.0))
            {
              mean_values[query_points.size()] = msr.compute_expectation(query_points.size());
              deviation_values[query_points.size()] = msr.compute_deviation(query_points.size());
            }  

            if(deviation_values[query_points.size()]==Number_type(0.0))
              *ot++ = single_sample_result - mean_values[query_points.size()];
            else
              *ot++ = (single_sample_result - mean_values[query_points.size()])/deviation_values[query_points.size()];
          }
          
        } // else of if( query_points.size() < 1 )

      } // if( in.good() )

    } // while( in.good() )

    in.close();
    return number_of_queries;

  } // csv_matrix_query( ... )

  
  template < class KernelType >  
  template < class Measure, class OutputIterator>
  int FunctionalMeasures::Measure_base<KernelType>::
  _matrix_query( std::vector<Point> &points, std::vector<std::string> &names, 
                 std::vector< std::vector<bool> > &matrix, 
                 Measure &msr, bool standardised, OutputIterator ot )
  {
    std::vector<Number_type> mean_values, deviation_values;
    std::vector<int> column_to_node_vec;

    mean_values.assign(points.size()+1, Number_type(-1.0));
    deviation_values.assign(points.size()+1, Number_type(-1.0));

    if(names.size() < points.size())
    {
      std::string warning(" Warning: the input matrix has fewer columns than the number of species in the point set.");
      Exception_functor().issue_warning(warning);
    }

    std::vector<bool> checked_names;
    checked_names.assign(points.size(),false);

    std::map<std::string,int> point_taxa;

    for(int i=0; i< points.size(); i++)
      point_taxa[points[i].taxon] = i;
 
    for( int i=0; i<names.size(); i++ )
    {
      typename std::map<std::string,int>::iterator pt_it = point_taxa.find(names[i]);

      if( pt_it == point_taxa.end() )
      {
        std::string exception_msg;
        exception_msg += " One of the species names in input the matrix was not found in the point set (";
        exception_msg += names[i];
        exception_msg += ") \n";  
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }
      else 
      {
        if(checked_names[(*pt_it).second] == true)
        {
          std::string exception_msg;
          exception_msg += " Two or more columns of the input matrix share the same species name (";
          exception_msg += (*pt_it).first;
          exception_msg += ") \n";  
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        } 
        else
          checked_names[(*pt_it).second] = true;
 
      } // if( pt_it == point_taxa.end() ) 

      column_to_node_vec.push_back((*pt_it).second);

    } // for( int i=0; i<names.size(); i++ )

    // Read the rest of the matrix, executing a query per row.

    for(int i=0; i<matrix.size(); i++)
    {
      std::vector<int> query_points; 

      for(int j=0; j<matrix[i].size(); j++)
        if(matrix[i][j] == true) 
          query_points.push_back( column_to_node_vec[j] );

      if( query_points.size() < 1 )
        *ot++ = 0.0;
      else
      {
        Number_type single_sample_result = msr(query_points.begin(), query_points.end());

        if(standardised == false)  
          *ot++ = single_sample_result;
        else
        {
          if(mean_values[query_points.size()] == Number_type(-1.0))
          {
            mean_values[query_points.size()] = msr.compute_expectation(query_points.size());
            deviation_values[query_points.size()] = msr.compute_deviation(query_points.size());
          }  

          if(deviation_values[query_points.size()]==Number_type(0.0))
            *ot++ = single_sample_result - mean_values[query_points.size()];
          else
            *ot++ = (single_sample_result - mean_values[query_points.size()])/deviation_values[query_points.size()];

        } // else of if(standardised == false)
          
      } // else of if( query_points.size() < 1 )

    } // for(int i=0; i<matrix.size(); i++)

    return matrix.size();

  } // _matrix_query( ... )

#endif //MEASURE_BASE_IMPL_H
