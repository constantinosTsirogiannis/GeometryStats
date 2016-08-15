#ifndef CONSTRUCT_POINT_SET_IMPL_H
#define CONSTRUCT_POINT_SET_IMPL_H

  // Input: A csv file storing a matrix where each row corresponds to a 
  // species and each column corresponds to a coordinate of the point representing 
  // this species (except the first column which stores the species names).
  // The second argument is an output iterator of the type std::back_insert_iterator<std::vector< Point > >.
  // The third argument is an integer d that indicates the number of coordinates that
  // we should extract for any point in the matrix e.g. for d=3 the three first coordinates
  // are extracted from every row. If d=-1 then all coordinates are extracted for each species.
  // Output: A vector of d-dimensional points, each point also storing the name of the corresponding species.
  
  template < class KernelType >  
  template<class OutputIterator>
  void FunctionalMeasures::Construct_point_set<KernelType>::
  operator()( char *filename, OutputIterator ot, int d)
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

    // Through out the first row

    std::string line;
    std::getline(in,line);

   // Read line-by-line the rest of the rows and 

    std::vector<std::string> names;
    std::vector<std::vector<double> > vals;

    while( in.good() )
    {
      line.clear();

      int start=0, count=0;
      std::getline(in,line);
     

      // First word in the row is a taxon name.

      names.push_back(std::string(""));

      while( start < line.size() && line[start] != ',' )
      {
        names.back().push_back(line[start]);
        start++;
      }

      if(start >= int(line.size())-1)
      {
        names.pop_back();
        break;
      }  

      vals.push_back(std::vector<double>());

      std::stringstream ss(line.substr(start+1));

      while(ss.good() && count < d || d==-1)
      {
        double value;
        char a;

        ss >> value;
        ss >> a;
   
        vals.back().push_back(value);

        if( a != ',' && a !='\r' && a !='\n' )
        {
          std::string exception_msg;
          exception_msg += " The matrix file has wrong syntax.\n";
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

        count++;

      } // while(ss.good() && count < d || d==-1)

      if(vals.size() > 1)
        if( vals.back().size() != vals[0].size() )
        {
          std::string exception_msg;
          exception_msg += " The matrix rows do not have the same number of elements.\n";
          Exception_type excp;
          excp.get_error_message(exception_msg);
          Exception_functor excf;
          excf(excp);
        }

    } // while( in.good() )

    in.close();

    this->operator()(names,vals,ot,-1);

  } // operator()( char *filename, OutputIterator ot, int d=-1)


  // Input: A vector with species names, and a matrix where each row corresponds to a 
  // species and each column corresponds to a coordinate of the point representing 
  // this species.
  // The third argument is an output iterator of the type std::back_insert_iterator<std::vector< Point > >.
  // The fourth argument is an integer d that indicates the number of coordinates that
  // we should extract for any point in the matrix e.g. for d=3 the three first coordinates
  // are extracted from every row. If d equals -1 then all coordinates are extracted for each species.
  // Output: A vector of d-dimensional points, each point also storing the name of the corresponding species.
  
  template < class KernelType >  
  template<class OutputIterator>
  void FunctionalMeasures::Construct_point_set<KernelType>::
  operator()( std::vector<std::string> &names, std::vector<std::vector<double> > &vals,  
               OutputIterator ot, int d)
  {
    if(names.size() != vals.size())
    {
      std::string exception_msg;
      exception_msg += " The matrix rows and the number of species names are not equal.\n";
      Exception_type excp;
      excp.get_error_message(exception_msg);
      Exception_functor excf;
      excf(excp);
    }

    if(vals.size() == 0)
      return;

    if(d==-1)
    {
      for(int i=0; i<vals.size(); i++)
      {
        Point p;

        p.taxon = names[i];

        for(int j=0; j<vals[i].size(); j++)
          p.push_back(vals[i][j]);

        *ot++ = p;

      } // for(int i=0; i<vals.size(); i++)      

    } 
    else // of if(d==-1)
    {
      if(d > vals[0].size())
      {
        std::string exception_msg;
        exception_msg += " The specified number of dimensions exceeds the number of columns in the matrix.\n";
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      } 

      for(int i=0; i<vals.size(); i++)
      {
        Point p;

        p.taxon = names[i];
 
        for(int j=0; j<d; j++)
          p.push_back(vals[i][j]);

        *ot++ = p;

      } // for(int i=0; i<vals.size(); i++)

    } // if(d==-1)

  } // operator()( std::vector<std::string> &names, std::vector<std::vector<double> > &vals, ... )


  // Similar as the previous function, except now matrix vals is represented by
  // an object of the Rcpp type NumericMatrix.
  
//  template < class KernelType >  
//  template<class OutputIterator>
//  void FunctionalMeasures::Construct_point_set<KernelType>::
//  operator()( std::vector<std::string> &names, NumericMatrix &vals,  
//               OutputIterator ot, int d=-1)
//  {
//    if(names.size() != vals.size())
//    {
//      std::string exception_msg;
//      exception_msg += " The matrix rows and the number of species names are not equal.\n";
//      Exception_type excp;
//      excp.get_error_message(exception_msg);
//      Exception_functor excf;
//      excf(excp);
//    }
//
//    if(vals.nrow() == 0 || vals.ncol() == 0)
//      return;
//
//    if(d==-1)
//    {
//      for(int i=0; i<vals.nrow(); i++)
//      {
//        Point p;
//
//        p.taxon = names[i];
//
//        for(int j=0; j<vals.ncol(); j++)
//          p.push_back(vals(i,j));
//
//        *ot++ = p;
//
//      } // for(int i=0; i<vals.nrow(); i++)      
//
//    } 
//    else // of if(d==-1)
//    {
//      if(d > vals.ncol())
//      {
//        std::string exception_msg;
//        exception_msg += " The specified number of dimensions exceeds the number of columns in the matrix.\n";
//        Exception_type excp;
//        excp.get_error_message(exception_msg);
//        Exception_functor excf;
//        excf(excp);
//      } 
//
//      for(int i=0; i<vals.nrow(); i++)
//      {
//        Point p;
//
//        p.taxon = names[i];
// 
//        for(int j=0; j<d; j++)
//          p.push_back(vals(i,j));
//
//        *ot++ = p;
//
//      } // for(int i=0; i<vals.nrow(); i++)
//
//    } // if(d==-1)
//
//  } // operator()( std::vector<std::string> &names, NumericMatrix &vals, ... )

  template < class KernelType >  
  void FunctionalMeasures::Construct_point_set<KernelType>::
        write_point_set_to_file(std::vector<Point> &points, char *filename)
  { 
    if(points.size() == 0)
    {
      std::ofstream of(filename);
      of.close();
      return;
    }  

    for(int i=0; i<points.size(); i++)
      if(points[i].dim() != points[0].dim())
      {
        std::string exception_msg;
        exception_msg += " Not all points in the input set have  the same dimension.\n";
        Exception_type excp;
        excp.get_error_message(exception_msg);
        Exception_functor excf;
        excf(excp);
      }

    std::ofstream of(filename);
 
    of << "TAXON"; 

    for(int i=0; i<points[0].dim(); i++)
      of << " , C" << i; 

    of << std::endl;

    for(int i=0; i<points.size(); i++)
    {
      of << points[i].taxon; 

      for(int j=0; j<points[i].dim(); j++)
        of << " , " << points[i][j]; 

      of << std::endl;
    }

    of.close(); 

  } // write_point_set_to_file(std::vector<Point> &points, char *filename)

#endif //CONSTRUCT_POINT_SET_IMPL_H
