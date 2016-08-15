#ifndef EXCEPTION_FUNCTOR_SILENT_H
#define EXCEPTION_FUNCTOR_SILENT_H

#include<string>
#include<cstdlib>
#include<exception>
#include<vector>
#include"Exception_type.h"


class Warning_list_type
{
 public: 

  int number_of_warnings()
  { return _warnings.size();}

  std::string operator[](int i)
  { return _warnings[i];}

  void push_warning(std::string &w)
  { _warnings.push_back(w);} 

  void clear()
  { _warnings.clear();}

 private:

  std::vector<std::string> _warnings;
}; 

Warning_list_type warning_list;

namespace ExceptionRelatedTypes
{
  class Exception_functor
  {
   public:

    Exception_functor(){}
    
    void operator()(Exception_type excp)
    { throw excp;}

    void issue_warning(std::string str)
    { warning_list.push_warning(str);}

  }; // class Exception_functor

} // ExceptionRelatedTypes

#endif // EXCEPTION_FUNCTOR_SILENT_H
