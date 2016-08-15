#ifndef EXCEPTION_FUNCTOR_LOUD_H
#define EXCEPTION_FUNCTOR_LOUD_H

#include<string>
#include<cstdlib>
#include<exception>
#include"Exception_type.h"

namespace ExceptionRelatedTypes
{
  class Exception_functor
  {
   public:

    Exception_functor(){}
    
    void operator()(Exception_type excp)
    {
      std::string msg = excp.return_error_message();
      std::cout << " Error:" << msg << std::endl;
      exit(1);
    }

    void issue_warning(std::string str)
    { std::cout << str << std::endl << std::endl;}

  }; // class Exception_functor

} // ExceptionRelatedTypes

#endif // EXCEPTION_FUNCTOR_LOUD_H
