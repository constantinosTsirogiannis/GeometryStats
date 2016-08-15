#ifndef EXCEPTION_TYPE_H
#define EXCEPTION_TYPE_H

#include<string>
#include<cstdlib>

namespace ExceptionRelatedTypes
{
  class Exception_type
  {    
   public:

    void get_error_message(std::string msg)
    { _msg = msg;}

    std::string return_error_message()
    { return _msg; }

   private:

    std::string _msg;

  }; // class Exception_type


} // ExceptionRelatedTypes

#endif // EXCEPTION_TYPE_H
