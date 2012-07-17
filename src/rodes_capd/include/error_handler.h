/*   File: error_handler.h 
 
     Contains exception classes that may 
     be called by the mathematical functions.

     Latest edit: Mon Feb 7 2000
*/

#ifndef ERROR_HANDLER_H
#define ERROR_HANDLER_H

#include <iostream.h>

class Error_Handler
{
 private: 
  char *message;
 public:
  Error_Handler() : message("ERROR: unknown origin") { }
  Error_Handler(char *msg) { message = msg; }
  void Print_Message() const { cout << message << endl; }
};

#endif // ERROR_HANDLER_H
