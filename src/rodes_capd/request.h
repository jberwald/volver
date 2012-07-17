/*   File: request.h 
 
     A set of functions used for requesting
     access to a shared file ('source_name')
     having multiple users.

     Latest edit: Sat Mar 18 2000
*/

#ifndef REQUEST_H
#define REQUEST_H

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdlib>
#include <string>

////////////////////////////////////////////////////////////////////

void get_file     (const char *, const char *);

void release_file (const char *, const char *);

////////////////////////////////////////////////////////////////////

#endif // REQUEST_H
