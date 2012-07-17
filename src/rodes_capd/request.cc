/*   File: request.cc
 
     A set of functions used for requesting
     access to a shared file ('source_name')
     having multiple users.

     Latest edit: Sat Mar 18 2000
*/

#include "request.h"

using namespace std;

////////////////////////////////////////////////////////////////////

static const unsigned WAIT_FOR_FILE = 1; // Wait 1 sec if the file is occupied.

extern "C" 
{
  unsigned int sleep ( unsigned );
  //  static void sleep  (unsigned);
  int  rename (const char *, const char *);
}

////////////////////////////////////////////////////////////////////

// Calls to : 'is_file_available', and 'sleep' (external C)
void get_file(const char *source_name, const char *image_name)
{ 
  while( 0 != rename(source_name, image_name) )
    sleep(WAIT_FOR_FILE);
}

////////////////////////////////////////////////////////////////////

// Calls to : 'exit(1)', 'rename' (external C)
void release_file(const char *source_name, const char *image_name)
{
  int system_result;

  system_result = rename(image_name, source_name);

  if ( system_result != 0 ) // Ops! This should never occur!
    {
      cout << "Error: " << "release_file(" << source_name 
	   << ", " << image_name << ") failed!" << endl;
      exit(1);
    }
}

////////////////////////////////////////////////////////////////////
