/*   File: rodes.cc 
 
     The main program for running multiple processes sharing
     a data file <shared_file>. As long as there are undone
     grids left, we keep on working on them. When there are 
     none left, we quit.

     Compilation: type 'make' to see the alternatives. 

     Latest edit: Fri Jul 14 2000
*/


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <ctime>

#include "2d_classes.h"
#include "classes.h"
#include "convert.h"
#include "error_handler.h"
#include "list.h"
#include "return_map.h"
#include "request.h"

// The factor by which we widen the cone openings.
static const double ANG_FACTOR = 1.5;

// The minimal cone opening allowed.
static const double MIN_CONE_OPENING = 5.0 * Sup(DEG_TO_RAD);
 
// Wait one minute for a fresh grid.
static const unsigned WAIT_FOR_GRID = 60;

// Enumerations used for search for grid.
enum find {NONE_LEFT, WAITING_FOR_ONE, GOT_ONE};

enum command {START_TIMING, SHOW_TIMING, STOP_TIMING};

////////////////////////////////////////////////////////////////////

extern "C" 
{
  
  unsigned int sleep             (unsigned);
}

static void   print_info        (const char *);
static void   clock             (const command &);
static void   get_the_flags     (iterate &, const int &,
				 char *argv[], char *, char *);
static bool   get_a_grid        (iterate &, const char *, const char *);
static find   find_fresh_grid   (iterate &, const char *);
static void   terminate_process (const char *);                   
static void   work_on_grid      (iterate &, const char *, const char *);
static void   insert_it_List    (List<iterate> &, const char *,
				 const char *, const bool &); 
static void   update            (iterate &, const iterate &);

////////////////////////////////////////////////////////////////////

//using namespace std;

// Called by: none 
// Calls to : 'Take_care_of_the_flags', 'get_a_grid', 
//            'terminate_process', and 'work_on_grid'
int main(int argc, char *argv[])
{ 
  int counter = 0;
  char mult_file[99] = ""; // Name of the file owned by all processes.
  char proc_file[99] = ""; // Name of the file owned by this process.
  iterate it;

  get_the_flags(it, argc, argv, mult_file, proc_file);
 
  print_info(proc_file);  
  clock(START_TIMING);  

  while ( get_a_grid(it, mult_file, proc_file) )
    {
      counter++;

      work_on_grid(it, mult_file, proc_file);
      clock(SHOW_TIMING);
    }
  clock(STOP_TIMING); 
  terminate_process(proc_file);
 
  return 0;
}

////////////////////////////////////////////////////////////////////

// Called by: 'main'
// Calls to : none
static void print_info(const char *proc_file)
{
  cout << endl
       << "******************************************************" << endl
       << "****** RODES - Rigorous ODE Solver, version 1.1 ******" << endl
       << endl << "Latest compilation: " << __DATE__ << ", " << __TIME__ ;
#ifdef __linux
  cout << " (linux).";
#endif
#ifdef __sparc
  cout << " (sparc).";
#endif
  cout << endl;
#ifdef COMPUTE_C1
  cout << "Computing C0/C1-information ";
#else
  cout << "Computing C0-information ";
#endif
  cout << endl;
  cout.precision(3);
  cout.setf(ios::showpos);
  cout.setf(ios::scientific);
  cout << "ANG_FACTOR = " << ANG_FACTOR 
       << "; MIN_CONE_OPENING = " << MIN_CONE_OPENING * Inf(RAD_TO_DEG)
       << endl << endl;
  cout.unsetf(ios::showpos);
  cout.unsetf(ios::scientific);
  cout << "proc_file = " << proc_file << endl;
}

////////////////////////////////////////////////////////////////////

// Called by: 'main'
// Calls to : none
// 'START_TIMING' starts the clock, 'STOP_TIMING' stops the clock,
// and prints the elapsed time in hh:mm:ss.
static void clock(const command &cmd)
{
  struct tm *ptr;              // Used for timing the program 
  static time_t time_start;    // start time 

  if ( cmd == START_TIMING )
    {
      time_start = time(NULL);
      ptr = localtime(&time_start); 
      cout << "Started computations: " << asctime(ptr);  
    }
  else if ( cmd == SHOW_TIMING || cmd == STOP_TIMING )
    {
      int weeks   = 0;
      int days    = 0;
      int hours   = 0;
      int minutes = 0;
      double seconds;
      time_t time_stop;  

      time_stop = time(NULL);
      ptr = localtime(&time_stop);
      if ( cmd == STOP_TIMING )
	cout << "Completed computations: " << asctime(ptr);
      seconds = difftime(time_stop, time_start);
      cout << "Elapsed run time: "; 
      // Convert seconds to vv:dd:hh:mm:ss
      while ( seconds >= 60 )
	{
	  seconds -= 60;
	  minutes++;
	}
      while ( minutes >= 60 )
	{
	  minutes -= 60;
	  hours++;
	}
      while ( hours >= 24 )
	{
	  hours -= 24;
	  days++;
	}
      while ( days >= 7 )
	{
	  days -= 7;
	  weeks++;
	}
      if ( weeks > 0)
	cout << weeks << " week" << (weeks != 1 ? "s, " : ", ");
      if ( days > 0)
	cout << days << " day" << (days != 1 ? "s, " : ", ");
      if ( hours > 0)
	cout << hours << " hour" << (hours != 1 ? "s, " : ", ");
      if ( minutes > 0)
	cout << minutes << " minute" << (minutes != 1 ? "s, " : ", "); 
      // jjb -- changed to 'unsetf(ios...)' from setf( 0, ios...)
      cout.unsetf( ios::floatfield );

      cout << seconds << " second" << (seconds != 1 ? "s." : ".") << endl; 
    }
}

////////////////////////////////////////////////////////////////////

// Called by: 'main'
// Calls to : 'work_on_grid'
static void get_the_flags(iterate &it, const int &argc, char *argv[],
			  char *mult_name, char *proc_name)
{
  int add = 0;
#ifdef COMPUTE_C1
  add = 2; // Two more arguments
#endif

  if ( argc != 3 && argc != 6 + add && argc != 8 + add )
    { 
      cout << endl
	   << "Usage: " << endl;
      cout << "(1) \t" << argv[0] << " [proc_nr] <shared_file>\n\n";
      cout << "(2) \t" << argv[0] 
	   << " [proc_nr] <shared_file> [u] [v] [P]\n\n"
	   << "\twhich loads the rectangle R given by \n\n"
	   << "\t  R = 2^-P * [u - 1, u + 1]x[v - 1, v + 1]\n\n";
      cout << "(3) \t" << argv[0]
	   << " [proc_nr] <shared_file> [x.lo] [x.hi] [y.lo] [y.hi] [P]\n\n"
	   << "\twhich covers the rectangle with grids having diameter 2^-P.\n\n";
      cout << "When computing in C0/C1-mode:\n";
      cout << "In case (2) and (3) one must also give the angles (in degrees)\n"
	   << "for the cone boundary: [ang.lo] [ang.hi]" << endl << endl;
      exit(0);
    }

  // From here on, argc == 3, 6 or 8 (+ add)

  strcpy(mult_name, argv[2]);
  strcpy(proc_name, argv[2]);
  strcat(proc_name, "_");
  strcat(proc_name, argv[1]);

  List<iterate> it_List;

  if ( argc == 6 + add ) // Load a single iterate.
    {
      it.ndl.grd.u = atoi(argv[3]); // u coordinate.
      it.ndl.grd.v = atoi(argv[4]); // v coordinate.
      it.ndl.grd.P = atoi(argv[5]); // diam = 2^-P.
#ifdef COMPUTE_C1
      // CAPD intervalHull() takes iv1 and iv2 as args, so we convert
      // the floats to intervals. Don't really need this since we're
      // not computing C1.
      it.ndl.ang = DEG_TO_RAD * intervalHull( interval( atof( argv[6]) ), 
					      interval( atof( argv[7]) ) );
      it.ndl.pre_exp = LARGE_NUMBER;
      it.ndl.min_exp = LARGE_NUMBER;
#endif
      it.ndl.c_stat = NOT_DONE;
      it.ndl.h_stat = NOT_HIT;
      it.inf_grd = NULL_GRID; 
      it.sup_grd = NULL_GRID; 

      it_List += it;
      insert_it_List(it_List, mult_name, proc_name, true);
    } 
  else if ( argc == 8 + add ) // Load a rectangle -> several iterates.
    {
      BOX rect ( 2 ); 
      //Resize(rect, 2);
      double dbl[4];
      int power;

      dbl[0] = atof(argv[3]);  dbl[1] = atof(argv[4]);
      dbl[2] = atof(argv[5]);  dbl[3] = atof(argv[6]);
      power = atoi(argv[7]);

      rect(1) = intervalHull( interval( dbl[0] ), 
			      interval( dbl[1] ) );
      rect(2) = intervalHull( interval( dbl[2] ), 
			      interval( dbl[3] ) );
      cout << "rect:" << endl << rect << endl;
      rect_to_it_List(rect, power, it_List); // rectangle -> several iterates.
#ifdef COMPUTE_C1
      First(it_List);
      while ( !Finished(it_List) )
	{
	  Current(it_List).ndl.ang = DEG_TO_RAD * Hull(atof(argv[8]), atof(argv[9]));
	  Current(it_List).ndl.pre_exp = LARGE_NUMBER;
	  Current(it_List).ndl.min_exp = LARGE_NUMBER;
	  Next(it_List);
	}
#endif
      cout << "it_List:" << endl << it_List << endl;
      cout << "Is this OK? (y/n) ";
      char answer;
      cin >> answer;
      if ( answer == 'y' )
	insert_it_List(it_List, mult_name, proc_name, true);
      else
	{
	  cout << "Quitting on request!" << endl;
	  exit(0);
	}
    } 
}

////////////////////////////////////////////////////////////////////

// Called by: 'main'
// Calls to : 'get_file', 'find_fresh_grid', 'make_file_available', and 'rest'
static bool get_a_grid(iterate &it, const char *mult_name, const char *proc_name)
{
  find result;

  while (1)
    { 
      get_file(mult_name, proc_name);
      result = find_fresh_grid(it, proc_name);

      #ifdef DEBUG
      cout << "  get_a_grid : result = " << result << endl;
      #endif

      release_file(mult_name, proc_name);
      if ( result == NONE_LEFT )
	return false;
      else if ( result == GOT_ONE )
	return true;
      else if ( result == WAITING_FOR_ONE )
	sleep(WAIT_FOR_GRID);
    }
}

////////////////////////////////////////////////////////////////////

// Called by: 'get_a_grid'
// Calls to : none
static find find_fresh_grid(iterate &it, const char *proc_name)
{ 
    find result = NONE_LEFT;
 
    std::ifstream InFile(proc_name, ios::in);
    iterate temp_it;
    List<iterate> File_List;
    while ( InFile >> temp_it )
    {    
      if ( result != GOT_ONE )
	{
	  if ( temp_it.ndl.c_stat == NOT_DONE ) 
	    {
	      temp_it.ndl.c_stat = BEING_DONE;
	      it = temp_it;
	      result = GOT_ONE;
	    }
	  else if ( temp_it.ndl.c_stat == BEING_DONE )
	    result = WAITING_FOR_ONE;
	}
      File_List += temp_it;  
    }
  InFile.close();

    std::ofstream OutFile(proc_name, ios::out);
    First(File_List);
    while( !Finished(File_List) ) 
      {     
	OutFile << Current(File_List) << endl;
	Next(File_List);
      }
    OutFile.close();

  return result;
};

////////////////////////////////////////////////////////////////////

// Called by: 'main'
// Calls to : none - 'exit(0)'
static void terminate_process(const char *proc_name)
{
  cout << "terminate_process(" << proc_name << ")" << endl;
  exit(0);
}

////////////////////////////////////////////////////////////////////

// Called by: 'main' and 'Take_care_of_the_flags' 
// Calls to : 'Compute_the_return'(extern), 'insert_it_List',
//            'iterate_to_parcel', 'Get_Hull', 'pcl_List_to_it_List',
//            'Get_Image_Hull',  
static void work_on_grid(iterate &it, const char *mult_name, const char *proc_name)
{
    parcel pcl; //Resize(pcl.box, SYSDIM);
    //pcl.box

    List<iterate> it_List;
    List<parcel> pcl_List;

  iterate_to_parcel(it, pcl);

  try
    {
  /***************************************************************/
  /*        HERE WE MAKE THE ONLY CALL TO THE INTEGRATOR         */
  /*                                                             */
      if ( it.ndl.c_stat != RESERVED )
	Compute_the_return(pcl, pcl_List);

      cout << "HERE in rodes.cc" << endl;
  /*                                                             */
  /*                                                             */
  /***************************************************************/
    }
  catch( Error_Handler error )
    {
      error.Print_Message();
      cout << "Thrown and caught for the iterate" << endl
	   << it << endl;
      it.ndl.c_stat = FAILED; 
    }

  if ( it.ndl.c_stat == FAILED || it.ndl.c_stat == RESERVED )
    it_List += it; // it_List contains it only.
  else
    { // Store the returns in the list
      parcel pcl_hull; 
      Get_Hull(pcl_hull, pcl_List);
      it.ndl.c_stat = DONE; 
#ifdef COMPUTE_C1
      it.ndl.min_exp = Inf(pcl_hull.expansion); // Inf(E_i), E_i = Hull(E_{i,j}).
#endif
  /***************************************************************/
  /*        HERE WE GENERATE THE PRE-EXPANSION ESTIMATES         */
  /*                                                             */
      pcl_List_to_it_List(pcl_List, it.ndl.grd.P, it_List);

  /*  AFTER: it.ndl.h_stat == HIT, it.ndl.c_stat == NOT_DONE,    */
  /*         it.ndl.min_exp == LARGE_NUMBER.                     */
  /***************************************************************/

      New_Get_Image_Hull(it, it_List); // Gets 'inf_grd' and 'sup_grd'.
      it_List += it;                   // Add the initial it to the end of the List.
    }

  insert_it_List(it_List, mult_name, proc_name, false);
}

////////////////////////////////////////////////////////////////////

// Called by: 'work_on_grid' 
// Calls to : 'update'
static void insert_it_List(List<iterate> &Add_List, const char *mult_name,
			   const char *proc_name, const bool &external_input)
{
    iterate temp_it;
    List<iterate> Old_List;

    get_file(mult_name, proc_name);
    std::ifstream InFile(proc_name, ios::in);

    while ( InFile >> temp_it )
      Old_List += temp_it;  
    InFile.close();

    iterate old_it;
    iterate add_it;
    List<iterate> New_List;
    
    if ( Length(Old_List) > 0 )
      { // Length(Old_List) > 0
	First(Old_List);
	while( !Finished(Old_List) ) 
	  { // !Finished(Old_List)                       
	    old_it = Current(Old_List);	  
	    if ( !IsEmpty(Add_List) )
	      {
		First(Add_List);
		while( !Finished(Add_List) ) 
		  {     
		    add_it = Current(Add_List);
		    if ( old_it.ndl.grd == add_it.ndl.grd )
		      {
			/********************************************************/
			/*        HERE WE UPDATE THE STATUS OF 'old.it'         */
			/*                                                      */
			update(old_it, add_it);
			/*                                                      */
			/*                                                      */ 
			/********************************************************/
			RemoveCurrent(Add_List);
			if ( IsEmpty(Add_List) ) 
			  break;
		      }
		    else                        
		      Next(Add_List);
		  }
	      }
	    New_List += old_it;  
	    Next(Old_List);
	  }
      }
    
    if ( Length(Add_List) > 0 )
      { // Append the remains of Add_List to New_List.
	First(Add_List);
	while( !Finished(Add_List) ) 
	  { 
	    add_it = Current(Add_List);
	    if ( !external_input ) // Only widen cone openings if they
	      {                    // come from internal computations.
		/********************************************************/
		/*    Widen the cone opening to minimize the risk of    */
		/*    having to recompute due to cone leakage.          */
		/*                                                      */
#ifdef COMPUTE_C1                             
		add_it.ndl.ang = Rescale(add_it.ndl.ang, ANG_FACTOR); 
		if ( diam(add_it.ndl.ang) < MIN_CONE_OPENING )
		  {
		    double mid = Mid(add_it.ndl.ang);
		    add_it.ndl.ang = mid + SymHull(MIN_CONE_OPENING / 2.0);
		  }
#endif                  
		/*   Now the cone opening is at least MIN_CONE_OPENING. */
		/*                                                      */ 
		/********************************************************/    
	      }
	    New_List += add_it;                     
	    Next(Add_List);                     
	  }                                     
      }
    std::ofstream OutFile(proc_name, ios::out);
    First(New_List);
    while( !Finished(New_List) ) 
      { // Write New_List to owned file.    
	OutFile << Current(New_List) << endl;;
	Next(New_List);
      }
    OutFile.close();
    release_file(mult_name, proc_name);
} 

////////////////////////////////////////////////////////////////////

// Called by: 'insert_it_List' 
// Calls to : none
static void update(iterate &old_it, const iterate &new_it)
{
  iterate result = old_it;

  if ( new_it.ndl.c_stat == NOT_DONE ) // 'new_it' is an image...
    { // 'new_it' is not the source it.
      result.ndl.h_stat = new_it.ndl.h_stat;// ...and so it is HIT.
#ifdef COMPUTE_C1
      // 'old_it' has already been updated, or is now
      // under computation by a different process.
      result.ndl.pre_exp = Min(old_it.ndl.pre_exp, new_it.ndl.pre_exp);
      if ( Subset(new_it.ndl.ang, old_it.ndl.ang) )
	{ /* No need to widen the cone */ }
      else
	{ // We have to recompute with a wider cone.
	  interval dummy;
	  if ( Intersection(dummy, new_it.ndl.ang, old_it.ndl.ang) )
	    { // Double the angular differences.
	      double upper_ang_diff = Sup(new_it.ndl.ang) - Sup(old_it.ndl.ang);
	      if ( upper_ang_diff < 0.0 )
		upper_ang_diff = 0.0;
	      double lower_ang_diff = Inf(new_it.ndl.ang) - Inf(old_it.ndl.ang);
	      if ( lower_ang_diff > 0.0 )
		lower_ang_diff = 0.0;
	      result.ndl.ang = old_it.ndl.ang + 
		ANG_FACTOR * Hull(lower_ang_diff, upper_ang_diff);
	    }
	  else // The cones don't even intersect.
	    result.ndl.ang = Rescale(Hull(new_it.ndl.ang, old_it.ndl.ang), ANG_FACTOR);
	  if ( old_it.ndl.c_stat == BEING_DONE )
	    result.ndl.c_stat = DO_AGAIN; 
	  if ( old_it.ndl.c_stat == DONE )
	    result.ndl.c_stat = NOT_DONE; 
	  // We do not alter c_stat == FAILED or RESERVED.
	}
#endif // COMPUTE_C1
    }
  else // 'new_it' is the source, i.e., 'old_it' is
    {  // now under computation by this process.
      if ( old_it.ndl.c_stat == DO_AGAIN ) // The cone has been widened
	{                                  // during this computation.
	  result.ndl.c_stat = NOT_DONE;    // Only this process may release 
	}                                  // the 'it' for re-computation. 
      else // An uninterrupted computation.
	{
	  result.ndl.c_stat = new_it.ndl.c_stat; // DONE, FAILED or RESERVED.
	  result.inf_grd = new_it.inf_grd;       // Save the hull of
	  result.sup_grd = new_it.sup_grd;       // the images.
#ifdef COMPUTE_C1
	  result.ndl.min_exp = new_it.ndl.min_exp; // Save the minimal expansion.
#endif // COMPUTE_C1
	}
    }
  old_it = result; // Pass the result by reference.
}

////////////////////////////////////////////////////////////////////
