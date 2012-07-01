/*   File: expansion.cc 
 
     Reads the data set from <infile>, and computes the (coarse) accumulated expansion
     estimates over several iterates. The fundamental domain is defined by the boundary
     u values, u_min and u_max. The main structure is as follows:

       (1) Input bounds for the fundamental domain F_0;
       (2) Compute <R>(F_0);
       (3) Define F_1 = <R>(F_0) \cap (N \ F_0);
       (4) Use N_i \in F_1 as starting elements;
       (5) Iterate each N_i until it hits F_0;
       (6) Signal if the accumulated expansion is less than 2.0;

     Usage: expansion <infile>
 
     Tips: The values <u,v> = <-128, 512> seem to do the trick.

     Compilation: make expansion

     Latest edit: Tue Jul 18 2000
*/

#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <string.h>

#include "2d_classes.h"
#include "classes.h"
#include "list.h"
#include "request.h"

////////////////////////////////////////////////////////////////////

static void Get_The_Flags        (const int &, char *argv[], char *, 
				  int &, int &);

static void Retrieve_Indata      (const char    *, List<iterate> &);

static void Coarsen              (List<iterate> &);

iterate     Hull                 (const iterate &, const iterate &);

static void Generate_F_0_List    (List<iterate> &, List<iterate> &, int, int);

static void Generate_F_1_List    (List<iterate> &, List<iterate> &, List<iterate> &,
				  List<iterate> &);
static void Find_F_0_Inv_Sets    (List<iterate> &);

static void Compute_Image        (List<iterate> &, List<iterate> &, List<iterate> &,
				  List<iterate> &);

static void Reflect_To_Upper     (List<iterate> &);    

static void Symmetrize           (List<iterate> &);

static void Remove_Intersections (List<iterate> &, List<iterate> &);

static void Get_Return_Grids     (List<iterate> &, List<iterate> &, List<iterate> &,
				  const iterate &, bool);

static void Remove_Redundancies  (List<iterate> &);

static void Append               (List<iterate> &, List<iterate> &);

static void Copy                 (List<iterate> &, List<iterate> &);

double      Flow_Along           (const iterate &, List<iterate> &, List<iterate> &,
				  List<iterate> &);

bool        Equal                (const    grid &, const    grid &);

bool        Intersect            (const iterate &, List<iterate> &);

static void Empty                (List<iterate> &);

static void Get_Exp_And_Image_Hull(grid &, grid &, double &, List<iterate> &);
 
const double SCALE = pow(2, - 8);

////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  char source_file[99] = "";
  bool use_symmetry    = true; // Use_Symmetry by default.
  int U_MIN, U_MAX;
  iterate f_1_it;
  List<iterate> Total_List, F_0_List, F_1_List, E_List, Empty_List;

  Get_The_Flags(argc, argv, source_file, U_MIN, U_MAX);

  Retrieve_Indata(source_file, Total_List);
  Reflect_To_Upper(Total_List);    

  // Joins together elements of the same branch that are
  // stacked vertically. Total_List becomes ordered in u.
  // To compensate for info loss, we widen the span of
  // of inf_grd.v and sup.grd.v.
  Coarsen(Total_List);

  //  cout << "Total_List: " << endl << Total_List << endl; exit(0);

  if ( use_symmetry )
    {
      Symmetrize(Total_List);
      cout << "Symmetrized Total_List has cardinality "
	   << Length(Total_List) << endl;
    }

  // Define F_0 (F_0_List) via the bounds U_MIN and U_MAX.
  // Also computes the smallest min_exp in F_0, and verifies
  // that it indeed contains a fundamental domain.          
  // F_0 is always taken to be on the upper branch.
  Generate_F_0_List(F_0_List, Total_List, U_MIN, U_MAX);

  // Finds all N_i in F_0 such that <R>(N_i) intersects
  // N_i^+ or  N_i^-. Checks that E_i > \sqrt{2}.
  Find_F_0_Inv_Sets(F_0_List);

  // Compute <R>(F_0) and refine the image of F_0 by
  // defining F_1 = <R>(F_0) \cap (N \ F_0).
  // F_1 is always taken to be on the upper branch.
  Generate_F_1_List(F_0_List, F_1_List, Total_List, Empty_List);

  Symmetrize(F_0_List);

  // Take F_1 (F_1_List) as starting elements and iterate each
  // element until entering F_0. If the accumulated expansion 
  // is less than two, we signal an error.
  int counter = 1;
  int final = Length(F_1_List);
  double local_acc_exp;
  double min_acc_exp = Machine::PosInfinity;
  cout << endl << "Starting the proper iterations..." << endl;
  while( !IsEmpty(F_1_List) ) // Empties F_1_List.
    {
      f_1_it = First(F_1_List);
      RemoveCurrent(F_1_List);   
      cout << "Grid #" << counter << "/" << final << "\t " << f_1_it.ndl.grd;
      cout.flush();
      local_acc_exp = Flow_Along(f_1_it, F_0_List, Total_List, E_List);  
      cout << "; -- " << local_acc_exp << endl;
      min_acc_exp = Min(min_acc_exp, local_acc_exp);
      counter++; 
    }
  cout << endl << "The smallest accumulated expansion was " << min_acc_exp << endl;
  cout << "Bye!" << endl;

  return 0;
}

////////////////////////////////////////////////////////////////////

static void Get_The_Flags(const int  &argc, char *argv[], char *source_name,
			  int &U_MIN, int &U_MAX)
{
  if ( argc != 2)
    {
      cout << " Usage: " << endl;
      cout << "(1) \t" << argv[0] << " <input_file>\n";
      exit(0);
    }

  // From here on, argc == 2.
  strcpy(source_name, argv[1]);
  cout << endl << endl
       << "********************** EXPANSION 1.0 **********************"
       << endl << endl;
  cout << "Enter the u values for the fundamental domain F_0." << endl;
  cout << "This domain is always taken in the upper branch." << endl;
  char answer = 'n';
  
  while ( answer == 'n' )
    {
      cout << "Enter u_min: "; cin >> U_MIN;
      cout << "Enter u_max: "; cin >> U_MAX;
      cout << "This corresponds to F_0|x_1 = [" << U_MIN * SCALE
	   << ", " << U_MAX * SCALE << "]. Is this OK? (y/n) ";
      cin >> answer;
    }
}

////////////////////////////////////////////////////////////////////

// Loads all data from source_file into Total_List. If use_symmetry is
// true, all elements are given a symmetric twin.
static void Retrieve_Indata(const char *source_file, List<iterate> &Total_List)
{
  char tmp_file[99] = "";

  strcpy(tmp_file, source_file);
  strcat(tmp_file, ".exp");

  // Open up.
  get_file(source_file, tmp_file);
  ifstream Input_File(tmp_file, ios::in);
  if ( !Input_File )
    {
      cout << "File " << tmp_file << " could not be opened!" << endl;
      exit(1);
    }
 
  // Get all input data.
  iterate it;
  while ( Input_File >> it )
    {
      if ( it.inf_grd == NULL_GRID )
	{ cout << "Not loading " << it << endl;} // Don't load this iterate. 
      else if ( it.ndl.h_stat == NOT_HIT )
	{} // Don't load this iterate. 
      else
	Total_List += it;
    }

  // Close down.
  Input_File.close();
  release_file(source_file, tmp_file);

  cout << endl << "Loaded " << Length(Total_List) << " elements all in all." << endl;
  if ( Length(Total_List) == 0 ) // Sanity check.
    exit(0);
}

////////////////////////////////////////////////////////////////////

// Joins together elements of the same branch that are stacked vertically. 
static void Coarsen(List<iterate> &Total_List)
{
  static const int BIG_INT = 50; // At most 50 pieces stacked on top of each other.
  int min_u, max_u;
  iterate it;
  List<iterate> Temp1_List, Temp2_List;

  cout << "Coarsening Total_List, please wait..." << endl;

  it = First(Total_List);
  min_u = it.ndl.grd.u;
  max_u = min_u;
  Next(Total_List);
  while( !Finished(Total_List) )
    {
      it = Current(Total_List);
      if ( min_u > it.ndl.grd.u )
	min_u = it.ndl.grd.u;
      if ( max_u < it.ndl.grd.u )
	max_u = it.ndl.grd.u;
      Next(Total_List);
    } // Now min_u and max_u enclose the u values of Total_List.

  for ( int u = max_u; u >= min_u; u-- )
    { // Loop through all possible u values.
      First(Total_List);
      while( !IsEmpty(Total_List) && !Finished(Total_List) )
	{
	  it = Current(Total_List);
	  if ( it.ndl.grd.u == u )
	    {
	      Temp1_List += it;
	      RemoveCurrent(Total_List);
	    }
	  else
	    Next(Total_List);
	} // Now Temp1_List contains all (if any) it with it.ndl.grd.u == u.
      if ( !IsEmpty(Temp1_List) )
	{
	  if ( Length(Temp1_List) > BIG_INT )
	    cout << "Warning (Coarsen): need to increase BIG_INT = " << BIG_INT << endl;
	  it = First(Temp1_List);
	  RemoveCurrent(Temp1_List);
	  while( !IsEmpty(Temp1_List) ) // Clean Temp1_List out for the next run.
	    {
	      it = Hull(it, Current(Temp1_List));
	      RemoveCurrent(Temp1_List);
	    } // Now it contains the hull of all it:s with it.ndl.grd.u == u.
	  it.inf_grd.v -= BIG_INT;          
	  it.sup_grd.v += BIG_INT;   
	  Temp2_List += it;
	}
    } // Now Temp2_List contains all coarsed it:s from right to left.
  if ( !IsEmpty(Total_List) )
    {
      cout << "Hey! Total_List isn't empty!" << endl;
      Empty(Total_List); // Clean out Total_List.
    }
  Append(Total_List, Temp2_List);
  cout << "Coarsened Total_List has cardinality " << Length(Total_List) << endl;
}

////////////////////////////////////////////////////////////////////

// Get the hull of two iterates. Now, we don't care about the location
// (<u,v,P>) or the cone angles.
iterate Hull(const iterate &it_1, const iterate &it_2)
{
  iterate result = it_1;

  result.ndl.min_exp = Min(it_1.ndl.min_exp, it_2.ndl.min_exp);
  result.ndl.pre_exp = Min(it_1.ndl.pre_exp, it_2.ndl.pre_exp);

  if ( it_1.inf_grd.u > it_2.inf_grd.u )
    result.inf_grd.u = it_2.inf_grd.u;
  if ( it_1.inf_grd.v > it_2.inf_grd.v )
    result.inf_grd.v = it_2.inf_grd.v;
  if ( it_1.sup_grd.u < it_2.sup_grd.u )
    result.sup_grd.u = it_2.sup_grd.u;
  if ( it_1.sup_grd.v < it_2.sup_grd.v )
    result.sup_grd.v = it_2.sup_grd.v;

  return result;
}

////////////////////////////////////////////////////////////////////

// Generates F_0_List = Total_List \cap [U_MIN, U_MAX] (upper branch),
// and computes the smallest min_exp in F_0.
static void Generate_F_0_List(List<iterate> &F_0_List, List<iterate> &Total_List,
			      int U_MIN, int U_MAX)
{
  char answer;
  double min_min_exp = Machine::PosInfinity;
  iterate it;

  First(Total_List);
  while( !Finished(Total_List) )
    {
      it = Current(Total_List);
      if ( 5 * it.ndl.grd.v > 2 * it.ndl.grd.u ) // Upper brach.
	if ( U_MIN <= it.ndl.grd.u && it.ndl.grd.u <= U_MAX )
	  {
	    F_0_List += it;
	    min_min_exp = Min(min_min_exp, it.ndl.min_exp);
	  }
      Next(Total_List);
    }
  cout << endl << "Loaded " << Length(F_0_List) << " elements into F_0." << endl;
  int min_u = Last(F_0_List).ndl.grd.u;
  int max_u = First(F_0_List).ndl.grd.u;
  cout << "This corresponds to F_0|x_1 = [" << (min_u - 1) * SCALE
       << ", " << (max_u + 1)* SCALE << "]." << endl;
  cout << "Smallest min_exp in F_0 = " << min_min_exp << endl;

  if ( min_min_exp < 1.0 )
    {
      cout << endl << "Wanna quit? (y/n) "; cin >> answer;
      if ( answer == 'y' )
	exit(0);
    }

  // Now check for a fundamental domain.
  List<iterate> Return_List, E_List; 
  bool first_time = true;
  // Maps to upper branch.
  Get_Return_Grids(Return_List, Total_List, E_List, First(F_0_List), true); 
  while ( !IsEmpty(Return_List) )
    {
      if ( !Intersect(First(Return_List), F_0_List) && first_time)
	{
	  cout << "F_0 does not contain a fundamental domain!" << endl;
	  cout << endl << "Wanna quit? (y/n) "; cin >> answer;
	  if ( answer == 'y' )
	    exit(0);
	  first_time = false;
	}
      RemoveCurrent(Return_List);
    }
  cout << "F_0 contains a fundamental domain." << endl;
}

////////////////////////////////////////////////////////////////////

// Finds all N_i in F_0 such that <R>(N_i) intersects N_i^+ or  N_i^-.
// Then checks that E_i > \sqrt{2}.
static void Find_F_0_Inv_Sets(List<iterate> &F_0_List)
{
  bool found_one    = false;
  double min_factor = Machine::PosInfinity;
  double factor;
  iterate it;
  List<iterate> Return_List, Space_List, E_List;

  cout << endl << "Finding F_0-invariant sets..." << endl;

  Copy(F_0_List, Space_List); // Copies F_0 into Space.
  Symmetrize(Space_List);

  First(F_0_List);
  while( !Finished(F_0_List) )
    {
      it = Current(F_0_List);

      Empty(Return_List);
      Get_Return_Grids(Return_List, Space_List, E_List, it, true); // Maps to upper branch. 

      if ( Intersect(it, Return_List) )
	{
	  found_one = true;
	  factor = Max(it.ndl.min_exp, it.ndl.pre_exp);
	  min_factor = Min(min_factor, factor);
	  cout << "  Iterate " << it.ndl.grd << " is mapped accross its twin image" << endl;
	  if ( factor < sqrt(2) )
	    cout << "  but only has" << factor << " in expansion." << endl; 
	  else
	    cout << "  and has " << factor << " in expansion." << endl; 
	}
      Next(F_0_List);
    }
  if ( found_one )
    cout << "Any orbit completely within F_0 satisfies " << endl
	 << "|DR^n(x)*v| > a^n*|v|, with a = " << min_factor << endl; 
  else
    cout << "No orbit completely within F_0 could be found." << endl;
}

////////////////////////////////////////////////////////////////////

// Loops through F_0_List and computes the image of each element. These are
// stored in F_1_List, which is also trimmed disjoint from F_0_List.
// F_1_List is always taken on the upper branch.
static void Generate_F_1_List(List<iterate> &F_0_List, List<iterate> &F_1_List,
			      List<iterate> &Total_List, List<iterate> &E_List)
{
  cout << endl << "Generating F_1_List, please wait..." << endl;

  Compute_Image(F_0_List, F_1_List, Total_List, E_List); // Generates F_1_List.

  Reflect_To_Upper(F_1_List);                            // F_1_List -> upper branch..

  Remove_Redundancies(F_1_List);                         // Removes Equal elements.
  
  Remove_Intersections(F_1_List, F_0_List);              // Trims F_1_List.

  cout << "Loaded " << Length(F_1_List) << " elements into F_1." << endl;

  int min_u = Last(F_1_List).ndl.grd.u;
  int max_u = First(F_1_List).ndl.grd.u;
  cout << "This corresponds to F_1|x_1 = [" << (min_u - 1)* SCALE
       << ", " << (max_u + 1)* SCALE << "]." << endl;
}

////////////////////////////////////////////////////////////////////

// Computes the iterates of f_1_it, and keeps track of the accumulated 
// expansion over its orbit. When an iterate enters the fundamental domain
// a warning is printed if the accumulated expansion is less than 2.0. 
// The iterates are simply the union of all images.
double Flow_Along(const iterate &f_1_it, List<iterate> &F_0_List, 
		  List<iterate> &Total_List, List<iterate> &E_List)
{
  bool show_info = false;
  bool failure   = false;
  double acc_exp = f_1_it.ndl.pre_exp; // We have an (invisible) iteration.
  double min_acc_exp = Machine::PosInfinity;
  double min_exp;
  grid inf_grd, sup_grd;
  BOX rect(2);
  iterate it;
  List<iterate> Iterate_List, Image_List;
 
  Iterate_List += f_1_it;
  while ( !IsEmpty(Iterate_List) )
    {
      // First, we get min_exp and image hull.
      Get_Exp_And_Image_Hull(inf_grd, sup_grd, min_exp, Iterate_List);
      acc_exp *= min_exp;

      // Find the elements of the iterate's image. Store in Image_List.
      Empty(Image_List); // Clear the list.
      Compute_Image(Iterate_List, Image_List, Total_List, E_List);
      Empty(Iterate_List); // Clear the list.
  
      if ( show_info )
	cout << "   |Image_List|: " << Length(Image_List) << endl << endl;

      // Check if any elements of Image_List belong to F_0_List (or its twin).
      if ( IsEmpty(Image_List) )
	{ cout << "; Empty image!"; cout.flush(); }
      else
	{
	  First(Image_List);
	  while( !IsEmpty(Image_List) && !Finished(Image_List) )
	    {
	      it = Current(Image_List);
	      if ( Intersect(it, F_0_List) ) // Uses symmetry.
		{ // If we enter the fundamental domain...
		  if ( acc_exp < 2.0 )
		    { // ..and if we have not enough expansion.
		      cout << "Iterate " << f_1_it.ndl.grd << " only accumulates "
			   << acc_exp << " in expansion." << endl << endl; 
		      getchar();
		      failure = true;
		      break;
		    }
		  else
		    { // ..and if we DO have enough expansion.
		      min_acc_exp = Min(min_acc_exp, acc_exp);
		      if ( show_info )
			cout << "Iterate " << it.ndl.grd << " accumulated "
			     << acc_exp << " in expansion." << endl; 
		    }
		}
	      else // If we do not enter the fundamental domain...
		Iterate_List += it;
	      Next(Image_List);
	    }
	}
    }

  return min_acc_exp;
}

////////////////////////////////////////////////////////////////////

// Loops through Source_List and computes the image of each element. 
// These are filtered for redundancies, and stored in Image_List. 
static void Compute_Image(List<iterate> &Source_List, List<iterate> &Image_List,
			  List<iterate> &Total_List, List<iterate> &E_List)
{
  bool hit;
  iterate source_it, return_it, image_it;
  List<iterate> Return_List, Add_List;

  First(Source_List);
  while( !Finished(Source_List) )
    {
      source_it = Current(Source_List);

      Get_Return_Grids(Return_List, Total_List, E_List, source_it, false); // Does NOT use symmetry.

      if ( IsEmpty(Image_List) )
	Append(Image_List, Return_List); // Empties Return_List

      // We only add new elements to Image_List.
      while( !IsEmpty(Return_List) ) 
	{
	  return_it = First(Return_List);
	  RemoveCurrent(Return_List);
	  hit = false;
	  First(Image_List);
	  while( !Finished(Image_List) )
	    {
	      image_it = Current(Image_List);
	      if ( Equal(image_it.ndl.grd, return_it.ndl.grd) ) // Does NOT use symmetry.
		{
		  hit = true;
		  break; // Get next return_it.
		}
	      Next(Image_List);
	    }
	  if ( hit == false )
	    Add_List += return_it;
	}
      Append(Image_List, Add_List);
      Next(Source_List);
    }
}

////////////////////////////////////////////////////////////////////

// Loops through Source_List and filteres it from redundancies.
static void Remove_Redundancies(List<iterate> &Red_List)
{
  bool hit;
  iterate image_it, curr_it;
  List<iterate> Image_List;

  if ( IsEmpty(Red_List) )
    return;

  Image_List += First(Red_List); // The first element can't be redundant.
  RemoveCurrent(Red_List);

  // We only add new elements to Image_List.
  while( !IsEmpty(Red_List) ) 
    {
      curr_it = First(Red_List);
      RemoveCurrent(Red_List);
      hit = false;
      First(Image_List);
      while( !Finished(Image_List) )
	{
	  image_it = Current(Image_List);
	  if ( Equal(image_it.ndl.grd, curr_it.ndl.grd) ) // Does NOT use symmetry.
	    {
	      hit = true;
	      break; // Get next return_it.
	    }
	  Next(Image_List);
	}
      if ( hit == false )
	Image_List += curr_it;
    }
  Append(Red_List, Image_List);
}

////////////////////////////////////////////////////////////////////

// Mirrors all lower elements of Vol_List to the upper branch.
static void Reflect_To_Upper(List<iterate> &Vol_List)    
{
  grid          temp_grd;
  iterate       curr_it;
  List<iterate> Upper_List;
  
  while( !IsEmpty(Vol_List) ) // Empties Vol_List.
    {
      curr_it = First(Vol_List);
      RemoveCurrent(Vol_List);
      if ( 5 * curr_it.ndl.grd.v < 2 * curr_it.ndl.grd.u )
	{
	  curr_it.ndl.grd.v = - curr_it.ndl.grd.v;
	  curr_it.ndl.grd.u = - curr_it.ndl.grd.u;
	  temp_grd.u = curr_it.inf_grd.u;
	  temp_grd.v = curr_it.inf_grd.v;
	  curr_it.inf_grd.u = - curr_it.sup_grd.u;
	  curr_it.inf_grd.v = - curr_it.sup_grd.v;
	  curr_it.sup_grd.u = - temp_grd.u;
	  curr_it.sup_grd.v = - temp_grd.v;
	}
      Upper_List += curr_it;
    }
  // Now that Vol_List is empty, we fill 
  // it with the contents of Upper_List.
  Append(Vol_List, Upper_List);
}

////////////////////////////////////////////////////////////////////

// Mirrors a copy of all elements of Vol_List to the lower branch.
static void Symmetrize(List<iterate> &Vol_List)    
{
  grid          temp_grd;
  iterate       curr_it;
  List<iterate> Symmetric_List;
  
  while( !IsEmpty(Vol_List) ) // Empties Vol_List.
    {
      curr_it = First(Vol_List);
      RemoveCurrent(Vol_List);
      Symmetric_List += curr_it;
      curr_it.ndl.grd.v = - curr_it.ndl.grd.v;
      curr_it.ndl.grd.u = - curr_it.ndl.grd.u;
      temp_grd.u = curr_it.inf_grd.u;
      temp_grd.v = curr_it.inf_grd.v;
      curr_it.inf_grd.u = - curr_it.sup_grd.u;
      curr_it.inf_grd.v = - curr_it.sup_grd.v;
      curr_it.sup_grd.u = - temp_grd.u;
      curr_it.sup_grd.v = - temp_grd.v;
      Symmetric_List += curr_it;
    }
  // Now that Vol_List is empty, we fill 
  // it with the contents of Upper_List.
  Append(Vol_List, Symmetric_List);
}

////////////////////////////////////////////////////////////////////

// Removes all elements from In_List that have a match in Fixed_List.
static void Remove_Intersections(List<iterate> &In_List, List<iterate> &Fixed_List)
{
  iterate fixed_it, in_it;
  
  First(Fixed_List);
  while( !Finished(Fixed_List) )
    {
      fixed_it = Current(Fixed_List);
      First(In_List);
      while( !IsEmpty(In_List) && !Finished(In_List) )
	{
	  in_it = Current(In_List);
	  if ( Equal(in_it.ndl.grd, fixed_it.ndl.grd) ) // Does NOT use symmetry.
	    RemoveCurrent(In_List);
	  else
	    Next(In_List);
	}
      Next(Fixed_List);
    }
}

////////////////////////////////////////////////////////////////////

// Loops through all of Total_List and singles out the elements that belong
// to the return of 'it'. These are passed by reference in Return_List. 
// If it \in E_List, we remove all its images that also are in E_List.
// Assumes that Return_List is empty at start.
static void Get_Return_Grids(List<iterate> &Return_List, List<iterate> &Total_List, 
			     List<iterate> &E_List, const iterate &it, bool upper)
{
  bool E_iterate = false;
  iterate curr_it;
  grid temp_grd = NULL_GRID;

  if( !IsEmpty(Return_List) )
    {
      cout << "Error (Get_Return_Grids): Return_List not empty!" << endl;
      exit(1);
    }

  if ( Intersect(it, E_List) )
    E_iterate = true;

  First(Total_List);
  while( !Finished(Total_List) )
    {
      curr_it = Current(Total_List);
      if ( it.inf_grd.u <= curr_it.ndl.grd.u && curr_it.ndl.grd.u <= it.sup_grd.u )
	if ( it.inf_grd.v <= curr_it.ndl.grd.v && curr_it.ndl.grd.v <= it.sup_grd.v )
	  {
	    if ( upper ) // Reflect up if needed.
	      if ( 5 * curr_it.ndl.grd.v < 2 * curr_it.ndl.grd.u )
		{
		  curr_it.ndl.grd.v = - curr_it.ndl.grd.v;
		  curr_it.ndl.grd.u = - curr_it.ndl.grd.u;
		  temp_grd.u = curr_it.inf_grd.u;
		  temp_grd.v = curr_it.inf_grd.v;
		  curr_it.inf_grd.u = - curr_it.sup_grd.u;
		  curr_it.inf_grd.v = - curr_it.sup_grd.v;
		  curr_it.sup_grd.u = - temp_grd.u;
		  curr_it.sup_grd.v = - temp_grd.v;
		}
	    if ( E_iterate && Intersect(curr_it, E_List) )
	      {} // Don't save the image.
	    else
	      Return_List += curr_it;
	  }
      Next(Total_List);
    }
}

////////////////////////////////////////////////////////////////////

// Appends the contents of Add_List to Out_List, and empties Add_List.
static void Append(List<iterate> &Out_List, List<iterate> &Add_List)
{
  while( !IsEmpty(Add_List) )
    {
      Out_List += First(Add_List);
      RemoveCurrent(Add_List);
    }
}

////////////////////////////////////////////////////////////////////

// Copies the contents of Original_List to Copy_List.
static void Copy(List<iterate> &Original_List, List<iterate> &Copy_List)
{
  if( IsEmpty(Original_List) )
    return;
  
  Empty(Copy_List);
  First(Original_List);
  while( !Finished(Original_List) )
    {
      Copy_List += Current(Original_List);
      Next(Original_List);
    }
}

////////////////////////////////////////////////////////////////////

bool Equal(const grid &grd_1, const grid &grd_2)
{
  if ( grd_1.u == grd_2.u && grd_1.v == grd_2.v )
    return true;
  else
    return false;
}

////////////////////////////////////////////////////////////////////

bool Intersect(const iterate &it, List<iterate> &it_List)
{
  if ( IsEmpty(it_List) )
    return false;

  First(it_List);
  while( !Finished(it_List) )
    {
      if ( it.ndl.grd == Current(it_List).ndl.grd )
	return true;
      Next(it_List);
    }
  return false;
}

////////////////////////////////////////////////////////////////////

static void Empty(List<iterate> &it_List)
{
  if ( IsEmpty(it_List) )
    return;

  First(it_List);
  while ( !IsEmpty(it_List) )
    RemoveCurrent(it_List);
}

////////////////////////////////////////////////////////////////////

static void Get_Exp_And_Image_Hull(grid &inf_grd, grid &sup_grd, double &min_exp, 
				   List<iterate> &Iterate_List)
{
  grid min_grd = First(Iterate_List).inf_grd;
  grid max_grd = First(Iterate_List).sup_grd;
  int min_u = min_grd.u;
  int max_u = max_grd.u;
  int min_v = min_grd.v;
  int max_v = max_grd.v;
  inf_grd.P = First(Iterate_List).ndl.grd.P;
  sup_grd.P = inf_grd.P;
  min_exp   = First(Iterate_List).ndl.min_exp;
  Next(Iterate_List);
  while( !Finished(Iterate_List) )
    {
      min_grd = Current(Iterate_List).inf_grd;
      max_grd = Current(Iterate_List).sup_grd;
      if ( min_u > min_grd.u )
	min_u = min_grd.u;
      if ( max_u < max_grd.u )
	max_u = max_grd.u;
      if ( min_v > min_grd.v )
	min_v = min_grd.v;
      if ( max_v < max_grd.v )
	max_v = max_grd.v;
      min_exp = Min(min_exp, Current(Iterate_List).ndl.min_exp);
      Next(Iterate_List);
    }
  inf_grd.u = min_u;
  inf_grd.v = min_v;
  sup_grd.u = max_u;
  sup_grd.v = max_v;
}

////////////////////////////////////////////////////////////////////










