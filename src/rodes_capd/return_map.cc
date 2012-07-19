/*   File: return_map.cc (BIAS)

     Contains all high level functions
     needed to compute a return. 

     Latest edit: Mon Apr 10 2000 
*/

#include "return_map.h"

////////////////////////////////////////////////////////////////////

static bool   Too_Large         (const parcel &, const double &);
static bool   Switching         (const parcel &,       short  &, const double &);
static void   Build_Switch_Box  (const parcel &, const short  &,       BOX    &, const double &);
static bool   Switch_Box_True   (const parcel &,       short  &, const BOX    &);
static bool   Stop              (const parcel &, const double &, const stop_parameters &);
static void   Flow              (      parcel &, const double &);
static void   Update_Transversal(      parcel &, const short  &, const double &,
				 const double &, const double &, const double &);
static void   Flow_The_Parcel   (const parcel &, List<parcel> &,
				 const stop_parameters &, const double &, const double &);
static void   Set_Max_Size      (      double &,       double &, const parcel &);

////////////////////////////////////////////////////////////////////

// Called by 'Flow_The_Parcel' and 'Multiple_Partition'.
// Checks if a box has grown too large.
static bool Too_Large(const parcel &pcl, const double &size)
{
  for (register short i = 1; i <= SYSDIM; i++)
    if ( i != pcl.trvl )
      if ( Sup(pcl.box(i)) - Inf(pcl.box(i)) > size )
	return true;

  return false;
}

////////////////////////////////////////////////////////////////////

// Called by 'Flow_The_Parcel'.
// Splits the box of pcl into
// pieces of side length in [0.5, 1] * size
void Multiple_Partition(const parcel &pcl, List<parcel> &PSinit,
			const double &size)
{
  bool no_split = false;    // true when no BOX was split
  int  n_in;                // Number of elements of PSinit
  parcel tmp_pcl;           // Temporary storage 
#if defined(__sparc)
  List<parcel> DummyList;
#endif

  PSinit += pcl;
  while ( !no_split )
    { // At least one BOX was split
      n_in = Length(PSinit);
      no_split = true;
      for(register int i = 0; i < n_in; i++)
	{ // Loop through PSinit
	  tmp_pcl = First(PSinit);
	  RemoveCurrent(PSinit);
	  if ( Too_Large(tmp_pcl, size) )
	    {
	      Single_Partition(tmp_pcl, PSinit, size);
	      no_split = false;
	    }
	  else // put it back at the end
	    PSinit += tmp_pcl;
	}
    }
}

////////////////////////////////////////////////////////////////////

// Checks if we have a potential switching situation.
static bool Switching(const parcel &pcl, short &trvl, const double &more_than_one)
{
  double factor = more_than_one; // Prevents flipping between leading directions;
  BOX vf(SYSDIM);  Vf_Range(vf, pcl.box);// Can be NR = non-rigorous.

  trvl = pcl.trvl;
  for ( register short i = 1; i <= SYSDIM; i++ )  
    if ( i != pcl.trvl )                       
      if ( Mig(vf(i)) > factor * Abs(vf(trvl)) )  // <-- abs ( interval )
	{ 
	  trvl = i;     // Check which component of the vector field 
	  factor = 1.0; // is dominating, and get the right sign
	}                 
  if ( pcl.trvl == trvl ) // If we are not switching, return 'false'
    return false;
  return true; 
}

////////////////////////////////////////////////////////////////////

// If we have a potential new transversal,
// we build a switch box - Switch_Box.
static void Build_Switch_Box(const parcel &pcl, const short &trvl, BOX &Switch_Box,
			     const double &more_than_one)
{
  double level = Sup(pcl.box(pcl.trvl)); // Inf == Sup.
  double diam  = Sup(pcl.box(trvl)) - Inf(pcl.box(trvl));

  Switch_Box = pcl.box;                       
  for ( register short i = 1; i <= SYSDIM; i++ ) 
    if ( (i != trvl) && (i != pcl.trvl) )       // Widen the non-essential sides
      Switch_Box(i) = Hull(Inf(Switch_Box(i)) - diam, Sup(Switch_Box(i)) + diam);

  diam /= more_than_one;
  if ( pcl.sign == 1 )
    Switch_Box(pcl.trvl) = Hull(level, level + diam);
  else
    Switch_Box(pcl.trvl) = Hull(level - diam, level);
}

////////////////////////////////////////////////////////////////////

// Return 'false' if we don't stay inside Switch_Box.
// Otherwise, return 'true', and pass the new direction
// via the reference to 'trvl'.
static bool Switch_Box_True(const parcel &pcl, short &trvl, const BOX &Switch_Box)
{
  BOX vf(SYSDIM); Vf_Range(vf, Switch_Box);              

  for ( register short i = 1; i <= SYSDIM; i++ ) 
    if ( i != trvl )
      if ( Abs(vf(i)) > Mig(vf(trvl)) )
	return false;
  // If we stay inside Switch_Box, we store the new direction
  // in trvl, which is passed by reference, and return the value 'true'.
  if ( Sign(vf(trvl)) == -1 )
    trvl = - trvl;             // trvl \in {+-1,..., +-SYSDIM}
  return true;                    
}

////////////////////////////////////////////////////////////////////

// Called by 'Flow_The_Parcel'.
// Checks if the parcel is near the stopping condition
static bool Stop(const parcel &pcl, const double &dist, const stop_parameters &sp)
{
  if ( pcl.trvl == sp.trvl )
    if ( pcl.sign == sp.sign )
      {
	double diff = Sup(pcl.box(pcl.trvl)) - sp.level; // Inf = Sup
	if ( fabs(diff) < dist )
	  {
	    if ( (sp.sign == - 1) && (diff > 0.0) )
	      return true;
	    if ( (sp.sign == + 1) && (diff < 0.0) )
	      return true;
	  }
      }
  return false;
}

////////////////////////////////////////////////////////////////////

// Called by 'Flow_The_Parcel'.
// Flows the parcel as far as possible, but
// not by more than 'trvl_dist' at a time.
static void Flow(parcel &pcl, const double &trvl_dist)
{
  parcel result = pcl;  // Pass on the unchanged pieces by copying.
  BOX Outer_Box(SYSDIM);

  double mid, rad;   
  for ( register short i = 1; i <= SYSDIM; i++ )
    if ( i == pcl.trvl ) // Widen trvl direction.
      {
	if ( pcl.sign == 1 ) // Inf == Sup
	  Outer_Box [i] = Hull ( Inf ( pcl.box[i] ), pcl.box[i] + trvl_dist );	
	else
	  Outer_Box[i] = Hull ( Inf( pcl.box[i] ) - trvl_dist, pcl.box[i]);	
      }
    else // i != pcl.trvl
      {
	rad = (Sup(pcl.box(i)) - Inf(pcl.box(i))) / 2.0;
	mid = Inf(pcl.box(i)) + rad;
	rad *= SCALE_FACTOR;
	Outer_Box(i) = Hull(mid - rad, mid + rad);
      }

  interval time;
  Get_Flow_Time(time, result, trvl_dist, Outer_Box);

  // Now, we tighten the enclosure...
  IMatrix DPi(SYSDIM, SYSDIM);
  BOX Tight_Box;

  BOX Image = pcl.box + time * Vf_Range(Outer_Box);
  if ( pcl.sign == 1 )
    Image(pcl.trvl) = Hull(Sup(Outer_Box(pcl.trvl)));
  else
    Image(pcl.trvl) = Hull(Inf(Outer_Box(pcl.trvl)));
  Get_DPi_Matrix(DPi, Outer_Box, pcl.trvl, time, Image);

  Flow_By_Corner_Method(Tight_Box, DPi, pcl, Outer_Box);

#ifdef COMPUTE_C1  // ...and flow the tangent vectors
  Flow_Tangent_Vectors(result, pcl.sign*pcl.trvl, pcl.sign*pcl.trvl, DPi);
#endif

  result.box = Tight_Box;  // Update the outgoing result
  result.time += time;
  pcl = result;

    cout << "pcl.box = " << pcl.box << endl;

}

////////////////////////////////////////////////////////////////////

// Called by 'Flow_The_Parcel'.
// Changes the local Poincare' planes according to the dominating component
// of the vector field, and flows the parcel to the first new plane.
static void Update_Transversal(parcel &pcl, const short &trvl, const double &max_size, 
			       const double &max_size_over_three, 
			       const double &max_size_over_ten, const double &more_than_one)
{
  parcel result  = pcl;
  short new_trvl = abs(trvl);  
  short new_sign = Sign(trvl);
  stop_parameters stop_pmtr;

  stop_pmtr.max_d_step = max_size_over_ten; // Set the stop parameters
  stop_pmtr.trvl  = new_trvl;
  stop_pmtr.sign  = new_sign;             
  if ( new_sign == 1 )
    stop_pmtr.level = Sup(pcl.box(new_trvl));
  else
    stop_pmtr.level = Inf(pcl.box(new_trvl));

  // Split pcl into several small pieces -> Start_List, 
  // whose elements are flowed separately to the plane.
  List<parcel> Start_List, Stop_List, Image_List;
  parcel current_pcl;
  BOX Temp_Switch_Box(SYSDIM);

  Multiple_Partition(pcl, Start_List, max_size_over_three);
  First(Start_List);
  while ( !Finished(Start_List) )
    { // Work through Start_List
      current_pcl = Current(Start_List);
      Build_Switch_Box(current_pcl, new_trvl, Temp_Switch_Box, more_than_one);	
      Switch_Transversal(current_pcl, trvl, Temp_Switch_Box);

      if ( stop_pmtr.level == Sup(current_pcl.box(new_trvl)) ) // Sup = Inf
	Stop_List += current_pcl; //No need to flow
      else   // Here we make a call to 'Local_Flow_The_Parcel', which flows
	{    // the parcel straight to the new trvl plane.  
	  Local_Flow_The_Parcel(current_pcl, Image_List, stop_pmtr, max_size); 
	  while ( !IsEmpty(Image_List) )                       
	    { // Work through Image_List, and delete it.    
	      Stop_List += First(Image_List);
	      RemoveCurrent(Image_List);
	    }
	}
      Next(Start_List);
    }
  Get_Hull(result, Stop_List);
  result.message = 0;

  pcl = result;
}

////////////////////////////////////////////////////////////////////

// Called by 'Compute_the_return'.
// Flows in_pcl until the stopping condition is met.
// All parcels constituting the pseudo-image are returned by reference
// in Return_List. Any error in a lower function is handled
// by catching (in 'rodes') the thrown exception.
static void Flow_The_Parcel(const parcel &in_pcl, List<parcel> &Return_List,
			    const stop_parameters &sp, const double &max_size,
			    const double &more_than_one)
{
  short        trvl;     // The new direction when switching.
  double       dist;     // The trvl distance we attempt to flow
  BOX          sw_box;   // The box used for switching transversals
  parcel       pcl;      // The parcel from In_List under computation
  List<parcel> In_List;  // All intermediate images of in_pcl

  // Special constants used in 'Update_Transversal'.
  double max_size_over_three = max_size / 3.00;
  double max_size_over_ten   = max_size / 10.0;

  In_List += in_pcl; 
  
  while( !IsEmpty(In_List) )
    {  // Loop through all of In_List  
      pcl = First(In_List);
      RemoveCurrent(In_List);
      dist = sp.max_d_step;

      while (1) // Enter the flow loop
	{ 
	  if ( Too_Large(pcl, max_size) ) // If the box is too large, we 
	    {                             // partition it sufficiently.
	      Multiple_Partition(pcl, In_List, max_size);
	      break;
	    }
	  if ( Switching(pcl, trvl, more_than_one) ) // If we have a possible switching
	    {                                        // situation, we attempt to switch.
	      Build_Switch_Box(pcl, trvl, sw_box, more_than_one);
	      if ( Switch_Box_True(pcl, trvl, sw_box) )	
		{
		  Update_Transversal(pcl, trvl, max_size, max_size_over_three,
				     max_size_over_ten, more_than_one);
		  In_List += pcl; 
		}
	      else 
		// jjb -- Need Max of { [a,a], [b,b], [c,c] } to get max radius
		Multiple_Partition ( pcl, In_List, Max ( diam ( pcl.box ) ) / 2.0 );
	      break;                             
	    }
	  if ( pcl.message == STOP ) // If we we have completed a full
	    {                        // return, we store the parcel.
	      Return_List += pcl;
	      break;
	    }
	  if ( Stop(pcl, sp.max_d_step, sp) ) // If we are close to a return
	    {                                 // we decrease the trvl_dist.
	      dist = fabs(Sup(pcl.box(pcl.trvl) - sp.level));
	      pcl.message = CLOSE_STOP;
	    }
	  if ( Inf(Norm2(pcl.box)) < 1.0 ) // Improves accuracy near the fixed point.       
	    dist = Min(sp.max_d_step, Inf(Norm2(pcl.box)) / 10.0);
	  if ( Cube_Entry(pcl) ) // If we enter the cube containing the origin
	    {                    // we explicitly compute the outgoing image(s).
	      Cube_Exit(pcl, In_List, max_size);
	      break;
	    }
	  Flow(pcl, dist); // If none of the situations        
	}                  // above occured, we flow along.
    }
}

////////////////////////////////////////////////////////////////////

// Called by 'Update_Transversal' and 'Flatten_Parcels'.
// Flows in_pcl until the stopping condition is met. We do NOT allow
// for trvl switching or cube entries during the flow time.
// All parcels constituting the pseudo-image are returned by reference
// in Return_List. Any error in a lower function is handled
// by catching (in 'rodes') the thrown exception.
void Local_Flow_The_Parcel(const parcel &in_pcl, List<parcel> &Return_List,
			   const stop_parameters &sp, const double &max_size)
{
  double       dist;     // The trvl distance we attempt to flow
  BOX          sw_box;   // The box used for switching transversals
  parcel       pcl;      // The parcel from In_List under computation
  List<parcel> In_List;  // All intermediate images of in_pcl

  In_List += in_pcl; 
  
  while( !IsEmpty(In_List) )
    {  // Loop through all of In_List  
      pcl = First(In_List);
      RemoveCurrent(In_List);
      dist = sp.max_d_step;

      while (1) // Enter the flow loop
	{ 
	  if ( Too_Large(pcl, max_size) ) // If the box is too large, we 
	    {                             // partition it sufficiently.
	      Multiple_Partition(pcl, In_List, max_size);
	      break;
	    }
	  if ( pcl.message == STOP ) // If we we have completed a full
	    {                        // return, we store the parcel.
	      Return_List += pcl;
	      break;
	    }
	  if ( Stop(pcl, sp.max_d_step, sp) ) // If we are close to a return
	    {                                 // we decrease the trvl_dist.
	      dist = fabs(Sup(pcl.box(pcl.trvl) - sp.level));
	      pcl.message = CLOSE_STOP;
	    }
	  Flow(pcl, dist); // If none of the situations        
	}                  // above occured, we flow along.
    }
}

////////////////////////////////////////////////////////////////////

// Sets the resolution to be used in 'Flow_the_parcel'. This number
// depends on the x_1-distance to the stable manifold.
// Also sets the variable 'more_than_one'.
static void Set_Max_Size(double &max_size, double &more_than_one, const parcel &pcl)
{
  double x = Mid(pcl.box(1));
  double y = Mid(pcl.box(2));

  more_than_one = 1.1; // The default value used almost everywhere,

  if ( 5 * y < 2 * x ) // Reflect onto upper branch.
    x = - x;

  if ( x > + 4.5 )        // Zone #1  [+4.5, +inf]
    max_size = 0.010; 
  else if ( x > + 4.0 )   // Zone #2  [+4.0, +4.5]
    max_size = 0.012;
  else if ( x > + 3.5 )   // Zone #3  [+3.5, +4.0]
    max_size = 0.015;
  else if ( x > + 3.0 )   // Zone #4  [+3.0, +3.5]
    max_size = 0.019; 
  else if ( x > + 2.5 )   // Zone #5  [+2.5, +3.0]
    max_size = 0.024; 
  else if ( x > + 2.0 )   // Zone #6  [+2.0, +2.5]
    max_size = 0.030; 
  else if ( x > + 1.5 )   // Zone #7  [+1.5, +2.0]
    max_size = 0.032; 
  else if ( x > + 1.0 )   // Zone #8  [+1.0, +1.5]
    max_size = 0.035; 
  else if ( x > + 0.5 )   // Zone #9  [+0.5, +1.0] - contains W^s(0).
    max_size = 0.040; 
  else if ( x > + 0.0 )   // Zone #10 [+0.0, +0.5]
    max_size = 0.035; 
  else if ( x > - 1.0 )   // Zone #11 [-1.0, +0.0]
    max_size = 0.032; 
  else if ( x > - 2.0 )   // Zone #12 [-2.0, -1.0]
    max_size = 0.029; 
  else if ( x > - 3.0 )   // Zone #13 [-3.0, -2.0]
    max_size = 0.026; 
  else if ( x > - 4.0 )   // Zone #14 [-4.0, -3.0]
    max_size = 0.023;    
  else                    // Zone #15 [-4.296875, -4.0]
    max_size = 0.020;  

  if ( x < - 1100 / 256.0 )   // Zone #16 [-4.36328125, -4.296875]     
    {
      more_than_one = 1.2; // Here we use a slightly larger value.
      max_size = 0.015; 
    }
  if ( x < - 1118 / 256.0 )   // Zone #17 [-inf, -4.36328125]   
    {
      more_than_one = 1.3; // Here we use a slightly larger value.
      max_size = 0.012; 
    } 

  cout << "max_size = " << max_size << "; "; 
  cout << "more_than_one = " << more_than_one << "; "; 
  cout.flush();
}

////////////////////////////////////////////////////////////////////

// Called by: 'Work_on_grid'.
// Calls to : 'Set_Max_Size' and 'Flow_the_parcel'.
// Computes the image set of the parcel w.r.t. the global stopping 
// parameters, "glob_stop_param". The image set is returned by reference
// via Return_List. If any error occurs, an exception is thrown to be 
// caught in the calling function.
void Compute_the_return(const parcel &current_pcl, List<parcel> &Return_List)
{
  stop_parameters glob_stop_param;
  double max_size, max_dist_step, more_than_one;

  // Set the resolution depending on where we are.
  Set_Max_Size(max_size, more_than_one, current_pcl);
  max_dist_step = 10.0 * max_size;

  // Set the global stopping parameters
  glob_stop_param.level      = STOP_DIST_LEVEL;
  glob_stop_param.trvl       = STOP_TRANSVERSAL;
  glob_stop_param.sign       = STOP_SIGN;
  glob_stop_param.max_d_step = max_dist_step;

  Flow_The_Parcel(current_pcl, Return_List, glob_stop_param, max_size, more_than_one);
}

////////////////////////////////////////////////////////////////////
