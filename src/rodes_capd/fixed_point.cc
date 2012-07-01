/*   File: new_fixed_point.cc (BIAS)
 
     Contains all functions needed to deal
     with parcels entering the cube containing
     the fixed point at the origin.

     Latest edit: Mon May 30 2000

*/

#include "fixed_point.h"

// Parameters defining the cube.
static const double CUBE_RADIUS = 0.1;
static const INTERVAL SYMM_CUBE_RADIUS = SymHull(CUBE_RADIUS);

// Parameters for computing the exit.
static const INTERVAL PM_ONE_IV = SymHull(1.0);
static const INTERVAL MU = - E2_IV / E1_IV;
static const INTERVAL NU = - E3_IV / E1_IV;

////////////////////////////////////////////////////////////////////

static void Divide               (const parcel &, List<parcel> &, const double &, int);
static void Divide               (List<parcel> &, List<parcel> &, const double &, int);
static void Deform               (      parcel &);
static void Widen_via_entrance   (List<parcel> &, List<parcel> &);
static void Check_for_splittings (List<parcel> &, List<parcel> &);
static void Compute_exits        (List<parcel> &, List<parcel> &);
static void Compute_single_exit  (const parcel &, List<parcel> &);
static void Widen_via_exit       (List<parcel> &, List<parcel> &);
static void Flatten_Parcels      (List<parcel> &, List<parcel> &, const double &);

////////////////////////////////////////////////////////////////////

// Called by: 'Flow_The_Parcel'.
// Calls to : none
// Checks if the parcel enters the cube containing the origin.
bool Cube_Entry(const parcel &pcl)
{
  if ( pcl.sign == - 1 )
    if ( pcl.trvl == 3 )
      if ( Sup(pcl.box(3)) < CUBE_RADIUS ) // Inf == Sup
	{
	  if ( Abs(pcl.box(1)) < CUBE_RADIUS )
	    {
	      if ( Abs(pcl.box(2)) < CUBE_RADIUS )    
		return true; // pcl.box is completely inside the cube.
	    }
	  else if ( Subset(0.0, pcl.box(1)) )
	    return true; // pcl.box straddles the stable manifold.
	}
  return false;
}
////////////////////////////////////////////////////////////////////

// Called by: 'Cube_Exit', 'Flatten_Parcels'.
// Calls to : 'Multiple_Partition'.
static void Divide(const parcel &pcl, List<parcel> &Result_List, 
		   const double &max_size, int N)
{
  Multiple_Partition(pcl, Result_List, max_size / N);
}

////////////////////////////////////////////////////////////////////

// Called by: 'Cube_Exit'.
// Calls to : 'Divide'.
static void Divide(List<parcel> &In_List, List<parcel> &Result_List, 
		   const double &max_size, int N)
{
  List<parcel> dummy_List;

  parcel cur_pcl;
  First(In_List);
  while ( !Finished(In_List) )
    {
      cur_pcl = Current(In_List);
      Divide(cur_pcl, Result_List, max_size, N);
      Next(In_List);
    }
}

////////////////////////////////////////////////////////////////////

// Called by: 'Cube_Exit'.
// Calls to : 'Deform'.
// Returns a list of enclosures of the deformed incoming parcels.
static void Widen_via_entrance(List<parcel> &In_List, List<parcel> &Result_List)
{
  parcel cur_pcl;

  First(In_List);
  while ( !Finished(In_List) )
    {
      cur_pcl = Current(In_List);
      Deform(cur_pcl);
      Result_List += cur_pcl;
      Next(In_List);
    }
}

////////////////////////////////////////////////////////////////////

// Called by: 'Cube_Exit'.
// Calls to : none.
static void Check_for_splittings(List<parcel> &In_List, List<parcel> &Result_List)
{
  parcel cur_pcl;

  First(In_List);
  while ( !Finished(In_List) )
    {
      cur_pcl = Current(In_List);
      if ( Subset(0.0, cur_pcl.box(1)) )
	{
	  parcel left_pcl  = cur_pcl;
	  parcel right_pcl = cur_pcl;
	  
	  left_pcl.box(1)  = Hull(Inf(cur_pcl.box(1)), 0.0);
	  right_pcl.box(1) = Hull(0.0, Sup(cur_pcl.box(1)));
	  Result_List += left_pcl;
	  Result_List += right_pcl;
	}
      else
	Result_List += cur_pcl;

      Next(In_List);
    }
}

////////////////////////////////////////////////////////////////////

// Called by: 'Cube_Exit'.
// Calls to : 'Flow_Tangent_Vectors'.
// Returns (by reference) an enclosure of the deformed parcel.
static void Deform(parcel &pcl)
{
  short trvl = pcl.trvl;
  INTERVAL distance = Hull(Abs(pcl.box(trvl)));

  if ( trvl == 1 ) // Use a larger distance for the inverse.
    distance = 1 - Sqrt(1 - 2 * distance); 

  INTERVAL C0_DEFORM_TERM = PM_ONE_IV * Sqr(distance) / 2.0; // [- r^2 / 2, + r^2 / 2]

  for ( short i = 1; i <= DIM; i++ )
    pcl.box(i) += C0_DEFORM_TERM;

#ifdef COMPUTE_C1
  INTERVAL C1_NORM = 2.0 * distance;
  if ( trvl == 1 ) // Use the formula for the inverse.
    C1_NORM /= (1 - C1_NORM);
  
  INTERVAL C1_DEFORM_FACTOR = 1 + PM_ONE_IV * C1_NORM; 
  INTERVAL_MATRIX D_M(DIM, DIM);  Clear(D_M); // Diagonal matrix D_M.

  for (register short i = 1; i <= DIM; i++) // D_M(i, i) == C1_DEFORM_FACTOR,
    D_M(i, i) = C1_DEFORM_FACTOR;           // D_M(i, j) == [0, 0] if i != j.

  Flow_Tangent_Vectors(pcl, trvl, trvl, D_M);
#endif // COMPUTE_C1
}

////////////////////////////////////////////////////////////////////

// Called by: 'Cube_Exit', 'Compute_multiple_exit'.
// Calls to : 'Flow_Tangent_Vectors'.
static void Compute_single_exit(const parcel &in_pcl, List<parcel> &Result_List)
{
  INTERVAL x = in_pcl.box(1);
  INTERVAL y = in_pcl.box(2);
  INTERVAL z = in_pcl.box(3);
  INTERVAL abs_x = Hull(fabs(Inf(x)), fabs(Sup(x)));
  int sign_x = (Mid(x) > 0.0 ? 1 : -1); // Distinguish a left and a right exit.
  double exit_rad = CUBE_RADIUS;
  INTERVAL abs_x_over_exit_rad =  abs_x / exit_rad;
  parcel result = in_pcl;

  result.box(1) = sign_x * exit_rad;
  result.box(2) = y * Power(abs_x_over_exit_rad, MU);
  result.box(3) = z * Power(abs_x_over_exit_rad, NU); 

  // Has the parcel been split over W^ss(0) or not?
  bool splitting = false;
  if ( Inf(in_pcl.box(1)) == 0.0 || Sup(in_pcl.box(1)) == 0.0 )
    splitting = true;

  if ( splitting )
    { // Since Log10 and Power can not handle zeroes very well,
      // we pre-compute these results, and insert them via Hull.
      INTERVAL max_x_over_r = Hull(Sup(abs_x_over_exit_rad));
      double small_time = Inf(- 1.0 / E1_IV * Log10(max_x_over_r)); 
      result.time += Hull(small_time, Machine::PosInfinity);
#ifdef COMPUTE_C1
      INTERVAL_MATRIX P_M(DIM, DIM);  Clear(P_M); // Poincare map matrix.

      P_M(2, 1) = Hull(0.0, MU / exit_rad * y * sign_x * Power(max_x_over_r, MU - 1)); 
      P_M(2, 2) = Hull(0.0, Power(max_x_over_r, MU));
      P_M(3, 1) = sign_x * Hull(NU / exit_rad * z * Power(max_x_over_r, NU - 1), Machine::PosInfinity); 
      P_M(3, 3) = Hull(0.0, Power(max_x_over_r, NU)); // This element is never used.

      Flow_Tangent_Vectors(result, result.trvl, sign_x, P_M);
#endif // COMPUTE_C1
    }
  else // No splitting.
    {    
      result.time += - 1.0 / E1_IV * Log10(abs_x_over_exit_rad);
#ifdef COMPUTE_C1
      INTERVAL_MATRIX P_M(DIM, DIM);  Clear(P_M); // Poincare map matrix.

      P_M(2, 1) = MU / exit_rad * y * sign_x * Power(abs_x_over_exit_rad, MU - 1); 
      P_M(2, 2) = Power(abs_x_over_exit_rad, MU);
      P_M(3, 1) = NU / exit_rad * z * sign_x * Power(abs_x_over_exit_rad, NU - 1); 
      P_M(3, 3) = Power(abs_x_over_exit_rad, NU); // This element is never used.

      Flow_Tangent_Vectors(result, result.trvl, sign_x, P_M);
#endif // COMPUTE_C1
    }
  result.trvl = 1;
  result.sign = sign_x;
  Result_List += result;
}

////////////////////////////////////////////////////////////////////

// Called by: 'Cube_Exit'.
// Calls to : 'Deform'.
// Computes enclosures of all deformed outgoing elements 
// of In_List, which are appended to Result_List.
static void Widen_via_exit(List<parcel> &In_List, List<parcel> &Result_List)
{
  parcel current_pcl;
  List<parcel> DummyList;            // Added to shut the compiler up!

  First(In_List);
  while ( !Finished(In_List) )
    {
      current_pcl = Current(In_List);
      Deform(current_pcl);
      Result_List += current_pcl;
      Next(In_List);
    }
}

////////////////////////////////////////////////////////////////////

// Called by 'Cube_Exit'.
// Calls to : 'Multiple_Partition', 'Flow_The_Parcel'.
// Flows the deformed parcels to the appropriate plane.
static void Flatten_Parcels(List<parcel> &Lumpy_List, List<parcel> &Flat_List, const double &max_size)
{
  List<parcel> Start_List, Image_List;
  parcel current_pcl;
  stop_parameters stop_pmtr; 

  First(Lumpy_List);
  while ( !Finished(Lumpy_List) )
    { // Work through Lumpy_List
      current_pcl = Current(Lumpy_List);
      stop_pmtr.max_d_step = 0.5 * Diam(current_pcl.box(1)); // Set the stop parameters
      stop_pmtr.trvl  = current_pcl.trvl;
      stop_pmtr.sign  = current_pcl.sign;             
      if ( current_pcl.sign == 1 )
	stop_pmtr.level = Sup(current_pcl.box(current_pcl.trvl));
      else
	stop_pmtr.level = Inf(current_pcl.box(current_pcl.trvl));
      
      stop_pmtr.level *= 1.1; // Flow a little bit longer.

      Divide(current_pcl, Start_List, max_size, 1);
      while ( !IsEmpty(Start_List) )
	{ // Work through Start_List, and delete it.
	  current_pcl = First(Start_List);
	  RemoveCurrent(Start_List);
	  Local_Flow_The_Parcel(current_pcl, Image_List, stop_pmtr, max_size); // ADDED JUNE 5, 2000
	  while ( !IsEmpty(Image_List) )                       
	    { // Work through Image_List, and delete it.    
	      Flat_List += First(Image_List);
	      RemoveCurrent(Image_List);
	    }
	}
      Next(Lumpy_List);
    }
}

////////////////////////////////////////////////////////////////////

// Called by: 'Cube_Exit'.
// Calls to : 'Compute_single_exit'. 
static void Compute_exits(List<parcel> &In_List, List<parcel> &Result_List)
{
  parcel cur_pcl;
  
  First(In_List);
  while ( !Finished(In_List) )
    {
      cur_pcl = Current(In_List);
      Compute_single_exit(cur_pcl, Result_List);
      Next(In_List);
    }
}

////////////////////////////////////////////////////////////////////

// Called by: 'Flow_The_Parcel'.
// Calls to : 'Divide', 'Widen_via_entrance', 'Check_for_splittings',
//            'Compute_exits', 'Widen_via_exit', 'Flatten_Parcels'.
// Explicitly computes the image(s) of the parcel entering the cube.
// Adds the image(s) to the end of the list Image_List.
void Cube_Exit(const parcel &pcl, List<parcel> &Image_List, const double &max_size)
{
  parcel hull_pcl; Resize(hull_pcl.box, DIM);
  List<parcel> Split_List, Widened_List;
  List<parcel> Pre_Enter_List, Enter_List, Exit_List;
  List<parcel> Lumpy_List, Flat_List;

  // 1;
  Divide(pcl, Split_List, max_size, 1);

  // 2;
  Widen_via_entrance(Split_List, Widened_List);

  // 3;
  Check_for_splittings(Widened_List, Pre_Enter_List);

  // 4;
  Divide(Pre_Enter_List, Enter_List, max_size, 1); 

  // 5;                                  
  Compute_exits(Enter_List, Exit_List);

  // 6;
  Widen_via_exit(Exit_List, Lumpy_List);

  // 7;
  Flatten_Parcels(Lumpy_List, Flat_List, max_size);

  // Now append Flat_List to Image_List.
  // After debugging remove the code below,
  // and call Flatten_Parcels(Lumpy_List, Image_List).
  parcel cur_pcl;
  First(Flat_List);
  while ( !IsEmpty(Flat_List) )
    { // Work through Flat_List and remove its contents.
      cur_pcl = Current(Flat_List);
      cur_pcl.message = 0;
      Image_List += cur_pcl;
      RemoveCurrent(Flat_List);
    }
}

////////////////////////////////////////////////////////////////////
