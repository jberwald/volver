/*   File: flow_functions.cc (BIAS)

     Contains all medium level functions
     needed to compute a return. 

     Latest edit: Mon Apr 10 2000
*/

#include "flow_functions.h"

static const short SWITCH_TRVL_ITERATES = 2;

////////////////////////////////////////////////////////////////////

// Called by 'Update_Transversal'.
// Changes the local Poincare' planes according to the dominating component
// of the vector field, and translates the parcel to the first new plane.
void Switch_Transversal(parcel &pcl, const short &trvl, const BOX &Switch_Box)
{
  short new_trvl = abs(trvl);
  short new_sign = Sign(trvl);
  INTERVAL time;
  BOX vf(DIM);
  BOX Euler_Box(DIM);
  BOX Temp_Switch_Box = Switch_Box;
  parcel result = pcl;         // Get the easy pieces of result

  result.trvl = new_trvl;
  result.sign = new_sign;

  for ( register short j = 1; j <= SWITCH_TRVL_ITERATES; j++ ) // No point in doing
    {                                                          // more than a few laps.
      Vf_Range(vf, Temp_Switch_Box);             // get the flow times
      if ( new_sign == + 1)
	time = Hull( 0, + Diam(pcl.box(new_trvl)) / vf(new_trvl) );
      else
	time = Hull( 0, - Diam(pcl.box(new_trvl)) / vf(new_trvl) );    
      
      for ( register short i = 1; i <= DIM; i++ ) // Compute the Euler step
	if ( i != new_trvl )
	  Euler_Box(i) = pcl.box(i) + time * vf(i);

      if ( new_sign == + 1)
	Euler_Box(new_trvl) = Hull(Sup(pcl.box(new_trvl)));
      else
	Euler_Box(new_trvl) = Hull(Inf(pcl.box(new_trvl)));

      if ( j != SWITCH_TRVL_ITERATES ) // We use this result to shrink
	{	                       // Temp_Switch_Box, and loop again.
	  Temp_Switch_Box = Euler_Box; // pcl.box(i) \subset Euler_Box(i), i != trvl.
	  Temp_Switch_Box(new_trvl) = Switch_Box(new_trvl); 
	}
    }
  result.box = Euler_Box;  // Save the new codim1

#ifdef COMPUTE_C1   // Compute the C1 info
  INTERVAL_MATRIX DPi(DIM, DIM);
  Get_DPi_Matrix(DPi, Temp_Switch_Box, new_trvl, time, result.box);

  if ( pcl.sign == 1 )
    {
      if ( new_sign == 1 )      
	Flow_Tangent_Vectors(result, pcl.trvl, new_trvl, DPi);
      else
	Flow_Tangent_Vectors(result, pcl.trvl, -new_trvl, DPi);
    }
  else
    {
      if ( new_sign == 1 )      
	Flow_Tangent_Vectors(result, -pcl.trvl, new_trvl, DPi);
      else
	Flow_Tangent_Vectors(result, -pcl.trvl, -new_trvl, DPi);
    }
#endif

  result.time = pcl.time + time;  // Add the flow time
  pcl = result;
}

////////////////////////////////////////////////////////////////////

// Called by 'Work_Through_PSin', 'Test_Switch_Transversal' or 'Flow'.
// Splits the 'BOX' of pcl into smaller pieces; a side is halved if
// it is longer than 'size'. The result is passed via 'PSin'.
void Single_Partition(const parcel &pcl, List<parcel> &PSin, const double &size)
{
  register int i, j;                    // Basic counters  
  int nr = 1;                           // Current number of BOXes
  double dx;                            // Diam/radius of pcl.box(i)
  BOX bx[POWER];                        // Storage for the BOXes

  bx[0] = pcl.box;  // Initialize 

  for(i = 1; i <= DIM; i++)
    if ( i != pcl.trvl )
      {
	dx = Sup(pcl.box(i)) - Inf(pcl.box(i));
	if ( dx > size )      // Check if splitting is necessary
	  {
	    dx /= 2.0;
	    for(j = 0; j < nr; j++)
	      { // Split in two
		bx[nr + j] = bx[j];
		bx[nr + j](i) = Hull(Inf(bx[nr + j](i)) + dx, Sup(bx[nr + j](i)));
		bx[j](i) = Hull(Inf(bx[j](i)), Inf(bx[nr + j](i)));
	      }
	    nr += nr; // There are twice as many boxes now.
	  }
      }

  parcel tmp_pcl = pcl;     // Temporary storage.
  for(j = 0; j < nr; j++)   // Now we return the partitioned boxes.
    {
      tmp_pcl.box = bx[j];
      PSin += tmp_pcl; 
    }
}

////////////////////////////////////////////////////////////////////

// Called by 'main' and 'Update_Transversal'.
// Returns the hull of all parcels in PSdone
void Get_Hull(parcel &hull_pcl, List<parcel> &pcl_List)
{
  hull_pcl = First(pcl_List);  // Initialize the local hull
  
  Next(pcl_List);
  while( !Finished(pcl_List) ) // Get the rectangular hull of the returned boxes,
    {                          // the flow times, and the tangent vector angles.
      hull_pcl = Hull(hull_pcl, Current(pcl_List));
      Next(pcl_List);
    }
}

////////////////////////////////////////////////////////////////////

// Called by 'Flow'.
// Computes DPi, the matrix of spatial derivatives of 
// the local Poincare map. Returns it by reference. By
// definition, DPi has DIM * (DIM - 1) non-zero elements.
// DPi(i, j) = DPhi(i, j) - vf(i) / vf(trvl) * DPhi(trvl, j).

// This version is better since it only evaluates Vf_Range over
// the codim-1 box Image.
void Get_DPi_Matrix(INTERVAL_MATRIX &DPi, const BOX &Outer_Box, 
		    const short &trvl, const INTERVAL &time, const BOX &Image)
{
  register short i, j;
  INTERVAL quotient;
  BOX vf(DIM); Vf_Range(vf, Image); // Vf_Range(vf, Outer_Box);
  INTERVAL_MATRIX DPhi;

  Get_DPhi_Matrix(DPhi, Outer_Box, time);
  DPi = DPhi;
  for (i = 1; i <= DIM; i++)
    if ( i != trvl )
      {
	quotient = vf(i) / vf(trvl);
	for (j = 1; j <= DIM; j++)  
	  DPi(i, j) -= quotient * DPhi(trvl, j);
      }
    else
      for (j = 1; j <= DIM; j++) // Not really neccessary, but
	DPi(i, j) = 0.0;         // it makes debugging clearer.
}

////////////////////////////////////////////////////////////////////

// Called by 'Flow'.
// We find the smallest time t such that
//   Inner_Box(i) + t * vf(Outer_Box)(i) \subset Outer_Box(i)
// is violated for i != trvl. Then we see if Inner_Box has completely 
// flowed out through the trvl side of Outer_Box during this time. If
// this is the case, we set time = [0, t]. Otherwise, let dx denote the
// shortest distance traveled in the trvl direction, and set
//   time = [0, dx/Inf(Abs(vf(Outer_Box)))(trvl)],
// where the trvl side of Outer_Box has been trimmed.
// Now set
//   Outer_Box = Inner_Box + time * vf(Outer_Box),
// and finally return (by reference) 
//   time = dx/Inf(Abs(vf(Outer_Box)))(trvl).
// If a close return becomes a true return, we modify 'pcl.message'.
void Get_Flow_Time(INTERVAL &time, parcel &pcl, const double &trvl_dist,
		   BOX &Outer_Box)
{
  register short i;
  double temp_time;
  double trvl_dx, dx;
  BOX    Inner_Box = pcl.box;              // Initialize.
  BOX    vf(DIM);  Vf_Range(vf, Outer_Box);
  double min_time  = Machine::PosInfinity; // A huge initial guess.

  for ( i = 1; i <= DIM; i++ )  // Find the first time any point of Inner_Box hits  
    if ( i != pcl.trvl )        // a non-transversal boundary of Outer_Box.
      {                                     // Uses the fact that the non-transversal
	dx = Inf(SubBounds(Sup(Outer_Box(i)), Sup(Inner_Box(i))));
	temp_time = dx / Abs(vf(i));        // coordinates of Outer_Box/Inner_Box are
	min_time = Min(min_time, temp_time);// symmetric, i.e., have the same centers.
      }

  if ( pcl.sign == 1 )                      // Now let's see how far this takes 
    trvl_dx = min_time * Inf(vf(pcl.trvl)); // us in the transversal direction.
  else
    trvl_dx = min_time * Inf(-vf(pcl.trvl));
  trvl_dx = Min(trvl_dx, trvl_dist);                 

  if ( trvl_dx < trvl_dist )            // If necessary, we shrink Outer_Box
    {                                   // in the transversal direction.
      double level = Inf(Inner_Box(pcl.trvl)); // Inf == Sup. 
      if ( pcl.sign == 1 )
	Outer_Box(pcl.trvl) = Hull(level, AddBounds(level, trvl_dx));
      else
	Outer_Box(pcl.trvl) = Hull(SubBounds(level, trvl_dx), level);
      Vf_Range(vf, Outer_Box);       // Recompute the vector field.
      time = trvl_dx / vf(pcl.trvl); // Compute the flow times required for 
      if ( pcl.sign == - 1 )         // all points in Inner_Box to flow through 
	time = - time;               // Outer_Box in the transveral direction.
      time = Hull(0.0, time);
    }
  else // if we flowed as far as we wanted to (i.e., out of Outer_Box(trvl))
    {    
      time = Hull(0.0, min_time);
      if ( pcl.message == CLOSE_STOP )
	pcl.message = STOP;
    }

  for ( i = 1; i <= DIM; i++ )
    if ( i != pcl.trvl )
      Outer_Box(i) = Inner_Box(i) + time * vf(i);  // Tighten Outer_Box.

  time = trvl_dx / Vf_Range(Outer_Box, pcl.trvl);  // Tighten time.
  if ( pcl.sign == - 1 )
    time = - time;
}

////////////////////////////////////////////////////////////////////

void Flow_Tangent_Vectors(parcel &pcl, const short &old_trvl, 
			  const short &new_trvl, const INTERVAL_MATRIX &DPi)
{
#ifdef COMPUTE_C1
  register short sin_i = 0, cos_i = 0, zip_i = 0;

  //cout << pcl << endl;
  switch( abs(old_trvl) )  // Get incoming coordinates.
    {
    case 1:
      zip_i = 1; sin_i = 2; cos_i = 3;
      break;
    case 2:
      zip_i = 2; sin_i = 1; cos_i = 3;
      break;
    case 3:
      zip_i = 3; sin_i = 2; cos_i = 1;
      break;
    default:  
      cout << endl << "abs(old_trvl) = " << abs(old_trvl) << endl;
      char *msg = "Error: 'Flow_Tangent_Vectors'. First switch failed!";
      throw Error_Handler(msg); 
      break;
    }

  INTERVAL old_theta = pcl.angles;  // Initialize the tangent vectors.
  BOX u_vect(DIM); 
  BOX v_vect(DIM); 

  u_vect(zip_i) = Hull(0.0);             v_vect(zip_i) = Hull(0.0); 
  u_vect(sin_i) = Sin(Inf(old_theta));   v_vect(sin_i) = Sin(Sup(old_theta)); 
  u_vect(cos_i) = Cos(Inf(old_theta));   v_vect(cos_i) = Cos(Sup(old_theta)); 

  BOX r_u_vect = DPi * u_vect;  // Compute the images
  BOX r_v_vect = DPi * v_vect;  // ot the vectors.

  if ( old_trvl != new_trvl )
    {      // Get outgoing coordinates.
      switch( abs(new_trvl) )
	{
	case 1:   // View in x(2)-x(3) coordinates.
	  zip_i = 1; sin_i = 2; cos_i = 3;
	  break;
	case 2:   // View in x(1)-x(2) coordinates.
	  zip_i = 2; sin_i = 1; cos_i = 3;
	  break;
	case 3:   // View in x(1)-x(2) coordinates.
	  zip_i = 3; sin_i = 2; cos_i = 1;
	  break;
	default:  
	  char *msg = "Error: 'Flow_Tangent_Vectors'. Second switch failed!";
	  throw Error_Handler(msg); 
	  break;
	}
    }

  INTERVAL u_angles;
  INTERVAL v_angles;

  if ( Inf(r_u_vect(cos_i)) > 0.0 && Inf(r_v_vect(cos_i)) > 0.0 )     
    {             
      u_angles = ArcTan(r_u_vect(sin_i) / r_u_vect(cos_i));
      v_angles = ArcTan(r_v_vect(sin_i) / r_v_vect(cos_i));
    }
  else if ( Sup(r_u_vect(cos_i)) < 0.0 && Sup(r_v_vect(cos_i)) < 0.0 )     
    {             
      r_u_vect = - r_u_vect;       // Reflect r_u and r_v into right-hand    
      r_v_vect = - r_v_vect;       // side of the unit circle.

      u_angles = ArcTan(r_u_vect(sin_i) / r_u_vect(cos_i));
      v_angles = ArcTan(r_v_vect(sin_i) / r_v_vect(cos_i));
    }
  else //  Subset(0.0, r_u_vect(cos_i)) or Subset(0.0, r_v_vect(cos_i)).
    {   
      // Debugging...
      cout << "Error: old_trvl = " << old_trvl << "; new_trvl = " << new_trvl << endl;
      cout << pcl << endl;
      cout << " u_vect = " << u_vect << endl;
      cout << " v_vect = " << v_vect << endl;
      cout << " DPi = " << endl << DPi << endl;
      cout << " r_u_vect = " << r_u_vect << endl;
      cout << " r_v_vect = " << r_v_vect << endl;
      // End debugging.
      char *msg = "Error: 'Flow_Tangent_Vectors'."
	" This case is not implemented yet.";
      throw Error_Handler(msg); 
    }

  INTERVAL new_diam = v_angles - u_angles;
  if ( Sup(new_diam) < 0.0 )        // Turn the vectors right.
    new_diam = - new_diam;
  else if ( Subset(0.0, new_diam) ) // Zero is an interior point.
    new_diam = Hull(0.0, Abs(new_diam));

  double   old_diam  = Diam(old_theta);
  INTERVAL factor    = Sqrt( (1 + Cos(new_diam)) / (1 + Cos(old_diam)) );
  INTERVAL prel_exp  = Hull(Norm(r_u_vect), Norm(r_v_vect));  

  pcl.expansion *= prel_exp * Hull(1.0, factor);
  pcl.angles = Hull(u_angles, v_angles);      // Update results.

#endif // COMPUTE_C1
}

////////////////////////////////////////////////////////////////////
