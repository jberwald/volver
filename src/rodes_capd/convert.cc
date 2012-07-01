/*   File: convert.cc 
 
     Converts data_line structures to
     parcel structures and vice versa.
     Also filters parcel returns into
     data_line grids.

     Latest edit: Thu Apr 20 2000
*/

#include "convert.h"

////////////////////////////////////////////////////////////////////
// Called by: 'work_on_grid'.
// Calls to : 'none'
void iterate_to_parcel(const iterate &it, parcel &pcl)
{
  BOX rect = ndl_to_box(it.ndl);

  pcl.box(1) = rect(1);
  pcl.box(2) = rect(2);
  pcl.box(3) = INTERVAL(27.0, 27.0);
  pcl.time   = INTERVAL(0.0, 0.0);
  pcl.trvl = 3; pcl.sign = -1; pcl.message = 0;
#ifdef COMPUTE_C1
  pcl.angles    = it.ndl.ang;         // pcl.angles    = it.dl.ang;
  pcl.expansion = INTERVAL(1.0, 1.0); // pcl.expansion = it.dl.exp;
#endif
}

////////////////////////////////////////////////////////////////////
// Called by: 'get_the_flags', 'pcl_List_to_it_List'.
// Calls to : 'none'
void rect_to_it_List(const BOX &box, const int &power, List<iterate> &it_List)
{
  int i, j;
  int inf[2], sup[2];
  iterate it;
  BOX rect = pow(2, power) * box; 

  for ( i = 0; i < 2; i++ )
    {
      inf[i] = (int) ceil (Inf(rect(i + 1)));
      sup[i] = (int) floor(Sup(rect(i + 1)));
     
      // Make the coordinates odd.
      if ( inf[i] % 2 == 0 )
	inf[i]--;
      if ( sup[i] % 2 == 0 )
	sup[i]++;
    }

  it.ndl.grd.P  = power;
  it.ndl.c_stat = NOT_DONE;
  it.ndl.h_stat = NOT_HIT;
  it.inf_grd = NULL_GRID;
  it.sup_grd = NULL_GRID;

  for ( i = inf[0]; i <= sup[0]; i += 2 )
    for ( j = inf[1]; j <= sup[1]; j += 2 )
      {
	it.ndl.grd.u = i;
	it.ndl.grd.v = j;
	it_List += it;
      }
}

////////////////////////////////////////////////////////////////////
// Called by: 'work_on_grid'
// Calls to : 'rect_to_it_List'
//
// The Parcel_List contains all Q_{i,j}'s. The iterate version
// (with pre_exp) is returned via Iterate_List.
void pcl_List_to_it_List(List<parcel> &Parcel_List, const int &power,
			 List<iterate> &Iterate_List)
{
  parcel pcl;
  BOX rect(2);
  iterate it, cmp_it;
  List<iterate> it_List, Redundant_List;  
  List<parcel>  Dummy_pcl_List;         // Just to shut the compiler up!
  
  First(Parcel_List);
  while( !Finished(Parcel_List) )
    { // Loop through all parcels.
      pcl = Current(Parcel_List);
      rect(1) = pcl.box(1);
      rect(2) = pcl.box(2);
      rect_to_it_List(rect, power, it_List); // Gives it.ndl.c_stat == NOT_DONE, and
      while ( !IsEmpty(it_List) )            // it.inf_grd == it.sup_grd == NULL_GRID.
	{
	  it = First(it_List);
#ifdef COMPUTE_C1
	  it.ndl.ang = pcl.angles;
	  it.ndl.pre_exp = Inf(pcl.expansion); // Inf(E_{i,j}).
	  it.ndl.min_exp = LARGE_NUMBER; // The iterates are not known to be computed yet.
#endif
	  it.ndl.h_stat = HIT;  // The iterates are known to be hit.
	  Redundant_List += it; 
	  RemoveCurrent(it_List);
	}
      Next(Parcel_List);
    }

  // Loop through Redundant_List, and remove redundancies.
  while( !IsEmpty(Redundant_List) )
    {  // Get the first it and remove it.
      it = First(Redundant_List);
      RemoveCurrent(Redundant_List);

// ADDED Jul 10, 2000 TO REMOVE ANNOYING OUTPUT.

      if( IsEmpty(Redundant_List) )
        {
          Iterate_List += it;        
          break;
        }

// END ADDED.
	
      First(Redundant_List);
      while( !IsEmpty(Redundant_List) && !Finished(Redundant_List) )
	{ 
	  cmp_it = Current(Redundant_List); // Get the next it for comparison.
	  if ( it.ndl.grd == cmp_it.ndl.grd ) // If they represent the same grid...
	    {
#ifdef COMPUTE_C1                             // ... we take their hull...
	      it.ndl.ang = Hull(it.ndl.ang, cmp_it.ndl.ang);
	      it.ndl.pre_exp = Min(it.ndl.pre_exp, cmp_it.ndl.pre_exp); 
	      it.ndl.min_exp = Min(it.ndl.min_exp, cmp_it.ndl.min_exp); 
#endif
	      RemoveCurrent(Redundant_List);  // ... and remove the redundant it.
	    }
	  else
	    Next(Redundant_List);
	}
      Iterate_List += it; // Save the hull of all iterates that represent the same grid.
    }
}

////////////////////////////////////////////////////////////////////

// Called by: 'work_on_grid'.
// Calls to : 'grid_to_box'.
// Modifies the iterate's members 'inf_grd' and 'sup_grd'.
void New_Get_Image_Hull(iterate &it, List<iterate> &Iterate_List)
{
  grid tmp_grd = First(Iterate_List).ndl.grd;
  int min_u = tmp_grd.u;
  int max_u = tmp_grd.u;
  int min_v = tmp_grd.v;
  int max_v = tmp_grd.v;
  
  Next(Iterate_List);
  while( !Finished(Iterate_List) )
    {
      tmp_grd = Current(Iterate_List).ndl.grd;
      if ( min_u > tmp_grd.u )
	min_u = tmp_grd.u;
      if ( max_u < tmp_grd.u )
	max_u = tmp_grd.u;
      if ( min_v > tmp_grd.v )
	min_v = tmp_grd.v;
      if ( max_v < tmp_grd.v )
	max_v = tmp_grd.v;
      Next(Iterate_List);
    }
  it.inf_grd.u = min_u;
  it.inf_grd.v = min_v;
  it.sup_grd.u = max_u;
  it.sup_grd.v = max_v;
}

////////////////////////////////////////////////////////////////////
