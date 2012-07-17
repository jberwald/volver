/*   File: return_map.h (BIAS)
 
     Contains all high level functions
     needed to compute a return.

     Latest edit: Mon Apr 10 2000
*/

#ifndef RETURN_MAP_H
#define RETURN_MAP_H

#include "classes.h"
#include "error_handler.h"
#include "fixed_point.h"
#include "flow_functions.h"
#include "list.h"
#include "low_functions.h"

////////////////////////////////////////////////////////////////////

// Integrator parameters
const double SCALE_FACTOR  = 1.1;          // Used in "Flow' for the coarse enclosure.

// Stopping parameters
const double STOP_DIST_LEVEL  =  27.0;
const short  STOP_SIGN        = -1;
const short  STOP_TRANSVERSAL =  3;

////////////////////////////////////////////////////////////////////

void Compute_the_return   (const parcel &, List<parcel> &);

void Local_Flow_The_Parcel(const parcel &, List<parcel> &,
			   const stop_parameters &, const double &);

void Multiple_Partition   (const parcel &, List<parcel> &, const double &);

////////////////////////////////////////////////////////////////////

#endif // RETURN_MAP_H
