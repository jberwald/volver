/*   File: fixed_point.h (BIAS)
 
     Contains all functions needed to deal
     with parcels entering the cube containing
     the fixed point at the origin.

     Latest edit: Fri Feb 25 2000
*/

#ifndef FIXED_POINT_H
#define FIXED_POINT_H

#include <cmath>

#include "classes.h"
#include "flow_functions.h"
#include "list.h"
#include "return_map.h"
#include "vector_field.h"

////////////////////////////////////////////////////////////////////

bool Cube_Entry (const parcel &);

void Cube_Exit  (const parcel &, List<parcel> &, const double &);


////////////////////////////////////////////////////////////////////

#endif // FIXED_POINT_H
