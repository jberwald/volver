/*   File: convert.h 
 
     A set of functions used for converting
     between data_line/parcel based classes.

     Latest edit: Mon Mar 27 2000
*/

#ifndef CONVERT_H
#define CONVERT_H

#include <math.h>

#include "2d_classes.h"
#include "classes.h"
#include "flow_functions.h"
#include "list.h"

////////////////////////////////////////////////////////////////////

void iterate_to_parcel   (const iterate &, parcel    &);
void rect_to_it_List     (const BOX     &, const int &, List<iterate> &);
void pcl_List_to_it_List (List<parcel>  &, const int &, List<iterate> &);
void New_Get_Image_Hull  (iterate       &, List<iterate> &); // Both u and v.

////////////////////////////////////////////////////////////////////

#endif // CONVERT_H
