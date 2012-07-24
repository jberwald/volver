/*   File: flow_functions.h (BIAS)
 
     Contains all medium level functions
     needed to compute a return.

     Latest edit: Mon Apr 10 2000
*/

#ifndef FLOW_FUNCTIONS_H
#define FLOW_FUNCTIONS_H

#include "classes.h"
#include "error_handler.h"
#include "list.h"
#include "vector_field.h"
#include "low_functions.h"

////////////////////////////////////////////////////////////////////

enum MESSAGE_CODE { OK, TOO_LARGE, CLOSE_STOP, STOP, PARTITION };

////////////////////////////////////////////////////////////////////

void Switch_Transversal   (parcel &, const short &, const BOX &);

void Single_Partition     (const parcel &, List<parcel> &, const double &);

void Get_Hull             (parcel &, List<parcel> &);

void Get_DPi_Matrix       (IMatrix &, const BOX &, const short &, 
			   const interval &, const BOX &);

void Get_Flow_Time        (interval &, parcel &, const double &, BOX &);

void Flow_Tangent_Vectors (parcel &, const short &, const short &, 
			   const IMatrix &);

////////////////////////////////////////////////////////////////////

#endif // FLOW_FUNCTIONS_H
