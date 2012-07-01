/*   File: low_functions.h (BIAS)
 
     Contains all low level functions
     needed to compute a return.

     Latest edit: Mon Apr 10 2000
*/

#ifndef LOW_FUNCTIONS_H
#define LOW_FUNCTIONS_H

#include "classes.h"
#include "error_handler.h"
#include "vector_field.h"

////////////////////////////////////////////////////////////////////

int  Sign                  (const double &);

int  Int_Power             (const int &, const int &);

// POWER declared as a fixed global variable to save time.
const int POWER = Int_Power(2, DIM - 1);

void Some_May_Vanish       (BOX &, const INTERVAL_MATRIX &, const parcel &,
			    const BOX &, const BOX &, const double &);

void None_May_Vanish       (BOX &, const parcel &, const BOX &,
			    const BOX &, const double &);

void Flow_By_Corner_Method (BOX &, const INTERVAL_MATRIX &,
			    const parcel &, const BOX &);

void Get_DPhi_Matrix       (INTERVAL_MATRIX &, const BOX &, const INTERVAL &);

////////////////////////////////////////////////////////////////////

#endif // LOW_FUNCTIONS_H
