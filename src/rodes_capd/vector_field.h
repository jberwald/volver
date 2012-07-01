/*   File: vector_field.h (BIAS)
 
     Contains the vector field information

     Latest edit: Thu Feb 10 2000
*/

#ifndef VECTOR_FIELD_H
#define VECTOR_FIELD_H

#include "classes.h"
#include "error_handler.h"

////////////////////////////////////////////////////////////////////

// Classical constants for the 'lorenz' system
const double S = 10.0;
const double R = 28.0;
const double B = 8./3;

// Parameters in the Jordan normal form
const double TEMP = sqrt((S + 1) * (S + 1) + 4 * S * (R - 1));
const double K1 = S / TEMP;
const double K2 = (S - 1 + TEMP) / (2 * S);
const double K3 = (S - 1 - TEMP) / (2 * S);

// Eigenvalues at the origin
const double E1 = (- (S + 1) + TEMP) / 2;
const double E2 = (- (S + 1) - TEMP) / 2;
const double E3 = - B;

////////////////////////////////////////////////////////////////////
// Rigorous Constants from here.

// Intervals constructed only once
const INTERVAL S_IV = Init_Interval(S, S);
const INTERVAL R_IV = Init_Interval(R, R);
const INTERVAL B_IV = Init_Interval(B, B);

// Parameters in the Jordan normal form
const INTERVAL TEMP_IV = Sqrt((S_IV + 1)*(S_IV + 1) + 4*S_IV*(R_IV - 1));
const INTERVAL K1_IV = S_IV/TEMP_IV;
const INTERVAL K2_IV = (S_IV - 1 + TEMP_IV)/(2*S_IV);
const INTERVAL K3_IV = (S_IV - 1 - TEMP_IV)/(2*S_IV);

// Eigenvalues at the origin
const INTERVAL E1_IV = (-(S_IV + 1) + TEMP_IV)/2;
const INTERVAL E2_IV = (-(S_IV + 1) - TEMP_IV)/2;
const INTERVAL E3_IV = -B_IV;

// Shortcuts
const INTERVAL K2_PLUS_K3_IV = K2_IV + K3_IV; 
const INTERVAL TWO_K2_IV = 2 * K2_IV; 
const INTERVAL TWO_K3_IV = 2 * K3_IV; 
const INTERVAL ZERO_IV(0.0);
const INTERVAL ONE_IV(1.0);

////////////////////////////////////////////////////////////////////

void     NR_Vf_Range(VECTOR &, const VECTOR &);

BOX      Vf_Range (const BOX &); 
void     Vf_Range (BOX &, const BOX &);

INTERVAL Vf_Range (const BOX &, const short &);
void     Vf_Range (INTERVAL &, const BOX &, const short &);

void     DVf_Range(INTERVAL_MATRIX &, const BOX &);

INTERVAL DVf_Range(const BOX &, const short &, const short &);
void     DVf_Range(INTERVAL &, const BOX &, const short &, const short &);

////////////////////////////////////////////////////////////////////

#endif // VECTOR_FIELD_H
