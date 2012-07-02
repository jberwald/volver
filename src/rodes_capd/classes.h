/*   File: classes.h (BIAS) --> CONVERTED TO CAPD
 
     Contains all class/struct definitions
     used in 'maintest'.

     Latest edit: Wed Feb 23 2000
*/

#ifndef CLASSES_H
#define CLASSES_H

#include <iomanip>
#include <iostream>
#include <cmath>

////////////////////////////////////////////////////////////////////

/* #include "Constants.h"                 // PROFIL/BIAS header */
/* #include "Error.h"                     // PROFIL/BIAS header */
/* #include "Functions.h"                 // PROFIL/BIAS header */
/* #include "Interval.h"                  // PROFIL/BIAS header */
/* #include "IntervalMatrix.h"            // PROFIL/BIAS header */
/* #include "IntervalVector.h"            // PROFIL/BIAS header */
/* #include "Matrix.h"                    // PROFIL/BIAS header */
/* #include "Vector.h"                    // PROFIL/BIAS header */

//#include "capd/capdlib.h"
//#include "capd/intervals/Interval.h"
#include "capd/rounding/DoubleRounding.h"
#include "capd/intervals/DoubleInterval.h"
#include "capd/intervals/IntervalError.h"
#include "capd/vectalg/vectalgLib.h"


////////////////////////////////////////////////////////////////////

#include "error_handler.h"

////////////////////////////////////////////////////////////////////

const int    NUMBER_OF_DIGITS = 6;
const short  SYSDIM           = 3;

#undef  COMPUTE_C1 // Comment out next line for topological mode.
//#define COMPUTE_C1 

// Define CAPD interval and interval vector of fixed dimension
//typedef capd::intervals::Interval< double > DInterval;
//typedef capd::intervals::Interval< double > INTERVAL;

typedef capd::vectalg::Vector< interval, SYSDIM > IVector;
////////////////////////////////////////////////////////////////////

const interval PI = interval::pi(); // jjb -- = Succ(Hull(Constant::Pi));
const interval DEG_TO_RAD = PI / 180.0;
const interval RAD_TO_DEG = 180.0 / PI;

interval Init_Interval  (const double   &, const double   &);
interval Center         (const interval &);
interval Radius         (const interval &);
interval Symm_Radius    (const interval &); 
void     Mid_And_SymRad (      interval &,       interval &, const interval &);
interval Rescale        (const interval &, const double   &); 
bool     Subset         (const double   &, const interval &);
bool     Subset         (const interval &, const interval &);
int      Sign           (const interval &);
void     Show_Interval  (const interval &);

////////////////////////////////////////////////////////////////////

#define BOX IVector    // Shorthand

BOX  Center         (const BOX &);
BOX  Radius         (const BOX &);
BOX  Symm_Radius    (const BOX &);
void Mid_And_SymRad (      BOX &,       BOX    &, const BOX &);
BOX  Rescale        (const BOX &, const double &);
bool Subset         (const BOX &, const BOX    &);

////////////////////////////////////////////////////////////////////

// Redefine some of the global CAPD functions to align with common functions in RODES
double Sup ( const interval & );
double Inf ( const interval & );

////////////////////////////////////////////////////////////////////
class parcel
{
public:
  BOX box;                   // The coordinates of all the variables 
#ifdef COMPUTE_C1
  interval angles;
  interval expansion;
#endif
  short trvl;                // the transversal variable: 1,...,DIM
  short sign;                // the direction of the flow: - 1 or + 1
  interval time;             // The "time" variable   
  short message;             // Any message that needs to be passed on 
  friend parcel Hull(const parcel &, const parcel &);
  friend ostream & operator << (ostream &, const parcel &);
};

////////////////////////////////////////////////////////////////////

typedef struct
{
  double R; // Red
  double G; // Green 
  double B; // Blue
} colour;

const colour BLACK     = {0, 0, 0};
const colour RED       = {1, 0, 0};
const colour GREEN     = {0, 1, 0};
const colour BLUE      = {0, 0, 1};
const colour PURPLE    = {1, 0, 1};
const colour YELLOW    = {1, 1, 0};
const colour TURQUOISE = {0, 1, 1};

////////////////////////////////////////////////////////////////////

typedef struct          // The structure defining
{                       // local/global stopping conds
  double level; 
  short  trvl; 
  short  sign; 
  double max_d_step; 
} stop_parameters;

////////////////////////////////////////////////////////////////////

#endif // CLASSES_H
