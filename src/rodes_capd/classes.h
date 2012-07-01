/*   File: classes.h (BIAS) --> CONVERTED TO CAPD
 
     Contains all class/struct definitions
     used in 'maintest'.

     Latest edit: Wed Feb 23 2000
*/

#ifndef CLASSES_H
#define CLASSES_H

#include <iomanip.h>
#include <iostream.h>
#include <math.h>

////////////////////////////////////////////////////////////////////

/* #include "Constants.h"                 // PROFIL/BIAS header */
/* #include "Error.h"                     // PROFIL/BIAS header */
/* #include "Functions.h"                 // PROFIL/BIAS header */
/* #include "Interval.h"                  // PROFIL/BIAS header */
/* #include "IntervalMatrix.h"            // PROFIL/BIAS header */
/* #include "IntervalVector.h"            // PROFIL/BIAS header */
/* #include "Matrix.h"                    // PROFIL/BIAS header */
/* #include "Vector.h"                    // PROFIL/BIAS header */

#include "capd/capdAlglib.h"
#include "capd/capdlib.h"
//#include "capd/intervals/Interval.h"

////////////////////////////////////////////////////////////////////

//#include "error_handler.h"

////////////////////////////////////////////////////////////////////

const int    NUMBER_OF_DIGITS = 6;
const short  DIM              = 3;

#undef  COMPUTE_C1 // Comment out next line for topological mode.
//#define COMPUTE_C1 

typedef PointBase< capd::intervals::Interval< double > > INTERVAL;
typedef capd::vectalg::Vector<DInterval,0> INTERVAL_VECTOR;
////////////////////////////////////////////////////////////////////

const INTERVAL INTERVAL::pi = Succ(Hull(Constant::Pi));
const INTERVAL DEG_TO_RAD = INTERVAL::pi / 180.0;
const INTERVAL RAD_TO_DEG = 180.0 / PI;

INTERVAL Init_Interval  (const double   &, const double   &);
INTERVAL Center         (const INTERVAL &);
INTERVAL Radius         (const INTERVAL &);
INTERVAL Symm_Radius    (const INTERVAL &); 
void     Mid_And_SymRad (      INTERVAL &,       INTERVAL &, const INTERVAL &);
INTERVAL Rescale        (const INTERVAL &, const double   &); 
bool     Subset         (const double   &, const INTERVAL &);
bool     Subset         (const INTERVAL &, const INTERVAL &);
int      Sign           (const INTERVAL &);
void     Show_Interval  (const INTERVAL &);

////////////////////////////////////////////////////////////////////

#define BOX INTERVAL_VECTOR    // Shorthand

BOX  Center         (const BOX &);
BOX  Radius         (const BOX &);
BOX  Symm_Radius    (const BOX &);
void Mid_And_SymRad (      BOX &,       BOX    &, const BOX &);
BOX  Rescale        (const BOX &, const double &);
bool Subset         (const BOX &, const BOX    &);

////////////////////////////////////////////////////////////////////

class parcel
{
public:
  BOX box;                   // The coordinates of all the variables 
#ifdef COMPUTE_C1
  INTERVAL angles;
  INTERVAL expansion;
#endif
  short trvl;                // the transversal variable: 1,...,DIM
  short sign;                // the direction of the flow: - 1 or + 1
  INTERVAL time;             // The "time" variable   
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
