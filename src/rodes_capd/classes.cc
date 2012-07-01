/*   File: classes.cc (BIAS)
 
     Contains all class/struct definitions
     used in 'maintest'.

     Latest edit: Wed Feb 23 2000
*/

#include "classes.h"

////////////////////////////////////////////////////////////////////

// Initializes an 'INTERVAL' out of 
// the two 'doubles' lo and hi. Used 
// for external input of 'lo' and 'hi'.
// **jjb -- Not sure why we don't just use INTERVAL here **
INTERVAL Init_Interval(const double &lo, const double &hi) 
{ return Succ(Hull(lo, hi)); }

// Returns an INTERVAL containing        Called mid in BIAS, and
// the center of an INTERVAL             there only returns a REAL
// jjb - since CAPD mid returns an interval, we change here to just
// using mid.
INTERVAL Center(const INTERVAL &iv)
{ return mid(iv); }

// Returns an INTERVAL containing
// the radius of an INTERVAL
// jjb - since CAPD diam returns an interval, we change here to just
// using diam.
INTERVAL Radius(const INTERVAL &iv)
{ return diam(iv) / 2.0; }

// Returns an INTERVAL containing
// the symmetric radius of an INTERVAL
INTERVAL Symm_Radius(const INTERVAL &iv)
{  
    INTERVAL _r = diam( iv ) / 2.0;
    double r = _r . leftBound ( );
    return INTERVAL( -r, r ); 
}

// Returns (by reference) intervals containing both
// The center and the symmetric radius of 'iv'.
void Mid_And_SymRad(INTERVAL &center, INTERVAL &symmrad, const INTERVAL &iv)
{
    symmrad = Symm_Radius( iv ); //SymHull(diam(iv) / 2.0);
    center = mid ( iv );
}

// Rescales the INTERVAL by a non-neg. factor, i.e.,
// radius([iv_out]) = factor*radius([iv_in]);
// center([iv_out]) = center([iv_in]);
INTERVAL Rescale(const INTERVAL &iv, const double &factor)
{
  return INTERVAL (Center(iv) + factor * Symm_Radius(iv));
}

// Checks if the 'double' dbl is
// contained in the 'INTERVAL' iv
// The <= seems to work the same in CAPD
bool Subset(const double &dbl, const INTERVAL &iv) 
{ return ( dbl <= iv ); }

// Checks if the 'INTERVAL' iv1 is
// a subset of the 'INTERVAL' iv2
bool Subset(const INTERVAL &iv1, const INTERVAL &iv2)  
{ return ( iv1.subset ( iv2 ) ); }

// Returns the sign of all numbers in iv
// If iv contain zero, a warning message 
// is printed, and zero is returned.
// jjb -- replaced Sup, Inf with CAPD equivs
int Sign(const INTERVAL &iv)
{
  if ( rightBound ( iv ) < 0 )
    return - 1;
  if ( leftBound ( iv ) > 0 )
    return + 1;

  throw Error_Handler("Sign: The interval contains zero!");
  return 0;
}

// Prints an 'INTERVAL' to the standard output
void Show_Interval(const INTERVAL &iv)
{
  cout.precision( NUMBER_OF_DIGITS );
  cout.setf( ios::showpos );
  cout.setf( ios::scientific );
  cout << "Showing interval (via fcn):" << endl
       << "[" << Inf(iv) << ", " << Sup(iv) << "];"
       << " dx = " << diam(iv) << endl;
  cout.unsetf( ios::showpos );
  cout.unsetf( ios::scientific );
}

////////////////////////////////////////////////////////////////////

// Returns a BOX containing
// the center of a BOX
// jjb -- grab the mid interval in each dimension
BOX Center(const BOX &bx)
{
  BOX temp(DIM);

  for (register short i = 1; i <= DIM; i++)
    temp(i) = mid ( bx (i ) ); 

  return temp;
}

// Returns a BOX containing
// the radius of a BOX
BOX Radius(const BOX &bx)
{
  BOX temp(DIM);

  for (register short i = 1; i <= DIM; i++) 
    temp(i) = diam ( bx ( i ) ) / 2.0;  

  return temp;
}

// Returns a BOX containing
// the symmetric radius of a BOX
BOX Symm_Radius(const BOX &bx)
{
  BOX temp(DIM);

  for (register short i = 1; i <= DIM; i++)
    temp(i) = Sym_Radius ( bx(i) ); // jjb -- SymHull(diam(bx(i)) / 2.0);

  return temp;
}

// Returns (by reference) boxes containing both
// The center and the symmetric radius of 'bx'.
void Mid_And_SymRad(BOX &center, BOX &symmrad, const BOX &bx)
{
  for (register short i = 1; i <= DIM; i++)
    {
      symmrad(i) = Sym_Radius ( bx ( i ) );// jjb -- SymHull(diam(bx(i)) / 2.0);
      center(i) = mid ( bx ( i ) );
    } 
}

// Rescales the BOX by a non-neg. factor, i.e.,
// radius([bx_out]) = factor*radius([bx_in]);
// center([bx_out]) = center([bx_in]);
BOX Rescale(const BOX &bx, const double &factor)
{
  BOX result(DIM);

  for (register short i = 1; i <= DIM; i++)
    result(i) = Center(bx(i)) + factor * Symm_Radius(bx(i));
      
  return result;
}

// Checks if the 'BOX' box1 is
// a subset of the 'BOX' box2
bool Subset(const BOX &box1, const BOX &box2)  
{ return ( box1.subset( box2 ) ); } //(box1 <= box2); }

////////////////////////////////////////////////////////////////////

parcel Hull(const parcel &pcl_1, const parcel &pcl_2)
{
  parcel result = pcl_1;

  result.box  = Hull(pcl_1.box,  pcl_2.box);
  result.time = Hull(pcl_1.time, pcl_2.time);
#ifdef COMPUTE_C1
  result.angles = Hull(pcl_1.angles, pcl_2.angles);
  result.expansion = Hull(pcl_1.expansion, pcl_2.expansion);
#endif  
  return result;
}

ostream & operator << (ostream &out, const parcel &pcl)
{   
  out.precision(NUMBER_OF_DIGITS);
  out.setf(ios::showpos);
  out.setf(ios::scientific);
  out << "Showing parcel (via class):" << endl
      << "  transversal = " << pcl.trvl
      << ";  sign = " << pcl.sign
      << ";  message = " << pcl.message << endl;
  out << "  time = " << pcl.time << "; dt = " << diam(pcl.time) << endl;
  out << "  box  = ";
  for ( int i = 1; i <= Dimension(pcl.box); i++ ) 
    {
      if ( i != 1 )
	out << "         ";
      out << pcl.box(i) << "; dx = " << diam(pcl.box(i)) << endl;
    }
#ifdef COMPUTE_C1
  INTERVAL angles = RAD_TO_DEG * pcl.angles;
  out << "  ang  = " << angles << "; da = " << diam(angles) << " (deg)\n"
      << "  exp  = " << pcl.expansion << "; de = " << diam(pcl.expansion)<< endl;
#endif
  out.unsetf(ios::showpos);
  out.unsetf(ios::scientific);

  return out; 
}

////////////////////////////////////////////////////////////////////
