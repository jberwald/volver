/*   File: classes.cc (BIAS)
 
     Contains all class/struct definitions
     used in 'maintest'.

     Latest edit: Wed Feb 23 2000
*/

#include "classes.h"

////////////////////////////////////////////////////////////////////

// Initializes an 'interval' out of 
// the two 'doubles' lo and hi. Used 
// for external input of 'lo' and 'hi'.
// **jjb -- Not sure why we don't just use interval here **
interval Init_Interval(const double &lo, const double &hi) 
{ return interval( lo, hi ); } //Succ(Hull(lo, hi)); }

// Returns an interval containing        Called mid in BIAS, and
// the center of an interval             there only returns a REAL
// jjb - since CAPD mid returns an interval, we change here to just
// using mid.
interval Center(const interval &iv)
{ return mid(iv); }

// Returns an interval containing
// the radius of an interval
// jjb - since CAPD diam returns an interval, we change here to just
// using diam.
interval Radius(const interval &iv)
{ return diam(iv) / 2.0; }


// Returns an interval containing
// the symmetric radius of a real number
interval Symm_Radius(const double &dbl)
{  
    return interval( -dbl, dbl ); 
}

// Returns an interval containing
// the symmetric radius of an interval
interval Symm_Radius(const interval &iv)
{  
    interval _r = diam( iv ) / 2.0;
    double r = _r . leftBound ( );
    return interval( -r, r ); 
}

// Returns (by reference) intervals containing both
// The center and the symmetric radius of 'iv'.
void Mid_And_SymRad(interval &center, interval &symmrad, const interval &iv)
{
    symmrad = Symm_Radius( iv ); //SymHull(diam(iv) / 2.0);
    center = mid ( iv );
}

// Rescales the interval by a non-neg. factor, i.e.,
// radius([iv_out]) = factor*radius([iv_in]);
// center([iv_out]) = center([iv_in]);
interval Rescale(const interval &iv, const double &factor)
{
  return interval (Center(iv) + factor * Symm_Radius(iv));
}

// Checks if the 'double' dbl is
// contained in the 'interval' iv
// jb -- Using iv.contains(..)
bool Subset(const double &dbl, const interval &iv) 
{ return iv.contains ( dbl ); }

// Checks if the 'interval' iv1 is
// a subset of the 'interval' iv2
bool Subset(const interval &iv1, const interval &iv2)  
{ return ( iv1.subset ( iv2 ) ); }

// Returns the sign of all numbers in iv
// If iv contain zero, a warning message 
// is printed, and zero is returned.
// jjb -- replaced Sup, Inf with CAPD equivs
int Sign(const interval &iv)
{
  if ( rightBound ( iv ) < 0 )
    return - 1;
  if ( leftBound ( iv ) > 0 )
    return + 1;

  throw Error_Handler("Sign: The interval contains zero!");
  return 0;
}

// Prints an 'interval' to the standard output
void Show_Interval(const interval &iv)
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

// Returns a BOX containing the center of a BOX
// jjb -- grab the mid
// interval in each dimension.(i.e. an IVector or 2 or 3 dimensions
// with singleton intervals as coordinates.
BOX Center(const BOX &bx)
{
  BOX temp(SYSDIM);

  for (register short i = 1; i <= SYSDIM; i++)
    temp(i) = mid ( bx (i ) ); 

  return temp;
}

// Returns a BOX containing
// the radius of a BOX
BOX Radius(const BOX &bx)
{
  BOX temp(SYSDIM);

  for (register short i = 1; i <= SYSDIM; i++) 
    temp(i) = diam ( bx ( i ) ) / 2.0;  

  return temp;
}

// Returns a BOX containing
// the symmetric radius of a BOX
BOX Symm_Radius(const BOX &bx)
{
  BOX temp(SYSDIM);

  for (register short i = 1; i <= SYSDIM; i++)
    temp(i) = Symm_Radius ( bx(i) ); // jjb -- SymHull(diam(bx(i)) / 2.0);

  return temp;
}

// Returns (by reference) boxes containing both
// The center and the symmetric radius of 'bx'.
void Mid_And_SymRad(BOX &center, BOX &symmrad, const BOX &bx)
{
  for (register short i = 1; i <= SYSDIM; i++)
    {
      symmrad(i) = Symm_Radius ( bx ( i ) );// jjb -- SymHull(diam(bx(i)) / 2.0);
      center(i) = mid ( bx ( i ) );
    } 
}

// Rescales the BOX by a non-neg. factor, i.e.,
// radius([bx_out]) = factor*radius([bx_in]);
// center([bx_out]) = center([bx_in]);
BOX Rescale(const BOX &bx, const double &factor)
{
  BOX result(SYSDIM);

  for (register short i = 1; i <= SYSDIM; i++)
    result(i) = Center(bx(i)) + factor * Symm_Radius(bx(i));
      
  return result;
}

// Checks if the 'BOX' box1 is
// a subset of the 'BOX' box2
bool Subset(const BOX &box1, const BOX &box2)  
{ return ( subset( box1, box2 ) ); } //(box1 <= box2); }

////////////////////////////////////////////////////////////////////

// jjb -- Hull used to overload PROFIL's Hull(interval,interval) by
// calling PROFIL's Hull where "intervalHull" is below. We modify it
// here to perform "hull of parcel" calculations.
parcel Hull(const parcel &pcl_1, const parcel &pcl_2)
{
  parcel result = pcl_1;
    
  // jjb -- CAPD global function 'intervalHull' replaces PROFIL's
  // 'Hull(interval,interval)'
  result.box  = intervalHull(pcl_1.box,  pcl_2.box);
  result.time = intervalHull(pcl_1.time, pcl_2.time);
#ifdef COMPUTE_C1
  result.angles = intervalHull(pcl_1.angles, pcl_2.angles);
  result.expansion = intervalHull(pcl_1.expansion, pcl_2.expansion);
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
    for ( int i = 1; i <= pcl.box.dimension( ); i++ ) 
    {
      if ( i != 1 )
	out << "         ";
      out << pcl.box(i) << "; dx = " << diam(pcl.box(i)) << endl;
    }
#ifdef COMPUTE_C1
  interval angles = RAD_TO_DEG * pcl.angles;
  out << "  ang  = " << angles << "; da = " << diam(angles) << " (deg)\n"
      << "  exp  = " << pcl.expansion << "; de = " << diam(pcl.expansion)<< endl;
#endif
  out.unsetf(ios::showpos);
  out.unsetf(ios::scientific);

  return out; 
}

////////////////////////////////////////////////////////////////////
// We recast common RODES functions to use CAPD global functions

double Sup ( const interval &iv )
{ 
    return iv . rightBound ();
}

double Inf ( const interval &iv )
{
    return iv . leftBound ();
}

// Returns an IVector containing the lower bounds of the IVector vec. as singletons. 
DVector Inf ( const IVector &vec )
{
  // fill infVec with the infima of vec. 
    DVector infVec ( vec.dimension() );
    for ( register short i = 0; i < vec.dimension(); i++ )
      {
	infVec [ i ] = vec [ i ] . leftBound(); 
      }
    return infVec;
}

// Returns an IVector containing the lower bounds of the IVector vec. as singletons. 
DVector Sup ( const IVector &vec )
{
  // fill infVec with the infima of vec. 
    DVector supVec ( vec.dimension() );
    for ( register short i = 0; i < vec.dimension(); i++ )
      {
	supVec [ i ] = vec [ i ] . rightBound(); 
      }
    return supVec;
}

interval AddBounds ( const double &d1, const double &d2 )
{
    interval iv1 ( d1 ), iv2 ( d2 );
    return iv1 + iv2;
}

interval SubBounds ( const double &d1, const double &d2 )
{
    interval iv1 ( d1 ), iv2 ( d2 );
    return iv2 - iv1;
}

IVector SubBounds ( const DVector &iv1, const DVector &iv2 )
{
    DVector result = iv2 - iv1;
    IVector sb ( result );
    return sb;
}
// To comply with overloaded PROFIL function Hull() which takes
// interval arguments. This is basically a wrapper around CAPD's
// intervalHull() which returns a BOX (==IVector)
interval Hull ( const double &d1, const double &d2 )
{
    return interval ( d1, d2 );
}

interval Hull ( const interval &iv1, const interval &iv2 )
{
    return capd::intervals::intervalHull ( iv1, iv2 );
}

interval Hull ( const double &dbl, const interval &iv )
{
    return capd::intervals::intervalHull ( interval ( dbl ), iv );
}

interval Hull ( const double &dbl )
{
    return interval ( dbl );
}

// END OVERLOADING OF HULL()

double Mig ( const interval &iv )
{
    if ( iv.leftBound() > 0 )
      return iv.leftBound();
    else if ( iv.rightBound() < 0 )
      return iv.rightBound();
    else
      return 0.;
}

double Abs ( const interval &iv )
{
    return iabs ( iv ) . leftBound();
}

// Max (Min) on a BOX assume that we are passing in an IVector of
// diameters obtained from diam(). Thus, we have SYSDIM intervals {
// [a,a], ..., [d,d] }, and we need to find the maximum (minimum). So
// we simply iterate through the intervals in the vector.
double Max ( const BOX &box )
{
    double theMax = box[0] . leftBound ( ); // singleton intervals
    for ( register int i = 0; i < box . dimension(); ++i )
      { 
	if ( box [ i ] . leftBound () > theMax )
	  theMax = box [ i ] . leftBound ( );
      }
    return theMax;
}

// This is a 'normal' min function
double Min ( const double &d1, const double &d2 )
{
    return ! ( d2 < d1 ) ? d1 : d2;    
}

// This is an interval norm
interval Norm2 ( const BOX &bx )
{ 
    return bx.euclNorm();
}

double Mid ( const interval &iv )
{
    return iv . mid ( ) . leftBound ( ) ;
}


// Intersection: IVector
bool Intersection ( IVector &overlap, const IVector &x, const IVector &y )
{
    try
      {
	overlap = intersection ( x, y );
	return true;
      }
    catch ( ... )
      {
	return false;
      }
}

// Intersection: IMatrix
bool Intersection ( IMatrix &overlap, const IMatrix &x, const IMatrix &y )
{
    return intersection ( x, y, overlap );
}

double Diam ( const interval &iv )
{
    return diam ( iv ) . leftBound();
}

interval DivBounds ( const double &d1, const double &d2 )
{
    interval iv ( d1 / d2 );
    return iv;
}

IVector DivBounds ( const DVector &dvec, const double &d )
{
    DVector vec_copy = dvec;
    vec_copy /= d;
    IVector ivec ( vec_copy ); 
    return ivec;
}

// Return column 'idx' of the matrix 'mat'
IVector Col ( const IMatrix &mat, const int &idx )
{
    IVector iv ( mat . column ( idx ) );
    return iv;
}

void SetCol ( IMatrix &mat, const int &col_num, const IVector &iv )
{
    for ( register int i = 1; i <= mat.numberOfRows(); ++i )
      mat ( i, col_num ) = iv [ i ];
}
////////////////////////////////////////////////////////////////////
