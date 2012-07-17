/*   File: 2d_classes.h

     Contains classes: grid, data_line, 
     and iterate. These are used only in
     the high level modules of RODES.

     Latest edit: Mon Jul 3 2000
*/

#ifndef TWO_D_CLASSES_H
#define TWO_D_CLASSES_H 

#include <iostream.h>
#include <iomanip.h>
#include <math.h>

#include "classes.h"

// Enumerations used for data_line status.
enum {NOT_DONE, BEING_DONE, DONE, DO_AGAIN, FAILED, RESERVED};
enum {NOT_HIT, HIT};

static const int    OUTPUT_PRECISION = 17; // string - double conversion exact.
static const int    GRID_WIDTH       =  5; // White padding for output.
static const double LARGE_NUMBER  = 1E+90; // Used for finding minimals.

////////////////////////////////////////////////////////////////////

class grid
{
 public:
  int u;         // Lower u coorninate's corner
  int v;         // Lower v coorninate's corner
  int P;         // Scale = 2^-P  

  bool operator == (const grid &rhs) const
    { 
      if ( P == rhs.P )
	{
	  if ( u == rhs.u && v == rhs.v )   // We use the symmetry
	    return true;                    // of the Lorenz eq's.
	  if ( u == -rhs.u && v == -rhs.v )
	    return true;
	}
      return false;
    }
  bool operator != (const grid &rhs) const
    {
      if ( *this == rhs )
	return false;
      return true;
    }
  friend BOX grid_to_box(const grid &g)
    {
      BOX result(2);
      result(1) = interval(g.u - 1, g.u + 1);
      result(2) = interval(g.v - 1, g.v + 1);
      result *= pow(2, -g.P);
      return result;
    }
  friend ostream & operator << (ostream &out, const grid &g)
    {
      out.setf(ios::showpos);
      out << setw(GRID_WIDTH) << g.u << " " << setw(GRID_WIDTH) << g.v;
      out.unsetf(ios::showpos);
      out << " " << g.P;

      return out;
    }
  friend istream & operator >> (istream &in, grid &g)
    {
      in >> g.u >> g.v >> g.P;
      return in;
    }
};

const grid NULL_GRID = {0, 0, 0};

////////////////////////////////////////////////////////////////////

class new_data_line
{
public:
  grid grd;
  int c_stat;     // Computational status 
  int h_stat;     // Hit status
#ifdef COMPUTE_C1
  interval ang;   // Cone angles
  double pre_exp; // Pre-expansion
  double min_exp; // Minimal expansion
#endif

  bool operator == (const new_data_line &rhs) const
    { 
      return ( grd == rhs.grd && c_stat == rhs.c_stat && h_stat == rhs.h_stat );
    }
  bool operator != (const new_data_line &rhs) const
    {
      return (grd != rhs.grd);
    }
  friend BOX ndl_to_box(const new_data_line &ndl)
    {
      return grid_to_box(ndl.grd);
    }
  friend ostream & operator << (ostream &out, const new_data_line &ndl)
    {
      out << ndl.grd << "   " << ndl.c_stat << " " << ndl.h_stat;
#ifdef COMPUTE_C1
      interval angles = RAD_TO_DEG * ndl.ang;

      out.precision(OUTPUT_PRECISION);
      out.setf(ios::showpos);
      out.setf(ios::scientific); 
      out << "   " << Inf(angles) << " " << Sup(angles)
	  << " " << ndl.pre_exp << " " << ndl.min_exp;
      out.unsetf(ios::showpos);
      out.unsetf(ios::scientific);
#endif // COMPUTE_C1
      return out;
    }
  friend istream & operator >> (istream &in, new_data_line &ndl)
    {
      in >> ndl.grd >> ndl.c_stat >> ndl.h_stat;
#ifdef COMPUTE_C1
      double lo, hi;
      in >> lo >> hi;
      ndl.ang = DEG_TO_RAD * Hull(lo, hi);
      in >> ndl.pre_exp >> ndl.min_exp;
#endif
      return in;
    }
};

////////////////////////////////////////////////////////////////////

class iterate
{
 public:
  new_data_line ndl;
  grid inf_grd;      // The left-most grid of the return of dl.grd
  grid sup_grd;      // The right-most grid of the return of dl.grd

  friend ostream & operator << (ostream &out, const iterate &it)
    {
      out << it.ndl << "   " << it.inf_grd << "   " << it.sup_grd;
      return out;
    }
  friend istream & operator >> (istream &in, iterate &it)
    {
      in >> it.ndl >> it.inf_grd >> it.sup_grd;
      return in;
    }
};

////////////////////////////////////////////////////////////////////

#endif // TWO_D_CLASSES_H

