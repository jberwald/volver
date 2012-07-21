/*   File: low_functions.cc (BIAS)

     Contains all low level functions
     needed to compute a return.

     Latest edit: Mon Apr 10 2000
*/

#include "low_functions.h"

////////////////////////////////////////////////////////////////////

static void LU_Decompose    (IMatrix &, const IMatrix &, int *);

static void LU_Backsub      (IVector &, const IMatrix &, int *,
			     const IVector &);

static void Invert_And_Mult (IMatrix &, const IMatrix &,
			     const IMatrix &);

////////////////////////////////////////////////////////////////////

static IMatrix Id (const short &dim)
{
    IMatrix I(dim, dim);

    I.clear(); 
    for (short i = 1; i <= dim; i++) 
      I(i,i) = 1.0;
    
    return I;
}

// ID declared as a fixed global variable to save time.
static const IMatrix ID = Id(SYSDIM); 

////////////////////////////////////////////////////////////////////

int Sign(const double &dbl)
{ return (dbl >= 0 ? 1 : - 1); }

////////////////////////////////////////////////////////////////////

int Int_Power(const int &base, const int &n)
{
  register int p = 1;
  for ( register int i = 1; i <= n; i++ )
    p *= base;
  
  return p;
}

////////////////////////////////////////////////////////////////////

// Called by 'Flow_By_Corner_Method'.
// Computes the flow using corners/sides when at least one
// partial derivative off the diagonal of DPi may vanish.
void Some_May_Vanish(BOX &result, const IMatrix &DPi, 
		     const parcel &pcl, const BOX &Outer_Box, 
		     const BOX &dx, const double &sign_trvl_dist)
{
    short trvl      = pcl.trvl;     // Shorthand
    BOX   Inner_Box = pcl.box;
    DVector dInf, dSup;  // for storage of Inf and Sup results

    // Array of two IVectors of length SYSDIM
    // IVector::setDefaultDimension( SYSDIM ); // This should set the
					    // dimension for the other
					    // IVectors below. ? What
					    // is the scope of this
					    // declaration ?
    BOX *Poincare = new ( SYSDIM ) IVector [ 2 ]; // we want two vectors with number of spatial coords = SYSDIM
  
    // trvl \in {1,2,...SYSDIM}, so shift trvl by -1
    if ( pcl.sign == - 1 )     // Set the trvl coordinates
      // we want the singleton, so we use an interval (Hull in PROFIL)
      Poincare [ 0 ][ trvl-1 ] = interval ( Inf ( Outer_Box [trvl-1] ) ); // probably need operator[]
    else
      Poincare [ 0 ][ trvl-1 ] = interval ( Sup ( Outer_Box [trvl-1] ) );
    Poincare [ 1 ][ trvl-1 ] = Poincare [ 0 ][ trvl-1 ]; 
      
    double   start;
    interval local_dist;
    BOX      vf;
    BOX *Side_Box = new ( SYSDIM )IVector [ 2 ];
  
  for (register short i = 1; i <= SYSDIM; i++)  
    if ( i != trvl )         
      { // Loop through i (i is the component-coordinate)

	Side_Box[0] = Outer_Box;	// Make Side_Box (lower)	
	Side_Box[0][i] = Inf ( Inner_Box [ i ] ) + dx [ i ];
	
	Side_Box[1] = Outer_Box;	// Make Side_Box (upper) 
	// apply necessary distortion to box
	Side_Box[1] [ i ] = Sup(Inner_Box[i]) + dx[i];	    
	
	for (register short j = 1; j <= SYSDIM; j++) 
	  { // Loop through j (j is the partial derivative-coordinate)
	    // Skip i = j (we know that the diagonal elements are non-zero) 
	    // Skip j = trvl (we already know what codim1 we will hit)
	    if ( (i != j) && (j != trvl) )  
	      for (register short k = 0; k < 2; k++)
		{ // Loop through k (k represents upper and lower box)	
		  if ( k == 0 )	 
		    start = Inf(Inner_Box[i]); // Get the initial value
		  else
		    start = Sup(Inner_Box[i]);
		  
		  if ( Subset(0.0, DPi(i, j)) == true )	
		    {    // compute min/max P[i] in the side boxes  
		      Vf_Range(vf, Side_Box[k]);
		      local_dist = sign_trvl_dist * vf[i] / vf(trvl);  // interval
		      Poincare[k][i] = start + local_dist; // double + interval
		    }   
		  else   // Compute min/max P[i] at the corners
		    { // Construct the small corner boxes 
		      //BOX Corner_Box ( 2 ); 
		      BOX *Corner_Box = new ( SYSDIM ) IVector [ 2 ];
		      interval iv[2];
  
		      Corner_Box[0] = Side_Box[k];
		      Corner_Box[0][ j ] = Inf(Inner_Box[ j ]) + dx[ j ];    
		      Corner_Box[1] = Side_Box[k];
		      Corner_Box[1][ j ] = Sup(Inner_Box[ j ]) + dx[ j ];	       

		      for (register short m = 0; m < 2; m++)
			{ // m indicates the Corner_Box in use
			  Vf_Range(vf, Corner_Box[m]);  
			  iv[m] = vf[i] / vf(trvl);	
			}	    
		      Poincare[k][i] = start + sign_trvl_dist 
			* intervalHull(iv[0], iv[1]);
		    } // Done with the corners
		}
	  }
      }
  // the hull of [a,a] and [b,b] 
  // Instead of adding function, just do the conversion here
    dInf = Inf ( Poincare[0] ); 
    dSup = Sup ( Poincare[1] );
    IVector pInf ( dInf );
    IVector pSup ( dSup );
    result = intervalHull ( pInf, pSup );
}

////////////////////////////////////////////////////////////////////

// Called by 'Flow_By_Corner_Method'.
// Computes the flow using corners/sides when none of the
// partial derivatives off the diagonal of DPi vanishes. 
void None_May_Vanish(BOX &result, const parcel &pcl, const BOX &Outer_Box, 
		     const BOX &dx, const double &sign_trvl_dist)
{
    register short i, j;
    double lo_coord, hi_coord;
    // allotting space for array of DVectors ( n = POWER )
    DVector *point_in = new ( SYSDIM ) DVector [ POWER ];
    DVector *point_out = new ( SYSDIM ) DVector [ POWER ];

    // DVector point_in[POWER];   // Storage for the corner points -- we need 2^(d-1) corners for the in points
    // DVector point_out[POWER];  // Storage for the corner points  << --- ***** THIS UTILIZES BUILT IN Int_Power() Will have to fils ASAP.
    DVector current_point(SYSDIM); // 

    BOX *corner_in = new ( SYSDIM ) BOX [ POWER ];
    BOX *corner_out = new ( SYSDIM ) BOX [ POWER ];

    // BOX corner_in[POWER];     // Storage for the corner boxes -- similar number of corner intervals needed as in point_in
    // BOX corner_out[POWER];    // Storage for the corner boxes 
    BOX current_box(SYSDIM);
    BOX lo_box(SYSDIM), hi_box(SYSDIM);
    int corner_in_cnt  = 0;
    int corner_out_cnt = 0;
    int point_in_cnt   = 0;
    int point_out_cnt  = 0;
    int loop_cnt;

  short trvl = pcl.trvl;    // Shorthand
  BOX Inner_Box = pcl.box;

  corner_out[corner_out_cnt++] = Outer_Box;
  point_out[point_out_cnt++] = Inf(Inner_Box); // Inf = Sup for trvl coordinate.
  
  for (i = 1; i <= SYSDIM; i++)
    if ( i != trvl )
      {	// Prepare the intersectors
	lo_box    = Outer_Box;
	lo_box[i] = Inf(Inner_Box[i]) + dx[i];
	hi_box    = Outer_Box;
	hi_box[i] = Sup(Inner_Box[i]) + dx[i];
	
	lo_coord = Inf(Inner_Box[i]);
	hi_coord = Sup(Inner_Box[i]);
	
	for (j = 0; j < corner_out_cnt; j++)  // Copy corner_out to corner_in
	  {	                              // and delete all in corner_out
	    corner_in[j] = corner_out[j];
	    point_in[j] = point_out[j];
	  }
	corner_in_cnt = corner_out_cnt;
	corner_out_cnt = 0;
	point_in_cnt = point_out_cnt;
	point_out_cnt = 0;
	
	loop_cnt = corner_in_cnt;
	for (j = 0; j < loop_cnt; j++) // loop_cnt = 2^(i - 1)
	  {
	    current_box = corner_in[j];
	    Intersection(corner_out[corner_out_cnt++], current_box, lo_box);
	    Intersection(corner_out[corner_out_cnt++], current_box, hi_box);
	    current_point = point_in[j];
	    point_out[point_out_cnt] = current_point;
	    point_out[point_out_cnt++][i] = lo_coord;
	    point_out[point_out_cnt] = current_point;
	    point_out[point_out_cnt++][i] = hi_coord;
	  }               // Now the corners are stored in 
      }                   // corner_out[i], i = 0,...,2^(SYSDIM - 1) - 1.

  i = 0; // Just for symmetric definitions
  BOX vf(SYSDIM);  Vf_Range(vf, corner_out[i]);
  interval corner_time = sign_trvl_dist / vf(trvl);

  // This is odd: DVector + ( interval * IVector )
  // From PROFIL doc: REAL or VECTOR operands may also be used instead of interval operands, as long as one operand is an INTERVAL type. 
    
  // Convert DVector to IVector (singletons)
    IVector dp ( point_out[i] );
    result = dp + ( corner_time * vf ); // A bit wasteful seeing that I don't use the trvl coord.
                                            // Probably faster to use add an inner loop:
    // Use dp here, too
  for (i = 1; i < POWER; i++)               // for (k = 1; k <= SYSDIM; k++)
    {                                       //   if ( k != trvl )
      Vf_Range(vf, corner_out[i]);          //     result(k) = point_out[i](k) + corner_time * vf(k);
      corner_time = sign_trvl_dist / vf(trvl);   
      result = intervalHull(result, dp + corner_time * vf);  // ** this is passing ivec and interval
    }

  if ( pcl.sign == 1 )
    result(trvl) = interval(Sup(Outer_Box(trvl)));
  else
    result(trvl) = interval(Inf(Outer_Box(trvl)));
}

////////////////////////////////////////////////////////////////////

// Called by 'Flow'.
// Flows pcl.box to the opposite transversal side of Outer_Box.
void Flow_By_Corner_Method(BOX &Result_Box, const IMatrix &DPi,
			   const parcel &pcl, const BOX &Outer_Box)
{
  // Compute the maximal change any point in Inner_Box can have while
  // flowing.  Store the convex hull of [inf,inf] and [sup,sup]. i.e.,
  // create an intervalHull() from the vectors of singleton intervals.
    BOX dx = capd::vectalg::intervalHull ( SubBounds( Inf( Outer_Box ), Inf( pcl.box ) ),
					   SubBounds( Sup( Outer_Box ), Sup( pcl.box ) ) );

  bool zero_indicator = false;     // Check for possible zeroes of the
  for (register short i = 1; i <= SYSDIM; i++)    
    {
      if ( i != pcl.trvl )           // partial (space) derivatives (DPi). 
	{
	  for (register short j = 1; j <= SYSDIM; j++)
	    {
	      if ( j != pcl.trvl )
		{
		  if ( Subset(0.0, DPi(i, j)) )
		    {
		      zero_indicator = true;
		      if ( i == j )
			{
			  char *msg = "Error: 'Flow_By_Corner_Method'"
			    " DPhi vanishes on the diagonal!";
			  throw Error_Handler(msg);
			}
		    }
		}
	    } // END  for j
	} 
    } // END for i
    double trvl_dist = Diam( dx( pcl.trvl ) );   // Get the transversal distance
    if ( pcl.sign == -1 )
      trvl_dist = - trvl_dist;

  if ( zero_indicator )
    Some_May_Vanish(Result_Box, DPi, pcl, Outer_Box, dx, trvl_dist);
  else
    None_May_Vanish(Result_Box, pcl, Outer_Box, dx, trvl_dist);    

}

////////////////////////////////////////////////////////////////////

static void LU_Decompose(IMatrix &R, const IMatrix &A, int *indx)
{
  register int i, j, k, imax;
  double big, dummy;
  static interval temp, sum;
  static IVector vv ( SYSDIM );

  R = A; // Copy A into R

  for (i = 1; i <= SYSDIM; i++) // Loop over all rows to get
    {                        // the implicit scaling info.
      big = 0.0;
      for (j = 1; j <= SYSDIM; j++)
	if ( (dummy = Abs(R(i,j))) > big ) 
	  big = dummy;
      if ( big == 0.0 )
	{
	  char *msg = "Error: 'LU_Decompose'. Singular Matrix";
	  throw Error_Handler(msg);
	}
      vv[i] = Sup ( DivBounds ( 1.0, big ) );  // Store the scaling
    }

  // Loop over columns for Crout's method.
  for (j = 1; j <= SYSDIM; j++)    
    {                      
      imax = j; // This row is missing in "Numerical Recipies in C".
      for (i = 1; i < j; i++)  // Solve for the upper elements.
	{
	  sum = R(i,j);
	  for (k = 1; k < i; k++)
	    sum -= R(i,k) * R(k,j); 
	  R(i,j) = sum;
	}
      big = 0.0;// Initialize the search for the biggest pivot element.
      for (i = j; i <= SYSDIM; i++)   // Solve for the diagonal and upper elements.
	{ 
	  sum = R(i,j);
	  for (k = 1; k < j; k++)
	    sum -= R(i,k) * R(k,j); 
	  R(i,j) = sum;
	  if ( (dummy = Abs(vv[i] * sum)) > big )
	    {                 // If the pivot is the best so far
	      big = dummy;    // we save it and its position.
	      imax = i;
	    }
	}
      if ( j != imax )       // Do we need to interchange rows?
	{
	  for (k = 1; k <= SYSDIM; k++) // If yes, we do so...
	    {
	      temp       = R(imax, k);
	      R(imax, k) = R(j, k);
	      R(j, k)    = temp;
	    }
	  vv [ imax ] = vv [ j ]; // ...and interchange the the scale factor
	}
      indx[j] = imax;
      if ( Subset(0.0, R(j, j)) ) // If the pivot element contains zero...
	{
	  char *msg = "Error: 'LU_Decompose'. Singular Matrix - "
	    "pivot element contains zero.";
	  throw Error_Handler(msg);
	}
      if (j != SYSDIM)
	{
	  temp = 1.0 / R(j, j);   // Now divide by the pivot element
	  for (i = j + 1; i <= SYSDIM; i++)
	    R(i, j) *= temp;
	}
    }                            // Get the next column in the decomposition.
}

////////////////////////////////////////////////////////////////////

static void LU_Backsub( IVector &r, const IMatrix &A, int *indx,
			const IVector &b )
{
  register int i, j, ip;
  static interval sum;

  r = b;
  for (i = 1; i <= SYSDIM; i++)  // Start with the forward substitution.
    {
      ip = indx[i];           // We unscramble the permutation  
      sum = r(ip);            // as we go...
      r(ip) = r [ i ];
      if ( i != 1 )
	for (j = 1; j <= i - 1; j++)
	  sum -= A(i, j) * r[ j ];
      r[i] = sum;
    }

  for (i = SYSDIM; i >= 1; i--) // Now we do the backsubstitution.
    {
     sum = r[i];  
     for (j = i + 1; j <= SYSDIM; j++)
       sum -= A(i, j) * r[ j ];
     r[i] = sum / A(i, i);  // Store the component of the solution vector.
    }
}

////////////////////////////////////////////////////////////////////

// Computes 'A^-1 * B' faster than 'Inverse(A) * B'.
static void Invert_And_Mult(IMatrix &Result, const IMatrix &A,
			    const IMatrix &B)
{
    static int indx[SYSDIM + 1];
    static IVector Result_IV;  // jjb -- this holds results in
			       // LU_Backsub(), so no need to fix a size
    static IMatrix LU_IM ( SYSDIM, SYSDIM );

  LU_Decompose(LU_IM, A, indx);
  for (register short j = 1; j <= SYSDIM; j++)
    {
      // *** Col ( ) defined in classes.h, wraps IMatrix method. ***
      LU_Backsub( Result_IV, LU_IM, indx, Col ( B, j ) );
      SetCol( Result, j, Result_IV );
    }
}

////////////////////////////////////////////////////////////////////
// Called by 'Get_DPi_Matrix'.
// Computes DPhi, the matrix of spatial derivatives of 
// the flow. Returns it by reference.
// Assumtion: DPhi(x, s) \subseteq ID + s * Delta_Matrix
// for all s \in [0,t]. This implies that 
//   DPhi - ID = \int_0^t{DVf * (ID + s * Delta_Matrix)}ds
//             = DVf * (t * ID + t^2 / 2 * Delta_Matrix)
//             = t * DVf * (ID + t / 2 * Delta_Matrix).
// Hence, we  just solve for Delta_Matrix:
//   Delta_Matrix = (ID - t / 2 * DVf)^-1 * DVf,
// and finally we set DPhi = ID + t * Delta_Matrix. Pronto!
void Get_DPhi_Matrix(IMatrix &DPhi, const BOX &Outer_Box, 
		     const interval &time)
{
  // Since Delta_Matrix is 'Result' in Invert_And_Mult, set it's size
  // here
    IMatrix Delta_Matrix ( SYSDIM, SYSDIM );  
    IMatrix DVf(SYSDIM, SYSDIM); 
    
    DVf_Range(DVf, Outer_Box);   
    IMatrix tDVf ( Hull( interval(0.0), time) * DVf );

    // make a copy of tDVF, since operator* not working with IMatrix,
    // but *= is working (??)
    //  Want: Invert_And_Mult( Delta_Matrix, (ID - 0.5 * tDVf), tDVf );
    // but we use tDVf_copy *= 0.5 instead.  
    IMatrix tDVf_copy ( tDVf );
    tDVf_copy *= 0.5;

    // jjb -- Delta_Matrix is not sized at this point
    Invert_And_Mult( Delta_Matrix, (ID - tDVf_copy), tDVf );
    IMatrix Exp_M ( ID + Delta_Matrix ); 
    IMatrix Pic ( ID + time * DVf * Exp_M ); 

  if ( !Intersection(DPhi, Pic,  Exp_M) )
    {
      cout << "Empty intersection!!!" << endl;
      DPhi = Pic;
    }
}

////////////////////////////////////////////////////////////////////
