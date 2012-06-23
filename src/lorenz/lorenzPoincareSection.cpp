#include <iostream>
#include "capd/capdlib.h"
#include "capd/poincare/PoincareMapJet.h"
using namespace capd;
using namespace std;

#define DEBUG

///////////////////////////////////////
//
// To compile: g++ -O2 lorenzPoincareSection.cpp -o lorenzPoincareSection `/Users/jberwald/src/capd/bin/capd-config --cflags --libs`
//
// Modified by: Jesse Berwald
//
///////////////////////////////////////
// The following computes the Poincare return map for the Lorenz map,
// with chaotic parameter values sigma=10, beta=8/3, rho=28

// First we define a class for easy computation of Poincare map and its derivative

class LorenzPoincareMap
{
public:
  typedef capd::poincare::PoincareMap<ITaylor> PoincareMap;

  IFunction section;
  IMap vectorField;
  ITaylor solver;
  PoincareMap pm;

  // initialize class with initializer list -- constructs derived objects
  // LorenzPoincareMap ( int order, interval _sigma, interval _rho, interval _beta )
  LorenzPoincareMap ( int order, interval _lambda1, interval _lambda2, interval _lambda3,
		      interval _k1, interval _k2, interval _k3 )
    : section("var:x,y,z;fun:z-0.10546875;"), //  the section is z-27
      //vectorField("par:b,r,s;var:x,y,z;fun:s*(y-x),x*(r-z)-y,x*y-b*z;"), // here is the vector field
      vectorField ( "par:a,b,c,j,k,l;var:x,y,z;fun:a*x-j*(x+y)*z,b*y+j*(x+y)*z,c*z+(x+y)*(k*x+l*y);"),
      solver(vectorField,order,0.1), // 0.1 is the time step, ignored since step control is turned on by default
      // MinusPlus means that the section is crossed with 'z' changing
      // sign from minus to plus Other acceptable values are PlusMinus
      // and None. The last means both directions are acceptable. This will probably have to be changed to None. 
      pm( solver, section, poincare::PlusMinus )  
                                  
  {
    // set parameter values b, r, s
    // vectorField.setParameter("b",_beta);
    // vectorField.setParameter("r",_rho);
    // vectorField.setParameter("s",_sigma);
    // set parameter values for the normal form of Lorenz equations
    vectorField.setParameter("a",_lambda1);
    vectorField.setParameter("b",_lambda2);
    vectorField.setParameter("c",_lambda3);
    vectorField.setParameter("j",_k1);
    vectorField.setParameter("k",_k2);
    vectorField.setParameter("l",_k3);
  }

  // this operator computes period-iteration of Poincare map
  IVector operator()(const IVector& u, int period)
  {
    // u is assumed to be on the section
    // so it is 2-dim, represented by coordinates (y,z).
    // We simply add 0 as the third coordinate i.e. we embed the vector
    // into the full 3d-space
    IVector px(3);
    px[0] = u[0];
    px[1] = u[1];
    px[2] = 0.10546875; // 27*2^{-8}

    cout << "before pm" << endl;
    cout << px << endl;

    // we define a doubleton representation of the set
    C0Rect2Set theSet( px );
    //for(int i=0;i<period;++i) // and compute period-iterations
    px = pm( theSet, period );

    cout << "after pm" << endl;
    cout << px << endl;

    // here we project the image 'x' onto 2-dimensional section
    //return IVector(2,x.begin()+1);

    // If you do not understand what is above, just forget and simply do:
    IVector result(2);
    result[0] = px[0];
    result[1] = px[1];
    return result;
    // What is above is faster and shorter.
  }

  // This operator computes derivative of the Poincare map.
  IMatrix dx(const IVector& u, int period)
  {
    // Again u is two dimensional, so embed it.
    IVector deriv(3);
    deriv[0] = u[0];
    deriv[1] = u[1];
    deriv[2] = 0.10546875; // 27*2^{-8} //interval( 27.0 );

    // for computing of derivative of PM we need an instance of logarithmic norm
    IEuclLNorm N;

    // We define a doubleton representation of the set and its derivative
    // constructor sets initial condition for variational equation to Identity
    C1Rect2Set theSet( deriv, N );

    // the matrix monodromyMatrix will store derivatives of the FLOW not Poincare map
    IMatrix monodromyMatrix(3,3);
    for(int i=0;i<period;++i)
      deriv = pm( theSet, monodromyMatrix );

    // This member function recomputes derivatives of the flow into derivatives of Poincare map
    IMatrix DP = pm.computeDP( deriv, monodromyMatrix );

    # ifdef DEBUG
      cout << "DP matrix " << DP << endl;
    #endif

    // as before, we extract from 3x3 matrix a 2x2 slice
    IMatrix result(2,2);
    result(1,1) = DP(1,1);
    result(1,2) = DP(1,2);
    result(2,1) = DP(2,1);
    result(2,2) = DP(2,2);
    // and return it
    return result;
  }

};


// ---------------------------------------------------------------------
// The following function computes the Interval Newton Operator for the Poincare map on a given set
// and verifies existence and uniqueness of a periodic point in a given set X.
// parameters are:
// pm - an instance of LorenzPoincareMap
// center - a very good approximation of a periodic point on section, center is assume dto be 2-dimensional
// X - interval vector in which we want to prove the existence of an unique periodic point.
// period - is a period of the point for the Poincare map
// Same as above function. Does not verify that N \cap X is
// nonempty. We do truncate the image of the Poincare map so
// computeOrbit() returns N \cap (domain)  -- not yet, but we will.
void computeOrbit( LorenzPoincareMap & pm, 
		   IVector center, 
		   IVector X, 
		   int period,
		   int order )
{
    cout.precision(15);
    // using previously defined class LorenzPoincareMap we compute
    // Poincare Map at the center. Here, the center is the center of a
    // box that we are mapping forward. We return the image
    // boxes. First, report the center that is being mapped.
    #ifdef DEBUG
    cout << endl;
    for ( int i=0; i<2; i++ )
      {
	cout << "interval " << i << " = " << fixed << center[i] << endl;
	cout << "  diam " << i << " = " << diam( center[i] ) << endl;
      }
    #endif

    // We flow forward until we cross the section from the same
    // side that we left from... ??
    IVector imCenter = pm( center, period );
    imCenter = pm( imCenter, period );
  
    for ( int i=0; i<2; i++ )
      cout << "  image[" << i << "] = " << fixed << imCenter[i] << endl;


    // derivative of PM on the set X,
    IMatrix DP = pm . dx( X, period );
    
    // and well known interval Newton operator
    IVector N = center - capd::matrixAlgorithms::gauss( DP - IMatrix::Identity(2),
    							imCenter - center );
  
}
// ----------------------------------- MAIN ----------------------------------------

int main(int argc, char* argv[])
{
  cout.precision(15);
  try{
    // Define an instance of LorenzPoincareMap
    int order = 10; //35;
    int period = 1;
    // interval sigma = interval ( 10. );
    // interval rho = interval ( 28. );
    // interval beta = interval ( 8. ) / interval( 3. );
    interval lambda1 = interval ( 118 ) / interval ( 10. );
    interval lambda2 = interval ( -228 ) / interval ( 10. );
    interval lambda3 = interval ( -267 ) / interval ( 100. );
    interval k1 = interval ( 29. ) / interval ( 100. );
    interval k2 = interval ( 22. ) / interval ( 10. );
    interval k3 = interval ( 13. ) / interval ( 10. );
    //LorenzPoincareMap pm( order, sigma, rho, beta );
    LorenzPoincareMap pm( order, lambda1, lambda2, lambda3, k1, k2, k3 );

    IVector center(2), X(2);
    // For all the orbits below the size of set in which we verify the
    // existence of an unique orbit is set to 7e-13 (up to the rounding:)
    X[0] = X[1] = 3.5e-13*interval(-1,1);


    // box centers
    center[0] = 1255.;
    center[1] = 727.;
    // center[0] = 1229.;
    // center[1] = 723.;
    //center += interval(-1,1) * 2e-1; 
    computeOrbit( pm, center, center, period, order );

   }catch(exception& e)
  {
    cout << "\n\nException caught: "<< e.what() << endl;
  }
  return 0;
} // END
