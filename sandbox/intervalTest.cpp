
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

// #include "capd/intervals/lib.h"
#include "capd/rounding/DoubleRounding.h"
#include "capd/intervals/DoubleInterval.h"
#include "capd/intervals/IntervalError.h"
#include "capd/vectalg/Vector.hpp"
#include "capd/vectalg/Matrix.hpp"


typedef capd::intervals::Interval< double > DInterval;
typedef capd::vectalg::Vector< interval, 0 > IVector;
typedef capd::vectalg::Vector< double, 0 > DVector;
typedef capd::vectalg::Matrix< double, 0, 0 > DMatrix;
typedef capd::vectalg::Matrix< DInterval, 0, 0 > IMatrix;

#define BOX IVector 

using namespace std;
namespace capd{
  namespace test{

class parcel 
{
public:
  BOX box;
};


    bool func ( interval &c )
    {
      try
	{
	  double num = 1.;
	  c = num / interval ( 2. );
	  return true;
	}
      catch ( ... ) 
	{
	  return false;
	}
    }

bool rodesTests()
{
    interval a( 2. );
    interval c( 0. );
  
    double dbl = 1.;
    double level = 27;
    double trvl_dx = 0.0011;
    interval result = interval( level - trvl_dx );

    cout << " result = " << level - trvl_dx << endl;
    

    interval arr[] = { interval( 1 ), interval( 2 ), interval( 3 ) };
    IVector ivec ( 3, arr );
   
    cout << "dbl / ivec[1] = " << dbl / ivec[1] << endl;
    
    interval d1[] = {interval(-1.,1.), interval(2.,2.), interval(3.,3.1), interval(4.,4.1)};
    IVector v1(4,d1);


    return true;
}

int main()
{
  try{

    rodesTests();

   } catch(capd::intervals::IntervalError<double> &e)
   {
      cout <<" interval Error : " << e.what();
   }

   catch(std::runtime_error &e)
   {
      cout << e.what();
   }

   return 0;
}


}}// namespace capd::test

int main()
{
   return capd::test::main();
}
