
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

//#include "capd/capdlib.h"
#include "capd/rounding/DoubleRounding.h"
#include "capd/intervals/DoubleInterval.h"
#include "capd/intervals/IntervalError.h"

using namespace std;

//typedef long INTERVAL;

//typedef capd::intervals::Interval< double > DInterval;
typedef int INTERVAL;
//typedef capd::vectalg::Vector< interval, 2> IVector;

namespace capd{
namespace test{


void basicsTest()
{
   interval a,
            b(1.0),
            c(2.0, 3.0),
            d(c),
            e(2.5, 3.0);
            
   cout <<"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n Base functions: constructors, acces to interval bounds, inclusions, relations";
   cout << " a() = " << a << "\n b(1.0) = " << b << " c(2.0,3.0) = " << c 
        << " d(c) = " << d << " e(2.5,3) = " << e 
        << "\n\n c.leftBound() = "  << c.leftBound() 
        << "   c.left() = "   << c.left()
        << "   left(c) = "  << left(c)
         << "\n c.rightBound() = "  << c.rightBound()       
        << "   c.right = "   << c.right()
        << "   right = "  << right(c)
        << "\n\n b==c : " << (b==c)  << "  b!=c : " << (b!=c) 
        << "  b>c : " << (b>c)    << "  b>=c : " << (b>=c) 
        << "  b<c : " << (b<c)    << "  b<=c : " << (b<=c) 
        << "\n\n b==1.0 : " << (b==1.0)  << "  b!=1.0 : " << (b!=1.0) 
        << "  b>1.0 : " << (b>1.0)    << "  b>=1.0 : " << (b>=1.0) 
        << "  b<1.0 : " << (b<1.0)    << "  b<=1.0 : " << (b<=1.0)
        << "\n\n " <<  c <<" contains " << e <<" :  " << c.contains(e)
        << "     "<<c <<" containsInInterior " << e << " : " << c.containsInInterior(e)
        << "\n "<<c << " contains 2.5 : " << c.contains(2.5)
        << "\n "<<d << " subset " << c << " : " << d.subset(c)
        << "   " << d << " subsetInterior " << c << " : " << d.subsetInterior(c)
         ;
   std::istringstream myStr("[3.21312312, 4.324324324]");
   myStr >> a;
   cout << "\n\n interval read from string \"[3.21312312, 4.324324324]\" = " << a;
   
   c.split(c,b);

   cout << "\n\n " << d << " = "  << c << "+" << b;
   double r;
   c = d;
   split(c, r);
   cout << "\n " << d << " = K(" << c << "," <<  r << ")";
   cout << endl;

}
 
void operatorsTest()
{
  interval a,
            b(1.0), 
            c(2.0, 3.0),
            d(c),
            e(-1.0, 4.0),
            f(4.0, 5.0);
   cout <<"\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
  
   cout << " a = " << a << "\n b = " << b << " c = " << c 
        << " d = " << d << " e = " << e << " f = " << f;
   
   cout << "\n c += e " << (c += e);
   cout << "   c -= e " << (c -= e);
   cout << "\n c *= e " << (c *= e);
   cout << "   c /= f " << (c /= f);
   cout << "\n b *(c+(-d)) / f - e = " << (b *(c+(-d))/ f  - e); 
   interval y = interval(-3, -2);
   cout << "\n y+3 = " << y+3 << " y-3 = " << y-3 << " y*3 = " << y*3<< " y/3 = " << y/3
        << "\n 3+y = " << 3+y << " 3-y = " << 3-y << " 3*y = " << y*3<< " 3/y = " << 3/y;

   cout << endl;
}



bool rodesTests()
{
    interval a ( 1.0, 3.0 );

    cout << "a = " << a << endl;

    IVAL _r = diam( a ) / 2.0;
    double r = _r . leftBound ( );
    IVAL symhull ( -r, r );
    cout << "symhull = " << symhull << endl;

    IVAL c( 1.0, 2.0 );
    cout << "c contained in a: " << ( c.subset( a ) ) << endl;

    // IVector Box; //(2);
    // Box[0] = a ;
    // Box[1] = c ;
    // cout << "the box " << Box << endl;

    cout << "PI " << IVAL::pi() << endl;

    return true;
}

// class BOX::IVector
// {
// public:
//   int DIM = 2;
//   BOX ( interval a, interval b );
// };

// BOX::BOX ( interval a, interval b ) 
// { 
    


int main()
{
  try{
    // basicsTest();
    //operatorsTest();
    //functionsTest();
   interval a(3, -1);
   //scalarTest();
   //multiplicationTest();
   //multiplicationTest2();

   //divisionTest();
   //divisionTest2();

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
