/**********************************************************************
 *
 * File:          smalldiv.cc
 *
 * Program:       smalldiv
 *
 * Author:        Warwick Tucker
 *
 * Date:          991124
 *
 * Purpose:       Calculates and expresses the smallest divisors
 *                up to a given maximal order and smoothness. 
 *                Uses the interval class from PROFIL/BIAS.
 *
 * Project:       The Lorenz Attractor Exists
 *
 * Compile:       make smalldiv
 *
 * Usage:         smalldiv                
 *
 * Input:         In my paper, I use smoothness = 10, and
 *                maximal order of resonance = 61. 
 *
 * Modified: 120701, by jjb
 *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

// #include "Interval.h"          /* PROFIL/BIAS header */
// #include "Functions.h"         /* PROFIL/BIAS header */

#include "capd/capdlib.h"       // jjb 

#define ANSW_SIZE 30
#define FALSE 0
#define TRUE 1

// CAPD interval typedef
typedef PointBase< capd::intervals::Interval< double > > interval;

interval lu, lss, ls;           /* The eigenvalues of the origin */      
interval R, S, B;               /* The parameter values */ 
interval uMin,ssMin,sMin;       /* The smallest divisors */
int total;                      /* Total number of iterates */ 
int MaxOrder, Smoothness;       /* Desired maximal order and smoothness */


/* Initialize the constants and ask for the  */
/* maximal order of resonance and smoothness */
void Init()        
{                   
  char answer[ANSW_SIZE];  
  double fl_R, fl_S, fl_B;

  // Convert to CAPD
    R = interval ( 28.0 ); 
    S = interval ( 10.0 ); 
    B = interval ( 8.0 ) / 3.0; 
  // R = Succ(Hull(28.0));
  // S = Succ(Hull(10.0));
  // B = Succ(Hull(8.0/3));
  printf(" \n");
  printf("**************************************************************\n\n");
  printf("Enter the desired smoothness (10): ");
  fgets(answer, sizeof(answer), stdin); 
  sscanf(answer, "%d", &Smoothness);
  printf("Enter the maximal order of resonance (61): ");
  fgets(answer, sizeof(answer), stdin); 
  sscanf(answer, "%d", &MaxOrder);
  printf("Do you want to change parameters? (y/n) ");
  fgets(answer, sizeof(answer), stdin); 
  if (answer[0] == 'y')
    {
      printf("Enter R: ");
      fgets(answer, sizeof(answer), stdin); 
      sscanf(answer, "%lf", &fl_R);   
      printf("Enter S: ");
      fgets(answer, sizeof(answer), stdin); 
      sscanf(answer, "%lf", &fl_S);  
      printf("Enter B: ");
      fgets(answer, sizeof(answer), stdin); 
      sscanf(answer, "%lf", &fl_B); 
      R = interval ( fl_R ); //Succ(Hull(fl_R));
      S = interval ( fl_S ); //Succ(Hull(fl_S));
      B = interval ( fl_B ); //Succ(Hull(fl_B));
    }
  lu  = ( - (S + 1) + Sqrt((S + 1)*(S + 1) + 4*S*(R - 1)))/2;
  lss = ( - (S + 1) - Sqrt((S + 1)*(S + 1) + 4*S*(R - 1)))/2;
  ls  = - B;
}

/* Compare all appearing divisors   */
/* and single out the smallest ones */
void Calculate()              
{                            
  register int n1, n2, n3;       /* The coefficients of resonance */
  int order;                     /* The present order of resonance */ 
  int counter1, counter2;        /* Two basic counters */
  int uMini=0,uMinj=0,uMink=0;   /* Minimizing coefficients wrt lu */
  int ssMini=0,ssMinj=0,ssMink=0;/* Minimizing coefficients wrt lss */
  int sMini=0,sMinj=0,sMink=0;   /* Minimizing coefficients wrt ls */
  interval divisor;              /* Current divisor */
           
  cout.precision(16);
  total = 0;
  for (order = 2; order <= MaxOrder; order++) /* We check all orders beween 2 and MaxOrder */   
    {
      uMin = -order*lss;          /* Guess a minimum resonance wrt lu */
      ssMin = -order*lss;         /* Guess a minimum resonance wrt lss */
      sMin = -order*lss;          /* Guess a minimum resonance wrt ls */
      for (counter1 = 0; counter1 <= order; counter1++) 
	{
	  n1 = order - counter1;            /* |n| = order */
	  for (counter2 = 0; counter2 <= counter1; counter2++)
	    {
	      n2 = order - n1 - counter2;   
	      n3 = order - n1 - n2;           
	      if ( (n1 < Smoothness) || (n2 + n3 < Smoothness) )  /* The filter */
		{
		  interval temp;
		  temp    = n1*lu + n2*lss + n3*ls - lu;
		  divisor = Hull(Mig(temp), Abs(temp));
		  if ( Inf(divisor) < Inf(uMin) ) 
		    {
		      uMin = divisor;
		      uMini = n1;
		      uMinj = n2;
		      uMink = n3;
		    }
		  temp    = n1*lu + n2*lss + n3*ls - lss;
		  divisor = Hull(Mig(temp), Abs(temp));
		  if ( Inf(divisor) < Inf(ssMin) )
		    {
		      ssMin = divisor;
		      ssMini = n1;
		      ssMinj = n2;
		      ssMink = n3;
		    }
		  temp    = n1*lu + n2*lss + n3*ls - ls;
		  divisor = Hull(Mig(temp), Abs(temp));
		  if ( Inf(divisor) < Inf(sMin) ) 
		    {
		      sMin = divisor;
		      sMini = n1;
		      sMinj = n2;
		      sMink = n3;
		    }
		  total++;
		}
	    }
	}
      divisor = uMin;
      if ( Inf(ssMin) < Inf(divisor) )
	divisor = ssMin;
      if ( Inf(sMin) < Inf(divisor) )
	divisor = sMin;     /* divisor = Omega(order) */	  
      printf("Smallest divisors of order %d:\n", order);
      printf("%d*lu + %d*lss + %d*ls - lu  = ", uMini,uMinj,uMink);
      cout << uMin << " diam = " << diam(uMin) << endl;
      printf("%d*lu + %d*lss + %d*ls - lss = ", ssMini,ssMinj,ssMink);
      cout << ssMin << " diam = " << diam(ssMin) << endl;
      printf("%d*lu + %d*lss + %d*ls - ls  = ", sMini,sMinj,sMink);
      cout << sMin << " diam = " << diam(sMin) << endl << endl;
    }
}

     
int main()
{
  Init();                
  cout << endl << "Parameters:" << endl;
  cout << "S = " << S << endl;
  cout << "R = " << R << endl;
  cout << "B = " << B << endl;
  cout << "lu  = " << lu << endl;
  cout << "lss = " << lss << endl;
  cout << "ls  = " << ls << endl;
  printf("Maximal order %d:\n", MaxOrder);
  printf("Smoothness: %d \n\n", Smoothness);
  Calculate();
  printf("Total amount of iterations: %d \n\n", total);
  printf("**************************************************************\n\n");
  printf("%c", 7);    /* Signal when done */
  return 0;
}

