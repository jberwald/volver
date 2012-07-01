/**********************************************************************
 *
 * File:          coeff.cc
 *
 * Program:       coeff
 *
 * Author:        Warwick Tucker
 *
 * Date:          991124 
 *
 * Purpose:       Calculates the coefficients a_{i,n} of the transformation
 *                in Prop.3.1, computes b_n = max_i |a_{i,n}|, and finally 
 *                constructs c_k = sum_{|n|=k} b_n. Uses the INTERVAL class
 *                from PROFIL/BIAS.
 *
 * Project:       The Lorenz Attractor Exists
 *
 * Compile:       make coeff
 *
 * Usage:         coeff
 *
 * Input:         In my paper, I use smoothness = 10, and 
 *                maximal order of resonance = 70. 
 *
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "Interval.h"             /* PROFIL/BIAS header */
#include "Functions.h"            /* PROFIL/BIAS header */

#define ANSW_SIZE 30              /* Maximal size for answers */
#define BIG 70                    /* Maximal computable order */
#define FALSE 0
#define TRUE 1


INTERVAL lu, lss, ls;             /* The eigenvalues of the origin */ 
INTERVAL k1, k2, k3;              /* Constants in the equations */     
INTERVAL R, S, B;                 /* The parameter values */ 
INTERVAL A[3][BIG][BIG][BIG];     /* The computed coefficients */
INTERVAL C[BIG];                  /* C[k] = sum_{|n|=k} max_i |a_{i,n}| */
INTERVAL T[2][BIG][BIG][BIG];     /* Two sums */
INTERVAL M[2][BIG][BIG][BIG];     /* Two products */
int MaxOrder;                     /* Desired maximal order */ 
int order;                        /* Accumulated order */
int Smoothness;                   /* Desired smoothness */
long int total;                   /* The number of computed coefficients */


/* Initializes the constants and gets the     */
/* desired orders of resonance and smoothness */
void Init()  
{
  char answer[ANSW_SIZE];       /* Temporary storage for answers */
  register int i, j, k;
  double fl_R, fl_S, fl_B;

  R = Succ(Hull(28.0));
  S = Succ(Hull(10.0));
  B = Succ(Hull(8.0/3));
  printf(" \n");
  printf("**************************************************************\n\n");
  printf("Enter the desired smoothness (10): ");
  fgets(answer, sizeof(answer), stdin); 
  sscanf(answer, "%d", &Smoothness);
  printf("Enter the maximal order of resonance (70): ");
  fgets(answer, sizeof(answer), stdin); 
  sscanf(answer, "%d", &MaxOrder);
  printf("Do you want to change parameters? (y/n)");
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
      R = Succ(Hull(fl_R));
      S = Succ(Hull(fl_S));
      B = Succ(Hull(fl_B));
    }
  lu  = ( - (S + 1) + Sqrt((S + 1)*(S + 1) + 4*S*(R - 1)))/2;
  lss = ( - (S + 1) - Sqrt((S + 1)*(S + 1) + 4*S*(R - 1)))/2;
  ls  = - B;
  k1 = S/(Sqrt((S + 1)*(S + 1) + 4*S*(R - 1)));
  k2 = (S - 1 + Sqrt((S + 1)*(S + 1) + 4*S*(R - 1)))/(2*S);
  k3 = (S - 1 - Sqrt((S + 1)*(S + 1) + 4*S*(R - 1)))/(2*S);
  total = 0;
  A[0][1][0][0] = Hull(1); A[0][0][1][0] = Hull(0); A[0][0][0][1] = Hull(0);
  A[1][1][0][0] = Hull(0); A[1][0][1][0] = Hull(1); A[1][0][0][1] = Hull(0);
  A[2][1][0][0] = Hull(0); A[2][0][1][0] = Hull(0); A[2][0][0][1] = Hull(1);
  for (i = 0; i < BIG; i++)
    {
      for (j = 0; j < BIG; j++)
	{
	  for (k = 0; k < BIG; k++)
	    {
	      M[0][i][j][k] = Hull(0.0);
	      M[1][i][j][k] = Hull(0.0);
	    }
	}
      C[i] = Hull(0.0);
    }     
  C[1] = Hull(1);
}


/* Computes the two sums [Phi_1 + Phi_2]_order */
/* and [k2*Phi_1 + k3*Phi_2]_order             */
void Add()     
{
  int counter, index;        /* Basic counters */ 
  register int n1, n2, n3;   /* The exponents  */
                             /* |n| = order    */
  for (counter = 0; counter <= order; counter++)
    {
      n1 = order - counter;     
      for (index = 0; index <= counter; index++)
	{
	  n2 = order - n1 - index;
	  n3 = order - n1 - n2;     
	  T[0][n1][n2][n3] = A[0][n1][n2][n3] + A[1][n1][n2][n3];
	  T[1][n1][n2][n3] = k2*A[0][n1][n2][n3] + k3*A[1][n1][n2][n3];
	}
    }
}


/* Computes the two products [(Phi_1 + Phi_2)*Phi_3]_order */
/* and [(Phi_1 + Phi_2)*(k2*Phi_1 + k3*Phi_2)]_order      */
void Multiply()     
{
  int counter1, counter2;    /* Basic counters */ 
  int index1, index2;        /* Basic counters */
  int level1, level2;        /* Basic counters */
  register int n1, n2, n3;   /* The exponents of the first factor */
  register int m1, m2, m3;   /* The exponents of the second factor */

  for (level1 = 1; level1 <= order; level1++)              /* level1 + level2 = order + 1 */ 
    {                                                      /* 1 <= level1,level2 <= order */
      for (counter1 = 0; counter1 <= level1; counter1++)
	{
	  n1 = level1 - counter1;                          /* |n| = level1 */
	  for (index1 = 0; index1 <= counter1; index1++)
	    {
	      n2 = level1 - n1 - index1;
	      n3 = level1 - n1 - n2;             	     
	      level2 = order + 1 - level1;	
	      for (counter2 = 0; counter2 <= level2; counter2++)
		{
		  m1 = level2 - counter2;	           /* |m| = level2 */	  
		  for (index2 = 0; index2 <= counter2; index2++)
		    {
		      m2 = level2 - m1 - index2;
		      m3 = level2 - m1 - m2;	   
		      M[0][n1 + m1][n2 + m2][n3 + m3] += T[0][n1][n2][n3]*A[2][m1][m2][m3];
		      M[1][n1 + m1][n2 + m2][n3 + m3] += T[0][n1][n2][n3]*T[1][m1][m2][m3];
		    }
		}
	    }
	}
    }
}


/* Updates the coefficients from */
/* degree order to order + 1     */
void Update()    
{
  int counter, index;            /* Basic counters */
  register int n1, n2, n3;       /* The exponents */
  INTERVAL predivisor;           /* n scalar lambda */ 
  INTERVAL b_n;                  /* temporary storage */

  for (counter = 0; counter <= order + 1; counter++)
    {
      n1 = order + 1 - counter;
      for (index = 0; index <= counter; index++)
	{
	  n2 = order + 1 - n1 - index;
	  n3 = order + 1 - n1 - n2;
	  if ( (n1 < Smoothness) || (n2 + n3 < Smoothness) )  /* The filter */ 
	    {
	      predivisor = n1*lu + n2*lss + n3*ls;
	      A[0][n1][n2][n3] = - k1*M[0][n1][n2][n3]/(predivisor - lu);
	      A[1][n1][n2][n3] = k1*M[0][n1][n2][n3]/(predivisor - lss);
	      A[2][n1][n2][n3] = M[1][n1][n2][n3]/(predivisor - ls);
	    }
	  else
	    {
	      A[0][n1][n2][n3] = Hull(0.0);
	      A[1][n1][n2][n3] = Hull(0.0);
	      A[2][n1][n2][n3] = Hull(0.0);
	    }
	  total++;
	  b_n = Abs(A[0][n1][n2][n3]);          /* b_n = max_i |a_{i,n}| */
	  if ( Inf(b_n) <= Abs(A[1][n1][n2][n3]) )
	    b_n = Abs(A[1][n1][n2][n3]);
	  if ( Inf(b_n) <= Abs(A[2][n1][n2][n3]) )
	    b_n = Abs(A[2][n1][n2][n3]);
	  C[n1 + n2 + n3] += b_n;               /* c_k = sum_{|n|=k} b_n */
	}
    }
}


int main()
{
  INTERVAL sum1;
  int i;
  int linebreak;     /* Counter for the print-out */

  linebreak = 0;
  Init();
  cout.precision(16);
  cout << endl << "Parameters:" << endl;
  cout << "S = " << S << Diam(S) << endl;
  cout << "R = " << R << Diam(R) << endl;
  cout << "B = " << B << Diam(B) << endl;
  cout << endl << "Coefficients:" << endl;
  cout << "lu  = " << lu << endl;
  cout << "lss = " << lss << endl;
  cout << "ls  = " << ls << endl;
  cout << "k1 = " << k1 << endl;
  cout << "k2 = " << k2 << endl;
  cout << "k3 = " << k3 << endl;
  printf("Maximal order: %d \n", MaxOrder);
  printf("Smoothness: %d \n\n", Smoothness);

  for (order = 1; order < MaxOrder; order++)
    {
      Add();       /* Adds terms of degree order */
      Multiply();  /* Gives products of degree order + 1 */
      Update();	   /* Gives coefficient of degree order + 1 */
    }
  for (order = 1; order <= MaxOrder; order++)
    {
      printf("C[%d] = ",order);
      cout << C[order] << endl;
    }
  if (order > 9)
    {
      printf(" \n\n");
      sum1 = Hull(0.0);
      for( i = 1; i < 11; i++ )
	sum1 += C[11-i]*pow(9.0/5,11-i);
      printf("The sum C[1]*(9/5)^1 +...+ C[10]*(9/5)^10 equals ");
      cout << sum1 << endl;
    }
  if (order > 19)
    {
      printf(" \n");
      sum1 = Hull(0.0);
      for( i = 1; i < 20; i++ )
	sum1 += C[i]*C[20-i];
      printf("The sum C[1]*C[19] +...+ C[19]*C[1] equals ");
      cout << sum1 << endl;
    }
  printf("\n");
  printf("Total number of computed coefficients: %ld \n", 3*total);
  printf("**************************************************************\n\n");
  printf("%c", 7);       /* Signal when done */

  return 0;
}

