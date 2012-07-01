/*   File: vector_field.cc (BIAS)

     Contains the vector field information
 
     Latest edit: Fri Mar 31 2000
*/

#include "vector_field.h"

////////////////////////////////////////////////////////////////////

static void Naive_Vf_Range(BOX      &, const BOX &);

static void Naive_Vf_Range(INTERVAL &, const BOX &, const short &);

static void DVf_Range     (BOX      &, const BOX &, const short &);

////////////////////////////////////////////////////////////////////

// Returns a VECTOR containing the range of 
// the 'lorenz' vector field w.r.t. the BOX bx.           
void NR_Vf_Range(VECTOR &result, const VECTOR &v)
{
  double C1 = v(1) + v(2);
  double C2 = K1 * C1 * v(3);

  result(1) = E1 * v(1) - C2;
  result(2) = E2 * v(2) + C2;
  result(3) = E3 * v(3) + C1 * (K2 * v(1) + K3 * v(2));
}

////////////////////////////////////////////////////////////////////

// Returns a BOX containing the range of 
// the 'lorenz' vector field w.r.t. the BOX bx.           
static void Naive_Vf_Range(BOX &result, const BOX &bx)
{
  INTERVAL C1 = bx(1) + bx(2);
  INTERVAL C2 = K1_IV * C1 * bx(3);

  result(1) = E1_IV * bx(1) - C2;
  result(2) = E2_IV * bx(2) + C2;
  result(3) = E3_IV * bx(3) + C1 * (K2_IV * bx(1) + K3_IV * bx(2));
}

////////////////////////////////////////////////////////////////////

// Its MVT approach
BOX Vf_Range(const BOX &bx)
{
  BOX center(DIM);
  BOX symmrad(DIM);
  BOX vf_box(DIM);
  INTERVAL_MATRIX DVf(DIM, DIM);

  Mid_And_SymRad(center, symmrad, bx);
  Naive_Vf_Range(vf_box, center);
  DVf_Range(DVf, bx);

  return vf_box + DVf * symmrad;
}

////////////////////////////////////////////////////////////////////

void Vf_Range(BOX &result, const BOX &bx)
{
  static BOX center(DIM);
  static BOX symmrad(DIM);
  static BOX vf_box(DIM);
  static INTERVAL_MATRIX DVf(DIM, DIM);

  Mid_And_SymRad(center, symmrad, bx);
  Naive_Vf_Range(vf_box, center);
  DVf_Range(DVf, bx);

  result = vf_box + DVf * symmrad;
}

////////////////////////////////////////////////////////////////////

// Returns an INTERVAL containing the range of 
// the i:th component of the 'lorenz' vector field 
// w.r.t. the BOX bx.
static void Naive_Vf_Range(INTERVAL &result, const BOX &bx, const short &i)
{
  if ( i == 1 )
    result = E1_IV * bx(1) - K1_IV * (bx(1) + bx(2)) * bx(3);
  else if ( i == 2 )
    result = E2_IV * bx(2) + K1_IV * (bx(1) + bx(2)) * bx(3);
  else if ( i == 3 )
    result = E3_IV * bx(3) + (bx(1) + bx(2)) * (K2_IV * bx(1) + K3_IV * bx(2));
  else
    {
      char *msg = "Error: 'Vf_Range'. Index i out of range.";
      throw Error_Handler(msg);
    }
}

////////////////////////////////////////////////////////////////////

// Its MVT approach
INTERVAL Vf_Range(const BOX &bx, const short &i)
{
  INTERVAL vf_center;
  BOX center(DIM);
  BOX symmrad(DIM);
  BOX DVf_row(DIM);

  Mid_And_SymRad(center, symmrad, bx);
  Naive_Vf_Range(vf_center, center, i);
  DVf_Range(DVf_row, bx, i);

  return vf_center + DVf_row * symmrad;
}

////////////////////////////////////////////////////////////////////

void Vf_Range(INTERVAL &result, const BOX &bx, const short &i)
{
  INTERVAL vf_center;
  BOX center(DIM);
  BOX symmrad(DIM);
  BOX DVf_row(DIM);

  Mid_And_SymRad(center, symmrad, bx);
  Naive_Vf_Range(vf_center, center, i);
  DVf_Range(DVf_row, bx, i);

  result = vf_center + DVf_row * symmrad;
}

////////////////////////////////////////////////////////////////////

// Returns the partial derivatives of
// the vector field Vf_Range wrt the BOX bx
void DVf_Range(INTERVAL_MATRIX &result, const BOX &bx)
{
  INTERVAL C1 = K1_IV * (bx(1) + bx(2));
  INTERVAL C2 = K1_IV * bx(3);

  result(1, 1) = E1_IV - C2;
  result(1, 2) = - C2;
  result(1, 3) = - C1;
  result(2, 1) = + C2;
  result(2, 2) = E2_IV + C2;
  result(2, 3) = + C1;
  result(3, 1) = TWO_K2_IV * bx(1) + K2_PLUS_K3_IV * bx(2);
  result(3, 2) = K2_PLUS_K3_IV * bx(1) + TWO_K3_IV * bx(2);
  result(3, 3) = E3_IV;
}

////////////////////////////////////////////////////////////////////

// Returns the i:th row of the partial derivatives
// of the vector field Vf_Range wrt the BOX bx
static void DVf_Range(BOX &result, const BOX &bx, const short &i)
{
  INTERVAL C1 = K1_IV * (bx(1) + bx(2));
  INTERVAL C2 = K1_IV * bx(3);

  if ( i == 1 )
    {
      result(1) = E1_IV - C2;
      result(2) = - C2;
      result(3) = - C1;
      return;
    }
  if ( i == 2)
    {
      result(1) = + C2;
      result(2) = E2_IV + C2;
      result(3) = + C1;
      return;
    }
  if ( i == 3 )
    { 
      result(1) = TWO_K2_IV * bx(1) + K2_PLUS_K3_IV * bx(2);
      result(2) = K2_PLUS_K3_IV * bx(1) + TWO_K3_IV * bx(2);
      result(3) = E3_IV;
      return;
    }
  char *msg = "Error: 'DVf_Range'. Index row out of range.";
  throw Error_Handler(msg);
}

////////////////////////////////////////////////////////////////////

// Returns the (i,j)-partial derivative of
// the vector field Vf_Range wrt the BOX bx
INTERVAL DVf_Range(const BOX &bx, const short &i, const short &j)
{
  if (i == 1)
    {
      if (j == 1)
	return E1_IV - K1_IV * bx(3);
      if (j == 2)
	return - K1_IV * bx(3);
      if (j == 3)
	return - K1_IV * (bx(1) + bx(2));
    }
  if (i == 2)
    {
      if (j == 1)
	return + K1_IV * bx(3);
      if (j == 2)
	return E2_IV + K1_IV * bx(3);
      if (j == 3)
	return + K1_IV * (bx(1) + bx(2));
    }
  if (i == 3)
    {
      if (j == 1)
	return TWO_K2_IV * bx(1) + (K2_PLUS_K3_IV) * bx(2);
      if (j == 2)
	return (K2_PLUS_K3_IV) * bx(1) + TWO_K3_IV * bx(2);
      if (j == 3)
	return E3_IV;
    }
  char *msg = "Error: 'DVf_Range'. Index (i,j) out of range.";
  throw Error_Handler(msg);
  return ZERO_IV;          // This line keeps the compiler with "-Wall" happy.
}

////////////////////////////////////////////////////////////////////

void DVf_Range(INTERVAL &result, const BOX &bx, const short &i, const short &j)
{
  if (i == 1)
    {
      if (j == 1)
	result = E1_IV - K1_IV * bx(3);
      else if (j == 2)
	result = - K1_IV * bx(3);
      else if (j == 3)
	result = - K1_IV * (bx(1) + bx(2));
    }
  else if (i == 2)
    {
      if (j == 1)
	result = + K1_IV * bx(3);
      else if (j == 2)
	result = E2_IV + K1_IV * bx(3);
      else if (j == 3)
	result = + K1_IV * (bx(1) + bx(2));
    }
  else if (i == 3)
    {
      if (j == 1)
	result = TWO_K2_IV * bx(1) + (K2_PLUS_K3_IV) * bx(2);
      else if (j == 2)
	result = (K2_PLUS_K3_IV) * bx(1) + TWO_K3_IV * bx(2);
      else if (j == 3)
	result = E3_IV;
    }
}

////////////////////////////////////////////////////////////////////


