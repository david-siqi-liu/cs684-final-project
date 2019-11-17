/*********************************************************************************************

  Copyright © 1998-2000 Massachusetts Institute of Technology 
  All Rights Reserved 

*********************************************************************************************/

#include "matrix.h"

int InverseMatrix3x3(float* M, float* invM )
{
        float m11 = M[4]*M[8] - M[5]*M[7];
        float m12 = M[2]*M[7] - M[1]*M[8];
        float m13 = M[1]*M[5] - M[2]*M[4];
        float d = M[0]*m11 + M[3]*m12 + M[6]*m13;
        if (d == 0)
		{
			
            return FALSE;
        
		}
		else 
		{
            d = 1/d;
            invM[0] = d*m11;
			invM[1] = d*m12;
			invM[2] = d*m13;
			invM[3] = d*(M[5]*M[6] - M[3]*M[8]);
			invM[4] = d*(M[0]*M[8] - M[2]*M[6]);
            invM[5] = d*(M[2]*M[3] - M[0]*M[5]);
			invM[6] = d*(M[3]*M[7] - M[4]*M[6]);
            invM[7] = d*(M[1]*M[6] - M[0]*M[7]);
			invM[8] = d*(M[0]*M[4] - M[1]*M[3]);
			return TRUE;
		}
}


///////////////////////////////////////////////////
// Compute Skew Symmetric Matrix 
///////////////////////////////////////////////////
//		M*V = e x V
//		[  0	-v_z	 v_y ]
//		[  v_z	0		-v_x ]
//		[ -v_y	v_x		0	 ]
///////////////////////////////////////////////////		 
void compute_skew_symmetric_matrix(const float* const V/*[3]*/, float* M/*[9]*/)
{
	float v_x = V[0];
	float v_y = V[1];
	float v_z = V[2];

	M[0] = 0;
	M[1] = -v_z;
	M[2] = v_y;
	M[3] = v_z;
	M[4] = 0;
	M[5] = -v_x;
	M[6] = -v_y;
	M[7] = v_x;
	M[8] = 0;
}

////////////////////////////////////////////////////////
//
//		M = A*B
//		A -- 3x3 matrix  
//		B -- 3x3 matrix
//
////////////////////////////////////////////////////////

void MultiplyMatrix3x3(const float* const A, const float* const B, float* M)
{
	M[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
	M[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
	M[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8];
	M[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
	M[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
	M[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8];
	M[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6];
	M[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7];
	M[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
}


////////////////////////////////////////////////////////
//
//		M = A*B
//		A - 3x3 matrix 
//		B - 3d vector
////////////////////////////////////////////////////////

void MultiplyVector(const float* const A, const float* const B, float* M)
{
	M[0] = A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
	M[1] = A[3]*B[0] + A[4]*B[1] + A[5]*B[2];
	M[2] = A[6]*B[0] + A[7]*B[1] + A[8]*B[2];
}

////////////////////////////////////////////////////////
//
//		M = scalar*A
//		A -- 3x3 matrix  
//
////////////////////////////////////////////////////////
void ScalarMultiplyMatrix3x3(float scalar, float* A, float* M)
{
	M[0] = scalar*A[0];
	M[1] = scalar*A[1];
	M[2] = scalar*A[2];
	M[3] = scalar*A[3];
	M[4] = scalar*A[4];
	M[5] = scalar*A[5];
	M[6] = scalar*A[6];
	M[7] = scalar*A[7];
	M[8] = scalar*A[8];

}

////////////////////////////////////////////////////////
//
//		C = AxB
//		A,B,C -- 3d Vectors  
//
////////////////////////////////////////////////////////
void cross(float* A, float* B, float* C)
{
	C[0]= A[1]*B[2] - A[2]*B[1];
	C[1]= A[2]*B[0] - A[0]*B[2];
	C[2]= A[0]*B[1] - A[1]*B[0];
}


float determinant(float* M)
{
    return (M[0]*(M[4]*M[8] - M[5]*M[7]) +
            M[1]*(M[5]*M[6] - M[3]*M[8]) +
            M[2]*(M[3]*M[7] - M[4]*M[6]));
}







