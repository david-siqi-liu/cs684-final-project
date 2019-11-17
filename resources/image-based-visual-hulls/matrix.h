#define SMALL 1e-20		/* Small enough to be considered zero */
//#define SWAP(A,B,TEMP)    (TEMP=(A), (A)=(B), (B)=(TEMP))
#define ABS(A)		((A) >  0   ? (A) : -(A))
#define MIN(A,B)	((A) <  (B) ? (A) :  (B))
#define MAX(A,B)	((A) >  (B) ? (A) :  (B))

#ifndef TRUE
#define TRUE  (0==0)
#endif
#ifndef FALSE
#define FALSE (0==1)
#endif

#include "stdafx.h"

int InverseMatrix3x3(float* M, float* invM );
void compute_skew_symmetric_matrix(const float* const V/*[3]*/, float* M/*[9]*/);
void MultiplyMatrix3x3(const float* const A, const float* const B, float* M);
void ScalarMultiplyMatrix3x3(const float scalar, const float* const A, float* M);
void MultiplyVector(const float* const A,const float* const B, float* M);
void cross(const float* const A,const  float* const B, float* C);
float determinant(float* M);


static int MAT3_invert3(
	register double source[3][3],	/* A 3 x 3 matrix to be inverted */
	register double inv[3][3] 		/* The inverse of the source matrix */   );
