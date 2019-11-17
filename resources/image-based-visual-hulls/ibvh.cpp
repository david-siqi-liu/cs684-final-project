
#include "stdafx.h"
#include "memory.h"
#include "math.h"
#include "matrix.h"
#include "stdlib.h"
#include "assert.h"
#include "float.h"



//#include <time.h>
//#include <windows.h>
//#include "mmsystem.h"

#define FOREGROUND			0
#define BACKGROUND			255
#define MAX_CAMERAS			16

#define DES_IMAGE_WIDTH		256
#define DES_IMAGE_HEIGHT	256
#define REF_IMAGE_WIDTH		256
#define REF_IMAGE_HEIGHT	256

#define DX_SAMPLE_RATE		4
#define DY_SAMPLE_RATE		4
#define INFINITY				1000

#define DEPTH_THRESH		2.5

// comment out the next line to get normal image
#define OUTPUT_DEPTH			0 

// comment out the next line to get sparse image
#define INTERPOLATE				0


#define INPUT_TEST_SCENE		"c:\\ibvh\\tri\\Tri%ds.raw"

#define OUTPUT_FILE				"c:\\ibvh\\tri\\test.raw"

///////////////////////////
// TYPES



typedef struct vector3D
{
	float x;
	float y;
	float z;
} VECTOR3D;


// this is ray data structure
typedef struct line
{
	int					count;
	int					maxminCamera; // this keeps the camera number that determines the edge
	float*				t;
	VECTOR3D*			normal;

} LINE; 


typedef struct slopeel
{
	float slope;
	float startX;
	float startY;
	float endX;
	float endY;
	bool  direction;
	int wedge_cache_index;
} SLOPEEL; 

typedef struct wedge
{
	int computed;
	int num_edges;
	int first_edge_off;
	int last_edge_off;
} WEDGE;

typedef struct wedge_cache
{
	int edge_count;
	WEDGE wedges[2*REF_IMAGE_WIDTH+2*REF_IMAGE_HEIGHT];
} WEDGE_CACHE;



typedef struct epipole
{
	int		u;
	int		v;
	float	ufl;
	float	vfl;
	bool	negative;
	bool	positive;
	bool	inside;
} EPIPOLE; 

typedef struct iterOrder
{
	float	cos;
	int		pos;
} ITERORDER;

typedef struct edge
{
	float		x1;
	float		x2;
	float		y1;
	float		y2;

} EDGE; 

typedef struct edgeel
{
	EDGE	e;
	float	distance;
	float	dxinv;
	float	dy_dx;
} EDGEEL; 




//////////////////////////////////
// VARS


WEDGE_CACHE cache[MAX_CAMERAS];
static			EDGEEL edgebuffer_[MAX_CAMERAS][50000];

float	V10[DES_IMAGE_WIDTH*DES_IMAGE_HEIGHT];
float	V11[DES_IMAGE_WIDTH*DES_IMAGE_HEIGHT];
float	V12[DES_IMAGE_WIDTH*DES_IMAGE_HEIGHT];
// Ps holds P matrices for all cameras
float Ps[MAX_CAMERAS][9];
// Cs holds Center of proj for all cameras
float Cs[MAX_CAMERAS][9];

// Silhouettes hold silhouettes for all cameras
unsigned char* Silhouettes[MAX_CAMERAS];

// Number of cameras
int ImageCount = 0;




// this is the output 
LINE** T;

float** FMs;
float** CDifs;
float** RefPInvDesP;
float** RefPInv;
EPIPOLE* E0s_des;
EPIPOLE* E0s_ref;
ITERORDER* iterations;


float*			tBuffer1;
float*			tBuffer2;
VECTOR3D*		normalBuffer1;
VECTOR3D*		normalBuffer2;
int				tnOffset = 0;
int				firstPasstnOffset = 0;
int				max_samples = 1000000;


// helpers
LINE* tempT;
LINE* temp2T;
float*  tB1;
VECTOR3D* nB1;






void setupFrameSet(float** &FMs, float** &CDifs, float** &RefPInvDesP, float** &RefPInv, float* DesP, float* DesC, EPIPOLE* &E0s_des, EPIPOLE* &E0s_ref, ITERORDER* &iterations);
int compute_minimum_epipolar_extent(
			float* M,
			float* T1,
			float* T2,
			int width,
			int height,
			EPIPOLE&  epi);
int compute_fundamental_matrix(
		float* P1 /* 3x3 matrix */,
		float* P2 /* 3x3 matrix */,
		float* T1 /* 3d Vector  */,
		float* T2 /* 3d Vector  */,
		float* FM /* fundamental Matrix */
		);

void clearWedgeCache(WEDGE_CACHE& cache);
void sampleCone(LINE** T, int ref_camera, int des_camera);
void sampleConeAndMerge(LINE** T, int ref_camera, int des_camera, int iteration);
void sampleConeSecondPass(LINE** T, int ref_camera, int des_camera);
void sampleConeAndMergeSecondPass(LINE** T, int ref_camera, int des_camera, int iteration);
bool ComputeSlopeElement(SLOPEEL& element, EPIPOLE& E0, float* Einf);
int compute_bbox_intersections(const EPIPOLE& e0, float dyOverdx, float* interU,float* interV);
int compute_effective_intersection1(float* Interval1Start,float* Interval1End,
								   float* Xo,float* Xinf,
								   float** Start,float** End);
int compute_effective_intersection2(float* Interval1Start,float* Interval1End,
								   float* Xo,float* Xinf,
								   float** Start,float** End);
int compute_effective_intersection3(float* Interval1Start,float* Interval1End,
								   float* Xo,float* Xinf,
								   float** Start,float** End);
void traverse_wedge(int wedge_index, EPIPOLE& e0, unsigned char* sil, WEDGE& w,int image_no);
void traverse_int_line(int x1, int y1, int x2, int y2,unsigned char* silhouette,int image_no);
inline void add_vertical_edge(int ix1, int iy1, unsigned char& val,int image_no);
inline void add_horizontal_edge(int ix1, int iy1, unsigned char& val,int image_no);
int compute_inner_bbox_intersections(float* point, float dyOverdx, float* interU,float* interV);
void switchToOrdered(float* vPoint, float* interU, float* interV);

void computeVHull(int cam_no);
inline bool compute_T(float& t, float* P2, float* C2minC1, float V1_0, float V1_1, float V1_2, float u1, float v1);
void compute3DIntervals(WEDGE& w, SLOPEEL& current_el, EPIPOLE& E0,
						float* RefP, float* RefCminusDesC, unsigned char *RefSilhouette,LINE* tempT,
						float* V1_0,float* V1_1,float* V1_2,int index, int order);

int  compute_wedge_cache_index(float flU, float flV);
int  compute_wedge_cache_index(float U1, float V1, float U2, float V2, float vp_U, float vp_V);
void merge2TLists(const LINE& ray1, const LINE& ray2, LINE& rayout);
void compute_normal(VECTOR3D& n,float* RefP, float x1,float y1,float x2,float y2);
void allocateMemory(int count);
void cleanBuffer();














int main(int argc, char* argv[])
{
	
	for (int i = 0; i < 6; i++)
	{
		char filename[256];
		sprintf(filename,INPUT_TEST_SCENE,i);
		FILE* fp = NULL;
		fp = fopen(filename,"rb");
		Silhouettes[i] = (unsigned char*) malloc(256*256*sizeof(unsigned char));
		fread(Silhouettes[i],sizeof(unsigned char),256*256,fp);
		fclose(fp);
		
	}
	
	Ps[0][0] =  (float)0.284352;
	Ps[0][1] =  (float)0.0;
	Ps[0][2] =  (float)-36.397;
	Ps[0][3] =  (float)0.0;
	Ps[0][4] =  (float)0.0;
	Ps[0][5] =  (float)100.0;
	Ps[0][6] =  (float)0.0;
	Ps[0][7] =  (float)-0.284352;
	Ps[0][8] =  (float)36.397; 
	Cs[0][0] = (float)0.0;
	Cs[0][1] = (float)-100.0;
	Cs[0][2] = (float)0.0; 


	Ps[1][0] =  (float)0.170611;
	Ps[1][1] =  (float)0.0;
	Ps[1][2] =  (float)58.1618;
	Ps[1][3] =  (float)-0.227481;
	Ps[1][4] =  (float)0.0;
	Ps[1][5] =  (float)89.1176;
	Ps[1][6] =  (float)0.0;
	Ps[1][7] =  (float)-0.284352;
	Ps[1][8] =  (float)36.397;
	Cs[1][0] =  (float)-80.0;
	Cs[1][1] =  (float)-60.0;
	Cs[1][2] =  (float)0.0; 

	Ps[2][0] =  (float)0.208955;
	Ps[2][1] =  (float)0.12064;
	Ps[2][2] =  (float)-102.188;
	Ps[2][3] =  (float)0.208955;
	Ps[2][4] =  (float)-0.12064;
	Ps[2][5] =  (float)48.6957;
	Ps[2][6] =  (float)0.0;
	Ps[2][7] =  (float)-0.24128;
	Ps[2][8] =  (float)-29.1161; 
	Cs[2][0] =  (float)60.0;
	Cs[2][1] =  (float)-60.0;
	Cs[2][2] =  (float)60.0;

	Ps[3][0] =  (float)0.242068;
	Ps[3][1] =  (float)-0.0206896;
	Ps[3][2] =  (float)-48.3364;
	Ps[3][3] =  (float)0.0605169;
	Ps[3][4] =  (float)0.0827585;
	Ps[3][5] =  (float)61.6607;
	Ps[3][6] =   (float)0;
	Ps[3][7] =  (float)-0.234482;
	Ps[3][8] =  (float)60.0138; 
	Cs[3][0] =  (float)20;
	Cs[3][1] =  (float)-80;
	Cs[3][2] =  (float)-30; 

	Ps[4][0] =  (float)0.233099;
	Ps[4][1] =  (float)0.0254332;
	Ps[4][2] =  (float)6.90787;
	Ps[4][3] =  (float)-0.11655;
	Ps[4][4] =  (float)0.0508664;
	Ps[4][5] =  (float)88.4074;
	Ps[4][6] =  (float)0;
	Ps[4][7] =  (float)-0.254332;
	Ps[4][8] =  (float)52.5545; 
	Cs[4][0] =  (float) -40;
	Cs[4][1] =  (float)-80;
	Cs[4][2] =  (float)-20; 


	Ps[5][0] =  (float)0.312787;
	Ps[5][1] =  (float)0.000000;
	Ps[5][2] =  (float)-40.036701;
	Ps[5][3] =  (float)0.000000;
	Ps[5][4] =  (float)-0.312787;
	Ps[5][5] =  (float)40.036701;
	Ps[5][6] =  (float)0.000000;
	Ps[5][7] =  (float)0.000000;
	Ps[5][8] =  (float)-110.000000; 
	Cs[5][0] =  (float)0;
	Cs[5][1] =  (float)0;
	Cs[5][2] =  (float)110; 



	ImageCount = 6;
	allocateMemory(ImageCount);
	cleanBuffer();


	float DesP[9];
	float DesC[3];

	int cam_no = 0;

	DesP[0] = Ps[cam_no][0];
	DesP[1] = Ps[cam_no][1];
	DesP[2] = Ps[cam_no][2];
	DesP[3] = Ps[cam_no][3];
	DesP[4] = Ps[cam_no][4];
	DesP[5] = Ps[cam_no][5];
	DesP[6] = Ps[cam_no][6];
	DesP[7] = Ps[cam_no][7];
	DesP[8] = Ps[cam_no][8];
	

	DesC[0] = Cs[cam_no][0];
	DesC[1] = Cs[cam_no][1];
	DesC[2] = Cs[cam_no][2];

	setupFrameSet(FMs, CDifs, RefPInvDesP,RefPInv, DesP, DesC,E0s_des, E0s_ref,iterations);
	

//	DWORD StartTime1 = timeGetTime();
	
//	for (int it = 0; it < 100; it++)
	{
		cleanBuffer();
		computeVHull(cam_no);
	}
//	DWORD EndTime1 = timeGetTime();
//	float time_per_frame = EndTime1-StartTime1;
//	time_per_frame =  time_per_frame / 100.0;
//	printf("time per frame: %f , fps: %f \n",time_per_frame/1000., 1000./time_per_frame);


	float minDepth = +100000.0f;
	float maxDepth = -100000.0f;


	int src_offset = 0;
	int dst_offset = 0;
	int x,y;

	for (y = 0; y < DES_IMAGE_HEIGHT; y++)
	{
		for (x = 0; x < DES_IMAGE_WIDTH; x++)
		{
			if (T[src_offset]->count > 0)
			{
				float depth = T[src_offset]->t[0];
				if (depth < minDepth) 
					minDepth = depth;
				if (depth > maxDepth) 
					maxDepth = depth;
			}
			src_offset++;
		}
	}


	float light_dir[3];

	light_dir[0]  = -DesP[0]*(DES_IMAGE_WIDTH/2.) - DesP[1]*(DES_IMAGE_HEIGHT/2.) - DesP[2];
	light_dir[1]  = -DesP[3]*(DES_IMAGE_WIDTH/2.) - DesP[4]*(DES_IMAGE_HEIGHT/2.) - DesP[5];
	light_dir[2]  = -DesP[6]*(DES_IMAGE_WIDTH/2.) - DesP[7]*(DES_IMAGE_HEIGHT/2.) - DesP[8];


	float len = light_dir[0]*light_dir[0]+ light_dir[1]*light_dir[1]+light_dir[2]*light_dir[2];
	len = sqrt(len);
	light_dir[0] = light_dir[0]/len;
	light_dir[1] = light_dir[1]/len;
	light_dir[2] = light_dir[2]/len;

	unsigned char* frameptr_ = (unsigned char*) malloc (DES_IMAGE_HEIGHT* DES_IMAGE_WIDTH);
	memset(frameptr_,0,256*256);
	src_offset = 0;
	float depthDiff = maxDepth - minDepth;
	for ( y = 0; y < DES_IMAGE_HEIGHT-4; y++)
	{
		dst_offset = DES_IMAGE_WIDTH*y;
		src_offset = DES_IMAGE_WIDTH*y;

		for (x = 0; x < DES_IMAGE_WIDTH-4; x++)
		{

			int ix = x - (x % 4);
			int iy = y - (y % 4);
			int off = ix + DES_IMAGE_WIDTH*iy;
			int off_right = 4 + ix + DES_IMAGE_WIDTH*iy;
			int off_down = ix + DES_IMAGE_WIDTH*(iy+4);
			int off_down_right = 4 + ix + DES_IMAGE_WIDTH*(iy+4);

#ifdef INTERPOLATE	

			if (T[off]->count > 0 && T[off_right]->count > 0 && T[off_down]->count > 0  && T[off_down_right]->count > 0)
			{
				if (T[src_offset]->count == -1)
				{
#ifdef OUTPUT_DEPTH	
					float t1 = T[off]->t[0];
					float t2 = T[off_right]->t[0];
					float t3 = T[off_down]->t[0];
					float t4 = T[off_down_right]->t[0];
					float du = ((float)(x % 4))/4.0f;
					float dv = ((float)(y % 4))/4.0f;
					
					float Pu = (1-du)*t1 + du*t2;
					float Pd = (1-du)*t3 + du*t4;
					float t = Pu*(1-dv) +Pd*dv;
					
					float depth = (t-minDepth)/(maxDepth-minDepth);
					depth = depth*255.0f;
					if (depth < 0) depth = 0;
					if (depth > 255) depth = 255;
					frameptr_[dst_offset] = 255-(int)depth;



#else
					float ray_x  = DesP[0]*x + DesP[1]*y + DesP[2];
					float ray_y  = DesP[3]*x + DesP[4]*y + DesP[5];
					float ray_z  = DesP[6]*x + DesP[7]*y + DesP[8];
					
					float ny = T[off]->normal[0].y;
					float nx = T[off]->normal[0].x;
					float nz = T[off]->normal[0].z;

					float len = sqrt(nx*nx + ny*ny + nz*nz);
					nx = nx / len;
					ny = ny / len;
					nz = nz / len;

					float cosA = ray_x*nx+ray_y*ny+ray_z*nz;
					if (cosA > 0)
					{
						nx = -nx;
						ny = -ny;
						nz = -nz;
					}

					cosA  = nx*light_dir[0] +ny*light_dir[1]+nz*light_dir[2];
					cosA = cosA*255.0f;
					if (cosA < 0) cosA = 0;
					if (cosA > 255) cosA = 255;
					frameptr_[dst_offset] = (int)cosA;
#endif


					dst_offset++;
					src_offset++;
					continue;
				}
			}
#endif
			if (T[src_offset]->count != -1)
			{
				if (T[src_offset]->count > 0)
				{
#ifdef OUTPUT_DEPTH	
					float depth = (T[src_offset]->t[0]-minDepth)/(maxDepth-minDepth);
					depth = depth*255.0f;
					if (depth < 0) depth = 0;
					if (depth > 255) depth = 255;
					frameptr_[dst_offset] = 255-(int)depth;
#else
				
					float ray_x  = DesP[0]*x + DesP[1]*y + DesP[2];
					float ray_y  = DesP[3]*x + DesP[4]*y + DesP[5];
					float ray_z  = DesP[6]*x + DesP[7]*y + DesP[8];


					float nx = T[src_offset]->normal[0].x;
					float ny = T[src_offset]->normal[0].y;
					float nz = T[src_offset]->normal[0].z;

					float len = sqrt(nx*nx + ny*ny + nz*nz);
					nx = nx / len;
					ny = ny / len;
					nz = nz / len;

					float cosA = ray_x*nx+ray_y*ny+ray_z*nz;
					if (cosA > 0)
					{
						nx = -nx;
						ny = -ny;
						nz = -nz;
					}

					cosA  = nx*light_dir[0] +ny*light_dir[1]+nz*light_dir[2];
					cosA = cosA*255.0f;
					if (cosA < 0) cosA = 0;
					if (cosA > 255) cosA = 255;
					frameptr_[dst_offset] = (int)cosA;
#endif

				}
				else
				{
					frameptr_[dst_offset] = 0;
				}
			}
			else
			{
				frameptr_[dst_offset] = 0;
			}
			dst_offset++;
			src_offset++;
		}
	}




	FILE* fp = NULL;
	fp = fopen(OUTPUT_FILE,"wb");
	fwrite(frameptr_,sizeof(unsigned char),256*256,fp);
	fclose(fp);

	return 0;
}





void computeVHull(int cam_no) // you compute the visual hull for camera_no camera
{

	int iter;
	unsigned char* 	RefSilhouette;
	unsigned char* 	DesSilhouette;

	int Count = ImageCount;
	int Conesprocessed;

	float* RefP;
	float* RefC;
	float* RefCminusDesC;
	float* FM;
	float* TempMatrix;


	EPIPOLE E0;

	// FIRST PASS
	Conesprocessed = 0;
	for (iter = 0; iter < Count; iter++)
	{
		int order = iterations[iter].pos;

		if (order == cam_no) continue;

		clearWedgeCache(cache[order]);
		
		if (Conesprocessed == 0)
			sampleCone(T,order, cam_no);
		else
			sampleConeAndMerge(T,order, cam_no,Conesprocessed);
		Conesprocessed++;
	}

	// SECOND PASS
	firstPasstnOffset = tnOffset;
	Conesprocessed = 0;
	for (iter = 0; iter < Count; iter++)
	{
		int order = iterations[iter].pos;

		if (order == cam_no) continue;

		RefC =  Cs[order];
		RefP =  Ps[order];
		RefSilhouette = Silhouettes[order];
		DesSilhouette = Silhouettes[cam_no];
		RefCminusDesC = CDifs[order];
		FM = FMs[order];
		TempMatrix = RefPInvDesP[order];
		E0 = E0s_des[order];	
		
		if (Conesprocessed == 0)
			sampleConeSecondPass(T,order, cam_no);
		else
			sampleConeAndMergeSecondPass(T,order, cam_no,Conesprocessed);
		Conesprocessed++;
	}
	
}






int compare( const void *arg1, const void *arg2 )
{
	ITERORDER el1 = *((ITERORDER* ) arg1);
	ITERORDER el2 = *((ITERORDER* ) arg2);

	if (el1.cos < el2.cos)
		return 1;
	else
	if (el1.cos > el2.cos)
		return -1;
	return 0;
}


void setupFrameSet( float** &FMs, float** &CDifs, float** &RefPInvDesP, float** &RefPInv, float* DesP, float* DesC, EPIPOLE* &E0s_des, EPIPOLE* &E0s_ref, ITERORDER* &iterations)
{
	float* P;
	float* C;
	int Count = ImageCount;

	
	int index1;
	int width;
	int height;
	// compute look vector for Desired
	float DesLook[3];
	DesLook[0] = DesP[0]*DES_IMAGE_WIDTH/2 + DesP[1]*DES_IMAGE_HEIGHT/2 +DesP[2];
	DesLook[1] = DesP[3]*DES_IMAGE_WIDTH/2 + DesP[4]*DES_IMAGE_HEIGHT/2 +DesP[5];
	DesLook[2] = DesP[6]*DES_IMAGE_WIDTH/2 + DesP[7]*DES_IMAGE_HEIGHT/2 +DesP[8];
	float DesLookLen = (DesLook[0]*DesLook[0] + DesLook[1]*DesLook[1] + DesLook[2]*DesLook[2]);
	DesLookLen = (float)sqrt(DesLookLen);
	DesLook[0] = DesLook[0] / DesLookLen;
	DesLook[1] = DesLook[1] / DesLookLen;
	DesLook[2] = DesLook[2] / DesLookLen;

	float RefLook[3];

	for (index1 = 0; index1 < Count; index1++)
	{
		P = Ps[index1];
		C = Cs[index1];
		width = REF_IMAGE_WIDTH;
		height = REF_IMAGE_HEIGHT;
		// epipole of desired on reference image
		compute_minimum_epipolar_extent(/* ref */ P,
										/* main */ DesC,
										/* ref */ C,
										REF_IMAGE_WIDTH,
										REF_IMAGE_HEIGHT,
										E0s_des[index1]);
		// epipole of reference on desired image
		compute_minimum_epipolar_extent(/* desired */ DesP,
										/* reference */ C,
										/* desired */ DesC,
										DES_IMAGE_WIDTH,
										DES_IMAGE_HEIGHT,
										E0s_ref[index1]);

		compute_fundamental_matrix( /* main */ DesP , /* ref */ P ,
									/* main */ DesC , /* ref */ C ,
   									FMs[index1] );
		CDifs[index1][0] = C[0]- DesC[0];
		CDifs[index1][1] = C[1]- DesC[1];
		CDifs[index1][2] = C[2]- DesC[2];
		if (InverseMatrix3x3(P, RefPInv[index1] ) == false)
			exit(4);
		MultiplyMatrix3x3(RefPInv[index1],DesP, RefPInvDesP[index1]);


		RefLook[0] = P[0]*width/2 + P[1]*height/2 +P[2];
		RefLook[1] = P[3]*width/2 + P[4]*height/2 +P[5];
		RefLook[2] = P[6]*width/2 + P[7]*height/2 +P[8];
		float RefLookLen = (RefLook[0]*RefLook[0] + RefLook[1]*RefLook[1] + RefLook[2]*RefLook[2]);
		RefLookLen = (float)sqrt(RefLookLen);
		RefLook[0] = RefLook[0] / RefLookLen;
		RefLook[1] = RefLook[1] / RefLookLen;
		RefLook[2] = RefLook[2] / RefLookLen;

		float cosRefLDesL = DesLook[0]*RefLook[0] + DesLook[1]*RefLook[1] + DesLook[2]*RefLook[2];
		cosRefLDesL = ABS(cosRefLDesL);

		iterations[index1].cos = cosRefLDesL;
		iterations[index1].pos = index1;
	}

	qsort( iterations, Count, sizeof( ITERORDER), compare  );

}

/////////////////////////////////////////////////
//
// Compute Minimum Extent of Epipolar Line Segment
//		Xo = P2inv*(T1-T2)
//
/////////////////////////////////////////////////
int compute_minimum_epipolar_extent(
			float* M,
			float* T1,
			float* T2,
			int width,
			int height,
			EPIPOLE&  epi)
{
	double Tdiff[3];
	double invM[9];
	double e0[3];

	// Compute P2 Inverse
//	if (InverseMatrix3x3(P2, P2inv ) == false)
//		return false;

    double m11 = M[4]*M[8] - M[5]*M[7];
    double m12 = M[2]*M[7] - M[1]*M[8];
    double m13 = M[1]*M[5] - M[2]*M[4];
    double d = M[0]*m11 + M[3]*m12 + M[6]*m13;

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


	// Compute T1 - T2
	Tdiff[0] = T1[0] - T2[0];
	Tdiff[1] = T1[1] - T2[1];
	Tdiff[2] = T1[2] - T2[2];

	// Compute Minimum Extent
	//MultiplyVector(P2inv, Tdiff,e0);

	e0[0] = invM[0]*Tdiff[0] + invM[1]*Tdiff[1] + invM[2]*Tdiff[2];
	e0[1] = invM[3]*Tdiff[0] + invM[4]*Tdiff[1] + invM[5]*Tdiff[2];
	e0[2] = invM[6]*Tdiff[0] + invM[7]*Tdiff[1] + invM[8]*Tdiff[2];


	epi.ufl = (float)(e0[0]/e0[2]);
	epi.u = (int)(epi.ufl+0.5f);

	epi.vfl = (float)(e0[1]/e0[2]);
	epi.v = (int)(epi.vfl+0.5f);

	if (e0[2] <= 0)
	{
		epi.negative = true;
		epi.positive = false;
	}
	else
	{
		epi.negative = false;
		epi.positive = true;
	}

	if  ( epi.u >= 0 && epi.v >= 0 &&
		  epi.u <= width-1 && epi.v <= height-1 )
	{
		epi.inside = true;
	}
	else
	{
		epi.inside = false;
	}
	return true;

}


// compute fundamental matrix
int compute_fundamental_matrix(
		float* P1 /* 3x3 matrix */,
		float* P2 /* 3x3 matrix */,
		float* T1 /* 3d Vector  */,
		float* T2 /* 3d Vector  */,
		float* FM /* fundamental Matrix */
		)
{
	float	P2inv[9];
	float	Tdiff[3];
	float	e1[3];
	float	ssM[9];
	float   tempM[9];

	// Compute P2 Inverse
	if (InverseMatrix3x3(P2, P2inv ) == false)
		return false;

	// Compute T1 - T2
	Tdiff[0] = T1[0] - T2[0];
	Tdiff[1] = T1[1] - T2[1];
	Tdiff[2] = T1[2] - T2[2];

	// Compute P2Inv*(T1-T2)
	MultiplyVector(P2inv, Tdiff,e1);

	// Compute Skew Symmetric Matrix 
	compute_skew_symmetric_matrix(e1, ssM);

	// Compute ssM*P2inv*P1
	MultiplyMatrix3x3(ssM, P2inv, tempM);
	MultiplyMatrix3x3(tempM,P1, FM);

	return true;
}

void clearWedgeCache(WEDGE_CACHE& cache)
{
	cache.edge_count = 0;
	for (int i = 0; i < 2*REF_IMAGE_WIDTH+2*REF_IMAGE_HEIGHT; i++)
	{
		cache.wedges[i].computed = 0;
	}
}



void sampleCone(LINE** T, int ref_camera, int des_camera)
{
	unsigned char* 	RefSilhouette = Silhouettes[ref_camera];
	unsigned char* 	DesSilhouette = Silhouettes[des_camera];

	float DesP[9];
	float DesC[3];


	DesP[0] = Ps[des_camera][0];
	DesP[1] = Ps[des_camera][1];
	DesP[2] = Ps[des_camera][2];
	DesP[3] = Ps[des_camera][3];
	DesP[4] = Ps[des_camera][4];
	DesP[5] = Ps[des_camera][5];
	DesP[6] = Ps[des_camera][6];
	DesP[7] = Ps[des_camera][7];
	DesP[8] = Ps[des_camera][8];
	

	DesC[0] = Cs[des_camera][0];
	DesC[1] = Cs[des_camera][1];
	DesC[2] = Cs[des_camera][2];


	float* RefP = Ps[ref_camera];
	float* RefC = Cs[ref_camera];
	
	float* RefCminusDesC = CDifs[ref_camera];
	float* FM = FMs[ref_camera];
	float* TempMatrix = RefPInvDesP[ref_camera];

	EPIPOLE E0 = E0s_des[ref_camera];	


	int index;
	float* tBuffer;
	vector3D* nBuffer;
	float flv, flu;
	tnOffset = 0;
	SLOPEEL current_el;
	int wedge_cache_index;
	float	scanline_infX[3];
	float	infX[3];

	
	tBuffer = tBuffer1;
	nBuffer = normalBuffer1;
	
	// first pass
	for (int v = 0; v < DES_IMAGE_HEIGHT; v+=DY_SAMPLE_RATE)
	{
		flv = (float)v;
		
		int scanline_index = v*DES_IMAGE_WIDTH;
		index = scanline_index;
		scanline_infX[0] = TempMatrix[1]*flv + TempMatrix[2];
		scanline_infX[1] = TempMatrix[4]*flv + TempMatrix[5];
		scanline_infX[2] = TempMatrix[7]*flv + TempMatrix[8];
	
		for (int u = 0; u < DES_IMAGE_WIDTH; u+=DX_SAMPLE_RATE)
		{

			T[index]->count = 0;
			T[index]->maxminCamera = ref_camera;

			if (DesSilhouette[index] == BACKGROUND)
			{
				index +=DX_SAMPLE_RATE;
				continue;
			}

			flu = (float)u;
			infX[2] = TempMatrix[6]*flu + scanline_infX[2];
			if ((E0.negative == true) && (infX[2] <= 0))
			{
				index +=DX_SAMPLE_RATE;
				continue;
			}
		
			infX[0] = TempMatrix[0]*flu + scanline_infX[0];
			infX[1] = TempMatrix[3]*flu + scanline_infX[1];
			if (ComputeSlopeElement(current_el,  E0, infX) == false)
			{			
				index +=DX_SAMPLE_RATE;
				continue;
			}
			
			wedge_cache_index = current_el.wedge_cache_index;
			if (cache[ref_camera].wedges[wedge_cache_index].computed == 0)
				traverse_wedge(wedge_cache_index,E0,RefSilhouette, cache[ref_camera].wedges[wedge_cache_index],ref_camera);
						
			
			T[index]->count = 0;
			T[index]->maxminCamera = ref_camera;
			T[index]->t = &(tBuffer[tnOffset]);
			T[index]->normal = &(nBuffer[tnOffset]);

			V10[index] = DesP[0]*u + DesP[1]*v + DesP[2];
			V11[index] = DesP[3]*u + DesP[4]*v + DesP[5];
			V12[index] = DesP[6]*u + DesP[7]*v + DesP[8];

			compute3DIntervals(cache[ref_camera].wedges[wedge_cache_index], current_el, E0,
						RefP, RefCminusDesC, RefSilhouette,T[index],
						 V10, V11,V12,index, ref_camera);

			tnOffset+= T[index]->count;
														
			assert(T[index]->count % 2 == 0);

			// update indexes
			index +=DX_SAMPLE_RATE;

		}
	}
			
}


void sampleConeAndMerge(LINE** T, int ref_camera, int des_camera, int iteration)
{

	unsigned char* 	RefSilhouette = Silhouettes[ref_camera];
	unsigned char* 	DesSilhouette = Silhouettes[des_camera];

	float DesP[9];
	float DesC[3];


	DesP[0] = Ps[des_camera][0];
	DesP[1] = Ps[des_camera][1];
	DesP[2] = Ps[des_camera][2];
	DesP[3] = Ps[des_camera][3];
	DesP[4] = Ps[des_camera][4];
	DesP[5] = Ps[des_camera][5];
	DesP[6] = Ps[des_camera][6];
	DesP[7] = Ps[des_camera][7];
	DesP[8] = Ps[des_camera][8];
	

	DesC[0] = Cs[des_camera][0];
	DesC[1] = Cs[des_camera][1];
	DesC[2] = Cs[des_camera][2];


	float* RefP = Ps[ref_camera];
	float* RefC = Cs[ref_camera];
	
	float* RefCminusDesC = CDifs[ref_camera];
	float* FM = FMs[ref_camera];
	float* TempMatrix = RefPInvDesP[ref_camera];

	EPIPOLE E0 = E0s_des[ref_camera];	

	int index;
	float* tBuffer;
	vector3D* nBuffer;
	float flv, flu;
	tnOffset = 0;
	SLOPEEL current_el;
	int wedge_cache_index;
	float	scanline_infX[3];
	float	infX[3];
	
	
	

	LINE* temp3T = NULL;
	tempT->count = 0;
	tempT->t = tB1;
	tempT->normal = nB1;
	temp2T->count = 0;

	
	if (iteration % 2 == 0)
	{
		// in	buf 2
		// out	buf 1
		tBuffer = tBuffer1;
		nBuffer = normalBuffer1;
	}
	else
	{
		// in	buf 1
		// out	buf 2
		tBuffer = tBuffer2;
		nBuffer = normalBuffer2;
	}


	for (int v = 0; v < DES_IMAGE_HEIGHT; v+=DY_SAMPLE_RATE)
	{
		flv = (float)v;
		
		int scanline_index = v*DES_IMAGE_WIDTH;
		index = scanline_index;
		scanline_infX[0] = TempMatrix[1]*flv + TempMatrix[2];
		scanline_infX[1] = TempMatrix[4]*flv + TempMatrix[5];
		scanline_infX[2] = TempMatrix[7]*flv + TempMatrix[8];
	
		for (int u = 0; u < DES_IMAGE_WIDTH; u+=DX_SAMPLE_RATE)
		{
			if (T[index]->count == 0)
			{
				index +=DX_SAMPLE_RATE;
				continue;
			}

			tempT->count = 0;
			tempT->maxminCamera = ref_camera;
			flu = (float)u;

			infX[2] = TempMatrix[6]*flu + scanline_infX[2];
			if ((E0.negative == true) && (infX[2] <= 0))
			{
				T[index]->count = 0;
				index +=DX_SAMPLE_RATE;
				continue;
			}
		
			infX[0] = TempMatrix[0]*flu + scanline_infX[0];
			infX[1] = TempMatrix[3]*flu + scanline_infX[1];
			if (ComputeSlopeElement(current_el,  E0, infX) == false)
			{
				T[index]->count = 0;				
				index +=DX_SAMPLE_RATE;
				continue;
			}
			
			wedge_cache_index = current_el.wedge_cache_index;

			if (cache[ref_camera].wedges[wedge_cache_index].computed == 0)
				traverse_wedge(wedge_cache_index,E0,RefSilhouette, cache[ref_camera].wedges[wedge_cache_index],ref_camera);


			V10[index] = DesP[0]*u + DesP[1]*v + DesP[2];
			V11[index] = DesP[3]*u + DesP[4]*v + DesP[5];
			V12[index] = DesP[6]*u + DesP[7]*v + DesP[8];

			compute3DIntervals(cache[ref_camera].wedges[wedge_cache_index], current_el, E0,
						RefP, RefCminusDesC, RefSilhouette,tempT,
						 V10, V11,V12,index, ref_camera);


			assert(tempT->count % 2 == 0);			
			if (tempT->count == 0)
			{
				T[index]->count = 0;
			}
			else
			{
				temp2T->t = &(tBuffer[tnOffset]);
				temp2T->normal = &(nBuffer[tnOffset]);
				merge2TLists(*tempT,*T[index],*temp2T);
				tnOffset += temp2T->count;
				temp3T = T[index];
				T[index] = temp2T;
				temp2T = temp3T;
				assert (T[index]->count % 2 == 0);
			}


			// update indexes
			index +=DX_SAMPLE_RATE;

		}
	}
}

void sampleConeSecondPass(LINE** T, int ref_camera, int des_camera)
{
	unsigned char* 	RefSilhouette = Silhouettes[ref_camera];
	unsigned char* 	DesSilhouette = Silhouettes[des_camera];

	float DesP[9];
	float DesC[3];


	DesP[0] = Ps[des_camera][0];
	DesP[1] = Ps[des_camera][1];
	DesP[2] = Ps[des_camera][2];
	DesP[3] = Ps[des_camera][3];
	DesP[4] = Ps[des_camera][4];
	DesP[5] = Ps[des_camera][5];
	DesP[6] = Ps[des_camera][6];
	DesP[7] = Ps[des_camera][7];
	DesP[8] = Ps[des_camera][8];
	

	DesC[0] = Cs[des_camera][0];
	DesC[1] = Cs[des_camera][1];
	DesC[2] = Cs[des_camera][2];


	float* RefP = Ps[ref_camera];
	float* RefC = Cs[ref_camera];
	
	float* RefCminusDesC = CDifs[ref_camera];
	float* FM = FMs[ref_camera];
	float* TempMatrix = RefPInvDesP[ref_camera];

	EPIPOLE E0 = E0s_des[ref_camera];	


	int index;
	float* tBuffer;
	vector3D* nBuffer;
	float flv, flu;
	tnOffset = firstPasstnOffset;
	SLOPEEL current_el;
	int wedge_cache_index;
	float	scanline_infX[3];
	float	infX[3];

	
	tBuffer = tBuffer1;
	nBuffer = normalBuffer1;
	
	// first pass
	for (int v = 0; v < DES_IMAGE_HEIGHT-4; v+=DY_SAMPLE_RATE)
	{
		flv = (float)v;
		
		int scanline_index = v*DES_IMAGE_WIDTH;
		index = scanline_index;	
		for (int u = 0; u < DES_IMAGE_WIDTH-4; u+=DX_SAMPLE_RATE)
		{

			flu = (float)u;
			int index_right = index + DX_SAMPLE_RATE;
			int index_down = index + DY_SAMPLE_RATE*DES_IMAGE_WIDTH;
			int index_down_right = index + DY_SAMPLE_RATE*DES_IMAGE_WIDTH + DX_SAMPLE_RATE;

			bool dothesquare = false;

			if (((T[index]->count > 0 || T[index_right]->count > 0 || T[index_down]->count >0 || T[index_down_right]->count > 0) &&
				(T[index]->count <= 0|| T[index_right]->count <= 0 || T[index_down]->count<= 0 || T[index_down_right]->count <= 0)))
			{

				dothesquare = true;
			}
			
			
			if (((DesSilhouette[index] == BACKGROUND || DesSilhouette[index_right]== BACKGROUND || DesSilhouette[index_down]== BACKGROUND || DesSilhouette[index_down_right]== BACKGROUND) &&
				(DesSilhouette[index]== FOREGROUND || DesSilhouette[index_right]== FOREGROUND || DesSilhouette[index_down]== FOREGROUND || DesSilhouette[index_down_right]== FOREGROUND)))
			{
				dothesquare = true;
			}




			if (T[index]->count > 0 && T[index_right]->count > 0 && T[index_down]->count > 0 && T[index_down_right]->count > 0)
			{

				float ray_x  = DesP[0]*flu + DesP[1]*flv + DesP[2];
				float ray_y  = DesP[3]*flu + DesP[4]*flv + DesP[5];
				float ray_z  = DesP[6]*flu + DesP[7]*flv + DesP[8];

				float nx = T[index]->normal[0].x;
				float ny = T[index]->normal[0].y;
				float nz = T[index]->normal[0].z;
				
				float len = sqrt(nx*nx + ny*ny + nz*nz);
				nx = nx / len;
				ny = ny / len;
				nz = nz / len;

				
				float t1 = T[index]->t[0] ;
				float x1 = t1 * ray_x;
				float y1 = t1 * ray_y;
				float z1 = t1 * ray_z;
				
				float t2 = T[index_right]->t[0];
				float x2 = t2 *(ray_x+ DesP[0]*DX_SAMPLE_RATE);
				float y2 = t2 *(ray_y+ DesP[3]*DX_SAMPLE_RATE);
				float z2 = t2 *(ray_z+ DesP[6]*DX_SAMPLE_RATE);
				
				float t3 = T[index_down]->t[0];
				float x3 = t3 * (ray_x+ DesP[1]*DY_SAMPLE_RATE);
				float y3 = t3 * (ray_y+ DesP[4]*DY_SAMPLE_RATE);
				float z3 = t3 * (ray_z+ DesP[7]*DY_SAMPLE_RATE);
				
				float t4 = T[index_down_right]->t[0];
				float x4 = t4 * (ray_x+DesP[0]*DX_SAMPLE_RATE+DesP[1]*DY_SAMPLE_RATE);
				float y4 = t4 * (ray_y+DesP[3]*DX_SAMPLE_RATE+DesP[4]*DY_SAMPLE_RATE);
				float z4 = t4 * (ray_z+DesP[6]*DX_SAMPLE_RATE+DesP[7]*DY_SAMPLE_RATE);
				
				
				float d1 = nx*x1 + ny*y1 + nz*z1;
				float d2 = nx*x2 + ny*y2 + nz*z2 - d1;
				float d3 = nx*x3 + ny*y3 + nz*z3 - d1;
				float d4 = nx*x4 + ny*y4 + nz*z4 - d1;
				
				
				if (ABS(d2) > DEPTH_THRESH || ABS(d3) > DEPTH_THRESH || ABS(d4) > DEPTH_THRESH)
				{
					dothesquare = true;
				}

			}

/*
			bool bfound = false;
			bool ffound = false;
			for (int iv = v; iv < v+4; iv++)
			{
				for (int iu = u; iu < u+4; iu++)
				{
					int i1 = iv*DES_IMAGE_WIDTH + iu;
					if (DesSilhouette[i1] == BACKGROUND)
						bfound = true;
					if (DesSilhouette[i1] == FOREGROUND)
						ffound = true;
				}
			}

			if (bfound && ffound)
				dothesquare = true;
*/

			if (dothesquare == true)
			{
				for (int iv = v; iv < v+DY_SAMPLE_RATE; iv++)
				{
					flv = (float)iv;
					int scanline_index = iv*DES_IMAGE_WIDTH + u;
					int iindex = scanline_index;
					scanline_infX[0] = TempMatrix[1]*flv + TempMatrix[2];
					scanline_infX[1] = TempMatrix[4]*flv + TempMatrix[5];
					scanline_infX[2] = TempMatrix[7]*flv + TempMatrix[8];
					
					for (int iu = u; iu < u+DX_SAMPLE_RATE; iu++)
					{

						if (iv == v && iu == u)
						{
							iindex++;
							continue;
						}

						if (DesSilhouette[iindex] == BACKGROUND)
						{
							T[iindex]->count = 0;
							iindex ++;
							continue;
						}

						flu = (float)iu;

						T[iindex]->count = 0;
						T[iindex]->maxminCamera = ref_camera;

						infX[2] = TempMatrix[6]*flu + scanline_infX[2];
						if ((E0.negative == true) && (infX[2] <= 0))
						{
							iindex++;
							continue;
						}
						
						infX[0] = TempMatrix[0]*flu + scanline_infX[0];
						infX[1] = TempMatrix[3]*flu + scanline_infX[1];
						if (ComputeSlopeElement(current_el,  E0, infX) == false)
						{			
							iindex++;
							continue;
						}
						
						wedge_cache_index = current_el.wedge_cache_index;
						if (cache[ref_camera].wedges[wedge_cache_index].computed == 0)
							traverse_wedge(wedge_cache_index,E0,RefSilhouette, cache[ref_camera].wedges[wedge_cache_index],ref_camera);
												
						
						T[iindex]->count = 0;
						T[iindex]->maxminCamera = ref_camera;
						T[iindex]->t = &(tBuffer[tnOffset]);
						T[iindex]->normal = &(nBuffer[tnOffset]);
						
						V10[iindex] = DesP[0]*flu + DesP[1]*flv + DesP[2];
						V11[iindex] = DesP[3]*flu + DesP[4]*flv + DesP[5];
						V12[iindex] = DesP[6]*flu + DesP[7]*flv + DesP[8];
						
						compute3DIntervals(cache[ref_camera].wedges[wedge_cache_index], current_el, E0,
							RefP, RefCminusDesC, RefSilhouette,T[iindex],
							V10, V11,V12,iindex, ref_camera);
						
						tnOffset+= T[iindex]->count;
						
						assert(T[iindex]->count % 2 == 0);
						
						// update indexes
						iindex ++;
						
					}
				}
			}

			index +=DX_SAMPLE_RATE;
		}
	}			
}

void sampleConeAndMergeSecondPass(LINE** T, int ref_camera, int des_camera, int iteration)
{
	unsigned char* 	RefSilhouette = Silhouettes[ref_camera];
	unsigned char* 	DesSilhouette = Silhouettes[des_camera];

	float DesP[9];
	float DesC[3];


	DesP[0] = Ps[des_camera][0];
	DesP[1] = Ps[des_camera][1];
	DesP[2] = Ps[des_camera][2];
	DesP[3] = Ps[des_camera][3];
	DesP[4] = Ps[des_camera][4];
	DesP[5] = Ps[des_camera][5];
	DesP[6] = Ps[des_camera][6];
	DesP[7] = Ps[des_camera][7];
	DesP[8] = Ps[des_camera][8];
	

	DesC[0] = Cs[des_camera][0];
	DesC[1] = Cs[des_camera][1];
	DesC[2] = Cs[des_camera][2];


	float* RefP = Ps[ref_camera];
	float* RefC = Cs[ref_camera];
	
	float* RefCminusDesC = CDifs[ref_camera];
	float* FM = FMs[ref_camera];
	float* TempMatrix = RefPInvDesP[ref_camera];

	EPIPOLE E0 = E0s_des[ref_camera];	

	int index;
	float* tBuffer;
	vector3D* nBuffer;
	float flv, flu;
	tnOffset = firstPasstnOffset;
	SLOPEEL current_el;
	int wedge_cache_index;
	float	scanline_infX[3];
	float	infX[3];
	
	
	LINE* temp3T = NULL;
	tempT->count = 0;
	tempT->t = tB1;
	tempT->normal = nB1;
	temp2T->count = 0;

	
	if (iteration % 2 == 0)
	{
		// in	buf 2
		// out	buf 1
		tBuffer = tBuffer1;
		nBuffer = normalBuffer1;
	}
	else
	{
		// in	buf 1
		// out	buf 2
		tBuffer = tBuffer2;
		nBuffer = normalBuffer2;
	}


	for (int v = 0; v < DES_IMAGE_HEIGHT-4; v+=DY_SAMPLE_RATE)
	{
		flv = (float)v;
		
		int scanline_index = v*DES_IMAGE_WIDTH;
		index = scanline_index;
	
		for (int u = 0; u < DES_IMAGE_WIDTH-4; u+=DX_SAMPLE_RATE)
		{

			int index_right = index + DX_SAMPLE_RATE;
			int index_down = index + DY_SAMPLE_RATE*DES_IMAGE_WIDTH;
			int index_down_right = index + DY_SAMPLE_RATE*DES_IMAGE_WIDTH + DX_SAMPLE_RATE;
			flu = (float)u;
			bool dothesquare = false;

			if (((T[index]->count > 0 || T[index_right]->count > 0 || T[index_down]->count >0 || T[index_down_right]->count > 0) &&
				(T[index]->count <= 0|| T[index_right]->count <= 0 || T[index_down]->count<= 0 || T[index_down_right]->count <= 0)))
			{

				dothesquare = true;
			}

			if (((DesSilhouette[index] == BACKGROUND || DesSilhouette[index_right]== BACKGROUND || DesSilhouette[index_down]== BACKGROUND || DesSilhouette[index_down_right]== BACKGROUND) &&
				(DesSilhouette[index]== FOREGROUND || DesSilhouette[index_right]== FOREGROUND || DesSilhouette[index_down]== FOREGROUND || DesSilhouette[index_down_right]== FOREGROUND)))
			{

				dothesquare = true;
			}

			if (T[index]->count > 0 && T[index_right]->count > 0 && T[index_down]->count > 0 && T[index_down_right]->count > 0)
			{

				float ray_x  = DesP[0]*flu + DesP[1]*flv + DesP[2];
				float ray_y  = DesP[3]*flu + DesP[4]*flv + DesP[5];
				float ray_z  = DesP[6]*flu + DesP[7]*flv + DesP[8];

				float nx = T[index]->normal[0].x;
				float ny = T[index]->normal[0].y;
				float nz = T[index]->normal[0].z;
				
				float len = sqrt(nx*nx + ny*ny + nz*nz);
				nx = nx / len;
				ny = ny / len;
				nz = nz / len;

				
				float t1 = T[index]->t[0] ;
				float x1 = t1 * ray_x;
				float y1 = t1 * ray_y;
				float z1 = t1 * ray_z;
				
				float t2 = T[index_right]->t[0];
				float x2 = t2 *(ray_x+ DesP[0]*DX_SAMPLE_RATE);
				float y2 = t2 *(ray_y+ DesP[3]*DX_SAMPLE_RATE);
				float z2 = t2 *(ray_z+ DesP[6]*DX_SAMPLE_RATE);
				
				float t3 = T[index_down]->t[0];
				float x3 = t3 * (ray_x+ DesP[1]*DY_SAMPLE_RATE);
				float y3 = t3 * (ray_y+ DesP[4]*DY_SAMPLE_RATE);
				float z3 = t3 * (ray_z+ DesP[7]*DY_SAMPLE_RATE);
				
				float t4 = T[index_down_right]->t[0];
				float x4 = t4 * (ray_x + DesP[0]*DX_SAMPLE_RATE + DesP[1]*DY_SAMPLE_RATE);
				float y4 = t4 * (ray_y + DesP[3]*DX_SAMPLE_RATE + DesP[4]*DY_SAMPLE_RATE);
				float z4 = t4 * (ray_z + DesP[6]*DX_SAMPLE_RATE + DesP[7]*DY_SAMPLE_RATE);
				
				
				float d1 = nx*x1 + ny*y1 + nz*z1;
				float d2 = nx*x2 + ny*y2 + nz*z2 - d1;
				float d3 = nx*x3 + ny*y3 + nz*z3 - d1;
				float d4 = nx*x4 + ny*y4 + nz*z4 - d1;
				
				if (ABS(d2) > DEPTH_THRESH || ABS(d3) > DEPTH_THRESH || ABS(d4) > DEPTH_THRESH)
				{
					dothesquare = true;
				}

			}

/*
			bool bfound = false;
			bool ffound = false;
			for (int iv = v; iv < v+4; iv++)
			{
				for (int iu = u; iu < u+4; iu++)
				{
					int i1 = iv*DES_IMAGE_WIDTH + iu;
					if (DesSilhouette[i1] == BACKGROUND)
						bfound = true;
					if (DesSilhouette[i1] == FOREGROUND)
						ffound = true;
				}
			}

			if (bfound && ffound)
				dothesquare = true;
*/

			if (dothesquare == true)
			{
				for (int iv = v; iv < v+DY_SAMPLE_RATE; iv++)
				{
					flv = (float)iv;
					int scanline_index = iv*DES_IMAGE_WIDTH + u;
					int iindex = scanline_index;
					scanline_infX[0] = TempMatrix[1]*flv + TempMatrix[2];
					scanline_infX[1] = TempMatrix[4]*flv + TempMatrix[5];
					scanline_infX[2] = TempMatrix[7]*flv + TempMatrix[8];
					
					for (int iu = u; iu < u+DX_SAMPLE_RATE; iu++)
					{

						if (iv == v && iu == u)
						{
							iindex++;
							continue;
						}
						
						if (T[iindex]->count == 0)
						{
							iindex ++;
							continue;
						}
						
						
						flu = (float)iu;
						
						tempT->count = 0;
						tempT->maxminCamera = ref_camera;

						infX[2] = TempMatrix[6]*flu + scanline_infX[2];
						if ((E0.negative == true) && (infX[2] <= 0))
						{
							iindex ++;
							continue;
						}
						
						infX[0] = TempMatrix[0]*flu + scanline_infX[0];
						infX[1] = TempMatrix[3]*flu + scanline_infX[1];
						if (ComputeSlopeElement(current_el,  E0, infX) == false)
						{			
							iindex ++;
							continue;
						}
						
						wedge_cache_index = current_el.wedge_cache_index;
						
						if (cache[ref_camera].wedges[wedge_cache_index].computed == 0)
							traverse_wedge(wedge_cache_index,E0,RefSilhouette, cache[ref_camera].wedges[wedge_cache_index],ref_camera);
						
						
						V10[iindex] = DesP[0]*flu + DesP[1]*flv + DesP[2];
						V11[iindex] = DesP[3]*flu + DesP[4]*flv + DesP[5];
						V12[iindex] = DesP[6]*flu + DesP[7]*flv + DesP[8];
						
						compute3DIntervals(cache[ref_camera].wedges[wedge_cache_index], current_el, E0,
							RefP, RefCminusDesC, RefSilhouette,tempT,
							V10, V11,V12,iindex, ref_camera);
						
						
						assert(tempT->count % 2 == 0);			
						if (tempT->count == 0)
						{
							T[iindex]->count = 0;
						}
						else
						{
							temp2T->t = &(tBuffer[tnOffset]);
							temp2T->normal = &(nBuffer[tnOffset]);
							merge2TLists(*tempT,*T[iindex],*temp2T);
							tnOffset += temp2T->count;
							temp3T = T[iindex];
							T[iindex] = temp2T;
							temp2T = temp3T;
							assert (T[iindex]->count % 2 == 0);
						}
						iindex++;

					}
				}
			}			

			// update indexes
			index +=DX_SAMPLE_RATE;

		}
	}

}







bool ComputeSlopeElement(SLOPEEL& element, EPIPOLE& E0, float* Einf)
{

	// Line Intersections
	float interU[2]; /* U coordinates */
	float interV[2]; /* V coordinates */

	float Interval1Start[2];
	float Interval1End[2];
    float Interval2Start[2];
	float Interval2End[2];
	float dx;
	float dy;
	float inv;
	float eInfUfl;
	float eInfVfl;
	float dyOverdx;

	float*	Start = NULL;
	float*	End = NULL;
	
	element.direction = true;
	element.wedge_cache_index = -1;


	if (Einf[2] == 0)
	{
		element.direction = E0.positive;
		//dx = Einf[0] - E0.ufl;
		//dy = Einf[1] - E0.vfl;
		dyOverdx = Einf[1] / Einf[0];

		int ret= compute_bbox_intersections(E0, dyOverdx,interU, interV);
		if ( ret != 2)
		{	
			if (ret == 0)
			{
				// line has no intersections with image bbox
				element.startX = -1;
				element.startY = -1;		
				element.endX = -1;
				element.endY = -1;
				if (Einf[0] == 0.0f)
					element.slope = -FLT_MAX +1;
				else
					element.slope = dyOverdx;
				return false;
			}
			else 
			{
				// line intersects only one corner
				interU[1] = interU[0];
				interV[1] = interV[0];
			}
		}

		element.wedge_cache_index = compute_wedge_cache_index(interU[0],interV[0],interU[1], interV[1],E0.ufl, E0.vfl);
		
		if (Einf[0] == 0.0f)
			element.slope = -FLT_MAX +1;
		else
			element.slope = dyOverdx;
		element.startX = interU[0];
		element.startY = interV[0];		
		element.endX = interU[1];
		element.endY = interV[1];


		return true;
	}

	if ((E0.positive == false) && (Einf[2] >= 0))
	{
		/* (E0[2] <= 0) && (Einf[2] >= 0) */
		/* Case 3 */
		element.direction = false;
		
		// compute intersection (u, v) -> E0_ref  with desired image

		inv = 1/Einf[2];
		eInfUfl = Einf[0]*inv;
		eInfVfl = Einf[1]*inv;


		// compute intersections
		dx = eInfUfl - E0.ufl;
		dy = eInfVfl - E0.vfl;
		dyOverdx = dy / dx;
		int ret= compute_bbox_intersections(E0, dyOverdx,interU, interV);
		if ( ret != 2)
		{	
			if (ret == 0)
			{
				// line has no intersections with image bbox
				element.startX = -1;
				element.startY = -1;		
				element.endX = -1;
				element.endY = -1;
				if (dx == 0.0f)
					element.slope = -FLT_MAX +1;
				else
					element.slope = dyOverdx;
				return false;
			}
			else 
			{
				// line intersects only one corner
				interU[1] = interU[0];
				interV[1] = interV[0];
			}
		}

		element.wedge_cache_index = compute_wedge_cache_index(interU[0],interV[0],interU[1], interV[1],E0.ufl, E0.vfl);

		Interval1Start[0] = interU[0];
		Interval1Start[1] = interV[0];

		Interval1End[0] = interU[1];
		Interval1End[1] = interV[1];

		Interval2Start[0] = E0.ufl;
		Interval2Start[1] = E0.vfl;

		Interval2End[0] = eInfUfl;
		Interval2End[1] = eInfVfl;

		// compute intersection with the line segment E0 -> Einf

		if (compute_effective_intersection3(Interval1Start, Interval1End,
									   Interval2Start, Interval2End, &Start, &End) == false)
		{
			element.startX = -1;
			element.startY = -1;		
			element.endX = -1;
			element.endY = -1;
			if (dx == 0.0f)
				element.slope = -FLT_MAX +1;
			else
				element.slope = dyOverdx;
			return false;
		}

	}
	else if ((E0.positive == true) && (Einf[2] >= 0))
	{
		/* case 1 */
		element.direction = true;
		// compute intersection (u, v) -> E0_ref  with desired image						
		inv = 1/Einf[2];
		eInfUfl = Einf[0]*inv;
		eInfVfl = Einf[1]*inv;


		// compute intersections
		dx = eInfUfl - E0.ufl;
		dy = eInfVfl - E0.vfl;

		dyOverdx = dy / dx;
		int ret= compute_bbox_intersections(E0, dyOverdx,interU, interV);
		if ( ret != 2)
		{	
			if (ret == 0)
			{
				element.startX = -1;
				element.startY = -1;		
				element.endX = -1;
				element.endY = -1;
				if (dx == 0.0f)
					element.slope = -FLT_MAX +1;
				else
					element.slope = dyOverdx;
				return false;
			}
			else
			{
				// line intersects only one corner
				interU[1] = interU[0];
				interV[1] = interV[0];
			}
		}

		element.wedge_cache_index = compute_wedge_cache_index(interU[0],interV[0],interU[1], interV[1],E0.ufl, E0.vfl);

		Interval1Start[0] = interU[0];
		Interval1Start[1] = interV[0];

		Interval1End[0] = interU[1];
		Interval1End[1] = interV[1];

		Interval2Start[0] = E0.ufl;
		Interval2Start[1] = E0.vfl;

		Interval2End[0] = eInfUfl;
		Interval2End[1] = eInfVfl;

		// compute intersection with the line segment E0 -> Einf
		if (compute_effective_intersection1(Interval1Start, Interval1End,
									   Interval2Start, Interval2End, &Start, &End) == false)
		{
			element.startX = -1;
			element.startY = -1;		
			element.endX = -1;
			element.endY = -1;
			if (dx == 0.0f)
				element.slope = -FLT_MAX +1;
			else
				element.slope = dyOverdx;
			return false;
		}
	}

	else //if ((E0.positive == true) && (Einf[2] <= 0))
	{
		/* Case 2 */
		element.direction = true;

		// compute intersection (u, v) -> E0_ref  with desired image						
		inv = 1/Einf[2];
		eInfUfl = Einf[0]*inv;
		eInfVfl = Einf[1]*inv;
		dx = eInfUfl - E0.ufl;
		dy = eInfVfl - E0.vfl;
		dyOverdx = dy / dx;
		int ret= compute_bbox_intersections(E0, dyOverdx,interU, interV);
		if ( ret != 2)
		{	
			if (ret == 0)
			{
				element.startX = -1;
				element.startY = -1;		
				element.endX = -1;
				element.endY = -1;
				if (dx == 0.0f)
					element.slope = -FLT_MAX +1;
				else
					element.slope = dyOverdx;
				return false;
			}else
			{
				// line intersects only one corner
				interU[1] = interU[0];
				interV[1] = interV[0];
			}
		}
		element.wedge_cache_index = compute_wedge_cache_index(interU[0],interV[0],interU[1], interV[1],E0.ufl, E0.vfl);

		Interval1Start[0] = interU[0];
		Interval1Start[1] = interV[0];

		Interval1End[0] = interU[1];
		Interval1End[1] = interV[1];

		Interval2Start[0] = E0.ufl;
		Interval2Start[1] = E0.vfl;

		Interval2End[0] = eInfUfl;
		Interval2End[1] = eInfVfl;

		// compute intersection with the line segment E0 -> Einf

		if (compute_effective_intersection2(Interval1Start, Interval1End,
									   Interval2Start, Interval2End, &Start, &End) == false)
		{
			element.startX = -1;
			element.startY = -1;		
			element.endX = -1;
			element.endY = -1;
			if (dx == 0.0f)
				element.slope = -FLT_MAX +1;
			else
				element.slope = dyOverdx;
			return false;
		}
	}

	if (dx == 0.0f)
		element.slope = -FLT_MAX +1;
	else
		element.slope = dyOverdx;
	element.startX = Start[0];
	element.startY = Start[1];		
	element.endX = End[0];
	element.endY = End[1];


	return true;
}

int compute_bbox_intersections(const EPIPOLE& e0, float dyOverdx, float* interU,float* interV)
{
	int		index = 0;

	// test line X3 = 0
	float Y3 = e0.vfl - (e0.ufl*dyOverdx);
	if ((Y3 > 0) && (Y3 < REF_IMAGE_HEIGHT-1))
	{
		interU[index] = 0;
		interV[index] = Y3;
		index++;
	}

	// test line X4 = IMAGE_WIDTH-1
	float Y4 = Y3 + ((REF_IMAGE_WIDTH-1)*dyOverdx);
	if ((Y4 >  0) && (Y4 < REF_IMAGE_HEIGHT-1))
	{
		interU[index] = REF_IMAGE_WIDTH-1;
		interV[index] = Y4;
		index++;
		if (index == 2) return index;
	}

	const float dxOverdy = 1.0f/dyOverdx;

	// test line Y1 = 0

	float X1 = e0.ufl - (e0.vfl * dxOverdy);
	if ((X1 >= 0) && (X1 <= REF_IMAGE_WIDTH-1))
	{
		interU[index] = X1;
		interV[index] = 0;
		index++;
		if (index == 2) return index;
	}

	// test line Y2 = (REF_IMAGE_HEIGHT-1)
	
	float X2 = X1 + ((REF_IMAGE_HEIGHT-1)* dxOverdy);
	if ((X2 >= 0) && (X2 <= REF_IMAGE_WIDTH-1))
	{
		interU[index] = X2;
		interV[index] = REF_IMAGE_HEIGHT-1;
		index++;
	}

	return index;
}

// CASE #1 
// Page 144
int compute_effective_intersection1(float* Interval1Start,float* Interval1End,
								   float* Xo,float* Xinf,
								   float** Start,float** End)
{
	float temp;
	bool SWAP = false;
	float dx = ABS(Xo[0]-Xinf[0]);
	float dy = ABS(Xo[1]-Xinf[1]);
	if (dy > dx)
	{
		SWAP = true;
		temp = Interval1Start[0];
		Interval1Start[0] = Interval1Start[1];
		Interval1Start[1] = temp;

		temp = Interval1End[0];
		Interval1End[0] = Interval1End[1];
		Interval1End[1] = temp;

		temp = Xo[0];
		Xo[0] = Xo[1];
		Xo[1] = temp;

		temp = Xinf[0];
		Xinf[0] = Xinf[1];
		Xinf[1] = temp;
	}

	if (Xo[0] > Xinf[0])
	{
		if (Xo[0] < MIN(Interval1Start[0], Interval1End[0]))
		{
			//				------
			//  Xinf<----Xo
			return false;
		}
		else
		if (Xinf[0] > MAX(Interval1Start[0], Interval1End[0]))
		{
			//				--------
			//				          Xinf<-----Xo
			return false;
		}
		else
		{
			//				------
			//		        Xinf<------Xo
			*Start =&MIN(Xo[0],MAX(Interval1Start[0], Interval1End[0]));
			*End = &MAX( Xinf[0],MIN(Interval1Start[0], Interval1End[0]));
		}
	}
	else
	{
		if (Xo[0] > MAX(Interval1Start[0], Interval1End[0]))
		{
			//				------
			//                      Xo----->Xinf
			return false;
		}
		else
		if (Xinf[0] < MIN(Interval1Start[0], Interval1End[0]))
		{
			//				------
			// Xo----->Xinf
			return false;
		}
		else
		{
			//				------
			//       Xo----->Xinf
			*Start =&MAX(Xo[0], MIN(Interval1Start[0], Interval1End[0]));
			*End = &MIN(Xinf[0],MAX(Interval1Start[0], Interval1End[0]));
		}
	}

	if (SWAP == true)
	{
		temp = Interval1Start[0];
		Interval1Start[0] = Interval1Start[1];
		Interval1Start[1] = temp;

		temp = Interval1End[0];
		Interval1End[0] = Interval1End[1];
		Interval1End[1] = temp;

		temp = Xo[0];
		Xo[0] = Xo[1];
		Xo[1] = temp;

		temp = Xinf[0];
		Xinf[0] = Xinf[1];
		Xinf[1] = temp;
	}
	return true;

}



// CASE #2 
// Page 144
int compute_effective_intersection2(float* Interval1Start,float* Interval1End,
								   float* Xo,float* Xinf,
								   float** Start,float** End)
{
	float temp;
	bool SWAP = false;
	float dx = ABS(Xo[0]-Xinf[0]);
	float dy = ABS(Xo[1]-Xinf[1]);
	if (dy > dx)
	{
		SWAP = true;
		temp = Interval1Start[0];
		Interval1Start[0] = Interval1Start[1];
		Interval1Start[1] = temp;

		temp = Interval1End[0];
		Interval1End[0] = Interval1End[1];
		Interval1End[1] = temp;

		temp = Xo[0];
		Xo[0] = Xo[1];
		Xo[1] = temp;

		temp = Xinf[0];
		Xinf[0] = Xinf[1];
		Xinf[1] = temp;
	}

	if (Xo[0] > Xinf[0])
	{
		if (Xo[0] < MIN(Interval1Start[0], Interval1End[0]))
		{
			//				------
			//  Xinf   Xo------------->
			*Start = &MIN(Interval1Start[0], Interval1End[0]);
			*End = &MAX(Interval1Start[0], Interval1End[0]);
		}
		else
		if (Xo[0] > MAX(Interval1Start[0], Interval1End[0]))
		{
			//				------
			//				 Xinf   Xo---------->
			return false;
		}
		else
		{
			//				------
			//		  Xinf    Xo------------->
			*Start = Xo;
			*End = &MAX(Interval1Start[0], Interval1End[0]);
		}
	}
	else
	{
		if (Xo[0] < MIN(Interval1Start[0], Interval1End[0]))
		{
			//				------
			//  <--------Xo    Xinf
			return false;
		}
		else
		if (Xo[0] > MAX(Interval1Start[0], Interval1End[0]))
		{
			//				------
			//	 <------------------Xo  Xinf
			*Start = &MAX(Interval1Start[0], Interval1End[0]);
			*End = &MIN(Interval1Start[0], Interval1End[0]);
		}
		else
		{
			//				------
			//  <--------------Xo   Xinf
			*Start = Xo;
			*End = &MIN(Interval1Start[0], Interval1End[0]);
		}
	}

	if (SWAP == true)
	{
		temp = Interval1Start[0];
		Interval1Start[0] = Interval1Start[1];
		Interval1Start[1] = temp;

		temp = Interval1End[0];
		Interval1End[0] = Interval1End[1];
		Interval1End[1] = temp;

		temp = Xo[0];
		Xo[0] = Xo[1];
		Xo[1] = temp;

		temp = Xinf[0];
		Xinf[0] = Xinf[1];
		Xinf[1] = temp;
	}
	return true;

}


// CASE #3 
// Page 144
int compute_effective_intersection3(float* Interval1Start,float* Interval1End,
								   float* Xo,float* Xinf,
								   float** Start,float** End)
{
	float temp;
	bool SWAP = false;
	float dx = ABS(Xo[0]-Xinf[0]);
	float dy = ABS(Xo[1]-Xinf[1]);
	if (dy > dx)
	{
		SWAP = true;
		temp = Interval1Start[0];
		Interval1Start[0] = Interval1Start[1];
		Interval1Start[1] = temp;

		temp = Interval1End[0];
		Interval1End[0] = Interval1End[1];
		Interval1End[1] = temp;

		temp = Xo[0];
		Xo[0] = Xo[1];
		Xo[1] = temp;

		temp = Xinf[0];
		Xinf[0] = Xinf[1];
		Xinf[1] = temp;
	}

	if (Xo[0] > Xinf[0])
	{
		if (Xinf[0] < MIN(Interval1Start[0], Interval1End[0]))
		{
			//				 ------
			//  -->---Xinf     Xo
			return false;
		}
		else
		if (Xinf[0] > MAX(Interval1Start[0], Interval1End[0]))
		{
			//				------
			//	--------------------->Xinf   Xo
			*Start = &MIN(Interval1Start[0], Interval1End[0]);
			*End = &MAX(Interval1Start[0], Interval1End[0]);
		}
		else
		{
			//				------
			//	-------------->Xinf   Xo
			*Start = &MIN(Interval1Start[0], Interval1End[0]);
			*End = Xinf;
		}
	}
	else
	{
		if (Xinf[0] < MIN(Interval1Start[0], Interval1End[0]))
		{
			//				------
			//  Xo  Xinf--<--------------------
			*Start = &MAX(Interval1Start[0], Interval1End[0]);
			*End = &MIN(Interval1Start[0], Interval1End[0]);
		}
		else
		if (Xinf[0] > MAX(Interval1Start[0], Interval1End[0]))
		{
			//				------
			//	Xo                  Xinf------<-----------
			return false;
		}
		else
		{
			//				------
			// Xo             Xinf-------------<---------
			*Start = &MAX(Interval1Start[0], Interval1End[0]);
			*End = Xinf;
		}
	}

	if (SWAP == true)
	{
		temp = Interval1Start[0];
		Interval1Start[0] = Interval1Start[1];
		Interval1Start[1] = temp;

		temp = Interval1End[0];
		Interval1End[0] = Interval1End[1];
		Interval1End[1] = temp;

		temp = Xo[0];
		Xo[0] = Xo[1];
		Xo[1] = temp;

		temp = Xinf[0];
		Xinf[0] = Xinf[1];
		Xinf[1] = temp;
	}
	return true;

}

void traverse_wedge(int wedge_index, EPIPOLE& e0, unsigned char* sil, WEDGE& w,int ref_cam)
{
	float x0;
	float y0;

	float mid_x;
	float mid_y;

	float vPoint[2];
	vPoint[0] = e0.ufl;
	vPoint[1] = e0.vfl;



	bool inside = true;
	// determine if center is inside the image
	if (vPoint[0] > REF_IMAGE_WIDTH || vPoint[0] < 0)
		inside = false;

	if (vPoint[1] > REF_IMAGE_HEIGHT || vPoint[1] < 0)
		inside = false;


	if (wedge_index < REF_IMAGE_WIDTH)
	{
		x0 = (float)wedge_index;
		y0 = 0;
		mid_x = x0 + 0.5f;
		mid_y = 0;
	}
	else if (wedge_index < 2*REF_IMAGE_WIDTH)
	{
		x0 = (float)(wedge_index-REF_IMAGE_WIDTH);
		y0 = REF_IMAGE_HEIGHT;
		mid_x = x0 + 0.5f;
		mid_y = REF_IMAGE_HEIGHT;
	} 
	else if (wedge_index < 2*REF_IMAGE_WIDTH+REF_IMAGE_HEIGHT)
	{
		x0 = 0;
		y0 = (float)(wedge_index-REF_IMAGE_WIDTH*2);
		mid_x = 0.0f;
		mid_y = y0+0.5f;
	} 
	else
	{
		x0 = REF_IMAGE_WIDTH;
		y0 = (float)(wedge_index-REF_IMAGE_WIDTH*2-REF_IMAGE_HEIGHT);
		mid_x = REF_IMAGE_WIDTH;
		mid_y = y0+0.5f;
	}
	
	float dx_mid = vPoint[0] - mid_x;
	float dy_mid = vPoint[1] - mid_y;

	float dyOverdx_mid = dy_mid / dx_mid;

	float interU0[2]; 
	float interV0[2]; 

	int current_wedge_count = cache[ref_cam].edge_count;

	if (inside == true)
	{
		traverse_int_line(vPoint[0],vPoint[1], mid_x, mid_y,sil,ref_cam);
	}
	else
	{
		int ret= compute_inner_bbox_intersections(vPoint, dyOverdx_mid,interU0, interV0);
//		assert(ret == 2);
		if (ret != 2)
		{
			w.num_edges = 0;
			w.first_edge_off = 0;
			w.last_edge_off = 0;
			return;
		}
		switchToOrdered(vPoint, interU0, interV0);
		traverse_int_line(interU0[0],interV0[0], interU0[1],interV0[1],sil,ref_cam);


	}
	w.computed = 1;
	if (cache[ref_cam].edge_count == current_wedge_count)
	{
		w.num_edges = 0;
		w.first_edge_off = 0;
		w.last_edge_off = 0;
	}
	else
	{
		w.first_edge_off = current_wedge_count;
		w.num_edges = cache[ref_cam].edge_count-current_wedge_count;
		w.last_edge_off = cache[ref_cam].edge_count-1;
	}	
}

void traverse_int_line(int x1, int y1, int x2, int y2,
					unsigned char* silhouette,int ref_cam)
{ 
	int offset = x1 + REF_IMAGE_WIDTH*y1;
	unsigned char val = silhouette[offset]; 

	int dx = (x1 - x2);
	dx = ABS(dx);
	int dy = (y1 - y2);
	dy = ABS(dy);
	int const1, const2, p, /*x, y,*/ step; 
	int off_step;
	if (dx >= dy) 
	{ 
		if (x1 < x2)
		{
			const1 = 2 * dy; 
			const2 = 2 * (dy - dx); 
			p = 2 * dy - dx; 
			step = (y1 > y2) ? (-1) : (1); 
			off_step = (y1 > y2) ? (-REF_IMAGE_WIDTH) : (REF_IMAGE_WIDTH); 
			if (silhouette[offset] != val)
				add_vertical_edge(x1, y1, val,ref_cam);
			while (x1 < x2) 
			{ 
				if (p < 0)
					p += const1; 
				else 
				{ 
					y1 += step;
					offset += off_step;
					p += const2; 
				}
				++x1;
				offset++;
				if (silhouette[offset] != val)
					add_vertical_edge(x1, y1, val,ref_cam);
			}
			if (val == FOREGROUND)
				add_vertical_edge(x1+1, y1, val,ref_cam);
		}
		else
		{
			const1 = 2 * dy; 
			const2 = 2 * (dy - dx); 
			p = 2 * dy - dx; 
			step = (y1 > y2) ? (-1) : (1);
			off_step = (y1 > y2) ? (-REF_IMAGE_WIDTH) : (REF_IMAGE_WIDTH); 
	
			if (silhouette[offset] != val)
				add_vertical_edge(x1, y1, val,ref_cam);
			while (x1 > x2) 
			{ 
				if (p < 0)
					p += const1; 
				else 
				{ 
					y1 += step; 
					offset += off_step;
					p += const2; 
				}
				--x1;
				offset--;
				if (silhouette[offset] != val)
					add_vertical_edge(x1, y1, val,ref_cam);
			}
			if (val == FOREGROUND)
				add_vertical_edge(x1-1, y1, val,ref_cam);

		}

	} 
	else 
	{ 
		if ( y1 < y2)
		{
			const1 = 2 * dx; 
			const2 = 2 * (dx - dy); 
			p = 2 * dx - dy; 
			step = (x1 > x2) ? (-1) : (1);
			off_step = (x1 > x2) ? (-1) : (1); 			
			if (silhouette[offset] != val)
				add_horizontal_edge(x1, y1, val,ref_cam);
			while (y1 < y2) 
			{	 
				if (p < 0)
				p += const1; 
				else 
				{ 
					x1 += step;
					offset+=off_step;
					p += const2; 
				}
				++y1;
				offset+=REF_IMAGE_WIDTH;
				if (silhouette[offset] != val)
					add_horizontal_edge(x1, y1, val,ref_cam);
			}
			if (val == FOREGROUND)
				add_horizontal_edge(x1, y1+1, val,ref_cam);
		}
		else
		{
			const1 = 2 * dx; 
			const2 = 2 * (dx - dy); 
			p = 2 * dx - dy; 
			step = (x1 > x2) ? (-1) : (1); 
			off_step = (x1 > x2) ? (-1) : (1); 
			if (silhouette[offset] != val)
				add_horizontal_edge(x1, y1, val,ref_cam);
			while (y1 > y2) 
			{	 
				if (p < 0)
				p += const1; 
				else 
				{ 
					x1 += step; 
					offset+=off_step;
					p += const2; 
				}
				--y1;
				offset-=REF_IMAGE_WIDTH;
				if (silhouette[offset] != val)
					add_horizontal_edge(x1, y1, val,ref_cam);
			}
			if (val == FOREGROUND)
				add_horizontal_edge(x1, y1-1, val,ref_cam);
		}
	}
		
} 

inline void add_vertical_edge(int ix1, int iy1, unsigned char& val,int image_no)
{
	edgebuffer_[image_no][cache[image_no].edge_count].e.x1 = ix1;
	edgebuffer_[image_no][cache[image_no].edge_count].e.x2 = ix1;
	edgebuffer_[image_no][cache[image_no].edge_count].e.y1 = iy1-100.0f;
	edgebuffer_[image_no][cache[image_no].edge_count].e.y2 = iy1+100.0f;
	cache[image_no].edge_count++;

	if (val == BACKGROUND)
		val = FOREGROUND;
	else
		val = BACKGROUND;
}

inline void add_horizontal_edge(int ix1, int iy1, unsigned char& val,int image_no)
{
	edgebuffer_[image_no][cache[image_no].edge_count].e.x1 = ix1-100;
	edgebuffer_[image_no][cache[image_no].edge_count].e.x2 = ix1+100;
	edgebuffer_[image_no][cache[image_no].edge_count].e.y1 = iy1;
	edgebuffer_[image_no][cache[image_no].edge_count].e.y2 = iy1;
	cache[image_no].edge_count++;

	if (val == BACKGROUND)
		val = FOREGROUND;
	else
		val = BACKGROUND;
}

int compute_inner_bbox_intersections(float* point, float dyOverdx, float* interU,float* interV)
{
	int		index = 0;

	// test line X3 = 0
	float Y3 = point[1] - (point[0]*dyOverdx)+dyOverdx;
	if ((Y3 >= 1) && (Y3 <= REF_IMAGE_HEIGHT-1))
	{
		interU[index] = 1;
		interV[index] = Y3;
		index++;
	}

	// test line X4 = IMAGE_WIDTH
	float Y4 = Y3 + ((REF_IMAGE_WIDTH-2)*dyOverdx);
	if ((Y4 >=  1) && (Y4 <= REF_IMAGE_HEIGHT-1))
	{
		interU[index] = REF_IMAGE_WIDTH-1;
		interV[index] = Y4;
		index++;
		if (index == 2) return index;
	}

	const float dxOverdy = 1.0f/dyOverdx;

	// test line Y1 = 0

	float X1 = point[0] - (point[1] * dxOverdy) +dxOverdy;
	if ((X1 >= 1) && (X1 <= REF_IMAGE_WIDTH-1))
	{
		interU[index] = X1;
		interV[index] = 1;
		index++;
		if (index == 2) return index;
	}

	// test line Y2 =	REF_IMAGE_HEIGHT
	
	float X2 = X1 + ((REF_IMAGE_HEIGHT-2)* dxOverdy);
	if ((X2 >= 1) && (X2 <= REF_IMAGE_WIDTH-1))
	{
		interU[index] = X2;
		interV[index] = REF_IMAGE_HEIGHT-1;
		index++;
	}

	return index;
}

void switchToOrdered(float* vPoint, float* interU, float* interV)
{
	float dist1;
	float dist2;
	float temp;

	temp = vPoint[0]-interU[0];
	dist1 = temp*temp;
	temp = vPoint[1]-interV[0];
	dist1 += temp*temp;

	temp = vPoint[0]-interU[1];
	dist2 = temp*temp;
	temp = vPoint[1]-interV[1];
	dist2 += temp*temp;

	if (dist1 > dist2)
	{
		// switch
		float tempU;
		float tempV;

		tempU = interU[0];
		tempV = interV[0];

		interU[0] = interU[1];
		interV[0] = interV[1];

		interU[1] = tempU;
		interV[1] = tempV;
	}
}


void compute3DIntervals(WEDGE& w, SLOPEEL& current_el, EPIPOLE& E0,
						float* RefP, float* RefCminusDesC, unsigned char *RefSilhouette,LINE* tempT,
						float* V1_0,float* V1_1,float* V1_2,int index, int order)

{

	int edge_off;
	int edge_inc;
	float startX;
	float startY;
	float endX;
	float endY;

	float x1;
	float x2;
	float y1;
	float y2;
	float x_inter;
	float y_inter;
	
	//iterate over all edges in the interval
	if (current_el.direction == true)
	{
		edge_off = w.first_edge_off;
		edge_inc = 1;
	}
	else
	{
		edge_off = w.last_edge_off;
		edge_inc = -1;
	}
	
	startX = current_el.startX;
	startY = current_el.startY;
	endX = current_el.endX;
	endY = current_el.endY;
	
	if (startX > 0 && startX < REF_IMAGE_WIDTH-1 && startY > 0 && startY <REF_IMAGE_HEIGHT-1)
	{
		// start is inside of the image
		int istartX = (int) startX;
		int istartXp1 = istartX+1;
		int istartY = (int) startY;
		int istartYp1 = istartY+1;
		int leftup = (REF_IMAGE_WIDTH*istartY + istartX);
		int rightup = leftup+1;
		int leftdown = leftup+REF_IMAGE_WIDTH;
		int rightdown = leftdown+1;
		if (RefSilhouette[leftup] == FOREGROUND && RefSilhouette[rightup] == FOREGROUND &&
			RefSilhouette[leftdown] == FOREGROUND && RefSilhouette[rightdown] == FOREGROUND)
		{
			tempT->t[0] = 0;
			tempT->normal[0].x = 0;
			tempT->normal[0].y = 0;
			tempT->normal[0].z = 1;
			tempT->count = 1;
		}
		else if (									  RefSilhouette[rightup] == FOREGROUND &&
			RefSilhouette[leftdown] == FOREGROUND && RefSilhouette[rightdown] == FOREGROUND)
		{
			float xdif = startX - istartX;
			float ydif = startY - istartY;
			if (ydif >= (1.0f-xdif))
			{
				tempT->t[0] = 0;
				tempT->normal[0].x = 0;
				tempT->normal[0].y = 0;
				tempT->normal[0].z = 1;
				tempT->count = 1;
			}
		}
		else if (RefSilhouette[leftup] == FOREGROUND && 
			RefSilhouette[leftdown] == FOREGROUND && RefSilhouette[rightdown] == FOREGROUND)
		{
			float xdif = startX - istartX;
			float ydif = startY - istartY;
			if (ydif >= xdif)
			{
				tempT->t[0] = 0;
				tempT->normal[0].x = 0;
				tempT->normal[0].y = 0;
				tempT->normal[0].z = 1;
				tempT->count = 1;
			}
		}
		else if (RefSilhouette[leftup] == FOREGROUND && RefSilhouette[rightup] == FOREGROUND &&
			RefSilhouette[rightdown] == FOREGROUND)
		{
			float xdif = startX - istartX;
			float ydif = startY - istartY;
			if (ydif <= xdif)
			{
				tempT->t[0] = 0;
				tempT->normal[0].x = 0;
				tempT->normal[0].y = 0;
				tempT->normal[0].z = 1;
				tempT->count = 1;
			}
		}
		else if (RefSilhouette[leftup] == FOREGROUND && RefSilhouette[rightup] == FOREGROUND &&
			RefSilhouette[leftdown] == FOREGROUND )
		{
			float xdif = startX - istartX;
			float ydif = startY - istartY;
			if (ydif <= (1.0f-xdif))
			{
				tempT->t[0] = 0;
				tempT->normal[0].x = 0;
				tempT->normal[0].y = 0;
				tempT->normal[0].z = 1;
				tempT->count = 1;
			}
		}
	}
	

	assert(w.num_edges % 2 == 0);
	for (int j = 0; j < w.num_edges; j++) 
	{	
		x1 = edgebuffer_[order][edge_off].e.x1;
		x2 = edgebuffer_[order][edge_off].e.x2;
		y1 = edgebuffer_[order][edge_off].e.y1;
		y2 = edgebuffer_[order][edge_off].e.y2;
		
		if (x1 == x2) 
		{
			
			// vertical edge
			x_inter = (float)x1;
			y_inter = E0.vfl + current_el.slope*(x_inter - E0.ufl);
			if (ABS(current_el.slope) > 1)
			{
				// use Y
				if  ( (y_inter >= startY && y_inter <= endY) || ( y_inter <= startY && y_inter >= endY) )
				{
					//assert ((x_inter >= startX && x_inter <= endX) || ( x_inter <= startX && x_inter >= endX));
					// calculate intersection point in 3D
					compute_T(tempT->t[tempT->count],RefP,RefCminusDesC,V1_0[index],V1_1[index],V1_2[index],x_inter,y_inter);
					compute_normal(tempT->normal[tempT->count],RefP,x1,y1,x2,y2);
					tempT->count++;
				}
			}
			else
			{
				if  ( (x_inter >= startX && x_inter <= endX) || ( x_inter <= startX && x_inter >= endX) )
				{
					//assert((y_inter >= startY && y_inter <= endY) || ( y_inter <= startY && y_inter >= endY));
					// calculate intersection point in 3D
					compute_T(tempT->t[tempT->count],RefP,RefCminusDesC,V1_0[index],V1_1[index],V1_2[index],x_inter,y_inter);
					compute_normal(tempT->normal[tempT->count],RefP,x1,y1,x2,y2);
					tempT->count++;
				}
			}
		}
		else if (y1 == y2)
		{
			// horizontal edge
			y_inter = (float)y1;
			x_inter = E0.ufl + (y_inter - E0.vfl)/current_el.slope;
			if (ABS(current_el.slope) > 1)
			{
				// use Y
				if  ( (y_inter >= startY && y_inter <= endY) || ( y_inter <= startY && y_inter >= endY) )
				{
					//assert ((x_inter >= startX && x_inter <= endX) || ( x_inter <= startX && x_inter >= endX));
					// calculate intersection point in 3D
					compute_T(tempT->t[tempT->count],RefP,RefCminusDesC,V1_0[index],V1_1[index],V1_2[index],x_inter,y_inter);
					compute_normal(tempT->normal[tempT->count],RefP,x1,y1,x2,y2);
					tempT->count++;
				}
			}
			else
			{
				if  ( (x_inter >= startX && x_inter <= endX) || ( x_inter <= startX && x_inter >= endX) )
				{
					//assert((y_inter >= startY && y_inter <= endY) || ( y_inter <= startY && y_inter >= endY));
					// calculate intersection point in 3D
					compute_T(tempT->t[tempT->count],RefP,RefCminusDesC,V1_0[index],V1_1[index],V1_2[index],x_inter,y_inter);
					compute_normal(tempT->normal[tempT->count],RefP,x1,y1,x2,y2);
					tempT->count++;
				}
			}
			
		} 
		else
		{
			
			float a1 = current_el.slope;
			float a2 = edgebuffer_[order][edge_off].dy_dx;// dy*dx_inv;
			float adif = a1-a2;
			float adif_inv = 1.0f/ adif;
			float dx_inv = edgebuffer_[order][edge_off].dxinv;// 1.0f / dx;
			float b1 = E0.vfl - a1*E0.ufl;
			float dy = y1 - y2;
			float b2 = (y2*x1 - y1*x2)*dx_inv;
			
			if (adif > -0.000001 && adif < 0.000001)
			{
				x_inter = (x1+x2) * 0.5f;
				y_inter = (y1+y2) * 0.5f;
			}
			else
			{
				x_inter = (b2-b1) * adif_inv;
				y_inter = (a1*b2 - b1*a2)*adif_inv;
			}
			if (ABS(a1) < 1)
			{
				if ( x_inter <= MAX(startX ,endX) && x_inter >= MIN(startX, endX))  
					//	&& x_inter <= MAX(x1, x2) && x_inter >= MIN(x1, x2) )
				{
					compute_T(tempT->t[tempT->count],RefP,RefCminusDesC,V1_0[index],V1_1[index],V1_2[index],x_inter,y_inter);
					compute_normal(tempT->normal[tempT->count],RefP,x1,y1,x2,y2);
					tempT->count++;
				}
			}
			else
			{
				if ( y_inter <= MAX(startY ,endY) && y_inter >= MIN(startY, endY))  
					//	&& x_inter <= MAX(x1, x2) && x_inter >= MIN(x1, x2) )
				{
					compute_T(tempT->t[tempT->count],RefP,RefCminusDesC,V1_0[index],V1_1[index],V1_2[index],x_inter,y_inter);
					compute_normal(tempT->normal[tempT->count],RefP,x1,y1,x2,y2);
					tempT->count++;
				}
			}
		}
		
		// get next element from the list
		edge_off += edge_inc;
	}
					
	if (endX > 0 && endX < REF_IMAGE_WIDTH-1 && endY > 0 && endY <REF_IMAGE_HEIGHT-1)
	{
		if (tempT->count % 2 != 0)
		{
			tempT->t[tempT->count] = FLT_MAX;
			tempT->normal[tempT->count].x = 0;
			tempT->normal[tempT->count].y = 0;
			tempT->normal[tempT->count].z = 1;
			tempT->count++;
		}
	}
	assert(tempT->count % 2 == 0);
}

inline bool compute_T(float& t, float* P2, float* C2minC1, float V1_0, float V1_1, float V1_2, float u1, float v1)
{

	float V2[3];
		
	//MultiplyVector(P2, X2, V2);
	V2[0] = P2[0]*u1 + P2[1]*v1 + P2[2];
	V2[1] = P2[3]*u1 + P2[4]*v1 + P2[5];
	V2[2] = P2[6]*u1 + P2[7]*v1 + P2[8];

	
	float V1xV2[3];
	//cross(V1,V2,V1xV2);
	V1xV2[0]= V1_1*V2[2] - V1_2*V2[1];
	V1xV2[1]= V1_2*V2[0] - V1_0*V2[2];
	V1xV2[2]= V1_0*V2[1] - V1_1*V2[0];	
	
	
	float lengthSq = (V1xV2[0]*V1xV2[0])+(V1xV2[1]*V1xV2[1])+(V1xV2[2]*V1xV2[2]);

	if (lengthSq == 0)
	{
		t =INFINITY;
		return true;
	}
	t =	(C2minC1[0]*(V2[1]*V1xV2[2] - V1xV2[1]*V2[2]) +
         V2[0]*(V1xV2[1]*C2minC1[2] - C2minC1[1]*V1xV2[2]) +
         V1xV2[0]*(C2minC1[1]*V2[2] - V2[1]*C2minC1[2]))/lengthSq;

	return true;


}

int compute_wedge_cache_index(float flU, float flV)
{
	if (flU == 1.0 || flV == 1.0)
		flU = 1.0;

	int U = (int) flU;
	int V = (int) flV;
	if (flV == 0.0f)
	{
		return U;
	}
	if (flV == REF_IMAGE_HEIGHT-1)
	{
		return REF_IMAGE_WIDTH + U;
	}
	if (flU == 0.0f)
	{
		return REF_IMAGE_WIDTH*2 + V;
	}
	if (flU == REF_IMAGE_WIDTH-1)
	{
		return REF_IMAGE_WIDTH*2 + REF_IMAGE_HEIGHT + V;
	}
	return -1;
}

int compute_wedge_cache_index(float U1, float V1, float U2, float V2, float vp_U, float vp_V)
{
	float dist1;
	float dist2;
	float temp;

	temp = vp_U - U1;
	dist1 = temp*temp;
	temp = vp_V - V1;
	dist1 += temp*temp;

	temp = vp_U - U2;
	dist2 = temp*temp;
	temp = vp_V - V2;
	dist2 += temp*temp;

	if (dist1 > dist2)
		return compute_wedge_cache_index(U1, V1);
	else
		return compute_wedge_cache_index(U2, V2);

}


void merge2TLists(const LINE& ray1, const LINE& ray2, LINE& rayout)
{
	rayout.count = 0;
	rayout.maxminCamera = -1;
	int ray1index = 0;
	int ray2index = 0;
	float start;
	float end; 

	while (ray1index < ray1.count && ray2index < ray2.count)
	{
		start = MAX(ray1.t[ray1index] , ray2.t[ray2index]);
		end = MIN(ray1.t[ray1index+1] , ray2.t[ray2index+1]);
		if (start <= end)
		{
			// add interval
			rayout.t[rayout.count] = start;
			if (start == ray1.t[ray1index])
				rayout.normal[rayout.count] = ray1.normal[ray1index];
			else
				rayout.normal[rayout.count] = ray2.normal[ray2index];
			rayout.count++;
			//if (rayout.count >= max_intervals) { MessageBox(NULL, "MAX INTERVALS too small","3D SCAN", MB_OK); exit(0);}
			rayout.t[rayout.count] = end;
			if (end == ray1.t[ray1index+1])
				rayout.normal[rayout.count] = ray1.normal[ray1index+1];
			else
				rayout.normal[rayout.count] = ray2.normal[ray2index+1];
			rayout.count++;
			//if (rayout.count >= max_intervals) { MessageBox(NULL, "MAX INTERVALS too small","3D SCAN", MB_OK); exit(0);}
			if (start == ray1.t[ray1index])
				rayout.maxminCamera = ray1.maxminCamera;
			else
				rayout.maxminCamera = ray2.maxminCamera;

			if (ray2.t[ray2index+1]> ray1.t[ray1index+1])
			{
				// ray1 ends first
				// drop ray1 keep ray2
				// move the pointer
				ray1index+=2;
			}
			else
			{
				// ray2 ends first
				// drop ray2 keep ray1
				// move the pointer
				ray2index+=2;
			}
			break;
		}
		if (ray2.t[ray2index+1]> ray1.t[ray1index+1])
		{
			// ray1 ends first
			// drop ray1 keep ray2
			// move the pointer
			ray1index+=2;
		}
		else
		{
			// ray2 ends first
			// drop ray2 keep ray1
			// move the pointer
			ray2index+=2;
		}
	}

	while (ray1index < ray1.count && ray2index < ray2.count)
	{
		start = MAX(ray1.t[ray1index] , ray2.t[ray2index]);
		end = MIN(ray1.t[ray1index+1] , ray2.t[ray2index+1]);
		if (start <= end)
		{
			// add interval
			rayout.t[rayout.count] = start;
			if (start == ray1.t[ray1index])
				rayout.normal[rayout.count] = ray1.normal[ray1index];
			else
				rayout.normal[rayout.count] = ray2.normal[ray2index];
			rayout.count++;
			//if (rayout.count >= max_intervals) { MessageBox(NULL, "MAX INTERVALS too small","3D SCAN", MB_OK); exit(0);}
			rayout.t[rayout.count] = end;
			if (end == ray1.t[ray1index+1])
				rayout.normal[rayout.count] = ray1.normal[ray1index+1];
			else
				rayout.normal[rayout.count] = ray2.normal[ray2index+1];
			rayout.count++;
			//if (rayout.count >= max_intervals) { MessageBox(NULL, "MAX INTERVALS too small","3D SCAN", MB_OK); exit(0);}
		}
		if (ray2.t[ray2index+1]> ray1.t[ray1index+1])
		{
			// ray1 ends first
			// drop ray1 keep ray2
			// move the pointer
			ray1index+=2;
		}
		else
		{
			// ray2 ends first
			// drop ray2 keep ray1
			// move the pointer
			ray2index+=2;
		}
	}
	//if (rayout.count >= max_intervals) { MessageBox(NULL, "MAX INTERVALS too small","3D SCAN", MB_OK); exit(0);}

}

void compute_normal(VECTOR3D& n,float* RefP, float x1,float y1,float x2,float y2)
{

	float vec1[3];
	float vec2[3];

	vec1[0] = RefP[0] *x1 + RefP[1]*y1 +RefP[2];
	vec1[1] = RefP[3] *x1 + RefP[4]*y1 +RefP[5];
	vec1[2] = RefP[6] *x1 + RefP[7]*y1 +RefP[8];
	
	vec2[0] = RefP[0] *x2 + RefP[1]*y2 +RefP[2];
	vec2[1] = RefP[3] *x2 + RefP[4]*y2 +RefP[5];
	vec2[2] = RefP[6] *x2 + RefP[7]*y2 +RefP[8];
	
	n.x = vec1[1]* vec2[2] - vec1[2]*vec2[1];
	n.y = vec1[2]* vec2[0] - vec1[0]*vec2[2];
	n.z = vec1[0]* vec2[1] - vec1[1]*vec2[0];
	
}

void allocateMemory(int Count)
{
	
	FMs = (float **) (new int[Count]);
	CDifs = (float **) (new int[Count]);
	RefPInv  = (float **) (new int[Count]);
	RefPInvDesP = (float **) (new int[Count]);
	E0s_des = new EPIPOLE[Count];
	E0s_ref = new EPIPOLE[Count];		
	iterations = new ITERORDER[Count];
	
	// allocate rays
	T = new LINE*[DES_IMAGE_WIDTH*DES_IMAGE_HEIGHT];
	
	// allocate samples and normals
	tBuffer1 = (float *)malloc(max_samples*sizeof(float));
	tBuffer2 = (float *)malloc(max_samples*sizeof(float));
	
	normalBuffer1 = (VECTOR3D *)malloc(max_samples*sizeof(VECTOR3D));
	normalBuffer2 = (VECTOR3D *)malloc(max_samples*sizeof(VECTOR3D));
	
	tempT = new LINE;
	temp2T = new LINE;
	tB1 = (float* )malloc(1000*sizeof (float));
	nB1 = (VECTOR3D*)malloc(1000*sizeof(VECTOR3D));
	
	for (int index1 = 0; index1 < Count; index1++)
	{
		FMs[index1] = new float[9];
		RefPInv[index1] = new float[9];
		RefPInvDesP[index1] = new float[9];
		CDifs[index1] = new float[3];
	}
}

void cleanBuffer()
{
	int off = 0;
	for (int v = 0; v < DES_IMAGE_HEIGHT; v++)
	{
		for (int u = 0; u < DES_IMAGE_WIDTH; u++)
		{
			T[off] = new LINE;
			T[off]->count = -1;
			off++;
		}
		
	}
}