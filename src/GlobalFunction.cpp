
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdlib.h>
#include <assert.h>

#include <gl/gl.h>
#include <gl/glu.h>

#include <fstream>
extern std::ofstream logout ;

#include "GLDrawer.h"


#include "GlobalFunction.h"

using namespace vcg;
using namespace std;
using namespace tri;

#include "DataMgr.h"




vector<int> GlobalFun::GetRandomCards(int Max)
{
	vector<int> nCard(Max, 0);
	srand(time(NULL));
	for(int i=0; i < Max; i++)
	{
		nCard[i] = i;
	}
	random_shuffle(nCard.begin(), nCard.begin() + Max);


	return nCard;
}

double GlobalFun::computeEulerDist(Point3f& p1, Point3f& p2)
{
	double dist2 = (p1-p2).SquaredNorm();
	if (dist2 < 1e-8 || dist2 > 1e8)
	{
		return 0;
	}
	return sqrt(dist2);
}

double GlobalFun::computeEulerDistSquare(Point3f& p1, Point3f& p2)
{
	return (p1-p2).SquaredNorm();
}

double GlobalFun::computeProjDist(Point3f& p1, Point3f& p2, Point3f& normal_of_p1)
{
	return (p2-p1) * normal_of_p1.Normalize();
}



double GlobalFun::computeProjDistSquare(Point3f& p1, Point3f& p2, Point3f& normal_of_p1)
{
	double proj_dist = computeProjDist(p1, p2, normal_of_p1);
	return proj_dist * proj_dist;
}

double GlobalFun::computePerpendicularDistSquare(Point3f& p1, Point3f& p2, Point3f& normal_of_p1)
{
	//Point3f v_p2_p1 = p1-p2;
	//double proj_dist = computeProjDist(p1, p2, normal_of_p1);
	//Point3f v_proj = /*p1 + */normal_of_p1 * proj_dist;
	//   return (v_p2_p1 + v_proj).SquaredNorm();
	double proj_dist = computeProjDist(p1, p2, normal_of_p1);
	Point3f proj_p = p1 + normal_of_p1 * proj_dist;
	return (proj_p - p2).SquaredNorm();
}



double GlobalFun::computePerpendicularDist(Point3f& p1, Point3f& p2, Point3f& normal_of_p1)
{
	return sqrt(computePerpendicularDistSquare(p1, p2, normal_of_p1));
}





bool GlobalFun::isTwoPoint3fTheSame(Point3f& v0, Point3f& v1)
{
	if (abs(v0[0] - v1[0]) < 1e-7 &&  
		abs(v0[1] - v1[1]) < 1e-7 && 
		abs(v0[2] - v1[2]) < 1e-7)
	{
		return true;
	}

	return false;

}

bool GlobalFun::isTwoPoint3fOpposite(Point3f& v0, Point3f& v1)
{
	if (abs(-v0[0] - v1[0]) < 1e-7 &&  
		abs(-v0[1] - v1[1]) < 1e-7 && 
		abs(-v0[2] - v1[2]) < 1e-7)
	{
		return true;
	}

	return false;
}


double GlobalFun::computeRealAngleOfTwoVertor(Point3f v0, Point3f v1)
{
	v0.Normalize();
	v1.Normalize();


	if (isTwoPoint3fTheSame(v0, v1))
	{
		return 0;
	}

	if (isTwoPoint3fOpposite(v0, v1))
	{
		return 180;
	}

	double angle_cos = v0 * v1;
	if (angle_cos > 1)
	{
		angle_cos = 0.99;
	}
	if (angle_cos < -1)
	{
		angle_cos = -0.99;
	}
	if (angle_cos > 0 && angle_cos < 1e-8)
	{
		return 90;
	}

	double angle = acos(angle_cos) * 180. / 3.1415926 ;

	if (angle < 0 || angle > 180)
	{
		cout << "compute angle wrong!!" << endl;
		//system("Pause");
		return -1;
	}



	return angle;
}

bool GlobalFun::gluInvertMatrix(const float m[16], float invOut[16])
{
	double inv[16], det;
	int i;

	inv[0] = m[5]  * m[10] * m[15] - 
		m[5]  * m[11] * m[14] - 
		m[9]  * m[6]  * m[15] + 
		m[9]  * m[7]  * m[14] +
		m[13] * m[6]  * m[11] - 
		m[13] * m[7]  * m[10];

	inv[4] = -m[4]  * m[10] * m[15] + 
		m[4]  * m[11] * m[14] + 
		m[8]  * m[6]  * m[15] - 
		m[8]  * m[7]  * m[14] - 
		m[12] * m[6]  * m[11] + 
		m[12] * m[7]  * m[10];

	inv[8] = m[4]  * m[9] * m[15] - 
		m[4]  * m[11] * m[13] - 
		m[8]  * m[5] * m[15] + 
		m[8]  * m[7] * m[13] + 
		m[12] * m[5] * m[11] - 
		m[12] * m[7] * m[9];

	inv[12] = -m[4]  * m[9] * m[14] + 
		m[4]  * m[10] * m[13] +
		m[8]  * m[5] * m[14] - 
		m[8]  * m[6] * m[13] - 
		m[12] * m[5] * m[10] + 
		m[12] * m[6] * m[9];

	inv[1] = -m[1]  * m[10] * m[15] + 
		m[1]  * m[11] * m[14] + 
		m[9]  * m[2] * m[15] - 
		m[9]  * m[3] * m[14] - 
		m[13] * m[2] * m[11] + 
		m[13] * m[3] * m[10];

	inv[5] = m[0]  * m[10] * m[15] - 
		m[0]  * m[11] * m[14] - 
		m[8]  * m[2] * m[15] + 
		m[8]  * m[3] * m[14] + 
		m[12] * m[2] * m[11] - 
		m[12] * m[3] * m[10];

	inv[9] = -m[0]  * m[9] * m[15] + 
		m[0]  * m[11] * m[13] + 
		m[8]  * m[1] * m[15] - 
		m[8]  * m[3] * m[13] - 
		m[12] * m[1] * m[11] + 
		m[12] * m[3] * m[9];

	inv[13] = m[0]  * m[9] * m[14] - 
		m[0]  * m[10] * m[13] - 
		m[8]  * m[1] * m[14] + 
		m[8]  * m[2] * m[13] + 
		m[12] * m[1] * m[10] - 
		m[12] * m[2] * m[9];

	inv[2] = m[1]  * m[6] * m[15] - 
		m[1]  * m[7] * m[14] - 
		m[5]  * m[2] * m[15] + 
		m[5]  * m[3] * m[14] + 
		m[13] * m[2] * m[7] - 
		m[13] * m[3] * m[6];

	inv[6] = -m[0]  * m[6] * m[15] + 
		m[0]  * m[7] * m[14] + 
		m[4]  * m[2] * m[15] - 
		m[4]  * m[3] * m[14] - 
		m[12] * m[2] * m[7] + 
		m[12] * m[3] * m[6];

	inv[10] = m[0]  * m[5] * m[15] - 
		m[0]  * m[7] * m[13] - 
		m[4]  * m[1] * m[15] + 
		m[4]  * m[3] * m[13] + 
		m[12] * m[1] * m[7] - 
		m[12] * m[3] * m[5];

	inv[14] = -m[0]  * m[5] * m[14] + 
		m[0]  * m[6] * m[13] + 
		m[4]  * m[1] * m[14] - 
		m[4]  * m[2] * m[13] - 
		m[12] * m[1] * m[6] + 
		m[12] * m[2] * m[5];

	inv[3] = -m[1] * m[6] * m[11] + 
		m[1] * m[7] * m[10] + 
		m[5] * m[2] * m[11] - 
		m[5] * m[3] * m[10] - 
		m[9] * m[2] * m[7] + 
		m[9] * m[3] * m[6];

	inv[7] = m[0] * m[6] * m[11] - 
		m[0] * m[7] * m[10] - 
		m[4] * m[2] * m[11] + 
		m[4] * m[3] * m[10] + 
		m[8] * m[2] * m[7] - 
		m[8] * m[3] * m[6];

	inv[11] = -m[0] * m[5] * m[11] + 
		m[0] * m[7] * m[9] + 
		m[4] * m[1] * m[11] - 
		m[4] * m[3] * m[9] - 
		m[8] * m[1] * m[7] + 
		m[8] * m[3] * m[5];

	inv[15] = m[0] * m[5] * m[10] - 
		m[0] * m[6] * m[9] - 
		m[4] * m[1] * m[10] + 
		m[4] * m[2] * m[9] + 
		m[8] * m[1] * m[6] - 
		m[8] * m[2] * m[5];

	det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

	if (det == 0)
		return false;

	det = 1.0 / det;

	for (i = 0; i < 16; i++)
		invOut[i] = inv[i] * det;

	return true;
}



Point3f GlobalFun::vecMulMatrix( Point3f point0, float invmvMatrix[16] ) {


	Point3f res(0,0,0) ;

	float vec[4] ;
	vec[0] = point0.X() ;
	vec[1] = point0.Y() ;
	vec[2] = point0.Z() ;
	vec[3] = 1.0 ;

	double res3 = 0.0 ;

	for( int i = 0; i<4; ++i )
		res.X() += vec[i] * invmvMatrix[i*4+0] ;

	for( int i = 0; i<4; ++i )
		res.Y() += vec[i] * invmvMatrix[i*4+1] ;

	for( int i = 0; i<4; ++i )
		res.Z() += vec[i] * invmvMatrix[i*4+2] ;

	for( int i = 0; i<4; ++i )
		res3 += vec[i] * invmvMatrix[i*4+3] ;

	return res/res3 ;

}


// 功能：判断点是否在多边形内 
// 方法：求解通过该点的水平线与多边形各边的交点 
// 结论：单边交点为奇数，成立!
//参数： 
// POINT p 指定的某个点 
// LPPOINT ptPolygon 多边形的各个顶点坐标（首末点可以不一致） 
// int nCount 多边形定点的个数

bool GlobalFun::PtInPolygon (Point2f p, std::vector<Point2f> &ptPolygon ) { 
	int nCross = 0;
	int nCount = ptPolygon.size() ;

	for (int i = 0; i < nCount; i++) 
	{ 
		Point2f p1 = ptPolygon[i]; 
		Point2f p2 = ptPolygon[(i + 1) % nCount];
		// 求解 y=p.y 与 p1p2 的交点
		if ( p1.Y() == p2.Y()) // p1p2 与 y=p0.y平行 
			continue;
		if ( p.Y()< min(p1.Y(), p2.Y()) ) // 交点在p1p2延长线上 
			continue; 
		if ( p.Y()>= max(p1.Y(), p2.Y()) ) // 交点在p1p2延长线上 
			continue;
		// 求交点的 X 坐标 -------------------------------------------------------------- 
		double x = (double)(p.Y()- p1.Y()) * (double)(p2.X() - p1.X()) / (double)(p2.Y() - p1.Y()) + p1.X();
		if ( x > p.X() ) 
			nCross++; // 只统计单边交点 
	}
	// 单边交点为偶数，点在多边形之外 --- 

	return (nCross % 2 == 1); 
}


double3 GlobalFun::angle2color( double angle ) {


	double dis = angle / (30.0 * 3.14/180.0) ;
	
	dis = (std::min)( (std::max)(dis, 0.0 ), 1.0 ) ;
	dis = 1.0 - dis ;

	double3 baseColor[9];
	baseColor[0] = double3(1.0, 0.0, 0.0) ;
	baseColor[1] = double3(1.0, 0.7, 0.0) ;
	baseColor[2] = double3(1.0, 1.0, 0.0) ;
	baseColor[3] = double3(0.7, 1.0, 0.0) ;
	baseColor[4] = double3(0.0, 1.0, 0.0) ;
	baseColor[5] = double3(0.0, 1.0, 0.7) ;
	baseColor[6] = double3(0.0, 1.0, 1.0) ;
	baseColor[7] = double3(0.0, 0.7, 1.0) ;
	baseColor[8] = double3(0.0, 0.0, 1.0) ;


	double step = 1.0 / 8.0 ;

	int baseID = dis/step;
	if( baseID == 8 )
		baseID = 7 ;

	double3 mixedColor =  baseColor[baseID] * (baseID*step+step- dis)  + baseColor[baseID+1] * (dis-baseID*step) ;
	mixedColor = (mixedColor/step)  ;
	return mixedColor ;

}



double3 GlobalFun::scalar2color( double scalar ) {


	double dis = scalar ;

	dis = (std::min)( (std::max)(dis, 0.0 ), 1.0 ) ;
	dis = dis ;

	double3 baseColor[9];
	baseColor[0] = double3(1.0, 0.0, 0.0) ;
	baseColor[1] = double3(1.0, 0.7, 0.0) ;
	baseColor[2] = double3(1.0, 1.0, 0.0) ;
	baseColor[3] = double3(0.7, 1.0, 0.0) ;
	baseColor[4] = double3(0.0, 1.0, 0.0) ;
	baseColor[5] = double3(0.0, 1.0, 0.7) ;
	baseColor[6] = double3(0.0, 1.0, 1.0) ;
	baseColor[7] = double3(0.0, 0.2, 1.0) ;
	baseColor[8] = double3(0.0, 0.0, 1.0) ;


	baseColor[0] = double3(1.0, 0.0, 0.0) ;
	baseColor[1] = double3(1.0, 1.0, 0.0) ;
	baseColor[2] = double3(0.3, 1.0, 0.0) ;
	baseColor[3] = double3(0.7, 1.0, 0.0) ;
	baseColor[4] = double3(0.0, 1.0, 0.0) ;
	baseColor[5] = double3(0.0, 1.0, 1.0) ;
	baseColor[6] = double3(0.0, 0.7, 1.0) ;
	baseColor[7] = double3(0.0, 0.3, 1.0) ;
	baseColor[8] = double3(0.0, 0.0, 1.0) ;

	double step = 1.0 / 8.0 ;

	int baseID = dis/step;
	if( baseID == 8 )
		baseID = 7 ;

	double3 mixedColor =  baseColor[baseID] * (baseID*step+step- dis)  + baseColor[baseID+1] * (dis-baseID*step) ;
	mixedColor = (mixedColor/step)  ;
	return mixedColor ;

}



void RotateArbitraryAxis(float m[4][4], Point3f axis, float theta){

	axis.Normalize() ;

	float u = axis[0];
	float v = axis[1];
	float w = axis[2];

	m[0][0] = cosf(theta) + (u * u) * (1 - cosf(theta));
	m[0][1] = u * v * (1 - cosf(theta)) + w * sinf(theta);
	m[0][2] = u * w * (1 - cosf(theta)) - v * sinf(theta);
	m[0][3] = 0;

	m[1][0] = u * v * (1 - cosf(theta)) - w * sinf(theta);
	m[1][1] = cosf(theta) + v * v * (1 - cosf(theta));
	m[1][2] = w * v * (1 - cosf(theta)) + u * sinf(theta);
	m[1][3] = 0;

	m[2][0] = u * w * (1 - cosf(theta)) + v * sinf(theta);
	m[2][1] = v * w * (1 - cosf(theta)) - u * sinf(theta);
	m[2][2] = cosf(theta) + w * w * (1 - cosf(theta));
	m[2][3] = 0;

	m[3][0] = 0;
	m[3][1] = 0;
	m[3][2] = 0;
	m[3][3] = 1;

}


void GlobalFun::Rotate_Point3D(float theta, Point3f normal, Point3f ptIn,Point3f &ptOut ){


	//double nx =normal.X();
	//double ny =normal.Y();
	//double nz =normal.Z();


	//double len = sqrtf(nx * nx + ny * ny + nz * nz); //单位化
	//nx /= len;        ny /= len;            nz /= len;   

	//ptOut[0] =  ptIn[0] * (cosf(theta) + nx * nx * (1 - cosf(theta))) +    //transform by matrix
	//	ptIn[1] * (nx * ny * (1 - cosf(theta)) - nz * sinf(theta)) + 
	//	ptIn[2] * (nx * nz * (1 - cosf(theta) + ny * sinf(theta)));
	//ptOut[1] = ptIn[0] * (nx * ny * (1 - cosf(theta)) + nz * sinf(theta)) +  
	//	ptIn[1] * (ny * ny * (1 - cosf(theta)) + cosf(theta)) + 
	//	ptIn[2] * (ny * nz * (1 - cosf(theta)) - nx * sinf(theta));
	//ptOut[2] = ptIn[0] * (nx * nz * (1 - cosf(theta) - ny * sinf(theta))) + 
	//	ptIn[1] * (ny * nz * (1 -cosf(theta)) + nx * sinf(theta)) + 
	//	ptIn[2] * (nz * nz * (1 - cosf(theta)) + cosf(theta));


	float m[4][4] ;
	RotateArbitraryAxis(m, normal, theta ) ;

	for( int i=0; i<3; ++i ){
		ptOut[i] = 0 ;
		for( int j=0; j<3; ++j )
			ptOut[i] += ptIn[j] * m[j][i] ;
	}


}


void GlobalFun::Rotate_Point3Ds(float theta, Point3f &normal, std::vector<Point3f> &ptIns,std::vector<Point3f>  &ptOuts ){
	ptOuts.clear();
	ptOuts.resize( ptIns.size() ) ;
	for( int i=0; i<ptIns.size(); ++i )
		 GlobalFun::Rotate_Point3D( theta, normal, ptIns[i], ptOuts[i] ) ;


}



bool GlobalFun::broken3dCurve( Curve3D c3,  CurveArray3D &c3s ){
// returned value indicate wether it was broken 
	c3s.clear() ;

	double radius = global_paraMgr.data.getDouble("CGrid Radius");

	bool noGap = true ;
	int n=c3.size();
	int startId = 0 ;
	for( int i=0; i<c3.size(); ++i )
		if( c3[i] == Point3f(1.0e10, 1.0e10, 1.0e10 )  ){
			noGap=false ;
			startId = (i+1)%n ;
			break ;
		}

	// there is no gap
	if( noGap ){
		c3s.push_back(c3) ;
		return false;
	}

	// there is gap
	c3s.resize(1) ;
	for( int i=0; i<n-1; ++i ){
		
		if( c3[ (startId+i)%n ] == Point3f(1.0e10, 1.0e10, 1.0e10 ) && c3s.back().size() != 0  )
			c3s.resize( c3s.size()+1) ;
		else if( c3[ (startId+i)%n ] != Point3f(1.0e10, 1.0e10, 1.0e10 ) ) {
			c3s.back().push_back( c3[ (startId+i)%n ] ) ;
		}
	}


	// remove too short curve
	for( int i=0; i<c3s.size(); ++i ){
		if( c3s[i].size() <= 3  ){
			c3s.erase( c3s.begin() + i) ;

			i-- ;
		}
	}



	return true ;
}



void GlobalFun::smooth3dCurve( CurveArray3D &c3s, bool closed ) {


	int winside = 7; 
	double sigma = 0.3 * ( winside / 2.0 - 1.0 )  + 0.8  ;

	CurveArray3D c3s_bk = c3s ;

	double e = 2.71828 ;
	for( int cid=0; cid<c3s.size(); ++cid ){

		int n = c3s[cid].size() ;
		for( int i=0; i<n; ++i ){


			double totalWeight = 0.0 ; 
			Point3f sum(0,0,0) ;
			for( int offset=-winside/2 ; offset<winside/2; ++offset ){
				if( !closed &&  i+offset >= c3s_bk[cid].size() )
					continue ;

				sum = sum + c3s_bk[cid][ ( i+offset+n) %n ] * pow(e, (-offset*offset/(2*sigma*sigma) ) ) ;
				totalWeight += pow(e, (-offset*offset/(2*sigma*sigma) ) ) ;
			}
			c3s[cid][i] = sum / totalWeight ;

		}
	}

}

Point3f GlobalFun::centerOfPoints( std::vector<Point3f> & points ){

	Point3f center(0,0,0 ) ;
	for( int i=0; i<points.size(); ++i )
		center += points[i] ;
	center /= points.size() ;

	return center;
}


Point2f GlobalFun::centerOfPoints( std::vector<Point2f> & points ){

	Point2f center(0,0 ) ;
	for( int i=0; i<points.size(); ++i )
		center += points[i] ;
	center /= points.size() ;

	return center;
}

void GlobalFun::computePCAplane( std::vector<Point3f> & points, Point3f &vector1, Point3f &vector2, Point3f & center){

	// calculate center
	center = centerOfPoints(points) ;

	// calculate covariance matrix
	Matrix33d covariance_matrix;
	Point3f diff;
	covariance_matrix.SetZero();
	for ( int i=0; i<points.size(); i++)
	{
		Point3f tP = points[i];
		diff = center - tP;

		for (int i=0; i<3; i++)
			for (int j=0; j<3; j++)
				covariance_matrix[i][j] += diff[i]*diff[j];
	
	}

	// calculate eigenvectors
	Point3f   eigenvalues;
	Matrix33d	eigenvectors;
	int required_rotations;
	vcg::Jacobi< Matrix33d, Point3f >(covariance_matrix, eigenvalues, eigenvectors, required_rotations);
	vcg::SortEigenvaluesAndEigenvectors< Matrix33d, Point3f >(eigenvalues, eigenvectors);

	for (int d=0; d<3; d++)
		vector1[d] = eigenvectors[d][0];
	for (int d=0; d<3; d++)
		vector2[d] = eigenvectors[d][1];
	//for (int d=0; d<3; d++)
	//	normal[d] = eigenvectors[d][2];


	vector1.Normalize() ;
	vector2.Normalize() ;
	//normal.Normalize() ;


}

//Point2d mapPoint3ToPoint2( Point3f p, Point3f &vector1, Point3f &vector2, Point3f &center ){  // will not use
//
//	// The formula of plane is Ax + By + Cz + D = 0
//	// compute A,B,C, D
//
//	Point3f normal = vector1^vector2 ;
//	double a = normal[0] ;
//	double b = normal[1] ;
//	double c = normal[2] ;
//	double d = -a*center[0] - b*center[1] - c*center[2] ;
//
//
//}



void GlobalFun::translatePoints(  std::vector<Point3f> & points , Point3f trans ){
	for( int i=0; i<points.size(); ++i)
		points[i] += trans ;
}


extern std::vector<Point3f > globalPointsToDraw1  ;  ;
extern std::vector<Point3f > globalPointsToDraw2  ;  ;


void GlobalFun::computeKNN( const std::vector<Point3f> &dataPoints, const std::vector<Point3f> &queryPoints, std::vector<std::vector<int>> &nearestIds, int nnk ) {
	//nearestIds.clear();
	//nearestIds.resize(queryPoints.size()) ;
	//for( int i=0;i<nearestIds.size(); ++i )
	//	nearestIds[i].push_back(0);
	//return ;

	// ------------------------------------------------------------ build kd-tree for points cloud ------------------------------------------

	int					nPts = dataPoints.size();					// actual number of data points
	ANNpointArray		dataPts;				// data points
	ANNpoint			queryPt;				// query point
	ANNidxArray			nnIdx;					// near neighbor indices
	ANNdistArray		dists;					// near neighbor distances
	ANNkd_tree*			kdTree;					// search structure

	int				k				= nnk;			// number of nearest neighbors
	int				dim				= 3;			// dimension
	double			eps				= 0.001;			// error bound
	int				maxPts			= 300000;			// maximum number of data points

	if (dataPoints.size() >= maxPts)
	{
		cout << "Too many data, fail: "<<__FILE__<<":"<<__LINE__ << endl;
		return;
	}


	queryPt = annAllocPt(dim);					// allocate query point
	dataPts = annAllocPts(nPts, dim);			// allocate data points
	nnIdx = new ANNidx[k];						// allocate near neigh indices
	dists = new ANNdist[k];						// allocate near neighbor dists

	nPts = dataPoints.size();									// read data points
	for( int i=0; i<dataPoints.size(); ++i)
	{
		for(int j = 0; j < 3; j++)
		{
			dataPts[i][j] = double(dataPoints[i][j]); 
		}
	}

	kdTree = new ANNkd_tree(					// build search structure
		dataPts,					// the data points
		nPts,						// number of points
		dim);						// dimension of space

	// ------------------------------------------  kd-tree built ------------------------------------

	nearestIds.clear();

	for( int i=0; i<queryPoints.size(); ++i ){

		for(int j = 0; j < 3; j++){
			queryPt[j] = double(queryPoints[i][j]); 
		}

		kdTree->annkSearch(						// search
			queryPt,						// query point
			k,								// number of near neighbors
			nnIdx,							// nearest neighbors (returned)
			dists,							// distance (returned)
			eps);							// error bound
	
		std::vector<int> tmp ;
		for( int i=0; i<k; ++i)
			tmp.push_back( nnIdx[i] ) ;

		nearestIds.push_back( tmp ) ;

	}

	
	// realease knn stuffs
	delete [] nnIdx;							// clean things up
	delete [] dists;
	annDeallocPt( queryPt ) ;
	annDeallocPts( dataPts ) ;
	delete kdTree;
	//annClose();									// done with ANN
}


void GlobalFun::computeKNN( const std::vector<Point2f> &dataPoints, const std::vector<Point2f> &queryPoints, std::vector<std::vector<int>> &nearestIds, int nnk ) {
	unsigned lastclock = clock() ;

	//nearestIds.clear();
	//nearestIds.resize(queryPoints.size()) ;
	//for( int i=0;i<nearestIds.size(); ++i )
	//	nearestIds[i].push_back(0);
	//return ;
	// ------------------------------------------------------------ build kd-tree for points cloud ------------------------------------------

	int					nPts = dataPoints.size();					// actual number of data points
	ANNpointArray		dataPts;				// data points
	ANNpoint			queryPt;				// query point
	ANNidxArray			nnIdx;					// near neighbor indices
	ANNdistArray		dists;					// near neighbor distances
	ANNkd_tree*			kdTree;					// search structure

	int				k				= nnk;			// number of nearest neighbors
	int				dim				= 2;			// dimension
	double			eps				= 0.001;			// error bound
	int				maxPts			= 300000;			// maximum number of data points

	if (dataPoints.size() >= maxPts)
	{
		cout << "Too many data, fail: "<<__FILE__<<":"<<__LINE__ << endl;
		return;
	}


	queryPt = annAllocPt(dim);					// allocate query point
	dataPts = annAllocPts(nPts, dim);			// allocate data points
	nnIdx = new ANNidx[k];						// allocate near neigh indices
	dists = new ANNdist[k];						// allocate near neighbor dists

	nPts = dataPoints.size();									// read data points
	for( int i=0; i<dataPoints.size(); ++i)
	{
		for(int j = 0; j < dim; j++)
		{
			dataPts[i][j] = double(dataPoints[i][j]); 
		}
	}

	kdTree = new ANNkd_tree(					// build search structure
		dataPts,					// the data points
		nPts,						// number of points
		dim);						// dimension of space

	// ------------------------------------------  kd-tree built ------------------------------------


	//std::cout<<"KNN: kdtree built time - " << (clock()-lastclock)<<"\n" ; lastclock = clock() ;

	nearestIds.clear();

	for( int i=0; i<queryPoints.size(); ++i ){

		for(int j = 0; j < dim; j++){
			queryPt[j] = double(queryPoints[i][j]); 
		}

		kdTree->annkSearch(						// search
			queryPt,						// query point
			k,								// number of near neighbors
			nnIdx,							// nearest neighbors (returned)
			dists,							// distance (returned)
			eps);							// error bound

		std::vector<int> tmp ;
		for( int i=0; i<k; ++i)
			tmp.push_back( nnIdx[i] ) ;

		nearestIds.push_back( tmp ) ;

	}


	// realease knn stuffs
	delete [] nnIdx;							// clean things up
	delete [] dists;
	annDeallocPt( queryPt ) ;
	annDeallocPts( dataPts ) ;
	delete kdTree;
	//annClose();	

	//std::cout<<"KNN: query time - " << (clock()-lastclock) <<std::endl;
	//system("pause") ;
}


void GlobalFun::projectCurveToPointSurface_bynormal( std::vector<Point3f> &surface,std::vector<Point3f> &normals, 
	std::vector<Point3f> &curve0, std::vector<Point3f> &projectedCurve, int nnk = 5) {


		projectedCurve.clear() ;

		std::vector<std::vector<int>> nearestIds ;
		GlobalFun::computeKNN(surface, curve0, nearestIds, nnk ) ;

		int n=curve0.size() ;

		// compute a projection position with each nearest point and its normal
		for( int i=0; i<n;++i ){

			std::vector<Point3f> positions( nnk );

			for( int k=0; k<nnk; ++k ){
				Point3f sp = curve0[i] - surface[ nearestIds[i][k] ];
				Point3f normal = normals[ nearestIds[i][k] ] ;
				double cosa =  normal*sp/(sp.Norm() * normal.Norm())   ;

				Point3f p_dst = - normal.Normalize()  * sp.Norm() * cosa ;

				Point3f dst = curve0[i] + p_dst ;

				positions[k] = dst ;
			}


			// compute means
			Point3f means(0,0,0) ;
			for( int k=0; k<nnk;++k)
				means += positions[k];	
			means /= nnk ;


			// compute variance
			double variance = 0 ;
			for( int k=0; k<nnk;++k)
				variance += (positions[k] - means ).SquaredNorm() ;	
			variance /= nnk  ;
			

			if( sqrt(variance) < 0.01 )
				projectedCurve.push_back( means ) ;


		}


}


void GlobalFun::removeKnot( std::vector<Point3f> &curveIpt, std::vector<Point3f> &afterRemoveKnot ) {

	afterRemoveKnot = curveIpt;

	for( int i=0; i<afterRemoveKnot.size(); ++i ){

		int n = afterRemoveKnot.size();

		Point3f v1 = afterRemoveKnot[i] - afterRemoveKnot[(i-1+n)%n]  ;
		Point3f v2 = afterRemoveKnot[(i+1)%n] - afterRemoveKnot[i]  ;

		if( v1*v2 < 0 ){
			afterRemoveKnot.erase( afterRemoveKnot.begin()+i) ;
			i-- ;
		}

	}

}

Box2f GlobalFun::computeBox( std::vector<Point2f> pts ){

	Box2f box ;
	box.min = Point2f(1e10, 1e10 ) ;
	box.max = Point2f(-1e10, -1e10 ) ;

	for( int i=0; i<pts.size(); ++i ){
		if( pts[i].X() < box.min.X() )  box.min.X() =  pts[i].X()  ;
		if( pts[i].Y() < box.min.Y() )  box.min.Y() =  pts[i].Y()  ;
		if( pts[i].X() > box.max.X() )  box.max.X() =  pts[i].X()  ;
		if( pts[i].Y() > box.max.Y() )  box.max.Y() =  pts[i].Y()  ;
	}

	return box ;
}


Box2f GlobalFun::computeBox( std::vector<std::vector<Point2f>> pts ){

	std::vector<Point2f> ptss ;
	for( int i=0; i<pts.size(); ++i )
		for( int j=0; j<pts[i].size(); ++j )
			ptss.push_back( pts[i][j] ) ;
	return computeBox( ptss ) ;
}


void GlobalFun::draw2dPointsOnScreen( std::vector<Point2f> pts, double3 color, int size ) {

	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	glDisable(GL_LIGHTING) ;

	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);
	int width = viewport[2];
	int height = viewport[3];

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0,width,height,0,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glColor4f(color.r, color.g, color.b, 0.7);

	glPointSize(size) ;

	

	glBegin(GL_POINTS);
	for( int i=0; i< pts.size(); ++ i )
		glVertex2f( pts[i].X(),viewport[3] - pts[i].Y()  ) ;
	glEnd();


	// Closing 2D
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);



	//glEnable(GL_LIGHTING) ;
}



void GlobalFun::draw2dCurveOnScreen( std::vector<Point2f> pts, double3 color, int size, bool closed  ) {

	draw2dCurveOnScreen( pts, GLColor(color.x, color.y, color.z, 1),size, closed ) ;
}

void GlobalFun::draw2dCurveOnScreen( std::vector<Point2f> pts, GLColor color, int size, bool closed  ) {

	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);


	glDisable(GL_LIGHTING) ;

	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);
	int width = viewport[2];
	int height = viewport[3];

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0,width,height,0,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glColor4f(color.r, color.g, color.b, color.a );

	glLineWidth(size) ;
	if( closed )
		glBegin(GL_LINE_LOOP);
	else
		glBegin(GL_LINE_STRIP);

	for( int i=0; i< pts.size(); ++ i )
		glVertex2f( pts[i].X(),viewport[3] - pts[i].Y()  ) ;
	glEnd();


	// Closing 2D
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);



	glEnable(GL_LIGHTING) ;
}


void GlobalFun::draw2dCurveOnScreen_fusiform( std::vector<Point2f> pts, double3 color, int size  ) {

	if( pts.size() < 5 )
		return ;

	//if( (pts[0] - pts[1]).Norm() < 1.0 )
	//	pts.erase( pts.begin() + 1) ;

	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	glDisable(GL_LIGHTING) ;

	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);
	int width = viewport[2];
	int height = viewport[3];

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0,width,height,0,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glColor4f(color.r, color.g, color.b, 0.7);

	glLineWidth(size) ;
	
	std::vector<Point2f>  smCurve ;
	ReconstructorUtility::getUniformCubicBSpline(pts,smCurve,200) ;

	glLineWidth(1.0) ;

	int mp = smCurve.size()/2 ;

	

	glBegin(GL_LINE_STRIP);

	for( int i=0; i< smCurve.size(); ++ i ){
		Point2f p = smCurve[i] ;
		if( i!=0 && i!=smCurve.size()-1 ){
			Point2f dir = smCurve[i+1] - smCurve[i-1] ;
			dir = Point2f( -dir.Y(), dir.X() ) ; dir.Normalize() ;
			p = p + dir * pow(2.718, -(i-mp)*(i-mp) / 50.0/50.0 ) * size  ;
		}
		glVertex2f( p.X(),viewport[3] - p.Y()  ) ;
	}

	for( int i=smCurve.size()-1;  i>=0;  --i ){
		Point2f p = smCurve[i] ;
		if( i!=0 && i!=smCurve.size()-1 ){
			Point2f dir = smCurve[i-1] - smCurve[i+1] ;
			dir = Point2f( -dir.Y(), dir.X() ) ; dir.Normalize() ;
			p = p + dir * pow(2.718, -(i-mp)*(i-mp) / 50.0/50.0  ) * size  ;
		}
		glVertex2f( p.X(),viewport[3] - p.Y()  ) ;
	}

	glEnd();

	//glBegin(GL_LINE_STRIP);

	//for( int i=0; i< smCurve.size(); ++ i )
	//	glVertex2f( smCurve[i].X(),viewport[3] - smCurve[i].Y()  ) ;
	//glEnd();


	// Closing 2D
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);



	glEnable(GL_LIGHTING) ;
}

void GlobalFun::draw3dCurves( 	CurveArray3D &curves3d_, GLColor color, bool closed, bool  randColor, double size , int step) {
	
	if( !curves3d_.size() )
		return ;

	CurveArray3D curves3d ;

	int lastId ;
	for( int i=0; i<curves3d_.size(); i+=step ){
		curves3d.push_back( curves3d_[i] ) ;
		lastId = i ;
	}

	if( (curves3d_.size()-1 -  lastId) > std::max(1, step/2) )
		curves3d.push_back( curves3d_.back() ) ;


	//GLfloat mat_ambient[4] = {0.6, 0.6, 0.6,1.0}; 
	//GLfloat mat_diffuse[4] = {0.0, 0.6, 0.6, 1.0 };
	//GLfloat mat_specular[] = {0, 1.0, 1.0, 1.0 };
	GLfloat mat_ambient[4] = {0.6, 0.6, 0.6,1.0}; 
	GLfloat mat_diffuse[4] = {0.6, 0.6, 0.6, 1.0 };
	GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0 };
	GLfloat shininess = 0.5*128;

	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse); 
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
	
	glEnable(GL_COLOR_MATERIAL) ;
	

	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	//glDisable( GL_BLEND ) ;
	
	glEnable(GL_LIGHTING) ;
	glEnable(GL_LIGHT0) ;
	glEnable(GL_LIGHT1) ;



	GLUquadricObj *objCylinder = gluNewQuadric();


	for( int i=0; i< curves3d.size(); ++ i ){

		if( randColor ){
			//color = glDrawer.randColor[i%32] ;
		}

		glColor4f( color.r, color.g, color.b, color.a);

		double radius = 0.002 * size ;
		int n = curves3d[i].size() ;
		for( int id = 0; id < curves3d[i].size(); ++id ){


			if( !closed && id == curves3d[i].size()-1 )
				break  ;
			//glLineWidth( 8.0 ) ;

			//glBegin(GL_LINES);
			//glVertex3f( curves3d[i][id].X(), curves3d[i][id].Y() , curves3d[i][id].Z()  ) ;
			//glVertex3f( curves3d[i][(id+1)%n ].X(),curves3d[i][(id+1)%n].Y(),curves3d[i][(id+1)%n].Z()  ) ;
			//glEnd();
			//glPointSize( 10 / 2) ;
			//if( id != 0){
			//	glBegin(GL_POINTS);
			//	glVertex3f( curves3d[i][id].X(), curves3d[i][id].Y() , curves3d[i][id].Z()  ) ;
			//	glEnd();

			//}


			Point3f p0 = curves3d[i][(id-1+n)%n ] ;
			Point3f p1 = curves3d[i][(id)%n ] ;
			Point3f p2 = curves3d[i][(id+1)%n ] ;
			Point3f p3 = curves3d[i][(id+2)%n ] ;

			Point3f dir = p2 - p1 ;

			Point3f zaxis(0,0,1) ;
			Point3f rAxis = zaxis ^  dir ;
			double angle = acos( zaxis*dir / dir.Norm() ) *180/3.14159;

			glPushMatrix() ;
			glTranslatef(curves3d[i][id] .X(), curves3d[i][id].Y(), curves3d[i][id].Z() ) ;
			glRotatef( angle, rAxis.X(), rAxis.Y(), rAxis.Z() ) ;
			gluCylinder(objCylinder, radius, radius, dir.Norm(), 20, 2);
			//gluDisk(objCylinder, 0, radius*0.99, 16, 3 ) ;
			glPopMatrix() ;


			//glPushMatrix() ;
			//glTranslatef(curves3d[i][id] .X(), curves3d[i][id].Y(), curves3d[i][id].Z() ) ;
			//glRotatef( angle, rAxis.X(), rAxis.Y(), rAxis.Z() ) ;
			//glTranslatef(0, 0,  dir.Norm()) ;
			////gluDisk(objCylinder, 0, radius*0.99, 32, 10 ) ;
			//glPopMatrix() ;

			//glPushMatrix() ;
			//glTranslatef(curves3d[i][id] .X(), curves3d[i][id].Y(), curves3d[i][id].Z() ) ;
			//gluSphere(objCylinder, radius*1.00, 20, 10 );
			//glPopMatrix() ;



		}

	}

	gluDeleteQuadric( objCylinder ) ;

	GLfloat mat_ambient_[4] = {0.2, 0.2, 0.2,1.0}; 
	GLfloat mat_diffuse_[4] = {0.8, 0.8, 0.8, 1.0 };
	GLfloat mat_specular_[] = {0.0, 0.0, 0.0, 1.0 };
	GLfloat shininess_ = 0;

	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_ambient_);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse_); 
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular_);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess_);

	glDisable(GL_LIGHT0) ;
	glDisable(GL_LIGHT1) ;
	glDisable(GL_LIGHTING) ;

}



void GlobalFun::draw3dCurves_nomat( 	CurveArray3D &curves3d, GLColor color, bool closed, bool  randColor, double size ) {



	glEnable(GL_LIGHTING) ;
	glEnable(GL_LIGHT0) ;
	glEnable(GL_LIGHT1) ;


	GLUquadricObj *objCylinder = gluNewQuadric();


	for( int i=0; i< curves3d.size(); ++ i ){

		if( randColor ){
			//color = glDrawer.randColor[i%32] ;
		}

		glColor4f( color.r, color.g, color.b, color.a);

		double radius = 0.002 * size ;
		int n = curves3d[i].size() ;
		for( int id = 0; id < curves3d[i].size(); ++id ){


			if( !closed && id == curves3d[i].size()-1 )
				break  ;


			Point3f p0 = curves3d[i][(id-1+n)%n ] ;
			Point3f p1 = curves3d[i][(id)%n ] ;
			Point3f p2 = curves3d[i][(id+1)%n ] ;
			Point3f p3 = curves3d[i][(id+2)%n ] ;

			Point3f dir = p2 - p1 ;

			Point3f zaxis(0,0,1) ;
			Point3f rAxis = zaxis ^  dir ;
			double angle = acos( zaxis*dir / dir.Norm() ) *180/3.14159;

			glPushMatrix() ;
			glTranslatef(curves3d[i][id] .X(), curves3d[i][id].Y(), curves3d[i][id].Z() ) ;
			glRotatef( angle, rAxis.X(), rAxis.Y(), rAxis.Z() ) ;
			gluCylinder(objCylinder, radius, radius, dir.Norm(), 20, 2);
			//gluDisk(objCylinder, 0, radius*0.99, 16, 3 ) ;
			glPopMatrix() ;


			//glPushMatrix() ;
			//glTranslatef(curves3d[i][id] .X(), curves3d[i][id].Y(), curves3d[i][id].Z() ) ;
			//glRotatef( angle, rAxis.X(), rAxis.Y(), rAxis.Z() ) ;
			//glTranslatef(0, 0,  dir.Norm()) ;
			////gluDisk(objCylinder, 0, radius*0.99, 32, 10 ) ;
			//glPopMatrix() ;

			//glPushMatrix() ;
			//glTranslatef(curves3d[i][id] .X(), curves3d[i][id].Y(), curves3d[i][id].Z() ) ;
			//gluSphere(objCylinder, radius*1.00, 20, 10 );
			//glPopMatrix() ;



		}

	}

	gluDeleteQuadric( objCylinder ) ;


}
//  [9/2/2014 CX]
bool GlobalFun::mycompfunc (std::pair<double,int> a,std::pair<double,int> b) {    return ( a.first < b.first ) ; }
bool GlobalFun::myintcompfunc (std::pair<int,int> a,std::pair<int,int> b) {    return ( a.first < b.first ) ; }
bool GlobalFun::myint2compfunc1 (std::pair<int,int> a,std::pair<int,int> b) {    return ( a.first < b.first ) ; }
bool GlobalFun::myint2compfunc2 (std::pair<int,int> a,std::pair<int,int> b) {    return ( a.second < b.second ) ; }