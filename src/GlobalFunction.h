

#pragma once
#include "Types.h"
#include <vector>
#include "CMesh.h"
#include <fstream>
#include <float.h>
#include <QString>
#include <iostream>
#include <time.h>
#include <string>
#include <ctime>
#include<algorithm>
#include <math.h>
#include "ANN/ANN.h"
#include "DataMgr.h"
#include "GLDrawer.h"

#define EIGEN_DEFAULT_TO_ROW_MAJOR
#define EIGEN_EXCEPTIONS

//#include <Eigen/Dense>

//class EllipseShape ;
//class Skeletonization ;
class GLColor;

using namespace std;
using namespace vcg;


//typedef Eigen::MatrixXd Matrix;

namespace GlobalFun
{
	
	double computeEulerDist(Point3f& p1, Point3f& p2);
	double computeEulerDistSquare(Point3f& p1, Point3f& p2);
	double computeProjDist(Point3f& p1, Point3f& p2, Point3f& normal_of_p1);
	double computeProjDistSquare(Point3f& p1, Point3f& p2, Point3f& normal_of_p1);
	double computePerpendicularDistSquare(Point3f& p1, Point3f& p2, Point3f& normal_of_p1);
	double computePerpendicularDist(Point3f& p1, Point3f& p2, Point3f& normal_of_p1);
	vector<int> GetRandomCards(int Max);

	double computeRealAngleOfTwoVertor(Point3f v0, Point3f v1);
	bool isTwoPoint3fTheSame(Point3f& v0, Point3f& v1);
	bool isTwoPoint3fOpposite(Point3f& v0, Point3f& v1);


	bool gluInvertMatrix(const float m[16], float invOut[16]) ;

	Point3f vecMulMatrix( Point3f point0, float invmvMatrix[16] ) ;

	bool PtInPolygon (Point2f x, vector<Point2f> &y) ;

	double3 angle2color( double angle ) ;//

	void translatePoints(  vector<Point3f> & points , Point3f trans ) ;
	void Rotate_Point3D(float theta, Point3f normal, Point3f ptIn,Point3f &ptOut ) ;
	void Rotate_Point3Ds(float theta, Point3f &normal, vector<Point3f> &ptIns,vector<Point3f>  &ptOuts ) ;

	bool broken3dCurve( vector<Point3f> c3,  vector<vector<Point3f>> &c3s );

	void smooth3dCurve(  vector<vector<Point3f>> &c3s, bool closed = true) ;
	double3 scalar2color( double scalar ) ;


	Point3f centerOfPoints( vector<Point3f> & points );
	Point2f centerOfPoints( std::vector<Point2f> & points ) ;

	void computePCAplane( vector<Point3f> & points, Point3f &vector1, Point3f &vector2, Point3f & center) ;


	void computeKNN( const vector<Point3f> &dataPts, const vector<Point3f> &queryPts, vector<vector<int>> &nearestIds, int nnk ) ;
	void computeKNN( const vector<Point2f> &dataPts,const  vector<Point2f> &queryPts, vector<vector<int>> &nearestIds, int nnk ) ;

	void projectCurveToPointSurface_bynormal( vector<Point3f> &surface,vector<Point3f> &normals, vector<Point3f> &curve0, vector<Point3f> &projectedCurve, int nnk ) ;
	void removeKnot( vector<Point3f> &afterProject, vector<Point3f> &afterRemoveKnot ) ;

	void draw2dPointsOnScreen( vector<Point2f> pts, double3 color, int size = 3 ) ;
	void draw2dCurveOnScreen( vector<Point2f> pts, double3 color, int size = 3 , bool closed = true) ;
	void draw2dCurveOnScreen( std::vector<Point2f> pts, GLColor color, int size, bool closed  = true ) ;

	void draw2dCurveOnScreen_fusiform( std::vector<Point2f> pts, double3 color, int size  );

	//void draw3dCurves( 	CurveArray3D &curves3d, GLColor color, bool closed, bool  randColor, double size ) ;
	void draw3dCurves( 	CurveArray3D &curves3d, GLColor color, bool closed, bool  randColor, double size , int step = 1) ;
	void draw3dCurves_nomat( 	CurveArray3D &curves3d, GLColor color, bool closed, bool  randColor, double size );

	Box2f computeBox( vector<Point2f> pts ) ;
	Box2f computeBox( vector<vector<Point2f>> pts );
	//  [9/2/2014 CX]
	bool mycompfunc (std::pair<double,int> a,std::pair<double,int> b);
	bool myintcompfunc (std::pair<int,int> a,std::pair<int,int> b);
	bool myint2compfunc1 (std::pair<int,int> a,std::pair<int,int> b);
	bool myint2compfunc2 (std::pair<int,int> a,std::pair<int,int> b);
} ;

class Timer
{
public:

	void start(const string& str)
	{
		cout << endl;
		starttime = clock();
		mid_start = clock();
		cout << "@@@@@ Time Count Strat For: " << str << endl;

		_str = str;
	}

	void insert(const string& str)
	{
		mid_end = clock();
		timeused = mid_end - mid_start;
		cout << "##" << str << "  time used:  " << timeused / double(CLOCKS_PER_SEC) << " seconds." << endl;
		mid_start = clock();
	}

	void end()
	{
		stoptime = clock();
		timeused = stoptime - starttime;
		cout << /*endl <<*/ "@@@@ finish	" << _str << "  time used:  " << timeused / double(CLOCKS_PER_SEC) << " seconds." << endl;
		cout << endl;
	}

private:
	int starttime, mid_start, mid_end, stoptime, timeused;
	string _str;
};




/* Useful code template

(1)
for(int i = 0; i < samples->vert.size(); i++)
{
CVertex& v = samples->vert[i];

for (int j = 0; j < v.neighbors.size(); j++)
{
CVertex& t = samples->vert[v.neighbors[j]];
}
}

(2)
int count = 0;
time.start("Test 2");
CMesh::VertexIterator vi;
Point3f p0 = Point3f(0,0,0);
for(vi = original->vert.begin(); vi != original->vert.end(); ++vi)
{
count += GlobalFun::computeEulerDistSquare(p0, vi->P());
}
cout << count << endl;
time.end();


time.start("Test 1");
for(int i = 0; i < original->vert.size(); i++)
{
CVertex& v = original->vert[i];
count += (p0 - v.P()).SquaredNorm();
}
cout << count << endl;
time.end();



*/

