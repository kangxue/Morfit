/*
	This file is part of the software described in the paper:
	Kangxue Yin, Hui Huang, Hao(Richard) Zhang, Minglun Gong, Daniel Cohen-or, Baoquan Chen. 
	Morfit: Interactive Surface Reconstruction from Incomplete Point Clouds with Curve-Driven Topology and Geometry Control. 
	ACM Transactions on Graphics 33(6)(Proc.SIGGRAPH ASIA 2014 )

	Copyright (C) <2013-2014>  <Kangxue Yin - yinkangxue@gmail.com>
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#ifndef _ReconstructorUtility_h_
#define  _ReconstructorUtility_h_

#include <pcl/surface/on_nurbs/fitting_curve_2d_pdm.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_tdm.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_sdm.h>
#include <pcl/surface/on_nurbs/triangulation.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>

//#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/surface/3rdparty/opennurbs/opennurbs_nurbscurve.h>

#include "smtf.h"
#include "skeleton_mul.h"
#include "skelpath.h"
//#include "..\nurbs++-3.0.11\nurbs\nurbs.h"

namespace ReconstructorUtility{

	template <class _Point>

	double vectorIntersectionAngle( _Point p1, _Point p2 ) ;

	std::vector<double> convertTransform2Vector ( const std::vector<ST> &stf ) ;

	std::vector<ST> convertVector2Transform (const std::vector<double> &x0 ) ;

	std::vector<Point2f> transformProfile( const std::vector<Point2f> &input ,  ST stf) ;

	std::vector<std::vector<Point2f>> transformProfiles( const std::vector<Point2f> &input ,  std::vector<ST> &stf) ;


	Profile2D mergeProfile2d( const Profile2D &temProf, const Profile2D &dataProf );

	Profile3D transformProfile3D( Profile3D input ,  ST stf, cutPlane cutpl_src, Point3f cter_src, cutPlane cutpl_dst, Point3f cter_dst ) ;

	Profile3D convert2dProfileTo3d( const Profile2D &prof2d, cutPlane cutpl, Point3f center ) ;
	Profile2D convert3dProfileTo2d( const Profile3D &prof3d, cutPlane cutpl, Point3f center ) ;
	ProfilePL convert2dProfileToPL( const Profile2D &prof2d );
	Profile2D convertPLProfileTo2d(   ProfilePL &profPL) ;
	std::vector<double> convertPLProfileTo360(   ProfilePL &profPL) ;

	//fit the nubs 
	ON_NurbsCurve fitNURBSforProfile2d( Profile2D &prof );

	ON_NurbsCurve fitAndUniformNURBSforProfile2d( Profile2D &prof );

	std::vector<Point2f> discretizeNurbs( ON_NurbsCurve nbs, int sampleNum, bool fast =false  ) ;

	ON_NurbsCurve deformNurbs( const ON_NurbsCurve nbs_in, NRTSF ngTrans, const std::vector<double> & weights ) ;

	double profileDis( const Profile2D &prof1,  const Profile2D &prof2 ) ;

	double profileDis_d( Profile2D prof1,  Profile2D prof2 ) ;
	double profileDis_d( Profile3D prof1,  Profile3D prof2 ) ;
	double ProfileDis_force(Profile2D prof1,  Profile2D prof2) ;
	double profileDis_360deg( Profile2D prof1,  Profile2D prof2 ) ;

	std::vector<double> convertNRTransform2Vector ( const std::vector<NRTSF> &nrtsf ) ;
	std::vector<NRTSF> convertVector2NRTransform (const std::vector<double> &x0, int cvNum ) ;

	std::vector<int> getSharpCVofNurbs(ON_NurbsCurve tmp) ;

	std::vector<int2> detectCorrespondenceForTmps( ON_NurbsCurve T1, ON_NurbsCurve T2, std::vector<int> sharpId1, std::vector<int> sharpId2  ) ;

	void saveTmpOf1Brch( std::vector<ON_NurbsCurve> templates, std::string filename ) ;
	void loadTmpOf1Brch( std::vector<ON_NurbsCurve> &templates, std::string filename ) ;
	void saveTemplates(  std::vector<std::vector<ON_NurbsCurve>> templates,  std::vector<std::vector<int>> templateIdss,  std::string dir ) ;
	void loadTemplates(  std::vector<std::vector<ON_NurbsCurve>> &templates,  std::vector<std::vector<int>> &templateIdss,  std::string dir ) ;

	std::vector<Point3f> upsampleResult( std::vector<Profile3D> profiles, int2 scope, bool firstIsTip, bool lastIsTip  ) ;

	bool vectorIntersect( Point2f p0, Point2f p1, Point2f p2, Point2f p3 ) ;

	void getBezier( std::vector<Point3f> ctrpoints,  std::vector<Point3f>  &output, int nump ) ;
	void getBezier( std::vector<Point2f> ctrpoints,  std::vector<Point2f>  &output, int nump ) ;

	void getUniformCubicBSpline( std::vector<Point3f> input, std::vector<Point3f> &output, int nump ) ;
	void getUniformCubicBSpline( std::vector<Point2f> input, std::vector<Point2f> &output, int nump ) ;

	void getNurbs( std::vector<Point3f> ctrpoints,  std::vector<Point3f>  &output, int nump, int degree ) ;
	void getPeriodicNurbs( std::vector<Point3f> ctrpoints,  std::vector<Point3f>  &output, int nump, int degree ) ;

	std::vector<Point3f>  cutPointsSet( std::vector<Point3f> points, Point3f center, Point3f normal ) ;
	
	std::vector<std::vector<Point3f>> segmentPointset( std::vector<Point3f> points, std::vector<std::vector<Point3f>> brches ) ;

	std::vector<Point3f> getSlice( std::vector<Point3f> &pointset, cutPlane cutp, Point3f cter) ;

	double evaluateConfident( Profile2D &points );

	double round(double r);

	std::vector<Point3f> getHermite(  Point3f p0_, Point3f p1_, Point3f tan0_, Point3f tan1_, double sampleStep ) ;
	Point3f get_hermite_value( Point3f p0_, Point3f p1_, Point3f tan0_, Point3f tan1_, double t ) ;

	void changeCoorOfTmps( ON_NurbsCurve &tmp, cutPlane srcpl, cutPlane dstpl, Point3f srccter, Point3f dstcter  ) ;


	void upsampleTip( std::vector<Point3f> &res, Point3f tip, std::vector<Point3f> neiprof,  std::vector<Point3f> neiprof_full, std::vector<Point3f> neineiprof_full ) ;

	std::vector<Point3f> convertLocal3dToWorld3d( std::vector<Point3f> points, float *mvmatrix ) ;
	std::vector<Point3f> convertLocal3dToWorld3d( std::vector<Point3f> points ) ;

	std::vector<Point3f> convertWorld3dToLocal3d( std::vector<Point3f> points, float *mvmatrix ) ;
	std::vector<Point3f> convertWorld3dToLocal3d( std::vector<Point3f> points ) ;
	
	std::vector<Point2f> project3dPointsOntoScreen( std::vector<Point3f> points, float *mvmatrix ) ;
	std::vector<Point2f> project3dPointsOntoScreen( std::vector<Point3f> points ) ;

	std::vector<Point2f> convertScreen2dToWorld2d( std::vector<Point2f> curve ) ;
	std::vector<Point2f> convertWorld2dToScreen2d( std::vector<Point2f> curve ) ;
	std::vector<Point3f> convertScreen2dToWorld3d( std::vector<Point2f> curve ) ;
	std::vector<Point3f> convertScreen2dToLocal3d( std::vector<Point2f> curve );

	int NearstPoint( std::vector<Point3f> &datapts, Point3f pt) ;
	int NearstPoint( std::vector<Point2f> &datapts, Point2f pt) ;

	void computeFlannKNN( const std::vector<Point2f> &dataPoints, const std::vector<Point2f> &queryPoints, std::vector<std::vector<int>> &nearestIds, int nnk ) ;
	void computeFlannKNN( const std::vector<Point3f> &dataPoints, const std::vector<Point3f> &queryPoints, std::vector<std::vector<int>> &nearestIds, int nnk ) ;
	void MycomputeFlannKNN(std::vector<Point3f> &dataPoints,std::vector<Point3f> &queryPoints,std::vector<Point3f> &nmls,std::vector<int> &nearestId);

	Point3f ComputeEigenvector( std::vector< Point3f > nmls);

	std::vector<Profile2D> cvtProf2Pixels( std::vector<Profile2D> profiles, Box2f box ) ;
	std::vector<Profile2D>  cvtPixels2Prof( std::vector<Profile2D> pixels, Box2f box  ) ;

	std::vector<Point2f> getCvPts(ON_NurbsCurve nbs) ;

	void setCvOfNurbs ( ON_NurbsCurve &nbs, std::vector<Point2f> CVs ) ;
	void setCvOfNurbs ( ON_NurbsCurve &nbs, std::vector<Point2f> CVs, std::vector<double> Weights ) ;

	bool isTip( skelpath skel, Point3f p ) ;

	bool isMidNode( std::vector<Curve3D>  impSkel, Point3f p ) ;

	bool CalPlaneLineIntersectPoint(Point3f &interPoint, Point3f planeVector, Point3f planePoint, Point3f lineVector, Point3f linePoint  ) ;

	void smoothVectorData( std::vector<double> &v, int winsize ) ;

	void CopyCmesh( CMesh & dst, CMesh &src ) ;

	void drawPointsAsSphere( std::vector<Point3f> pts, GLColor color, int size ) ;
	void draw3dCurves( 	Curve3D &curves3d, GLColor color, bool closed, bool  randColor, double size ) ;
	bool loadSkeletonAndProfiles( skelpath &skel, int ogSkelId ) ;

	void updateSkelDepth( skeleton_mul &src, skeleton_mul dst ) ;
	void updateSkelDepth( skelpath &src, skelpath dst ) ;

	template< class T>
	bool inVerctor( std::vector<T> vec, T ele){
		for( int i=0; i<vec.size(); ++i )
			if( vec[i] == ele )
				return true ;
		return false ;
	}

	template< class T>
	int FindInVerctor( std::vector<T> vec, T ele){
		for( int i=0; i<vec.size(); ++i )
			if( vec[i] == ele )
				return i ;
		return -1 ;
	}

	template <class T>
	int getMinId( std::vector<T> vals ) {

		T minValue = vals[0];
		int minId = 0 ;

		for( int i=0; i<vals.size(); ++i ){
			if( vals[i] < minValue ){
				minValue = vals[i]  ;
				minId = i ;
			}
		}
		return minId ;
	}


}


#endif