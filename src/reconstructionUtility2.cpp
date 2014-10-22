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



#include "GlobalFunction.h"
#include "reconstructionUtility.h"
#include "smtf.h"
#include "skeleton_mul.h"
#include "reconstructorPara.h"
//#include "declarations.h"
#include "appstate.h"
#include "gl/GL.h"
#include<numeric>

extern GLdouble globalModelviewMatrix[16] ;

extern double glareaFov  ;

namespace ReconstructorUtility{

	std::vector<Point3f>  cutPointsSet( std::vector<Point3f> points, Point3f center, Point3f normal ) {
	// preserve points on positive direction	
		std::vector<Point3f> res ;
		for( int i=0; i<points.size(); ++i )
			if( (points[i]-center) * normal >= 0 )
				res.push_back( points[i]) ;

		return res ;
	}


	std::vector<std::vector<Point3f>> segmentPointset( std::vector<Point3f> points, std::vector<std::vector<Point3f>> brches ) {

		//divide the point cloud by the branched
		std::vector<std::vector<Point3f>> & smoothedBrchPts = brches ;

		std::vector<Point3f> branchPoints ;
		std::vector<int2> index ;
		for( int i=0; i<smoothedBrchPts.size(); ++i ){
			for( int j=0; j<smoothedBrchPts[i].size(); ++j ){
				branchPoints.push_back(smoothedBrchPts[i][j] ) ;
				index.push_back( int2(i,j) ) ;
			}
		}

		std::vector<std::vector<int>> nearestNeibs ;
		GlobalFun::computeKNN( branchPoints,points, nearestNeibs, 1 ) ;

		std::vector<std::vector<Point3f>>  fragments(brches.size()) ;
		for( int i=0; i<nearestNeibs.size(); ++i )
			fragments[ index[nearestNeibs[i][0]].x ].push_back( points[i]) ;

		return fragments ;
	}


	std::vector<Point3f> getSlice( std::vector<Point3f> &pointset, cutPlane cutp, Point3f cter){
		//get the slice of the branch ,
		//the center is the third parameter , the cut plane at the center is the second parameter
		double left_scope = ReconstructorPara::branchSampleStep / 2 ;
		double right_scope = ReconstructorPara::branchSampleStep / 2 ;

		std::vector<Point3f> slice  ;
		for( int pid =0; pid<pointset.size(); ++pid ){
			double r = (pointset[pid]-cter).Norm() ;
			double dis = (pointset[pid]-cter) * cutp.planeNormal  ;
			if(  dis > -left_scope && dis <right_scope   )
				slice.push_back( pointset[pid] ) ;
		}

		return slice ;

	}

	double evaluateConfident( Profile2D &points ){
	
		// evaluate confidence with the distance of neighbor points
		if( points.size() < 5 )
			return 0 ;

		// evaluate confidence with max span degree
		std::vector<std::pair<double,int>> angles ; 
		Point2f axisy = Point2f(0.0, 1.0) ;

		for( int pid=0; pid<points.size(); ++pid ){

			double alpha = ReconstructorUtility::vectorIntersectionAngle( points[pid], axisy ) ;

			if( points[pid][0] < 0.0 )
				alpha = 3.14159*2 - alpha  ;
			angles.push_back(std::pair<double,int>(alpha, pid) ) ;
		}

		bool GlobalFun::mycompfunc (std::pair<double,int> a,std::pair<double,int> b) ;
		std::sort( angles.begin(), angles.end(), GlobalFun::mycompfunc ) ;

		double sum_span_degree=0.0;
		std::vector<double> span_degree ;
		for( int i=1; i<angles.size(); ++i )
			span_degree.push_back(angles[i].first - angles[i-1].first);
		span_degree.push_back(3.14159 * 2 - ( angles.back().first-angles[0].first ) );
		std::sort( span_degree.begin(), span_degree.end()) ;
		for(int i=1;i<=5;i++)
			sum_span_degree += span_degree[span_degree.size()-i];
		//sum_span_degree /= span_degree.size();
		return 1/sum_span_degree;

	}


	ProfilePL convert2dProfileToPL( const Profile2D &prof2d ) {

		ProfilePL ppl ;
		Point2f axisy = Point2f(0.0, 1.0) ;
		for( int pid=0; pid<prof2d.size(); ++pid ){

			double alpha = ReconstructorUtility::vectorIntersectionAngle( prof2d[pid], axisy ) ;
			if( prof2d[pid][0] < 0.0 )
				alpha = 3.14159*2 - alpha  ;

			ppl.push_back( PointPL(alpha, prof2d[pid].Norm() ) ) ;
		}

		return ppl ;
	}



	Profile2D convertPLProfileTo2d(   ProfilePL &profPL) {

		Profile2D prof2d ;

		Point2f axisy = Point2f(0.0, 1.0) ;
		for( int pid=0; pid<profPL.size(); ++pid ){


			Point2f axy0 = axisy * profPL[pid].rd() ;
			double theta = - profPL[pid].th() ;

			double x = axy0.X() * cos(theta) - axy0.Y() * sin(theta) ;
			double y = axy0.X() * sin(theta) + axy0.Y() * cos(theta) ;

			prof2d.push_back( Point2f(x,y) ) ;

		}

		return prof2d ;
	}



	std::vector<double> convertPLProfileTo360(   ProfilePL &profPL) {

		std::vector<std::vector<double>> radii ( 360 ) ;

		for( int pid=0; pid<profPL.size(); ++pid ){
			int d = round( profPL[pid].th() * 180.0 / 3.14159 ) ;
			d = (d+360)%360 ;
			radii[d].push_back( profPL[pid].rd() ) ;
		}

		
		std::vector<double> radius(360,0) ;

		for( int d=0; d<360; ++d ){

			if( radii[d].size()!=0 ){
				for( int i=0; i<radii[d].size(); ++i)
					radius[d] += radii[d][i] ;
				radius[d] /= radii[d].size();
				continue;
			}

			// then radii[d] is empty


			int id0= d;
			int id1 = d ;
			double r0 = 0 ;
			double r1 = 0 ;

			// search left
			for( int ofs=0; ofs<180; ++ofs ){
				int id = (d - ofs + 360)%360 ;
				if( radii[id].size() ){
					id0 = id ;
					for( int i=0; i<radii[id].size(); ++i)
						r0 += radii[id][i] ;
					r0 /= radii[id].size();
					break;
				}
			}

			// search right
			for( int ofs=0; ofs<180; ++ofs ){
				int id = (d + ofs + 360)%360 ;
				if( radii[id].size() ){
					id1 = id ;
					for( int i=0; i<radii[id].size(); ++i)
						r1 += radii[id][i] ;
					r1 /= radii[id].size();
					break;
				}
			}


			// interplate radius linearly
			radius[d] = (r1-r0) * (d - id0) / (id1-id0)  + r0;

		}
		return radius ;

	}


	std::vector<Point3f> getHermite(  Point3f p0_, Point3f p1_, Point3f tan0_, Point3f tan1_, double sampleStep ){

		tan0_.Normalize() ;
		tan1_.Normalize() ;

		std::vector<Point3f> hermite ;
		for( int i=0; i<= 1000; ++i ){
			double t = i/1000.0 ;
			hermite.push_back( get_hermite_value( p0_, p1_, tan0_, tan1_, t ) ) ;
		}


		double length = 0;
		for( int i=1; i<hermite.size(); ++i)
			length += ( hermite[i] - hermite[i-1]).Norm() ; 

		double tstep = sampleStep/length ;

		hermite.clear() ;
		for( double t=0; t<1.0; t+= tstep ){
			hermite.push_back( get_hermite_value( p0_, p1_, tan0_, tan1_, t ) ) ;
		}
		hermite.push_back( get_hermite_value( p0_, p1_, tan0_, tan1_, 1.0 ) ) ;


		return hermite ;

	}


	Point3f get_hermite_value( Point3f p0_, Point3f p1_, Point3f tan0_, Point3f tan1_, double t ){


		double3 p0 = double3( p0_.X(), p0_.Y() , p0_.Z() ) ;
		double3 p1 = double3( p1_.X(), p1_.Y(), p1_.Z() ) ;
		double3 tan0 = double3( tan0_.X(), tan0_.Y() , tan0_.Z() ) ;
		double3 tan1 = double3( tan1_.X(), tan1_.Y() , tan1_.Z() ) ;


		tan0.normalize() ;
		tan1.normalize() ;

		double dis = (p0-p1).norm();

		tan0 = tan0 * dis;
		tan1 = tan1 * dis ;



		double m[4][4] = { {2,-2,1,1},{-3,3,-2,-1},{0,0,1,0},{1,0,0,0} } ;

		double T[4] ;
		double temp[4] ;
		double3 result ;

		T[0] = t*t*t ;
		T[1] = t*t;
		T[2] = t ;
		T[3] = 1 ;
		for( int x = 0; x<4; ++x ){
			temp[x]  = T[0] * m[0][x] ;
			temp[x] += T[1] * m[1][x] ;
			temp[x] += T[2] * m[2][x] ;
			temp[x] += T[3] * m[3][x] ;
		}

		result.x = temp[0] * p0.x + temp[1] * p1.x + temp[2] * tan0.x + temp[3] * tan1.x ;
		result.y = temp[0] * p0.y + temp[1] * p1.y + temp[2] * tan0.y + temp[3] * tan1.y ;
		result.z = temp[0] * p0.z + temp[1] * p1.z + temp[2] * tan0.z + temp[3] * tan1.z ;

		return Point3f(result.x, result.y, result.z ) ;

	}


	//void morphingBetween2BoundaryProfiles( on_nurbs nbs1,on_nurbs nbs2, std::vector<Point3f> Brch1,std::vector<Point3f> Brch2,  ) {


	//}

	std::vector<Point2f> getCvOfNurbs( ON_NurbsCurve nbs ){

		int cvnum = nbs.CVCount() - nbs.Degree() ;

		std::vector<Point2f> CVs ;
		for( int i=0; i<cvnum; ++i ){

			ON_4dPoint p;
			nbs.GetCV( i, p) ; p.x/=p.w; p.y/=p.w ;
			CVs.push_back( Point2f(p.x, p.y) ) ;
		}

		return CVs ;
	}

	std::vector<double> getCVWeightsOfNurbs( ON_NurbsCurve nbs ) {
	
		int cvnum = nbs.CVCount() - nbs.Degree() ;

		std::vector<double> ws ;
		for( int i=0; i<cvnum; ++i ){

			ON_4dPoint p;
			nbs.GetCV( i, p) ; p.x/=p.w; p.y/=p.w ;
			ws.push_back( p.w ) ;
		}

		return ws ;
	}


	void setCvOfNurbs ( ON_NurbsCurve &nbs, std::vector<Point2f> CVs ){

		int cvnum = nbs.CVCount() - nbs.Degree() ;

		for( int i=0; i<cvnum; ++i ){

			ON_4dPoint p;
			nbs.GetCV( i, p) ; 

			p.x = CVs[i].X() * p.w ;
			p.y = CVs[i].Y() * p.w ;

			nbs.SetCV( i, p) ;


			if( i< nbs.Degree())
				nbs.SetCV( cvnum+i,p ) ;
		}

	}

	void setCvOfNurbs ( ON_NurbsCurve &nbs, std::vector<Point2f> CVs, std::vector<double> Weights ){

		int cvnum = nbs.CVCount() - nbs.Degree() ;

		for( int i=0; i<cvnum; ++i ){

			ON_4dPoint p;

			p.x = CVs[i].X() * Weights[i] ;
			p.y = CVs[i].Y() * Weights[i]  ;
			p.z = 0;
			p.w = Weights[i] ;

			nbs.SetCV( i, p) ;

			if( i< nbs.Degree())
				nbs.SetCV( cvnum+i,p ) ;
		}

	}

	std::vector<Point2f> rotate2dPoints( std::vector<Point2f> input, Point2f center, double angle ){

		std::vector<Point2f> output ;

		for( int i=0; i<input.size(); ++i ){

			Point2f p = input[i] - center ;

			// p = Ap

			double x = p.X() * cos(angle) - p.Y() * sin(angle) ;
			double y = p.X() * sin(angle) + p.Y() * cos(angle) ;

			output.push_back( Point2f(x,y) + center) ;
		}
		return output ;
	}


	void changeCoorOfTmps( ON_NurbsCurve &tmp, cutPlane srcpl, cutPlane dstpl, Point3f srccter, Point3f dstcter  ) {
		
		//change the coordinate of the template
		int cvnum = tmp.CVCount() - tmp.Degree() ;

		std::vector<Point2f> CVs = getCvOfNurbs( tmp ) ;
		std::vector<double> Ws = getCVWeightsOfNurbs( tmp ) ;

		std::vector<Point3f> CVs3d = convert2dProfileTo3d( CVs, srcpl, srccter ) ;

		CVs = convert3dProfileTo2d(CVs3d, dstpl,dstcter ) ;

		// check direction of cut plane normal
		if( srcpl.planeNormal * dstpl.planeNormal < 0 ){
			// inverse Control Points
			std::vector<Point2f> newCVs ;
			std::vector<double> newWs ;
			for( int i = cvnum -1; i>=0; i-- ){
				newCVs.push_back( CVs[i] ) ;
				newWs.push_back( Ws[i] ) ;
			}
			CVs = newCVs ;
			Ws = newWs ;
		}
		setCvOfNurbs( tmp, CVs, Ws ) ;

		return ;
	}


	std::vector<Point3f> convertLocal3dToWorld3d( std::vector<Point3f> points, float *mvmatrix ) {

		std::vector<Point3f> c_w3d ;
		for( int i=0; i<points.size(); ++i )
			c_w3d.push_back( GlobalFun::vecMulMatrix(points[i], mvmatrix) ) ;

		return c_w3d ;
	}
	std::vector<Point3f> convertLocal3dToWorld3d( std::vector<Point3f> points ) {

		float mvmatrix[16] ;
		for(int i=0; i<16; ++i )
			mvmatrix[i] = globalModelviewMatrix[i] ;
		return convertLocal3dToWorld3d(points, mvmatrix) ;
	}


	std::vector<Point3f> convertWorld3dToLocal3d( std::vector<Point3f> points, float *mvmatrix ) {

		float invmvmatrix[16] ;
		GlobalFun::gluInvertMatrix(mvmatrix,invmvmatrix ) ;

		std::vector<Point3f> c_l3d ;
		for( int i=0; i<points.size(); ++i )
			c_l3d.push_back( GlobalFun::vecMulMatrix(points[i], invmvmatrix) ) ;

		return c_l3d ;
	}

	std::vector<Point3f> convertWorld3dToLocal3d( std::vector<Point3f> points ) {

		float mvmatrix[16] ;
		for(int i=0; i<16; ++i )
			mvmatrix[i] = globalModelviewMatrix[i] ;

		float invmvmatrix[16] ;
		GlobalFun::gluInvertMatrix(mvmatrix,invmvmatrix ) ;

		std::vector<Point3f> c_l3d ;
		for( int i=0; i<points.size(); ++i )
			c_l3d.push_back( GlobalFun::vecMulMatrix(points[i], invmvmatrix) ) ;

		return c_l3d ;
	}


	std::vector<Point2f> project3dPointsOntoScreen( std::vector<Point3f> points, float *mvmatrix ){


		double fov = glareaFov * 3.1415926 / 180;
		double nearPlane = -0.1 ;
		double farPlane = -10.0 ;
		double camaraDistForDrawCurve = 0.0 ;

		GLint viewport[4];
		glGetIntegerv (GL_VIEWPORT, viewport);


		// width and height of screen
		double width = viewport[2] ;    
		double height = viewport[3] ;   
		double height_world = tan(fov*0.5) * nearPlane  * 2;
		double width_world = height_world *  width / height ;;

		Point3f origin(0,0,-camaraDistForDrawCurve ) ;	

		// get 2d points on screen
		double zpos_scn = -camaraDistForDrawCurve + nearPlane  ;
		std::vector<Point2f> points2d ;
		for( int i=0; i< points.size(); ++i ){
			Point3f point = points[i] ;
			Point3f point_inworld = GlobalFun::vecMulMatrix(point, mvmatrix) ;
			Point3f vec = point_inworld - origin ;
			Point3f screenDot = ( vec / vec.Z() ) * ( zpos_scn - origin.Z() ) ; 
			Point2f screenDot_noz( screenDot.X(), screenDot.Y() ) ;
			points2d.push_back( screenDot_noz ) ;
		}

		return points2d ;

	}

	std::vector<Point2f> project3dPointsOntoScreen( std::vector<Point3f> points ){

		float mvmatrix[16] ;
		for(int i=0; i<16; ++i )
			mvmatrix[i] = globalModelviewMatrix[i] ;
	
		return project3dPointsOntoScreen( points,mvmatrix) ;

	}

	std::vector<Point2f> convertScreen2dToWorld2d( std::vector<Point2f> curve ){
	
		Curve2D &c2 = curve ;
		int c2size = c2.size() ;

		double fov = glareaFov * 3.1415926 / 180;
		double nearPlane = -0.1 ;
		double farPlane = -10.0 ;
		double camaraDistForDrawCurve = 0.0 ;

		GLint viewport[4];
		glGetIntegerv (GL_VIEWPORT, viewport);


		// width and height of screen
		double width = viewport[2] ;    
		double height = viewport[3] ;   
		double height_world = tan(fov*0.5) * nearPlane  * 2;
		double width_world = height_world *  width / height ;;

		// get 2d curves coordinate in world
		Point3f origin(0,0,-camaraDistForDrawCurve ) ;	
		Curve2D curve2d ;
		for( int i=0; i<c2size; ++i )
			curve2d.push_back( Point2f( (1.0 - c2[i].X()/width) * width_world - width_world/2,  (1.0 - c2[i].Y()/height) * height_world - height_world/2 ) ) ; // on screen

		return curve2d ;
	}

	std::vector<Point2f> convertWorld2dToScreen2d( std::vector<Point2f> curve ){  // mark

		Curve2D &c2 = curve ;
		int c2size = c2.size() ;

		double fov = glareaFov * 3.1415926 / 180;
		double nearPlane = -0.1 ;
		double farPlane = -10.0 ;
		double camaraDistForDrawCurve = 0.0 ;

		GLint viewport[4];
		glGetIntegerv (GL_VIEWPORT, viewport);

		// width and height of screen
		double width = viewport[2] ;    
		double height = viewport[3] ;   
		double height_world = tan(fov*0.5) * nearPlane  * 2;
		double width_world = height_world *  width / height ;;

		Point3f origin(0,0,-camaraDistForDrawCurve ) ;	
		Curve2D curve_s2d ;
		for( int i=0; i<c2size; ++i )
			curve_s2d.push_back( Point2f( width - (c2[i].X() + width_world/2 ) * width/width_world,   height - (c2[i].Y() + height_world/2 ) * height/height_world) ) ;

		return curve_s2d ;
	}


	std::vector<Point3f> convertScreen2dToWorld3d( std::vector<Point2f> curve ){

		std::vector<Point2f> c_w2d = convertScreen2dToWorld2d(curve) ;

		std::vector<Point3f> c_w3d ;
		for( int i=0; i<c_w2d.size(); ++i )
			c_w3d.push_back( Point3f( c_w2d[i][0],c_w2d[i][1],-0.1) ) ;
		return c_w3d ;
	}


	std::vector<Point3f> convertScreen2dToLocal3d( std::vector<Point2f> curve ){

		std::vector<Point3f> c_w3d = convertScreen2dToWorld3d(curve) ;

		float mvmatrix[16] ;
		float invmvmatrix[16] ;
		for(int i=0; i<16; ++i )
			mvmatrix[i] = globalModelviewMatrix[i] ;

		GlobalFun::gluInvertMatrix(mvmatrix,invmvmatrix ) ;

		std::vector<Point3f> c_l3d ;
		for( int i=0; i<c_w3d.size(); ++i )
			c_l3d.push_back( GlobalFun::vecMulMatrix(c_w3d[i], invmvmatrix) ) ;

		return c_l3d ;
	}

	int NearstPoint( std::vector<Point3f> &datapts, Point3f pt){
		//(3D) find the nearest point to the second parameter "pt" among the first parameter datapts
		double mindis = 1e10 ;
		int nstId =0;

		for( int i=0;i<datapts.size(); ++i )
			if( (datapts[i]-pt).Norm() < mindis ){
				mindis = (datapts[i]-pt).Norm() ;
				nstId = i ;
			}
		return nstId ;
	}

	int NearstPoint( std::vector<Point2f> &datapts, Point2f pt){
		//(2D) find the nearest point to the second parameter "pt" among the first parameter datapts
		double mindis = 1e10 ;
		int nstId =0;

		for( int i=0;i<datapts.size(); ++i )
			if( (datapts[i]-pt).Norm() < mindis ){
				mindis = (datapts[i]-pt).Norm() ;
				nstId = i ;
			}
		return nstId ;
	}

	//ON_NurbsCurve buildNurbs( std::vector<Point2f> cvpts, std::vector<double> Weights  ){

	//	unsigned order = ReconstructorPara::NurbsDegree + 1;
	//	unsigned n_control_points ( ReconstructorPara::cvNum + ReconstructorPara::NurbsDegree);
	//	ON_NurbsCurve *pnbs = ON_NurbsCurve::New(3,true,order, n_control_points ) ;
	//	ON_NurbsCurve nbs = *pnbs ;  delete pnbs ;


	//	setCvOfNurbs( nbs, cvpts,Weights ) ;

	//	return nbs ;

	//}

	std::vector<Point2f> getCvPts(ON_NurbsCurve nbs){

		Profile2D  ctrlpoints ; 
		int cvnum = nbs.CVCount() - nbs.Degree() ;
		for( int cpid = 0; cpid<cvnum ; ++cpid ){
			ON_4dPoint p ;
			nbs.GetCV(cpid, p) ;
			Point2f cvp = Point2f(p[0]/p[3], p[1]/p[3]) ;
			ctrlpoints.push_back( cvp  ) ;
		}

		return ctrlpoints ;
	}

	std::vector<Profile2D> cvtProf2Pixels( std::vector<Profile2D> profiles, Box2f box ){

		int n=profiles.size() ;

		double halfwdth = (std::max)(  (std::max)( fabs( box.max.X() ), fabs( box.min.X() ) ) , (std::max)( fabs( box.max.Y() ), fabs( box.min.Y() ) )    ) ;  
		box.min = Point2f(-halfwdth,-halfwdth) ;
		box.max = Point2f(halfwdth,halfwdth) ;

		//box.min -= Point2f(0.02,0.02) ;
		//box.max += Point2f(0.02,0.02) ;

		box.min *= 1.3 ;
		box.max *= 1.3 ;

		float dim = (box.max - box.min).X() > (box.max - box.min).Y() ? (box.max - box.min).X() : (box.max - box.min).Y() ;


		//convert profiles to pixels
		GLint viewport[4];
		glGetIntegerv (GL_VIEWPORT, viewport);
		int height = viewport[3] ;

		std::vector<Profile2D> pixelss(n) ;
		for( int i=0; i<n; ++i  ){

			int start_y ;
			if( n<=3  )  start_y = ( height - ReconstructorPara::profileOnLeftScreenWidth * n ) /2 +  ReconstructorPara::profileOnLeftScreenWidth * i;
			else if( i>2) start_y = ( height - ReconstructorPara::profileOnLeftScreenWidth * (n-3) ) /2 +  ReconstructorPara::profileOnLeftScreenWidth * (i-3);
			else start_y = ( height - ReconstructorPara::profileOnLeftScreenWidth * 3 ) /2 +  ReconstructorPara::profileOnLeftScreenWidth * i;

			int start_x = 50;
			if( i>=3)
				start_x += ReconstructorPara::profileOnLeftScreenWidth ;

			for( int j=0; j<profiles[i].size(); ++j ){
				pixelss[i].push_back (  Point2f(start_x, start_y) + (profiles[i][j] - box.min ) * ReconstructorPara::profileOnLeftScreenWidth / dim );
			}

		}
		return pixelss ;
	}


	std::vector<Profile2D>  cvtPixels2Prof( std::vector<Profile2D> pixels, Box2f box  ){

		int n = pixels.size() ;

		double halfwdth = (std::max)(  (std::max)( fabs( box.max.X() ), fabs( box.min.X() ) ) , (std::max)( fabs( box.max.Y() ), fabs( box.min.Y() ) )    ) ;  
		box.min = Point2f(-halfwdth,-halfwdth) ;
		box.max = Point2f(halfwdth,halfwdth) ;

		//box.min -= Point2f(0.02,0.02) ;
		//box.max += Point2f(0.02,0.02) ;

		box.min *= 1.3 ;
		box.max *= 1.3 ;

		float dim = (box.max - box.min).X() > (box.max - box.min).Y() ? (box.max - box.min).X() : (box.max - box.min).Y() ;


		// convert pixels to control points
		std::vector<Profile2D> profs ;

		profs.resize(n) ;

		GLint viewport[4];
		glGetIntegerv (GL_VIEWPORT, viewport);
		int height = viewport[3] ;


		for( int i=0; i<n; ++i  ){

			//int start_y = ( height - ReconstructorPara::profileOnLeftScreenWidth * n ) /2 +  ReconstructorPara::profileOnLeftScreenWidth * i;
			//int start_y = ( height - ReconstructorPara::profileOnLeftScreenWidth * 3 ) /2 +  ReconstructorPara::profileOnLeftScreenWidth * (i%3);	
			int start_y ;
			if( n<=3  )  start_y = ( height - ReconstructorPara::profileOnLeftScreenWidth * n ) /2 +  ReconstructorPara::profileOnLeftScreenWidth * i;
			else if( i>2) start_y = ( height - ReconstructorPara::profileOnLeftScreenWidth * (n-3) ) /2 +  ReconstructorPara::profileOnLeftScreenWidth * (i-3);
			else start_y = ( height - ReconstructorPara::profileOnLeftScreenWidth * 3 ) /2 +  ReconstructorPara::profileOnLeftScreenWidth * i;


			int start_x = 50;
			if( i>=3)
				start_x += ReconstructorPara::profileOnLeftScreenWidth ;


			for( int j=0; j<pixels[i].size(); ++j ){
				profs[i].push_back(  ( pixels[i][j] - Point2f(start_x, start_y) ) *  (dim / ReconstructorPara::profileOnLeftScreenWidth) +  box.min  );
			}
		}

		return profs ;
	}



	bool isTip( skelpath skel, Point3f p ) {

		int count = 0 ;
		int2 NodeId ;
		std::vector<std::vector<Point3f>> &impSkel = skel.smoothedBrchPts ;
		for( int i=0; i<impSkel.size(); ++i ){
			if( (impSkel[i][0] - p).Norm() < 0.05 ){
				count++ ;
				NodeId = int2(i,0) ;
			}
			if( (impSkel[i].back() - p).Norm() < 0.05 ){
				count++ ;
				NodeId = int2(i, impSkel[i].size()-1 ) ;
			}
		}
		if( count != 1 )
			return false ;
		int Bid = NodeId.x ;
		int Nid = NodeId.y ;
		std::vector<Profile2D> &prof2d = skel.profiles2d[Bid] ;
		while( prof2d[Nid].size() < 10 ){
			if( Nid < prof2d.size()/2 ) Nid++ ; 
			else Nid-- ;
			if( Nid == prof2d.size()/2 )
				return true ;

		}

		Box2f box = GlobalFun::computeBox( prof2d[Nid] ) ;

		if( (box.max - box.min).Norm() < 0.1 )
			return true ;
		else
			return false ;


	}





	bool isMidNode( std::vector<Curve3D>  impSkel, Point3f p ) {

		int count = 0 ;
		int2 NodeId ;
		for( int i=0; i<impSkel.size(); ++i ){
			if( (impSkel[i][0] - p).Norm() < 0.05 ){
				count++ ;
				NodeId = int2(i,0) ;
			}
			if( (impSkel[i].back() - p).Norm() < 0.05 ){
				count++ ;
				NodeId = int2(i, impSkel[i].size()-1 ) ;
			}
		}

		if( count > 1 )
			return true ;
		else
			return false ;


	}




	bool CalPlaneLineIntersectPoint(Point3f &interPoint, Point3f planeVector, Point3f planePoint, Point3f lineVector, Point3f linePoint  ){

		//calculate if the plane and the line is intersect ,and return the bool result
		//if do intersect , calculate the intersect point , and save it in the first parameter
		float vp1, vp2, vp3, n1, n2, n3, v1, v2, v3, m1, m2, m3, t,vpt;  
		vp1 = planeVector[0];  
		vp2 = planeVector[1];  
		vp3 = planeVector[2];  
		n1 = planePoint[0];  
		n2 = planePoint[1];  
		n3 = planePoint[2];  
		v1 = lineVector[0];  
		v2 = lineVector[1];  
		v3 = lineVector[2];  
		m1 = linePoint[0];  
		m2 = linePoint[1];  
		m3 = linePoint[2];  
		vpt = v1 * vp1 + v2 * vp2 + v3 * vp3;  


		if (vpt == 0) {  
			// the plane and the line is parallel or overlapping
			return false ;
		}  
		else  
		{  
			// the plane and the line is intersect
			t = ((n1 - m1) * vp1 + (n2 - m2) * vp2 + (n3 - m3) * vp3) / vpt;  
			interPoint[0] = m1 + v1 * t;  
			interPoint[1] = m2 + v2 * t;  
			interPoint[2] = m3 + v3 * t;  
		}  
		return true;  
	}

	void smoothVectorData( std::vector<double> &v, int winsize ){


		double sigma = 0.3 * ( winsize / 2.0 - 1.0 )  + 0.8  ;
		 std::vector<double> u = v ;
		for( int i=0; i<v.size(); ++i ){
			
			double totalWeight = 0 ;
			double totalValues = 0 ;
			for( int j=0; j<u.size(); ++j ){
				int offset = abs(i-j) ;
				double weight = pow(2.71828, -offset*offset/(sigma*sigma) ) ; 
				totalWeight += weight ;
				totalValues += u[j] * weight ;
			}
			v[i] = totalValues / totalWeight ;
		}
	}


	void CopyCmesh( CMesh & dst, CMesh &src ){

		triangleList  tl(src) ;

		int nvertices = tl.nvertices ;
		int nfaces = tl.nfaces ;
		float *vertices = tl.vertices ;
		float *normals = tl.normals ;
		unsigned *indices = tl.indices ;
		
		CMesh &m = dst ;

		m.Clear() ;

		for( int i=0; i<nvertices; ++i ){
			CMesh::VertexIterator vi = vcg::tri::Allocator<CMesh>::AddVertices(m,1);
			vi->P()=CMesh::CoordType ( vertices[i* 3 + 0] , vertices[ i * 3 + 1] , vertices[ i * 3 + 2] ); 
			vi->N()=CMesh::CoordType ( normals [i* 3 + 0] , normals[ i * 3 + 1] , normals[ i * 3 + 2] ); 
		}

		for( int i=0; i<nfaces; ++i ){

			CMesh::FaceIterator fi = vcg::tri::Allocator<CMesh>::AddFaces(m,1);
			CMesh::VertexPointer ivp[3];
			ivp[2]= &( m.vert[indices[i*3+0]] ) ;
			ivp[1]= &( m.vert[indices[i*3+1]] ) ;
			ivp[0]= &( m.vert[indices[i*3+2]] ) ;
			fi->V(0)=ivp[0];
			fi->V(1)=ivp[1];
			fi->V(2)=ivp[2];
		}

		m.vn = m.vert.size() ;

	}



	void drawPointsAsSphere( std::vector<Point3f> pts, GLColor color, int size ){


		glEnable(GL_LIGHTING) ;
		glEnable( GL_LIGHT0 ) ;
		glEnable( GL_LIGHT1 ) ;

		double radius = 0.002 * size ;

		GLUquadricObj *objCylinder = gluNewQuadric();

		glColor4f( color.r, color.g, color.b, color.a ) ;

		for( int i=0; i<pts.size(); ++i){

			glPushMatrix() ;
			glTranslatef(pts[i].X(), pts[i].Y(), pts[i].Z() ) ;
			gluSphere(objCylinder, radius, 20, 20 );
			glPopMatrix() ;
		}
		
		gluDeleteQuadric( objCylinder ) ;

		glDisable(GL_LIGHTING) ;
	}

	void draw3dCurves( 	Curve3D &curves3d, GLColor color, bool closed, bool  randColor, double size ){

		glColor4f( color.r,color.g,color.b, color.a) ;
		glLineWidth(size) ;

		if( closed )
			glBegin(GL_LINE_LOOP) ;
		else
			glBegin(GL_LINE_STRIP) ;

		for( int j=0; j<curves3d.size(); ++j  )
			glVertex3f( curves3d[j].X(),  curves3d[j].Y(),  curves3d[j].Z() ) ;
		glEnd() ;


	}


	bool loadSkeletonAndProfiles( skelpath &skel, int ogSkelId ){

		char buffer[1024] ;

		std::string fname1 = appstate.path + "\/ogskel_" + itoa( ogSkelId, buffer,10 ) + ".txt" ;
		
		std::ifstream ifs(fname1) ;

		if( ifs.fail() ){
			std::cout << "didn't find " << fname1 <<std::endl;
			return false ;
		}


		std::vector<int> tmpIdx ;
		std::vector<ON_NurbsCurve> nbss;
		Curve3D  brch ;
		std::vector<double>  conf ;

		// load skeleton
		int brchSize ;
		ifs >> brchSize ;
		for( int i=0; i<brchSize; ++i ){
			double x,y,z ;  
			ifs >> x >> y >> z;  
			brch.push_back( Point3f(x,y,z) ) ;
		}


		// load confidence
		int confSize ;
		ifs >> confSize ;
		for( int i=0; i<confSize; ++i ){
			double t;  ifs >>t ;   conf.push_back( t) ;
		}


		// load tmpId 
		int tmpNum ;
		ifs >> tmpNum ;
		for( int i=0; i<tmpNum; ++i ){
			int id;
			ifs >> id ;
			tmpIdx.push_back( id ) ;
		}


		// load nurbs

		int nbsNum ;
		ifs >> nbsNum ;

		for( int i=0; i<nbsNum; ++i ){
			Curve2D cvs  ;
				
			for( int j=0; j< ReconstructorPara::cvNum; ++j ){
				Point2f p ;
				ifs >> p.X() >> p.Y() ;
				cvs.push_back(p) ;
			}

			ON_NurbsCurve nbs = ReconstructorUtility::fitNURBSforProfile2d(  std::vector<Point2f>(10, Point2f(0,0) ) ) ;
			ReconstructorUtility::setCvOfNurbs(nbs, cvs) ;
			nbss.push_back( nbs ) ;
		}

		ifs.close() ;


		if( tmpIdx.size() != nbss.size() ){
			std::cout <<__FILE__<<__LINE__<<std::endl;
			system("error") ;
		}

		// rank tmpId and nbs
		std::vector<int> newTmpId = tmpIdx ;
		std::sort( newTmpId.begin(), newTmpId.end() ) ;

		std::vector<ON_NurbsCurve> newNbss;
		for( int i=0; i<newTmpId.size(); ++i )
			newNbss.push_back(  nbss[ReconstructorUtility::FindInVerctor(tmpIdx, newTmpId[i] ) ] ) ;
		
		tmpIdx = newTmpId ;
		nbss = newNbss ;




		skel.smoothedBrchPts = std::vector<Curve3D>(1, brch) ;

		skel.calculateCutplanes() ;
		skel.tempalteIdss.resize(1) ;
		skel.Nurbsss.resize(1) ;
		skel.tempalteIdss[0] = tmpIdx;
		skel.Nurbsss[0] = nbss;
		skel.discretizeNurbs() ;




		
		// detect loop
		int n=brch.size() ;
		if( (brch[0] - brch.back()).Norm() < 0.05  && (brch[0] - brch[1] ) *(brch[n-1] - brch[n-2] ) < 0  ){
			skel.brch0IsLoop = true ;
		}else
			skel.brch0IsLoop = false ;



	}
	void updateSkelDepth( skeleton_mul &src, skeleton_mul dst ){ // not rubost engough

		CurveArray3D brches = src.smoothedBrchPts ;

		CurveArray2D brches_w2d ;  
		CurveArray3D brches_w3d ;
		for( int i=0; i<brches.size(); ++i ){
			brches_w2d.push_back(  project3dPointsOntoScreen( brches[i] ) ) ;
			brches_w3d.push_back( convertLocal3dToWorld3d( brches[i] )) ;
		}

		Curve2D refpts_w2d ;
		Curve3D refpts_w3d ;

		for( int i=0; i<dst.smoothedBrchPts.size(); ++i ){
			if( refpts_w2d.size() == 0 ) {
				refpts_w2d = project3dPointsOntoScreen( dst.smoothedBrchPts[i] ) ;
				refpts_w3d = convertLocal3dToWorld3d( dst.smoothedBrchPts[i] ) ;
			}else{
				Curve2D brchiw2d = project3dPointsOntoScreen( dst.smoothedBrchPts[i] ) ;
				refpts_w2d.insert( refpts_w2d.begin(), brchiw2d.begin(), brchiw2d.end()  ) ;

				Curve3D brchiw3d = convertLocal3dToWorld3d( dst.smoothedBrchPts[i] ) ;
				refpts_w3d.insert( refpts_w3d.begin(), brchiw3d.begin(), brchiw3d.end() ) ;
			}
		}


		//Curve2D refJoints_w2d ;
		//Curve3D refJoints_w3d ;
		//Curve3D joints ;  std::vector<std::vector<int>> abjbid ;
		//dst.getJoints( joints, abjbid ) ;
		//refJoints_w2d = project3dPointsOntoScreen(joints ) ;
		//refJoints_w3d = convertLocal3dToWorld3d(joints ) ;


		for( int bid = 0; bid<brches_w2d.size(); ++bid ){

			std::vector<std::vector<int>> nstId ,nstIdJoints;

			ReconstructorUtility::computeFlannKNN( refpts_w2d, brches_w2d[bid],nstId,1 ) ;
			//ReconstructorUtility::computeFlannKNN( refJoints_w2d, brches_w2d[bid],nstIdJoints,1 ) ;

			for( int nid=0; nid<brches_w2d[ bid ].size() ; ++nid ){

				//if( nid!=0 && nid != brches_w2d[ bid ].size()-1 ){
				double z = refpts_w3d[ nstId[nid][0] ].Z() ;
				brches_w3d[bid][nid] = brches_w3d[bid][nid] * ( z/brches_w3d[bid][nid].Z() ) ;
				//}else{
				//	//double z = refJoints_w3d[ nstIdJoints[nid][0] ].Z() ;
				//	//brches_w3d[bid][nid] = brches_w3d[bid][nid] * ( z/brches_w3d[bid][nid].Z() ) ;
				//}
			}


		}


		Curve3D joints ;  std::vector<std::vector<int>> abjbid ;

		src.getJoints( joints, abjbid ) ;

		for( int i=0; i<abjbid.size(); ++i ){

			Point3f sum(0,0,0 ) ;
			for( int j=0; j<abjbid[i].size(); ++j ){
				if( (src.smoothedBrchPts[abjbid[i][j]][0] - joints[i] ).Norm() < (src.smoothedBrchPts[abjbid[i][j]].back() - joints[i] ).Norm() )
					sum = sum + brches_w3d[abjbid[i][j]][0] ;
				else
					sum = sum + brches_w3d[abjbid[i][j]].back() ;
			}

			Point3f jpos = sum / abjbid[i].size() ;
			for( int j=0; j<abjbid[i].size(); ++j ){
				if( (src.smoothedBrchPts[abjbid[i][j]][0] - joints[i] ).Norm() < (src.smoothedBrchPts[abjbid[i][j]].back() - joints[i] ).Norm() )
					brches_w3d[abjbid[i][j]][0] = jpos;
				else
					brches_w3d[abjbid[i][j]].back() = jpos;
			}

		}


		for( int i=0; i<brches.size(); ++i ){
			brches[i] = convertWorld3dToLocal3d( brches_w3d[i] ) ;
			getBezier( brches[i], brches[i], 100 ) ;
		}


		for( int i=0; i<abjbid.size(); ++i ){

			Point3f sum(0,0,0 ) ;
			for( int j=0; j<abjbid[i].size(); ++j ){
				if( (src.smoothedBrchPts[abjbid[i][j]][0] - joints[i] ).Norm() < (src.smoothedBrchPts[abjbid[i][j]].back() - joints[i] ).Norm() )
					sum = sum + brches[abjbid[i][j]][0] ;
				else
					sum = sum + brches[abjbid[i][j]].back() ;
			}

			Point3f jpos = sum / abjbid[i].size() ;
			for( int j=0; j<abjbid[i].size(); ++j ){
				if( (src.smoothedBrchPts[abjbid[i][j]][0] - joints[i] ).Norm() < (src.smoothedBrchPts[abjbid[i][j]].back() - joints[i] ).Norm() )
					brches[abjbid[i][j]][0] = jpos;
				else
					brches[abjbid[i][j]].back() = jpos;
			}

		}

		src = skeleton_mul( brches , false) ;
	}

	void updateSkelDepth( skelpath &src, skelpath dst ){ // not rubost engough

		CurveArray3D brches = src.smoothedBrchPts ;

		CurveArray2D brches_w2d ;  
		CurveArray3D brches_w3d ;
		for( int i=0; i<brches.size(); ++i ){
			brches_w2d.push_back(  project3dPointsOntoScreen( brches[i] ) ) ;
			brches_w3d.push_back( convertLocal3dToWorld3d( brches[i] )) ;
		}

		Curve2D refpts_w2d ;
		Curve3D refpts_w3d ;
		
		for( int i=0; i<dst.smoothedBrchPts.size(); ++i ){
			if( refpts_w2d.size() == 0 ) {
				refpts_w2d = project3dPointsOntoScreen( dst.smoothedBrchPts[i] ) ;
				refpts_w3d = convertLocal3dToWorld3d( dst.smoothedBrchPts[i] ) ;
			}else{
				Curve2D brchiw2d = project3dPointsOntoScreen( dst.smoothedBrchPts[i] ) ;
				refpts_w2d.insert( refpts_w2d.begin(), brchiw2d.begin(), brchiw2d.end()  ) ;

				Curve3D brchiw3d = convertLocal3dToWorld3d( dst.smoothedBrchPts[i] ) ;
				refpts_w3d.insert( refpts_w3d.begin(), brchiw3d.begin(), brchiw3d.end() ) ;
			}
		}

		
		//Curve2D refJoints_w2d ;
		//Curve3D refJoints_w3d ;
		//Curve3D joints ;  std::vector<std::vector<int>> abjbid ;
		//dst.getJoints( joints, abjbid ) ;
		//refJoints_w2d = project3dPointsOntoScreen(joints ) ;
		//refJoints_w3d = convertLocal3dToWorld3d(joints ) ;


		for( int bid = 0; bid<brches_w2d.size(); ++bid ){

			std::vector<std::vector<int>> nstId ,nstIdJoints;
			
			ReconstructorUtility::computeFlannKNN( refpts_w2d, brches_w2d[bid],nstId,1 ) ;
			//ReconstructorUtility::computeFlannKNN( refJoints_w2d, brches_w2d[bid],nstIdJoints,1 ) ;

			for( int nid=0; nid<brches_w2d[ bid ].size() ; ++nid ){
				
				//if( nid!=0 && nid != brches_w2d[ bid ].size()-1 ){
					double z = refpts_w3d[ nstId[nid][0] ].Z() ;
					brches_w3d[bid][nid] = brches_w3d[bid][nid] * ( z/brches_w3d[bid][nid].Z() ) ;
				//}else{
				//	//double z = refJoints_w3d[ nstIdJoints[nid][0] ].Z() ;
				//	//brches_w3d[bid][nid] = brches_w3d[bid][nid] * ( z/brches_w3d[bid][nid].Z() ) ;
				//}
			}
			

		}

		
		Curve3D joints ;  std::vector<std::vector<int>> abjbid ;
		
		src.getJoints( joints, abjbid ) ;

		for( int i=0; i<abjbid.size(); ++i ){

			Point3f sum(0,0,0 ) ;
			for( int j=0; j<abjbid[i].size(); ++j ){
				if( (src.smoothedBrchPts[abjbid[i][j]][0] - joints[i] ).Norm() < (src.smoothedBrchPts[abjbid[i][j]].back() - joints[i] ).Norm() )
					sum = sum + brches_w3d[abjbid[i][j]][0] ;
				else
					sum = sum + brches_w3d[abjbid[i][j]].back() ;
			}

			Point3f jpos = sum / abjbid[i].size() ;
			for( int j=0; j<abjbid[i].size(); ++j ){
				if( (src.smoothedBrchPts[abjbid[i][j]][0] - joints[i] ).Norm() < (src.smoothedBrchPts[abjbid[i][j]].back() - joints[i] ).Norm() )
					brches_w3d[abjbid[i][j]][0] = jpos;
				else
					brches_w3d[abjbid[i][j]].back() = jpos;
			}

		}


		for( int i=0; i<brches.size(); ++i ){
			brches[i] = convertWorld3dToLocal3d( brches_w3d[i] ) ;
			getBezier( brches[i], brches[i], 100 ) ;
		}


		for( int i=0; i<abjbid.size(); ++i ){

			Point3f sum(0,0,0 ) ;
			for( int j=0; j<abjbid[i].size(); ++j ){
				if( (src.smoothedBrchPts[abjbid[i][j]][0] - joints[i] ).Norm() < (src.smoothedBrchPts[abjbid[i][j]].back() - joints[i] ).Norm() )
					sum = sum + brches[abjbid[i][j]][0] ;
				else
					sum = sum + brches[abjbid[i][j]].back() ;
			}

			Point3f jpos = sum / abjbid[i].size() ;
			for( int j=0; j<abjbid[i].size(); ++j ){
				if( (src.smoothedBrchPts[abjbid[i][j]][0] - joints[i] ).Norm() < (src.smoothedBrchPts[abjbid[i][j]].back() - joints[i] ).Norm() )
					brches[abjbid[i][j]][0] = jpos;
				else
					brches[abjbid[i][j]].back() = jpos;
			}

		}

		src = skelpath( brches , false) ;
	}

	//added by huajie,2014/8/21
	using namespace Eigen;
	//************************************
	// 函数名称: ComputeEigenvector
	// 函数说明：ComputeEigenvector
	// 作    者：huajie.Shi
	// 日    期：2014.8.21
	// 返 回 值: Point3f
	// 参    数: std::vector< Point3f > nmls
	//************************************
	Point3f ComputeEigenvector( std::vector< Point3f > nmls){
		int size = nmls.size();
		double x, y, z, x2, y2, z2, xy, xz, yz;
		
		x = y = z = x2 = y2 = z2 = xy = xz = yz = 0;
		for(int i=0; i<size; i++){
			x += nmls[i].X();
			y += nmls[i].Y();
			z += nmls[i].Z();
			x2 += nmls[i].X() * nmls[i].X();
			y2 += nmls[i].Y() * nmls[i].Y();
			z2 += nmls[i].Z() * nmls[i].Z();
			xy += nmls[i].X() * nmls[i].Y();
			xz += nmls[i].X() * nmls[i].Z();
			yz += nmls[i].Y() * nmls[i].Z();
		}
		x /= size; x2 /= size;
		y /= size; y2 /= size;
		z /= size; z2 /= size;
		xy /= size; xz /= size; yz /= size;
		MatrixXf M(3,3); 
		M(0,0) = x2 - x*x; M(0,1) = xy - x*y; M(0,2) = xz - x*z;
		M(1,0) = xy - x*y; M(1,1) = y2 - y*y; M(1,2) = yz - y*z;
		M(2,0) = xz - x*z; M(2,1) = yz - y*z; M(2,2) = z2 - z*z;
		EigenSolver<MatrixXf> ces(M);
		//cout << "The first eigenvector of the 3x3 matrix of ones is:"
		//	<< endl << ces.eigenvectors().col(0) << endl;
		//cout << "The first eigenvalue of the 3x3 matrix of ones is:"
		//	<< endl << ces.eigenvalues() << endl;
		Point3f normal;
		double a[3];
		int maxCol=0, minCol=0;
		a[0] = ces.eigenvalues()(0,0).real();
		a[1] = ces.eigenvalues()(1,0).real();
		a[2] = ces.eigenvalues()(2,0).real();
		for(int i=0; i<3; i++){
			if( fabs(a[i] - max( a[0], max( a[1], a[2] ) )) < 0.00001 )
				maxCol = i;
			if( fabs(a[i] - min( a[0], min( a[1], a[2] ) )) < 0.00001 )
				minCol = i;
		}
		normal.X() = ces.eigenvectors()(0,minCol).real();
		normal.Y() = ces.eigenvectors()(1,minCol).real();
		normal.Z() = ces.eigenvectors()(2,minCol).real();
		return normal;
	}
}