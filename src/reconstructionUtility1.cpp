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
#include <cmath>
namespace ReconstructorUtility{

	template <class _Point>

	//return the degree of the angle between two vector
	double vectorIntersectionAngle( _Point p1, _Point p2 ){
		double degree = acos( (std::min)( (std::max)( -1.0, (double)( p1 * p2) /( p1.Norm()*p2.Norm())  ) , 1.0 )  ) ;
		return degree ;
	}

	std::vector<double> convertTransform2Vector ( const std::vector<ST> &stf ) {
		std::vector<double> x0 ;
		for( int i=0; i<stf.size(); ++i ){
			x0.push_back( stf[i].rotateAngle ) ;
			x0.push_back( stf[i].scaling ) ;
			x0.push_back( stf[i].translate[0] ); 
			x0.push_back( stf[i].translate[1] ) ;
		}
		return x0 ;
	}

	std::vector<ST> convertVector2Transform (const std::vector<double> &x0 ) {
		std::vector<ST> stf( x0.size()/4 ) ;
		for( int i=0; i<x0.size()/4; ++i ){
			stf[i].rotateAngle  = x0[i*4] ;
			stf[i].scaling      = x0[i*4+1] ;
			stf[i].translate[0] = x0[i*4+2] ;
			stf[i].translate[1] = x0[i*4+3] ;
		}
		return stf ;
	}

	std::vector<double> convertNRTransform2Vector ( const std::vector<NRTSF> &nrtsf ) {
		std::vector<double> x0 ;
		for( int i=0; i<nrtsf.size(); ++i ){

			for( int j=0; j<nrtsf[i].t.size(); ++j ){
				x0.push_back( nrtsf[i].t[j].X() ) ;
				x0.push_back( nrtsf[i].t[j].Y() ) ;
			}

			x0.push_back( nrtsf[i].A ) ;

			x0.push_back( nrtsf[i].b.X() ) ;
			x0.push_back( nrtsf[i].b.Y() ) ;
		}

		return x0 ;
	}

	std::vector<NRTSF> convertVector2NRTransform (const std::vector<double> &x0, int cvNum ) {
		
		std::vector<NRTSF> res( x0.size()/(cvNum*2+3) );

		for( int i=0; i<res.size(); ++i ){
			
			res[i].t.resize(cvNum) ;

			int offset = i * (cvNum*2+3) ;

			for( int j=0; j<cvNum*2; ++j ){
				if( j%2 == 0 )
					res[i].t[j/2].X() = x0[offset+j] ; 
				else
					res[i].t[j/2].Y() = x0[offset+j] ; 
			}

			res[i].A = x0[offset+cvNum*2] ; 

			res[i].b.X() = x0[offset+cvNum*2+1] ; 
			res[i].b.Y() = x0[offset+cvNum*2+2] ; 
		}
		
		return res ;
	}

	Profile2D transformProfile( const Profile2D &input ,  ST stf) {

		//transform the profile(the first parameter input) with ST(the second parameter stf) ,
		//and return the result
		Profile2D res( input.size() ) ;

		double theta = stf.rotateAngle ;
		for( int i=0; i< input.size(); ++i ){

			double x = input[i][0] ;
			double y = input[i][1] ;

			res[i] = Point2f( x*cos(theta) - y*sin(theta) ,  x*sin(theta) + y*cos(theta) )  * stf.scaling + stf.translate ;
		}
		return res ;
	}

	Profile3D transformProfile3D( Profile3D input ,  ST stf, cutPlane cutpl_src, Point3f cter_src, cutPlane cutpl_dst, Point3f cter_dst ) {

		Profile3D tmp  ;
		// ------------------------ locally transformation -----------------------
		// calculating translation
		Point3f translation ; 
		//double degree = acos( std::min( std::max( -1.0, (double)(  Point2f(0,1) * stf.translate ) / stf.translate.Norm() ) , 1.0 )  ) ;
		double degree = vectorIntersectionAngle( Point2f(0,1) , stf.translate ) ;


		std::cout << "trans degree = " << degree <<std::endl;
		Point3f axis_ =  Point3f(0,1,0) ^ Point3f( stf.translate[0],stf.translate[1], 0 ) ;
		if( axis_[2] > 0 ){

			std::cout<< "axis_[2] > 0" <<std::endl;
			GlobalFun::Rotate_Point3D(degree, cutpl_src.planeNormal, cutpl_src.cutRadius, translation ) ;  
		}
		else{
			std::cout<< "axis_[2] < 0" <<std::endl;
			GlobalFun::Rotate_Point3D( -degree, cutpl_src.planeNormal, cutpl_src.cutRadius, translation ) ;  
		}
		translation = translation.Normalize() * stf.translate.Norm() ;   // mark
		//translation.Y() *= -1 ;

		std::cout << "2d translat: "<<stf.translate.X() <<","<< stf.translate.Y() <<std::endl;
		std::cout << "translation: "<<translation.X() <<","<<translation.Y() <<","<<translation.Z() <<"\n"<<std::endl;

		GlobalFun::translatePoints(input, -cter_src) ; 
		// rotate
		GlobalFun::Rotate_Point3Ds(stf.rotateAngle, cutpl_src.planeNormal, input, tmp ) ; input = tmp ;
		// scaling
		for( int i=0; i<input.size(); ++i ){
			Point3f norDirComp = cutpl_src.planeNormal  * (input[i] * cutpl_src.planeNormal ) ;
			input[i] = (input[i] - norDirComp) * stf.scaling + norDirComp;
		}
		//translate 
		GlobalFun::translatePoints(input, cter_src + translation ) ;

		// -----------------------------------------------------------------


		// transform it to dst plane
		GlobalFun::translatePoints(input, cter_dst-cter_src ) ;

		GlobalFun::translatePoints(input, -cter_dst) ; 

		//  match plane normal
		degree = vectorIntersectionAngle( cutpl_src.planeNormal, cutpl_dst.planeNormal) ;
		Point3f raxis = cutpl_src.planeNormal ^ cutpl_dst.planeNormal ;
		GlobalFun::Rotate_Point3Ds(degree, raxis, input, tmp ) ;  input = tmp ;

		// match raidus
		Point3f srcRadius ;
		GlobalFun::Rotate_Point3D( degree, raxis, cutpl_src.cutRadius, srcRadius ) ;
		degree = vectorIntersectionAngle( srcRadius, cutpl_dst.cutRadius ) ;
		GlobalFun::Rotate_Point3Ds(degree, cutpl_dst.planeNormal, input, tmp ) ;  input = tmp ;


		GlobalFun::translatePoints(input, cter_dst) ; 


		return input ;
	}
	
	std::vector<std::vector<Point2f>> transformProfiles( const std::vector<Point2f> &input ,  std::vector<ST> &stfs){
		std::vector<std::vector<Point2f>> ress( stfs.size() ) ;

		for( int i=0;i<stfs.size(); ++i)
			ress[i] = transformProfile( input, stfs[i] ) ;

		return ress ;
	}


	Profile2D mergeProfile2d( const Profile2D &temProf, const Profile2D &dataProf ){
		Profile2D merged ;
		merged = temProf ;
		merged.insert(merged.end(), dataProf.begin(), dataProf.end() ) ;

		return merged ;
	}

	Profile3D convert2dProfileTo3d( const Profile2D &prof2d, cutPlane cutpl, Point3f center ) {
		
		
		// first compute inter-media cutTadius
		std::vector<Point3f>   radiusVector ;
		radiusVector.push_back(  cutpl.cutRadius ) ;

		std::vector<Point3f>  tmp3d ;
		Point3f roAxis0 = (cutpl.planeNormal ^ Point3f( 0,0,1)).Normalize() ;
		double deg0 = vectorIntersectionAngle(cutpl.planeNormal,  Point3f( 0,0,1) ) ;

		if( fabs( deg0 ) > 3.13 ){
			roAxis0 = Point3f(0,1,0) ;
			GlobalFun::Rotate_Point3Ds( deg0, roAxis0, radiusVector, tmp3d ) ;  radiusVector = tmp3d ;
		}else if(  fabs( deg0 ) > 0.01 ){
			GlobalFun::Rotate_Point3Ds( deg0, roAxis0, radiusVector, tmp3d ) ;  radiusVector = tmp3d ;
		}

		//--------------- align y axis
		
		std::vector<Point3f>  dstpts3d ;
		for( int pid=0;pid<prof2d.size() ; ++pid)
			dstpts3d.push_back( Point3f(prof2d[pid][0],prof2d[pid][1], 0 )   ) ;

		std::vector<Point3f>  tmp ;
		Point3f roAxis = (radiusVector[0]^ Point3f( 0,1,0) ).Normalize();
		double deg = vectorIntersectionAngle(radiusVector[0],  Point3f( 0,1,0) ) ;

		 if( fabs( deg ) > 3.13 ){
			 roAxis = Point3f(0,0,1) ;
			 GlobalFun::Rotate_Point3Ds( -deg, roAxis, dstpts3d, tmp ) ;  dstpts3d = tmp ;
		 }else if(  fabs( deg ) > 0.01 ){
			 GlobalFun::Rotate_Point3Ds( -deg, roAxis, dstpts3d, tmp ) ;  dstpts3d = tmp ;
		 }

		//--------------- align plane normal axis
		//GlobalFun::Rotate_Point3Ds( -deg0, roAxis0, dstpts3d, tmp ) ;  dstpts3d = tmp ;
		 if( fabs( deg0 ) > 3.13 ){
			 roAxis0 = Point3f(0,1,0) ;
			 GlobalFun::Rotate_Point3Ds( -deg0, roAxis0, dstpts3d, tmp ) ;  dstpts3d = tmp ;
		}else if(  fabs( deg0 ) > 0.01 ){
			GlobalFun::Rotate_Point3Ds( -deg0, roAxis0, dstpts3d, tmp ) ;  dstpts3d = tmp ;
		}


		GlobalFun::translatePoints(dstpts3d, center ) ;

		if( dstpts3d.size() > 5 && 0){

			Point3f center = GlobalFun::centerOfPoints(dstpts3d) ;
			bool onePoint = true ;
			for( int i=0; i<dstpts3d.size(); ++i )
				if( (dstpts3d[i]-center).Norm() > 0.01  )
					onePoint = false ;
			if( !onePoint ){
				Point3f v1,v2 ,ct;
				GlobalFun::computePCAplane( dstpts3d, v1, v2, ct ) ;
				double angle = MyMin( ReconstructorUtility::vectorIntersectionAngle(( v1^v2), cutpl.planeNormal ),ReconstructorUtility::vectorIntersectionAngle(( v1^v2), -cutpl.planeNormal ))  ;
				if(  angle > 3.14/20 ){
					std::cout << "convert2dProfileTo3d Error! "<<__FILE__<<__LINE__<<std::endl;
					std::cout << "pca angle = " <<angle <<std::endl;
					std::cout<<"deg0 = "<<deg0<<std::endl;
					std::cout<<"roAxis0 = "<<roAxis0.X()<<", "<<roAxis0.Y()<<", "<<roAxis0.Z()<<std::endl;
					std::cout<<"deg = "<<deg<<std::endl;
					std::cout<<"roAxis = "<<roAxis.X()<<", "<<roAxis.Y()<<", "<<roAxis.Z()<<std::endl;
					std::cout<<"normal = " << cutpl.planeNormal.X() << ", "<<cutpl.planeNormal.Y() << ", "<<cutpl.planeNormal.Z() << std::endl;

					system("pause") ;
				}
			}
		}

		return dstpts3d ;

	}
	Profile2D convert3dProfileTo2d( const Profile3D &prof3d, cutPlane cutpl, Point3f center ) {

		std::vector<Point3f>  dstpts =  prof3d ;

		std::vector<Point3f>  tmp ;

		// add an extra points, for the sack of drawing radius
		dstpts.push_back(center +  cutpl.cutRadius ) ;

		// transform all points to make the cutRadius to be axis y
		GlobalFun::translatePoints(dstpts, -center ) ;

		Point3f roAxis = (cutpl.planeNormal ^ Point3f( 0,0,1)).Normalize() ;
		double deg = vectorIntersectionAngle(cutpl.planeNormal,  Point3f( 0,0,1) ) ;
		if( fabs( deg ) > 3.13 ){
			roAxis = Point3f(0,1,0) ;
			GlobalFun::Rotate_Point3Ds( deg, roAxis, dstpts, tmp ) ;  dstpts = tmp ;
		}else if(  fabs( deg ) > 0.01 ){
			GlobalFun::Rotate_Point3Ds( deg, roAxis, dstpts, tmp ) ;  dstpts = tmp ;
		}


		roAxis = dstpts.back()^ Point3f( 0,1,0) ;
		deg = vectorIntersectionAngle(dstpts.back(),  Point3f( 0,1,0) ) ;
		if( fabs( deg ) > 3.13 ){
			roAxis = Point3f(0,0,1) ;
			GlobalFun::Rotate_Point3Ds( deg, roAxis, dstpts, tmp ) ;  dstpts = tmp ;
		}else if(  fabs( deg ) > 0.05 ){
			GlobalFun::Rotate_Point3Ds( deg, roAxis, dstpts, tmp ) ;  dstpts = tmp ;
		}



		std::vector<Point2f>  dstpts2d ;
		for( int pid=0;pid<prof3d.size() ; ++pid)
			dstpts2d.push_back( Point2f(dstpts[pid][0],dstpts[pid][1] )   ) ;
			
		return dstpts2d ;
	}


	ON_NurbsCurve fitNURBSforProfile2d( Profile2D &prof ) {
		
		// convert to NURBS data structure
		pcl::on_nurbs::NurbsDataCurve2d data;
		for( int i=0; i<prof.size(); ++i )
			data.interior.push_back( Eigen::Vector2d (prof[i].X(), prof[i].Y()) );


		// CURVE PARAMETERS
		unsigned order = ReconstructorPara::NurbsDegree+1 ;  
		unsigned n_control_points ( ReconstructorPara::cvNum + ReconstructorPara::NurbsDegree);  

		// fit a 2 degree nurbs
		order = 3 ;
		n_control_points = ReconstructorPara::cvNum + 2 ;


		pcl::on_nurbs::FittingCurve2dSDM::Parameter curve_params;
		curve_params.smoothness = 0.1;
		curve_params.rScale = 1;

		//  CURVE FITTING
		ON_NurbsCurve curve = pcl::on_nurbs::FittingCurve2dPDM::initNurbsCurve2D (order, data.interior, n_control_points);

		pcl::on_nurbs::FittingCurve2dSDM fit_pdm (&data, curve);
		fit_pdm.assemble (curve_params);
		fit_pdm.solve ();

		int ncv = ReconstructorPara::cvNum + ReconstructorPara::NurbsDegree ;

		ON_NurbsCurve *pnbs = ON_NurbsCurve::New(3,true,ReconstructorPara::NurbsDegree+1, ncv ) ;
		ON_NurbsCurve nbs = *pnbs ;  delete pnbs ;
		for( int i=0; i<nbs.CVCount() - nbs.Degree(); ++i){
			ON_3dPoint p ;
			fit_pdm.m_nurbs.GetCV(i, p);
			nbs.SetCV( i, ON_4dPoint(p[0],p[1],p[2],1.0 ) ) ;

			if( i < nbs.Degree() )
				nbs.SetCV( i+ReconstructorPara::cvNum, ON_4dPoint(p[0],p[1],p[2],1.0 ) ) ;

		}

		{
			ON_NurbsCurve curveXX = pcl::on_nurbs::FittingCurve2dPDM::initNurbsCurve2D (ReconstructorPara::NurbsDegree+1, data.interior, ncv);
			pcl::on_nurbs::FittingCurve2dSDM fit_pdmXX (&data, curveXX);
			fit_pdmXX.assemble (curve_params);
			fit_pdmXX.solve ();

			for( int i=0;i<nbs.KnotCount();++i )
				nbs.SetKnot( i,fit_pdmXX.m_nurbs.Knot(i) ) ;
		}

		return nbs ;

		return fit_pdm.m_nurbs ;


	}



	//************************************
	// 函数名称: fitAndUniformNURBSforProfile2d
	// 函数说明：fit a NURBS and adjust the control point
	// 作    者：huajie.Shi
	// 日    期：2014.8.10
	// 返 回 值: ON_NurbsCurve
	// 参    数: Profile2D &prof
	//************************************
	ON_NurbsCurve fitAndUniformNURBSforProfile2d( Profile2D &prof ) {

		// convert to NURBS data structure
		pcl::on_nurbs::NurbsDataCurve2d data;
		for( int i=0; i<prof.size(); ++i )
			data.interior.push_back( Eigen::Vector2d (prof[i].X(), prof[i].Y()) );


		//  CURVE PARAMETERS 
		unsigned order = ReconstructorPara::NurbsDegree+1 ;  
		unsigned n_control_points ( ReconstructorPara::cvNum + ReconstructorPara::NurbsDegree);  

		// fit a 2 degree nurbs
		order = 3 ;
		n_control_points = ReconstructorPara::cvNum + 2 ;


		pcl::on_nurbs::FittingCurve2dSDM::Parameter curve_params;
		curve_params.smoothness = 0.1;
		curve_params.rScale = 1;

		//  CURVE FITTING 
		ON_NurbsCurve curve = pcl::on_nurbs::FittingCurve2dPDM::initNurbsCurve2D (order, data.interior, n_control_points);

		pcl::on_nurbs::FittingCurve2dSDM fit_pdm (&data, curve);
		fit_pdm.assemble (curve_params);
		fit_pdm.solve ();


		int ncv = ReconstructorPara::cvNum + ReconstructorPara::NurbsDegree ;

		ON_NurbsCurve *pnbs = ON_NurbsCurve::New(3,true,ReconstructorPara::NurbsDegree+1, ncv ) ;
		ON_NurbsCurve nbs = *pnbs ;  delete pnbs ;

		//added by huajie 
		std::vector<ON_3dPoint> pa ;
		pa.clear();
		for( int i=0; i<nbs.CVCount() - nbs.Degree(); ++i){
			ON_3dPoint p ;
			fit_pdm.m_nurbs.GetCV( i, p ) ;
			pa.push_back(p);

		}
		
		Point2f d1,d2;
		int s=pa.size();
		for(int j=0;j<s;j++){
			d1.X() = pa[   j   ].x - pa[(j-1+s)%s].x;
			d1.Y() = pa[   j   ].y - pa[(j-1+s)%s].y;
			d2.X() = pa[(j+1)%s].x - pa[(j-1+s)%s].x;
			d2.Y() = pa[(j+1)%s].y - pa[(j-1+s)%s].y;
			if( (d1.X() * d2.Y() - d1.Y() * d2.X()) > 0 ){
				pa[j].x = (pa[(j-1+s)%s].x + pa[(j+1)%s].x) / 2;
				pa[j].y = (pa[(j-1+s)%s].y + pa[(j+1)%s].y) / 2;
			}

		}

		std::vector<Point2f> profile;
		profile.clear();
		for(int i=0;i<pa.size();i++){
			profile.push_back(Point2f(pa[i].x,pa[i].y));
		}
		std::vector<Point2f> curves = ReconstructorUtility::discretizeNurbs( ReconstructorUtility::fitNURBSforProfile2d( profile ),  100 );
		Profile2D pro;
		pro.clear();
		int cvnum=nbs.CVCount() - nbs.Degree();
		int size = curves.size() / cvnum;
		for(int j=0; j< cvnum; ++j){
			pro.push_back(curves[j*size]);
		}
		//profile= ReconstructorUtility::discretizeNurbs( ReconstructorUtility::fitNURBSforProfile2d( pro ), cvnum);
		for(int i=0;i<pa.size();i++){
			pa[i].x=pro[i].X();
			pa[i].y=pro[i].Y();
		}

		for( int i=0; i<nbs.CVCount() - nbs.Degree(); ++i){
			ON_3dPoint p = pa[i] ;
			nbs.SetCV( i, ON_4dPoint(p[0],p[1],p[2],1.0 ) ) ;

			if( i < nbs.Degree() )
				nbs.SetCV( i+ReconstructorPara::cvNum, ON_4dPoint(p[0],p[1],p[2],1.0 ) ) ;

		}

		//////////////////////////////

		{
			ON_NurbsCurve curveXX = pcl::on_nurbs::FittingCurve2dPDM::initNurbsCurve2D (ReconstructorPara::NurbsDegree+1, data.interior, ncv);
			pcl::on_nurbs::FittingCurve2dSDM fit_pdmXX (&data, curveXX);
			fit_pdmXX.assemble (curve_params);
			fit_pdmXX.solve ();

			for( int i=0;i<nbs.KnotCount();++i )
				nbs.SetKnot( i,fit_pdmXX.m_nurbs.Knot(i) ) ;

		}
		return nbs ;
	}


	std::vector<Point2f> discretizeNurbs( ON_NurbsCurve nbs, int sampleNum, bool fast  ){

		//discretize nurbs ,according the sample number , and return the discretized nubs
		if( fast || 1 ){
			std::vector<Point2f> disNurbs ;
			for( int i=0; i<sampleNum;  ++i ){

				double t = i/(double)sampleNum ;

				ON_4dPoint p;
				nbs.Evaluate (t, 0, 4, p);

				disNurbs.push_back( Point2f(p.x, p.y ) ) ;
			}
			return disNurbs ;
		}


		double step = 1.0 / (sampleNum*5) ;
		std::vector<Point2f> disNurbs ;
		for( double t=0; t<1.0;  t+= step ){

			ON_4dPoint points;
			nbs.Evaluate (t, 0, 4, points);

			disNurbs.push_back( Point2f(points[0], points[1] ) ) ;
		}


		// calculate length
		double length = 0.0 ;
		for( int i=0; i<disNurbs.size()-1; ++i )
			length += ( disNurbs[i] - disNurbs[i+1]).Norm() ;
		length += (disNurbs[0]-disNurbs.back()).Norm() ;
		double lstep = length /sampleNum ;

		if( lstep < 2 )
			return disNurbs ;
		
		// down sample large step
		std::vector<Point2f> newdisNurbs ;
		int n = disNurbs.size();
		for( int i=0; i<disNurbs.size(); ++i ){
			int numstep = (disNurbs[(i+1)%n] - disNurbs[i]).Norm() /lstep ;
			if( numstep > 0){
				Point2f vec = (disNurbs[(i+1)%n] - disNurbs[i])/numstep ;
				for( int ns = 0; ns<numstep; ++ns )
					newdisNurbs.push_back(disNurbs[i]+vec*ns) ;
			}
			else
				newdisNurbs.push_back(disNurbs[i]) ;

		}
		disNurbs = newdisNurbs ;
		
		//delete too close adjacent points
		Point2f prepoint = disNurbs[0] ;
		for( int i=1; i< disNurbs.size()-1; ++i){
			if( (disNurbs[i] - prepoint).Norm() > lstep  )
				prepoint = disNurbs[i] ;
			else{
				disNurbs.erase( disNurbs.begin() + i) ;	
				i-- ;
			}
		}

		//std::cout<<disNurbs.size()<<std::endl;system("pause");

		if( disNurbs.size() == 0){

			ON_4dPoint points;
			nbs.GetCV(0, points) ;
			disNurbs.push_back( Point2f(points[0]/points[3], points[1]/points[3] ) ) ;

		}
		return disNurbs ;

	}

	ON_NurbsCurve deformNurbs( const ON_NurbsCurve nbs_in, NRTSF ngTrans, const std::vector<double>  &weights ){

		ON_NurbsCurve nbs = nbs_in ;

		if( nbs.CVCount() - nbs.Degree() != ngTrans.t.size()  ){
			std::cout<< "nurbs do not have same cvpointsnum with the transformation!!"<<std::endl;
			system("pause") ;
		}

		int cvnum = nbs.CVCount() - nbs.Degree() ;

		//change weights
		for( int i=0; i< cvnum; ++i ){
			ON_4dPoint p;
			nbs.GetCV( i, p) ;
			p.x *= weights[i] / p.w ;
			p.y *= weights[i] / p.w ;
			p.z *= weights[i] / p.w ;
			p.w *= weights[i] / p.w ;
			nbs.SetCV(i,p) ;

		}

		// transform
		ON_NurbsCurve res = nbs ;
		for( int cvid=0; cvid<cvnum; ++cvid){
			ON_4dPoint p;
			nbs.GetCV( cvid, p) ; p.x/=p.w ; p.y/=p.w ; p.z/=p.w ;

			p.x += ngTrans.t[cvid].X() ;
			p.y += ngTrans.t[cvid].Y() ;

			double x = p.x * cos( ngTrans.A ) - p.y * sin( ngTrans.A ) ;
			double y = p.x * sin( ngTrans.A ) + p.y * cos( ngTrans.A ) ;

			p.x = x + ngTrans.b.X() ;
			p.y = y + ngTrans.b.Y() ;


			p.x*=p.w ; p.y*=p.w ; p.z*=p.w ;

			res.SetCV(cvid, p) ;

			if( cvid <  nbs.Degree())
				res.SetCV(cvid + cvnum, p) ;
		}

		return res ;
	}


	double profileDis( const Profile2D &prof1,  const Profile2D &prof2 ){
		return ProfileDis_force(prof1, prof2);

		int size1 = prof1.size() ;
		int size2 = prof2.size() ;

		std::vector<std::vector<int>> neiId ;

		std::vector<double> dis1to2( size1 ) ;

		//GlobalFun::computeKNN(prof2, prof1, neiId, 1 ) ;
		ReconstructorUtility::computeFlannKNN(prof2, prof1, neiId, 1 ) ;
		for( int i=0; i<size1; ++i )
			dis1to2[i] = ( prof1[i] - prof2[ neiId[i][0] ] ).SquaredNorm() ;


		double dis = 0;
		for( int i=0; i<dis1to2.size(); ++i )
			dis += dis1to2[i];

		return sqrtf(dis) ;
	} 

	double ProfileDis_force(Profile2D prof1,  Profile2D prof2){

		double dis = 0;
		for( int i=0; i<prof1.size(); ++i ){
			double mindis = 10000;
			for( int j=0;j<prof2.size();++j ){
				double d = (prof1[i]-prof2[j]).SquaredNorm() ;
				if( d <mindis)
					mindis = d ;
			}
			dis+= mindis ;
		}

		return dis ;

	}

	double profileDis_d( Profile2D prof1,  Profile2D prof2 ){
		int size1 = prof1.size() ;
		int size2 = prof2.size() ;

		std::vector<std::vector<int>> neiId ;
		double dis = 0;


		std::vector<double> dis1to2( size1 ) ;

		GlobalFun::computeKNN(prof2, prof1, neiId, 1 ) ;
		for( int i=0; i<size1; ++i )
			dis1to2[i] = ( prof1[i] - prof2[ neiId[i][0] ] ).Norm() ;

		for( int i=0; i<dis1to2.size(); ++i )
			dis += dis1to2[i];

		std::vector<double> dis2to1( size2 ) ;


		GlobalFun::computeKNN(prof1, prof2, neiId, 1 ) ;
		for( int i=0; i<size2; ++i )
			dis2to1[i] = ( prof2[i] - prof1[ neiId[i][0] ] ).Norm() ;

		for( int i=0; i<dis2to1.size(); ++i )
			dis += dis2to1[i];

		return dis ;
	} 

	double profileDis_d( Profile3D prof1,  Profile3D prof2 ){
		int size1 = prof1.size() ;
		int size2 = prof2.size() ;

		std::vector<std::vector<int>> neiId ;
		double dis = 0;


		std::vector<double> dis1to2( size1 ) ;
		GlobalFun::computeKNN(prof2, prof1, neiId, 1 ) ;
		for( int i=0; i<size1; ++i )
			dis1to2[i] = ( prof1[i] - prof2[ neiId[i][0] ] ).Norm() ;

		for( int i=0; i<dis1to2.size(); ++i )
			dis += dis1to2[i];

		std::vector<double> dis2to1( size2 ) ;
		GlobalFun::computeKNN(prof1, prof2, neiId, 1 ) ;
		for( int i=0; i<size2; ++i )
			dis2to1[i] = ( prof2[i] - prof1[ neiId[i][0] ] ).Norm() ;

		for( int i=0; i<dis2to1.size(); ++i )
			dis += dis2to1[i];

		return dis ;
	} 



	double profileDis_360deg( Profile2D prof1,  Profile2D prof2 ){

		int size1 = prof1.size() ;
		int size2 = prof2.size() ;

		int degree_num = prof2.size()/2 ;
		int spanAngle = 360/degree_num ;

		std::vector< std::vector<Point2f> > prof11 (degree_num);
		std::vector< std::vector<Point2f> > prof22 (degree_num);

		for( int pid=0; pid<prof1.size(); ++pid ){

			Point2f p = prof1[pid] ;
			Point2f yaxis = Point2f( 0,1) ;
			double degree = acos( (std::min)( (std::max)( -1.0, (double)(yaxis* p)/p.Norm() ) , 1.0 )  ) ;
			if( p.X() < 0 )
				degree = 3.14159*2 - degree ;

			int deg = degree * 180.0 / 3.14159 ;

			deg = ((deg%360)/spanAngle)%degree_num ;

			prof11[deg].push_back( prof1[pid] ) ;
		}

		for( int pid=0; pid<prof2.size(); ++pid ){

			Point2f p = prof2[pid] ;
			Point2f yaxis = Point2f( 0,1) ;
			double degree = acos( (std::min)( (std::max)( -1.0, (double)(yaxis* p)/p.Norm() ) , 1.0 )  ) ;
			if( p.X() < 0 )
				degree = 3.14159*2 - degree ;

			int deg = degree * 180.0 / 3.14159 ;

			deg = ((deg%360)/spanAngle)%degree_num ;

			prof22[deg].push_back( prof2[pid] ) ;
		}

		double sumdis = 0; 
		
		for( int i=0; i<degree_num; ++i ){

			if( prof11[i].size() ){

				std::vector<Point2f> pwPoints;

				if( prof22[i].size() )
					pwPoints = prof22[i] ;
				else
					for( int offset = 1; offset < degree_num/2; ++offset ){
						if( prof22[(i+offset)%degree_num].size() ||  prof22[(i-offset)%degree_num].size() ){
							pwPoints.insert( pwPoints.begin(), prof22[(i+offset)%degree_num].begin(), prof22[(i+offset)%degree_num].end() ) ;
							pwPoints.insert( pwPoints.begin(), prof22[(i-offset)%degree_num].begin(), prof22[(i-offset)%degree_num].end() ) ;
							break;
						}
					}


				if( pwPoints.size() ){

					for( int pid1=0; pid1<prof11[i].size(); ++pid1 ){
						double mindis = 10000.0 ;
						for( int pid2=0; pid2<pwPoints.size(); ++pid2 ){
							double d = (prof11[i][pid1] - pwPoints[pid2]).SquaredNorm() ;
							if( d < mindis ) mindis = d ;
						}
						sumdis += mindis;
					}
				}


			}

		}

		return sumdis ;
	}

	std::vector<int> getSharpCVofNurbs(ON_NurbsCurve tmp){

		// set sharpIds

		std::vector<int> sharpIds ;

		int cvNum = ReconstructorPara::cvNum ;		
		int degree = ReconstructorPara::NurbsDegree ;

		for( int i=0; i<cvNum; ++i ){
			ON_4dPoint p;
			tmp.GetCV( i, p) ;
			if( p.w > 1.1 )
				sharpIds.push_back(i) ;		
		}
		return sharpIds ;
	}

	void resampleTemplates( std::vector<ON_NurbsCurve> &templates ){


	}


	std::vector<int2> detectCorrespondenceForTmps( ON_NurbsCurve T1, ON_NurbsCurve T2, std::vector<int> sharpId1, std::vector<int> sharpId2  ){
		
		
		std::vector<int2> corres ;

		if( sharpId1.size() == sharpId2.size() && sharpId1.size()!=0){
			// detect correspondences for sharp features

			std::vector<Point2f> sharpoint1 ;
			std::vector<Point2f> sharpoint2 ;

			for( int i=0; i<sharpId1.size(); ++i ){
				ON_4dPoint p;
				T1.GetCV( sharpId1[i], p) ;
				sharpoint1.push_back( Point2f(p.x/p.w, p.y/p.w) ) ;
				T2.GetCV( sharpId2[i], p) ;
				sharpoint2.push_back( Point2f(p.x/p.w, p.y/p.w) ) ;
			}

			int n=sharpId1.size() ;
			int bestOffset = 0;
			double minDis = 1000.0 ;

			for( int offset = 0; offset<n;++offset ){
				double sumdis = 0.0 ;
				for( int i=0;i<n; ++i )
					sumdis += (sharpoint1[i] - sharpoint2[(i+offset)%n]).SquaredNorm() ;

				if( sumdis < minDis){
					bestOffset = offset;
					minDis = sumdis ;
				}
			}

			for( int i=0;i<n;++i)
				corres.push_back( int2( sharpId1[i], sharpId2[(i+bestOffset)%n] ) ) ;

			// check if sharp feature distribute well
			n =  T1.CVCount()-T1.Degree() ;
			bool GoodDistribution = true ;
			for( int i=1; i<corres.size(); ++i )
				if( (corres[i].first-corres[i].second  +n )%n != (corres[0].first-corres[0].second  +n )%n  )
					GoodDistribution = false ;

			if( GoodDistribution ){
				int bestOffset = (corres[0].second - corres[0].first +n)%n;
				corres.clear() ;
				for( int i=0;i<n;++i)
					corres.push_back( int2( i, (i+bestOffset+n)%n ) ) ;
			}

		}else{
			// detect correspondences for all control points

			if( T1.CVCount() != T2.CVCount() )
				return corres ;

			std::vector<Point2f> cvp1 ;
			std::vector<Point2f> cvp2 ;

			int n =  T1.CVCount()-T1.Degree() ;

			for( int i=0; i<n; ++i ){
				ON_4dPoint p;
				T1.GetCV( i, p) ;
				cvp1.push_back( Point2f(p.x/p.w, p.y/p.w) ) ;
				T2.GetCV( i, p) ;
				cvp2.push_back( Point2f(p.x/p.w, p.y/p.w) ) ;
			}

			int bestOffset = 0;
			double minDis = 1000.0 ;

			for( int offset = 0; offset<n;++offset ){
				double sumdis = 0.0 ;
				for( int i=0;i<n; ++i )
					sumdis += (cvp1[i] - cvp2[(i+offset)%n]).SquaredNorm() ;

				if( sumdis < minDis){
					bestOffset = offset;
					minDis = sumdis ;
				}
			}

			for( int i=0;i<n;++i)
				corres.push_back( int2( i, (i+bestOffset+n)%n ) ) ;

		}

		return corres ;
	}


	void saveTmpOf1Brch( std::vector<ON_NurbsCurve> templates, std::string filename ){

		std::ofstream ofs( filename ) ;

		for( int i=0; i<templates.size(); ++i ){
			int cvnum = templates[i].CVCount() ;
			ofs << cvnum <<std::endl;

			// save cvs
			for( int j=0; j<cvnum; ++j ){
				ON_4dPoint p;
				templates[i].GetCV( j, p) ;
				ofs << p.x << " " <<p.y<< " " <<p.z<< " " <<p.w <<std::endl;
			}


		}

		ofs.close() ;

	}

	void loadTmpOf1Brch( std::vector<ON_NurbsCurve> &templates, std::string filename ){



		std::ifstream ifs( filename ) ;

		int cvnum ;
		int count = -1;
		while( ifs>> cvnum ){
		
			if( ++count >= templates.size() )
				break ;

			// load cvs
			for( int j=0; j<cvnum; ++j ){

				ON_4dPoint p;
				ifs >> p.x >> p.y>> p.z >> p.w ;

				templates[count].SetCV( j, p) ;
			}

		}

		ifs.close() ;

	}

	void saveTemplates(  std::vector<std::vector<ON_NurbsCurve>> templates,  std::vector<std::vector<int>> templateIdss, std::string dir ){

		std::cout << "saveTemplates "<<std::endl; 
		std::cout << "dir = "<<dir<<std::endl;


		std::string tidfname = dir + "\\Tids.txt" ;

		std::ofstream ofs( tidfname) ;

		ofs << templateIdss.size() << std::endl;

		for( int i = 0; i<templateIdss.size(); ++i ){
			ofs << templateIdss[i].size() << std::endl;
			for( int j=0; j< templateIdss[i].size(); ++j  )
				ofs << templateIdss[i][j] << std::endl;
		}

		char buff[1024] ;
		for( int i=0; i<templates.size(); ++i ){
			std::string filename = dir + std::string("\\Tmps.") + std::string( itoa(i,buff,10 ) ) + ".txt" ;
			saveTmpOf1Brch( templates[i], filename  ) ;
		}

		return ;
	} 

	void loadTemplates(  std::vector<std::vector<ON_NurbsCurve>> &templates,  std::vector<std::vector<int>> &templateIdss,  std::string dir ){

		std::cout << "loadTemplates "<<std::endl; 

		std::cout << "dir = "<<dir<<std::endl;

		std::string tidfname = dir + "\\Tids.txt" ;
		std::ifstream ifs( tidfname ) ;

		int numbrch ;
		ifs >> numbrch ;

		templateIdss.resize( numbrch ) ;

		for( int i=0; i< numbrch; ++i ){
			int tnum ;
			ifs>>tnum ;

			templateIdss[i].resize( tnum ) ;

			for( int j=0; j<tnum; ++j )
				ifs >> templateIdss[i][j] ;
		}


		templates.resize( templateIdss.size() ) ;
		char buff[1024] ;
		for( int i=0; i<templates.size(); ++i ){

			templates[i].resize( templateIdss[i].size() ) ;
			for( int id = 1; id < templates[i].size(); ++id )
				templates[i][id] = templates[i][0] ;


			std::string filename = dir + std::string("\\Tmps.") + std::string( itoa(i,buff,10 ) ) + ".txt" ;
			loadTmpOf1Brch( templates[i], filename  ) ;
		}

		return ;
	} 




	double lengthOfDisNurbs( Profile3D disNbs ){
		double len = 0;
		for( int i=0; i<disNbs.size()-1; ++i )
			len += ( disNbs[i] - disNbs[i+1] ).Norm() ;
		len += ( disNbs[0] - disNbs.back() ).Norm() ;

		return len ;
	}

	double round(double r)	{
		//r>0,return floor (r+0.5)
		//r<0,return ceil (r-0.5)
		return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
	}

	std::vector<Point3f> upsampleResult( std::vector<Profile3D> profiles, int2 scope, bool firstIsTip, bool lastIsTip  ){
		
		 std::vector<Profile3D> profiles_bk = profiles ;

		// check input profiles
		for( int i=0;i<profiles.size(); ++i ){
			std::vector<Point3f> &pts = profiles[i] ;
			for( int pid=0;pid<pts.size(); ++pid)
				if( _isnan (pts[pid][0]) || !_finite(pts[pid][0]) ||_isnan (pts[pid][1]) || !_finite(pts[pid][1]) ||_isnan (pts[pid][2]) || !_finite(pts[pid][2])  ){

					std::cout<<" profiles[i].size()=="<< profiles[i].size()<<std::endl;
					std::cout << pts[pid][0] <<","<<pts[pid][1]<<","<<pts[pid][2]<<std::endl;
					std::cout<<"profiles["<<i<<"]["<<pid<<"] is NAN!\n"<<__FILE__<<__LINE__<<std::endl;;
					system("pause") ;
				}


		}

		double upstep = ReconstructorPara::upsampleStep ;

		// make uniform distribution
		for( int i=0; i<profiles.size(); ++i ){

			if( i==0&&firstIsTip){
				Point3f center = GlobalFun::centerOfPoints( profiles[0] ) ;
				profiles[0].clear() ;
				profiles[0].push_back( center ) ;
				continue;
			}

			if( i==profiles.size()-1 && lastIsTip ){
				Point3f center = GlobalFun::centerOfPoints( profiles.back() ) ;
				profiles.back().clear() ;
				profiles.back().push_back( center ) ;
				continue;
			}

			int j0 = 0; //rand()%profiles[i].size() ;

			for( int inc=0; inc<=profiles[i].size(); ++inc ){

				if( profiles[i].size()==1 )
					break ;
				int j1 = (j0 + inc)%profiles[i].size() ;
				int j2 = (j0 + inc+1)%profiles[i].size() ;

				if( ( profiles[i][j1]-profiles[i][j2]).Norm() <upstep ){
					profiles[i].erase( profiles[i].begin()+j2) ;
					inc-- ;
				}

			}

		
		}
		

		std::vector<Point3f> points ;

	
#define  ThreadNum 6

		std::vector<Profile3D> subpoints( profiles.size()-1 ) ;

		// upsample
		for( int pid=scope.x; pid<=scope.y-1; ++pid){

			Profile3D &prof1 = profiles[pid] ;
			Profile3D &prof2 = profiles[pid+1] ;

			std::vector<std::vector<int>> nstIds ;
			GlobalFun::computeKNN( prof2, prof1, nstIds , 1)  ;

			if( pid == 0 && firstIsTip ){
				upsampleTip( subpoints[pid], profiles[pid][0],profiles[pid+1],profiles_bk[pid+1], profiles_bk[pid+2] ) ;
				continue;
			}

			if( pid == profiles.size()-2 && lastIsTip ){
				upsampleTip( subpoints[pid], profiles[pid][0],profiles[pid-1],profiles_bk[pid-1],profiles_bk[pid-2] ) ;
				continue;
			}

			// interpolate linearly
			for( int i=0; i<prof1.size(); ++i ){

				Point3f p1 = prof1[i]  ;
				Point3f p2 = prof2[nstIds[i][0]] ;

				int nstep = round( (p1-p2).Norm()/ upstep );

				Point3f v = (p2-p1) / nstep ;


				for( int j=0; j<nstep; ++j )
					subpoints[pid].push_back( p1 + v*j ) ;
			}

		}

		for( int i=0; i<subpoints.size(); ++i )
			if( subpoints[i].size() )
				points.insert( points.begin(),subpoints[i].begin(), subpoints[i].end() ) ;

		//points.insert( points.begin(), profiles.back().begin()+scope.x, profiles.back().end()+scope.y+1 ) ;

		for( int i=0;i<subpoints.size(); ++i ){
			std::vector<Point3f> &pts = subpoints[i] ;

			for( int pid=0;pid<pts.size(); ++pid)
				if( _isnan (pts[pid][0]) || !_finite(pts[pid][0]) ||_isnan (pts[pid][1]) || !_finite(pts[pid][1]) ||_isnan (pts[pid][2]) || !_finite(pts[pid][2])  ){

					std::cout<<" subpoints[i].size()=="<< subpoints[i].size()<<std::endl;
					std::cout << pts[pid][0] <<","<<pts[pid][1]<<","<<pts[pid][2]<<std::endl;
					std::cout<<"subpoints["<<i<<"]["<<pid<<"] is NAN!\n"<<__FILE__<<__LINE__<<std::endl;;
					system("pause") ;
				}

		}
		return points ;
	}

	void upsampleTip( std::vector<Point3f> &res, Point3f tip, std::vector<Point3f> neiprof,  std::vector<Point3f> neiprof_full, std::vector<Point3f> neineiprof_full ) {

		//up sample  the tip by Hermite , and save the result in the first parameter 
		double upstep = ReconstructorPara::upsampleStep ;

		res.clear() ;

		for( int i=0; i<neiprof.size(); ++i ){
			Point3f p = neiprof[i] ;
			
			int nnid = NearstPoint(neineiprof_full,p) ;
			Point3f normal = p-neineiprof_full[nnid] ;

			Point3f rad = p-GlobalFun::centerOfPoints(neiprof_full) ;

			std::vector<Point3f> sub = getHermite(p, tip, normal, rad, upstep) ;

			res.insert( res.begin(), sub.begin(), sub.end() ) ;
		}
		res.push_back( tip ) ;
	}



	bool vectorIntersect( Point2f p0, Point2f p1, Point2f p2, Point2f p3 ){
		double x,y,x0,y0,x1,y1,x2,y2,x3,y3,k1,k2;
		x0 = p0.X() ;
		y0 = p0.Y() ;
		x1 = p1.X() ;
		y1 = p1.Y() ;
		x2 = p2.X() ;
		y2 = p2.Y() ;
		x3 = p3.X() ;
		y3 = p3.Y() ;

		k1=(y0-y1)/(x0-x1);
		k2=(y2-y3)/(x2-x3);
		x=-(y0-k1*x0-(y2-k2*x2))/(k1-k2);
		y=y0+(x-x0)*k1;

		if( (Point2f(x,y)-Point2f(x3,y3)) * (Point2f(x,y)-Point2f(x2,y2)) <= 0 && (Point2f(x,y)-Point2f(x0,y0)) * (Point2f(x,y)-Point2f(x1,y1)) <=0)
			return true ;
		else
			return false ;
			
	}

	void getBernsteinCoefficient( std::vector<double> &bc, int n   ){

		bc.resize(n);

		n=n-1;

		int j,k;
		for (k=0;k<=n;k++)
		{ //compute n! / (k!*(n-k)!)
			bc[k] = 1;
			for (j = n;j>=k+1;j--){
				bc[k] *= (double)j / (double)(j-k);
			}
		}


	}

	void getBezier( std::vector<Point3f> ctrpoints,  std::vector<Point3f>  &output, int nump ){


			std::vector<double> bc ;

			const int n = ctrpoints.size() ;

			getBernsteinCoefficient(bc, n ) ;

			// store the Bezier curve
			output.clear();

			double ustep = 1.0/nump ;
			for( double u=0; u<=1.0; u+=ustep ){
				Point3f p(0,0,0) ;
				for( int k = 0; k<n; ++k ){
					p.X() += bc[k] * pow(u, k) * pow( 1-u, n-1-k) * ctrpoints[k].X() ;
					p.Y() += bc[k] * pow(u, k) * pow( 1-u, n-1-k) * ctrpoints[k].Y() ;
					p.Z() += bc[k] * pow(u, k) * pow( 1-u, n-1-k) * ctrpoints[k].Z() ;
				}	

				output.push_back( p ) ;
			}
	}

	void getBezier( std::vector<Point2f> ctrpoints,  std::vector<Point2f>  &output, int nump ){


		std::vector<double> bc ;

		const int n = ctrpoints.size() ;

		getBernsteinCoefficient(bc, n ) ;

		// store the Bezier curve
		output.clear();

		double ustep = 1.0/nump ;
		for( double u=0; u<=1.0; u+=ustep ){
			Point2f p(0,0) ;
			for( int k = 0; k<n; ++k ){
				p.X() += bc[k] * pow(u, k) * pow( 1-u, n-1-k) * ctrpoints[k].X() ;
				p.Y() += bc[k] * pow(u, k) * pow( 1-u, n-1-k) * ctrpoints[k].Y() ;
			}	

			output.push_back( p ) ;
		}
	}



	double M[4][4] = {
		-1, 3, -3, 1,
		3,-6,3,0,
		-3,0,3,0,
		1,4,1,0
	} ;

	void getUniformCubicBSpline( std::vector<Point3f> input, std::vector<Point3f> &output, int nump ){
		//  refer http://en.wikipedia.org/wiki/B-spline

		output.clear() ;

		int n=input.size() ;
		input.insert( input.begin(), input[0]) ; input[0] = input[1] ;
		input.push_back( input.back() ) ;
		input.push_back( input.back() ) ;

		for( int i=1; i+2<input.size(); ++i  ){
			double tstep = 1.0/nump * n ;
			for( double t = 0; t<1.0; t+= tstep ){

				// index of control points


				double T[4] ; 
				for( int id = 0 ; id< 3; ++id )
					T[id] = pow( t, 3.0-id ) /6.0 ;
				T[3] = 1.0/6.0 ;


				double TM[4] ;
				for( int id = 0 ; id< 4; ++id ){
					TM[id] = 0;
					for( int j = 0; j<4; ++j)
						TM[id] += T[j] * M[j][id] ;
				}

				Point3f St(0,0,0);
				for( int id = 0 ; id< 4; ++id ){
					St = St +  input[i-1 + id ] * TM[id] ;
				}

				output.push_back(St) ;

			}

		}


	}


	void getUniformCubicBSpline( std::vector<Point2f> input, std::vector<Point2f> &output, int nump ){
		//  refer http://en.wikipedia.org/wiki/B-spline

		output.clear() ;

		int n=input.size() ;
		input.insert( input.begin(), input[0]) ; input[0] = input[1] ;
		input.push_back( input.back() ) ;
		input.push_back( input.back() ) ;

		for( int i=1; i+2<input.size(); ++i  ){
			double tstep = 1.0/nump * n ;
			for( double t = 0; t<1.0; t+= tstep ){

				// index of control points


				double T[4] ; 
				for( int id = 0 ; id< 3; ++id )
					T[id] = pow( t, 3.0-id ) /6.0 ;
				T[3] = 1.0/6.0 ;


				double TM[4] ;
				for( int id = 0 ; id< 4; ++id ){
					TM[id] = 0;
					for( int j = 0; j<4; ++j)
						TM[id] += T[j] * M[j][id] ;
				}

				Point2f St(0,0);
				for( int id = 0 ; id< 4; ++id ){
					St = St +  input[i-1 + id ] * TM[id] ;
				}

				output.push_back(St) ;

			}

		}


	}

	void getNurbs( std::vector<Point3f> ctrpoints,  std::vector<Point3f>  &output, int nump, int degree ){

		int ncv = ctrpoints.size() ;

		if( degree > ncv/2 )
			degree = ncv/2 ;
		
		ON_3dPoint * cvs = new ON_3dPoint[ncv] ;
		for( int i=0; i<ncv; ++i ){
			ON_3dPoint p;
			p.x = ctrpoints[i].X() ;
			p.y = ctrpoints[i].Y() ;
			p.z = ctrpoints[i].Z() ;
			cvs[i] =  p ;
		}

		ON_NurbsCurve nbs;
		nbs.CreateClampedUniformNurbs( 3, degree+1, ncv, cvs, 1.0/(ncv-degree) ) ;

		double firstKnot = nbs.Knot(0) ;
		double lastKnot = nbs.Knot( nbs.KnotCount() - 1) ;


		output.clear() ;
		for( int i=0; i<nump; ++i ){

			double t = i/(nump-1.0) ;

			ON_4dPoint p;
			nbs.Evaluate ( t, 0, 4, p);

			output.push_back( Point3f(p.x, p.y, p.z ) ) ;
		}



	}



	void getPeriodicNurbs( std::vector<Point3f> ctrpoints,  std::vector<Point3f>  &output, int nump, int degree ){

		int ncv = ctrpoints.size() ;

		if( degree > ncv/2 )
			degree = ncv/2 ;

		ON_3dPoint * cvs = new ON_3dPoint[ncv] ;
		for( int i=0; i<ncv; ++i ){
			ON_3dPoint p;
			p.x = ctrpoints[i].X() ;
			p.y = ctrpoints[i].Y() ;
			p.z = ctrpoints[i].Z() ;
			cvs[i] =  p ;
		}

		ON_NurbsCurve nbs;
		nbs.CreatePeriodicUniformNurbs( 3, degree+1, ncv, cvs, 1.0/(ncv) ) ;


		output.clear() ;
		for( int i=0; i<nump; ++i ){

			double t = i/(nump-1.0) ;

			ON_4dPoint p;
			nbs.Evaluate ( t, 0, 4, p);

			output.push_back( Point3f(p.x, p.y, p.z ) ) ;
		}



	}







}

