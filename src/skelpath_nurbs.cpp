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
#include "skelpath.h"
#include "bspline.h"
#include "Types.h"

#include <fstream>
extern std::ofstream logout ;

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <algorithm>    // std::sort

#include <GL/GL.h>

//#include <cv.h>
//#include <highgui.h>

#include <string> ;

#include "reconstructionUtility.h"
#include "reconstructorPara.h"

#include "declarations.h"


#include "reconstructorPara.h" 

void skelpath::normalize( Box3f box ) {

	//normalize to fit the box
	std::cout<<" skelpath::normalize"<<std::endl;
	std::cout << "bbox = "<<box.min.X() << " " <<box.min.Y() << " "<<box.min.Z() << "\n"
		<<box.max.X() << " "<<box.max.Y() << " "<<box.max.Z() << std::endl;

	float max_x = abs((box.min - box.max).X());
	float max_y = abs((box.min - box.max).Y());
	float max_z = abs((box.min - box.max).Z());
	float max_length = max_x > max_y ? max_x : max_y;
	max_length = max_length > max_z ? max_length : max_z;

	for(int i = 0; i < nodes.size(); i++){

		Point3f& p = nodes[i];

		p -= box.min;
		p /= max_length;

		p = (p - Point3f( 0.5*max_x/max_length,  0.5*max_y/max_length,  0.5*max_z/max_length));
	}


	convertVerticesListToGraph() ;
	convertVeticesListToBranchPoints() ;
	smoothEachBranch() ;
	calculateCutplanes() ;

}


void skelpath::fitInitialNurbs(){

	//calculate the nurbsss and save it in the member variable "Nurbsss"
	//using the function "ReconstructorUtility::fitAndUniformNURBSforProfile2d()"
	Nurbsss.clear() ;
	Nurbsss.resize( tempalteIdss.size() ) ;
	for( int i=0; i<tempalteIdss.size(); ++i ){
		Nurbsss[i].resize( tempalteIdss[i].size() ) ;

		for( int j=0; j<tempalteIdss[i].size(); ++j ){
			//modified by huajie,2014.8.10
			Nurbsss[i][j] = ReconstructorUtility::fitAndUniformNURBSforProfile2d( profiles2d[i][tempalteIdss[i][j]] ) ;
			//Nurbsss[i][j].MakePeriodicUniformKnotVector() ;
			//Nurbsss[i][j] = ReconstructorUtility::fitNURBSforProfile2d( profiles2d[i][tempalteIdss[i][j]] ) ;
		}
	}
}

void skelpath::initSharpLable(){

	//initialization the member variable "NurbsCVPixelsSharpLable"
	NurbsCVPixelsSharpLable.clear();
	NurbsCVPixelsSharpLable.resize( tempalteIdss[branchIdofTemplatetoDraw].size() ) ;

	int cvnum = Nurbsss[0][0].CVCount() ;
	for( int i=0; i<NurbsCVPixelsSharpLable.size(); ++i ){

		NurbsCVPixelsSharpLable[i] = std::vector<bool>( cvnum, false) ;

	}
}

void skelpath::discretizeNurbs() {

	disNurbsss.clear() ;
	disNurbsss.resize( Nurbsss.size() ) ;

	disNurbsss3d.clear() ;
	disNurbsss3d.resize( Nurbsss.size() ) ;


	for( int i=0; i<Nurbsss.size(); ++i ){

		disNurbsss[i].resize( Nurbsss[i].size() ) ;
		disNurbsss3d[i].resize( Nurbsss[i].size() ) ;


		for( int j=0; j<Nurbsss[i].size(); ++j ){

			// discretize nurbs
			disNurbsss[i][j] = ReconstructorUtility::discretizeNurbs( Nurbsss[i][j],200 ) ;

			// convert it to 3d curve
			disNurbsss3d[i][j] = ReconstructorUtility::convert2dProfileTo3d(disNurbsss[i][j],cutplanes[i][tempalteIdss[i][j]], smoothedBrchPts[i][tempalteIdss[i][j]] ) ;

		}


	}



}



void skelpath::convertTemplate2Pixles(){

	//get the profile of the template and convert it to pixels with " ReconstructorUtility::cvtProf2Pixels()"
	int bid = branchIdofTemplatetoDraw ;
	int n=tempalteIdss[bid].size() ;

	std::vector<Profile2D> profiles ;
	for( int i=0; i<n; ++i ){
		Profile2D profile =  profiles2d[bid][tempalteIdss[bid][i] ] ;
		profiles.push_back(profile) ;
	}

	Box2f box = GlobalFun::computeBox( profiles ) ;
	temPixels = ReconstructorUtility::cvtProf2Pixels( profiles ,box) ;
}


void skelpath::convertDisNurbs2Pixles(){

	int bid = branchIdofTemplatetoDraw ;
	int n=tempalteIdss[bid].size() ;

	std::vector<Profile2D> profiles ;
	for( int i=0; i<n; ++i ){
		Profile2D profile =  profiles2d[bid][tempalteIdss[bid][i] ] ;
		profiles.push_back(profile) ;
	}


	Box2f box = GlobalFun::computeBox( profiles ) ;


	// convert discretized nurbs to pixels
	std::vector<Profile2D> curves ;
	for( int i=0; i<disNurbsss[bid].size(); ++i ){
		Profile2D curve =  disNurbsss[bid][i] ;
		curves.push_back(curve) ;
	}
	disNurbsssPixels = ReconstructorUtility::cvtProf2Pixels(curves, box) ;


	double halfwdth = (std::max)(  (std::max)( fabs( box.max.X() ), fabs( box.min.X() ) ) , (std::max)( fabs( box.max.Y() ), fabs( box.min.Y() ) )    ) ;  
	Box2f box2  ;
	box2.min = Point2f(-halfwdth,-halfwdth) * 1.3;
	box2.max = Point2f(halfwdth,halfwdth) * 1.3;


	// convert control points to pixels
	std::vector<Profile2D> ctrlpoinss ;
	for( int i=0; i<Nurbsss[bid].size(); ++i ){

		Profile2D  ctrlpoints ; 
		int cvnum = Nurbsss[bid][i].CVCount() - Nurbsss[bid][i].Degree() ;
		for( int cpid = 0; cpid<cvnum + Nurbsss[bid][i].Degree(); ++cpid ){
			ON_4dPoint p ;
			Nurbsss[bid][i].GetCV(cpid, p) ;
			Point2f cvp = Point2f(p[0]/p[3], p[1]/p[3]) ;

			// correct out out boundary control points
			if( cvp.X() < box2.min.X() ){ cvp.X() = box2.min.X() ; p[0] =  box2.min.X()* p[3] ; }
			if( cvp.X() > box2.max.X() ){ cvp.X() = box2.max.X() ; p[0] =  box2.max.X()* p[3] ; }

			if( cvp.Y() < box2.min.Y() ){ cvp.Y() = box2.min.Y() ; p[1] =  box2.min.Y() * p[3] ; }
			if( cvp.Y() > box2.max.Y() ){ cvp.Y() = box2.max.Y() ; p[1] =  box2.max.Y() * p[3] ; }

			Nurbsss[bid][i].SetCV(cpid, p) ;
			if( cpid < Nurbsss[bid][i].Degree() )
				Nurbsss[bid][i].GetCV(cpid+cvnum, p) ;

			ctrlpoints.push_back( cvp  ) ;
		}

		ctrlpoinss.push_back(ctrlpoints) ;
	}
	NurbsCVPixels = ReconstructorUtility::cvtProf2Pixels(ctrlpoinss, box) ;


}

void skelpath::setSharpNurbsCV( int NurbsId, int CVId) {

	std::cout<<"void skelpath::setSharpNurbsCV( int NurbsId, int CVId)"<<std::endl;
	std::cout<<"NurbsId="<<NurbsId<<std::endl;
	std::cout<<"CVId="<<CVId<<std::endl;

	// set sharp feature by weighting
	NurbsCVPixelsSharpLable[NurbsId][CVId] = true ;

	double w = ReconstructorPara::weight_sharpCV ;
	ON_4dPoint p ;
	Nurbsss[branchIdofTemplatetoDraw][NurbsId].GetCV(CVId, p);  p[0] *= w ; p[1] *= w ; p[2] *= w ; p[3] *= w ;
	Nurbsss[branchIdofTemplatetoDraw][NurbsId].SetCV( CVId, p) ;
	if( CVId < ReconstructorPara::NurbsDegree )
		Nurbsss[branchIdofTemplatetoDraw][NurbsId].SetCV( CVId + ReconstructorPara::cvNum, p) ;

	std::cout << "weight: " ;
	for( int i=0; i<Nurbsss[branchIdofTemplatetoDraw].size(); ++i )
		for( int j=0; j<Nurbsss[branchIdofTemplatetoDraw][i].CVCount(); ++j ){
			ON_4dPoint p ;
			Nurbsss[branchIdofTemplatetoDraw][i].GetCV(j, p); 
			std::cout << p.w <<" " ;

		}
		std::cout<<std::endl;

		discretizeNurbs() ;
		convertDisNurbs2Pixles();

		return ;



		// set sharp feature by knot multipilcity

		// set lable 
		NurbsCVPixelsSharpLable[NurbsId][CVId] = true ;

		std::cout<<"void skelpath::setSharpNurbsCV( int NurbsId, int CVId)"<<std::endl;
		std::cout<<"NurbsId="<<NurbsId<<std::endl;
		std::cout<<"CVId="<<CVId<<std::endl;


		ON_NurbsCurve & nurbs0 = Nurbsss[branchIdofTemplatetoDraw][NurbsId] ;

		std::cout<<"degree = "<<nurbs0.Degree() <<std::endl;

		std::cout<<"before Set sharp feature:"<<std::endl;
		for( int k = 0; k<nurbs0.KnotCount(); ++k )
			std::cout <<  nurbs0.Knot(k) <<std::endl ;


		int cvnum = nurbs0.CVCount() ;

		if( CVId >= 3 && CVId < cvnum - 3){
			double midKnotValue = nurbs0.Knot(CVId+1) ;
			nurbs0.SetKnot(CVId+0, midKnotValue ) ;
			nurbs0.SetKnot(CVId+1, midKnotValue ) ;
			nurbs0.SetKnot(CVId+2, midKnotValue ) ;
		}
		else
		{
			if( CVId == 0 || CVId== cvnum-3 ){
				nurbs0.SetKnot(0, 0 ) ;
				nurbs0.SetKnot(1, 0 ) ;
				nurbs0.SetKnot(2, 0 ) ;
				nurbs0.SetKnot(cvnum-3, 1  ) ;
				nurbs0.SetKnot(cvnum-2, 1  ) ;
				nurbs0.SetKnot(cvnum-1, 1  ) ;
			}

			if( CVId == 1 || CVId== cvnum-2 ){
				nurbs0.SetKnot(1, 0 ) ;
				nurbs0.SetKnot(2, 0 ) ;
				nurbs0.SetKnot(3, 0 ) ;
				nurbs0.SetKnot(cvnum-2, 1  ) ;
				nurbs0.SetKnot(cvnum-1, 1  ) ;
				nurbs0.SetKnot(cvnum+0, 1  ) ;
			}
			if( CVId == 2 || CVId== cvnum-1 ){
				nurbs0.SetKnot(2, 0 ) ;
				nurbs0.SetKnot(3, 0 ) ;
				nurbs0.SetKnot(4, 0 ) ;
				nurbs0.SetKnot(cvnum-1, 1  ) ;
				nurbs0.SetKnot(cvnum+0, 1  ) ;
				nurbs0.SetKnot(cvnum+1, 1  ) ;
			}

		}




		discretizeNurbs() ;
		convertDisNurbs2Pixles();

		std::cout<<"\nafter Set sharp feature:"<<std::endl;
		for( int k = 0; k<nurbs0.KnotCount(); ++k )
			std::cout <<  nurbs0.Knot(k) <<std::endl ;

}



void skelpath::cvtPixels2CV(){


	int bid = branchIdofTemplatetoDraw ;
	int n=tempalteIdss[bid].size() ;

	std::vector<Profile2D> profiles ;
	for( int i=0; i<n; ++i ){
		Profile2D profile =  profiles2d[bid][tempalteIdss[bid][i] ] ;
		profiles.push_back(profile) ;
	}
	Box2f box = GlobalFun::computeBox( profiles ) ;


	// convert pixels to control points
	std::vector<Profile2D> ctrlpoinss= ReconstructorUtility::cvtPixels2Prof( NurbsCVPixels, box ) ;

	for( int i=0; i<Nurbsss[bid].size(); ++i ){

		ON::point_style cvtype = Nurbsss[bid][i].CVStyle() ;

		int cvnum = Nurbsss[bid][i].CVCount() ;
		for( int cpid = 0; cpid<cvnum; ++cpid ){

			ON_4dPoint p ;
			Nurbsss[bid][i].GetCV(cpid, p);

			p[0] = ctrlpoinss[i][cpid].X()*p[3] ;
			p[1] = ctrlpoinss[i][cpid].Y()*p[3] ;
			p[2] = 0 ;

			Nurbsss[bid][i].SetCV( cpid,cvtype, p ) ;
		}
	}

}

void skelpath::drawTemplates() {

	//draw templates using the member variable "temPixels"
	for( int i=0; i<temPixels.size(); ++i)
		GlobalFun::draw2dPointsOnScreen( temPixels[i], double3(0.5, 0.5, 1.0), 4 );

}

void skelpath::draw3dNurbs(){

	//draw 3d nurbs using member variable"disNurbsss3d"
	if( disNurbsss3d.size() > 0){

		for( int i=0; i<disNurbsss3d[0].size(); ++i )
			if( i != activeTmpId )
				GlobalFun::draw3dCurves(std::vector<Curve3D>(1,disNurbsss3d[0][i]) , GLColor( 0.8,0, 0, 1), true, false, ReconstructorPara::nurbsDisplaySize ) ;

	}

}

void skelpath::drawNurbs() {

	//draw sharp control point, control points and nurbs
	int bid = branchIdofTemplatetoDraw ;

	// draw sharp control point
	std::vector<Point2f> sharpcvpoint ;
	for( int i=0; i<NurbsCVPixels.size(); ++i){
		for( int j=0; j<NurbsCVPixels[i].size(); ++j )
			if( NurbsCVPixelsSharpLable[i][j] )
				sharpcvpoint.push_back( NurbsCVPixels[i][j] ) ;
	}

	GlobalFun::draw2dPointsOnScreen( sharpcvpoint, double3(0.0, 0.0, 0.0 ), 10 );

	// draw control points
	for( int i=0; i<NurbsCVPixels.size(); ++i){
		GlobalFun::draw2dPointsOnScreen( NurbsCVPixels[i], double3(1.0, 0.0, 0.0 ), 10 );
		GlobalFun::draw2dCurveOnScreen( NurbsCVPixels[i], double3(0.0, 0.0, 0.0 ), 1, false );
	}

	// draw nurbs
	for( int i=0; i<disNurbsssPixels.size(); ++i)
		GlobalFun::draw2dCurveOnScreen( disNurbsssPixels[i], double3(0.0, 1.0, 1.0 ), 4 );

}


void skelpath::setActivBranch(int bid ){ 	
	branchIdofTemplatetoDraw = bid  ;


	convertTemplate2Pixles() ;

	convertDisNurbs2Pixles() ;
	initSharpLable() ;

}
