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


#include "GLArea.h"
extern GLArea *glarea_global_ptr ;

#include "GlobalFunction.h"
#include "skelpath.h"
#include "Types.h"
#include <GL/GL.h>

#include <string> ;

#include "reconstructionUtility.h"
#include "reconstructorPara.h"

#include "reconstructorPara.h" 
#include "appstate.h"
double cornerWidth =400 ;
double planeWidth =0.6;
double edgeWidth =100;
Point2f profCenter ;

void skelpath::selectActiveTmp( Point2f px ) {

	//select template

	std::vector<Point3f> tmpCters ;
	for( int i=0; i<tempalteIdss[0].size(); ++i )
		tmpCters.push_back( smoothedBrchPts[0][ tempalteIdss[0][i] ] ) ;

	std::vector<Point2f> tmpCters_s2d = ReconstructorUtility::convertWorld2dToScreen2d( ReconstructorUtility::project3dPointsOntoScreen( tmpCters )) ;

	int nstId = ReconstructorUtility::NearstPoint( tmpCters_s2d, px ) ;
	if( (tmpCters_s2d[nstId]-px).Norm() < 40 ){
		if( activeTmpId == nstId ){
			activeTmpId = -1;
			activeCVID = -1 ;
		}
		else {
			activeTmpId = nstId ;
			activeCVID = -1 ;
		}
	}



	// update planeWidth & profileCter
	if( activeTmpId >= 0 ){

		if( profiles2d[0][tempalteIdss[0][activeTmpId]].size() >  5){
			Box2f box = GlobalFun::computeBox( profiles2d[0][tempalteIdss[0][activeTmpId]] ) ;
			planeWidth = 2 * std::max(  std::max( fabs(box.min.X()),  fabs(box.min.Y() ) ), std::max( fabs(box.max.X()),  fabs(box.max.Y() ) ) );
			planeWidth *= 1.3 ;

			profCenter = ( box.max + box.min ) * 0.5 ;
		}
		else {
			planeWidth = 0.3 ;
			profCenter *= 0 ;
		}
	}



	std::cout << "skelpath::selectActiveTmp:\n" << "activeTmpId = " <<activeTmpId << "\nactiveCVID = "<<activeCVID<<std::endl;

}

void skelpath::selectActiveTmpCV( Point2f px ){

	//select control point of template

	if( activeTmpId < 0 )
		return ;

	cutPlane ctpl = cutplanes[0][tempalteIdss[0][activeTmpId]] ;
	Point3f cter = smoothedBrchPts[0][ tempalteIdss[0][activeTmpId] ] ;

	std::vector<Point2f> cvs_l2d = ReconstructorUtility::getCvPts( Nurbsss[0][activeTmpId] ) ;
	std::vector<Point3f> cvs = ReconstructorUtility::convert2dProfileTo3d( cvs_l2d, ctpl, cter );

	std::vector<Point2f> cvs_s2d = ReconstructorUtility::convertWorld2dToScreen2d( ReconstructorUtility::project3dPointsOntoScreen( cvs )) ;

	int nstId = ReconstructorUtility::NearstPoint( cvs_s2d, px ) ;
	if( (cvs_s2d[nstId]-px).Norm() < 40 ){
		activeCVID = nstId ;
		cvSelectionType = 1 ;

	}
	//else
	//	activeCVID = -1 ;


	if( activeCVID == -1 ){ // select from corner

		cvs_s2d.clear();
		for( int i=0; i<cvs_l2d.size(); ++i )
			cvs_s2d.push_back(  (cvs_l2d[i] - profCenter) * cornerWidth/planeWidth + Point2f( 1,1) * (cornerWidth/2 + edgeWidth)  ) ;

		int nstId = ReconstructorUtility::NearstPoint( cvs_s2d, px ) ;
		if( (cvs_s2d[nstId]-px).Norm() < 40 ){
			activeCVID = nstId ;
			cvSelectionType = 2 ;
		}
	}


	std::cout << "skelpath::selectActiveTmpCV:\n" << "activeTmpId = " <<activeTmpId << "\nactiveCVID = "<<activeCVID<<std::endl;
}

void skelpath::moveActiveTmpCV( Point2f pixelPos ) {

	//move the control point

	std::cout << " skelpath::moveActiveTmpCV( "<<pixelPos.X()<<", "<<pixelPos.Y()<<")"<<std::endl;

	if( activeCVID<0)
		return ;

	cutPlane ctpl = cutplanes[0][tempalteIdss[0][activeTmpId]] ;
	Point3f cter = smoothedBrchPts[0][ tempalteIdss[0][activeTmpId] ] ;

	ON_NurbsCurve &nbs = Nurbsss[0][activeTmpId]  ;


	if( cvSelectionType == 1 ){

		Point3f pixelPos_l3d = ReconstructorUtility::convertWorld3dToLocal3d( ReconstructorUtility::convertScreen2dToWorld3d( std::vector<Point2f>(1,pixelPos) ) ) [0];
		Point3f original_l3d = ReconstructorUtility::convertWorld3dToLocal3d( std::vector<Point3f>(1, Point3f(0,0,0)) )[0] ;


		Point3f intersectPoint ;
		if( ReconstructorUtility::CalPlaneLineIntersectPoint(intersectPoint, ctpl.planeNormal, cter, original_l3d-pixelPos_l3d, pixelPos_l3d) ){

			Point2f dst = ReconstructorUtility::convert3dProfileTo2d( std::vector<Point3f>(1,intersectPoint),ctpl, cter )[0] ;

			std::vector<Point2f> cvs = ReconstructorUtility::getCvPts( Nurbsss[0][activeTmpId] ) ;
			cvs[activeCVID] = dst ;

			ReconstructorUtility::setCvOfNurbs( nbs, cvs) ;
			discretizeNurbs() ;
		}
	}else if( cvSelectionType == 2 ){


		std::vector<Point2f> cvs = ReconstructorUtility::getCvPts( Nurbsss[0][activeTmpId] ) ;

		Point2f dst =( pixelPos - Point2f( 1,1) * (cornerWidth/2 + edgeWidth) )  * ( planeWidth/cornerWidth)  +  profCenter ;
		cvs[activeCVID] = dst ;

		ReconstructorUtility::setCvOfNurbs( nbs, cvs) ;
		discretizeNurbs() ;
	}


}

void skelpath::deselectActiveTmpCV() {

	//release the control point 

	activeCVID = -1 ;
}
void skelpath::deselectActiveTmp() {

	//release the template

	activeTmpId = -1 ;
}

void skelpath::drawActiveTmp() {

	//draw active template

	if( activeTmpId<0)
		return ;

	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);




	cutPlane ctpl = cutplanes[0][tempalteIdss[0][activeTmpId]] ;
	Point3f cter = smoothedBrchPts[0][ tempalteIdss[0][activeTmpId] ] ;

	std::vector<Point3f> cvs = ReconstructorUtility::convert2dProfileTo3d( ReconstructorUtility::getCvPts( Nurbsss[0][activeTmpId] ), ctpl, cter );




	// draw control points
	for( int i=0; i<cvs.size(); ++i){
		if( i==activeCVID ){
			glColor4f(1,0,0, 1) ;
			glPointSize(30) ;
		}
		else{
			glColor4f(0,0,1, 1) ;
			glPointSize(20) ;
		}
		glBegin( GL_POINTS) ;
		glVertex3f( cvs[i][0],  cvs[i][1],  cvs[i][2] ) ;
		glEnd() ;
	}

	// draw line segment between control points
	glColor3f(0,0,0) ;
	glLineWidth(5) ;
	glBegin(GL_LINE_LOOP) ;
	for( int i=0; i<cvs.size(); ++i)
		glVertex3f( cvs[i][0],  cvs[i][1],  cvs[i][2] ) ;
	glEnd() ;




	//if( activeTmpId >= 0 )
	//	GlobalFun::draw3dCurves(std::vector<Curve3D>(1,disNurbsss3d[0][activeTmpId]) ,GLColor( 0.8,0, 0, 0.7), true, false, ReconstructorPara::nurbsDisplaySize ) ;

	if( activeTmpId >= 0 ){

		glColor4f(0.8,0,0, 0.7) ;
		glLineWidth(10) ;
		glBegin(GL_LINE_LOOP) ;
		cvs = disNurbsss3d[0][activeTmpId] ;
		for( int i=0; i<cvs.size(); ++i)
			glVertex3f( cvs[i][0],  cvs[i][1],  cvs[i][2] ) ;

		glEnd() ;
	}

}

void skelpath::drawActiveTmp_cutplane(){

	if( activeTmpId<0)
		return ;

	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);




	cutPlane ctpl = cutplanes[0][tempalteIdss[0][activeTmpId]] ;
	Point3f cter = smoothedBrchPts[0][ tempalteIdss[0][activeTmpId] ] ;


	// draw plane
	std::vector<Point2f> vertices;
	vertices.push_back( Point2f(1,1) ) ;
	vertices.push_back( Point2f(1,-1) ) ;
	vertices.push_back( Point2f(-1,-1) ) ;
	vertices.push_back( Point2f(-1,1) ) ;

	for( int i=0; i<vertices.size(); ++i )
		vertices[i] *= planeWidth/2 ;

	std::vector<Point3f> v3d = ReconstructorUtility::convert2dProfileTo3d(vertices, ctpl, cter );


	glDisable( GL_CULL_FACE ) ;
	glDisable(GL_LIGHTING) ;
	glColor4f(1,0,1,0.3) ;




	glBegin( GL_POLYGON ) ;
	for( int i=0; i<v3d.size(); ++i)
		glVertex3f( v3d[i][0],  v3d[i][1],  v3d[i][2] ) ;
	glEnd() ;
}

void skelpath::drawActiveTmpOnCorner(){



	if( activeTmpId<0)
		return ;


	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);


	double scaleRatio = cornerWidth / planeWidth ;


	glMatrixMode(GL_PROJECTION);

	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, glarea_global_ptr->width() ,0,glarea_global_ptr->height(), -1,1);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef( cornerWidth/2+edgeWidth, cornerWidth/2+edgeWidth, 0) ;

	glScalef(scaleRatio,scaleRatio, 1 ) ;

	glPushMatrix() ;
	glTranslatef( -profCenter.X(), -profCenter.Y(), 0) ;


	std::vector<Point2f> cvs2d =  ReconstructorUtility::getCvPts( Nurbsss[0][activeTmpId] ) ;
	std::vector<Point3f> cvs3d ; 
	for( int i=0; i<cvs2d.size(); ++i) 
		cvs3d.push_back( Point3f( cvs2d[i].X(), cvs2d[i].Y(), 0) ) ;

	// draw data points
	glColor3f(0,0,0) ;
	glPointSize(7) ;
	glBegin(GL_POINTS) ;
	for( int i=0; i<profiles2d[0][tempalteIdss[0][ activeTmpId] ].size(); ++i )
		glVertex3f( profiles2d[0][tempalteIdss[0][ activeTmpId]][i].X(), profiles2d[0][tempalteIdss[0][ activeTmpId]][i].Y(), 0 ) ;
	glEnd() ;



	// draw control points
	glDisable(GL_LIGHTING) ;
	for( int i=0; i<cvs2d.size(); ++i ){
		if( i==activeCVID )	glColor4f(1.0,0,0, 1) ;  
		else glColor4f(0,0,1,1) ;

		glPointSize( 20*(1 + 0.5*(int)( i==activeCVID )) ) ;
		glBegin(GL_POINTS);
		glVertex3f(cvs2d[i] .X(), cvs2d[i].Y(), 0.02  ) ;
		glEnd() ;
	}

	glDisable(GL_LIGHTING) ;
	glColor4f(0,0,1, 0.7);  glLineWidth(3);  
	glBegin(GL_LINE_LOOP) ;
	for( int i=0; i<cvs2d.size(); ++i )  glVertex3f(cvs2d[i].X(),cvs2d[i].Y(), 0.02 ) ;
	glEnd() ;


	// draw curves
	Curve2D  curve2d = ReconstructorUtility::discretizeNurbs(Nurbsss[0][activeTmpId],100) ;
	std::vector<Point3f> curve3d ; 
	for( int i=0; i<curve2d.size(); ++i) 
		curve3d.push_back( Point3f( curve2d[i].X(), curve2d[i].Y(), 0) ) ;

	//GlobalFun::draw3dCurves( std::vector<Curve3D>(1, curve3d ) ,GLColor( 0.8,0, 0, 0.7), true, false, ReconstructorPara::nurbsDisplaySize ) ;
	glDisable(GL_LIGHTING) ;
	glColor4f(0.8,0, 0, 0.7);  glLineWidth(20);  
	glBegin(GL_LINE_LOOP) ;
	for( int i=0; i<curve2d.size(); ++i )  glVertex3f(curve2d[i].X(),curve2d[i].Y(),0.01 ) ;
	glEnd() ;




	glPopMatrix() ;

	// draw plane
	std::vector<Point2f> vertices;
	vertices.push_back( Point2f(1,1) ) ;
	vertices.push_back( Point2f(1,-1) ) ;
	vertices.push_back( Point2f(-1,-1) ) ;
	vertices.push_back( Point2f(-1,1) ) ;

	for( int i=0; i<vertices.size(); ++i )
		vertices[i] *= planeWidth/2;

	std::vector<Point3f> v3d ;
	for( int i=0; i<vertices.size(); ++i) 
		v3d.push_back( Point3f( vertices[i].X(), vertices[i].Y(), 0) ) ;


	glDisable( GL_CULL_FACE ) ;
	glDisable(GL_LIGHTING) ;
	glColor4f(0,1,1,0.2) ;


	glBegin( GL_POLYGON ) ;
	for( int i=0; i<v3d.size(); ++i)
		glVertex3f( v3d[i][0],  v3d[i][1],  v3d[i][2] ) ;
	glEnd() ;


	//glPointSize(10) ;
	//glBegin( GL_POINTS);
	//glVertex2f(100/scaleRatio,100/scaleRatio) ;
	//glEnd() ;

	// Closing 2D
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}



void skelpath::drawSpecialView(){

}
