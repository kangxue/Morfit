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

#include "skeleton_mul.h"
#include "sweeper.h"

#include "appstate.h"

extern double cornerWidth;
extern double planeWidth;
extern double edgeWidth ;
extern Point2f profCenter ;


void sweeper::selectCtrlProf( Curve2D stroke){

	//select the control profile , save in data member"ctrlSkelId" and "ctrlProfId"
	std::cout << "sweeper::selectCtrlProf()" <<std::endl;
	// upsample the stroke
	for( int i=0; i<stroke.size()-1; ++i )
		if( (stroke[i] - stroke[i+1]).Norm() > 2 )
			stroke.insert( stroke.begin()+i+1,stroke[i] +(stroke[i+1] - stroke[i]).Normalize() ) ;


	int SkelId = 0;
	int ProfId = 0 ;
	double minDis = 1e10 ;
	for( int sid =0; sid<settings.size(); ++sid){

		settings[sid].smoothedBrchPts[0] ;
		for( int pid = 0 ;pid <settings[sid].profiles2d[0].size(); ++pid ){


			cutPlane ctpl = settings[sid].cutplanes[0][pid] ;
			Point3f cter = settings[sid].smoothedBrchPts[0][pid] ;

			ON_NurbsCurve nbs = settings[sid].solvers[0].finalNurbs[pid] ;
			
			Curve2D nbs_s2d = ReconstructorUtility::convertWorld2dToScreen2d( ReconstructorUtility::project3dPointsOntoScreen( 
								ReconstructorUtility::convert2dProfileTo3d( ReconstructorUtility::discretizeNurbs(nbs,200), ctpl, cter ) ) ) ;
			
			std::vector<std::vector<int>> nstIds ;
			ReconstructorUtility::computeFlannKNN(nbs_s2d, stroke, nstIds, 1) ;

			double dis = 0 ;
			for( int i=0; i<stroke.size(); ++i )
				dis += (nbs_s2d[ nstIds[i][0] ] - stroke[i] ).Norm() ;


			if( dis < minDis ){
				minDis = dis ;
				SkelId = sid ;
				ProfId = pid ;
			}
		}
	}

	if( minDis/stroke.size() < 40 ){
		
		if( ctrlSkelId == SkelId && ctrlProfId == ProfId ){
			ctrlSkelId = -1 ;
			ctrlProfId = -1 ;
		}else{
			ctrlSkelId = SkelId ;
			ctrlProfId = ProfId ;
			ctrNurbs = settings[ctrlSkelId].solvers[0].finalNurbs[ctrlProfId] ;
		}
	}

	std::cout << "sweeper::selectCtrlProf:\n" << "ctrlSkelId = " <<ctrlSkelId << "\ctrlProfId = "<<ctrlProfId<<std::endl;

		// update planeWidth & profileCter
	if( ctrlSkelId >= 0 && ctrlProfId>= 0){

		Profile2D prof = settings[ctrlSkelId].solvers[0].finalProfiles[ctrlProfId] ;

		if( prof.size() >  5){
			Box2f box = GlobalFun::computeBox( prof ) ;
			planeWidth = 2 * std::max(  std::max( fabs(box.min.X()),  fabs(box.min.Y() ) ), std::max( fabs(box.max.X()),  fabs(box.max.Y() ) ) );
			planeWidth *= 1.3 ;

			profCenter = ( box.max + box.min ) * 0.5 ;
		}
		else {
			planeWidth = 0.3 ;
			profCenter *= 0 ;
		}


		planeWidth = 0.25 ;
	}
}

void sweeper::selectCtrlProfCV( Point2f px ) {
	
	//select the control cv , save in data member"ctrlCVId" and "ctrlCVSelectionType"
	std::cout << "sweeper::selectCtrlProfCV( "<<px.X() <<", "<<px.Y() << ")" <<std::endl;

	if( ctrlSkelId<0 || ctrlProfId<0 )
		return ;

	cutPlane ctpl = settings[ctrlSkelId].cutplanes[0][ctrlProfId] ;
	Point3f cter = settings[ctrlSkelId].smoothedBrchPts[0][ctrlProfId] ;

	std::vector<Point2f> cvs_l2d = ReconstructorUtility::getCvPts(ctrNurbs) ;
	std::vector<Point2f> cvs = ReconstructorUtility::convertWorld2dToScreen2d( ReconstructorUtility::project3dPointsOntoScreen( 
										ReconstructorUtility::convert2dProfileTo3d( cvs_l2d, ctpl, cter )   )) ;

	int nstId = ReconstructorUtility::NearstPoint( cvs, px ) ;
	if( (cvs[nstId] - px).Norm() < 20 ){
		ctrlCVId = nstId ;
		ctrlCVSelectionType = 1 ;
	}


	if( ctrlCVId == -1 ){ // select from corner
		
		std::vector<Point2f> cvs_s2d;
		for( int i=0; i<cvs_l2d.size(); ++i )
			cvs_s2d.push_back(  (cvs_l2d[i] - profCenter) * cornerWidth/planeWidth + Point2f( 1,1) * (cornerWidth/2 + edgeWidth)  ) ;
		
		int nstId = ReconstructorUtility::NearstPoint( cvs_s2d, px ) ;
		if( (cvs_s2d[nstId]-px).Norm() < 40 ){
			ctrlCVId = nstId ;
			ctrlCVSelectionType = 2 ;
		}
	}


	std::cout << "ctrlCVId = " << ctrlCVId <<std::endl;
}

void sweeper::moveCtrlProfCV(  Point2f pixelPos ) {

	std::cout << "sweeper::moveCtrlProfCV( "<<pixelPos.X() <<", "<<pixelPos.Y() << ")" <<std::endl;

	if( ctrlSkelId<0 || ctrlProfId<0 || ctrlCVId<0 )
		return ;

	cutPlane ctpl = settings[ctrlSkelId].cutplanes[0][ctrlProfId] ;
	Point3f cter = settings[ctrlSkelId].smoothedBrchPts[0][ctrlProfId] ;


	Point3f pixelPos_l3d = ReconstructorUtility::convertWorld3dToLocal3d( ReconstructorUtility::convertScreen2dToWorld3d( std::vector<Point2f>(1,pixelPos) ) ) [0];
	Point3f original_l3d = ReconstructorUtility::convertWorld3dToLocal3d( std::vector<Point3f>(1, Point3f(0,0,0)) )[0] ;


	if( ctrlCVSelectionType == 1 ){
		Point3f intersectPoint ;
		if( ReconstructorUtility::CalPlaneLineIntersectPoint(intersectPoint, ctpl.planeNormal, cter, original_l3d-pixelPos_l3d, pixelPos_l3d) ){

			Point2f dst = ReconstructorUtility::convert3dProfileTo2d( std::vector<Point3f>(1,intersectPoint),ctpl, cter )[0] ;

			std::vector<Point2f> cvs =  ReconstructorUtility::getCvPts(ctrNurbs) ;

			cvs[ctrlCVId] = dst ;

			ReconstructorUtility::setCvOfNurbs( ctrNurbs, cvs) ;
		}
	}else if(  ctrlCVSelectionType == 2 ){
		
		std::vector<Point2f> cvs = ReconstructorUtility::getCvPts(ctrNurbs ) ;

		Point2f dst =( pixelPos - Point2f( 1,1) * (cornerWidth/2 + edgeWidth) )  * ( planeWidth/cornerWidth)  +  profCenter ;
		cvs[ctrlCVId] = dst ;
		ReconstructorUtility::setCvOfNurbs( ctrNurbs, cvs) ;
	}

}

void sweeper::drawCtrlNurbs() {
	//std::cout << "sweeper::drawCtrlNurbs()" <<std::endl;


	if( ctrlSkelId<0 || ctrlProfId<0 )
		return ;


	cutPlane ctpl = settings[ctrlSkelId].cutplanes[0][ctrlProfId] ;
	Point3f cter = settings[ctrlSkelId].smoothedBrchPts[0][ctrlProfId] ;

	std::vector<Point3f> cvs = ReconstructorUtility::convert2dProfileTo3d( ReconstructorUtility::getCvPts(ctrNurbs), ctpl, cter ) ;




	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);




	// draw control points
	glEnable( GL_LIGHTING ) ;
	glEnable( GL_LIGHT0) ;
	glEnable( GL_LIGHT1) ;
	glDisable(GL_LIGHTING) ;
	for( int i=0; i<cvs.size(); ++i ){
		if( i==ctrlCVId )	glColor4f(1.0,0,0, 1) ;  
		else glColor4f(0,0,1,1) ;

		GLUquadricObj *objCylinder = gluNewQuadric();
		glPushMatrix() ;
		glTranslatef(cvs[i] .X(), cvs[i].Y(), cvs[i].Z() ) ;
		gluSphere(objCylinder, 0.006*(1 + 0.5*(int)( i==ctrlCVId )), 20, 10 );
		glPopMatrix() ;
		gluDeleteQuadric( objCylinder ) ;
	}

	GlobalFun::draw3dCurves( std::vector<Curve3D>(1,cvs),GLColor(0,0,1), true, false, 0.5) ;
	glDisable(GL_LIGHTING) ;


	// draw 3d nurbs
	GlobalFun::draw3dCurves(std::vector<Curve3D>(1,ReconstructorUtility::convert2dProfileTo3d( ReconstructorUtility::discretizeNurbs(ctrNurbs,100) , ctpl, cter)), GLColor( 0.8,0, 0, 1), true, false, ReconstructorPara::nurbsDisplaySize) ;


	
}


void sweeper::drawCtrlNurbs_cutplane() {


	if( ctrlSkelId<0 || ctrlProfId<0 )
		return ;


	cutPlane ctpl = settings[ctrlSkelId].cutplanes[0][ctrlProfId] ;
	Point3f cter = settings[ctrlSkelId].smoothedBrchPts[0][ctrlProfId] ;

	std::vector<Point3f> cvs = ReconstructorUtility::convert2dProfileTo3d( ReconstructorUtility::getCvPts(ctrNurbs), ctpl, cter ) ;

	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);



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


void sweeper::drawCtrlNurbsOnCorner(){



	if( ctrlSkelId<0 || ctrlProfId<0 )
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

	//// draw


	std::vector<Point2f> cvs2d =  ReconstructorUtility::getCvPts(ctrNurbs);

	

	// draw data points
	glColor3f(0,0,0) ;
	glPointSize(7) ;
	glBegin(GL_POINTS) ;
	
	Curve2D prof2d = settings[ctrlSkelId].profiles2d[0][ctrlProfId] ;
	
	for( int i=0; i<prof2d.size(); ++i )
		glVertex3f( prof2d[i].X(), prof2d[i].Y(), 0 ) ;
	glEnd() ;



	// draw control points
	glDisable(GL_LIGHTING) ;
	for( int i=0; i<cvs2d.size(); ++i ){
		if( i==ctrlCVId )	glColor4f(1.0,0,0, 1) ;  
		else glColor4f(0,0,1,1) ;

		glPointSize( 20*(1 + 0.5*(int)( i==ctrlCVId )) ) ;
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
	Curve2D  curve2d = ReconstructorUtility::discretizeNurbs( ctrNurbs,100) ;
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


	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

}

void sweeper::startEditResultStroke( Point2f p ) {


	std::cout << "sweeper::startEditResultStroke()" <<std::endl;

	editResultStroke.clear() ;
	editResultStroke.push_back( p ) ;
}

void sweeper::moveEditResultStroke( Point2f p ) {

	std::cout << "sweeper::moveEditResultStroke()" <<std::endl;

	if( (editResultStroke.back() - p).Norm() > 5 )
		editResultStroke.push_back( p ) ;
}

void sweeper::EndEditResultStroke( ) {
	std::cout << "sweeper::EndEditResultStroke()" <<std::endl;


	if( editResultStroke.size() > 2 )
		selectCtrlProf(editResultStroke) ;

	editResultStroke.clear() ;
}


void sweeper::updateResult() {

	std::cout << "sweeper::updateResult()" <<std::endl;

	settings[ctrlSkelId].updateResult( ctrlProfId, ctrNurbs) ;

	meshes[ctrlSkelId] = settings[ctrlSkelId].resultMesh[0] ;

	ongoingSkel = settings[ctrlSkelId] ;


	mergeMesh() ;

} 