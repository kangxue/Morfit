
#include "GLArea.h"
#include <GL/glu.h>
#include <stack>
#include <QCursor>
#include <QGuiApplication>
#include <QMessageBox>
#include "appstate.h"
#include "reconstructorPara.h"
#include <fstream>
#include <QColorDialog>
#include <QMenu>
extern std::ofstream logout ;
extern AppState appstate ;

#include "reconstructionUtility.h"
#include "declarations.h"
#include "skeleton_mul.h"
std::vector<std::pair<Point3f, Point3f> > globalVectorsToDraw  ; 

std::vector<Point3f > globalPointsToDraw1  ;  ;
std::vector<Point3f > globalPointsToDraw2  ;  ;


GLdouble globalProjectMatrix[16] ;
GLdouble globalModelviewMatrix[16] ;

std::vector<std::pair<Point3f, GLColor> >  normaDiffFieldToDraw ;



void  drawPoints( std::vector<Point3f >& points, GLColor color, int size = 10 ){
	glDisable(GL_LIGHTING) ;

	glColor4f(color.r, color.g, color.b, 1 ) ;

	glPointSize(size) ;

	glBegin(GL_POINTS ) ;

	for( int i=0; i<points.size(); ++i )
		glVertex3f(points[i][0],points[i][1],points[i][2] ) ;
	glEnd() ;
} 

void  drawPoint( Point3f & point, GLColor color ){
	glDisable(GL_LIGHTING) ;

	glColor4f(color.r, color.g, color.b, 1 ) ;

	glPointSize(8) ;

	glBegin(GL_POINTS ) ;

	glVertex3f(point[0], point[1], point[2] ) ;

	glEnd() ;
} 


void drawAxises(){

	glLineWidth(2) ;

	glColor4f(1,0,0, 1 ) ;
	glBegin( GL_LINES ) ;
	glVertex3f(0,0,0) ;
	glVertex3f(1,0,0) ;
	glEnd() ;

	glColor4f(0,1,0, 1 ) ;
	glBegin( GL_LINES ) ;
	glVertex3f(0,0,0) ;
	glVertex3f(0,1,0) ;
	glEnd() ;

	glColor4f(0,0,1, 1 ) ;
	glBegin( GL_LINES ) ;
	glVertex3f(0,0,0) ;
	glVertex3f(0,0,1) ;
	glEnd() ;
}
//added by huajie
extern	std::vector<Point3f> dsts;
extern  std::vector<Point3f> fraAroundTip;
extern std::vector<Point3f> centerOfGravitys;
extern std::vector<Point3f> centerTip;
extern std::vector<Point3f> allNmls;
extern std::vector<Point3f> edges;
void GLArea::paintGL() 
{
	bool orthProject = false ;

	//define some colors
	GLColor randColor[10] ;
	randColor[0] = GLColor( 0.5,0, 0.5) ;
	randColor[1] = GLColor( 0,  1,  0) ;
	randColor[2] = GLColor( 0,  0,1.0) ;
	randColor[3] = GLColor( 1.0,0.7,0.5) ;
	randColor[4] = GLColor( 1.0,1.0,0.0) ;
	randColor[5] = GLColor( 0.0,1.0,1.0) ;
	randColor[6] = GLColor( 1,0,0) ;
	randColor[7] = GLColor( 0.5,1.0,0.5) ;
	randColor[8] = GLColor( 1,0.5,1.0) ;
	randColor[9] = GLColor( 0.7,0.5,0.7) ;

	glEnable(GL_DEPTH_TEST) ;

	if( appstate.cullFace )
		glEnable(GL_CULL_FACE) ;
	else
		glDisable(GL_CULL_FACE) ;

	glCullFace(GL_BACK) ;

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	paintMutex.lock();{

		if (is_paintGL_locked)
		{
			goto PAINT_RETURN;
		}

		lightOnOff(para->getBool("Light On or Off"));

		GLColor color(global_paraMgr.drawer.getColor("Background Color"));
		//glClearColor(		color.r, color.g, color.b, 1); 

		Point3f lightPos = para->getPoint3f("Light Position");
		float lpos[4];
		lpos[0] = lightPos[0];
		lpos[1] = lightPos[1];
		lpos[2] = lightPos[2];
		lpos[3] = 0;
		glLightfv(GL_LIGHT0, GL_POSITION, lpos);       
		lpos[0] = -lightPos[0];
		lpos[1] = -lightPos[1];
		lpos[2] = -lightPos[2];
		lpos[3] = 0;
		glLightfv(GL_LIGHT1, GL_POSITION, lpos);     


		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);   

		double SnapResolutionScale = global_paraMgr.glarea.getDouble("Snapshot Resolution");
		if (takeSnapTile)
		{
			double normal_width = global_paraMgr.drawer.getDouble("Normal Line Width");
			double dot_size = global_paraMgr.drawer.getDouble("Sample Dot Size");
			double original_dot_size = global_paraMgr.drawer.getDouble("Original Dot Size");

			global_paraMgr.drawer.setValue("Normal Line Width", DoubleValue(normal_width * SnapResolutionScale * SnapResolutionScale * snapDrawScal));
			global_paraMgr.drawer.setValue("Sample Dot Size", DoubleValue(dot_size * SnapResolutionScale * SnapResolutionScale * snapDrawScal));
			global_paraMgr.drawer.setValue("Original Dot Size", DoubleValue(original_dot_size * SnapResolutionScale * SnapResolutionScale * snapDrawScal));
		}

		glLoadIdentity();
		gluLookAt(0, 0, -3, 0, 0, 0, 0, 1, 0); 

		camaraDistForDrawCurve = -3 ;

		if (takeSnapTile)
		{
			setView();// high resolution snapshot
		}


		//Drawing the scene 
		glPushMatrix();
		trackball.GetView();

		Point3f c = -gl_box.Center();
		double radius = 2.0f/gl_box.Diag();

		trackball.Apply(false);

		glPushMatrix();
		glScalef(radius, radius, radius);
		glTranslatef(c[0], c[1], c[2]);

		View<float> view;
		view.GetView();
		Point3f viewpoint = view.ViewPoint();
		glDrawer.setViewPoint(viewpoint);


		if( mvmLoaded )
			glLoadMatrixf(mvmatrix_loaded);

		if( appstate.cameraRotating ){

			// axis set for  dancers
			glTranslatef( -0.015, -0.432, +0.006) ;
			glRotatef( appstate.cameraRotateDegree *180/3.14, -0.03646, 0.9993, 0.0064 ) ;
			glTranslatef( 0.015, 0.432, -0.006) ;


		}

		float oldMVMatrix[16] ; for( int i=0;i<16; ++i ) oldMVMatrix[i] = modelviewMatrix[i] ;
		glGetFloatv(GL_MODELVIEW_MATRIX, modelviewMatrix);
		glGetFloatv(GL_PROJECTION_MATRIX, projectMatrix);

		glGetDoublev(GL_MODELVIEW_MATRIX, globalModelviewMatrix);
		glGetDoublev(GL_PROJECTION_MATRIX, globalProjectMatrix);


#ifndef _DEBUG

//		if( !appstate.cameraRotating )
			if( dataMgr.swp.mergedMesh.nfaces > 0)
				for( int i=0; i<16; ++i )
					if( oldMVMatrix[i]!= modelviewMatrix[i])
						dataMgr.swp.mergedMesh.sortDepth() ; 
#endif

		//// draw points
		//bindShader("Oren-Nayar") ;
		//glDrawer.drawPolygonOnScreen( dataMgr.getCurrentOriginal(), global_paraMgr.drawer.getDouble("Original Dot Size")*0.02  * (takeSnapTile?0.25:1), 4 ) ;
		//releaseShader() ;

		////added by huajie
		//	drawPoints( dsts, GLColor(255,0,0), 10 );
		//	drawPoints( fraAroundTip, GLColor(0,0,255), 3 );
		//	drawPoints(centerOfGravitys,GLColor(0,255,255),10);
		//	drawPoints(centerTip,GLColor(255,0,0),15);
		//	glLineWidth(2) ;
		//	for(int i=0; i<fraAroundTip.size(); i++){
		//		glColor4f(0,1,0, 1 ) ;
		//		glBegin( GL_LINES ) ;
		//		glVertex3f( fraAroundTip[i].X(), fraAroundTip[i].Y(), fraAroundTip[i].Z() ) ;
		//		glVertex3f( fraAroundTip[i].X() + 0.01*allNmls[i].X(), fraAroundTip[i].Y() + 0.01*allNmls[i].Y(), fraAroundTip[i].Z() + 0.01*allNmls[i].Z() ) ;
		//		glEnd() ;
		//		glPointSize(5);
		//		glColor4f(1,0,0, 1 ) ;
		//		glBegin(GL_POINTS); 
		//		glVertex3f( fraAroundTip[i].X() + 0.01*allNmls[i].X(), fraAroundTip[i].Y() + 0.01*allNmls[i].Y(), fraAroundTip[i].Z() + 0.01*allNmls[i].Z() ) ;
		//		glEnd();
		//	}
		//	drawPoints(edges,GLColor(255,0,0),20);

		//// draw curves
		//for( int i=0; i<greenCurves.size(); ++i )
		//	GlobalFun::draw2dCurveOnScreen( greenCurves[i], double3(0,1,0) ,6,false ) ;
		//for( int i=0; i<blueCurves.size(); ++i )
		//	GlobalFun::draw2dCurveOnScreen( blueCurves[i],  double3(0,0.5,1) ,6,false ) ;



		if( (get_operating_state()!=onesweepReady && get_operating_state()!=onesweepUpdateResult) ||  
			(get_operating_state()==onesweepReady&&dataMgr.swp.ongoingSkel.activeTmpId < 0) || 
			(get_operating_state()==onesweepUpdateResult&&dataMgr.swp.ctrlSkelId<0) ){

			if( 0 || !appstate.displayProfileConfidence ){

				std::vector<Point3f> pset1, pset2 ;
				for( int i=0; i<dataMgr.swp.points.size(); ++i )
					if( dataMgr.swp.ptsSelected[i] )  
						pset2.push_back( dataMgr.swp.points[i] ) ;
					else 
						pset1.push_back( dataMgr.swp.points[i] ) ;
				
				drawPoints(pset1, GLColor(0,0,0), global_paraMgr.drawer.getDouble("Original Dot Size") ) ;	
				if( get_operating_state() == onesweepSegment )
					drawPoints(pset2, GLColor(0,0.1,0), global_paraMgr.drawer.getDouble("Original Dot Size")  * 3 ) ;
				else
					drawPoints(pset2, GLColor(0,0.1,0), global_paraMgr.drawer.getDouble("Original Dot Size")  * 3 ) ;

			}else if( dataMgr.swp.ongoingSkel.profiles2d.size()==1){
				//std::vector<Profile3D> profs3d =  dataMgr.swp.ongoingSkel.profiles3d[0] ;
				//std::vector<double> profconf =  dataMgr.swp.ongoingSkel.profile_conf[0] ;
				//for( int i=0; i<profs3d.size(); ++i ){
				//	double3 color = GlobalFun::scalar2color(1.0 - profconf[i]) ;
				//	drawPoints(profs3d[i], GLColor(color.x, color.y, color.z),global_paraMgr.drawer.getDouble("Original Dot Size")  ) ;
				//}
			}

		}else if(get_operating_state()==onesweepReady&&dataMgr.swp.ongoingSkel.activeTmpId>=0){

			glDisable( GL_LIGHTING) ;

			skelpath &skel = dataMgr.swp.ongoingSkel ;

			for( int i=0; i<skel.profiles3d[0].size(); ++i ){
				if( i==skel.tempalteIdss[0][skel.activeTmpId] ){
					glColor4f(0,0,0,1.0) ;
					glPointSize( global_paraMgr.drawer.getDouble("Original Dot Size")*4) ;
				}
				else{
					glColor4f(0.0,0.0,0.0,0.7) ;
					glPointSize( global_paraMgr.drawer.getDouble("Original Dot Size")) ;

				}

				for( int j=0; j<skel.profiles3d[0][i].size();++j ){
					Point3f p = skel.profiles3d[0][i][j] ;
					glBegin( GL_POINTS ) ;
					glVertex3f(p.X(), p.Y(), p.Z() ) ;
					glEnd() ;
				}
			}

			if( !appstate.displayProfileConfidence )
				drawPoints(dataMgr.swp.points, GLColor(0,0,0), 1) ;

		}else if(get_operating_state()==onesweepUpdateResult&&dataMgr.swp.ctrlSkelId>=0 &&  !appstate.displayProfileConfidence ){
			glDisable( GL_LIGHTING) ;
			skelpath &skeltoprint = dataMgr.swp.settings[dataMgr.swp.ctrlSkelId] ;
			for( int i=0; i<skeltoprint.profiles3d[0].size(); ++i ){
				if( i==dataMgr.swp.ctrlProfId ){
					glColor4f(0,0,0,1.0) ;
					glPointSize( global_paraMgr.drawer.getDouble("Original Dot Size")*4 ) ;
				}
				else{
					glColor4f(0.0,0.0,0.0,0.7) ;
					glPointSize( global_paraMgr.drawer.getDouble("Original Dot Size")) ;


				}

				for( int j=0; j<skeltoprint.profiles3d[0][i].size();++j ){
					Point3f p = skeltoprint.profiles3d[0][i][j] ;
					glBegin( GL_POINTS ) ;
					glVertex3f(p.X(), p.Y(), p.Z() ) ;
					glEnd() ;
				}
			}

			drawPoints(dataMgr.swp.points, GLColor(0,0,0), 1 ) ;
		}


		// draw skeleton
		if( get_operating_state() != onesweepReady && get_operating_state()!=onesweepCorrespondence){
			//bindShader("xray") ;

			if( appstate.displaySkeleton ){

				GLColor color1 = GLColor(101/255.0, 152/255.0, 1, 1) ;
				GLColor color2 = GLColor(0.8,0,0,1) ;
				GLColor color3 = GLColor(0xcd/255.0, 0x7f/255.0 ,  0x32/255.0,1) ;
				GLColor color4 = GLColor(0x8e/255.0, 0x23/255.0 ,  0x6b/255.0,1) ;
				GLColor color5 = GLColor( 0, 0.8,0,1) ;
				GLColor color6 = GLColor( 0.5, 0,0.5,1) ;

				GlobalFun::draw3dCurves( dataMgr.skel.smoothedBrchPts,color2 , false, false, ReconstructorPara::skelDisplaySize * 0.7 ) ;

				std::vector<Point3f> joints;
				std::vector<std::vector<int>> ajBrchId ;
				dataMgr.skel.getJoints(joints, ajBrchId)  ;

				ReconstructorUtility::drawPointsAsSphere( joints, color2, 8* 0.7 ) ;

			}

			//if( appstate.displaySkeleton )
			//	for( int i=0; i< dataMgr.skel.smoothedBrchPts.size(); ++i )
			//		GlobalFun::draw2dCurveOnScreen( ReconstructorUtility::convertWorld2dToScreen2d( ReconstructorUtility::project3dPointsOntoScreen(dataMgr.skel.smoothedBrchPts[i])),
			//											GLColor(101/255.0, 152/255.0, 1, 0.6), 6, false ) ;



			//dataMgr.skel.drawBrchesInDiffColor() ;
			//releaseShader() ;
		}


		if(   appstate.displaySettings ){

			for( int i=0; i<appstate.skelProfs.size(); ++i ){

				GLColor color = randColor[i%10] ;
				GlobalFun::draw3dCurves( dataMgr.skel.smoothedBrchPts,GLColor(0,0.8,0) , false, false, ReconstructorPara::skelDisplaySize*0.6  ) ;
				std::vector<Point3f> joints;
				std::vector<std::vector<int>> ajBrchId ;
				dataMgr.skel.getJoints(joints, ajBrchId)  ;
				ReconstructorUtility::drawPointsAsSphere( joints, GLColor(0,0.8,0), 8* 0.7 ) ;


				color = GLColor(0.8,0,0) ;

				GlobalFun::draw3dCurves(appstate.skelProfs[i].disNurbsss3d[0]  , color, true, false, ReconstructorPara::nurbsDisplaySize*0.6) ;


			}


		}


		if(  appstate.displayMeshDiffColor ){

			
			appstate.resultMesh.drawItself( appstate.MeshDiffColors ) ;

		}



		if( get_operating_state() == onesweepReady || get_operating_state() == onesweepCorrespondence ){

			if( appstate.displaySkeleton )
				dataMgr.swp.ongoingSkel.drawBrchesInDiffColor( GLColor(0, 0.8,0,1) );
			
			if( appstate.displayCutRadii)
				dataMgr.swp.ongoingSkel.drawCutRadii() ;


			if( appstate.displaySeedProfile )
				dataMgr.swp.ongoingSkel.draw3dNurbs() ;
			
			dataMgr.swp.ongoingSkel.drawActiveTmp() ;
			dataMgr.swp.ongoingSkel.drawActiveTmpOnCorner () ;
			

			//dataMgr.swp.ongoingSkel.drawTemplates() ;
			//dataMgr.swp.ongoingSkel.drawNurbs() ;
			//GlobalFun::draw3dCurves( std::vector<Curve3D>(1,dataMgr.swp.featureStroke3d), GLColor(0.98,0.5,0.04), false, false, 1.0 ) ;
			//GlobalFun::draw2dCurveOnScreen(ReconstructorUtility::convertWorld2dToScreen2d( ReconstructorUtility::project3dPointsOntoScreen(dataMgr.swp.ongoingSkel.smoothedBrchPts[0])),double3(0,0,0), 3, false  ) ;

		}

		if( get_operating_state() == onesweep ){
			GlobalFun::draw2dCurveOnScreen( dataMgr.swp.ongoingSkelStroke, double3(0,1,0) ,6,false ) ;
			GlobalFun::draw2dCurveOnScreen( dataMgr.swp.ongoingFeatureStroke, double3(0.98,0.5,0.04) ,6,false ) ;
			//GlobalFun::draw3dCurves( std::vector<Curve3D>(1,dataMgr.swp.ongoingSkelStroke3d), GLColor(1.0,0,0), false, false, 1.0 ) ;
		}



		if( get_operating_state() == onesweepUpdateResult ){
			dataMgr.swp.drawCtrlNurbs() ;
			dataMgr.swp.drawCtrlNurbsOnCorner() ;
			
			// draw stroke
			GlobalFun::draw2dCurveOnScreen( dataMgr.swp.editResultStroke, double3(0,0,1), 3, false ) ;
		}


		if( get_operating_state() == onesweep || get_operating_state() == onesweepReady || get_operating_state() == onesweepUpdateResult ){
/*
			if( appstate.displayTrajectory && dataMgr.swp.ongoingSkel.solvers.size() )
				dataMgr.swp.ongoingSkel.solvers[0].drawFinalCVcorres() ;*/


			if( appstate.displayTrajectory  )
				for( int i=0; i<dataMgr.swp.settings.size(); ++i )
					dataMgr.swp.settings[i].solvers[0].drawFinalCVcorres() ;
		}

		if( get_operating_state() == onesweepCorrespondence ){
			if( appstate.displayTrajectory )
				GlobalFun::draw3dCurves( dataMgr.swp.corresTrajectory, GLColor(0,0,0), false, false, ReconstructorPara::trajDisplaySize ) ;
		}


		if( appstate.displayResultProfile ) {

			// draw profiles

			appstate.sweepingAnimation_meshNum = 0 ;

			if( get_operating_state() != onesweepCorrespondence ){

				int profileIdCount = 0 ;


				int n = 0 ;

				for( int sid=0; sid<dataMgr.swp.settings.size(); ++sid ){
					std::vector<std::vector<Profile3D>> &prof3ds = dataMgr.swp.settings[sid].resultProf3d ; 
					for( int i=0; i<prof3ds.size(); ++i )
						n += prof3ds[i].size() ;
				}

				for( int sid=0; sid<dataMgr.swp.settings.size(); ++sid ){
					std::vector<std::vector<Profile3D>> &prof3ds = dataMgr.swp.settings[sid].resultProf3d ; 


					std::vector<Profile3D> p3ds;

					for( int i=prof3ds[0].size()-1; i>=0; i--  )
						p3ds.push_back( prof3ds[0][i] ) ;

					//for( int i=0; i<p3ds.size(); ++i ){
					//	Profile3D newprof ;
					//	for( int j=0; j<p3ds[i].size(); j+=2 )
					//		newprof.push_back( p3ds[i][j] ) ;
					//	p3ds[i] = newprof ;
					//}


					// invert sequence of profiles for dinosaur
					int pnum = p3ds.size() ;


					if( appstate.sweepAnimation &&  profileIdCount < appstate.sweepingId &&  profileIdCount + p3ds.size()> appstate.sweepingId )
						p3ds.erase( p3ds.begin() + (appstate.sweepingId -profileIdCount),  p3ds.end()   ) ;





					//profileIdCount += p3ds.size() ;
					
					
					// add for dispear together -----------------------
					n = 0;
					for( int sid=0; sid<dataMgr.swp.settings.size(); ++sid ){
						std::vector<std::vector<Profile3D>> &prof3ds = dataMgr.swp.settings[sid].resultProf3d ; 
						for( int i=0; i<prof3ds.size(); ++i )
							 if( n< prof3ds[i].size() )
								 n = prof3ds[i].size() ;
					}
					// --------------------------------------------



					double ridius = 1.0 ;
					
					
					if( appstate.sweepingId > n )
						ridius = std::max( 0.0, 1.0 - (appstate.sweepingId-n)/10.0  ) ;



					

					if( radius > 0.001 )
						GlobalFun::draw3dCurves( p3ds, GLColor(0.1,0.1,0.1,1), true, false, ridius, 2 ) ;

					if( sid!= 0 )
						appstate.sweepingAnimation_meshNum += 1 ;
					if( sid == dataMgr.swp.settings.size()-1 && profileIdCount < appstate.sweepingId )
						appstate.sweepingAnimation_meshNum += 1 ;



					// add for dispear together ---------------------
					if(  appstate.sweepingId > n)
						appstate.sweepingAnimation_meshNum = dataMgr.swp.settings.size() ;
					else
						appstate.sweepingAnimation_meshNum = 0 ;
					//std::cout << "radius = " << radius <<std::endl;
					//std::cout << "sweepingAnimation_meshNum = " << appstate.sweepingAnimation_meshNum <<std::endl;
					//std::cout << "dataMgr.swp.settings.size() = " << dataMgr.swp.settings.size() <<std::endl;
					// --------------------------------------------

					if( appstate.sweepAnimation && profileIdCount >= appstate.sweepingId )
						break ;
				}



			}else{


				int profileIdCount = 0 ;

				std::vector<Profile3D> boundaryProfiles ;
				for( int sid=0; sid<dataMgr.swp.settings.size(); ++sid ){

					std::vector<std::vector<Profile3D>> &prof3ds = dataMgr.swp.settings[sid].resultProf3d ; 

					 int TmpIdx;
					 
					 if( sid == 0 )
						 TmpIdx = dataMgr.swp.ongoingSkel.tempalteIdss[0][sid] ;
					 else
						 TmpIdx = dataMgr.swp.ongoingSkel.tempalteIdss[0][sid] - (dataMgr.swp.ongoingSkel.tempalteIdss[0][sid] + dataMgr.swp.ongoingSkel.tempalteIdss[0][sid-1])/2;


					std::vector<Profile3D> p3ds_1, p3ds_2;

					for( int i=TmpIdx; i< prof3ds[0].size(); ++i )
						p3ds_1.push_back( prof3ds[0][i]  ) ;

					for( int i=TmpIdx; i>= 0; --i )
						p3ds_2.push_back( prof3ds[0][i]  ) ;



					if( appstate.sweepAnimation &&  profileIdCount < appstate.sweepingId &&  profileIdCount + p3ds_1.size()> appstate.sweepingId )
						p3ds_1.erase( p3ds_1.begin() + (appstate.sweepingId -profileIdCount),  p3ds_1.end()   ) ;

					if( appstate.sweepAnimation &&  profileIdCount < appstate.sweepingId &&  profileIdCount + p3ds_2.size()> appstate.sweepingId )
						p3ds_2.erase( p3ds_2.begin() + (appstate.sweepingId -profileIdCount),  p3ds_2.end()   ) ;

					if( sid != dataMgr.swp.settings.size()-1 && ( p3ds_1.size() < (appstate.sweepingId -profileIdCount)|| !appstate.sweepAnimation)  )
						boundaryProfiles.push_back( p3ds_1.back() ) ;
					if( sid != 0 && ( p3ds_2.size() < (appstate.sweepingId -profileIdCount)  ||  appstate.sweepAnimation)  )
						boundaryProfiles.push_back( p3ds_2.back() ) ;


					GlobalFun::draw3dCurves( p3ds_1, GLColor(0.1,0.1,0.1,1), true, false, 1, 1 ) ;
					GlobalFun::draw3dCurves( p3ds_2, GLColor(0.1,0.1,0.1,1), true, false, 1, 1 ) ;


					profileIdCount += std::max( p3ds_1.size(), p3ds_2.size() ) ;


					if( sid!= 0 )
						appstate.sweepingAnimation_meshNum += 1 ;
					if( sid == dataMgr.swp.settings.size()-1 && profileIdCount < appstate.sweepingId )
						appstate.sweepingAnimation_meshNum += 1 ;


					if( appstate.sweepAnimation && profileIdCount >= appstate.sweepingId )
						break ;
				}


				// draw boundary profiles

				if( boundaryProfiles.size() )
					draw3dCurves(boundaryProfiles, GLColor(0/255.0, 0/255.0, 255/255.0), true, false, 1.5 ) ;

			}
		}



		glCullFace(GL_FRONT) ;

		if( appstate.cullFace )
			glEnable(GL_CULL_FACE) ;
		else
			glDisable(GL_CULL_FACE) ;


	

		bindShader("Oren-Nayar") ;


		if( appstate.displayMesh  && !appstate.displayMeshVertOnly){
		
			if( appstate.sweepAnimation && appstate.sweepingAnimation_meshNum <dataMgr.swp.meshes.size()){
				for( int i=0; i<dataMgr.swp.meshes.size(); ++i ){

					if( i< appstate.sweepingAnimation_meshNum ){

						dataMgr.swp.meshes[i].sortDepth() ;
						dataMgr.swp.meshes[i].drawItself() ;
					}
				}
			}
			else{

					dataMgr.swp.mergedMesh.drawItself() ;


			}
		}

		else if( appstate.displayMesh  && appstate.displayMeshVertOnly)
			glDrawer.drawPolygonOnScreen( &dataMgr.swp.mergedMesh.cmeshFormat, global_paraMgr.drawer.getDouble("Original Dot Size")*0.02  * (takeSnapTile?0.25:1), (takeSnapTile?24:6 )) ;

		releaseShader();
		glCullFace(GL_BACK) ;





		if( get_operating_state() == onesweepUpdateResult ){
			dataMgr.swp.drawCtrlNurbs_cutplane() ;
		}
		if( get_operating_state() == onesweepReady ){
			dataMgr.swp.ongoingSkel.drawActiveTmp_cutplane() ;
		}

		//drawPoints( globalPointsToDraw1, GLColor(0,1,0 ) ) ;

		if( dataMgr.curves3d.size()!=0 )
			draw3dCurves(dataMgr.curves3d, GLColor(0.0, 0.5, 0.5) ) ;

		if( dataMgr.snappedCurves3d.size()!=0 )
			draw3dCurves(dataMgr.snappedCurves3d, GLColor(0.5, 0.5, 0.0) ) ;

		if( appstate.displayControlPoints ){
			////drawPoints( globalPointsToDraw2 ,  GLColor(0,0,1 ) ) ;
			////drawPoints( globalPointsToDraw1, GLColor(0,1,0 ) ) ;

			//for( int i=0;i<globalVectorsToDraw.size(); ++i ){
			//	Point3f p1 = globalVectorsToDraw[i].first ;
			//	Point3f p2 = globalVectorsToDraw[i].second ;
			//	glBegin(GL_LINES) ;
			//	glVertex3f(p1[0], p1[1], p1[2] ) ;
			//	glVertex3f(p2[0], p2[1], p2[2] ) ;
			//	glEnd();
			//}

		}





		//if( get_operating_state() == editResult ){
		//	dataMgr.skel.drawCtrlTmp() ;
		//	for( int i=0; i<dataMgr.skel.resultMesh.size(); ++i )
		//		dataMgr.skel.resultMesh[i].drawItself() ;
		//}
		if( curveForRefineSgement.size() ){
		
			if(isLeftPressed  )
				GlobalFun::draw2dCurveOnScreen( curveForRefineSgement, double3(0,1,1) ,4,false) ;
			else if( isRightPressed )
				GlobalFun::draw2dCurveOnScreen( curveForRefineSgement, double3(0,0.5,0.5) ,4,false) ;

		}

		/*********************************added begin*****************************************/
		if(!zcOp.mouseTrack.empty())
			draw2dCurves(zcOp.mouseTrack);

		//draw3dMapTo2d();

		/*********************************added end*****************************************/

		//drawAxises() ;




		glPopMatrix();
		glPopMatrix();

		if (takeSnapTile){
			cout << "snap shot!" << endl;
			pasteTile();

			double normal_width = global_paraMgr.drawer.getDouble("Normal Line Width");
			double dot_size = global_paraMgr.drawer.getDouble("Sample Dot Size");
			double original_dot_size = global_paraMgr.drawer.getDouble("Original Dot Size");

			global_paraMgr.drawer.setValue("Normal Line Width", DoubleValue(normal_width / (SnapResolutionScale * SnapResolutionScale * snapDrawScal)));
			global_paraMgr.drawer.setValue("Sample Dot Size", DoubleValue(dot_size / (SnapResolutionScale * SnapResolutionScale * snapDrawScal)));
			global_paraMgr.drawer.setValue("Original Dot Size", DoubleValue(original_dot_size / (SnapResolutionScale * SnapResolutionScale * snapDrawScal)));
		}


		if( dataMgr.curves2d.size()!=0 )
			draw2dCurves() ;



	}

PAINT_RETURN:
	paintMutex.unlock();
}
