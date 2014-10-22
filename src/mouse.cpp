
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

extern double glareaFov  ;

extern std::vector<Curve2D > greenCurves  ;  
extern std::vector<Curve2D > blueCurves  ;  


void GLArea::wheelEvent(QWheelEvent *e) 
{
	const int WHEEL_STEP = 120 * 5;
	double change_rate = 0.1;
	double change = (e->delta() < 0) ? (1 + change_rate) : (1 - change_rate);
	double size_temp = 0.0;


	if( (e->modifiers() & Qt::AltModifier) && (e->modifiers() & Qt::ControlModifier) )
	{
		size_temp = global_paraMgr.drawer.getDouble("Normal Line Length");
		global_paraMgr.drawer.setValue("Normal Line Length", DoubleValue(size_temp * change));
	}
	else if( (e->modifiers() & Qt::ShiftModifier) && (e->modifiers() & Qt::ControlModifier) )
	{
	}
	else if((e->modifiers() & Qt::ShiftModifier) && (e->modifiers() & Qt::AltModifier))
	{
	}
	else
	{
		switch(e->modifiers())
		{
		case Qt::ControlModifier:

			if(para->getBool("Show Samples") &&  para->getBool("Show Samples Dot")
				&&para->getBool("Show Original") && para->getBool("Show Original Dot") )
			{
				size_temp = global_paraMgr.drawer.getDouble("Sample Dot Size") * change;
				if(size_temp < 1)
				{
					size_temp = 1;
				}
				global_paraMgr.drawer.setValue("Sample Dot Size", DoubleValue(size_temp));

				size_temp = global_paraMgr.drawer.getDouble("Original Dot Size") * change;
				if(size_temp < 1)
				{
					size_temp = 1;
				}

				global_paraMgr.drawer.setValue("Original Dot Size", DoubleValue(size_temp));
			}
			else if(para->getBool("Show Samples") &&  para->getBool("Show Samples Dot") )
			{
				size_temp = global_paraMgr.drawer.getDouble("Sample Dot Size") * change;
				if(size_temp < 1)
				{
					size_temp = 1;
				}
				global_paraMgr.drawer.setValue("Sample Dot Size", DoubleValue(size_temp));
			}
			else if(para->getBool("Show Original") && para->getBool("Show Original Dot") )
			{
				size_temp = global_paraMgr.drawer.getDouble("Original Dot Size") * change;
				if(size_temp < 1)
				{
					size_temp = 1;
				}
				global_paraMgr.drawer.setValue("Original Dot Size", DoubleValue(size_temp));
			}
			else
			{
				size_temp = global_paraMgr.drawer.getDouble("Sample Draw Width");
				global_paraMgr.drawer.setValue("Sample Draw Width", DoubleValue(size_temp * change));

				size_temp = global_paraMgr.drawer.getDouble("Original Draw Width");
				global_paraMgr.drawer.setValue("Original Draw Width", DoubleValue(size_temp * change));

				size_temp = global_paraMgr.drawer.getDouble("Quade Draw Width");
				global_paraMgr.drawer.setValue("Quade Draw Width", DoubleValue(size_temp * change));


			}
			emit needUpdateStatus();
			break;

		case Qt::ShiftModifier:
			size_temp = global_paraMgr.data.getDouble("Down Sample Num");
			global_paraMgr.setGlobalParameter("Down Sample Num", DoubleValue(size_temp * change));
			emit needUpdateStatus();
			break;

		case  Qt::AltModifier:
			size_temp = global_paraMgr.data.getDouble("CGrid Radius");
			global_paraMgr.setGlobalParameter("CGrid Radius", DoubleValue(size_temp * change));
			initSetting();
			break;

		default:
			trackball.MouseWheel( e->delta()/ float(WHEEL_STEP));
			break;
		}
	}

	updateGL();
}





void GLArea::mouseDoubleClickEvent(QMouseEvent *e){

	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);
	x2 = e->x();
	y2 = viewport[3] - e->y();
	Point2f mpos = Point2f(x2, y2) ;
	
	if( get_operating_state() == onesweepReady ){
		dataMgr.swp.ongoingSkel.selectActiveTmp(mpos) ; 
	}

	if( get_operating_state() == onesweepUpdateResult )
		dataMgr.swp.ctrlSkelId = dataMgr.swp.ctrlProfId = -1;

}

void GLArea::mouseMoveEvent(QMouseEvent *e)
{


	Qt::KeyboardModifiers keyMod = QGuiApplication::keyboardModifiers();
	bool ctrlPressed = keyMod.testFlag(Qt::ControlModifier);
	bool altPressed = keyMod.testFlag(Qt::AltModifier);
	bool shiftPressed = keyMod.testFlag(Qt::ShiftModifier) ;

	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);

	x2 = e->x();
	y2 = viewport[3] - e->y();


	if( get_operating_state() == onesweep && altPressed  ){
		if( isLeftPressed )
			greenCurves.back().push_back( Point2f(x2, y2) ) ;//draw colorful line for fun
		else if( isRightPressed )
			blueCurves.back().push_back( Point2f(x2, y2) ) ;//draw colorful line for fun
		return ;
	}

	//onesweep state
	if( get_operating_state() == onesweep ){


		if( isRightPressed && ctrlPressed ){
			if( (Point2f(x2,y2) - dataMgr.swp.ongoingSkelStroke.back()).Norm() > 5 )
				dataMgr.swp.moveSkelStroke( Point2f(x2,y2) ) ;//select skeleton
			updateGL() ;
			return ;
		}

		if( isRightPressed && shiftPressed ){
			if( (Point2f(x2,y2) - dataMgr.swp.ongoingFeatureStroke.back()).Norm() > 5 )
				dataMgr.swp.moveFeatureStroke( Point2f(x2,y2) ) ;
			updateGL() ;
			return ;
		}

	}

	//onesweepReady state,move control points
	if( get_operating_state() == onesweepReady ){

		if( dataMgr.swp.ongoingSkel.activeTmpId >= 0 &&  dataMgr.swp.ongoingSkel.activeCVID >= 0 && isLeftPressed ){
			dataMgr.swp.ongoingSkel.moveActiveTmpCV(Point2f(x2,y2)) ; 
			updateGL() ;
			return;
		}
	}
	//onesweepSegment 
	if( get_operating_state() == onesweepSegment ){		

		if( isLeftPressed || isRightPressed ){
			Point2f currPos = Point2f(e->x(), viewport[3] - e->y() )  ;
			if(curveForRefineSgement.size()==0 || (curveForRefineSgement.back()-currPos).Norm() > 15 )
				curveForRefineSgement.push_back( currPos ) ;
			updateGL() ;
			return ;

		}
	}
	//U
	if( get_operating_state() == onesweepUpdateResult ){

		if( dataMgr.swp.ctrlSkelId >= 0 &&  dataMgr.swp.ctrlProfId >= 0 &&  dataMgr.swp.ctrlCVId >= 0 && isLeftPressed ){
			dataMgr.swp.moveCtrlProfCV(Point2f(x2,y2)) ; 
			updateGL() ;
			return;
		}

		if( ctrlPressed && isRightPressed && dataMgr.swp.editResultStroke.size() ){
			dataMgr.swp.moveEditResultStroke( Point2f(x2,y2) ) ;
			updateGL() ;
			return ;
		}



	}


	
	if( get_operating_state() == profEditMode && isLeftPressed && SelectedNurbsId>=0 && SelectedCtrPointId >= 0){

		int ndeg = ReconstructorPara::NurbsDegree  ;
		//暂时
		std::vector<Point2f> &ctrpts = dataMgr.swp.ongoingSkel.NurbsCVPixels[SelectedNurbsId] ; 
		ctrpts[SelectedCtrPointId] = Point2f(x2, y2) ;
		if( SelectedCtrPointId < ndeg ) ctrpts[SelectedCtrPointId+ctrpts.size()-ndeg] = Point2f(x2, y2) ;
		if( SelectedCtrPointId >= ctrpts.size()-ndeg ) ctrpts[SelectedCtrPointId-(ctrpts.size()-ndeg)] = Point2f(x2, y2) ;

		dataMgr.swp.ongoingSkel.cvtPixels2CV();
		dataMgr.swp.ongoingSkel.discretizeNurbs();
		dataMgr.swp.ongoingSkel.convertDisNurbs2Pixles();
		updateGL() ;
		return;
	}


	/*if( get_operating_state() == editResult && isLeftPressed && SelectedCtrPointId >= 0 && dataMgr.skel.ctrTmpValid ){

	dataMgr.skel.moveCtrTmpCV(SelectedCtrPointId, Point2f(x2, y2) ) ;
	updateGL() ;
	return;
	}*/

	//??
	if( get_operating_state() == segmentMode && isLeftPressed ){
		Point2f currPos = Point2f(e->x(), viewport[3] - e->y() )  ;
		if(curveForRefineSgement.size()==0 || (curveForRefineSgement.back()-currPos).Norm() > 15 )
			curveForRefineSgement.push_back( currPos ) ;
		updateGL() ;
		return ;
	}




	/*********************************added begin*****************************************/
	//E
	if( get_operating_state() == skelEditState && isLeftPressed ){
		Point2f currPos = Point2f(e->x(), viewport[3] - e->y() )  ;
		if(zcOp.mouseTrack.size()==0 || (zcOp.mouseTrack.back()-currPos).Norm() > 10 )
			zcOp.mouseTrack.push_back( currPos ) ;

		updateGL() ;
		return ;
	}
	if( get_operating_state() == skelEditState && isRightPressed ){
		Point2f currPos = Point2f(e->x(), viewport[3] - e->y() )  ;
		zcOp.mouseTrack.resize(1);
		zcOp.mouseTrack.push_back(currPos);

		updateGL() ;
		return ;
	}
	/*********************************added end*****************************************/

	// if the state is freeLassoState, then store the position of cursor into 2d curves
	//L
	if( get_operating_state() == freeLassoState && isLeftPressed ){
		Point2f currPos = Point2f(e->x(), viewport[3] - e->y() )  ;
		if(dataMgr.curves2d.back().size()==0 || (dataMgr.curves2d.back().back()-currPos).Norm() > 5 )
			dataMgr.curves2d.back().push_back( currPos ) ;

		updateGL() ;
		return ;
	}


	if (isDefaultTrackBall())
	{
		trackball.MouseMove(e->x(),height()-e->y());
	}
	else 
	{
		trackball_light.MouseMove(e->x(),height()-e->y());
	}

	updateGL();
}

vcg::Trackball::Button QT2VCG(Qt::MouseButton qtbt,  Qt::KeyboardModifiers modifiers)
{
	int vcgbt= vcg::Trackball::BUTTON_NONE;
	if(qtbt & Qt::LeftButton		) vcgbt |= vcg::Trackball::BUTTON_LEFT;
	if(qtbt & Qt::RightButton		) vcgbt |= vcg::Trackball::BUTTON_RIGHT;
	if(qtbt & Qt::MidButton			) vcgbt |= vcg::Trackball::BUTTON_MIDDLE;
	if(modifiers & Qt::ShiftModifier		)	vcgbt |= vcg::Trackball::KEY_SHIFT;
	if(modifiers & Qt::ControlModifier ) vcgbt |= vcg::Trackball::KEY_CTRL;
	if(modifiers & Qt::AltModifier     ) vcgbt |= vcg::Trackball::KEY_ALT;
	return vcg::Trackball::Button(vcgbt);
}

void GLArea::mousePressEvent(QMouseEvent *e)
{
	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);

	x2 = e->x();
	y2 = viewport[3] - e->y();


	Qt::KeyboardModifiers keyMod = QGuiApplication::keyboardModifiers();
	bool ctrlPressed = keyMod.testFlag(Qt::ControlModifier);
	bool altPressed = keyMod.testFlag(Qt::AltModifier);
	bool shiftPressed = keyMod.testFlag(Qt::ShiftModifier) ;

	//std::cout << "mousePress" <<std::endl;



	if(e->button() == Qt::LeftButton){

		isLeftPressed = true ;

		//for fun
		if( get_operating_state() == onesweep && altPressed ){
			greenCurves.resize( greenCurves.size() + 1 ) ;
			return ;
		}

		if( get_operating_state() == onesweepReady ){

			if( !ctrlPressed ){
				// select control points

				Point2f mpos = Point2f(x2, y2) ;

				dataMgr.swp.ongoingSkel.selectActiveTmpCV(mpos) ; 

			}else{ 

				std::cout << " set sharp feature" <<std::endl;
				// set sharp feature
				Point2f mpos = Point2f(x2, y2) ;
				int nearstNurbsId = 0 ;
				int nearstCVId = 0 ;
				double mindis = 1e10 ;
				for( int i=0; i<dataMgr.swp.ongoingSkel.NurbsCVPixels.size(); ++i ){
					for( int j=0; j<dataMgr.swp.ongoingSkel.NurbsCVPixels[i].size(); ++j ){
						double dis = (mpos - dataMgr.swp.ongoingSkel.NurbsCVPixels[i][j]).Norm()  ;
						if( dis < mindis ){
							mindis = dis ;
							nearstNurbsId = i ;
							nearstCVId = j ;
						}
					}
				}

				if( mindis < 20 ){
					dataMgr.swp.ongoingSkel.setSharpNurbsCV(nearstNurbsId, nearstCVId) ;

					updateGL() ;
					return;
				}

			}
		}

		if( get_operating_state() == onesweepUpdateResult ){

			dataMgr.swp.selectCtrlProfCV( Point2f(x2, y2)) ;

		}

		if( get_operating_state() == segmentMode  ){
			std::cout<<"segmentMode,已被我注释掉！！！"<<endl;
			//暂时
			if( altPressed ){
				dataMgr.swp.ongoingSkel.addTipSeed(  Point2f(e->x(), viewport[3] - e->y() ) ) ;
				return ;
			}

			if( shiftPressed ){
				dataMgr.swp.ongoingSkel.removeTips(Point2f(e->x(), viewport[3] - e->y() ) ) ;
				return ;
			}

			curveForRefineSgement.clear() ;
			GLint viewport[4];
			glGetIntegerv (GL_VIEWPORT, viewport);
			curveForRefineSgement.push_back( Point2f(e->x(), viewport[3] - e->y() ) ) ;

			return ;
		}



		if( get_operating_state() == onesweepSegment ){			

			curveForRefineSgement.clear() ;
			GLint viewport[4];
			glGetIntegerv (GL_VIEWPORT, viewport);
			curveForRefineSgement.push_back( Point2f(e->x(), viewport[3] - e->y() ) ) ;
			return ;
		}

		/*********************************added begin*****************************************/
		zcOp.mouseTrack.clear();
		if(get_operating_state()==skelEditState)
		{
			GLint viewport[4];
			glGetIntegerv(GL_VIEWPORT,viewport);
			zcOp.mouseTrack.push_back( Point2f(e->x(), viewport[3] - e->y() ) );
			
			updateGL() ;
			return;
		}
		/*********************************added end*****************************************/


		dataMgr.curves2d.clear();
		if( get_operating_state() == freeLassoState ){
			dataMgr.curves2d.resize( 1 ) ;
			GLint viewport[4];
			glGetIntegerv (GL_VIEWPORT, viewport);
			dataMgr.curves2d.back().push_back( Point2f(e->x(), viewport[3] - e->y() ) ) ;

			return ;
		}

	}
	/*********************************added begin*****************************************/
	if(e->button()==Qt::RightButton)
	{

		isRightPressed = true ;

		//for fun
		if( get_operating_state() == onesweep && altPressed ){
			blueCurves.resize( blueCurves.size() + 1 ) ;
			return ;
		}

		if( get_operating_state() == onesweep || get_operating_state() == onesweepReady ){
			if( ctrlPressed ){
				set_operating_state(onesweep) ;
				dataMgr.swp.startSkelStroke( Point2f(x2,y2) ) ;
				return ;
			}

			if( shiftPressed ){
				set_operating_state(onesweep) ;
				dataMgr.swp.startFeatureStroke( Point2f(x2,y2) ) ;
				return ;
			}
		}

		if( get_operating_state() == onesweepUpdateResult ){

			if( ctrlPressed ){
				dataMgr.swp.startEditResultStroke( Point2f(x2, y2)) ;
			}
		}

		zcOp.mouseTrack.clear();
		if(get_operating_state()==skelEditState)
		{
			isRightPressed=true;
			GLint viewport[4];
			glGetIntegerv(GL_VIEWPORT,viewport);
			zcOp.mouseTrack.push_back( Point2f(e->x(), viewport[3] - e->y() ) );

			updateGL() ;
			return;
		}


		if( get_operating_state() == segmentMode  ){

			//暂时
			if( ctrlPressed ){
				dataMgr.swp.ongoingSkel.switchTemplate(  Point2f(e->x(), viewport[3] - e->y() ) ) ;
				return ;
			}
		}

		if( get_operating_state() == onesweepReady  ){
			dataMgr.swp.ongoingSkel.switchTemplate_online(  Point2f(e->x(), viewport[3] - e->y() ) ) ;
			return ;
		}

		if( get_operating_state() == onesweepSegment ){			
			curveForRefineSgement.clear() ;
			GLint viewport[4];
			glGetIntegerv (GL_VIEWPORT, viewport);
			curveForRefineSgement.push_back( Point2f(e->x(), viewport[3] - e->y() ) ) ;
			return ;
		}

	}
	/*********************************added end*****************************************/

	if ((e->modifiers() & Qt::ShiftModifier) && (e->modifiers() & Qt::ControlModifier) &&
		(e->button()==Qt::LeftButton) )
		activeDefaultTrackball=false;
	else activeDefaultTrackball=true;

	if (isDefaultTrackBall())
	{
		if(e->button() == Qt::LeftButton)
			trackball.MouseDown(e->x(), height() - e->y(), QT2VCG(e->button(), e->modifiers() ) );     

	}
	else trackball_light.MouseDown(e->x(),height()-e->y(), QT2VCG(e->button(), Qt::NoModifier ) );

	//isDragging = true;
	update();
	updateGL();
}
void GLArea::mouseReleaseEvent(QMouseEvent *e)
{
	Qt::KeyboardModifiers keyMod = QGuiApplication ::keyboardModifiers ();
	bool ctrlPressed = keyMod.testFlag(Qt::ControlModifier);
	bool altPressed = keyMod.testFlag(Qt::AltModifier);
	bool shiftPressed = keyMod.testFlag(Qt::ShiftModifier) ;

	//std::cout << "mouseRelease" <<std::endl;

	if( get_operating_state() == onesweep && altPressed  ){
		if( e->button() == Qt::LeftButton )
			if( greenCurves.back().size() < 5 )  greenCurves.resize( greenCurves.size()-1 ) ;
			else if( e->button() == Qt::RightButton )
				if( blueCurves.back().size() < 5 )  blueCurves.resize( blueCurves.size()-1 ) ;

		isLeftPressed = isRightPressed = false ;

		return ;
	}


	if(e->button() == Qt::LeftButton){

		isLeftPressed = false ;
		if( get_operating_state() == profEditMode || get_operating_state() == editResult || get_operating_state() == onesweepReady ){


			SelectedNurbsId = -1 ;
			SelectedCtrPointId = -1 ;


			dataMgr.swp.ongoingSkel.deselectActiveTmpCV() ;

			updateGL() ;
			//return;
		}


		if( get_operating_state() == onesweepUpdateResult  ){

			dataMgr.swp.ctrlCVId = -1 ;
		}


		if( get_operating_state() == segmentMode  ){

			if( altPressed || shiftPressed)
				return ;
			//暂时
			if( curveForRefineSgement.size() > 2 )
			dataMgr.swp.ongoingSkel.refineSegment( curveForRefineSgement,ctrlPressed ) ;  //,ctrlPressed很奇怪啊，为什么没return
			curveForRefineSgement.clear() ;
			updateGL() ;
			return ;
		}



		if( get_operating_state() == onesweepSegment ){	
			if( curveForRefineSgement.size() > 3 )
				dataMgr.swp.segment( curveForRefineSgement, shiftPressed, ctrlPressed, true ) ;
			curveForRefineSgement.clear() ;
			updateGL() ;
			return ;
		}
		/*********************************added begin*****************************************/

		if(get_operating_state()==skelEditState)
		{
			zcOp.doOp_Left(dataMgr.skel);
			zcOp.mouseTrack.clear();
			updateGL();
			return;
		}
		zcOp.mouseTrack.clear();
		/*********************************added end*****************************************/
	}

	/*********************************added begin*****************************************/
	if(e->button()==Qt::RightButton)
	{
		isRightPressed = false ;

		if( get_operating_state() == onesweep  ){

			if( ctrlPressed ){
				float mvmatrix[16];  for( int i=0; i<16;++i ) mvmatrix[i]=globalModelviewMatrix[i] ;
				bool settingBuilt = dataMgr.swp.endSkelStroke( Point2f(x2,y2),mvmatrix ) ;
				if( settingBuilt )
					set_operating_state(onesweepReady) ;
				updateGL() ;
			}
			if( shiftPressed ){
				float mvmatrix[16];  for( int i=0; i<16;++i ) mvmatrix[i]=globalModelviewMatrix[i] ;
				dataMgr.swp.endFeatureStroke( Point2f(x2,y2),mvmatrix ) ;
				updateGL() ;
			}

		}
		if( get_operating_state() == onesweepUpdateResult  ){

			if( ctrlPressed ){
				dataMgr.swp.EndEditResultStroke() ;
				updateGL() ;
			}
		}

		if(get_operating_state()==skelEditState) {
			zcOp.doOp_Right(dataMgr.skel);
			zcOp.mouseTrack.clear();
			updateGL();
			return;
		}
		zcOp.mouseTrack.clear();


		if( get_operating_state() == onesweepSegment ){	
			if( curveForRefineSgement.size() > 3 )
				dataMgr.swp.segment( curveForRefineSgement, shiftPressed, ctrlPressed,false ) ;
			curveForRefineSgement.clear() ;
			updateGL() ;
			return ;
		}

	}
	zcOp.mouseTrack.clear();
	/*********************************added end*****************************************/

	if(e->button() == Qt::LeftButton)
		trackball.MouseUp(e->x(),height()-e->y(), QT2VCG(e->button(), e->modifiers() ) );

	updateGL();
}




void GLArea::keyReleaseEvent ( QKeyEvent * e ) {

	Qt::KeyboardModifiers keyMod = QGuiApplication ::keyboardModifiers ();
	bool ctrlPressed = keyMod.testFlag(Qt::ControlModifier);
	bool altPressed = keyMod.testFlag(Qt::AltModifier);
	bool shiftPressed = keyMod.testFlag(Qt::ShiftModifier) ;


	if( e->key() == Qt::Key_I  ){


		if(dataMgr.swp.meshes.size()==0 ){

			// initialize sweeper
			std::vector<Point3f> pointset ;
			for( int i=0; i<dataMgr.original.vert.size(); ++i )
				pointset.push_back(dataMgr.original.vert[i].P() ) ;
			std::vector<Point3f> normals ;
			for( int i=0; i<dataMgr.original.vert.size(); ++i )
				normals.push_back(dataMgr.original.vert[i].N() ) ;

			dataMgr.skel = skeleton_mul(dataMgr.skel.smoothedBrchPts, true ) ;

			std::vector<Point3f> allJoints;
			std::vector<std::vector<int>> ajBrchId ;
			dataMgr.skel.getJoints(allJoints,ajBrchId);
			cout<<"points number"<<allJoints.size()<<endl;
			dataMgr.swp = sweeper(pointset,normals, dataMgr.skel.smoothedBrchPts,allJoints,dataMgr.theBranches ) ;
		}


		set_operating_state(onesweep);
	}

	if( e->key() == Qt::Key_S  ){
		if( get_operating_state() == onesweepReady )
			set_operating_state(onesweepSegment) ;
		else if( get_operating_state() == onesweepSegment )
			set_operating_state(onesweepReady) ;
	}


	if( e->key() == Qt::Key_R ){
		if( get_operating_state() == onesweepReady ){
			// solve
			
			dataMgr.swp.reconLastStroke() ;
		}

		if( get_operating_state() == onesweepUpdateResult )
			if( dataMgr.swp.ctrlSkelId>=0 && dataMgr.swp.ctrlProfId>=0)
				dataMgr.swp.updateResult() ;

	}

	if( e->key() == Qt::Key_U ){
		if( get_operating_state() == onesweepReady || get_operating_state()==onesweep){
			set_operating_state(onesweepUpdateResult) ;
			dataMgr.swp.ptsSelected = std::vector<bool>(dataMgr.swp.points.size(), false) ;
		}
	}

	


	if( e->key() == Qt::Key_M ){   // store modelview
		std::ofstream ofs(appstate.path + "\/mvmatrix.txt") ;
		for( int i=0; i<16;++i ) ofs << globalModelviewMatrix[i] << std::endl;
		ofs.close() ;
	}
	if( e->key() == Qt::Key_L ){  // recover modelview

		if( !mvmLoaded ){

			std::ifstream ifs(appstate.path + "\/mvmatrix.txt") ;
			for( int i=0; i<16;++i ) ifs >> mvmatrix_loaded[i] ;
			if( !ifs.fail() )
				mvmLoaded = true ;
			ifs.close() ;
		}else
			mvmLoaded = false ;
	}


	//if( e->key() == Qt::Key_Return ){  
	//	std::cout << "input cmd please:"<<std::endl;
	//	std::string cmd ; cin >> cmd ;
	//	runCommand(cmd) ;
	//}



	if( e->key() == Qt::Key_Escape  ){

		if( get_operating_state() == onesweepReady ){
			set_operating_state(onesweep) ;
			dataMgr.swp.ptsSelected = std::vector<bool>(dataMgr.swp.points.size(), false) ;
		}else if( get_operating_state() == onesweep ){
			if( dataMgr.swp.meshes.size() ){
				dataMgr.swp.meshes.erase( dataMgr.swp.meshes.begin() + dataMgr.swp.meshes.size() - 1 ) ;
				dataMgr.swp.meshLable.erase( dataMgr.swp.meshLable.begin() + dataMgr.swp.meshLable.size() - 1 ) ;
				dataMgr.swp.settings.erase( dataMgr.swp.settings.begin() + dataMgr.swp.settings.size() - 1 ) ;

				dataMgr.swp.mergeMesh() ;

			}
		}
	}


	if( e->key() == Qt::Key_F1 ){
		initShader() ;
	}

	extern 	double planeWidth  ;
	if( e->key() == Qt::Key_Left ){
		planeWidth -= 0.1 ;
		if( planeWidth <0.1 )
			planeWidth = 0.05 ;
		if( planeWidth <0.05 )
			planeWidth = 0.025 ;


		std::cout << "planeWidth = " << planeWidth <<std::endl; 
	}

	if( e->key() == Qt::Key_Right ){
		planeWidth += 0.1 ;
		if( planeWidth > 1 )
			planeWidth = 1 ;
		std::cout << "planeWidth = " << planeWidth <<std::endl; 
	}




	if(e->key()==Qt::Key_Control) trackball.MouseUp(0,0, QT2VCG(Qt::NoButton, Qt::ControlModifier ) );
	if(e->key()==Qt::Key_Shift) trackball.MouseUp(0,0, QT2VCG(Qt::NoButton, Qt::ShiftModifier ) );
	if(e->key()==Qt::Key_Alt) trackball.MouseUp(0,0, QT2VCG(Qt::NoButton, Qt::AltModifier ) );



	updateGL() ;




}



void GLArea::contextMenuEvent ( QContextMenuEvent * e )
{
	//QMenu *menu = new QMenu(this);

	//QFont font;
	//font.setPointSize(14);
	//menu->setFont(font);
	////menu->addAction(lassoAction);
	//for( int i=0; i<contextActions.size(); ++i )
	//	menu->addAction( contextActions[i] ) ;

	//menu->exec(QCursor::pos());


}
