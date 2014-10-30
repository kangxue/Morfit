
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

double glareaFov  ;
std::vector<Curve2D > greenCurves  ;  ;
std::vector<Curve2D > blueCurves  ;  ;

GLArea::GLArea(QWidget *parent): QGLWidget(/*QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer |QGL::SampleBuffers),*/ parent),
	                             para(global_paraMgr.getGlareaParameterSet()),
								 glDrawer(global_paraMgr.getDrawerParameterSet()),
								 dataMgr(global_paraMgr.getDataParameterSet()),
								 paintMutex(QMutex::NonRecursive)
{
	setMouseTracking(true); 
	isRightPressed = false;

	trackball_light.center=Point3f(0, 0, 0);
	trackball_light.radius= 1;
	activeDefaultTrackball=true;

	glareaFov = fov = 60;
	clipRatioFar = 1;
	clipRatioNear = 1;
	nearPlane = .2f;
	farPlane = 5.f;
	takeSnapTile = false;
	default_snap_path = QString( ".\\snapshot\\");

	current_snap_path = default_snap_path;

	QDir dir ;
	dir.mkpath ( dir.currentPath() + QString("\/snapshot\/") );

	std::cout <<"shot dir = " << (dir.currentPath() + QString("\/snapshot\/")).toStdString()<<std::endl;

	snapDrawScal = 1;
	is_paintGL_locked = false;

	cout << "GLArea constructed" << endl;

	CVertex v;
	cout << sizeof(v) << endl;

	state = startState ;

	isLeftPressed = false ;
	
	SelectedNurbsId = -1 ;
	SelectedCtrPointId = -1 ;

	mvmLoaded = false ;
}

GLArea::~GLArea(void)
{
}


void GLArea::initializeGL()
{
	//initializing the GL
	cout << "initializeGL" << endl;

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0); 
	glEnable(GL_LIGHT1);
	//glEnable(GL_COLOR_MATERIAL);  

	glShadeModel (GL_SMOOTH);
	glShadeModel(GL_FLAT);

	//glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST ); 
	GLfloat mat_ambient[4] = {0.1745, 0.01175, 0.01175,1.0}; 
	GLfloat mat_diffuse[4] = {0.61424, 0.04136, 0.04136, 1.0 };
	GLfloat mat_specular[] = {0.727811, 0.626959, 0.626959, 1.0 };
	GLfloat shininess = 0.6*128;

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_ambient);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse); 
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, shininess);

	glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT );
	glColorMaterial( GL_FRONT_AND_BACK, GL_DIFFUSE );
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	glEnable(GL_NORMALIZE);

	glDisable(GL_CULL_FACE);
	glColor4f(1, 1, 1, 1);   

	glEnable(GL_LIGHTING);

	initLight();


	trackball.center = Point3f(0, 0, 0);
	trackball.radius = 1;

	// force to open anti
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_POLYGON_SMOOTH);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	glLoadIdentity();
	initShader() ;
	glClearColor(1,1,1,0) ;

}

void GLArea::initLight()
{
	GLColor color_ambient =  QColor(150, 150, 150) ;
	GLColor color_diffuse =  QColor(241, 241, 241) ;
	GLColor color_specular = QColor(255, 255, 255) ;
	float ambient[4] = {color_ambient.r, color_ambient.g, color_ambient.b, 1.0};
	float diffuse[4] = {color_diffuse.r, color_diffuse.g, color_diffuse.b, 1.0};
	float specular[4] = {color_specular.r, color_specular.g, color_specular.b, 1.0};

	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT2, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT3, GL_AMBIENT, ambient);

	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT3, GL_DIFFUSE, diffuse);

	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT2, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT3, GL_SPECULAR, specular);


	float position0[4] = {0,0,-1,0} ; 
	float position1[4] = {0,0,1,0} ; 
	float position2[4] = {1,0,0,0} ; 
	float position3[4] = {-1,0,0,0} ; 
	glLightfv(GL_LIGHT0, GL_POSITION, position0);
	glLightfv(GL_LIGHT1, GL_POSITION, position1);
	glLightfv(GL_LIGHT2, GL_POSITION, position2);
	glLightfv(GL_LIGHT3, GL_POSITION, position3);

}

void GLArea::resizeGL(int w, int h)
{
	//cout << "resizeGL" << endl;

	glViewport(0, 0, (GLint)w, (GLint)h);  

	//Applies subsequent matrix operations to the projection matrix stack.
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();  

	float r = w/(float)h;
	//set up a perspective projection matrix
	gluPerspective(fov, r, 0.1, 10);

	//Applies subsequent matrix operations to the modelview matrix stack.
	glMatrixMode(GL_MODELVIEW);  
}


void GLArea::lightOnOff(bool _val)
{
	if(_val)
	{
		glEnable(GL_LIGHTING);
	}
	else
	{
		glDisable(GL_LIGHTING);
	}
}

void GLArea::initAfterOpenFile()
{
	zcOp.init();
	initView();
	initSetting();
}

void GLArea::initSetting()
{
	emit needUpdateStatus();
}

void GLArea::initView()
{
	dataMgr.recomputeBox();

	if (!dataMgr.isOriginalEmpty())
	{
		gl_box = dataMgr.getCurrentOriginal()->bbox;
	}
	else if(!dataMgr.isSamplesEmpty())
	{
		gl_box = dataMgr.getCurrentSamples()->bbox;
	}
}

//void GLArea::runPointCloudAlgorithm(PointCloudAlgorithm& algorithm)
//{
//	//input the Algorithm name ,and the run the algorithm
//	paintMutex.lock();
//
//	QString name = algorithm.getParameterSet()->getString("Algorithm Name");
//	cout << "*********************************** Start  " << name.toStdString() << "  ***********************************" << endl;
//	int starttime, stoptime, timeused;
//	starttime = clock();
//
//	algorithm.setInput(&dataMgr);
//	algorithm.run();
//	algorithm.clear();
//
//	stoptime = clock();
//	timeused = stoptime - starttime;
//
//	int currentUsedTime = timeused/double(CLOCKS_PER_SEC);
//
//	cout << "time used:  " << timeused/double(CLOCKS_PER_SEC) << " seconds." << endl;
//	cout << "*********************************** End  " << name.toStdString() << "  ***********************************" << endl;
//	cout << endl << endl << endl;
//
//	paintMutex.unlock();
//}
//

void GLArea::openByDrop(QString fileName)
{
	//user can load a file by drop(.ply)
	if(fileName.endsWith("ply"))
	{
		if (fileName.contains("original"))
		{
			dataMgr.loadPlyToOriginal(fileName);
		}
		dataMgr.loadPlyToSample(fileName);
	}

	initAfterOpenFile();
	updateGL();
}

void GLArea::loadDefaultModel()
{
	dataMgr.loadPlyToSample("default.ply");
	dataMgr.loadPlyToOriginal("default_original.ply");
	initAfterOpenFile();
	updateGL();
}

void GLArea::drawPickRect()
{
	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0,width(),height(),0,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glColor3f(1,1,1);

	glBegin(GL_LINE_LOOP);
	glVertex2f(x1,viewport[3] - y1);
	glVertex2f(x2,viewport[3] - y1);
	glVertex2f(x2,viewport[3] - y2);
	glVertex2f(x1,viewport[3] - y2);
	glEnd();
	glDisable(GL_LOGIC_OP);

	// Closing 2D
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}


void GLArea::draw2dCurves()
{
	//glEnable(GL_BLEND) ;
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;

	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

	glDisable(GL_LIGHTING) ;

	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0,width(),height(),0,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glColor4f(1,0,1, 0.7);
	for( int i=0; i< dataMgr.curves2d.size(); ++ i ){
		int n = dataMgr.curves2d[i].size();
		for( int id = 0; id < dataMgr.curves2d[i].size(); ++id ){

			glLineWidth( 6.0 ) ;

			if( (id != dataMgr.curves2d[i].size() - 1) || ! isLeftPressed  ){
				glBegin(GL_LINES);
				glVertex2f( dataMgr.curves2d[i][id].X(),viewport[3] - dataMgr.curves2d[i][id].Y()  ) ;
				glVertex2f( dataMgr.curves2d[i][ (id+1)%n].X(),viewport[3] - dataMgr.curves2d[i][(id+1)%n].Y()  ) ;
				glEnd();
			}

			glPointSize( 8 / 2) ;
			glBegin(GL_POINTS);
			glVertex2f( dataMgr.curves2d[i][id].X(),viewport[3] - dataMgr.curves2d[i][id].Y()  ) ;
			glEnd();

		}

	}


	// Closing 2D
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);



	glEnable(GL_LIGHTING) ;
}
/*********************************added begin*****************************************/
void GLArea::draw2dCurves(vector<Point2f> track)
{
	//glEnable(GL_BLEND) ;
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;

	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);


	glDisable(GL_LIGHTING) ;

	GLint viewport[4];
	glGetIntegerv (GL_VIEWPORT, viewport);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0,width(),height(),0,-1,1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	if(isLeftPressed)
	{
		glColor4f(1,0,1, 0.7);
		int n = track.size();
		for( int id = 0; id < track.size(); ++id ){

			glLineWidth( 6.0 ) ;

			if( (id != track.size() - 1) && isLeftPressed  ){
				glBegin(GL_LINES);
				glVertex2f( track[id].X(),viewport[3] - track[id].Y()  ) ;
				glVertex2f( track[ (id+1)%n].X(),viewport[3] - track[(id+1)%n].Y()  ) ;
				glEnd();
			}

			glPointSize( 8 / 2) ;
			glBegin(GL_POINTS);
			glVertex2f( track[id].X(),viewport[3] - track[id].Y()  ) ;
			glEnd();

		}
	}
	else if(isRightPressed)
	{
		glColor4f(0,1,0.7,1);
		int n = track.size();
		for( int id = 0; id < track.size(); ++id ){

			glLineWidth( 3.0 ) ;

			if( (id != track.size() - 1) && isRightPressed  ){
				glBegin(GL_LINES);
				glVertex2f( track[id].X(),viewport[3] - track[id].Y()  ) ;
				glVertex2f( track[ (id+1)%n].X(),viewport[3] - track[(id+1)%n].Y()  ) ;
				glEnd();
			}

			glPointSize( 8 / 2) ;
			glBegin(GL_POINTS);
			glVertex2f( track[id].X(),viewport[3] - track[id].Y()  ) ;
			glEnd();

		}

		if(track.size()==2)
		{
			Point2f v=(track[1]-track[0]).Normalize();
			GLfloat arrowLenth=20.0f;

			glPushMatrix();
			glTranslatef(track[1].X(),viewport[3]-track[1].Y(),0.0f);
			glRotatef(20.0f,0.0f,0.0f,1.0f);
			glTranslatef(-track[1].X(),-viewport[3]+track[1].Y(),0.0f);
	
			glBegin(GL_LINES);
			glVertex2f(track[1].X(),viewport[3]-track[1].Y());
			glVertex2f(track[1].X()-arrowLenth*v.X(),viewport[3]-(track[1].Y()-arrowLenth*v.Y()));
			glEnd();
			glPopMatrix();

			glPushMatrix();
			glTranslatef(track[1].X(),viewport[3]-track[1].Y(),0.0f);
			glRotatef(-20.0f,0.0f,0.0f,1.0f);
			glTranslatef(-track[1].X(),-viewport[3]+track[1].Y(),0.0f);
			glBegin(GL_LINES);
			glVertex2f(track[1].X(),viewport[3]-track[1].Y());
			glVertex2f(track[1].X()-arrowLenth*v.X(),viewport[3]-(track[1].Y()-arrowLenth*v.Y()));
			glEnd();
			glPopMatrix();
		}
	}

	// Closing 2D
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

	glEnable(GL_LIGHTING) ;
}

void GLArea::draw3dMapTo2d()
{
	typedef lemon::ListGraph::Node Node ;
	typedef lemon::ListGraph::NodeIt NodeIt ;
	typedef lemon::ListGraph::Edge Edge ;
	typedef lemon::ListGraph::EdgeIt EdgeIt ;
	typedef lemon::ListGraph::IncEdgeIt IncEdgeIt ;

	typedef lemon::ListGraph::NodeMap<Point3f> P3fNodeMap ;

	lemon::ListGraph & g = *(dataMgr.skel.skeletonGraph);
	P3fNodeMap &nodeMap=*(dataMgr.skel.graphNodeP3f);

	for(EdgeIt e=EdgeIt(g);e!=lemon::INVALID;e++)
	{
		Point3f p1=nodeMap[g.u(e)],p2=nodeMap[g.v(e)];
		GLint viewport[4];   
		GLdouble modelview[16];    
		GLdouble projection[16];   
		GLdouble winX[2], winY[2], winZ[2];   
		glGetIntegerv(GL_VIEWPORT, viewport);   
		glGetDoublev(GL_MODELVIEW_MATRIX, modelview);    
		glGetDoublev(GL_PROJECTION_MATRIX, projection);    
		gluProject(p1.X(),p1.Y(),p1.Z(),modelview,projection,viewport,winX,winY,winZ);
		gluProject(p2.X(),p2.Y(),p2.Z(),modelview,projection,viewport,winX+1,winY+1,winZ+1);

		//winY[0]=viewport[3]-winY[0];
		//winY[1]=viewport[3]-winY[1];

		glEnable(GL_POINT_SMOOTH);   
		glEnable(GL_LINE_SMOOTH);   
		glEnable(GL_POLYGON_SMOOTH);

		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);


		glDisable(GL_LIGHTING) ;

		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glOrtho(0,width(),0,height(),-1,1);
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		glColor4f(0,0,1, 0.7);

		glLineWidth( 6.0 ) ;

	

		glBegin(GL_LINES);
		glVertex2f(winX[0],winY[0]);
		glVertex2f(winX[1],winY[1]);
		glEnd();

		glPopMatrix(); // restore modelview
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);

		glEnable(GL_LIGHTING) ;
	}
}
/*********************************added end*****************************************/

void GLArea::draw3dCurves( 	CurveArray3D &curves3d, GLColor color, bool closed, bool  randColor, double size )
{
	//using global function draw3dCurves()
	GlobalFun::draw3dCurves(curves3d,color, closed, randColor, size ) ;
	return ;
}


void GLArea::drawLightBall()
{
	// ============== LIGHT TRACKBALL ==============
	// Apply the trackball for the light direction

	glPushMatrix();
	trackball_light.GetView();
	trackball_light.Apply(!(isDefaultTrackBall()));

	static float lightPosF[]={0.0,0.0,1.0,0.0};
	glLightfv(GL_LIGHT0,GL_POSITION,lightPosF);
	static float lightPosB[]={0.0,0.0,-1.0,0.0};
	glLightfv(GL_LIGHT1,GL_POSITION,lightPosB);

	if (!(isDefaultTrackBall()))
	{
		glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);
		glColor3f(1,1,0);
		glDisable(GL_LIGHTING);
		const unsigned int lineNum=3;
		glBegin(GL_LINES);
		for(unsigned int i=0;i<=lineNum;++i)
			for(unsigned int j=0;j<=lineNum;++j) {
				glVertex3f(-1.0f+i*2.0/lineNum,-1.0f+j*2.0/lineNum,-2);
				glVertex3f(-1.0f+i*2.0/lineNum,-1.0f+j*2.0/lineNum, 2);
			}
			glEnd();
			glPopAttrib();
	}
	glPopMatrix();
}


void GLArea::changeColor(QString paraName)
{
	QColor qcolor ;
	

	if(!paraName.contains("Light"))
		qcolor = global_paraMgr.drawer.getColor(paraName);
	else
		qcolor = global_paraMgr.glarea.getColor(paraName);

	qcolor = QColorDialog::getColor(qcolor);

	if(qcolor.isValid()){

		if(paraName.contains("Light"))
		{
			GLColor color(qcolor);
			float light_col[4] = {color.r, color.g, color.b, 1.0};

			if(paraName.contains("Ambient"))
			{
				glLightfv(GL_LIGHT0, GL_AMBIENT, light_col);
				glLightfv(GL_LIGHT1, GL_AMBIENT, light_col);

			}
			else if(paraName.contains("Diffuse"))
			{
				glLightfv(GL_LIGHT0, GL_DIFFUSE, light_col);
				glLightfv(GL_LIGHT1, GL_DIFFUSE, light_col);

			}
			else if(paraName.contains("Specular"))
			{
				glLightfv(GL_LIGHT0, GL_SPECULAR, light_col);
				glLightfv(GL_LIGHT1, GL_SPECULAR, light_col);
			}
			para->setValue(paraName, ColorValue(qcolor));
		}

		if(!paraName.contains("Light"))
			global_paraMgr.drawer.setValue(paraName, ColorValue(qcolor));
		else
			global_paraMgr.glarea.setValue(paraName, ColorValue(qcolor));

	} 
}

int GLArea::pickPoint(int x, int y, vector<int> &result, int width, int height,bool only_one)
{
	//pick points,and save in vector result, return the number of points
	if(dataMgr.isSamplesEmpty())
		return -1;

	result.clear();

	if(width==0 ||height==0) return 0; 
	long hits;	
	CMesh* samples = dataMgr.getCurrentSamples();

	if (global_paraMgr.drawer.getBool("Use Pick Original"))
	{
		samples = dataMgr.getCurrentOriginal();
	}

	int sz= samples->vert.size();

	//if (global_paraMgr.drawer.getBool("Use Pick Skeleton"))
	//{
	//	sz = dataMgr.getCurrentSkeleton()->size;
	//}
	GLuint *selectBuf =new GLuint[sz*5];

	//  static unsigned int selectBuf[16384];
	glSelectBuffer(sz * 5, selectBuf);
	glRenderMode(GL_SELECT);
	glInitNames();

	/* Because LoadName() won't work with no names on the stack */
	glPushName(-1);

	double mp[16];

	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);
	glMatrixMode(GL_PROJECTION);
	glGetDoublev(GL_PROJECTION_MATRIX ,mp);
	glPushMatrix();
	glLoadIdentity();
	gluPickMatrix(x, y, width, height, viewport);
	glMultMatrixd(mp);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	//doPick = true;
	global_paraMgr.drawer.setValue("Doing Pick", BoolValue(true));


	//if (global_paraMgr.drawer.getBool("Use Pick Skeleton"))
	//{
	//	just_draw_skel_for_pick = true;
	//}
	paintGL();
	//just_draw_skel_for_pick = false;


	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	hits = glRenderMode(GL_RENDER);

	//xstring buf;
	if (hits <= 0)     return 0;

	vector<int> H;

	if (hits > 1 && global_paraMgr.drawer.getBool("Use Pick Mode2"))
	{
		hits--;
	}

	for(long ii=0;ii<hits;ii++)
	{
		//H.push_back( std::pair<double,unsigned int>(selectBuf[ii*4+1]/4294967295.0,selectBuf[ii*4+3]));
		H.push_back(selectBuf[ii*4+3]);
	}
	sort(H.begin(),H.end());


	//Only One Pick
	for(long i = 0 ;i < H.size();i++)
	{
		if(H[i] >= 0 && H[i] < sz)
		{
			result.push_back(H[i]);
			if(only_one)
				break;
		}
	}

	delete [] selectBuf;

	global_paraMgr.drawer.setValue("Doing Pick", BoolValue(false));

	if(only_one)
	{
		if(result.empty())
			return -1;
		else
		{
			return result[0];
		}
	}

	return result.size();
}


void GLArea::setView()
{
	glViewport(0,0, this->width(),this->height());
	curSiz.setWidth(this->width());
	curSiz.setHeight(this->height());

	GLfloat fAspect = (GLfloat)curSiz.width()/ curSiz.height();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// This parameter is the one that controls:
	// HOW LARGE IS THE TRACKBALL ICON ON THE SCREEN.
	float viewRatio = 1.75f;
	float cameraDist = viewRatio / tanf(math::ToRad(fov*.5f));

	if(fov==5)
		cameraDist = 1000; // small hack for orthographic projection where camera distance is rather meaningless...
	nearPlane = cameraDist - 2.f*clipRatioNear;
	farPlane =  cameraDist + 10.f*clipRatioFar;
	if(nearPlane<=cameraDist*.1f) nearPlane=cameraDist*.1f;

	if (!takeSnapTile)
	{
		if(fov==5)	glOrtho( -viewRatio*fAspect, viewRatio*fAspect, -viewRatio, viewRatio, cameraDist - 2.f*clipRatioNear, cameraDist+2.f*clipRatioFar);
		else    		gluPerspective(fov, fAspect, nearPlane, farPlane);
	}
	else	setTiledView(fov, viewRatio, fAspect, nearPlane, farPlane, cameraDist);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0, 0, -cameraDist,0, 0, 0, 0, 1, 0);
	camaraDistForDrawCurve = cameraDist ;

}

void GLArea::pasteTile()
{
	glPushAttrib(GL_ENABLE_BIT);
	QImage tileBuffer=grabFrameBuffer(true).mirrored(false,true);

	if (snapBuffer.isNull())
		snapBuffer = QImage(tileBuffer.width() * ss.resolution, tileBuffer.height() * ss.resolution, tileBuffer.format());

	uchar *snapPtr = snapBuffer.bits() + (tileBuffer.bytesPerLine() * tileCol) + ((totalCols * tileRow) * tileBuffer.byteCount());   //mark
	uchar *tilePtr = tileBuffer.bits();

	for (int y=0; y < tileBuffer.height(); y++)
	{
		memcpy((void*) snapPtr, (void*) tilePtr, tileBuffer.bytesPerLine());
		snapPtr+=tileBuffer.bytesPerLine() * totalCols;
		tilePtr+=tileBuffer.bytesPerLine();
	}

	tileCol++;

	if (tileCol >= totalCols)
	{
		tileCol=0;
		tileRow++;

		if (tileRow >= totalRows)
		{
			QString outfile=QString("%1/%2%3.png")
				.arg(ss.outdir)
				.arg(ss.basename)
				.arg("");
			//.arg(ss.counter++,2,10,QChar('0'));
			bool ret = (snapBuffer.mirrored(false,true)).save(outfile,"PNG");

			takeSnapTile=false;
			recoverView();

			snapBuffer=QImage();
		}
	}
	update();
	glPopAttrib();
}

void GLArea::recoverView()
{
	int w = width();
	int h = height();

	glViewport(0, 0, (GLint)w, (GLint)h);  
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();  

	float r = w/(float)h;
	gluPerspective(fov, r, 0.1, 10);
	glMatrixMode(GL_MODELVIEW);  
}

void GLArea::setTiledView(GLdouble fovY, float viewRatio, float fAspect, GLdouble zNear, GLdouble zFar,  float cameraDist)
{
	if(fovY<=5)
	{
		GLdouble fLeft   = -viewRatio*fAspect;
		GLdouble fRight  =  viewRatio*fAspect;
		GLdouble fBottom = -viewRatio;
		GLdouble fTop    =  viewRatio;

		GLdouble tDimX = fabs(fRight-fLeft) / totalCols;
		GLdouble tDimY = fabs(fTop-fBottom) / totalRows;


		glOrtho(fLeft   + tDimX * tileCol, fLeft   + tDimX * (tileCol+1),     /* left, right */
			fBottom + tDimY * tileRow, fBottom + tDimY * (tileRow+1),     /* bottom, top */
			cameraDist - 2.f*clipRatioNear, cameraDist+2.f*clipRatioFar);
	}
	else
	{
		GLdouble fTop    = zNear * tan(math::ToRad(fovY/2.0));
		GLdouble fRight  = fTop * fAspect;
		GLdouble fBottom = -fTop;
		GLdouble fLeft   = -fRight;

		// tile Dimension
		GLdouble tDimX = fabs(fRight-fLeft) / totalCols;
		GLdouble tDimY = fabs(fTop-fBottom) / totalRows;

		glFrustum(fLeft   + tDimX * tileCol, fLeft   + tDimX * (tileCol+1),
			fBottom + tDimY * tileRow, fBottom + tDimY * (tileRow+1), zNear, zFar);
	}
}

void GLArea::saveSnapshot()
{
	is_paintGL_locked = true;

	double SnapResolutionScale = para->getDouble("Snapshot Resolution");
	ss.resolution = SnapResolutionScale;
	snapDrawScal = 1. / SnapResolutionScale;
	
	totalCols = totalRows = ss.resolution;
	tileRow = tileCol = 0;
	ss.setOutDir(current_snap_path);
	ss.GetSysTime();

	takeSnapTile = true;
	updateGL();

	if (SnapResolutionScale > 1)
	{
		for (int i = 0; i < 4; i++)
		{
			is_paintGL_locked = false;
			updateGL();
			is_paintGL_locked = true;
		}
	}
	else
	{
		is_paintGL_locked = false;
		updateGL();
		is_paintGL_locked = true;
	}

	is_paintGL_locked = false;
}

void GLArea::setLassoCursor( ){

	QPixmap lassoCursorMap;

	if( lassoCursorMap.load("Icons/lasso2.png")){
		QCursor lassoCursor( lassoCursorMap , 1, 1) ;
		setCursor( lassoCursor ) ;
	}

}


void GLArea::resetLassoCursor( ){


	setCursor(Qt::ArrowCursor);
}

void GLArea::setContextAction(std::vector<QAction*> & acts) {

	contextActions = acts ;

}



//  [8/27/2014 CX]

void GLArea::set_operating_state(OperatingStates state_value){
	state=state_value;
	string operating_state;
	switch(state_value){
	case freeLassoState:
		operating_state="freeLassoState";
		break;
	case startState:
		operating_state="startState";
		break;
	case skelEditState:
		operating_state="skelEditState";
		break;
	case profEditMode:
		operating_state="profEditMode";
		break;
	case segmentMode:
		operating_state="segmentMode";
		break;
	case editResult:
		operating_state="editResult";
		break;
	case onesweep:
		operating_state="onesweep";
		break;
	case onesweepSegment:
		operating_state="onesweepSegment";
		break;
	case onesweepReady:
		operating_state="onesweepReady";
		break;
	case onesweepUpdateResult:
		operating_state="onesweepUpdateResult";
		break;
	case onesweepCorrespondence:
		operating_state="onesweepCorrespondence";
		break;
	}
	
	std::cout<<"set_operating_state:"<<operating_state<<endl;
}

OperatingStates GLArea::get_operating_state(){
	return state;
}

QString GLArea::get_operating_state_of_qstring(){
	QString statesofstring;
	switch(state){
	case freeLassoState:
		statesofstring=".....";
		break;
	case startState:
		statesofstring="[Start]: Press \'I\' to intialize morfit sweeper. Press \'E\' to edit the skeleton... ";
		break;
	case skelEditState:
		statesofstring="[Skeleton]: Press \'E\' again to exit. Ctrl+Z to undo. Ctrl+Y to redo.";
		break;
	case profEditMode:
		statesofstring=".....";
		break;
	case segmentMode:
		statesofstring=".....";
		break;
	case editResult:
		statesofstring=".....";
		break;
	case onesweep:
		statesofstring="[Stroke]: Draw strokes with your right mouse button (with ctrl/shift pressed). Green stroke to select sweeping path. Orange stroke to indicate sharp features.";
		break;
	case onesweepSegment:
		statesofstring="[Segment]: Draw one circle with your left button to refine the segmentation, with ctrl/shift pressed, just like photoshop. Right button is for enable/disable tips.  ";
		break;
	case onesweepReady:
		statesofstring="[Ready]: Press \'R\' to reconstruct the shape around the selected path. You can also double click the red profiles to edit it. Or, Press \'S\' to refine segmentation.  ";
		break;
	case onesweepUpdateResult:
		statesofstring="[Update]: draw a blue curve to select one profile. Edit it by drag the control points. Press \'R\' to reconstruct the shape again. Press \'I\' to go back to stroke state.";
		break;
	case onesweepCorrespondence:
		statesofstring=".....";
		break;
	}
	return statesofstring;

}