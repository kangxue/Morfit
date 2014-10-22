

#pragma once
#include "gl/glew.h"
#pragma comment(lib,"glew32.lib")
#pragma comment(lib,"opengl32.lib")
#pragma comment(lib,"glu32.lib")

//#pragma comment(lib,"glut64.lib") // for 64bit
#pragma comment(lib,"glut32.lib") // for 32bit


#include <QWidget>
#include <QColor>
#include <QImage>
#include <QPoint>
#include <QtGui>
#include <QtOpenGL/QGLWidget>
#include<QString>

#include <QOpenGLFunctions>
#include <QOpenGLShader>
#include <QOpenGLShaderProgram>

#include <wrap/gl/trimesh.h>
#include <wrap/gui/trackball.h>
#include <vcg\space\point3.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <iostream>

//#include "Algorithm/anistropicPCA_Normal.h"
#include "DataMgr.h"
#include "GLDrawer.h"
#include "CMesh.h"
#include "ParameterMgr.h"
//#include "Algorithm/PointCloudAlgorithm.h"
#include "Operator_zc.h"

using std::cout;
using std::endl;
using vcg::Point3f;
using namespace vcg;
using std::vector;
class MainWindow ;

enum OperatingStates{
	startState,          
	freeLassoState,             
	skelEditState,            
	profEditMode,             
	segmentMode,         
	editResult,          
	onesweep,              
	onesweepSegment,    
	onesweepReady,     
	onesweepUpdateResult,     
	onesweepCorrespondence   
};

class GLArea : public QGLWidget
{
public:
	Q_OBJECT

public:
	GLArea(QWidget *parent = 0);
	~GLArea(void);

	void initializeGL();
	void resizeGL(int w, int h);

	void paintGL(); 
	
	void loadDefaultModel();

	void openByDrop(QString fileName);

	void initAfterOpenFile();
	void initView();
	void initSetting();


	void setContextAction(std::vector<QAction*> &) ;

	void setLassoCursor() ;
	void resetLassoCursor() ;

private:
	//void runPointCloudAlgorithm(PointCloudAlgorithm& algorithm);

signals:
	void needUpdateStatus();

private:
	void initLight();
	void lightOnOff(bool _val);

private: 
	bool activeDefaultTrackball; 
	bool isDefaultTrackBall()   { return activeDefaultTrackball; }
	vcg::Trackball trackball_light;
	void drawLightBall();

private: 
	int x1, y1, x2, y2;
	bool doPick;
	vector<int> pickList;

	int pickPoint(int x, int y, vector<int> &result, int width=4, int height=4, bool only_one = true);

public:
	bool isLeftPressed;
	bool isRightPressed;

	void drawPickRect();
	void draw2dCurves();

	void draw3dCurves( 	CurveArray3D &curves3d, GLColor color, bool closed=true, bool randColor=false, double size = 1.0 ) ;
	
	


private:
	int tileCol, tileRow, totalCols, totalRows;
	QImage snapBuffer;
	bool takeSnapTile;
	SnapshotSetting ss;

	void pasteTile();
	void setView(); 
	void recoverView();

	void setTiledView(GLdouble fovY, float viewRatio, float fAspect, GLdouble zNear, GLdouble zFar,  float cameraDist);
	QSize curSiz;
	float fov;
	float clipRatioFar;
	float clipRatioNear;
	float nearPlane;
	float farPlane;

	QString default_snap_path;
	QString current_snap_path;
	double snapDrawScal;
	bool is_paintGL_locked;

public:
	void saveSnapshot();
	void changeColor(QString paraName);
	void contextMenuEvent(QContextMenuEvent *);


private:
	vcg::Trackball trackball;
	vcg::Box3f gl_box;
	/*********************************controler begin*****************************************/
private:
	void wheelEvent(QWheelEvent *e);
	void mouseDoubleClickEvent(QMouseEvent * event);
	void mouseMoveEvent(QMouseEvent *e);
	void mousePressEvent(QMouseEvent *e);
	void mouseReleaseEvent(QMouseEvent *e);
public:
	void keyReleaseEvent ( QKeyEvent *e);
	/*********************************controler end*****************************************/
private:
	QMutex paintMutex;

	std::vector<QAction*> contextActions ;

public:
	DataMgr dataMgr;
	GLDrawer glDrawer;	
	RichParameterSet* para;

	OperatingStates state ;
	OperatingStates get_operating_state();
	void set_operating_state(OperatingStates state_value);
	QString get_operating_state_of_qstring();

	float projectMatrix[16] ;
	float modelviewMatrix[16] ;
	double camaraDistForDrawCurve ;

	Point3f screenRec[2] ;
	std::vector<bool> sampleIsSelected ;  
	std::vector<bool> sampleIsBoudary ; 
	std::vector<int> sampleClusterId ; 

	Operator_zc zcOp;
	void GLArea::draw2dCurves(vector<Point2f> track);
	void GLArea::draw3dMapTo2d();

	int SelectedNurbsId ;
	int SelectedCtrPointId ;

	std::vector<Point2f> curveForRefineSgement ;
	std::vector<Point2f> curveForFeatureSweep ;
	std::vector<Point2f> onesweepCurve ;

	// for shader
	void initShader() ;
	void SetupShaders( QString shaderName, std::vector<uniformValue> univalues );
	void bindShader( QString shaderName ) ;
	void releaseShader() ;
	std::vector<QOpenGLShaderProgram*>        shaderPrograms;
	std::vector<QString>        shaderNames;

	// for recover model view matrix
	float mvmatrix_loaded[16] ;
	bool mvmLoaded ;

	// command
	void runCommand( std::string cmd  ) ;

};//GLAERA_GLARERA_H

