

#include "mainwindow.h"
#include <QFileDialog>
#include "appstate.h"
#include "para.h"

//this is a Global pointer for GLArea
GLArea *glarea_global_ptr ;


//std::stringstream buffer; 

MainWindow::MainWindow(QWidget *parent, Qt::WindowFlags flags)
: QMainWindow(parent, flags)
{

	
	
		//std::streambuf * old = std::cout.rdbuf(buffer.rdbuf()); 

	cout << "MainWindow constructed" << endl;
	ui.setupUi(this);
	area = new GLArea(this);
	setCentralWidget(area);

	init();
	initWidgets();
	iniStatusBar();
	initConnect();


	ReconstructorPara::cvNum = 10 ;
	ReconstructorPara::weight_sharpCV = 5 ;
	QString file = QString("gecko\/gecko.ply");


	if(!file.size()) return;
	
	appstate.plyfilename = file.toStdString() ;
	std::string path = appstate.plyfilename ;
	path.erase( path.begin()+path.find_last_of("\/"), path.end() ) ;
	appstate.path = path ;

	clearData() ;
	area->dataMgr.loadPlyToOriginal(file);
	area->dataMgr.loadSkeletonWithPly(file);
	normalizeData();
	area->initAfterOpenFile();
	area->updateGL();


	CMesh* samples = area->dataMgr.getCurrentOriginal();

	area->dataMgr.intializeSwp() ;

	glarea_global_ptr = area ;

}

MainWindow::~MainWindow()
{
	if(area) delete area;
	area = NULL;
}

void MainWindow::initWidgets()
{
	ui.actionShow_Samples->setChecked(global_paraMgr.glarea.getBool("Show Samples"));
	ui.actionShow_Original->setChecked(global_paraMgr.glarea.getBool("Show Original"));
	ui.actionShow_Normal_diff->setChecked(global_paraMgr.glarea.getBool("Show Normal Diff"));
	ui.actionShow_Tangent->setChecked(global_paraMgr.glarea.getBool("Show Tangent"));
	ui.actionShow_Vectors->setChecked(global_paraMgr.glarea.getBool("Show Vectors"));
	ui.actionShow_Neighborhood_Ball->setChecked(global_paraMgr.glarea.getBool("Show Radius"));
	ui.actionCull_Points->setChecked(global_paraMgr.drawer.getBool("Need Cull Points"));
}

void MainWindow::initConnect()
{
	if (!connect(area,SIGNAL(needUpdateStatus()),this,SLOT(updateStatusBar())))
	{
		cout << "can not connect signal" << endl;
	}
	//import the file of ply
	connect(ui.actionImport_Ply, SIGNAL(triggered()), this, SLOT(openFile()));
	//save the result of ply
	connect(ui.actionSave_Ply, SIGNAL(triggered()), this, SLOT(saveFile()));

	//snapshot
	connect(ui.actionSnapShot_2, SIGNAL(triggered()), this, SLOT(saveSnapshot()));


	
	connect( ui.actionConfig, SIGNAL(toggled(bool)), this, SLOT( showConfig(bool) ) ) ;



	connect(ui.actionShow_cutRadii,SIGNAL(toggled(bool)),this,SLOT(showCutRadii(bool)));
	connect(ui.actionShow_trajectory_2,SIGNAL(toggled(bool)),this,SLOT(showTrajectory(bool)));
	connect(ui.actionShow_seedProfile,SIGNAL(toggled(bool)),this,SLOT(showSeedProfile(bool)));
	connect(ui.actionShow_skeleton,SIGNAL(toggled(bool)),this,SLOT(showSkeleton(bool)));
	connect(ui.actionShow_mesh,SIGNAL(toggled(bool)),this,SLOT(showMesh(bool)));
	connect(ui.actionShow_resultProfile,SIGNAL(toggled(bool)),this,SLOT(showResultProfile(bool)));
	connect(ui.actionShow_skeletonConfidence,SIGNAL(toggled(bool)),this,SLOT(showSkeletonConfidence(bool)));
	connect(ui.actionShow_points,SIGNAL(toggled(bool)),this,SLOT(showPoints(bool)));
	connect(ui.actionShow_cullFace,SIGNAL(toggled(bool)),this,SLOT(showCullFace(bool)));
	connect(ui.actionShow_meshVertOnly,SIGNAL(toggled(bool)),this,SLOT(showMeshVertOnly(bool)));
	connect(ui.actionShow_settings,SIGNAL(toggled(bool)),this,SLOT(showSettings(bool)));
	
	connect(ui.actionShowCv,SIGNAL(toggled(bool)),this,SLOT(showControlPoints(bool)));


	// skeleton edit
	connect( ui.actionEditMode, SIGNAL(toggled(bool)), this, SLOT( flipEditMode(bool) ) ) ; // mark
	connect( ui.actionUndo, SIGNAL(triggered(bool)), this, SLOT( undoOp_zc(bool) ) ) ;
	connect( ui.actionRedo, SIGNAL(triggered(bool)), this, SLOT( redoOp_zc(bool) ) ) ;
	connect(ui.actionSave_Skel,SIGNAL(triggered(bool)),this,SLOT(saveSkel_zc(bool)));


	connect( ui.actionRecenter, SIGNAL(triggered(bool) ), this, SLOT( skelRecenter(bool) ) ) ;



	// load and save ongoing nurbs
	connect( ui.actionSave_og_nbs, SIGNAL(triggered(bool) ), this, SLOT( save_og_nbs(bool) ) ) ;
	connect( ui.actionLoad_og_nbs, SIGNAL(triggered(bool) ), this, SLOT( load_og_nbs(bool) ) ) ;


	// load and save mesh
	connect( ui.actionSave_mesh, SIGNAL(triggered(bool) ), this, SLOT( save_mesh(bool) ) ) ;
	connect( ui.actionLoad_mesh, SIGNAL(triggered(bool) ), this, SLOT( load_mesh(bool) ) ) ;
	connect( ui.actionSave_brch_mesh, SIGNAL(triggered(bool) ), this, SLOT( save_brch_mesh(bool) ) ) ;

	// load and save settings
	connect( ui.actionSave_settings, SIGNAL(triggered(bool) ), this, SLOT( save_settings(bool) ) ) ;
	connect( ui.actionLoad_settings, SIGNAL(triggered(bool) ), this, SLOT( load_settings(bool) ) ) ;


}

void MainWindow::iniStatusBar()
{
	//Initializing the status bar
	QStatusBar *status_bar = statusBar();

	state_lable = new QLabel;
	state_lable->setMinimumSize(200,30);
	state_lable->setFrameShape(QFrame::NoFrame);
	state_lable->setFrameShadow(QFrame::Plain);

	radius_label = new QLabel;
	radius_label->setMinimumSize(200,30);
	radius_label->setFrameShape(QFrame::NoFrame);
	radius_label->setFrameShadow(QFrame::Plain);
	
	sample_size_lable = new QLabel;
	sample_size_lable->setMinimumSize(200,30);
	sample_size_lable->setFrameShape(QFrame::NoFrame);
	sample_size_lable->setFrameShadow(QFrame::Plain);

	error_label = new QLabel;
	error_label->setMinimumSize(200,30);
	error_label->setFrameShape(QFrame::NoFrame);
	error_label->setFrameShadow(QFrame::Plain);

	state_lable = new QLabel;
	state_lable->setMinimumSize(200,30);
	state_lable->setFrameShape(QFrame::NoFrame);
	state_lable->setFrameShadow(QFrame::Plain);

	updateStatusBar();
	status_bar->addWidget(state_lable);
	
}

void MainWindow::updateStatusBar()
{
	//update the information of the status bar
	QString title = strTitle +  " - " + area->dataMgr.curr_file_name;
	if( area->dataMgr.curr_file_name.isEmpty() )
		setWindowTitle(strTitle);
	else
		setWindowTitle(title);

	int s_size = 0;
	if(!area->dataMgr.isSamplesEmpty())
		s_size = area->dataMgr.getCurrentSamples()->vert.size();

	QString str = "Sample points: " + QString::number(s_size);
	sample_size_lable->setText(str);

	double radius  = global_paraMgr.data.getDouble("CGrid Radius");
	QString strRadius = "Radius: " + QString::number(radius);
	radius_label->setText(strRadius);

	double error = 0;
	QString running_name = global_paraMgr.glarea.getString("Running Algorithm Name");

	QString strError = "Movement: " + QString::number(error);
	error_label->setText(strError);

	QString stateofop = area->get_operating_state_of_qstring();
	state_lable->setText(stateofop);

	update();
	repaint();
}

void MainWindow::init()
{
	strTitle = "Morfit_v1.0";

	paraDlg_Config = NULL;
	paras = &global_paraMgr;
}

void MainWindow::showConfig( bool t)
{

	if(paraDlg_Config != 0)
	{
		paraDlg_Config->close();
		delete paraDlg_Config;
		paraDlg_Config = 0 ;
		cout<<"delete paraDlg_Config"<<endl;
	}


	if( !t)
		return ;

	paras->initConfigParameter();
	paraDlg_Config = new StdParaDlg(paras, area, this);
	paraDlg_Config->setAllowedAreas(Qt::LeftDockWidgetArea
		| Qt::RightDockWidgetArea);
	addDockWidget(Qt::RightDockWidgetArea,paraDlg_Config);

	paraDlg_Config->setFloating(false);
	paraDlg_Config->hide();

	paraDlg_Config->showConfigParaDlg();
}



void MainWindow::openFile()
{
	QString file = QFileDialog::getOpenFileName(this, "Select a ply file", "", "*.ply");
	if(!file.size()) return;

	appstate.plyfilename = file.toStdString() ;
	std::string path = appstate.plyfilename ;
	path.erase( path.begin()+path.find_last_of("\/"), path.end() ) ;
	appstate.path = path ;

	clearData() ;

	// get possion filename
	std::string possionfile = file.toStdString() ;
	possionfile.insert( possionfile.size()-4, "_possion") ;

	area->dataMgr.loadPlyToOriginal(file);
	area->dataMgr.loadSkeletonWithPly(file);

	normalizeData();

	area->initAfterOpenFile();
	area->updateGL();


	CMesh* samples = area->dataMgr.getCurrentOriginal();
	area->dataMgr.intializeSwp() ;
}

void MainWindow::normalizeData()
{
	area->dataMgr.normalizeAllMesh();
	area->initView();
	area->updateGL();
}

void MainWindow::clearData()
{
	area->dataMgr.clearData();
	area->sampleIsSelected.clear() ;
	area->sampleIsBoudary.clear() ;
	area->sampleClusterId.clear();


	area->updateGL();

}

void MainWindow::saveFile()
{
	QString file = QFileDialog::getSaveFileName(this, "Save samples as", "", "*.ply");
	if(!file.size()) return;

	area->dataMgr.savePly(file, *area->dataMgr.getCurrentOriginal());

	return ;
}



void MainWindow::saveSnapshot()
{
	area->saveSnapshot();
}





void MainWindow::keyReleaseEvent( QKeyEvent * e){

	area->keyReleaseEvent( e ) ;

	if( e->key() == Qt::Key_Escape &&  area->get_operating_state() == freeLassoState && !area->isLeftPressed ) {
		
		ui.actionLasso->setChecked(false) ;

	}else if( e->key() == Qt::Key_Escape &&  area->get_operating_state() != freeLassoState ) {
		area->dataMgr.curves3d.clear() ;
		for( int i=0; i<area->sampleIsSelected.size(); ++i )
			area->sampleIsSelected[i] = false ;
	}

	extern AppState appstate ;
	if( e->key() == Qt::Key_J)
		appstate.displayJoint = !appstate.displayJoint ;

	updateStatusBar();


}



void MainWindow::showCutRadii(bool val) {

	extern AppState appstate ;
	appstate.displayCutRadii = val ;
}

//add

void MainWindow::showTrajectory(bool val) {

	extern AppState appstate ;
	appstate.displayTrajectory = val ;
}


void MainWindow::showSeedProfile(bool val) {

	extern AppState appstate ;
	appstate.displaySeedProfile = val ;
}


void MainWindow::showSkeleton(bool val) {

	extern AppState appstate ;
	appstate.displaySkeleton = val ;
}


void MainWindow::showMesh(bool val) {

	extern AppState appstate ;
	appstate.displayMesh = val ;
}


void MainWindow::showResultProfile(bool val) {

	extern AppState appstate ;
	appstate.displayResultProfile = val ;
}


void MainWindow::showSkeletonConfidence(bool val) {

	extern AppState appstate ;
	appstate.displaySkeletonConfidence = val ;
}


void MainWindow::showPoints(bool val) {

	extern AppState appstate ;
	appstate.displayPoints = val ;
}


void MainWindow::showCullFace(bool val) {

	extern AppState appstate ;
	appstate.cullFace = val ;
}


void MainWindow::showMeshVertOnly(bool val) {

	extern AppState appstate ;
	appstate.displayMeshVertOnly = val ;
}


void MainWindow::showSettings(bool val) {

	extern AppState appstate ;
	appstate.displaySettings = val ;
}
//end

void MainWindow::showControlPoints(bool val) {

	extern AppState appstate ;

	appstate.displaySettings = !appstate.displaySettings ;

	int skelID = 0;
	skelpath skel = area->dataMgr.swp.ongoingSkel ;
	appstate.skelProfs.clear() ;
	while( ReconstructorUtility::loadSkeletonAndProfiles(skel, skelID)) {
		appstate.skelProfs.push_back( skel) ;
		skelID++ ;

	}
}



void MainWindow::save_og_nbs(bool trigged) {
	area->runCommand("save og nbs") ;
}
void MainWindow::load_og_nbs(bool trigged) {
	area->runCommand("load og nbs") ;
}

void MainWindow::save_mesh(bool trigged) {
	area->runCommand("save mesh") ;
}
void MainWindow::load_mesh(bool trigged) {
	area->runCommand("load mesh") ;
}



void MainWindow::save_settings(bool trigged) {
	area->runCommand("save settings") ;
}
void MainWindow::load_settings(bool trigged) {
	area->runCommand("load settings") ;
}




void MainWindow::save_brch_mesh(bool trigged) {
	area->runCommand("save brch mesh") ;
}
void MainWindow::skelRecenter( bool trigged){

	area->runCommand("recenter") ;
}



void MainWindow::undoOp_zc(bool val)
{
	area->zcOp.undo(area->dataMgr.skel);
}
void MainWindow::redoOp_zc(bool val)
{
	area->zcOp.redo(area->dataMgr.skel);
}
void MainWindow::saveSkel_zc(bool val)
{
	QString file = QFileDialog::getSaveFileName(this, "Save samples as", "", "*.skel");
	if(!file.size()) return;

	area->zcOp.saveSkeletonAsSkel(area->dataMgr.skel,file);

	area->updateGL();
}

