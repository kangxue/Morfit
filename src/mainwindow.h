


#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "GLArea.h"
#include <QMainWindow>
//#include "Algorithm/normal_extrapolation.h"
#include "ui_mainwindow.h"
#include "UI/std_para_dlg.h"
#include "UI/dlg_config_para.h"
#include "ParameterMgr.h"


class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = 0, Qt::WindowFlags flags = 0);
	~MainWindow();


private:
	void init();

	void initWidgets();
	void initConnect();
	void iniStatusBar();


private slots:
	void updateStatusBar();


	void openFile();

	void saveFile();

	
	void clearData();
	void normalizeData();

	void saveSnapshot();

	void showConfig(bool);

	//show things you want
	void showCutRadii(bool) ;

	void showTrajectory(bool) ;

	void showSeedProfile(bool) ;

	void showSkeleton(bool) ;

	void showMesh(bool) ;

	void showResultProfile(bool) ;

	void showSkeletonConfidence(bool) ;

	void showPoints(bool) ;

	void showCullFace(bool) ;

	void showMeshVertOnly(bool) ;

	void showSettings(bool) ;


	void showControlPoints(bool) ;
	


	//  the skeleton recenter
	void skelRecenter( bool);

	//call the function of zcop
	void undoOp_zc(bool);
	void redoOp_zc(bool);

	//save skeleton end with ".skel"
	void saveSkel_zc(bool);


	//execute the following by area->runCommand(...)
	void save_og_nbs(bool) ;
	void load_og_nbs(bool) ;
	void save_mesh(bool) ;
	void load_mesh(bool) ;
	void save_brch_mesh(bool) ;
	void save_settings(bool) ;
	void load_settings(bool) ;


	//call function keyReleaseEvent of area to deal with the key release event
	void keyReleaseEvent( QKeyEvent * e);

private:
	//data member  for class GLarea
	GLArea* area;

	//the title
	QString strTitle;

	//labels for status bar
	QLabel * sample_size_lable;
	QLabel * radius_label;
	QLabel * error_label;
	QLabel * state_lable;

	//data member for parameter manager
	ParameterMgr * paras;

	//call the parameter dialog
	StdParaDlg * paraDlg_Config;

	//ui for main windows
	Ui::mainwindowClass ui;
};

#endif // MAINWINDOW_H
