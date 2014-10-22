
#include "ParameterMgr.h"
#include <iostream>
#include <fstream>
#include "appstate.h"


int ParameterMgr::init_time = 0;
ParameterMgr global_paraMgr;

ParameterMgr::ParameterMgr(void)
{
	init_time++;
	if(init_time > 1)
	{
		std::cout << "can not init ParameterMgr twice!" << endl;
		return;
	}

	grid_r = 0.05;

	initDataMgrParameter();
	initDrawerParameter();		
	initGlareaParameter();
	initSkeletonParameter();
}

ParameterMgr::~ParameterMgr(void)
{
}





void ParameterMgr::setGlobalParameter(QString paraName,Value& val)
{
	if(glarea.hasParameter(paraName))
		glarea.setValue(paraName, val);
	if(data.hasParameter(paraName))
		data.setValue(paraName, val);
	if(drawer.hasParameter(paraName))
		drawer.setValue(paraName, val);
	if (skeleton.hasParameter(paraName))
		skeleton.setValue(paraName, val);
	if (config.hasParameter(paraName))
		config.setValue(paraName,val);
}

//initialization for the parameter of "dataMgr"  
void ParameterMgr::initDataMgrParameter()
{
	//data.addParam(new RichInt("Sub Sample", 1));
	data.addParam(new RichDouble("Init Radius Para", 1.0));
	data.addParam(new RichDouble("Down Sample Num", 10000));
	data.addParam(new RichDouble("CGrid Radius", grid_r));
}


//initialization for the parameter of GLarea 
void ParameterMgr::initGlareaParameter()
{
	glarea.addParam(new RichString("Running Algorithm Name", "") );

	glarea.addParam(new RichBool("Light On or Off", false) );

	glarea.addParam(new RichBool("Show Normal", false) );
	glarea.addParam(new RichBool("Show Normal Diff", false) );
	glarea.addParam(new RichBool("Show Tangent", false) );
	glarea.addParam(new RichBool("Show Vectors", false) );
	
	glarea.addParam(new RichBool("Show Samples", true) );
	glarea.addParam(new RichBool("Show Samples Quad", false) );
	glarea.addParam(new RichBool("Show Samples Dot", true) );
	glarea.addParam(new RichBool("Show Samples Circle", false) );
	glarea.addParam(new RichBool("Show Samples 2D", false) );
	glarea.addParam(new RichBool("Render Sample Disks", false) );

	glarea.addParam(new RichBool("Show Original", true) );
	glarea.addParam(new RichBool("Show Original Quad", false) );
	glarea.addParam(new RichBool("Show Original Dot", true) );
	glarea.addParam(new RichBool("Show Original Circle", false) );
	glarea.addParam(new RichBool("Show Skeleton", true));



	glarea.addParam(new RichBool("Show Radius", false));
	glarea.addParam(new RichBool("Show Radius Use Pick", true));
	glarea.addParam(new RichBool("Show Red Radius Line", true));
	glarea.addParam(new RichBool("Multiply Pick Point", true) );

	glarea.addParam(new RichBool("GLarea Busying", false) );


	glarea.addParam(new RichPoint3f("Light Position", vcg::Point3f(-4.0, -4.0, -4.0)));
	glarea.addParam(new RichColor("Light Ambient Color", QColor(85, 85, 85)));
	glarea.addParam(new RichColor("Light Diffuse Color", QColor(241, 241, 241)));
	glarea.addParam(new RichColor("Light Specular Color", QColor(255, 255, 255)));



	//glarea.addParam(new RichPoint3f("Light Position", vcg::Point3f(2.0, 4.0, -4.0)));
	//glarea.addParam(new RichColor("Light Ambient Color", QColor(32, 32, 32)));  
	//glarea.addParam(new RichColor("Light Diffuse Color", QColor(204, 204, 204)));
	//glarea.addParam(new RichColor("Light Specular Color", QColor(255, 255, 255)));


	glarea.addParam(new RichDouble("Snapshot Resolution", 4));
}

//initialization for the parameter of GLdrawer
void ParameterMgr::initDrawerParameter()
{
	drawer.addParam(new RichBool("Doing Pick", false));
	drawer.addParam(new RichBool("Need Cull Points", false) );
	drawer.addParam(new RichBool("Use Pick Original", false));
	drawer.addParam(new RichBool("Use Pick Mode2", true) );
	drawer.addParam(new RichBool("Skeleton Light", true));

	drawer.addParam(new RichDouble("Original Draw Width", 0.0015));
	drawer.addParam(new RichDouble("Sample Draw Width", 0.0030));
	drawer.addParam(new RichDouble("Quade Draw Width", 0.0030));
	drawer.addParam(new RichDouble("Sample Dot Size", 2));
	drawer.addParam(new RichDouble("Original Dot Size", 1));
	drawer.addParam(new RichDouble("Normal Line Width", 2));
	drawer.addParam(new RichDouble("Normal Line Length", 0.05));

	drawer.addParam(new RichColor("Background Color", QColor(255, 255, 255) ));
	drawer.addParam(new RichColor("Normal Line Color", QColor(0, 0, 255) ));
	drawer.addParam(new RichColor("Sample Point Color", QColor(255, 0, 0) ));
	drawer.addParam(new RichColor("Original Point Color", QColor(48, 48, 48) ));
	drawer.addParam(new RichColor("Feature Color", QColor(0, 0, 255) ));
	drawer.addParam(new RichColor("Pick Point Color", QColor(128, 128, 0) ));
	drawer.addParam(new RichColor("Pick Point DNN Color", QColor(0, 0, 155) ));

	drawer.addParam(new RichColor("Skeleton Bone Color", QColor(200, 0, 0) ));
	drawer.addParam(new RichColor("Skeleton Node Color", QColor(50, 250, 50) ));
	drawer.addParam(new RichColor("Skeleton Branch Color", QColor(0, 0, 0)));
	drawer.addParam(new RichDouble("Skeleton Bone Width", 100)); // ./10000
	drawer.addParam(new RichDouble("Skeleton Node Size", 120)); // ./10000
	drawer.addParam(new RichDouble("Skeleton Branch Size", 30)); // ./10000
}

//initialization for the parameter of skeleton 
void ParameterMgr::initSkeletonParameter()
{
	skeleton.addParam(new RichString("Algorithm Name", "Skeletonization") );
	skeleton.addParam(new RichDouble("Num Of Iterate Time", 1));

	skeleton.addParam(new RichDouble("CGrid Radius", grid_r));
	skeleton.addParam(new RichDouble("H Gaussian Para", 4));
	skeleton.addParam(new RichDouble("Repulsion Power", 1.0));
	skeleton.addParam(new RichDouble("Average Power", 2.0));
	skeleton.addParam(new RichBool("Need Compute Density", true));
	
	
	skeleton.addParam(new RichDouble("Current Movement Error", 0.0));
	skeleton.addParam(new RichBool("Run Auto Wlop One Step", false));
	skeleton.addParam(new RichBool("Run Auto Wlop One Stage", false));
	skeleton.addParam(new RichBool("Process Should Stop", false));

	skeleton.addParam(new RichBool("Step1 Detect Skeleton Feature", false));
	skeleton.addParam(new RichBool("Step2 Run Search New Branchs", false));
	skeleton.addParam(new RichBool("Step3 Clean And Update Radius", false));


	//init
	skeleton.addParam(new RichDouble("Max Iterate Time", 25));
	skeleton.addParam(new RichDouble("Stop And Grow Error", 0.0001));
	skeleton.addParam(new RichDouble("Initial Radius", -1.));
	skeleton.addParam(new RichDouble("Radius Update Speed", 0.5));

	//step0
	skeleton.addParam(new RichDouble("Repulsion Mu", 0.35));
	skeleton.addParam(new RichDouble("Repulsion Mu2", 0.15));
	skeleton.addParam(new RichDouble("Follow Sample Radius", 0.5));
	skeleton.addParam(new RichDouble("Follow Sample Max Angle", 80));// should add to UI
	skeleton.addParam(new RichDouble("Inactive And Keep Virtual Angle", 60)); // should add to UI
	skeleton.addParam(new RichDouble("Save Virtual Angle", 30)); // should add to UI

	skeleton.addParam(new RichDouble("Grow Accept Sigma", 0.8));// should add to UI
	skeleton.addParam(new RichDouble("Bad Virtual Angle", 110));// 1-15


	//step1
	skeleton.addParam(new RichDouble("Combine Too Close Threshold", 0.01));
	skeleton.addParam(new RichDouble("Sigma KNN", 6));//this one is hart to choose, should be small for narrow region, but will lead to unnecessary small branches
	skeleton.addParam(new RichDouble("Eigen Feature Identification Threshold", 0.901));

	//step2
	skeleton.addParam(new RichDouble("Branches Search Angle", 25));
	skeleton.addParam(new RichDouble("Virtual Head Accecpt Angle", 25));
	skeleton.addParam(new RichDouble("Snake Search Max Dist Blue", 0.4));
	skeleton.addParam(new RichDouble("Accept Branch Size", 4)); // important
	skeleton.addParam(new RichDouble("Branch Search Max Dist Yellow", 0.1));

	skeleton.addParam(new RichDouble("Branches Merge Max Dist", 0.06));

	skeleton.addParam(new RichDouble("Two Virtual Very Close Dist", 0.005));
	skeleton.addParam(new RichDouble("Branch Search KNN", 10));
	skeleton.addParam(new RichDouble("Combine Similar Angle", 150));
	skeleton.addParam(new RichDouble("Grow Search Radius", 0.15));

	skeleton.addParam(new RichDouble("Add Accept Branch Size", 0));


	//step3
	skeleton.addParam(new RichDouble("Clean Near Branches Dist", 0.05));
	skeleton.addParam(new RichDouble("Fix Original Weight", 0.91));
	skeleton.addParam(new RichDouble("Curve Segment Length", 0.08));
	skeleton.addParam(new RichInt("Fix Original Mode", 4)); // 1 for noisy , 4 for clean

	//strategy
	skeleton.addParam(new RichBool("Use Nearby Combine Strategy", true));
	skeleton.addParam(new RichBool("Use Go Through Strategy", false));
	skeleton.addParam(new RichBool("Use Aggresive Growth Strategy", false));
	skeleton.addParam(new RichBool("Use Clean Points When Following Strategy", true));
	skeleton.addParam(new RichBool("Use All Connect Strategy", true));
	skeleton.addParam(new RichBool("Use Plus Perpendicular Dist Strategy", false));
	skeleton.addParam(new RichBool("Use Kill Too Close Strategy", false));
	skeleton.addParam(new RichBool("Use Compute Eigen Ignore Branch Strategy", false));
	skeleton.addParam(new RichBool("Use Virtual Group Merge Strategy", false));
	skeleton.addParam(new RichBool("Use Final Merge Strategy", true));
	skeleton.addParam(new RichBool("Use Search New Twice Strategy", false));
	skeleton.addParam(new RichBool("Inactive Overlap Strategy", false));
	skeleton.addParam(new RichBool("Move Overlap Strategy", false));
	skeleton.addParam(new RichBool("Use Virtual Near Body Stop Strategy", false));
	skeleton.addParam(new RichBool("Need To Keep Big Bug", false));
	skeleton.addParam(new RichBool("2D Leaf", false));




	skeleton.addParam(new RichDouble("Change Strategy Radius", 0.45));
	skeleton.addParam(new RichBool("Is Radius Big", false));


	//wLop.addParam(new RichBool("Run Grow Pick Curve", false));
	//wLop.addParam(new RichBool("Run Grow Pick Curve", false));


	skeleton.addParam(new RichDouble("All Connect Dist", 0.25)); //important
	skeleton.addParam(new RichDouble("All Connect Angle", 30)); //important


	//skeleton.addParam(new RichBool("Need Compute Density", true));
	skeleton.addParam(new RichBool("Need Recentering", true));
	skeleton.addParam(new RichBool("Run Auto Skeleton", false));
	skeleton.addParam(new RichBool("Not Load Parameter", false));
	skeleton.addParam(new RichBool("Have To Keep Virtual", false));




}

//initialization for the parameter of config 
void ParameterMgr::initConfigParameter()
{
	double t_data=1000;
	double t_rigid=100;
	double t_smooth=1;
	double t_tmpconstraint=0;
	if(appstate.plyfilename=="")
	{
		cout<<"no ply name!"<<endl;
		return;
	}

	//std::string dir = appstate.plyfilename ;
	//dir.erase( dir.begin() + dir.find_last_of("\/"), dir.end()  );
	//std::string parafname = dir + "\\config.txt" ;

	//std::ifstream paraifs(parafname) ;
	//if( paraifs.good() ) paraifs >> t_data ;
	//if( paraifs.good() ) paraifs >> t_rigid ;
	//if( paraifs.good() ) paraifs >> t_smooth;
	//if( paraifs.good() ) paraifs >> t_tmpconstraint ;

	if(config.hasParameter("Algorithm Name")==false)
	{
		config.addParam(new RichString("Algorithm Name", "config") );
	}
	
	if(config.hasParameter("data")==false)
	{
		config.addParam(new RichDouble("data", t_data ));
		cout<<"t_data"<<t_data<<endl;
	}

	if(config.hasParameter("rigid")==false)
	{
		config.addParam(new RichDouble("rigid", t_rigid ));
		cout<<"t_rigid"<<t_rigid<<endl;

	}
	if(config.hasParameter("smooth")==false)
	{
		config.addParam(new RichDouble("smooth", t_smooth ));
		cout<<"t_smooth"<<t_smooth<<endl;

	}
	if(config.hasParameter("tmpconstraint")==false)
	{
		config.addParam(new RichDouble("tmpconstraint", t_tmpconstraint ));
		cout<<"t_tmpconstraint"<<t_tmpconstraint<<endl;

	}


	if(config.hasParameter("Tolerance")==false)
	{
		config.addParam(new RichDouble("Tolerance", 0.01 ));
		cout<<"Tolerance"<<0.01<<endl;

	}
}