#include "UI/dlg_config_para.h"
#include "appstate.h"
#include <fstream>
#include "reconstructionUtility.h"
#include "reconstructorPara.h"
#include "para.h"
#include <iostream>
extern AppState appstate ;
ConfigParaDlg::ConfigParaDlg(QWidget *p, ParameterMgr * _paras, GLArea * _area) : QFrame(p)
{
	ui = new Ui::para_config;
	ConfigParaDlg::ui->setupUi(this);
	m_paras = _paras;
	area = _area;
	if(!initWidgets())
	{
		cerr << "Warning:  ConfigParaDlg::initWidgets failed!" << endl;
		return ;
	}
	initConnects();
}


void ConfigParaDlg::initConnects()
{
	if (!connect(area,SIGNAL(needUpdateStatus()),this,SLOT(initWidgets())))
	{
		cout << "can not connect signal" << endl;
	}
	if(!connect((ui->data),SIGNAL(valueChanged(double)),this,SLOT(getData(double))))
	{
		cout << "cannot connect ConfigParaDlg::getData(double)." << endl;
	}
	if(!connect((ui->rigid),SIGNAL(valueChanged(double)),this,SLOT(getRigid(double))))
	{
		cout << "cannot connect ConfigParaDlg::getRigid(double)." << endl;
	}
	if(!connect((ui->smooth),SIGNAL(valueChanged(double)),this,SLOT(getSmooth(double))))
	{
		cout << "cannot connect ConfigParaDlg::getSmooth(double)." << endl;
	}
	if(!connect((ui->tolerance),SIGNAL(valueChanged(double)),this,SLOT(getTolerance(double))))
	{
		cout << "cannot connect ConfigParaDlg::getTolerance(double)." << endl;
	}

	//if(!connect(ui->radius,SIGNAL(valueChanged(double)),this,SLOT(getRadiusValues(double))))
	//{
	//	cerr << "cannot connect WlopParaDlg::getDoubleValues(double)." << endl;
	//}

	//
	////
	//if(!connect(ui->wlop_apply,SIGNAL(clicked()),this,SLOT(applyWlop())))
	//{
	//	cerr << "cannot connect WlopParaDlg::applyWlop()." << endl;
	//}
	if(!connect(ui->okButton,SIGNAL(clicked()),this,SLOT(applyChanges())))
	{
		cout<<"can not connect ConfigParaDlg::applyChanges()"<<endl;
	}

}

bool ConfigParaDlg::initWidgets()
{
	ui->data->setValue(m_paras->config.getDouble("data"));
	ui->rigid->setValue(m_paras->config.getDouble("rigid"));
	ui->smooth->setValue(m_paras->config.getDouble("smooth"));
	ui->tolerance->setValue(m_paras->config.getDouble("Tolerance"));
	return true;
}


void ConfigParaDlg::setFrameConent()
{
	if(layout()) delete layout();
	showNormal();
	adjustSize();
}

void ConfigParaDlg::getData(double _val)
{
	m_paras->config.setValue("data",DoubleValue(_val));
	applyChanges();
	area->updateGL();
}

void ConfigParaDlg::getRigid(double _val)
{
	m_paras->config.setValue("rigid",DoubleValue(_val));
	applyChanges();
	area->updateGL();
}

void ConfigParaDlg::getSmooth(double _val)
{
	m_paras->config.setValue("smooth",DoubleValue(_val));
	applyChanges();
	area->updateGL();
}

void ConfigParaDlg::getTolerance(double _val)
{
	m_paras->config.setValue("Tolerance",DoubleValue(_val));
	applyChanges();
	area->updateGL();
}

void ConfigParaDlg::applyChanges()
{
	/*std::string dir = appstate.plyfilename ;
	dir.erase( dir.begin() + dir.find_last_of("\/"), dir.end()  );
	std::string parafname = dir + "\\config.txt" ;

	std::ofstream paraofs(parafname) ;
	if(paraofs.good())
	{
		paraofs<<m_paras->config.getDouble("data")<<endl
			<<m_paras->config.getDouble("rigid")<<endl
			<<m_paras->config.getDouble("smooth")<<endl
			<<m_paras->config.getDouble("tmpconstraint")<<endl
			<<"1e5"<<endl;
		cout<<"success ofstream"<<endl;
		cout<<parafname<<endl;
	}
	paraofs.clear();
	paraofs.close();
*/

	completionSolverPara::w_data=m_paras->config.getDouble("data");
	completionSolverPara::w_rigid=m_paras->config.getDouble("rigid");
	completionSolverPara::w_smooth=m_paras->config.getDouble("smooth");
	ReconstructorPara::MMA_FTOL=m_paras->config.getDouble("Tolerance");


	std::cout<<"-------------------------------------------------"<<std::endl;
	std::cout<< "applyChanges called!" 
		<<"\nw_data = " << completionSolverPara::w_data 
		<<"\nw_rigid = " << completionSolverPara::w_rigid 
		<<"\nw_smooth = " << completionSolverPara::w_smooth 
		<<"\nMMA_FTOL = " << ReconstructorPara::MMA_FTOL  
		<<std::endl;
	std::cout<<"-------------------------------------------------"<<std::endl;

}

ConfigParaDlg::~ConfigParaDlg()
{
	cout << "De-construct ConfigParaDlg Frame." << endl;
	delete ui;
	ui = NULL;
}