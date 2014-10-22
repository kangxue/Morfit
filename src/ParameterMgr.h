
#pragma once
#include "Parameter.h"
#include "CMesh.h"

class ParameterMgr
{
public:
	//parameter manager
	ParameterMgr(void);
	~ParameterMgr(void);

	//get the data members
	RichParameterSet* getDataParameterSet(){ return &data; }
	RichParameterSet* getDrawerParameterSet(){ return &drawer; }
	RichParameterSet* getGlareaParameterSet(){ return &glarea; }
	RichParameterSet* getConfigParameterSet(){ return &config;}
	RichParameterSet* getSkeletonParameterSet(){ return &skeleton; }	

	void setGlobalParameter(QString paraName,Value& val);

	//initialization for the parameter of config 
	void initConfigParameter();
private:
	//initialization for the parameter of "dataMgr"  
	void initDataMgrParameter();

	//initialization for the parameter of GLdrawer 
	void initDrawerParameter();

	//initialization for the parameter of GLarea 
	void initGlareaParameter();

	//initialization for the parameter of skeleton 
	void initSkeletonParameter();

public:
	//data member in parameter set
	RichParameterSet glarea;
	RichParameterSet data;
	RichParameterSet drawer;
	RichParameterSet skeleton;
	RichParameterSet config;
private:
	static int init_time;
	double grid_r;
};

extern ParameterMgr global_paraMgr;