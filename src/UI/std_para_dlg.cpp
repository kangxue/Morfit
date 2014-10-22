#include "UI/std_para_dlg.h"

StdParaDlg::StdParaDlg(ParameterMgr* _paras, GLArea * _area, QWidget* parent /* = 0 */) 
  : QDockWidget(QString("Plugin"), parent)
{
	paras = _paras;
	gla = _area;

	init(); // it's important
}

void StdParaDlg::init()
{
	para_config = NULL;
	mainFrame = NULL;
}

bool StdParaDlg::showConfigParaDlg()
{	
	createFrame();
	loadConfigParaDlg();
	return true;
}

void StdParaDlg::createFrame()
{
	if(mainFrame) delete mainFrame;
	QFrame *newFrame = new QFrame;
	setWidget(newFrame);
	setSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum);
	mainFrame = newFrame;	
	
}

void StdParaDlg::loadConfigParaDlg()
{
	assert(mainFrame);
	//mainFrame->hide();
	QGridLayout *gridLayout = new QGridLayout(mainFrame);
	//QVBoxLayout *vLayout = new QVBoxLayout(mainFrame);
	mainFrame->setLayout(gridLayout);
	setWindowTitle("Config");
	para_config = new ConfigParaDlg(this,paras,gla);
	para_config->setFrameConent();

	gridLayout->setRowMinimumHeight(2,620);
	gridLayout->setColumnMinimumWidth(1,315);
	gridLayout->addWidget(para_config,1,0,7,10);
	mainFrame->showNormal();
	mainFrame->adjustSize();
	//set the minimum size so it will shrink down to the right size
	this->setMinimumSize(mainFrame->sizeHint());

	this->showNormal();
	this->adjustSize();
}

void StdParaDlg::closeClick()
{
	cout << "close." << endl;
	close();
}

StdParaDlg::~StdParaDlg()
{
	cout << "De-construct StdParaDlg." << endl;
	// just set it to NULL
	gla = NULL;
	if(para_config)
		delete para_config;
	para_config=NULL;

	if(mainFrame)
		delete mainFrame;
	mainFrame = NULL;
}


