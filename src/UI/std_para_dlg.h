#pragma once

#include "glarea.h"
#include <iostream>
#include <QDialog>
#include <QtGui>
#include <QtWidgets/QFrame>
#include <QtWidgets/QWidget>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QGridLayout>


#include "ParameterMgr.h"
#include "UI/dlg_config_para.h"

using namespace std;

class StdParaDlg : public QDockWidget
{
	Q_OBJECT
	public:
		StdParaDlg(ParameterMgr* _paras, GLArea * _area, QWidget* parent = 0);
		~StdParaDlg();
		bool showConfigParaDlg();

	private:
		void init();
		void createFrame();
		void loadConfigParaDlg();

	private slots:
		 void closeClick();
		 
	private:
		ConfigParaDlg* para_config;
		ParameterMgr * paras;
		QFrame * mainFrame;
		GLArea * gla;
};
