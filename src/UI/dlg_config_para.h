#pragma once

#include "glarea.h"
#include <QtGui>
#include <QtWidgets/QFrame>
#include <QtWidgets/QWidget>
#include <iostream>

#include "ui_para_config.h"
#include "ParameterMgr.h"

using namespace std;


class ConfigParaDlg : public QFrame
{
	Q_OBJECT
public:
	ConfigParaDlg(QWidget *p, ParameterMgr * _paras, GLArea * _area);
	~ConfigParaDlg();
	void initConnects();
	void setFrameConent();
	public slots:	
		bool initWidgets();


signals:
	void parameterChanged();

	private slots:
		//get the value from the dialog
		void getData(double _val);
		void getRigid(double _val);
		void getSmooth(double _val);
		void getTolerance(double _val);

		//save the changes
		void applyChanges();
private:
	Ui::para_config * ui;
	ParameterMgr * m_paras;
	GLArea * area;
};