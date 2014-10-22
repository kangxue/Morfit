/*
	This file is part of the software described in the paper:
	Kangxue Yin, Hui Huang, Hao(Richard) Zhang, Minglun Gong, Daniel Cohen-or, Baoquan Chen. 
	Morfit: Interactive Surface Reconstruction from Incomplete Point Clouds with Curve-Driven Topology and Geometry Control. 
	ACM Transactions on Graphics 33(6)(Proc.SIGGRAPH ASIA 2014 )

	Copyright (C) <2013-2014>  <Kangxue Yin - yinkangxue@gmail.com>
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/



#ifndef _nrigidtransform_solver_
#define  _nrigidtransform_solver_

//#define  _LargeScope_ 

#include "skelpath.h"
#include "GlobalFunction.h"
#include "reconstructionUtility.h"
#include "reconstructorPara.h"
#include <algorithm>



class completionSolver{
public:
	completionSolver() {}

	completionSolver(const skelpath &skel, int branId, bool isTip_first, bool IsTip_last, std::vector<int2> corres_cues );

	// for save and load
	std::string cvtToString() ;
	completionSolver( std::string str) ;


	std::vector<double> solve(  std::vector<double> x0_from_caller ) ;
	std::vector<double> resolve( int newTmpId, ON_NurbsCurve newTmp ) ;


	std::vector<double> finalX ;
	std::vector<NRTSF> FinalNRTrans ;                      // non-rigid transform
	std::vector<Profile2D> finalProfiles ;				   // final profiles
	std::vector<Profile3D> finalCtrPoints ;				   // final profiles
	std::vector<ON_NurbsCurve> finalNurbs ;


	std::vector<Profile2D> profiles2d ;                    // 2d profiles
	std::vector<cutPlane> cutplanes ;                      //profiling planes of each profile
	std::vector<Point3f> centers ; 


	std::vector<int> templatesIds ;                  
	std::vector<ON_NurbsCurve> Tmps ;                      // templates

	std::vector<int2> TScope ;
	std::vector< std::vector<int> > sharpIds ;
	std::vector< std::vector<int> > cvIdx ;
	bool firstIsTip ;
	bool lastIsTip ;

	std::vector<double> tmpConsWeights ;                 

	// for solving the obj piecewise
	int LeftProfId ;
	std::vector<double> lastX ;
	bool oneIterationSolverCalled ;


	// for loop 
	bool isLoop ;

	bool wdataScaled ;

	//std::vector<Point2f> transformProfile( int proId,  ST stf) ;

	double E_data( std::vector<NRTSF> &x , int pid );
	double E_rigid( std::vector<NRTSF> &x,  int pid );
	double E_smooth( std::vector<NRTSF> &x ,  int pid, int changedCvID = -1 ) ;
	double E_featrueStroke( std::vector<NRTSF> &x, int profId  ) ;
	double E_template_constrain( std::vector<NRTSF> &x, int pid ) ;
	double E_tip_constrain( std::vector<NRTSF> &x, int pid ) ;
	double E_conti(  NRTSF nrt1,NRTSF nrt2 , NRTSF nrt1l,NRTSF nrt2r , int Tid1,int Tid2 ) ;
	double objFunc( std::vector<double> x ) ;
	double gradient( const std::vector<double> &x , int xdiff_id, double oriProfEnergy ) ;
	std::vector<double> oneIterationSolver(  std::vector<double> x0   ) ;
	std::vector<double> morphing( ) ;
	std::vector<double> init_afterFeatureStroke( ) ;

	std::vector<ON_NurbsCurve> consolidateFinalNurbs( std::vector<ON_NurbsCurve> finalNurbs ) ;
	std::vector<ON_NurbsCurve> smoothFinalTrajectory( std::vector<ON_NurbsCurve> finalNurbs ) ;

	std::vector<std::vector<double>>  WeightOfCV ;

	std::vector<double>  solveOneSegment(  std::vector<double> x0, int leftProfId  ) ;

	// for sharp feature stroke
	void setFeatureStrokeCVId( int CVId, Curve2D featureStroke_w2d,Curve3D featureStroke_l3d, int2 sharpScope ) ;
	std::vector<int>  featureStrokeCVIds ;
	std::vector<MVMatrix>  mvmtrices ;
	std::vector<Curve2D>  featureStroke_w2d ;
	std::vector<Curve3D>  featureStroke_l3d ;
	std::vector<int2> sharpStrokeScope ;

	// for detecting correspondence
	std::vector<int2> corres_cues ;
	std::vector< std::vector<int2> > corres ;
	void detectCorrespondence () ;

	// for debug
	void drawFinalCVcorres() ;

	//added by huajie
	void initMofit(std::vector<NRTSF> &x);
	void adjustCtrlPoints(std::vector<NRTSF> &x);
	double xmult(Point2f p1,Point2f p2,Point2f p0);
	void MiniDisWith2pointss(Point2f p,Point2f q,int n);
	void MiniDisWithpointss(Point2f pi,int n);

public:

	double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data ) ;


	~completionSolver() {}
};


#else

class completionSolver ;

#endif//COMPLETIONSOLVER_COMPLETIONSOLVER_V2_H