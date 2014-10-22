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




namespace ReconstructorPara{

	double LSWforTsolver = 1 ;
	//double skeletonSampleLength = 0.04;
	double MMA_FTOL = 1.0e-2 ;

	//int profileNumEachBranch = 30;

	double branchSampleStep = 0.02;

	int profileOnLeftScreenWidth = 300 ;

	int profileDownSampleNum = 50 ;

	int NurbsDegree = 3 ;

	int cvNum = 10 ;

	int maxTemplateNum = 6 ;

	int weight_sharpCV = 10 ;


	int minimumBranchPointNum = 2;

	double upsampleStep = 0.005 ;

	int finalProfilesSampleNum = 200 ;

	int winside_skel_deform = 10 ;


	double skelDisplaySize = 3.0 ;
	double nurbsDisplaySize = 2.0 ;
	double trajDisplaySize = 0.5 ;
}


namespace completionSolverPara{
	//int winside = 7; 
	int winside_mscurvature = 7; 

	int smooth_winside = 11 ;

	int pwCurveSampleNum = 50 ;

	double w_data = 1000;
	double w_rigid = 1000;
	double w_smooth = 1;
	double w_tmpconstraint = 0 ;
	double w_tipconstraint = 1e5 ;
	double w_featureconstraint = 0 ;
	double one = 1.0 ;


	int edit_winside = 5 ;



	//double w_smooth_Scal_end = 1e-3;
	//double w_smooth_Scal_mid = 1;
	double w_continuity ;
	double w_cMs =  1;

	//double init_regular = 0.001 ;
	//double regular_incspeed = 2.0 ;

	double tStep = 0.004 ;
	double AStep = 0.05 ;
	double bStep = 0.004 ;

	double t_ub = 0.5 ;
	double t_lb = -0.5 ;
	double A_ub = 3.14 ;
	double A_lb = -3.14 ;
	double b_ub = 0.1 ;
	double b_lb = -0.1 ;

	bool diagnosticMessage = true ;

	bool smoothTrajactory = false ;

	int trajectoryNurbsDegree = 40;

}