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


namespace completionSolverPara{
	//int winside = 7; 
	extern int winside_mscurvature; 

	extern int smooth_winside;

	extern int pwCurveSampleNum ;

	extern double w_data ;
	extern double w_rigid ;
	extern double w_smooth ;
	extern double w_tmpconstraint ;
	extern double w_tipconstraint ;
	extern double w_featureconstraint ;
	extern double one ;

	extern int edit_winside  ;

	//double w_smooth_Scal_end = 1e-3;
	//double w_smooth_Scal_mid = 1;
	extern double w_continuity ;
	extern double w_cMs;

	//double init_regular = 0.001 ;
	//double regular_incspeed = 2.0 ;

	extern double tStep ;
	extern double AStep  ;
	extern double bStep  ;

	extern double t_ub  ;
	extern double t_lb  ;
	extern double A_ub  ;
	extern double A_lb  ;
	extern double b_ub  ;
	extern double b_lb  ;

	extern bool diagnosticMessage ;
	extern bool smoothTrajactory ;

	extern int trajectoryNurbsDegree ;

}