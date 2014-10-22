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


#ifndef _ReconstructorPara_
#define _ReconstructorPara_


namespace ReconstructorPara{

	extern double LSWforTsolver ;
	//extern double skeletonSampleLength ;
	extern double MMA_FTOL ;

	//extern int profileNumEachBranch ;

	extern double branchSampleStep;

	extern int profileOnLeftScreenWidth ;

	extern int profileDownSampleNum ;

	extern int NurbsDegree ;

	extern int cvNum ;

	extern int maxTemplateNum ;

	extern int weight_sharpCV;

	extern int minimumBranchPointNum ;

	extern 	double upsampleStep ;

	extern int finalProfilesSampleNum ;

	extern int winside_skel_deform  ;

	extern double skelDisplaySize  ;
	extern double nurbsDisplaySize  ;
	extern double trajDisplaySize  ;
}

#endif