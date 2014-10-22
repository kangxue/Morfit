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



#include "completionSolver_v2.h"
#include "reconstructionUtility.h"



void completionSolver::drawFinalCVcorres() {

	int cvnum = ReconstructorPara::cvNum ;

	std::vector<Curve3D> cvarrays;

	for( int cvi=0; cvi<cvnum;++cvi ){
		
		Curve3D curve ;
		for( int i=0; i<finalCtrPoints.size();++i )
			curve.push_back( finalCtrPoints[i][cvi] ) ;
		
		cvarrays.push_back( curve ) ;
	}

	for( int i=0; i<cvarrays.size(); ++i ){
		//GlobalFun::draw3dCurves( cvarrays, GLColor(0,0,0), isLoop, false, ReconstructorPara::trajDisplaySize ) ;
	//	GlobalFun::draw3dCurves( cvarrays, GLColor(0,0,0), isLoop, false, ReconstructorPara::trajDisplaySize ) ;
		
		glColor3f( 0,0,0) ;
		glLineWidth(2) ;
		for( int i=0; i<cvarrays.size(); ++i ){
			glBegin(GL_LINE_STRIP) ;
			for( int j=0; j<cvarrays[i].size(); ++j  )
				glVertex3f( cvarrays[i][j].X(),  cvarrays[i][j].Y(),  cvarrays[i][j].Z() ) ;
			glEnd() ;
		}
	}

}