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



#include "skelpath.h"
#include "bspline.h"
#include "Types.h"
#include "completionSolver_v2.h"

#include <fstream>
extern std::ofstream logout ;

#include <Eigen/Core>
#include <Eigen/Eigen>
#include <algorithm>
#include <GL/GL.h>
//#include <cv.h>
//#include <highgui.h>
#include <string> ;
#include "reconstructionUtility.h"
#include "reconstructorPara.h"
#include "declarations.h"
#include "reconstructorPara.h" 
#include "appstate.h"

std::vector<Point3f> skelpath::reconstruct() {
	// return point cloud
	std::vector<Point3f> points ;


	if( solvers.size() == 0 ){
		if( tempalteIdss.size() == 0)
			return points ;
		solvers = generateJointSetting() ;
	}

	jointSettingNum = solvers.size() ;
	cutplanes_bp.clear() ;
	cutplanes_bp.resize( solvers.size() ) ;
	branches_bp.clear() ;
	branches_bp.resize( solvers.size() ) ;
	subpoints.clear() ;
	resultProf3d.clear() ;
	resultProf2d.clear() ;
	for( int iter =0; iter<solvers.size(); ++iter ){
		std::cout<<"iter="<<iter<<std::endl;

		completionSolver &comSolver = solvers[iter];

		char tmp[10] ;
		std::string xfilename = std::string("x") + std::string( itoa(iter, tmp,10) )+".txt" ;

		std::ifstream xifs(xfilename) ;
		std::vector<double> inx,outx ;
		//double tt ;
		//while( xifs >> tt )
		//	inx.push_back(tt);

		outx = comSolver.solve( inx ) ;

		std::ofstream xofs(xfilename) ;
		for( int i=0;i<outx.size(); ++i )
			xofs << outx[i] << std::endl;
		xofs.close() ;
		xifs.close() ;


		std::cout<<"mark--------1" <<std::endl;
		std::vector<cutPlane> cutpls = comSolver.cutplanes ;
		std::vector<Point3f> brch = comSolver.centers ;



		std::vector<std::vector<Point3f>> profiles3d ;
		std::vector<std::vector<Point2f>> profiles2d ;
		for( int i=0; i<comSolver.finalProfiles.size(); ++i ){

			Profile3D p3d =  ReconstructorUtility::convert2dProfileTo3d( comSolver.finalProfiles[i], cutpls[i], brch[i] ) ;
			profiles3d.push_back( p3d ) ;
			profiles2d.push_back( comSolver.finalProfiles[i] ) ;

			// store cutplane and plnormal of branch pairs for debuging
			cutplanes_bp[iter].push_back(cutpls[i] );
			branches_bp[iter].push_back(brch[i] );

		}

		resultProf3d.push_back( profiles3d ) ;
		resultProf2d.push_back( profiles2d ) ;


		std::cout<<"mark--------2" <<std::endl;

		std::vector<Point3f> subpoint = ReconstructorUtility::upsampleResult( profiles3d, int2(0, brch.size()-1 ), comSolver.firstIsTip, comSolver.lastIsTip ) ;
		points.insert(points.begin(), subpoint.begin(), subpoint.end() ) ;

		subpoints.push_back( subpoint ) ; 	// for blend

		std::cout<<"mark--------3" <<std::endl;

		// store control points
		extern std::vector<Point3f > globalPointsToDraw1  ; 
		for( int i=0; i<comSolver.finalCtrPoints.size(); ++i )
			for( int j=0; j<comSolver.finalCtrPoints[i].size(); ++j )
				globalPointsToDraw1.push_back( comSolver.finalCtrPoints[i][j] ) ;


		// store sharp features
		extern std::vector<Point3f > globalPointsToDraw2  ;

		for( int i=0; i<comSolver.finalCtrPoints.size(); ++i )
			for( int k=0; k<comSolver.sharpIds[0].size(); ++k)
				globalPointsToDraw2.push_back( comSolver.finalCtrPoints[i][comSolver.sharpIds[0][k]] ) ;

	}


	generateMesh() ;

	return points ;
}


bool GlobalFun::mycompfunc (std::pair<double,int> a,std::pair<double,int> b) ;

std::vector<completionSolver> skelpath::generateJointSetting(){


	std::vector<completionSolver> solvers ;
	// if there is only one branch
	if( smoothedBrchPts.size() == 1 ){
		extern std::vector<Curve3D> glonalOriSkel ;
		bool firstIsTip=( ReconstructorUtility::isTip(*this, smoothedBrchPts[0][0]  )  || ReconstructorUtility::isMidNode( glonalOriSkel, smoothedBrchPts[0][0]) );
		bool lastIsTip =( ReconstructorUtility::isTip(*this, smoothedBrchPts[0].back()  ) || ReconstructorUtility::isMidNode(glonalOriSkel, smoothedBrchPts[0].back()));
		if( firstMustBeTip && !firstMustBeNotTip ) firstIsTip = true ;
		if( !firstMustBeTip && firstMustBeNotTip ) firstIsTip = false ;
		if( lastMustBeTip && !lastMustBeNotTip ) lastIsTip = true ;
		if( !lastMustBeTip && lastMustBeNotTip ) lastIsTip = false ;
		if( brch0IsLoop )
			firstIsTip = lastIsTip = false ;
		solvers.push_back(completionSolver(*this, 0, firstIsTip, lastIsTip, corres_cues)  ) ;
		return solvers ;
	}else{

		std::cout << "multi-branch skeleton!" <<std::endl;
		system("pause") ;

	}

	return solvers ;
}