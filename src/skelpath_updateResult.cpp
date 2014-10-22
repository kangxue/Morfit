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


#include "GlobalFunction.h"
#include "skelpath.h"
#include "bspline.h"


#include <fstream>
extern std::ofstream logout ;

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <algorithm>    // std::sort

#include <GL/GL.h>

//#include <cv.h>
//#include <highgui.h>

#include <string> ;

#include "reconstructionUtility.h"
#include "reconstructorPara.h"

#include "declarations.h"
#include "reconstructorPara.h" 



std::vector<Point3f> skelpath::updateResult(int newProfId, ON_NurbsCurve newProf ){

	std::cout << "skelpath::updateResult" <<std::endl;

	std::vector<Point3f> points_;


	solvers[0].resolve(newProfId, newProf ) ;


	subpointsLast = subpoints ;
	points.clear();
	subpoints.clear() ;
	resultProf3d.clear() ;
	resultProf2d.clear() ;
	for( int iter =0; iter<solvers.size(); ++iter ){
		completionSolver &comSolver = solvers[iter];


		std::vector<cutPlane> cutpls = comSolver.cutplanes ;
		std::vector<Point3f> brch = comSolver.centers ;

		std::vector<std::vector<Point3f>> profiles3d ;
		std::vector<std::vector<Point2f>> profiles2d ;
		for( int i=0; i<comSolver.finalProfiles.size(); ++i ){

			Profile3D p3d =  ReconstructorUtility::convert2dProfileTo3d( comSolver.finalProfiles[i], cutpls[i], brch[i] ) ;
			profiles3d.push_back( p3d ) ;
			profiles2d.push_back( comSolver.finalProfiles[i] ) ;

		}

		resultProf3d.push_back( profiles3d ) ;
		resultProf2d.push_back( profiles2d ) ;

		std::vector<Point3f> subpoint = ReconstructorUtility::upsampleResult( profiles3d, int2(0, brch.size()-1 ), comSolver.firstIsTip, comSolver.lastIsTip ) ;
		points.insert(points.begin(), subpoint.begin(), subpoint.end() ) ;

		subpoints.push_back( subpoint ) ; 	// for blend

	}




	generateMesh() ;

	return points_ ;
}
