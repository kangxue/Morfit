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


// ---

#include "GlobalFunction.h"
#include "skeleton_mul.h"
#include "bspline.h"
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

void skeleton_mul::normalize( Box3f box ) {


	std::cout<<" skeleton_mul::normalize"<<std::endl;
	std::cout << "bbox = "<<box.min.X() << " " <<box.min.Y() << " "<<box.min.Z() << "\n"
		<<box.max.X() << " "<<box.max.Y() << " "<<box.max.Z() << std::endl;

	float max_x = abs((box.min - box.max).X());
	float max_y = abs((box.min - box.max).Y());
	float max_z = abs((box.min - box.max).Z());
	float max_length = max_x > max_y ? max_x : max_y;
	max_length = max_length > max_z ? max_length : max_z;

	for(int i = 0; i < nodes.size(); i++){

		Point3f& p = nodes[i];

		p -= box.min;
		p /= max_length;

		p = (p - Point3f( 0.5*max_x/max_length,  0.5*max_y/max_length,  0.5*max_z/max_length));
	}

	convertVerticesListToGraph() ;
	convertVeticesListToBranchPoints() ;
	smoothEachBranch() ;
	calculateCutplanes() ;
}