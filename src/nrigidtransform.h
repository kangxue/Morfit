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



#ifndef _nrigidtransform_h_
#define _nrigidtransform_h_

#include "CMesh.h"


namespace completionSolverPara{

	extern double angleGraStep ;
	extern double scalingGraStep ;
	extern double transGraStep  ;
}

typedef class nrigidTransform{
public:

	nrigidTransform( int n){ 
		t.clear();
		for( int i=0; i<n;++i){
			t.push_back(Point2f(0,0)) ;
		}

		A = 0 ;
		b = Point2f(0,0) ;
	}

	nrigidTransform(){}

	std::vector<Point2f> t ;  // translations of each control point
	double A ; // rotation
	Point2f b ; // translation


}NRTSF;



#endif