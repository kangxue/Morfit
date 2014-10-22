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


#ifndef _skeleton_mul_h_ 
#define _skeleton_mul_h_

#include "Types.h"
#include<vector>
#include "CMesh.h"

#include <lemon/list_graph.h>

#include <pcl/surface/on_nurbs/fitting_curve_2d_pdm.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_tdm.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_sdm.h>
#include <pcl/surface/on_nurbs/triangulation.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
//#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/surface/3rdparty/opennurbs/opennurbs_nurbscurve.h>

#include "smtf.h"
#include "nrigidtransform.h"
#include "Types.h"
#include "triangleList.h"

class completionSolver ;

class skeleton_mul {

public:
	skeleton_mul(){  iditifier = clock(); }
	skeleton_mul( std::vector<std::vector<Point3f>> & brchPts, bool dsample ) ;
	skeleton_mul( completionSolver solver ) ;

	unsigned iditifier ;

	std::vector<Point3f> reconstruct() ;
	
	// data structure utilities
	void convertVerticesListToGraph() ;
	void convertGraphToVerticeList() ;
	void convertVeticesListToBranchPoints() ;
	void convertBranchPointsToVeticesList(  double disToMerge ) ;

	//calculate and get the member variable "smoothedBrchPts"
	void smoothEachBranch() ;

	//get the joints , buffer joints in the first parameter "joints"
	void getJoints( std::vector<Point3f> &joints, std::vector<std::vector<int>> &ajBrchId ) ;

	//get the tips ,buffer tips in the "endpoints"  
	void getTips( std::vector<int2> &endpoints) ;

	//calculate and get the member variable "cutplanes"
	void calculateCutplanes() ;

	//  branch points representation
	std::vector<std::vector<Point3f>>  brchPts ;

	// graph representation 
	lemon::ListGraph  *skeletonGraph ;
	lemon::ListGraph::NodeMap<Point3f> *graphNodeP3f;

	// vertices list representation 
	std::vector<Point3f> nodes ;  
	std::vector<std::vector<int>> branches ;

	// degree of each points
	std::vector<std::vector<int>>  brchNodeDegree ;

	// smoothed branch points, cannot used for calculus purpose 
	std::vector<std::vector<Point3f>>  smoothedBrchPts ;

	// profiling planes 
	std::vector<std::vector<cutPlane>> cutplanes ;

	void normalize( Box3f box ) ;/**/
	
	// for loop
	bool brch0IsLoop ;

};

#include "completionSolver_v2.h" 

#else

class sweeper ;

#endif