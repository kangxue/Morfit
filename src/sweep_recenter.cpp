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


#include <QMessageBox>

#include "skeleton_mul.h"
#include "sweeper.h"
	

void sweeper::recenterizeSkel( skeleton_mul &skel ) {

	if( settings.size() == 0 )
		return ;



	// get skeletal Points
	std::vector<Point3f> skeletalPoints ;
	for( int i=0; i<settings.size(); ++i )
		for( int j=0; j<settings[i].resultProf3d[0].size();++j ){
			skeletalPoints.push_back(  GlobalFun::centerOfPoints( settings[i].resultProf3d[0][j] ) ) ;
		}

		// get new branches by finding nearest skeletal points
	std::vector<std::vector<Point3f>> oldBrches = skel.smoothedBrchPts ;
	std::vector<std::vector<Point3f>> newBrches;
	for( int i=0; i<oldBrches.size(); ++i ){
		Curve3D brch ;
		for( int j=0; j<oldBrches[i].size(); ++j )
			brch.push_back( skeletalPoints[ReconstructorUtility::NearstPoint(skeletalPoints, oldBrches[i][j] )]  ) ;
		newBrches.push_back( brch ) ;
	}

	// merge joints
	std::vector<Point3f> joints;
	std::vector<std::vector<int>> ajBrchId ;
	skel.getJoints(joints, ajBrchId)  ;
	for( int i=0; i<joints.size(); ++i ){

		Point3f sum(0,0,0) ;
		for( int j=0; j<ajBrchId[i].size(); ++j ){
			int bid = ajBrchId[i][j] ;
			if( ( newBrches[bid][0]-joints[i] ).Norm() < ( newBrches[bid].back() - joints[i] ).Norm()  )
				sum += newBrches[bid][0] ;
			else
				sum += newBrches[bid].back() ;
		}

		Point3f Jnt = sum / ajBrchId[i].size() ;

		for( int j=0; j<ajBrchId[i].size(); ++j ){
			int bid = ajBrchId[i][j] ;
			if( ( newBrches[bid][0]-joints[i] ).Norm() < ( newBrches[bid].back() - joints[i] ).Norm()  )
				newBrches[bid][0]  = Jnt;
			else
				newBrches[bid].back() = Jnt ;
		}

	}

	// smooth with bezier
	for( int i=0; i<newBrches.size(); ++i  ){
		Curve3D brch ;
		ReconstructorUtility::getBezier( newBrches[i], brch,50 ) ;
		brch[0] = newBrches[i][0] ;
		brch.back() = newBrches[i].back() ;

		newBrches[i] = brch ;

	}


	//check validity
	double sumDis = 0;
	int count = 0;
	for( int i=0; i<skel.smoothedBrchPts.size(); ++i )
		for( int j=0; j<skel.smoothedBrchPts[i].size(); ++j ){
			count++ ;
			sumDis += (skel.smoothedBrchPts[i][j]-newBrches[i][j]).Norm() ;
		}

	if(sumDis/count > 0.05 ){
		QMessageBox msgBox;
		msgBox.setText("Please do recenter only after all branches are reconstructed!");
		msgBox.exec();

		return;
	}

	skel.smoothedBrchPts = newBrches ;
	skel.brchPts = newBrches ;
	
	skel.convertBranchPointsToVeticesList( 0.0002 ) ;
	skel.convertVerticesListToGraph() ;

	



}