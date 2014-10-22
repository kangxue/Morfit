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


#ifndef _appstate_h_
#define  _appstate_h_
#include <string>
#include "Types.h"
#include "triangleList.h"

#include "skeleton_mul.h"

class AppState{

	// AppState are some switches ,determining some function on or off
public:
	AppState(){
		//   default states
		displayCutRadii = false ;
		displayPoints = true ;
		displayControlPoints = false ;
		display3dJoints = true ;
		displayJoint = true ;
		displayJointBrchId = 0 ;
		displayTrajectory = false ;
		displaySeedProfile = true; 
		displaySkeleton = true;   
		displayMesh = true;      
		displayResultProfile = false ;
		displaySkeletonConfidence = false ;
		cullFace = false ;
		displayMeshVertOnly = false ;
		displayMeshDiffColor = false ;
		displayGeckoBoundary = false ;

		vedioMode = false ;

		cameraRotating = false ;
		cameraRotateDegree = 0 ;

		sweepAnimation = false ;
		sweepingId = 0 ;
		sweepingAnimation_meshNum = 0 ;

	}
	
	std::string plyfilename ;
	std::string path ;
	
	bool displayCutRadii ; 


	bool displayControlPoints ; 

	bool display3dJoints ;

	bool displayJoint ;

	int displayJointBrchId ;
	 
	bool displayTrajectory;             // 1
	bool displaySeedProfile;            // 2
	bool displaySkeleton;               // 3
	bool displayMesh;                   // 4
	bool displayResultProfile;          // 5
	bool displaySkeletonConfidence;     // 6
	bool displayPoints;					// 7
	bool cullFace;					    // 8
	bool displayMeshVertOnly;		    // 9

	bool displaySettings;               // 0
	std::vector<skelpath> skelProfs ;

	bool displayGeckoBoundary ;           

	bool displayMeshDiffColor;            
	std::vector<double3> MeshDiffColors ;
	triangleList resultMesh ;

	bool vedioMode ;

	bool cameraRotating ;
	double cameraRotateDegree ;

	bool sweepAnimation ;
	int sweepingId ;
	int sweepingAnimation_meshNum ;

};


extern AppState appstate ;

#endif // APPSTATE_APPSTATE_H