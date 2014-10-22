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




#include "skeleton_mul.h"
#include "sweeper.h"


void sweeper::generateCorrespondenceInterface() {
	
	// backup data
	settings_bk = settings ;
	meshes_bk = meshes;
	meshLable_bk = meshLable ;

	// get settings
	int n = ongoingSkel.Nurbsss[0].size() ;
	assert( n>0 ) ;

	std::vector<skelpath> corresSkels ;
	for( int i=0; i<n; ++i){

		int nid0,nid1,newId ;
		if( i==0 ) nid0 = 0;
		else nid0 = ( ongoingSkel.tempalteIdss[0][i-1] + ongoingSkel.tempalteIdss[0][i] )/2 ;
		if( i==n-1 )  nid1 = ongoingSkel.profiles2d[0].size() - 1;
		else nid1 = ( ongoingSkel.tempalteIdss[0][i+1] + ongoingSkel.tempalteIdss[0][i] )/2 - 1 ;

		if( i==0 ) newId = ongoingSkel.tempalteIdss[0][0] ;
		else  newId = ongoingSkel.tempalteIdss[0][i-1] - nid0;

		std::vector<Point3f> brch ;
		for( int i=nid0; i<=nid1; ++i )
			brch.push_back( ongoingSkel.smoothedBrchPts[0][i] ) ;

		std::vector<Profile2D> prof2d ;
		for( int i=nid0; i<=nid1; ++i )
			prof2d.push_back( ongoingSkel.profiles2d[0][i] ) ;

		std::vector<cutPlane> ctpl ;
		for( int i=nid0; i<=nid1; ++i )
			ctpl.push_back( ongoingSkel.cutplanes[0][i] ) ;


		skelpath skel(std::vector<Curve3D>(1,brch),false ) ;
		skel.profiles2d.resize(1);  skel.profiles2d[0] = prof2d ;
		skel.cutplanes.resize(1);  skel.cutplanes[0] = ctpl ;
		skel.smoothedBrchPts.resize(1);  skel.smoothedBrchPts[0] = brch ;
		skel.tempalteIdss.resize(1);  skel.tempalteIdss[0] = std::vector<int>(1,newId) ;
		skel.Nurbsss.resize(1);  skel.Nurbsss[0] = std::vector<ON_NurbsCurve>(1, ongoingSkel.Nurbsss[0][i]) ;

		skel.firstMustBeTip = false ;
		skel.firstMustBeNotTip = true ;
		skel.lastMustBeTip = false ;
		skel.lastMustBeNotTip = true ;

		corresSkels.push_back( skel ) ;
		
	}


	// solve surface
	for( int i=0; i<corresSkels.size(); ++i)
		corresSkels[i].reconstruct() ;

	// write meshes
	settings.clear() ;
	meshes.clear() ;
	meshLable.clear() ;
	for( int i=0; i<corresSkels.size(); ++i){
		settings.push_back( corresSkels[i] );
		meshes.push_back(  corresSkels[i].resultMesh[0] ) ;
		meshLable.push_back( corresSkels[i].iditifier ) ;
	}
	mergeMesh() ;

	// write trajectory
	corresTrajectory.clear() ;

	if( trajSegsToDisplay.size()==0 ){
		for( int cvid = 0; cvid < ReconstructorPara::cvNum; ++cvid){

			for( int i=0; i<corresSkels.size(); ++i ){
				Curve3D t ;
				for( int j=0; j<corresSkels[i].solvers[0].finalCtrPoints.size(); ++j )
					t.push_back( corresSkels[i].solvers[0].finalCtrPoints[j][cvid] ) ;
				corresTrajectory.push_back(t) ;
			}
		}

	}else for( int id=0; id<trajSegsToDisplay.size(); ++id  ){
		int cvid = trajSegsToDisplay[id] ;

		for( int i=0; i<corresSkels.size(); ++i ){
			Curve3D t ;
			for( int j=0; j<corresSkels[i].solvers[0].finalCtrPoints.size(); ++j )
				t.push_back( corresSkels[i].solvers[0].finalCtrPoints[j][cvid] ) ;
			corresTrajectory.push_back(t) ;
		}
			

	}

}



void sweeper::recoverAfterCorrespondenceInterface() {

	settings = settings_bk ;
	meshes = meshes_bk;
	meshLable = meshLable_bk ;
	mergeMesh() ;

}
