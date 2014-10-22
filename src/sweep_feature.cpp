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

void sweeper::startFeatureStroke(Point2f p) {

	//start feature stroke

	ptsSelected = std::vector<bool>(points.size(), false) ;
	std::cout<<" sweeper::startFeatureStroke"<<std::endl;
	ongoingFeatureStroke.clear() ;
	ongoingFeatureStroke.push_back(p) ;

}

void sweeper::moveFeatureStroke(Point2f p) {

	//during stroke

	std::cout<<" sweeper::moveFeatureStroke"<<std::endl;

	if( (ongoingFeatureStroke.back() - p).Norm() > 10  )
		ongoingFeatureStroke.push_back(p) ;
	else
		return ;

}

bool sweeper::endFeatureStroke(Point2f p, float *mvmatrix) {

	//end stroke , and call function "generateSetting()"

	std::cout<<" sweeper::endFeatureStroke"<<std::endl;

	if( ongoingFeatureStroke.size()>5){

		if( (ongoingFeatureStroke.back() - p).Norm() > 10  )
			ongoingFeatureStroke.push_back(p) ;

		strokes.push_back( ongoingFeatureStroke ) ;
		MVMatrix mvm ;
		for( int i=0; i<16; ++i )
			mvm.m[i] = mvmatrix[i] ;

		mvmatrices.push_back( mvm ) ;
	}
	else{
		ongoingFeatureStroke.clear() ;
		return false ;
	}

	featureStrokeSkelId = getFeatureStrokeSkelId() ;
	if( featureStrokeSkelId >= 0 )
		performFeatureStroke() ;

	ongoingFeatureStroke.clear() ;
	return true;
}

int sweeper::getFeatureStrokeSkelId() {


	// upsample the stroke
	Curve2D stroke = ongoingFeatureStroke ;
	for( int i=0; i<stroke.size()-1; ++i )
		if( (stroke[i] - stroke[i+1]).Norm() > 2 )
			stroke.insert( stroke.begin()+i+1,stroke[i] +(stroke[i+1] - stroke[i]).Normalize() ) ;


	int SkelId = 0;
	double minDis = 1e10 ;
	for( int sid =0; sid<settings.size(); ++sid){

		Curve2D skel_s2d = ReconstructorUtility::convertWorld2dToScreen2d( ReconstructorUtility::project3dPointsOntoScreen( settings[sid].smoothedBrchPts[0] ) ) ;

		std::vector<std::vector<int>> nstIds ;
		ReconstructorUtility::computeFlannKNN(skel_s2d, stroke, nstIds, 1) ;

		double dis = 0 ;
		for( int i=0; i<stroke.size(); ++i )
			dis += (skel_s2d[ nstIds[i][0] ] - stroke[i] ).Norm() ;

		if( dis < minDis ){
			minDis = dis ;
			SkelId = sid ;
		}
	}

	if( minDis/stroke.size() < 100 ){

		std::cout << "sweeper::getFeatureStrokeSkelId: featureStrokeSkelId = " << SkelId <<std::endl ;
		return SkelId ;
	}else
		return -1;


}

double round(double r)
{
	return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}
void sweeper::performFeatureStroke(){

	std::cout<<" sweeper::performFeatureStroke{"<<std::endl;
	int cvnum = ReconstructorPara::cvNum ;

	skelpath &strokedSkel = settings[featureStrokeSkelId] ;

	// project the stroke to the surface
	std::vector<Profile3D> prof3d = strokedSkel.resultProf3d[0];
	std::vector<Profile2D> prof2d =  strokedSkel.resultProf2d[0];
	std::vector<ON_NurbsCurve> nbs = strokedSkel.solvers[0].finalNurbs ;
	for( int i=0; i<nbs.size(); ++i ){
		prof2d[i] = ReconstructorUtility::discretizeNurbs(nbs[i],500) ;
		prof3d[i] = ReconstructorUtility::convert2dProfileTo3d( prof2d[i], strokedSkel.solvers[0].cutplanes[i], strokedSkel.solvers[0].centers[i] ) ;
	}

	std::vector<Profile2D> prof3d_w2d ;
	std::vector<Profile3D> prof3d_w3d ;
	for( int i=0; i<prof3d.size(); ++i ){
		prof3d_w2d.push_back(  ReconstructorUtility::project3dPointsOntoScreen( prof3d[i]) ) ;
		prof3d_w3d.push_back(  ReconstructorUtility::convertLocal3dToWorld3d( prof3d[i]) ) ;
	}

	Curve2D ogFeaStr_w2d = ReconstructorUtility::convertScreen2dToWorld2d( ongoingFeatureStroke ) ;
	ReconstructorUtility::getUniformCubicBSpline( ogFeaStr_w2d, ogFeaStr_w2d, 100) ;

	std::vector<int> interId(prof3d.size(),-1 ) ;
	for( int i=0; i<prof3d.size(); ++i){

		std::vector<int> candidateNode ;
		for( int j=0; j<prof3d[i].size(); ++j ){
			int n = prof3d[i].size() ;
			for( int k=0; k<ogFeaStr_w2d.size()-1; ++k )
				if( ReconstructorUtility::vectorIntersect( prof3d_w2d[i][j], prof3d_w2d[i][(j+1)%n], ogFeaStr_w2d[k], ogFeaStr_w2d[k+1] )  ){

					// candidate node must locate in front
					int back_count = 0;
					int front_count = 0;
					for( int id=0; id<prof3d_w3d[i].size(); ++id )
						if( prof3d_w3d[i][id].Z() < prof3d_w3d[i][j].Z() )
							back_count++ ;
						else
							front_count++ ;

					if( back_count > front_count ){
						candidateNode.push_back(j) ;
						candidateNode.push_back((j+1)%n) ;
					}

					//break;
				}
		}

		if( candidateNode.size() ){
			int nstNodeId = 0;
			double depth = -1e10 ;
			for( int id=0; id<candidateNode.size(); ++id ){
				if( prof3d_w3d[i][candidateNode[id]].Z() > depth ){ nstNodeId=candidateNode[id]; depth=prof3d_w3d[i][candidateNode[id]].Z() ; }
			}
			interId[i] = nstNodeId ;
		}
	}


	// output featureStroke3d for displaying
	featureStroke3d.clear() ;
	for( int i=0; i<interId.size(); ++i  ){
		if( interId[i]>0 ){
			featureStroke3d.push_back( prof3d[i][interId[i]] ) ;
		}
	}

	// find sharp control point
	int2 sharpScope;
	bool first = true ;
	for( int i=0; i<interId.size(); ++i ) {
		if( first&&interId[i]>=0 ){
			first = false ;
			sharpScope.x = i ;
		}

		if( !first && interId[i]>=0  )
			sharpScope.y = i ;
	}

	std::vector<Curve3D> cvpts ;
	for( int i=0 ; i<nbs.size(); ++i )
		cvpts.push_back( ReconstructorUtility::convert2dProfileTo3d( ReconstructorUtility::getCvPts(nbs[i]),  strokedSkel.solvers[0].cutplanes[i], strokedSkel.solvers[0].centers[i] ) ) ;

	std::vector<int> sharpCvId(prof3d.size(),-1 ) ;
	for( int i=0; i<interId.size(); ++i ) {
		if( interId[i] < 0 )
			continue ;

		//std::vector<double> distances ;
		//for( int j=0; j<cvpts[i].size(); ++j )
		//	distances.push_back( ( featureStroke3d[ReconstructorUtility::NearstPoint( featureStroke3d,cvpts[i][j])] - cvpts[i][j]).Norm() ) ;
		//sharpCvId[i] = ReconstructorUtility::getMinId(distances) ;

		sharpCvId[i] = ( (int)round( ( interId[i] / (double)prof3d[i].size() ) / (1.0/cvnum) + 1.0) +  cvnum )%cvnum;
	}

	// compute nstCVid
	std::vector<int> cvidDiff( ReconstructorPara::cvNum, 0) ;
	for( int i=0; i<sharpCvId.size(); ++i )
		if( sharpCvId[i] >= 0 )
			for( int j=0; j<cvidDiff.size(); ++j )
				cvidDiff[j] += std::min( abs(j-sharpCvId[i]), abs( ReconstructorPara::cvNum - abs(j-sharpCvId[i]) ) ) ;
	int nstCVID = ReconstructorUtility::getMinId(cvidDiff) ;

	strokedSkel.solvers[0].setFeatureStrokeCVId( nstCVID, ogFeaStr_w2d, featureStroke3d, sharpScope ) ;


	for( int i=0; i<sharpCvId.size(); ++i )
		std::cout<< sharpCvId[i] <<" " ;
	std::cout <<"\nnstCVID = " << nstCVID<<std::endl;

	// reconstruction
	strokedSkel.reconstruct() ;
	meshes[featureStrokeSkelId] = settings[featureStrokeSkelId].resultMesh[0] ;
	
	ongoingSkel = strokedSkel ;
	ongoingSkel=skelpath(strokedSkel.smoothedBrchPts,true);
	mergeMesh() ;

	std::cout <<"}"<<std::endl;

}