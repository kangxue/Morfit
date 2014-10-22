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


#include "GLArea.h"
extern GLArea *glarea_global_ptr ;


#include "skeleton_mul.h"
#include "sweeper.h"

#include "appstate.h"

std::vector<int> ptsubsetBetween2Line( std::vector<Point2f> points2d, Point2f p0, Point2f p1, Point2f t0){

	t0.Normalize() ;
	std::vector<int> res ;
	for( int i=0; i<points2d.size(); ++i ){
		Point2f p0p1 = p1-p0 ;
		if(  p0p1 * (points2d[i]-p1) <= 0 &&  !ReconstructorUtility::vectorIntersect(p0-t0*100, p0+t0*100, points2d[i], p1 ) )
			res.push_back( i ) ;
	}
	return res ;
}


std::vector<Curve3D> glonalOriSkel ;
sweeper::sweeper( std::vector<Point3f> &ptcld, std::vector<Point3f> &nmls,std::vector<Curve3D> &skel ) {
	points = ptcld ;
	normals = nmls ;
	glonalOriSkel = impSkel = skel ;

	std::vector<std::vector<int>> ajBrchId ;
	oriSkel=skeleton_mul(skel,true);
	oriSkel.getJoints(allJoints,ajBrchId);
	ptsSelected = std::vector<bool>(points.size(), false) ;

	//ReconstructorUtility::computeFlannKNN( points, points, kneibs, 50 * points.size() / 5000.0 ) ;
	
	//segment
	std::vector<Point3f> skeletalPts ;
	std::vector<int> bid ;
	for( int i=0; i<impSkel.size(); ++i ){
		for( int j=0; j<impSkel[i].size(); ++j ){
			skeletalPts.push_back( impSkel[i][j] ) ;
			bid.push_back(i) ;
		}
	}

	std::vector<std::vector<int>> nsts ;
	ReconstructorUtility::computeFlannKNN(skeletalPts, points,nsts,1 ) ;

	for( int i=0; i<points.size(); ++i )
		segLable.push_back( bid[ nsts[i][0] ] ) ;

	

	ctrlSkelId = -1;
	ctrlProfId = -1;
	ctrlCVId = -1;
	ctrlCVSelectionType = 1 ;
}

sweeper::sweeper( std::vector<Point3f> &ptcld, std::vector<Point3f> &nmls,std::vector<Curve3D> &skel, std::vector<Point3f> &allTheJoints ,std::vector< std::vector<Point3f> > &branches) {
	/*allJoints=std::vector<Point3f>(allTheJoints.size(),0);
	for(int jid=0; jid<allTheJoints.size();jid++ ){
		allJoints[jid]=allTheJoints[jid];
	}*/

	nml=nmls;
	oriSkel = skeleton_mul(branches, true) ;
	cout<<"branches.size()"<<branches.size()<<endl;
	allJoints=allTheJoints;
	cout<<"allJoints "<<allJoints.size()<<std::endl;
	points = ptcld ;
	normals = nmls ;
	glonalOriSkel = impSkel = skel ;

	ptsSelected = std::vector<bool>(points.size(), false) ;

	//ReconstructorUtility::computeFlannKNN( points, points, kneibs, 50 * points.size() / 5000.0 ) ;

	//segment
	std::vector<Point3f> skeletalPts ;
	std::vector<int> bid ;
	for( int i=0; i<impSkel.size(); ++i ){
		for( int j=0; j<impSkel[i].size(); ++j ){
			skeletalPts.push_back( impSkel[i][j] ) ;
			bid.push_back(i) ;
		}
	}
	cout<<"skeletalPts.size()"<<skeletalPts.size()<<endl;


	//std::vector<std::vector<int>> nsts ;
	//ReconstructorUtility::computeFlannKNN(skeletalPts, points,nsts,1 ) ;
	//for( int i=0; i<points.size(); ++i )
	//	segLable.push_back( bid[ nsts[i][0] ] ) ;

	//  [8/18/2014 CX]
	std::vector<int> nsts ;
	ReconstructorUtility::MycomputeFlannKNN(skeletalPts, points,nmls,nsts ) ;
	for( int i=0; i<points.size(); ++i )
		segLable.push_back( bid[ nsts[i] ] ) ;

	//--

	cout<<"nmls.size()"<<nmls.size()<<endl;
	ctrlSkelId = -1;
	ctrlProfId = -1;
	ctrlCVId = -1;
	ctrlCVSelectionType = 1 ;
}


void sweeper::startSkelStroke(Point2f p) {
	
	//start skeleton stroke
	ptsSelected = std::vector<bool>(points.size(), false) ;
	std::cout<<" sweeper::startSkelStroke"<<std::endl;
	ongoingSkelStroke.clear() ;
	ongoingSkelStroke.push_back(p) ;

}

void sweeper::moveSkelStroke(Point2f p) {

	//during stroke
	std::cout<<" sweeper::moveSkelStroke"<<std::endl;

	if( (ongoingSkelStroke.back() - p).Norm() > 10  )
		ongoingSkelStroke.push_back(p) ;
	else
		return ;

}

bool sweeper::endSkelStroke(Point2f p, float *mvmatrix) {

	//end stroke , and call function "generateSetting()"
	std::cout<<" sweeper::endSkelStroke"<<std::endl;
	
	if( ongoingSkelStroke.size()>5){

		if( (ongoingSkelStroke.back() - p).Norm() > 10  )
			ongoingSkelStroke.push_back(p) ;

		strokes.push_back( ongoingSkelStroke ) ;
		MVMatrix mvm ;
		for( int i=0; i<16; ++i )
			mvm.m[i] = mvmatrix[i] ;

		mvmatrices.push_back( mvm ) ;
		ongoingSkelStroke.clear() ;
	}
	else{
		ongoingSkelStroke.clear() ;
		return false ;
	}

	
	return generateSetting() ;
}

bool sweeper::generateSetting(){
	
	std::cout<< "sweeper::generateSetting" <<std::endl;

	// find associated impSkel branches
	std::vector<Curve2D>  impSkel_w2d ;
	for( int i=0; i<impSkel.size(); ++i )
		impSkel_w2d.push_back( ReconstructorUtility::project3dPointsOntoScreen(impSkel[i]) );
	Curve2D stroke_w2d = ReconstructorUtility::convertScreen2dToWorld2d( strokes.back() ) ;


	std::vector<Point2f> skeletalPts_w2d ;   // by find nearest skeletal point
	std::vector<int2> idLable ;
	for( int i=0; i<impSkel_w2d.size(); ++i  )
		for( int j=0; j<impSkel_w2d[i].size() ; ++ j){
			skeletalPts_w2d.push_back( impSkel_w2d[i][j] ) ;
			idLable.push_back( int2(i,j) ) ;
		}

	std::vector<int> asscoBrchId ;
	std::vector<std::vector<int>> nsts ;
	ReconstructorUtility::computeFlannKNN(skeletalPts_w2d, stroke_w2d, nsts,1 ) ; 
	for( int i=0; i<stroke_w2d.size(); ++i ){
		int2 BNid = idLable[nsts[i][0]] ;
		if( BNid[1] > impSkel[BNid[0]].size()/5 && BNid[1] < impSkel[BNid[0]].size()*4/5  )
			asscoBrchId.push_back( BNid[0] ) ;
	}
	for( int i=0; i<asscoBrchId.size(); ++i ){
		for( int j=i+1; j<asscoBrchId.size(); ++j )
			if( asscoBrchId[i]==asscoBrchId[j]){ asscoBrchId.erase(asscoBrchId.begin()+j); j--; }
	}

	if( asscoBrchId.size() == 0 )
		return false;

	// fit a bezier curve as the new branch
	Curve3D ctrpts = impSkel[asscoBrchId[0]];

	if( asscoBrchId.size() > 1  ){  // determin direction of first branch


		Curve3D &b0 = impSkel[asscoBrchId[0]] ;
		Curve3D &b1 = impSkel[asscoBrchId[1]] ;

		double d0 = std::min( (b1[0] - b0[0]).Norm() ,  (b1.back() - b0[0] ).Norm() );
		double d1 = std::min( (b1[0] - b0.back()).Norm() ,  (b1.back() - b0.back()).Norm() );

		if( d0 < d1 ){
			ctrpts.clear() ;
			for( int id=b0.size()-1; id>=0; --id )
				ctrpts.push_back(b0[id] ) ;
		}else{

			ctrpts.clear() ;
			for( int id=0; id<b0.size(); ++id )
				ctrpts.push_back(b0[id] ) ;
		}
	}

	for( int i=1; i<asscoBrchId.size(); ++i ){ // determine direction of the rest branch

		Curve3D &bi = impSkel[asscoBrchId[i]] ;

		if( (bi[0] - ctrpts.back()).Norm() <  (bi.back() - ctrpts.back()).Norm() )
			for( int id=0; id<bi.size(); ++id )
				ctrpts.push_back(bi[id] ) ;
		else
			for( int id=bi.size()-1; id>=0; --id )
				ctrpts.push_back(bi[id] ) ;

	}

	Curve3D branch;
	//ReconstructorUtility::getBezier(ctrpts, branch, ctrpts.size() ) ;
	ReconstructorUtility::getNurbs(ctrpts, branch,100,20 ) ;

	if( branches.size() == strokes.size() )
		branches.back() = branch ;
	else
		branches.push_back( branch ) ;


	// detect loop
	bool isLoop ;
	int n=branch.size() ;
	if( (branch[0] - branch.back()).Norm() < 0.05  && (branch[0] - branch[1] ) *(branch[n-1] - branch[n-2] ) < 0  ){
		isLoop = true ;
		std::cout << "skeleton is loop " <<std::endl;
		ReconstructorUtility::getPeriodicNurbs(branch, branch, 100,20 ) ;
	}
	else{
		isLoop = false ;
		std::cout << "skeleton is not loop " <<std::endl;
	}


	// ------------------ prepare the settings --------------
	std::vector<Point3f> subpts ;
	std::vector<int> subptsMapId;  //added by huajie,2014/8/19
	subptsMapId.clear();
	subpts.clear();
	ptsSelected = std::vector<bool>(points.size(), false) ;
	for( int i=0; i<segLable.size(); ++i ){
		bool found = false ;
		for( int j=0; j<asscoBrchId.size(); ++j ) 
			if( asscoBrchId[j] == segLable[i] )
				found = true ;

		if( found ){
			subpts.push_back(points[i]) ;
			ptsSelected[i] = true ;
		}
	}
	
	ongoingSkel = skelpath( std::vector<Curve3D>(1,branch),false ) ;
	//  [8/14/2014 Cao]
	ongoingSkel.brch0IsLoop = isLoop ;
	ongoingSkel.intialSegment( subpts ) ;
	ongoingSkel.calculateCutplanes() ;
	ongoingSkel.calculateProfiles() ;
	std::vector<Point3f> subpointsmp=subpts;
	ongoingSkel.addPointInJoints(points,subpts,ptsSelected,allJoints,nml);
	if(subpts.size() == subpointsmp.size()){
		cout<<"not add in fact!~~~~~~~~~~~~~~~~~"<<endl;
	}else{
		cout<<subpts.size()-subpointsmp.size()<<"  Points added!"<<endl;
		cout<<points.size()<<" points at all!"<<endl;
	}


	ongoingSkel.brch0IsLoop = isLoop ;

	ongoingSkel.intialSegment( subpts) ;

	ongoingSkel.calculateCutplanes() ;
	ongoingSkel.calculateProfiles() ;

	//delete by huajie
	//if( !ReconstructorUtility::isTip(ongoingSkel, ctrpts[0]) || ReconstructorUtility::isMidNode(impSkel, ctrpts[0]) ){
	//	ongoingSkel.removeTips( strokes.back()[0] ) ;  
	//	std::cout << "first end is not tip" <<std::endl;
	//}
	//if( !ReconstructorUtility::isTip(ongoingSkel, ctrpts.back()) || ReconstructorUtility::isMidNode(impSkel, ctrpts.back()) ){
	//	ongoingSkel.removeTips( strokes.back().back() ) ;
	//	std::cout << "last end is not tip" <<std::endl;
	//}




	if( !isLoop )
		ongoingSkel.snapTips(points, normals);



	ongoingSkel.intialSegment( subpts ) ;

	ongoingSkel.calculateCutplanes() ;
	ongoingSkel.calculateProfiles() ;
	ongoingSkel.initAfterClaculateProfiles() ;


	return true ;


}



void sweeper::reconLastStroke( ) {

	if( appstate.vedioMode ){
		int brchNum = settings.size() ;
		loadSettings() ;

		if( settings.size() >= brchNum+1 )
			settings.resize( brchNum+1 ) ;



		// write iditifier
		meshLable.resize(brchNum+1) ;
		meshes.resize(brchNum+1) ;

		mergeMesh() ;
		
		return ;
	}

	ongoingSkel.reconstruct() ;

	bool found = false;
	for( int i=0; i<meshes.size(); ++i ){
		if( meshLable[i] == ongoingSkel.iditifier){
			found = true ;
			meshes[i] = ongoingSkel.resultMesh[0] ;
		}
	}
	

	if( !found ){
		meshes.push_back( ongoingSkel.resultMesh[0] );
		meshLable.push_back( ongoingSkel.iditifier ) ;
		settings.push_back( ongoingSkel ) ;
	}

	
	mergeMesh() ;
}





void sweeper::segment(std::vector<Point2f> segStroke, bool shiftPressed, bool ctrPressed, bool leftbutton ) {

	std::cout << "sweeper::segment called" <<std::endl;

	if( leftbutton ){

		std::cout << "leftbutton"<<std::endl;

		std::vector<Point2f> stroke_w2d = ReconstructorUtility::convertScreen2dToWorld2d(  segStroke ) ;
		std::vector<Point2f> points_w2d = ReconstructorUtility::project3dPointsOntoScreen(  points ) ;

		for( int i=0; i<points_w2d.size(); ++i ){
			if( GlobalFun::PtInPolygon( points_w2d[i], stroke_w2d )  ){

				if( shiftPressed )
					ptsSelected[i] = false ;
				if( ctrPressed )
					ptsSelected[i] = true ;
			}
		}

	} else {

		std::cout << "rightButton"<<std::endl;

		std::vector<Point2f> stroke_w2d = ReconstructorUtility::convertScreen2dToWorld2d(  segStroke ) ;
		std::vector<Point2f> brch_w2d = ReconstructorUtility::project3dPointsOntoScreen(  ongoingSkel.smoothedBrchPts[0] ) ;

		if( GlobalFun::PtInPolygon( brch_w2d[0], stroke_w2d ) ){
			if( ctrPressed ){ ongoingSkel.firstMustBeTip = true ; ongoingSkel.firstMustBeNotTip = false ;  std::cout << "make first node be tip"<<std::endl;}
			if( shiftPressed ){	ongoingSkel.firstMustBeTip = false ; ongoingSkel.firstMustBeNotTip = true ;  std::cout << "make first node not tip"<<std::endl; }

		}
			
		if( GlobalFun::PtInPolygon( brch_w2d.back(), stroke_w2d ) ){
			if( ctrPressed ){ ongoingSkel.lastMustBeTip = true ; ongoingSkel.lastMustBeNotTip = false ; std::cout << "make last node be tip"<<std::endl; }
			if( shiftPressed ){	ongoingSkel.lastMustBeTip = false ; ongoingSkel.lastMustBeNotTip = true ; std::cout << "make last node not tip"<<std::endl;}
		}

	}

	// write the segment result to ongoingSkel
	std::vector<Point3f> subpts ;
	subpts.clear();
	for( int i=0; i<ptsSelected.size(); ++i )
		if( ptsSelected[i] ){
			subpts.push_back( points[i] ) ;
		}

	ongoingSkel.intialSegment(subpts) ;


	//delete by huajie
	//if( !ReconstructorUtility::isTip(ongoingSkel, ongoingSkel.smoothedBrchPts[0][0]) || ReconstructorUtility::isMidNode(impSkel, ongoingSkel.smoothedBrchPts[0][0]) ){
	//	ongoingSkel.removeTips( strokes.back()[0] ) ;
	//	std::cout << "first end is not tip" <<std::endl;
	//}
	//if( !ReconstructorUtility::isTip(ongoingSkel, ongoingSkel.smoothedBrchPts[0].back()) || ReconstructorUtility::isMidNode(impSkel, ongoingSkel.smoothedBrchPts[0].back()) ){
	//	ongoingSkel.removeTips( strokes.back().back() ) ;
	//	std::cout << "last end is not tip" <<std::endl;
	//}
	if( !ongoingSkel.brch0IsLoop )
		ongoingSkel.snapTips(points, normals);

	ongoingSkel.intialSegment( subpts ) ;
	ongoingSkel.calculateCutplanes() ;
	ongoingSkel.calculateProfiles() ;
	ongoingSkel.initAfterClaculateProfiles() ;

}

//
//void sweeper::estBranchAndNurbs() {
//
//	// ------------------ estimate 3d branch --------------
//
//	std::vector<Point3f> subpoints ;
//	for( int i=0; i<points.size(); ++i) 
//		if( ptsSelected[i] )  subpoints.push_back( points[i] ) ;
//	std::vector<Point2f> stroke_w2d = ReconstructorUtility::convertScreen2dToWorld2d(  strokes.back() ) ;
//	std::vector<Point2f> points_w2d = ReconstructorUtility::project3dPointsOntoScreen(  subpoints, mvmatrices.back().m  ) ;
//	std::vector<Point3f> points_w3d = ReconstructorUtility::convertLocal3dToWorld3d(  subpoints, mvmatrices.back().m  ) ;
//
//	std::vector<std::vector<int>> nsts ;
//	ReconstructorUtility::computeFlannKNN( stroke_w2d, points_w2d, nsts,1) ;
//
//	std::vector<Point3f> stroke_w3d ;
//	for( int i=0; i<stroke_w2d.size(); ++i ){
//		float min = 1e10 ; 
//		float max = -1e10 ;
//
//		int count = 0 ;
//		for( int j=0; j<points_w3d.size(); ++j ){
//
//			if( nsts[j][0] == i ){
//				if( points_w3d[j].Z() > max ) max = points_w3d[j].Z() ;
//				if( points_w3d[j].Z() < min ) min = points_w3d[j].Z()  ;
//				count++ ;
//			}
//
//		}
//
//		if( count > 5 )
//			stroke_w3d.push_back(  Point3f( stroke_w2d[i][0], stroke_w2d[i][1], -0.1 )/(-0.1)  * (min+max)/2 )  ;
//	}
//	//for( int i=0; i<stroke_w3d.size(); ++i )
//	//	if( fabs(stroke_w3d[i].Z()) > 10 ){ stroke_w3d.erase( stroke_w3d.begin()+i ) ; i-- ;}
//
//
//	Curve3D brch = ReconstructorUtility::convertWorld3dToLocal3d(stroke_w3d, mvmatrices.back().m  ) ;
//	ReconstructorUtility::getBezier(brch, brch, 100) ;
//
//	if( branches.size() == stroke_w2d.size() )
//		branches.back() = brch ;
//	else
//		branches.push_back( brch ) ;
//
//
//	// ------------------ prepare the settings --------------
//	ongoingSkel = skeleton_mul( std::vector<Curve3D>(1,brch),false ) ;
//	ongoingSkel.intialSegment( subpoints ) ;
//
//	ongoingSkel.removeTips( strokes.back().back() ) ;
//	ongoingSkel.snapTips();
//
//	ongoingSkel.calculateCutplanes() ;
//	ongoingSkel.calculateProfiles() ;
//	ongoingSkel.initAfterClaculateProfiles() ;
//
//
//	if( settings.size() == stroke_w2d.size() )
//		settings.back() = ongoingSkel ;
//	else
//		settings.push_back( ongoingSkel ) ;
//
//
//
//}




void sweeper::mergeMesh() {

	// merge all meshes

	mergedMesh.destroy() ;

	for( int i=0; i<meshes.size(); ++i )
		mergedMesh.addMesh( meshes[i]) ;

	mergedMesh.sortDepth() ;
}