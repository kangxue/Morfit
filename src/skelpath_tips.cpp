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
#include "skelpath.h"
#include <GL/GL.h>
#include "GlobalFunction.h"
#include "reconstructionUtility.h"
#include "appstate.h"
#include "Types.h"
//draw
std::vector<Point3f> dsts;
std::vector<Point3f> fraAroundTip;
std::vector<Point3f> centerOfGravitys;
std::vector<Point3f> centerTip;
std::vector<Point3f> allNmls;
std::vector<Point3f> edges;
////
extern std::vector<Curve3D> glonalOriSkel ;
extern std::vector<Point3f> BranchEnds ;


void skelpath::snapTips(std::vector<Point3f> ptcld , std::vector<Point3f> nmls) {

	//return ;
	edges.clear();
	std::vector<int> nsts ;
	std::vector<Point3f> skeletalPts ;
	std::vector<Point3f> ptsAroundTip;
	std::vector<Point3f> nmlAroundTip;
	std::vector<Point3f> normalAroundTip;
	allNmls.clear();
	nsts.clear();
	skeletalPts.clear();
	ptsAroundTip.clear();
	nmlAroundTip.clear();

	for( int i=0; i<glonalOriSkel.size(); ++i ){
		for( int j=0; j<glonalOriSkel[i].size(); ++j ){
			skeletalPts.push_back( glonalOriSkel[i][j] ) ;
		}
	}

	ReconstructorUtility::MycomputeFlannKNN( skeletalPts, ptcld, nmls, nsts );
	for( int eid = 0; eid<BranchEnds.size(); ++eid ){
		for(int sid=0; sid<skeletalPts.size(); sid++){
			if(BranchEnds[eid] == skeletalPts[sid]){
				for(int pid=0; pid<ptcld.size(); pid++){
					if(nsts[pid] == sid){
						ptsAroundTip.push_back(ptcld[pid]);
						nmlAroundTip.push_back(nmls[pid]);
					}
				}
			}
		}
	}

	// tips are initialized in intialSegment (){}
	// modify tips
	if( smoothedBrchPts.size() == 1 ) {

		if( firstMustBeTip && !firstMustBeNotTip  ){
			tips.push_back( int2(0, 0 )) ;
		}
		if( lastMustBeTip && !lastMustBeNotTip  ){
			tips.push_back( int2(0, smoothedBrchPts[0].size()-1 )) ;
		}

		if( !firstMustBeTip && firstMustBeNotTip  ){
			for( int i=0; i<tips.size(); ++i)
				if( tips[i].y == 0 ){
					tips.erase( tips.begin() + i ) ; --i ;
				}
		}
		if( !lastMustBeTip && lastMustBeNotTip  ){
			for( int i=0; i<tips.size(); ++i)
				if( tips[i].y == smoothedBrchPts[0].size()-1 ){
					tips.erase( tips.begin() + i ) ; --i ;
				}
		}

	}
	// find dst position of tips
	dsts.clear();
	centerOfGravitys.clear();
	fraAroundTip.clear();
	centerTip.clear();

	std::vector< Point3f > fra;
	std::vector< Point3f > farts;
	std::vector< Point3f > tipNormal;
	std::vector< Point3f > tipTangent;
	farts.clear();
	tipNormal.clear();
	tipTangent.clear();
	for( int i=0; i<tips.size(); ++i ){
		fra.clear();
		normalAroundTip.clear();
		int bid = tips[i].x ;
		int nid = tips[i].y ;
		//std::cout<<i<<"::  "<<"bid=="<<bid<<" nid=="<<nid<<endl;
		Point3f p = smoothedBrchPts[bid][nid] ;
		Point3f normal ;

		if( nid == 0 )
			normal =  smoothedBrchPts[bid][0] -  smoothedBrchPts[bid][1] ;
		else 
			normal =  smoothedBrchPts[bid][nid] -  smoothedBrchPts[bid][nid-1] ;
		normal.Normalize() ;

		Point3f centerOfGravity = Point3f(0.0,0.0,0.0);
		for( int pid=0; pid<ptsAroundTip.size(); ++pid ){

			Point3f v = ptsAroundTip[pid] - p ;
			double dis = v*normal ;
			double alpha1 = ReconstructorUtility::vectorIntersectionAngle( v, normal ) ;
			//double alpha2 = ReconstructorUtility::vectorIntersectionAngle( v, nmls[fragMapId[bid][pid]] ) ;

			if(alpha1 < 3.1415926/2 ) {
				bool flag = true;
				for(int j=0; j < smoothedBrchPts[bid].size(); ++j){
					if(j!=nid && (ptsAroundTip[pid] - smoothedBrchPts[bid][j]).Norm() <= v.Norm()){
						flag = false; break;
					}
				}
				if(flag ){
					fra.push_back(ptsAroundTip[pid]);
					normalAroundTip.push_back(nmlAroundTip[pid]);
				}
			}
		}

		int num = 0;
		for( int fid=0; fid<fra.size(); fid++){
			centerOfGravity += fra[fid];
			num++;
		}
		if(num == 0){
			double dis = 1000;
			Point3f nearFra;
			for(int fid=0; fid<fragments[bid].size(); fid++){
				if( dis > (fragments[bid][fid] - p).Norm() ){
					dis = (fragments[bid][fid] - p).Norm();
					nearFra = fragments[bid][fid];
				}
			}
			dis = 1000;
			int n;
			for(int sbid=0; sbid<smoothedBrchPts[bid].size(); sbid++){
				if( dis > (nearFra - smoothedBrchPts[bid][sbid]).Norm() ){
					dis = (nearFra - smoothedBrchPts[bid][sbid]).Norm();
					n = sbid;
				}
			}
			if(nid == 0){
				for(int sbid=nid; sbid<n; sbid++)
					smoothedBrchPts[bid][sbid] = smoothedBrchPts[bid][n];
			}else{
				for(int sbid=nid; sbid>n; sbid--)
					smoothedBrchPts[bid][sbid] = smoothedBrchPts[bid][n];
			}

			farts.push_back( Point3f(0,0,0) );
			tipTangent.push_back(normal);
			tipNormal.push_back(Point3f(0,0,0));
			continue;
		}
		centerOfGravity /= num;
		centerOfGravitys.push_back(centerOfGravity);
		Point3f totalNormal;
		Point3f farthest ;
		double maxdis = -1000 ;

		if( normalAroundTip.size() > 4 )	
			totalNormal = ReconstructorUtility::ComputeEigenvector( normalAroundTip );
		else
			totalNormal = normalAroundTip[0] ;

		if( ReconstructorUtility::vectorIntersectionAngle( totalNormal, normal ) > 3.1415926/2 )
			totalNormal *= -1;
		totalNormal.Normalize();

		tipTangent.push_back(normal);
		tipNormal.push_back(totalNormal);
		for( int pid=0; pid<fra.size(); ++pid ){
			fraAroundTip.push_back(fra[pid]);
			allNmls.push_back(normalAroundTip[pid]);
			Point3f v = fra[pid] - p ;
			double dis = v*totalNormal ;
			if(  dis > maxdis  &&  sqrtf(v*v-dis*dis) < ReconstructorPara::branchSampleStep){
				maxdis = dis ;
				farthest = fra[pid] ;
			}
		}

		if( maxdis == -1000 ){

			for( int pid=0; pid<fra.size(); ++pid ){
				Point3f v = fra[pid] - p ;
				double dis = v.Norm() ;
				if(  dis > maxdis ){
					maxdis = dis ;
					farthest = fra[pid] ;
				}
			}

			float len = fabs((farthest - p) * totalNormal);
			if( len > 1.5*fabs((farthest - p)*normal))
				farthest = p;
			else{
				farthest = p + totalNormal*len;
				totalNormal *= 100*len;
			}
			fraAroundTip.push_back( p );
			allNmls.push_back( totalNormal );
		}
		farts.push_back( farthest );
		dsts.push_back( farthest );

	}

	// dst position defined by user
	for( int i=0; i<tipSeed.size(); ++i ){
		int nstTipId = 0 ;
		double mindis = 1e10;
		for( int j=0; j<tips.size(); ++j ){
			int bid = tips[j].x ;
			int nid = tips[j].y ;
			Point3f p = smoothedBrchPts[bid][nid] ;
			if( (p-tipSeed[i]).Norm() < mindis ){
				mindis = (p-tipSeed[i]).Norm() ;
				nstTipId = j ;
			}
		}
		if( mindis < 0.1)
			dsts[nstTipId] = tipSeed[i] ;
	}

	// move tips
	for( int i=0; i<tips.size(); ++i ){

		int bid = tips[i].x ;
		int nid = tips[i].y ;
		if( farts[i].Norm() < 0.00001 ) continue;
		Point3f trans = farts[i] - smoothedBrchPts[bid][nid] ;

		if((smoothedBrchPts[bid][nid] - farts[i]).Norm() < ReconstructorPara::branchSampleStep*2){
			//if(nid == 0) smoothedBrchPts[bid].insert( smoothedBrchPts[bid].begin(), farts[i] );
			for(int j=0; trans.Norm()>0.00001; j++){
				smoothedBrchPts[bid][abs(nid - j)] += trans;
				trans = trans / 2;
			}
			continue;
		}
		std::vector<Point3f> hermitePts = ReconstructorUtility::getHermite( smoothedBrchPts[bid][nid], farts[i], tipTangent[i], tipNormal[i], ReconstructorPara::branchSampleStep );
		std::cout<<"hermitePts.size=="<<hermitePts.size()<<"nid=="<<nid<<endl;
		if(nid == 0){
			for(int hid=0; hid<hermitePts.size(); hid++){
				smoothedBrchPts[bid].insert( smoothedBrchPts[bid].begin(), hermitePts[hid] );
				if( tips.size() > 1 )
					tips[abs(i-1)].y++;
			}
		}else{
			for(int hid=0; hid<hermitePts.size(); hid++){
				smoothedBrchPts[bid].push_back(hermitePts[hid]);
				tips[i].y++;
			}
		}

	}

	brchPts = smoothedBrchPts ;
	smoothEachBranch() ;

}

void skelpath::snapFaceTips( std::vector<Curve3D>  impSkel) {

	// tips are initialized in intialSegment (){}
	if( smoothedBrchPts.size() > 1  ) 
		return ;
	skelpath skeletontmp=skelpath(brchPts, true) ;
	std::vector<int2> faceTips ;
	if( !ReconstructorUtility::isTip(skeletontmp, smoothedBrchPts[0][0] ) &&!ReconstructorUtility::isMidNode(impSkel, smoothedBrchPts[0][0] )  )
		faceTips.push_back(int2(0,0)) ;
	if( !ReconstructorUtility::isTip(skeletontmp, smoothedBrchPts[0].back() )&& !ReconstructorUtility::isMidNode(impSkel, smoothedBrchPts[0].back() ) )
		faceTips.push_back(int2(0,smoothedBrchPts[0].size()-1 )) ;



	// find dst position of tips
	std::vector<Point3f> dsts;
	for( int i=0; i<faceTips.size(); ++i ){

		int bid = faceTips[i].x ;
		int nid = faceTips[i].y ;

		Point3f p = smoothedBrchPts[bid][nid] ;
		Point3f normal ;

		if( nid == 0 )
			normal =  smoothedBrchPts[bid][0] -  smoothedBrchPts[bid][1] ;
		else 
			normal =  smoothedBrchPts[bid][nid] -  smoothedBrchPts[bid][nid-1] ;
		normal.Normalize() ;


		Point3f furthest ;
		double maxdis = -1000 ;
		for( int pid=0; pid<fragments[bid].size(); ++pid ){

			Point3f v = fragments[bid][pid] - p ;
			double dis = v*normal ;

			if(  dis > maxdis  &&  sqrtf(v*v-dis*dis) < 0.1 ){
				maxdis = (fragments[bid][pid] - p)*normal ;
				furthest = fragments[bid][pid] ;
			}

		}

		if( maxdis == -1000 ){
			for( int pid=0; pid<fragments[bid].size(); ++pid ){
				Point3f v = fragments[bid][pid] - p ;
				double dis = v*normal ;
				if(  dis > maxdis ){
					maxdis = (fragments[bid][pid] - p)*normal ;
					furthest = fragments[bid][pid] ;
				}
			}

		}


		dsts.push_back( p + normal * maxdis ) ;

	}



	// move tips
	for( int i=0; i<faceTips.size(); ++i ){

		int bid = faceTips[i].x ;
		int nid = faceTips[i].y ;

		Point3f trans = dsts[i] - smoothedBrchPts[bid][nid] ;


		double bsize = smoothedBrchPts[bid].size()-1 ;

		for( int j=0; j<smoothedBrchPts[bid].size(); ++j )
			smoothedBrchPts[bid][j] += trans * (bsize-abs(nid-j))/bsize ; 
	}

	brchPts = smoothedBrchPts ;
	smoothEachBranch() ;
	intialSegment( points ) ;
	calculateCutplanes();
	calculateProfiles() ;




	faceTips.clear() ;
	if( !ReconstructorUtility::isTip(skeletontmp, smoothedBrchPts[0][0] ) &&!ReconstructorUtility::isMidNode(impSkel, smoothedBrchPts[0][0] )  )
		faceTips.push_back(int2(0,0)) ;
	if( !ReconstructorUtility::isTip(skeletontmp, smoothedBrchPts[0].back() )&& !ReconstructorUtility::isMidNode(impSkel, smoothedBrchPts[0].back() ) )
		faceTips.push_back(int2(0,smoothedBrchPts[0].size()-1 )) ;

	int averageProfSize = 0;
	for( int i=0; i<profiles2d[0].size(); ++i)
		averageProfSize+=profiles2d[0].size() ;
	averageProfSize /= profiles2d[0].size() ;

	// find dst position of tips
	dsts.clear();
	for( int i=0; i<faceTips.size(); ++i ){

		int Bid = faceTips[i].x ;
		int Nid = faceTips[i].y ;

		std::vector<Profile2D> &prof2d = profiles2d[Bid] ;
		while( prof2d[Nid].size() < averageProfSize/2 )
			if( Nid < prof2d.size()/2 ) Nid++ ; 
			else Nid-- ;

			dsts.push_back(smoothedBrchPts[Bid][Nid] ) ;
	}


	// move tips
	for( int i=0; i<faceTips.size(); ++i ){

		int bid = faceTips[i].x ;
		int nid = faceTips[i].y ;

		Point3f trans = dsts[i] - smoothedBrchPts[bid][nid] ;


		double bsize = smoothedBrchPts[bid].size()-1 ;

		for( int j=0; j<smoothedBrchPts[bid].size(); ++j )
			smoothedBrchPts[bid][j] += trans * (bsize-abs(nid-j))/bsize ; 
	}

	brchPts = smoothedBrchPts ;
	smoothEachBranch() ;

}

void skelpath::addTipSeed( Point2f seed){

	// convert seed to world coordinate
	std::vector<Point2f> seedarray(1, seed) ;
	seedarray = ReconstructorUtility::convertScreen2dToWorld2d( seedarray) ;

	std::vector<Point2f> points2d = ReconstructorUtility::project3dPointsOntoScreen(points) ;

	std::vector<std::vector<int>> nstIds ;
	GlobalFun::computeKNN( points2d, seedarray, nstIds, 1 ) ;


	Point3f seed3d = points[nstIds[0][0]] ;

	bool areadyExist = false ;
	for( int i=0; i< tipSeed.size(); ++i ){
		if( (tipSeed[i]- seed3d).Norm() < 0.1 ){
			tipSeed[i] = seed3d ;
			areadyExist = true ;
		}
	}

	if( !areadyExist )
		tipSeed.push_back( seed3d ) ;

}


void skelpath::removeTips( Point2f seed){

	// convert seed to world coordinate
	std::vector<Point2f> seedarray(1, seed) ;
	seedarray = ReconstructorUtility::convertScreen2dToWorld2d( seedarray) ;
	Point2f seed2d = seedarray[0] ;

	// traverse all tips
	int nstTipId = 0 ;
	double minDis = 1e10 ;
	for( int i=0; i<tips.size(); ++i ){
		int bid = tips[i].x ;
		int nid = tips[i].y ;
		Point3f p = smoothedBrchPts[bid][nid] ;

		std::vector<Point3f> parray(1, p) ;
		std::vector<Point2f> parray2d ;
		parray2d = ReconstructorUtility::project3dPointsOntoScreen( parray) ;
		Point2f p2d = parray2d[0] ;

		if( (p2d - seed2d).Norm() < minDis ){
			minDis = (p2d - seed2d).Norm() ;
			nstTipId = i ;
		}

	}

	if( minDis < 0.05 ){

		tips.erase( tips.begin() + nstTipId ) ;
	}


}