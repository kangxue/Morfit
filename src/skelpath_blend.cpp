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

#include "skelpath.h"
#include "GlobalFunction.h"
#include "reconstructionUtility.h"
#include "gl/GL.h"
#include "appstate.h"


void skelpath::removeInnerPoints(){

	//std::vector<std::vector<Point3f>> subpoints ;
	//std::vector<std::vector<Point3f>> pwBrches ;
	//std::vector<std::vector<cutPlane>> pwBrchesCtpl ;
	//std::vector<int2> pwBrchesCmpnt ;
	//std::vector<int> pwBrcheJid ;
	// resultProf3d
	std::cout<<"skelpath::removeInnerPoints"<<std::endl;


	assert(subpoints.size()==pwBrches.size() && pwBrches.size()==pwBrchesCmpnt.size()&&pwBrchesCmpnt.size()==pwBrcheJid.size() &&pwBrcheJid.size()==resultProf3d.size() ) ;


	for( int i=0; i<pwBrcheJid.size(); ++i ){

		std::vector<bool> toRemove( subpoints[i].size(), false) ;

		for( int j=0;j<pwBrcheJid.size(); ++j ){

			if( i==j||pwBrcheJid[i]!=pwBrcheJid[j])
				continue ;

			bool needUpdate = true ;
			if( changedPwbranId.size()!=0  ){
				needUpdate = false ;
				for( int id=0; id<changedPwbranId.size(); ++id)
					if(changedPwbranId[id] == i || changedPwbranId[id] == j )
						needUpdate = true ;
			}

			if( !needUpdate )
				continue ;

			std::cout <<"###"<<std::endl;

			std::vector<Point3f> pwbrchJ_us;  // upsampled pwbranch,j
			for( int nid=0; nid<pwBrches[j].size()-1; ++nid )
				for( int n=0; n<10; ++n )
					pwbrchJ_us.push_back( pwBrches[j][nid] + (pwBrches[j][nid+1]-pwBrches[j][nid])*n/10.0 ) ;
			pwbrchJ_us.push_back( pwBrches[j].back() ) ;


			// remove points of cylinder i that located in cylinder j
			std::vector<std::vector<int>> nstId ;


			ReconstructorUtility::computeFlannKNN( subpoints[j], subpoints[i], nstId, 1 ) ;



			for( int pid=0; pid<subpoints[i].size(); ++pid ){

				//if( toRemove[pid] )
				//	continue ;

				int qid = nstId[pid][0] ;
				Point3f p = subpoints[i][pid] ;
				Point3f q = subpoints[j][qid] ;

				int bJid_p = ReconstructorUtility::NearstPoint( pwbrchJ_us, p) ;
				int bJid_q = ReconstructorUtility::NearstPoint( pwbrchJ_us, q) ;

				if(  (pwbrchJ_us[bJid_q]-q).Norm() - (pwbrchJ_us[bJid_p]-p).Norm() > 0.002  )
					toRemove[pid] = true ;


			}
		}

		for( int pid=0; pid<subpoints[i].size(); ++pid )
			if( toRemove[pid] ) {
				toRemove.erase( toRemove.begin() + pid ) ;
				subpoints[i].erase( subpoints[i].begin() + pid ) ;
				pid--;
			}

	}

}

std::vector<Point3f> skelpath::blendBranchBetweenJoints(){

	//std::vector<Point3f> points_;
	//for( int i=0; i<subpoints.size(); ++i )
	//	points_.insert(points_.begin(), subpoints[i].begin(), subpoints[i].end() ) ;
	//return points_ ;


	int Jnum = 1 ;
	for( int i=0; i<pwBrcheJid.size();  ++i)
		if( pwBrcheJid[i]+1>Jnum ) Jnum = pwBrcheJid[i]+1 ;

	std::vector<std::vector<Point3f>> jointPoints(Jnum) ;// points set for each joint
	for( int i=0; i<pwBrcheJid.size();  ++i)
		jointPoints[ pwBrcheJid[i]].insert( jointPoints[ pwBrcheJid[i] ].begin(), subpoints[i].begin(),subpoints[i].end() ) ;

	std::cout << "mark 1"<<std::endl;


	// --------------- generate points for blending
	std::vector<std::vector<Point3f>> PointSet1 ;
	std::vector<std::vector<Point3f>> PointSet2 ;
	//std::vector<int2> Jids ;
	std::vector<int> brchIds ;
	for( int i=0; i<pwBrcheJid.size(); ++i ){
		for( int j=i+1;j<pwBrcheJid.size(); ++j ){

			if( pwBrcheJid[i] == pwBrcheJid[j] )
				continue ;

			// common branch of 2 piecewise branch
			std::vector<int> cmnBrch ; // common branch
			if( pwBrchesCmpnt[i].x == pwBrchesCmpnt[j].x || pwBrchesCmpnt[i].x == pwBrchesCmpnt[j].y )
				cmnBrch.push_back( pwBrchesCmpnt[i].x );
			if( pwBrchesCmpnt[i].y == pwBrchesCmpnt[j].x || pwBrchesCmpnt[i].y == pwBrchesCmpnt[j].y )
				cmnBrch.push_back( pwBrchesCmpnt[i].y );

			if( cmnBrch.empty() )
				continue ;

			for( int cbid=0; cbid<cmnBrch.size(); ++cbid ){

				int bid = cmnBrch[cbid] ;

				// branch bid
				std::vector<Point3f> brch ;
				for( int id=brchScopes[bid].x; id<=brchScopes[bid].y; ++id )
					brch.push_back( smoothedBrchPts[bid][id] ) ;

				if( brch.size() <= 3 )
					continue ;

				std::vector<Point3f> brch_us ;  // upsampled branch bid
				for( int nid=0; nid<brch.size()-1; ++nid )
					for( int n=0; n<10; ++n )
						brch_us.push_back( brch[nid] + (brch[nid+1]-brch[nid])*n/10.0 ) ;
				brch_us.push_back( brch.back() ) ;

				// start and end index of piecewise branch segment that overlap with branch bid
				int pwSgmId1_i = ReconstructorUtility::NearstPoint( pwBrches[i], brch[0] ) ;
				int	pwSgmId2_i = ReconstructorUtility::NearstPoint( pwBrches[i], brch.back() ) ;
				int pwSgmId1_j = ReconstructorUtility::NearstPoint( pwBrches[j], brch[0] ) ;
				int pwSgmId2_j = ReconstructorUtility::NearstPoint( pwBrches[j], brch.back() ) ;
				if( pwSgmId1_i> pwSgmId2_i ) { int tmp = pwSgmId2_i; pwSgmId2_i=pwSgmId1_i ; pwSgmId1_i=tmp; }
				if( pwSgmId1_j> pwSgmId2_j ) { int tmp = pwSgmId2_j; pwSgmId2_j=pwSgmId1_j ; pwSgmId1_j=tmp; }
				if( pwSgmId1_i < brch.size()/2 )  pwSgmId1_i = 0;
				if( pwSgmId2_i > pwBrches[i].size() -  brch.size()/2 ) pwSgmId2_i = pwBrches[i].size()-1 ;
				if( pwSgmId1_j < brch.size()/2 )  pwSgmId1_j = 0;
				if( pwSgmId2_j > pwBrches[j].size() -  brch.size()/2 ) pwSgmId2_j = pwBrches[j].size()-1 ;


				int brchEnd_i = 0;  // which end of branch bid locates in pwBranch i 
				int brchEnd_j ;     // which end of branch bid locates in pwBranch j
				if( MyMin( (brch_us[0] - pwBrches[i][0]).Norm(), (brch_us[0] - pwBrches[i].back() ).Norm() ) 
					<  MyMin( (brch_us.back() - pwBrches[i][0]).Norm(), (brch_us.back() - pwBrches[i].back() ).Norm() )  )
					brchEnd_i = brch_us.size()-1 ;
				brchEnd_j = brch_us.size()-1 - brchEnd_i ;


				std::vector<Point3f> subpi, subpj ;  // points subset associated with branch bid

				// compute points subset associated with branch bid
				// and delete them from subpoints arrat
				std::vector<std::vector<int>> nstId ;
				GlobalFun::computeKNN( pwBrches[i], subpoints[i], nstId, 1 ) ;
				for( int pid=0; pid<subpoints[i].size(); ++pid )
					if( nstId[pid][0] >= pwSgmId1_i && nstId[pid][0]<= pwSgmId2_i){
						subpi.push_back( subpoints[i][pid] ) ;
						subpoints[i][pid] = Point3f(0,0,0) ;
					}
					for( int pid=0; pid<subpoints[i].size(); ++pid )
						if( subpoints[i][pid] == Point3f(0,0,0) ){
							subpoints[i].erase( subpoints[i].begin() + pid ) ; 
							pid--;
						}

						GlobalFun::computeKNN( pwBrches[j], subpoints[j], nstId, 1 ) ;
						for( int pid=0; pid<subpoints[j].size(); ++pid )
							if( nstId[pid][0] >= pwSgmId1_j && nstId[pid][0]<= pwSgmId2_j){
								subpj.push_back( subpoints[j][pid] ) ;
								subpoints[j][pid] = Point3f(0,0,0) ;
							}
							for( int pid=0; pid<subpoints[j].size(); ++pid )
								if( subpoints[j][pid] == Point3f(0,0,0) ){
									subpoints[j].erase( subpoints[j].begin() + pid ) ; 
									pid--;
								}


								bool found = false ;
								for( int id=0; id<PointSet1.size(); ++id )
									if( brchIds[id]==bid ){
										found = true ;
										if( brchEnd_i < brchEnd_j ){
											PointSet1[id].insert( PointSet1[id].begin(), subpi.begin(), subpi.end() ) ;
											PointSet2[id].insert( PointSet2[id].begin(), subpj.begin(), subpj.end() ) ;
										} else {
											PointSet1[id].insert( PointSet1[id].begin(), subpj.begin(), subpj.end() ) ;
											PointSet2[id].insert( PointSet2[id].begin(), subpi.begin(), subpi.end() ) ;
										}
									}

									if( !found ){
										if( brchEnd_i < brchEnd_j ){
											PointSet1.push_back( subpi ) ;
											PointSet2.push_back( subpj ) ;
										}else{
											PointSet1.push_back( subpj ) ;
											PointSet2.push_back( subpi ) ;
										}
										//Jids.push_back( int2( pwBrcheJid[i],pwBrcheJid[j] ) ) ;
										brchIds.push_back( bid ) ;
									}

			}		
		}
	}

	// ---------------- blend
	for( int id=0; id<PointSet1.size(); ++id ){

		std::vector<Point3f> subpi = PointSet1[id] ;
		std::vector<Point3f> subpj = PointSet2[id] ;
		int bid = brchIds[id] ;

		// branch bid
		std::vector<Point3f> brch ;
		for( int id=brchScopes[bid].x; id<=brchScopes[bid].y; ++id )
			brch.push_back( smoothedBrchPts[bid][id] ) ;
		if( brch.size() <= 3 )
			continue ;
		std::vector<Point3f> brch_us ;  // upsampled branch bid
		for( int nid=0; nid<brch.size()-1; ++nid )
			for( int n=0; n<10; ++n )
				brch_us.push_back( brch[nid] + (brch[nid+1]-brch[nid])*n/10.0 ) ;
		brch_us.push_back( brch.back() ) ;

		// weights
		std::vector<double> wpi,wpj ;

		std::vector<std::vector<int>> nstId ;
		GlobalFun::computeKNN( brch_us, subpi, nstId, 1 ) ;
		for( int pid=0; pid<subpi.size(); ++pid )
			wpi.push_back( 1.0 - nstId[pid][0] /(double)(brch_us.size()-1) ) ;

		GlobalFun::computeKNN( brch_us, subpj, nstId, 1 ) ;
		for( int pid=0; pid<subpj.size(); ++pid )
			wpj.push_back( nstId[pid][0]/(double)(brch_us.size()-1) ) ;

		// blend
		std::vector<Point3f> newsubpi(subpi.size() ) ;
		std::vector<Point3f> newsubpj(subpj.size() ) ;

		GlobalFun::computeKNN( subpj, subpi, nstId, 1 ) ;
		for( int pid=0; pid<subpi.size(); ++pid ){
			int qid = nstId[pid][0] ;
			newsubpi[pid] = (subpi[pid] * wpi[pid] + subpj[qid]*wpj[qid])/(wpi[pid]+wpj[qid]) ;
		}
		GlobalFun::computeKNN( subpi, subpj, nstId, 1 ) ;
		for( int pid=0; pid<subpj.size(); ++pid ){
			int qid = nstId[pid][0] ;
			newsubpj[pid] = (subpj[pid] * wpj[pid] + subpi[qid]*wpi[qid])/(wpj[pid]+wpi[qid]) ;
		}
		PointSet1[id] = newsubpi ;
		PointSet2[id] = newsubpj ;
	}


	//extern std::vector<Point3f > globalPointsToDraw1  ;
	//globalPointsToDraw1.clear();
	//for( int i=0; i<PointSet1.size(); ++i )
	//	globalPointsToDraw1.insert(globalPointsToDraw1.begin(), PointSet1[i].begin(), PointSet1[i].end() ) ;
	//for( int i=0; i<PointSet2.size(); ++i )
	//	globalPointsToDraw1.insert(globalPointsToDraw1.begin(), PointSet2[i].begin(), PointSet2[i].end() ) ;


	std::vector<Point3f> points;
	for( int i=0; i<subpoints.size(); ++i )
		points.insert(points.begin(), subpoints[i].begin(), subpoints[i].end() ) ;


	for( int i=0; i<PointSet1.size(); ++i )
		points.insert(points.begin(), PointSet1[i].begin(), PointSet1[i].end() ) ;
	for( int i=0; i<PointSet2.size(); ++i )
		points.insert(points.begin(), PointSet2[i].begin(), PointSet2[i].end() ) ;

	return points ;
}