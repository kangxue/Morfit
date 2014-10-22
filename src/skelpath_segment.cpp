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
#include "Types.h"
extern AppState appstate ;
std::vector<int2> globalBrchNodeInJoint ;
void skelpath::intialSegment( std::vector<Point3f> pts ) {

	// find tips
	getTips( tips) ;

	refineSeg_curves.clear() ;

	points = pts ;

	std::cout<<" skelpath::intialSegment(){}" <<std::endl;
	std::cout << "branchesNum = "<<smoothedBrchPts.size()<<std::endl;


	ptMap2BrchNode.clear() ;
	if( !points.size() || !smoothedBrchPts.size()|| !smoothedBrchPts[0].size() )
		return ;

	std::vector<Point3f> branchPoints ;
	std::vector<int2> index ;
	for( int i=0; i<smoothedBrchPts.size(); ++i ){
		for( int j=0; j<smoothedBrchPts[i].size(); ++j ){
			branchPoints.push_back(smoothedBrchPts[i][j] ) ;
			index.push_back( int2(i,j) ) ;
		}
	}

	std::vector<std::vector<int>> nearestNeibs ;
	GlobalFun::computeKNN( branchPoints,points, nearestNeibs, 1 ) ;


	for( int i=0; i<nearestNeibs.size(); ++i )
		ptMap2BrchNode.push_back( index[nearestNeibs[i][0]] ) ;


	fragments.clear() ;  	fragments_nid.clear() ;
	fragments.resize( smoothedBrchPts.size() ) ;  fragments_nid.resize( smoothedBrchPts.size() ) ;
	for( int i=0; i<points.size(); ++i ){
		fragments[ ptMap2BrchNode[i].x ].push_back( points[i]) ;
		fragments_nid[ ptMap2BrchNode[i].x ].push_back( ptMap2BrchNode[i].y ) ;
	}

	isJoint = std::vector<bool>( pts.size(), false) ;




	if( smoothedBrchPts.size() == 1 ){
		brchScopes = std::vector<int2>(1, int2( 0, smoothedBrchPts[0].size()-1) ) ;
		return ;
	}

	// compute BrchNodeInJoint
	computeBrchNodeInJoint() ;


	// initialize brchScopes
	brchScopes.clear();

	for( int Bid=0; Bid<smoothedBrchPts.size(); ++Bid ){

		// detect scope excluding JOINT
		int2 scope(-1, -1)  ;
		std::vector<int> brchNodeInJoint_bid ;
		for( int i=0; i<brchNodeInJoint.size(); ++i )
			if( brchNodeInJoint[i].x == Bid )
				brchNodeInJoint_bid.push_back( brchNodeInJoint[i].y ) ;

		std::sort( brchNodeInJoint_bid.begin(), brchNodeInJoint_bid.end() ) ;

		if( brchNodeInJoint_bid.size() && brchNodeInJoint_bid[0] == 0 ){
			for( int i=0; i< brchNodeInJoint_bid.size()-1; ++i )
				if(brchNodeInJoint_bid[i]+1!= brchNodeInJoint_bid[i+1])
					scope.x = i+1 ;
			if( scope.x == -1 )
				scope.x = brchNodeInJoint_bid.size() ;
		}else scope.x = 0 ;

		int lastId = smoothedBrchPts[Bid].size()-1 ;
		if( brchNodeInJoint_bid.size() && brchNodeInJoint_bid.back() == lastId ){
			for( int i=brchNodeInJoint_bid.size()-1; i>0 ; --i )
				if(brchNodeInJoint_bid[i]-1!= brchNodeInJoint_bid[i-1])
					scope.y = lastId - ( brchNodeInJoint_bid.size() - i ) ;
			if( scope.y == -1 )
				scope.y = lastId - brchNodeInJoint_bid.size() ;
		}else scope.y = lastId;

		brchScopes.push_back( scope ) ;
	}

}





void skelpath::refineSegment( std::vector<Point2f> curve , bool ctrlPressed) {
	std::cout<<" skelpath::refineSegment(){}" <<std::endl;

	if( !ptMap2BrchNode.size() )
		return ;


	Curve2D curve2d = ReconstructorUtility::convertScreen2dToWorld2d(curve ) ;
	std::vector<Point2f> points2d = ReconstructorUtility::project3dPointsOntoScreen(points);
	CurveArray2D brchPts2d( smoothedBrchPts.size() ) ;
	for( int i=0; i<smoothedBrchPts.size(); ++i )
		brchPts2d[i] =  ReconstructorUtility::project3dPointsOntoScreen(smoothedBrchPts[i]);




	// if the curve is a loop, update brchScopes
	if( (curve[0] - curve.back()).Norm() < 20 ){

		// find brchPts both in joints and in circle
		std::vector<int> foundIds ;
		for( int i=0; i<brchPts2d.size(); ++i ){
			for( int j=0; j<brchPts2d[i].size(); ++j )
				if(GlobalFun::PtInPolygon(brchPts2d[i][j], curve2d) ){
					for( int id=0; id<brchNodeInJoint.size(); ++id )
						if( brchNodeInJoint[id] ==  int2(i,j) )
							foundIds.push_back(id );
				}
		}

		if( foundIds.size() && ctrlPressed ){
			// erase them

			for( int i=0; i<foundIds.size(); ++i )
				brchNodeInJoint[foundIds[i]].x = -1 ;

			for( int i=0; i<brchNodeInJoint.size(); ++i )
				if( brchNodeInJoint[i].x == -1){
					brchNodeInJoint.erase( brchNodeInJoint.begin() + i) ;
					i-- ;
				}

		} else if( foundIds.size() ){

			// put all brchPts in circle into brchNodeInJoint
			for( int i=0; i<brchPts2d.size(); ++i )
				for( int j=0; j<brchPts2d[i].size(); ++j )
					if(GlobalFun::PtInPolygon(brchPts2d[i][j], curve2d) ){
						bool found = false;
						for( int fid=0; fid<foundIds.size(); ++fid )
							if( brchNodeInJoint[ foundIds[fid] ] == int2(i,j) )
								found = true ;
						if( !found )
							brchNodeInJoint.push_back( int2(i,j) ) ;
					}

		}





		brchScopes.clear() ;
		for( int Bid=0; Bid<smoothedBrchPts.size(); ++Bid ){

			// detect scope excluding JOINT
			int2 scope(-1, -1)  ;
			std::vector<int> brchNodeInJoint_bid ;
			for( int i=0; i<brchNodeInJoint.size(); ++i )
				if( brchNodeInJoint[i].x == Bid )
					brchNodeInJoint_bid.push_back( brchNodeInJoint[i].y ) ;

			std::sort( brchNodeInJoint_bid.begin(), brchNodeInJoint_bid.end() ) ;

			if( brchNodeInJoint_bid.size() && brchNodeInJoint_bid[0] == 0 ){
				for( int i=0; i< brchNodeInJoint_bid.size()-1; ++i )
					if(brchNodeInJoint_bid[i]+1!= brchNodeInJoint_bid[i+1])
						scope.x = i+1 ;
				if( scope.x == -1 )
					scope.x = brchNodeInJoint_bid.size() ;
			}else scope.x = 0 ;

			int lastId = smoothedBrchPts[Bid].size()-1 ;
			if( brchNodeInJoint_bid.size() && brchNodeInJoint_bid.back() == lastId ){
				for( int i=brchNodeInJoint_bid.size()-1; i>0 ; --i )
					if(brchNodeInJoint_bid[i]-1!= brchNodeInJoint_bid[i-1])
						scope.y = lastId - ( brchNodeInJoint_bid.size() - i ) ;
				if( scope.y == -1 )
					scope.y = lastId - brchNodeInJoint_bid.size() ;
			}else scope.y = lastId;

			brchScopes.push_back( scope ) ;
		}


	}else{

		//  do refine
		refineSeg_curves.push_back( curve ) ;

		// find points need to be re-segment
		std::vector<int> relatedPoints ;
		for( int i=0; i<points2d.size(); ++i ){

			bool intersect = false ;
			for( int j=0; j<curve2d.size()-1; ++j ){
				if( ReconstructorUtility::vectorIntersect( brchPts2d[ptMap2BrchNode[i].x][ptMap2BrchNode[i].y], points2d[i],curve2d[j], curve2d[j+1] ) ){
					intersect = true ;
					break;
				}
			}

			if( intersect )
				relatedPoints.push_back( i ) ;
		}

		// resegment relatedPoints
		for( int rpid=0; rpid<relatedPoints.size(); ++rpid ){
			int pid = relatedPoints[rpid] ;

			// find nearest branch points
			double minDis = 1.0e10 ;
			int2 minId = int2(-1,-1) ;

			for( int bid=0; bid<brchPts2d.size(); ++bid ){
				for( int bpid=0; bpid<brchPts2d[bid].size(); ++bpid ){

					bool intersect = false ;
					for( int j=0; j<curve2d.size()-1; ++j ){
						if( ReconstructorUtility::vectorIntersect(brchPts2d[bid][bpid], points2d[pid], curve2d[j], curve2d[j+1] ) ){
							intersect = true ; 	break;
						}
					}

					if( intersect ) continue ;

					double dis = (smoothedBrchPts[bid][bpid] - points[pid]).Norm() ;
					if( dis < minDis ){ minDis = dis;  minId = int2(bid, bpid); }

				}
			}

			if( minId.x >= 0 && minId.y>=0)
				ptMap2BrchNode[pid] = minId ;

		}
	}

	// write the result to fragments
	fragments.clear() ;
	fragments.resize( smoothedBrchPts.size() ) ;  fragments_nid.resize( smoothedBrchPts.size() ) ;
	for( int i=0; i<points.size(); ++i ){
		if( !isJoint[i] ){
			fragments[ ptMap2BrchNode[i].x ].push_back( points[i]) ;
			fragments_nid[ ptMap2BrchNode[i].x ].push_back( ptMap2BrchNode[i].y ) ;
		}
	}

}

//  [8/14/2014 XUCAO]

void skelpath::addPointInJoints(std::vector<Point3f> ptscld,std::vector<Point3f> &spts,std::vector<bool> &ptsSelected ,std::vector<Point3f> allJoints,std::vector<Point3f> nml){
	std::cout<<" skelpath::addPointInJoints(){}" <<std::endl;
	std::cout << "branchesNum = "<<smoothedBrchPts.size()<<std::endl;
	cout<<allJoints.size()<<" -joints"<<endl;
	float mindis=0.01;
	isJoint = std::vector<bool>( points.size(), false) ;
	for(int pid=0; pid<points.size() ; pid++){
		for(int jid=0;jid<allJoints.size();jid++){
			if((smoothedBrchPts[ptMap2BrchNode[pid].x][ptMap2BrchNode[pid].y] - allJoints[jid]).Norm()<=mindis){
				isJoint[pid]=true;
			}
		}
	}

	if( !points.size() || !smoothedBrchPts.size()|| !smoothedBrchPts[0].size() )
		return ;


	//calculate the map relationship between the point cloud and the branch points
	std::vector<Point3f> allPoints = ptscld;
	std::vector<Point3f> branchPoints ;
	std::vector<int2> index ;
	std::vector<int2> aptMap2BrchNode ;
	for( int i=0; i<smoothedBrchPts.size(); ++i ){
		for( int j=0; j<smoothedBrchPts[i].size(); ++j ){
			branchPoints.push_back(smoothedBrchPts[i][j] ) ;
			index.push_back( int2(i,j) ) ;
		}
	}

	std::vector<std::vector<int>> nearestNeibs ;
	GlobalFun::computeKNN( branchPoints,allPoints, nearestNeibs, 1 ) ;


	for( int i=0; i<nearestNeibs.size(); ++i )
		aptMap2BrchNode.push_back( index[nearestNeibs[i][0]] ) ;

	//adjust if the joint if a cross
	for(int jid=0;jid<allJoints.size();jid++){

		cout<<"deal with  the joint :"<<jid<<endl;
		double radius = 0;

		int tempbid=0;
		int tempbpid=0;
		bool found=false;
		float mindisTojoint=0.03;
		for(int bid=0;bid<smoothedBrchPts.size() && found==false;bid++){
			for(int bpid=0;bpid<smoothedBrchPts[bid].size() && found == false;bpid++){
				if((smoothedBrchPts[bid][bpid]-allJoints[jid]).Norm()<mindisTojoint){
					tempbid=bid;
					tempbpid=bpid;
					found=true;
				}
			}
		}

		for(int pid=0 ; pid<points.size() ; pid++){
			if(ptMap2BrchNode[pid].x==tempbid && ptMap2BrchNode[pid].y==tempbpid ){
				double dis = (points[pid]-smoothedBrchPts[tempbid][tempbpid]).Norm() ;
				if( dis > radius ){
					radius = dis ;
				}
			}
		}

		bool reachJoint=false;
		bool found1=false;
		bool found2=false;
		double radius1 =0;
		double radius2 =0;
		int2 rap1;
		int2 rap2;

		for(int bid=0;bid<smoothedBrchPts.size() && (found1==false || found2==false);bid++){
			for(int bpid=0;bpid<smoothedBrchPts[bid].size() && (found1==false || found2==false);bpid++){
				if(found1==false && found2 == false && reachJoint ==false){
					double dis1=(smoothedBrchPts[bid][bpid]-allJoints[jid]).Norm();
					if(fabs(dis1-1.5*radius)<=mindis){
						found1=true;
						rap1.x=bid;
						rap1.y=bpid;
					}
				}

				if(reachJoint==true && found2 == false){
					double dis2=(smoothedBrchPts[bid][bpid]-allJoints[jid]).Norm();
					if(fabs(dis2-1.5*radius)<=mindis){
						found2=true;
						rap2.x=bid;
						rap2.y=bpid;
					}
				}
				if((smoothedBrchPts[bid][bpid]-allJoints[jid]).Norm()<mindis){
					tempbid=bid;
					tempbpid=bpid;
					reachJoint=true;
				}
			}
		}
		if(found1==true && found2==true){

			//if yes for a cross
			cout<<"Found 2 points near Joints !!!!!~~~~~"<<endl;
			cout<<allJoints[jid].X()<<allJoints[jid].Y()<<allJoints[jid].Z()<<endl;;
			cout<<smoothedBrchPts[rap1.x][rap1.y].X()<<smoothedBrchPts[rap1.x][rap1.y].Y()<<smoothedBrchPts[rap1.x][rap1.y].Z()<<endl;
			cout<<smoothedBrchPts[rap2.x][rap2.y].X()<<smoothedBrchPts[rap2.x][rap2.y].Y()<<smoothedBrchPts[rap2.x][rap2.y].Z()<<endl;
			for(int pid =0 ;pid < points.size() ; pid ++){
				if(ptMap2BrchNode[pid].x==rap1.x && ptMap2BrchNode[pid].y==rap1.y ){
					double dis = (points[pid]-smoothedBrchPts[rap1.x][rap1.y]).Norm() ;
					if( dis > radius1 ){
						radius1 = dis ;
					}
				}
				if(ptMap2BrchNode[pid].x==rap2.x && ptMap2BrchNode[pid].y==rap2.y ){
					double dis = (points[pid]-smoothedBrchPts[rap2.x][rap2.y]).Norm() ;
					if( dis > radius2 ){
						radius2 = dis ;
					}
				}
			}
			cout<<"radius"<<radius<<endl;
			cout<<"radius1"<<radius1<<endl;
			cout<<"radius2"<<radius2<<endl;

			Point3f p2_p1 =smoothedBrchPts[rap2.x][rap2.y]-smoothedBrchPts[rap1.x][rap1.y];
			Point3f vct =p2_p1/(p2_p1.Norm());
			cout<<p2_p1.Norm()<<" p2_p1.norm()  "<<radius<<endl;


			for(int pid=0; pid < allPoints.size() ; pid++){
				Point3f q_p1 = allPoints[pid]-smoothedBrchPts[rap1.x][rap1.y];

				double t=(q_p1.X()*vct.X()+q_p1.Y()*vct.Y()+q_p1.Z()*vct.Z());
				if(t>p2_p1.Norm() || t<0){
					continue;
				}

				//add points
				double dis_sq=q_p1.SquaredNorm()-t*t;
				double distmp=radius1-t*(radius1-radius2)/p2_p1.Norm();
				if(dis_sq<=1.5*distmp*distmp){
					if(ptsSelected[pid]==false){
						float tmp=((allPoints[pid]-smoothedBrchPts[aptMap2BrchNode[pid].x][aptMap2BrchNode[pid].y])*nml[pid])/
							((allPoints[pid]-smoothedBrchPts[aptMap2BrchNode[pid].x][aptMap2BrchNode[pid].y]).Norm()*nml[pid].Norm());
						if(tmp>0 && tmp<0.7){
							spts.push_back(allPoints[pid]);
							ptsSelected[pid]=true;
						}

					}

				}

			}
		}
	}
}


void skelpath::drawSegmentation( int pointsize ){

	//draw the segment according to data member "fragments"
	GLColor randColor[10] ;

	randColor[0] = GLColor( 1.0,0, 0) ;
	randColor[1] = GLColor( 0,  1,  0) ;
	randColor[2] = GLColor( 0,  0,1.0) ;
	randColor[3] = GLColor( 1.0,0.0,1.0) ;
	randColor[4] = GLColor( 1.0,1.0,0.0) ;
	randColor[5] = GLColor( 0.0,1.0,1.0) ;
	randColor[6] = GLColor( 0.5,0.5,1.0) ;
	randColor[7] = GLColor( 0.5,1.0,0.5) ;
	randColor[8] = GLColor( 0.5,0.5,1.0) ;
	randColor[9] = GLColor( 0.7,0.5,0.7) ;

	glPointSize(pointsize) ;
	for( int i=0; i<fragments.size(); ++i ){

		glColor3f( randColor[i%10].r,randColor[i%10].g,randColor[i%10].b ) ;
		glBegin( GL_POINTS) ;
		for( int j=0; j<fragments[i].size(); ++j )
			glVertex3f(fragments[i][j].X(),fragments[i][j].Y(),fragments[i][j].Z() ) ;
		glEnd() ;

	}

}



void skelpath::drawBrchesInDiffColor( GLColor skelColor){

	//draw branches in different colors
	GLColor randColor[10] ;
	randColor[0] = GLColor( 1.0,0, 0) ;
	randColor[1] = GLColor( 0,  1,  0) ;
	randColor[2] = GLColor( 0,  0,1.0) ;
	randColor[3] = GLColor( 1.0,0.0,1.0) ;
	randColor[4] = GLColor( 1.0,1.0,0.0) ;
	randColor[5] = GLColor( 0.0,1.0,1.0) ;
	randColor[6] = GLColor( 0.5,0.5,1.0) ;
	randColor[7] = GLColor( 0.5,1.0,0.5) ;
	randColor[8] = GLColor( 0.5,0.5,1.0) ;
	randColor[9] = GLColor( 0.7,0.5,0.7) ;

	static std::vector<Point3f> brch ;
	static std::vector<double3> color ;
	static bool  initialized = false ;


	if( profiles2d.size()==1 && appstate.displaySkeletonConfidence ){

		// draw skeleton
		for( int i=0; i<smoothedBrchPts[0].size()-1; ++i ){
			double3 color0 = GlobalFun::scalar2color(0.8 - profile_conf[0][i]) ;
			double3 color1 = GlobalFun::scalar2color(0.8 - profile_conf[0][i+1]) ;
			Point3f p0 = smoothedBrchPts[0][i]  ;
			Point3f p1 = smoothedBrchPts[0][i+1]  ;

			for( int j=0; j<10; ++j ){

				Point3f p = p0 + (p1-p0) * j / 10.0 ;
				double3 c = color0 + (color1 - color0)* j / 10.0 ;
				Point3f pp = p0 + (p1-p0) * (j+1) / 10.0 ;
				//double3 cc = color0 + (color1 - color0)* (j+1) / 10.0 ;

				Curve3D curve ; curve.push_back(p) ;curve.push_back(pp) ;

				GlobalFun::draw3dCurves_nomat( std::vector<Curve3D>(1,curve ) , GLColor(c.x, c.y, c.z), false, false,ReconstructorPara::skelDisplaySize ) ;

				if( !initialized ){
					brch.push_back( p) ;
					color.push_back( c ) ;
				}
			}
		}

		initialized = true ;


		if( 0){

			// draw points
			glDisable(GL_LIGHTING) ;
			std::vector<std::vector<int>> nstId ;
			GlobalFun::computeKNN( brch, points, nstId, 1) ;
			for( int i=0; i< points.size(); ++i){
				double3 c = color[ nstId[i][0] ] ;
				glColor3f(c.x,c.y,c.z) ;
				glPointSize( global_paraMgr.drawer.getDouble("Original Dot Size")) ;
				glBegin(GL_POINTS) ;
				glVertex3f( points[i].X(), points[i].Y(), points[i].Z() ) ;
				glEnd() ;
			}
		}

	} else{
		for( int i=0; i<smoothedBrchPts.size(); ++i ){
			CurveArray3D cvs (1) ;
			cvs[0] = smoothedBrchPts[i] ;
			GlobalFun::draw3dCurves_nomat( cvs,skelColor , brch0IsLoop, false,ReconstructorPara::skelDisplaySize ) ;

		}
	}


	return ;

	//following not done
	// draw branch in joint 
	if( ! appstate.displayJoint )
		return ;

	if( this->brchNodeInJoint.size() == 0)
		return ;


	std::vector<std::vector<int>> nodes (smoothedBrchPts.size());
	for( int i=0; i<brchNodeInJoint.size(); ++i )
		nodes[brchNodeInJoint[i].x].push_back(brchNodeInJoint[i].y) ;

	for( int i=0; i<nodes.size(); ++i )
		std::sort( nodes[i].begin(), nodes[i].end() ) ;

	std::vector<std::pair<int,int>> brchNodeInJoint ;
	for( int i=0; i<nodes.size(); ++i )
		for( int j=0; j<nodes[i].size(); ++j )
			brchNodeInJoint.push_back(  std::pair<int,int>(i, nodes[i][j] ) );



	std::vector<Point3f> brchToDraw(1) ;
	brchToDraw[0] = smoothedBrchPts[ brchNodeInJoint[0].first ][ brchNodeInJoint[0].second ] ;
	for( int i=1; i<brchNodeInJoint.size(); ++i ){
		if( brchNodeInJoint[i].first == brchNodeInJoint[i-1].first && brchNodeInJoint[i].second == brchNodeInJoint[i-1].second+1 )
			brchToDraw.push_back( smoothedBrchPts[ brchNodeInJoint[i].first ][ brchNodeInJoint[i].second ] ) ;
		else{
			CurveArray3D cvs (1) ;
			cvs[0] = brchToDraw ;
			GlobalFun::draw3dCurves( cvs, GLColor(0.8,0.8,0.8) , false, false, 4.0 ) ;
			brchToDraw.clear() ;
			brchToDraw.push_back( smoothedBrchPts[ brchNodeInJoint[i].first ][ brchNodeInJoint[i].second ] ) ;
		}
		CurveArray3D cvs (1) ;
		cvs[0] = brchToDraw ;
		GlobalFun::draw3dCurves( cvs, GLColor(0.8,0.8,0.8) , false, false, 4.0 ) ;
		brchToDraw.clear() ;
		brchToDraw.push_back( smoothedBrchPts[ brchNodeInJoint[i].first ][ brchNodeInJoint[i].second ] ) ;

	}


}


void skelpath::recenter(){

	if( !profiles3d.size() )
		return ;

	for( int i=0; i<profiles3d.size(); ++i)
		if( !profiles3d[i].size() ) return ;


	brchPts.clear(); brchPts.resize( profiles3d.size()) ;
	for( int i=0; i<profiles3d.size(); ++i ){

		std::vector<Point3f> newBranchPts ;
		//newBranchPts.push_back( smoothedBrchPts[i][0] ) ;

		for( int j=0; j<profiles3d[i].size(); ++j ){

			if( profile_conf[i][j] < 0.1 )
				continue ;

			if( profiles3d[i][j].size() <10 )
				continue ;

			Point3f vector1, vector2, center  ;
			GlobalFun::computePCAplane( profiles3d[i][j],  vector1, vector2, center ) ;

			//Point3f  center(0,0,0) ;
			//for( int id=0; id<profiles3d[i][j].size(); ++id )
			//	center += profiles3d[i][j][id] ;
			//center/=profiles3d[i][j].size() ;

			newBranchPts.push_back( center) ;

		}

		//newBranchPts.push_back( smoothedBrchPts[i].back() ) ;

		//brchPts[i] = newBranchPts ;

		ReconstructorUtility::getBezier( newBranchPts, smoothedBrchPts[i], 100 ) ;
		//ReconstructorUtility::getUniformCubicBSpline( newBranchPts, smoothedBrchPts[i], 100 ) ;
	}

	std::vector<Point3f> joints;
	std::vector<std::vector<int>> ajBrchId ;
	getJoints( joints, ajBrchId ) ;



	std::vector<Point3f> newJoints(joints.size(), Point3f(0,0,0) );
	for( int i=0; i<joints.size(); ++i ){
		for( int j=0; j<ajBrchId[i].size(); ++j ){
			int bid = ajBrchId[i][j] ;

			if( (smoothedBrchPts[bid][0]-joints[i]).Norm() <  (smoothedBrchPts[bid].back()-joints[i]).Norm() )
				newJoints[i] += smoothedBrchPts[bid][0] ;
			else
				newJoints[i] += smoothedBrchPts[bid].back();
		}
	}

	for( int i=0; i<joints.size(); ++i )
		joints[i] = newJoints[i]/ajBrchId[i].size() ;


	for( int i=0; i<joints.size(); ++i ){
		for( int j=0; j<ajBrchId[i].size(); ++j ){
			int bid = ajBrchId[i][j] ;

			if( (smoothedBrchPts[bid][0]-joints[i]).Norm() <  (smoothedBrchPts[bid].back()-joints[i]).Norm() )
				smoothedBrchPts[bid].insert( smoothedBrchPts[bid].begin(), joints[i]) ;
			else
				smoothedBrchPts[bid].push_back(joints[i]) ;
		}
	}

	brchPts = smoothedBrchPts ;
	smoothEachBranch();

}




void skelpath::computeBrchNodeInJoint() {


	brchNodeInJoint.clear() ;

	calculateProfiles() ;
	calculateProfileConfidence() ;

	// build a simple profile confidence array for branches

	std::vector<Point3f> joints;
	std::vector<std::vector<int>> ajBrchId ;

	getJoints( joints, ajBrchId ) ;

	for( int Jid = 0; Jid<ajBrchId.size(); ++Jid ){
		for( int jbid = 0; jbid<ajBrchId[Jid].size(); ++jbid ){

			int Bid = ajBrchId[Jid][jbid] ;
			int Nid = 0 ;
			bool start0 = true ;

			if( ( smoothedBrchPts[Bid].back()-joints[Jid]).Norm() <  ( smoothedBrchPts[Bid][0]-joints[Jid] ).Norm()  ){
				Nid = smoothedBrchPts[Bid].size() - 1 ;
				start0 = false ;
			}

			for( ; start0? Nid<smoothedBrchPts[Bid].size()/2:Nid>smoothedBrchPts[Bid].size()/2;   start0? ++Nid : --Nid ){

				std::vector<Point3f> slice = ReconstructorUtility::getSlice( fragments[Bid], cutplanes[Bid][Nid], smoothedBrchPts[Bid][Nid]) ;
				std::vector<Point3f> slice_big = ReconstructorUtility::getSlice( points, cutplanes[Bid][Nid], smoothedBrchPts[Bid][Nid]) ;					 


				if( slice.size() < 10 ){
					brchNodeInJoint.push_back( int2(Bid,Nid) ) ;
					continue;
				}


				double confi =  ReconstructorUtility::evaluateConfident(  ReconstructorUtility::convert3dProfileTo2d(slice, cutplanes[Bid][Nid], smoothedBrchPts[Bid][Nid]) ) ;
				if( confi!=0)
					confi = (confi - profconf_maxmin[Bid].y)/(profconf_maxmin[Bid].x-profconf_maxmin[Bid].y) ;
				if( confi < 0.1 ){
					brchNodeInJoint.push_back( int2(Bid,Nid) ) ;
					continue;
				}

				// put all points in slice ahead
				for( int i=0; i<slice_big.size(); ++i )
					for( int j=0; j<slice.size(); ++j )
						if( slice[j] == slice_big[i]){
							slice_big.erase(slice_big.begin()+i) ;
							--i ;
							break;
						}
						slice_big.insert(slice_big.begin(), slice.begin(),  slice.begin()+slice.size()-1 ) ;


						std::vector<std::vector<int>> nstId ;
						GlobalFun::computeKNN( slice_big, slice,nstId,5);

						// if number of nstId in other fragments are below a value, then its the boundary slices
						int count = 0;
						for( int i=0; i<nstId.size(); ++i )
							for( int j=0;j<nstId[i].size(); ++j )
								if( nstId[i][j] > slice.size() )
									count ++ ;

						if( count < nstId.size() * nstId[0].size() * 0.01 )
							break ;
						else
							brchNodeInJoint.push_back( int2(Bid,Nid) ) ;
			}


		} // for( int jbid = 0; jbid<ajBrchId.size(); ++jbid )


	}//for( int Jid = 0; Jid<ajBrchId.size(); ++Jid )


}


void skelpath::compute3dJoints( ) {

	// calculate the joints and save it in the member variable "joints3d"
	std::cout<<"skelpath::compute3dJoints"<<std::endl;

	if( fragments.size() == 0 )
		return ;

	// get Joints and their adjacent branches
	std::vector<Point3f> joints;
	std::vector<std::vector<int>> ajBrchId ;
	getJoints( joints, ajBrchId ) ;

	// build the relationship between joints and brachNode 
	brchEndForJoints.clear() ;
	for( int i=0; i<ajBrchId.size(); ++i ){
		brchEndForJoints.resize(brchEndForJoints.size()+1) ;
		for( int j=0; j<ajBrchId[i].size(); ++j  ){

			int bid = ajBrchId[i][j] ;
			if( ( smoothedBrchPts[bid][0] -  joints[i] ).Norm() < 0.05 )
				brchEndForJoints.back().push_back( int2(bid, 0 )) ;

			if( ( smoothedBrchPts[bid].back() -  joints[i] ).Norm() < 0.05 )
				brchEndForJoints.back().push_back( int2(bid, smoothedBrchPts[bid].size()-1 )) ;

		}
	}

	//  compute joints 
	calculateCutplanes() ;

	joints3d.clear() ;
	for( int Jid=0; Jid < brchEndForJoints.size(); ++Jid ){
		std::vector<Point3f> jointPoints ;

		for( int id = 0; id<brchEndForJoints[Jid].size(); ++id ){
			int Bid = brchEndForJoints[Jid][id].x ;

			int Nid = brchScopes[Bid].x;  
			Point3f cutplaneNormal = -cutplanes[Bid][Nid].planeNormal ;

			if( brchEndForJoints[Jid][id].y > 0 ){
				Nid = brchScopes[Bid].y;
				cutplaneNormal = cutplanes[Bid][Nid].planeNormal ;
			}


			std::vector<Point3f> points = ReconstructorUtility::cutPointsSet( fragments[Bid],smoothedBrchPts[Bid][Nid], cutplaneNormal  ) ;
			if( points.size()> 0)
				jointPoints.insert( jointPoints.begin(), points.begin(), points.begin() + points.size()-1 ) ;


		}

		joints3d.push_back( jointPoints ) ;
	}

}
