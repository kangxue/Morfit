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


#include "GlobalFunction.h"
#include "skelpath.h"
#include "bspline.h"


#include <fstream>
extern std::ofstream logout ;

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <algorithm>

#include <GL/GL.h>

#include <string> ;

#include "reconstructionUtility.h"
#include "reconstructorPara.h"

#include "declarations.h"


//added by huajie,2014/8/15
extern std::vector<Point3f> allSkeletonEnds;
extern std::vector<Point3f> BranchEnds ;

skelpath::skelpath( std::vector<std::vector<Point3f>> & brchPts, bool dsample = true) {
	iditifier = clock(); 

	// filter each branch
	for( int bid=0;bid<brchPts.size(); ++bid ){

		for( int i=0; i< brchPts[bid].size(); ++i )
			for( int j=i+1; j< brchPts[bid].size(); ++j )
				if( (brchPts[bid][i] - brchPts[bid][j]).Norm() < 1e-3 ){
					brchPts[bid].erase(  brchPts[bid].begin() + j ) ;
					j-- ;
				}

	}



	// convert vertices list to graph
	skeletonGraph = new lemon::ListGraph ;
	graphNodeP3f = new lemon::ListGraph::NodeMap<Point3f>(*skeletonGraph) ;

	this->brchPts = brchPts ;


	convertBranchPointsToVeticesList( 0.001 ) ;
	convertVerticesListToGraph() ;
	convertGraphToVerticeList() ;
	convertVeticesListToBranchPoints() ;

	smoothEachBranch() ;
	calculateCutplanes() ;

	branchIdofTemplatetoDraw = 0 ;

	jointSettingNum = 0 ;
	//points = pointset ;


	//	SelectedTemplate = int2(-1,-1) ;

	//ctrTmpValid = false ;

	activeCVID = -1;
	activeTmpId = -1;




	firstMustBeTip = false;
	firstMustBeNotTip  = false;
	lastMustBeTip  = false;
	lastMustBeNotTip  = false;


	brch0IsLoop = false ;

	cvSelectionType = 1 ;

}




skelpath::skelpath( completionSolver solver ) {

	iditifier = clock(); 

	brchPts = std::vector<Curve3D>(1, solver.centers ) ;

	// convert vertices list to graph
	skeletonGraph = new lemon::ListGraph ;
	graphNodeP3f = new lemon::ListGraph::NodeMap<Point3f>(*skeletonGraph) ;

	this->brchPts = brchPts ;


	convertBranchPointsToVeticesList( 0.001 ) ;
	convertVerticesListToGraph() ;
	convertGraphToVerticeList() ;
	convertVeticesListToBranchPoints() ;


	smoothedBrchPts = std::vector<Curve3D>(1, solver.centers ) ;
	cutplanes = std::vector<std::vector<cutPlane>>(1, solver.cutplanes) ;
	profiles2d = std::vector<std::vector<Curve2D>>(1, solver.profiles2d) ;

	tempalteIdss.resize(1) ; tempalteIdss[0] = solver.templatesIds;  
	Nurbsss.resize(1) ; 	Nurbsss[0] = solver.Tmps ;

	discretizeNurbs() ;


	branchIdofTemplatetoDraw = 0 ;
	jointSettingNum = 0 ;
	activeCVID = -1;
	activeTmpId = -1;

	firstMustBeTip = false;
	firstMustBeNotTip  = false;
	lastMustBeTip  = false;
	lastMustBeNotTip  = false;


	brch0IsLoop = solver.isLoop ;

	solvers = std::vector<completionSolver>(1,solver) ;



	std::vector<std::vector<Point3f>> profiles3d ;
	std::vector<std::vector<Point2f>> profiles2d ;
	for( int i=0; i<solver.finalProfiles.size(); ++i ){

		Profile3D p3d =  ReconstructorUtility::convert2dProfileTo3d( solver.finalProfiles[i], solver.cutplanes[i], solver.centers[i] ) ;
		profiles3d.push_back( p3d ) ;
		profiles2d.push_back( solver.finalProfiles[i] ) ;

	}
	resultProf3d.clear() ;
	resultProf2d.clear() ;
	resultProf3d.push_back( profiles3d ) ;
	resultProf2d.push_back( profiles2d ) ;

}

void skelpath::convertVerticesListToGraph(){


	skeletonGraph->clear() ;

	typedef lemon::ListGraph::Node Node ;
	typedef lemon::ListGraph::NodeIt NodeIt ;

	std::vector< Node > nodeList ;
	for( int i=0; i< nodes.size(); ++i ){
		nodeList.push_back( skeletonGraph->addNode() ) ;
		(*graphNodeP3f) [ nodeList[i] ] = nodes[i] ;
	}

	for( int i=0; i<branches.size(); ++i ){
		for( int j=0; j<branches[i].size()-1; ++j)
			skeletonGraph->addEdge(nodeList[ branches[i][j] ], nodeList[ branches[i][j+1]] ) ;
	}

}

void skelpath::convertGraphToVerticeList(){

	nodes.clear() ;
	branches.clear() ;
	brchNodeDegree.clear() ;

	typedef lemon::ListGraph::Node Node ;
	typedef lemon::ListGraph::NodeIt NodeIt ;
	typedef lemon::ListGraph::Edge Edge ;
	typedef lemon::ListGraph::EdgeIt EdgeIt ;
	typedef lemon::ListGraph::IncEdgeIt IncEdgeIt ;

	typedef lemon::ListGraph::NodeMap<bool> BoolNodeMap ;

	lemon::ListGraph & g = *skeletonGraph ;

	// get nodeList
	std::vector< Node > nodeList ;
	for( NodeIt n(g) ; n!=lemon::INVALID; n++  ){
		nodeList.push_back( n ) ;
		nodes.push_back( (*graphNodeP3f ) [n] );

	}

	// start from node with degree <= 2, search utill reach end or joint.  
	bool finshed = false ;
	//std::vector<bool> nodeHasBeenSelected(nodeList.size(), false) ;
	BoolNodeMap nodeHasBeenSelected( g ) ;
	for( NodeIt n(g); n!=lemon::INVALID; n++ )
		nodeHasBeenSelected[n] = false ; 

	while( !finshed ){

		// test if there are seedNodes rest
		Node seedNode ;
		for( int i=0; i<nodeList.size(); ++i ){
			Node node = nodeList[i] ;
			if( !nodeHasBeenSelected[node] && (lemon::countIncEdges(g, node) == 2  || lemon::countIncEdges(g, node) == 1 )){ seedNode = node ;	break; }

			if( i==nodeList.size()-1 )
				finshed = true ;
		}

		// if there is
		if( !finshed ){

			std::vector< Node > nodeSubList ;
			nodeSubList.push_back(seedNode) ;
			nodeHasBeenSelected[seedNode] = true ;

			bool first = true ;
			for( IncEdgeIt e(g,seedNode );  e!=lemon::INVALID; e++  ){

				if( first ){
					// search forward, until the degree is not 2
					Node fwdNode =  g.oppositeNode(seedNode, e ) ;
					Edge lastEdge = e ;
					while( 1){

						if( lemon::countIncEdges(g, fwdNode) == 1 ||  lemon::countIncEdges(g, fwdNode) > 2 ||nodeHasBeenSelected[fwdNode] == true ){
							nodeSubList.push_back( fwdNode );
							nodeHasBeenSelected[fwdNode] = true ;
							break ;
						}else{

							nodeSubList.push_back( fwdNode );
							nodeHasBeenSelected[fwdNode] = true ;
						}

						for( IncEdgeIt iei(g, fwdNode) ;  iei!=lemon::INVALID; iei++  ){
							if( Edge(iei)!=lastEdge ){
								fwdNode = g.oppositeNode(fwdNode, iei ) ;
								lastEdge = iei ;
								break ;
							}
						}



					}

					first = false ;

				}else{		

					Edge lastEdge = e ;
					// search backward, until the degree is not 2
					Node fwdNode =  g.oppositeNode(seedNode, e ) ;
					while(1){

						if( lemon::countIncEdges(g, fwdNode) == 1 ||  lemon::countIncEdges(g, fwdNode) > 2 ||  nodeHasBeenSelected[fwdNode] == true ){
							nodeSubList.insert(nodeSubList.begin(), fwdNode );
							nodeHasBeenSelected[fwdNode] = true ;
							break ;
						}else{

							nodeSubList.insert(nodeSubList.begin(), fwdNode );
							nodeHasBeenSelected[fwdNode] = true ;
						}

						for( IncEdgeIt iei(g, fwdNode) ;  iei!=lemon::INVALID; iei++  ){
							if( Edge(iei)!=lastEdge ){
								fwdNode = g.oppositeNode(fwdNode, iei ) ;
								lastEdge = iei ;
								break ;
							}
						}

					}

				}
			}


			// convert nodeSubList
			branches.resize( branches.size()+1) ;
			brchNodeDegree.resize( brchNodeDegree.size()+1 ) ;
			for( int i=0; i<nodeSubList.size(); ++i ){
				for( int j=0; j<nodeList.size(); ++j ){
					if( nodeSubList[i] == nodeList[j] ){
						branches.back().push_back( j ) ;
						brchNodeDegree.back().push_back( lemon::countIncEdges(g,nodeSubList[i]) ) ;
						break ;
					}

				}
			}
		}


	}

	// store the edge that connect 2 joints as a branch
	for( EdgeIt e = EdgeIt(g); e!= lemon::INVALID; e++ ){

		if( lemon::countIncEdges(g, g.u(e) ) > 2 &&  lemon::countIncEdges(g, g.v(e) ) > 2 ){

			std::vector< Node > nodeSubList ;
			nodeSubList.push_back(g.u(e)) ;
			nodeSubList.push_back(g.v(e)) ;

			branches.resize( branches.size()+1) ;
			brchNodeDegree.resize( brchNodeDegree.size()+1 ) ;
			for( int i=0; i<nodeSubList.size(); ++i ){
				for( int j=0; j<nodeList.size(); ++j ){
					if( nodeSubList[i] == nodeList[j] ){
						branches.back().push_back( j ) ;
						brchNodeDegree.back().push_back( lemon::countIncEdges(g,nodeSubList[i]) ) ;
						break ;
					}

				}
			}

		}

	}
}


void skelpath::convertVeticesListToBranchPoints() {

	CurveArray3D curves ;
	for(int i = 0; i < branches.size(); i++)			{
		std::vector<Point3f> curve ;
		for( int pbid=0; pbid<branches[i].size(); ++pbid  )
			curve.push_back( nodes[ branches[i][pbid] ] ) ;
		curves.push_back(curve) ;
	}

	brchPts = curves ;

}


void skelpath::convertBranchPointsToVeticesList( double disToMerge ) {
	// convert points on each branch to vertices lists

	nodes.clear() ;
	branches.resize(brchPts.size()) ;

	std::vector<int> nodes_mutiplicity;

	for( int i=0;i<branches.size(); ++i )
		branches[i].resize(brchPts[i].size() );

	for( int i=0; i<brchPts.size(); ++i  ){
		for( int j=0; j<brchPts[i].size(); ++j ){

			bool hasBeenStored = false ;
			int foundId = 0 ;
			for( int nid=0; nid<nodes.size(); ++nid )
				if( ( brchPts[i][j] - nodes[nid]).Norm() < disToMerge ){

					nodes[nid] = ( nodes[nid] * nodes_mutiplicity[nid] + brchPts[i][j] ) / ( nodes_mutiplicity[nid] + 1) ;
					nodes_mutiplicity[nid]++ ;

					hasBeenStored = true ;
					foundId = nid ;
				}

				if( hasBeenStored )
					branches[i][j] = foundId ;
				else{
					nodes.push_back( brchPts[i][j] ) ;
					nodes_mutiplicity.push_back(1) ;
					branches[i][j] = nodes.size() - 1 ; 
				}

		}

	}
}



#include "cubicInterpolation.h"
void skelpath::smoothEachBranch() {

	smoothedBrchPts.clear() ;
	smoothedBrchPts.resize(brchPts.size()) ;
	for( int i=0; i<brchPts.size(); ++i ){
		if( brchPts[i].size()< 4 ){
			while( brchPts[i].size() < 4 ){
				std::vector<Point3f> tmp ;
				for( int j=0; j< brchPts[i].size()-1; ++j ){
					tmp.push_back( brchPts[i][j] ) ;
					tmp.push_back( (brchPts[i][j]+brchPts[i][j+1])/2 ) ;

				}
				tmp.push_back( brchPts[i].back() ) ;

				brchPts[i] = tmp ;
			}
		}

		std::vector<double> x, y, z, t ;
		for( int j=0; j<brchPts[i].size(); ++j){
			x.push_back(brchPts[i][j].X() );
			y.push_back(brchPts[i][j].Y() );
			z.push_back(brchPts[i][j].Z() );
			t.push_back( (double)j/(brchPts[i].size()-1) ) ;
		}

		Spline3Interp poly1( t, x);
		poly1.calcCoefs();
		std::vector<double> px,py, pz ;
		for( double t=0.0; t<1.0; t+=0.0005 )
			px.push_back(poly1.evaluate(t)) ;
		px.push_back(poly1.evaluate(1.0)) ;


		Spline3Interp poly2( t, y);
		poly2.calcCoefs();
		for( double t=0.0; t<1.0; t+=0.0005 )
			py.push_back(poly2.evaluate(t)) ;
		py.push_back(poly2.evaluate(1.0)) ;

		Spline3Interp poly3( t, z);
		poly3.calcCoefs();
		for( double t=0.0; t<1.0; t+=0.0005 )
			pz.push_back(poly3.evaluate(t)) ;
		pz.push_back(poly3.evaluate(1.0)) ;

		for( int id=0; id<px.size(); ++id )
			smoothedBrchPts[i].push_back( Point3f(px[id], py[id], pz[id]) ) ;

	}

	for( int i=0; i<smoothedBrchPts.size() ; ++i ){
		for( int j=1; j<smoothedBrchPts[i].size()-1; ++j ){
			if( (smoothedBrchPts[i][j-1] - smoothedBrchPts[i][j]).Norm() < ReconstructorPara::branchSampleStep ){
				smoothedBrchPts[i].erase( smoothedBrchPts[i].begin()+j) ;
				j-- ;
			}

			if( (smoothedBrchPts[i][smoothedBrchPts[i].size()-1] - smoothedBrchPts[i][smoothedBrchPts[i].size()-2]).Norm() <   ReconstructorPara::branchSampleStep * 0.5)
				smoothedBrchPts[i].erase( smoothedBrchPts[i].begin()+smoothedBrchPts[i].size()-2) ;
		}
	}

}

void skelpath::getJoints( std::vector<Point3f> &joints, std::vector<std::vector<int>> &ajBrchId ){

	typedef lemon::ListGraph::Node Node ;
	typedef lemon::ListGraph::NodeIt NodeIt ;
	typedef lemon::ListGraph::Edge Edge ;
	typedef lemon::ListGraph::EdgeIt EdgeIt ;
	typedef lemon::ListGraph::IncEdgeIt IncEdgeIt ;

	typedef lemon::ListGraph::NodeMap<bool> BoolNodeMap ;

	lemon::ListGraph & g = *skeletonGraph ;

	joints.clear();
	ajBrchId.clear();

	for( NodeIt n(g); n!=lemon::INVALID; ++n ){
		if( lemon::countIncEdges(g, n) > 2){
			joints.push_back( (*graphNodeP3f)[n] ) ;

			ajBrchId.resize( ajBrchId.size()+1 ) ;
			ajBrchId.back().resize( lemon::countIncEdges(g, n) ) ;
		}
	}

	//added by huajie,2014/8/15
	//////////
	allSkeletonEnds.clear();
	for( NodeIt n(g); n!=lemon::INVALID; ++n ){
		if( lemon::countIncEdges(g, n) == 1){
			allSkeletonEnds.push_back( (*graphNodeP3f)[n] ) ;
		}
	}
	//////////

	//std::cout<<"////joints.size() = "<<joints.size()<<endl;
	for( int jid=0; jid<joints.size(); ++jid ){


		std::vector<bool> hasBeenAdded( smoothedBrchPts.size(), false) ;

		for( int i=0; i<ajBrchId[jid].size(); ++i ){

			int minBid = 0 ;
			double mindis = 1e10 ;
			for( int j=0; j<smoothedBrchPts.size(); ++j ){
				if(  !hasBeenAdded[j] && (smoothedBrchPts[j][0]- joints[jid]).Norm() < mindis ){
					mindis = (smoothedBrchPts[j][0]- joints[jid]).Norm() ;
					minBid = j ;
				}
				if(  !hasBeenAdded[j] && (smoothedBrchPts[j].back()- joints[jid]).Norm() < mindis ){
					mindis = (smoothedBrchPts[j].back()- joints[jid]).Norm() ;
					minBid = j ;
				}
			}

			ajBrchId[jid][i] = minBid ;
			hasBeenAdded[minBid] = true ;

		}
	}


}


void skelpath::getTips( std::vector<int2> &endpoints) {

	typedef lemon::ListGraph::Node Node ;
	typedef lemon::ListGraph::NodeIt NodeIt ;
	typedef lemon::ListGraph::Edge Edge ;
	typedef lemon::ListGraph::EdgeIt EdgeIt ;
	typedef lemon::ListGraph::IncEdgeIt IncEdgeIt ;

	typedef lemon::ListGraph::NodeMap<bool> BoolNodeMap ;

	lemon::ListGraph & g = *skeletonGraph ;

	endpoints.clear() ;
	BranchEnds.clear();
	for( NodeIt n(g); n!=lemon::INVALID; ++n ){
		if( lemon::countIncEdges(g, n) == 1
			&& find(allSkeletonEnds.begin(), allSkeletonEnds.end(), (*graphNodeP3f)[n]) != allSkeletonEnds.end() ){ //modified by huajie,2014/8/15
				BranchEnds.push_back( (*graphNodeP3f)[n] ) ;
		}
	}
	////added by huajie 2014/8/14
	//for( int eid = 0; eid<allSkeletonEnds.size(); ++eid )
	//	std::cout<<"////~~allEnds["<<eid<<"] = "<<allSkeletonEnds[eid].X()<<","<<allSkeletonEnds[eid].Y()<<","<<allSkeletonEnds[eid].Z()<<endl;

	//std::cout<<"////BranchEnds.size() = "<<BranchEnds.size()<<endl;
	for( int eid = 0; eid<BranchEnds.size(); ++eid ){
		/*std::cout<<"////BranchEnds["<<eid<<"] = "<<BranchEnds[eid].X()<<","<<BranchEnds[eid].Y()<<","<<BranchEnds[eid].Z()<<endl;*/
		double minDis = 1e10;
		int2 minId(-1,-1) ;
		for( int bid = 0; bid<smoothedBrchPts.size(); ++bid ){

			if(  (smoothedBrchPts[bid][0]- BranchEnds[eid]).Norm() < minDis ){
				minDis = (smoothedBrchPts[bid][0]- BranchEnds[eid]).Norm() ;
				minId = int2(bid,0) ;
			}
			if(  (smoothedBrchPts[bid].back()- BranchEnds[eid]).Norm() < minDis ){
				minDis = (smoothedBrchPts[bid].back()- BranchEnds[eid]).Norm()  ;
				minId =  int2(bid,smoothedBrchPts[bid].size()-1) ;
			}

		}

		endpoints.push_back(minId) ;
	}
}