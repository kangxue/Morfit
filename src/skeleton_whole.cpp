/************************************************************************/
/*                       skel_datastructureUtility.cpp                  */
/************************************************************************/
#include "GlobalFunction.h"
#include "skeleton_whole.h"
#include "bspline.h"
#include <fstream>
extern std::ofstream logout ;
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <algorithm>    // std::sort
#include <GL/GL.h>
#include <string> ;
#include "reconstructionUtility.h"
#include "reconstructorPara.h"
#include "declarations.h"

#include "ellipse.h"
#include "appstate.h"
extern AppState appstate ;
#include <fstream>
extern std::ofstream logout ;

#include "GLArea.h"
extern GLArea *glarea_global_ptr ;
#include "completionSolver_v2.h"

//added by huajie,2014/8/15
std::vector<Point3f> allSkeletonEnds;
std::vector<Point3f> BranchEnds ;

skeleton_whole::skeleton_whole( std::vector<std::vector<Point3f>> & brchPts, bool dsample = true) {
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

	activeCVID = -1;
	activeTmpId = -1;




	firstMustBeTip = false;
	firstMustBeNotTip  = false;
	lastMustBeTip  = false;
	lastMustBeNotTip  = false;


	brch0IsLoop = false ;

	cvSelectionType = 1 ;

}



skeleton_whole::skeleton_whole( completionSolver solver ) {

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
	Nurbsss.resize(1) ; 	Nurbsss[0] = solver.Tmps ;  discretizeNurbs() ;


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


void skeleton_whole::convertVerticesListToGraph(){


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

void skeleton_whole::convertGraphToVerticeList(){

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



void skeleton_whole::convertVeticesListToBranchPoints() {

	CurveArray3D curves ;
	for(int i = 0; i < branches.size(); i++)			{
		std::vector<Point3f> curve ;
		for( int pbid=0; pbid<branches[i].size(); ++pbid  )
			curve.push_back( nodes[ branches[i][pbid] ] ) ;
		curves.push_back(curve) ;
	}

	brchPts = curves ;

}


void skeleton_whole::convertBranchPointsToVeticesList( double disToMerge ) {
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
void skeleton_whole::smoothEachBranch() {


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



void skeleton_whole::getJoints( std::vector<Point3f> &joints, std::vector<std::vector<int>> &ajBrchId ){

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


void skeleton_whole::getTips( std::vector<int2> &endpoints) {

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





/************************************************************************/
/*                             skel_blend.cpp                           */
/************************************************************************/


void skeleton_whole::removeInnerPoints(){

	std::cout<<"skeleton_whole::removeInnerPoints"<<std::endl;


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

std::vector<Point3f> skeleton_whole::blendBranchBetweenJoints(){
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



	std::vector<Point3f> points;
	for( int i=0; i<subpoints.size(); ++i )
		points.insert(points.begin(), subpoints[i].begin(), subpoints[i].end() ) ;


	for( int i=0; i<PointSet1.size(); ++i )
		points.insert(points.begin(), PointSet1[i].begin(), PointSet1[i].end() ) ;
	for( int i=0; i<PointSet2.size(); ++i )
		points.insert(points.begin(), PointSet2[i].begin(), PointSet2[i].end() ) ;

	return points ;
}



/************************************************************************/
/*        skel_getprofiles.cpp                                          */
/************************************************************************/

void skeleton_whole::calculateCutplanes(){


	cutplanes.clear() ;
	cutplanes.resize(smoothedBrchPts.size()) ;
	for( int i=0; i<smoothedBrchPts.size(); ++i ){

		cutplanes[i].resize(smoothedBrchPts[i].size() ) ;

		for( int j=0; j<smoothedBrchPts[i].size(); ++j ){

			if( smoothedBrchPts[i].size() < 3 )
				continue ;

			Point3f p = smoothedBrchPts[i][j] ;

			int id0 = j-1 ;
			int id1 = j ;
			int id2 = j+1 ;

			if( j== smoothedBrchPts[i].size() - 1 ) { id0--; id1-- ; id2-- ;}
			if( j== 0 ) { id0++; id1++ ; id2++;}

			Point3f v1 = smoothedBrchPts[i][id2] - smoothedBrchPts[i][id1] ;
			Point3f v2 = smoothedBrchPts[i][id1] - smoothedBrchPts[i][id0] ;

			v1.Normalize() ;
			v2.Normalize() ;

			Point3f NormOfPlane = (v1+v2).Normalize() ;   // tangent of the curve
			if( j== smoothedBrchPts[i].size() - 1 ) 
				NormOfPlane = v1 ;
			if( j==0 ) 
				NormOfPlane = v2 ;

			Point3f radius ;


			if( j==0){
				// choose the  intercept with minimal absolute value to compute a vector that can represent this plane
				// -- compute intercept by point-normal equation, A(x-x0)+B(y-y0)+C(z-z0)=0
				double tx = ( NormOfPlane[1] *  p[1] +  NormOfPlane[2] *  p[2] ) / NormOfPlane[0] + p[0] ;
				double ty = ( NormOfPlane[0] *  p[0] +  NormOfPlane[2] *  p[2] ) / NormOfPlane[1] + p[1] ;
				double tz = ( NormOfPlane[0] *  p[0] +  NormOfPlane[1] *  p[1] ) / NormOfPlane[2] + p[2] ;

				if( tx <= ty && tx <= tz )
					radius = Point3f(tx,0,0) - p ;
				else if( ty <= tz )
					radius = Point3f(0,ty,0) - p ;
				else
					radius = Point3f(0,0,tz) - p ;

			}else{
				// compute the projection of last radius on this plane, by rotation

				Point3f axis = cutplanes[i][j-1].planeNormal^NormOfPlane ;
				double degree = ReconstructorUtility::vectorIntersectionAngle( cutplanes[i][j-1].planeNormal,  NormOfPlane ) ;
				GlobalFun::Rotate_Point3D(degree, axis, cutplanes[i][j-1].cutRadius, radius ) ;

			}


			cutplanes[i][j].cutRadius = radius  ;
			cutplanes[i][j].planeNormal = NormOfPlane  ;

		}

	}

	// mark,  smooth cutplane for loop
	/// 


	if( brch0IsLoop ){

		std::vector<cutPlane> &ctpl0s = cutplanes[0] ;

		Point3f rad0 = ctpl0s[0].cutRadius ;
		Point3f rade = ctpl0s.back().cutRadius ;
		Point3f norm0 = ctpl0s[0].planeNormal ;
		Point3f norme = ctpl0s.back().planeNormal ;


		Point3f rade_proj ;

		Point3f axis =norme^norm0 ;
		double degree = ReconstructorUtility::vectorIntersectionAngle( norme, norm0 ) ;
		GlobalFun::Rotate_Point3D(degree, axis,rade, rade_proj ) ;


		degree = ReconstructorUtility::vectorIntersectionAngle(  rad0, rade_proj ) ;

		if( (rad0^rade_proj) * norm0 < 0  ) degree = -degree ;


		double delta_degree = degree / (ctpl0s.size()-1) ;

		for( int i=0; i<ctpl0s.size(); ++i ){

			double angle = delta_degree * (ctpl0s.size()-1 - i) ;

			GlobalFun::Rotate_Point3D( angle, ctpl0s[i].planeNormal, ctpl0s[i].cutRadius, ctpl0s[i].cutRadius ) ;

		}

	}

	///


}



void skeleton_whole::calculateProfiles(  ) {


	// mark , need to be rewrite

	std::vector<Point3f>  &pts = points ;
	double radius = 1;

	profiles3d.clear() ;  
	profiles3d.resize(smoothedBrchPts.size()) ;
	for( int i=0; i<smoothedBrchPts.size(); ++i ){

		profiles3d[i].resize(smoothedBrchPts[i].size() ) ;

		// calculate profiles from the segment result
		std::vector<Point3f> frag;
		if( fragments.size() ){
			frag = fragments[i] ;
			double maxdis = 0;
			for( int pid=0; pid<points.size(); ++pid ){
				if( ptMap2BrchNode[pid].x == i ){
					double dis = (points[pid]-smoothedBrchPts[i][ptMap2BrchNode[pid].y]).Norm() ;
					if( dis > maxdis )
						maxdis = dis ;
				}
			}
			radius = maxdis ;
		}
		else{
			frag = pts ;

			std::cout<<"error: no fragments!"<<std::endl;
			system("pause") ;
		}


		for( int pid = 0; pid<frag.size(); ++pid ){

			int j= fragments_nid[i][pid] ;
			if( fabs( ( frag[pid] - smoothedBrchPts[i][j] )* cutplanes[i][j].planeNormal)  < ReconstructorPara::branchSampleStep/2 )
				profiles3d[i][j].push_back( frag[pid] ) ;
		}
	}


	// downsample
	profiles3d_full = profiles3d ;
	for( int i=0; i<profiles3d.size(); ++i )
		for( int j=0; j<profiles3d[i].size(); ++j ){
			if( profiles3d[i][j].size() > ReconstructorPara::profileDownSampleNum ){
				std::random_shuffle( profiles3d[i][j].begin(),  profiles3d[i][j].end() ) ;
				profiles3d[i][j].erase( profiles3d[i][j].begin() + ReconstructorPara::profileDownSampleNum, profiles3d[i][j].end() ) ;
			}
		}

	profiles2d.clear() ;
	profiles2d.resize(profiles3d.size()) ;
	for( int i=0; i<profiles3d.size(); ++i ){

		profiles2d[i].resize(profiles3d[i].size()) ;
		for( int j=0; j<profiles3d[i].size(); ++j )
			profiles2d[i][j]  = ReconstructorUtility::convert3dProfileTo2d( profiles3d[i][j], cutplanes[i][j], smoothedBrchPts[i][j] ) ;

	}


}

void skeleton_whole::initAfterClaculateProfiles(){
	

	calculateSharpness() ;

	calculateProfileConfidence() ;

	selectTemplates() ;

	convertTemplate2Pixles() ;

	fitInitialNurbs() ;
	discretizeNurbs() ;

	convertDisNurbs2Pixles() ;
	initSharpLable() ;
}

bool mycompfunc (std::pair<double,int> a,std::pair<double,int> b) {    return ( a.first < b.first ) ; }


void skeleton_whole::calculateProfileConfidence(){

	profile_conf.clear() ;
	profile_conf.resize( profiles2d.size()) ;
	for( int i=0; i<profiles2d.size(); ++i ){

		profile_conf[i].resize( profiles2d[i].size() ) ;

		for( int j=0; j<profiles2d[i].size(); ++j ){

			if( profiles2d[i][j].size() <= 10 ){
				profile_conf[i][j] = 0.0 ;
				continue;
			}

			profile_conf[i][j] =  ReconstructorUtility::evaluateConfident(profiles2d[i][j]);

		}

	}


	// smooth and normalize confidence
	for( int i=0; i<profile_conf.size(); ++i ){

		ReconstructorUtility::smoothVectorData( profile_conf[i],5) ;

		std::vector<std::pair<double,int>> confIdx ;
		for( int j=0; j<profile_conf[i].size(); ++j )
			confIdx.push_back( std::pair<double,int>( profile_conf[i][j], j) ) ;
		std::sort(confIdx.begin(), confIdx.end() ) ;

		for( int j=0; j<confIdx.size(); ++j )
			profile_conf[i][confIdx[j].second] =  j/(double)(confIdx.size()-1) ;
	}
}


bool myintcompfunc (std::pair<int,int> a,std::pair<int,int> b) {    return ( a.first < b.first ) ; }
void skeleton_whole::selectTemplates() {

	int maxTnum = ReconstructorPara::maxTemplateNum ;
	tempalteIdss.clear() ;
	tempalteIdss.resize( profile_conf.size() ) ;
	for( int i=0; i<profile_conf.size(); ++i ){

		std::vector<double> profconf_i = profile_conf[i];

		int windsize = std::max( 10, (int)(15 * profconf_i.size()/100) );
		//int windsize = std::max( ReconstructorPara::profileNumEachBranch / 5, 7 ) ;
		double threshold = 0.3 ;

		for( int j=0; j<profconf_i.size(); ++j ){
			bool isTemplate = true ;
			if( profconf_i[j]<threshold ||profconf_i[j] == 0 )
				isTemplate = false ;
			for( int k = j-windsize/2; k<j+windsize/2; ++k  ){
				if( k< 0  || k>profconf_i.size()-1)
					continue ;
				if( profconf_i[k] > profconf_i[j] )
					isTemplate = false ;
			}

			if( isTemplate )
				tempalteIdss[i].push_back( j ) ;
		}


		// flip node in flipTemplatesId  
		for( int ftid=0; ftid<flipTemplatesId.size(); ++ftid )
			if( flipTemplatesId[ftid].x == i  ){
				int nid = flipTemplatesId[ftid].y ;
				int nstid = 0 ;
				int mindis = 1000000 ;
				for( int j=0; j<tempalteIdss[i].size(); ++j )
					if( (tempalteIdss[i][j] - nid) * (tempalteIdss[i][j] - nid) < mindis ){
						mindis = (tempalteIdss[i][j] - nid) * (tempalteIdss[i][j] - nid) ;
						nstid = j ;
					}


					if( abs( tempalteIdss[i][nstid]- nid ) < 2 )
						tempalteIdss[i].erase( tempalteIdss[i].begin() + nstid ) ;
					else
						tempalteIdss[i].push_back( nid ) ;

			}

		if(tempalteIdss[i].size() >  ReconstructorPara::maxTemplateNum && (appstate.plyfilename.find("interwinedPillar" ) == string::npos )){

			// preserve the max/min/middle 3 templates

			int maxTnum = ReconstructorPara::maxTemplateNum  ;

			int tn = tempalteIdss[i].size() ;
			std::vector< std::pair<int,int> > tids(tn) ;
			for( int id =0; id<tn; ++id ){
				tids[id].first =  tempalteIdss[i][id] ;
				tids[id].second = id ;
			}
			std::sort (tids.begin(), tids.end(), myintcompfunc );
			tempalteIdss[i].clear() ;
		

			for( int id=0; id<maxTnum; ++id )
				tempalteIdss[i].push_back( tids[tids.size()*id/maxTnum ].first ) ;

			// remove repeated templates
			for( int id=0; id<tempalteIdss[i].size(); ++id )
				for( int id1=id+1; id1<tempalteIdss[i].size(); ++id1 )
					if( tempalteIdss[i][id] == tempalteIdss[i][id1] ){
						tempalteIdss[i].erase( tempalteIdss[i].begin() + id1 ) ;
						id1-- ;
					}



		}
		if( tempalteIdss[i].size() == 0 ){
			int maxId =0;
			double maxConf = -1e10 ;
			for ( int j=profconf_i.size()/5; j<profconf_i.size()*4/5; ++j ){




				if( profconf_i[j] > maxConf  ){
					maxConf =  profconf_i[j] ;
					maxId = j ;
				}
			}

			tempalteIdss[i].push_back( maxId ) ;
		}

	}

	if( appstate.plyfilename.find("relief" ) != string::npos  ){
		tempalteIdss[0].clear() ;
		tempalteIdss[0].push_back( profile_conf[0].size()/2 ) ;

	}

}

void skeleton_whole::drawCutRadii() {

	glBegin( GL_LINES ) ;

	glColor3f(0.0,1.0,0 ) ;
	for( int i=0; i<smoothedBrchPts.size(); ++i ){
		for( int j=0; j<smoothedBrchPts[i].size(); ++j ){

			Point3f p1 = smoothedBrchPts[i][j] ;
			Point3f p2 = p1 + cutplanes[i][j].cutRadius * 0.1 ;
			Point3f p3 = p1 + cutplanes[i][j].planeNormal * 0.01 ;

			glVertex3f( p1[0], p1[1], p1[2] ) ;
			glVertex3f( p2[0], p2[1], p2[2] ) ;
			glVertex3f( p1[0], p1[1], p1[2] ) ;
			glVertex3f( p3[0], p3[1], p3[2] ) ;

		}
	}


	glEnd() ;


	glPointSize( 10.0 ) ;
	glBegin( GL_POINTS ) ;
	if( profiles3d.size() >=1  && profiles3d[0].size()>=1 )
		for( int i=0; i<profiles3d[0][0].size(); ++i  ){

			Point3f p1 = profiles3d[0][0][i] ;
			glVertex3f( p1[0], p1[1], p1[2] ) ;
		}


		glEnd() ;

}


void skeleton_whole::calculateSharpness(){

	std::cout<<"skeleton_whole::calculateSharpness() begin"<<std::endl;


	using namespace Eigen;

	sharpness.clear() ;

	sharpness.resize( profiles2d.size() ) ;
	for( int pro_i=0; pro_i<profiles2d.size(); ++pro_i ){


		sharpness[pro_i].resize( profiles2d[pro_i].size() ) ;

		for( int pro_j=0; pro_j<profiles2d[pro_i].size(); ++pro_j ){


			int pn = profiles2d[pro_i][pro_j].size() ;

			sharpness[pro_i][pro_j].resize( pn ) ;

			////-------------------- DO NOT COMPUTE SHARPNESS ----------------------

			//continue;

			for( int k=0; k<pn; ++k ){


				Point2f p =  profiles2d[pro_i][pro_j][k] ;
				std::vector<Point2f> neis ;

				// find 1/5 most nearest points to p
				int neiNum = pn / 5 ;
				std::vector< std::pair<double,int> > dis(pn) ;
				for( int id =0; id<pn; ++id ){
					dis[id].first =  ( profiles2d[pro_i][pro_j][id] - p ).Norm() ;
					dis[id].second = id ;
				}
				std::sort (dis.begin(), dis.end(), mycompfunc );   
				for( int id=0; id<neiNum; ++id)
					neis.push_back(profiles2d[pro_i][pro_j][dis[id].second ] ) ;


				//  --------------- begin PCA -------------



				unsigned int dimSpace = 2; // dimension space
				unsigned int m = 2;   // dimension of each point
				unsigned int n = neis.size();  // number of points 

				MatrixXd DataPoints(m, n);  // matrix (m x n)
				for( int id=0; id<n; ++id ){
					DataPoints(0,id) = neis[id][0] ;
					DataPoints(1,id) = neis[id][1] ;
				}

				if( k==0 ){
					for( int i=0; i<n;++i )
						logout <<"[ "<<DataPoints(0,i) <<" " << DataPoints(1,i) <<" ] " <<"," ;
					logout<<std::endl;


				}
				typedef std::pair<double, int> myPair;
				typedef std::vector<myPair> PermutationIndices;	

		
				double mean; VectorXd meanVector;
				for (int i = 0; i < DataPoints.rows(); i++){
					mean = (DataPoints.row(i).sum())/n;		 //compute mean
					meanVector  = VectorXd::Constant(n,mean); // create a vector with constant value = mean
					DataPoints.row(i) -= meanVector;
					// std::cout << meanVector.transpose() << "\n" << DataPoints.col(i).transpose() << "\n\n";
				}


				if( k==0 ){

					logout<<"after mean"<<std::endl;
					for( int i=0; i<n;++i )
						logout <<"[ "<<DataPoints(0,i) <<" " << DataPoints(1,i) <<" ] " <<"," ;
					logout<<std::endl;

				}

				// get the covariance matrix
				MatrixXd Covariance = MatrixXd::Zero(m, m);
				Covariance = (1 / (double) n) * DataPoints * DataPoints.transpose();
				//  std::cout << Covariance ;	

				if( k==0 ){

					logout<<"Covariance"<<std::endl;
					for( int i=0; i<2;++i )
						for( int j=0; j<2;++j )
							logout << Covariance(i,j) <<", ";
					logout<<std::endl;

				}


				// compute the eigenvalue on the Cov Matrix
				EigenSolver<MatrixXd> m_solve(Covariance);
				//std::cout << "Done";

				VectorXd eigenvalues = VectorXd::Zero(m);
				eigenvalues = m_solve.eigenvalues().real();

				MatrixXd eigenVectors = MatrixXd::Zero(n, m);  // matrix (n x m) (points, dims)
				eigenVectors = m_solve.eigenvectors().real();	


				if( k==0)
					for (unsigned int i = 0; i < m ; i++){
						logout <<  i << "-eigenvalue " << eigenvalues(i) << std::endl;

						logout << i << "-eigenvector " << eigenVectors.row(i) << std::endl;

						logout << std::endl ;
					}


					// save sharpness
					sharpness[pro_i][pro_j][k] = std::min( eigenvalues[0],eigenvalues[1] ) / ( eigenvalues[0] + eigenvalues[1] ) ;
			}
		}
	}


	// normalize sharpness
	double maxsharp = -100 ;
	double minsharp = 100 ;
	for( int i=0; i<sharpness.size(); ++i )
		for( int j=0; j<sharpness[i].size(); ++j )
			for( int k=0; k<sharpness[i][j].size(); ++k ){
				if( sharpness[i][j][k] > maxsharp ) maxsharp =  sharpness[i][j][k] ;
				if( sharpness[i][j][k] < minsharp ) minsharp =  sharpness[i][j][k] ;
			}

	double len = maxsharp - minsharp ;

	for( int i=0; i<sharpness.size(); ++i )
		for( int j=0; j<sharpness[i].size(); ++j )
			for( int k=0; k<sharpness[i][j].size(); ++k ){
				sharpness[i][j][k] = (sharpness[i][j][k] - minsharp) / len ;
			}

	std::cout<<"skeleton_whole::calculateSharpness() end"<<std::endl;
}

void skeleton_whole::switchTemplate( Point2f seed ){


	std::cout << "skel::switchTemplate"<<std::endl;
	std::cout << "seed = " << seed.X() <<" "<< seed.Y() << std::endl;

	// convert seed to world coordinate
	std::vector<Point2f> seedarray(1, seed) ;
	seedarray = ReconstructorUtility::convertScreen2dToWorld2d( seedarray) ;
	seed = seedarray[0] ;

	CurveArray2D brchPts2d( smoothedBrchPts.size() ) ;
	for( int i=0; i<smoothedBrchPts.size(); ++i )
		brchPts2d[i] =  ReconstructorUtility::project3dPointsOntoScreen(smoothedBrchPts[i]);

	// traverse all nodes
	int2 nstId(0,0) ;
	double minDis = 1e10 ;

	for( int i=0; i<brchPts2d.size(); ++i ){
		for( int j=0; j<brchPts2d[i].size(); ++j )
			if( (brchPts2d[i][j]-seed).Norm() < minDis  ){
				minDis = (brchPts2d[i][j]-seed).Norm() ;
				nstId = int2(i,j) ;
			}
	}

	if( minDis < 0.1 ){
		for( int i=0; i<flipTemplatesId.size(); ++i)
			if( flipTemplatesId[i] == nstId ){
				flipTemplatesId.erase( flipTemplatesId.begin() + i ) ;
				return ;
			}

		flipTemplatesId.push_back( nstId ) ;
		std::cout<<"flipTemplatesId addd <"<< nstId.X() <<","<<nstId.Y()<<">"<<std::endl;
	}

}


void skeleton_whole::switchTemplate_online(Point2f p) {

	switchTemplate(p) ;
	initAfterClaculateProfiles() ;
}



/************************************************************************/
/*                      skel_mesh.cpp                                   */
/************************************************************************/

void skeleton_whole::generateMesh() {

	resultMesh.clear() ;

	int psize = ReconstructorPara::finalProfilesSampleNum ;
	for( int pwbid = 0; pwbid<resultProf3d.size(); ++pwbid ) {

		//if( pwbid!=0 )
		//	continue ;

		std::vector<Point3f> vertices ;
		std::vector<Point3f> normals ;
		std::vector<int3> faces ;

		std::vector<Profile3D> &profs = resultProf3d[pwbid] ;

		for( int i=0; i<psize; ++i )
			vertices.push_back( profs[0][i] ) ;


		bool firstIsTip = solvers[pwbid].firstIsTip ;
		bool lastIsTip = solvers[pwbid].lastIsTip ;

		for( int i=0; i<profs.size()-1; ++i  ){

			// estimate offset
			std::vector<int> offsetCand ;
			for( int id=0; id<psize; id+=psize/10 ){
				int nstId = ReconstructorUtility::NearstPoint(profs[i+1], profs[i][id] ) ;
				offsetCand.push_back( (nstId - id + psize ) % psize ) ;
			}
			std::vector<int> diffSums ;
			for( int i=0; i<psize; ++i){
				double sum = 0;
				for( int i=0; i<offsetCand.size(); ++i )
					sum += (std::min)( (offsetCand[i] - i + psize) % psize , psize - (offsetCand[i] - i + psize) % psize ) ;
				diffSums.push_back(sum) ;
			}
			int maxSum = 0 ;  int offset=0;
			for( int i=0; i<diffSums.size(); ++i ){
				if( diffSums[i]>maxSum){ maxSum= diffSums[i]; offset = i;}
			}


			offset = 0 ;

			// generates triangles
			//  ------------ add vertices
			for( int id=0; id<psize; ++id )
				vertices.push_back( profs[i+1][id] ) ;

			while( normals.size()<vertices.size() )
				normals.push_back(Point3f(0,0,0)) ;



			for( int id=0; id<psize; ++id ){
				int id1 = (id+1)%psize ;

				int startId = i * psize ; 
				int startId1 = i * psize + psize ; 

				int3 t1 = int3(startId+id, startId + id1,                startId1 + (id1+offset)%psize ) ;
				int3 t2 = int3(startId+id, startId1 + (id1+offset)%psize, startId1 + (id+offset)%psize ) ;
				faces.push_back( t1 ) ;
				faces.push_back( t2 ) ;

				Point3f n1 = -(vertices[t1[1]] - vertices[t1[0]])^(vertices[t1[2]] - vertices[t1[1]]) ;
				Point3f n2 = -(vertices[t2[1]] - vertices[t2[0]])^(vertices[t2[2]] - vertices[t2[1]]) ;
				n1.Normalize() ;
				n2.Normalize() ;

				if( !(i==0&&firstIsTip) ){
					normals[t1[0]] += n1 ;     
					normals[t1[1]] += n1 ;
					normals[t1[2]] += n1 ;
				}

				if( !(i==profs.size()-1&&lastIsTip) ){
					normals[t2[0]] += n2 ;
					normals[t2[1]] += n2 ;
					normals[t2[2]] += n2 ;
				}


			}


		}

		if( brch0IsLoop ){
			// connect the first and last profile

			// estimate offset
			std::vector<int> offsetCand ;
			for( int id=0; id<psize; id+=psize/10 ){
				int nstId = ReconstructorUtility::NearstPoint(profs[0], profs.back()[id] ) ;
				offsetCand.push_back( (nstId - id + psize ) % psize ) ;
			}
			std::vector<int> diffSums ;
			for( int i=0; i<psize; ++i){
				double sum = 0;
				for( int i=0; i<offsetCand.size(); ++i )
					sum += (std::min)( (offsetCand[i] - i + psize) % psize , psize - (offsetCand[i] - i + psize) % psize ) ;
				diffSums.push_back(sum) ;
			}
			int maxSum = 0 ;  int offset=0;
			for( int i=0; i<diffSums.size(); ++i ){
				if( diffSums[i]>maxSum){ maxSum= diffSums[i]; offset = i;}
			}


			// generates triangles
			for( int id=0; id<psize; ++id ){
				int id1 = (id+1)%psize ;

				int startId  = (profs.size()-1) * psize ; 
				int startId1 = 0; 

				int3 t1 = int3(startId+id, startId + id1,                startId1 + (id1+offset)%psize ) ;
				int3 t2 = int3(startId+id, startId1 + (id1+offset)%psize, startId1 + (id+offset)%psize ) ;
				faces.push_back( t1 ) ;
				faces.push_back( t2 ) ;

			}


		}

		// calculate normal
		for( int i=0; i<profs.size(); ++i  ){
			for( int j=0; j<profs[i].size(); ++j ){

				if( i>0 && i<profs.size()-1 ){
					Point3f v1 = ( profs[i][(j+1)%psize] - profs[i][(j-1+psize)%psize] ).Normalize();
					Point3f v2 = ( profs[i+1][j] - profs[i-1][j] ).Normalize();

					normals[i * psize + j] = -v1^v2 ;
				}

				if( i==0&&firstIsTip ){
					Point3f n = solvers[pwbid].centers[0] -  solvers[pwbid].centers[1] ;
					normals[i * psize + j] = n ;
				}

				if( i==profs.size()-1&&lastIsTip ){
					Point3f n = solvers[pwbid].centers.back() -  solvers[pwbid].centers[ solvers[pwbid].centers.size()-2 ] ;
					normals[i * psize + j] = n ;
				}

			}

		}


		// update tip normal
		if( firstIsTip ){
			Point3f n = solvers[pwbid].centers[0] -  solvers[pwbid].centers[1] ;
			for( int i=0; i<psize; ++i)
				normals[i] = n ;
		}

		if( lastIsTip ){
			Point3f n = solvers[pwbid].centers.back() -  solvers[pwbid].centers[ solvers[pwbid].centers.size()-2 ] ;
			for( int i=0; i<psize; ++i)
				normals[ normals.size()-i-1 ] = n ;
		}


		// delete redundent tip vertices
		if( firstIsTip ){
			for( int i=0; i<psize-1; ++i) { normals.erase(normals.begin()); vertices.erase( vertices.begin()); }

			for( int i=0; i< faces.size(); ++i){
				if( faces[i].X() >= psize )  faces[i].X() -= psize-1 ;
				else  faces[i].X() = 0 ;

				if( faces[i].Y() >= psize )  faces[i].Y() -= psize-1 ;
				else  faces[i].Y() = 0 ;

				if( faces[i].Z() >= psize )  faces[i].Z() -= psize-1 ;
				else  faces[i].Z() = 0 ;
			}
		}
		if( lastIsTip ){
			for( int i=0; i<psize-1; ++i) { normals.erase(normals.begin() + normals.size()-1); vertices.erase(vertices.begin() + vertices.size()-1); }
			for( int i=0; i< faces.size(); ++i){
				if( faces[i].X() >= normals.size()-1 )  faces[i].X() =  normals.size()-1 ;
				if( faces[i].Y() >= normals.size()-1 )  faces[i].Y() =  normals.size()-1 ;
				if( faces[i].Z() >= normals.size()-1 )  faces[i].Z() =  normals.size()-1 ;
			}
		}


		// add endface ,,  do not share vertices with the extisting surface 
		if( !firstIsTip && !brch0IsLoop){
			Point3f cter = GlobalFun::centerOfPoints( profs[0] ) ;
			Point3f n = solvers[pwbid].centers[0] -  solvers[pwbid].centers[1] ;
			vertices.push_back( cter ) ;  int posmark = vertices.size();
			normals.push_back( n ) ;
			for( int i=0; i<psize; ++i ){
				vertices.push_back( vertices[i] ) ;
				normals.push_back( n ) ;
			}
			for( int i=0; i<psize; ++i )
				faces.push_back( int3(posmark-1, posmark+(i+1)%psize , posmark+i ) ) ;

		}
		if( !lastIsTip &&!brch0IsLoop ){
			Point3f cter = GlobalFun::centerOfPoints( profs.back() ) ;
			Point3f n = solvers[pwbid].centers.back() -  solvers[pwbid].centers[ solvers[pwbid].centers.size()-2 ] ;
			vertices.push_back( cter ) ;  int posmark = vertices.size();
			normals.push_back( n ) ;
			int pid = profs.size()-1 ;
			for( int i=0; i<psize; ++i ){
				vertices.push_back( vertices[pid*psize + i] ) ;
				normals.push_back( n ) ;
			}
			for( int i=0; i<psize; ++i )
				faces.push_back( int3(posmark-1 , posmark+i , posmark+(i+1)%psize) ) ;

		}


		// normalize normals
		for( int i=0; i<normals.size(); ++i ){
			normals[i].Normalize() ;
			//normals[i] *= -1 ;
		}

		// convert the vertices and faces to  c-style array

		triangleList TL( vertices, normals,faces ) ;

		resultMesh.push_back(TL) ;
	}

}


/************************************************************************/
/*                  skel_nurbs.cpp                                      */
/************************************************************************/

void skeleton_whole::normalize( Box3f box ) {


	std::cout<<" skeleton_whole::normalize"<<std::endl;
	std::cout << "bbox = "<<box.min.X() << " " <<box.min.Y() << " "<<box.min.Z() << "\n"
		<<box.max.X() << " "<<box.max.Y() << " "<<box.max.Z() << std::endl;

	float max_x = abs((box.min - box.max).X());
	float max_y = abs((box.min - box.max).Y());
	float max_z = abs((box.min - box.max).Z());
	float max_length = max_x > max_y ? max_x : max_y;
	max_length = max_length > max_z ? max_length : max_z;

	for(int i = 0; i < nodes.size(); i++){

		Point3f& p = nodes[i];

		p -= box.min;
		p /= max_length;

		p = (p - Point3f( 0.5*max_x/max_length,  0.5*max_y/max_length,  0.5*max_z/max_length));
	}


	convertVerticesListToGraph() ;
	convertVeticesListToBranchPoints() ;
	smoothEachBranch() ;
	calculateCutplanes() ;

}




void skeleton_whole::fitInitialNurbs(){

	Nurbsss.clear() ;
	Nurbsss.resize( tempalteIdss.size() ) ;
	for( int i=0; i<tempalteIdss.size(); ++i ){
		Nurbsss[i].resize( tempalteIdss[i].size() ) ;

		for( int j=0; j<tempalteIdss[i].size(); ++j ){
			//modified by huajie,2014.8.10
			Nurbsss[i][j] = ReconstructorUtility::fitAndUniformNURBSforProfile2d( profiles2d[i][tempalteIdss[i][j]] ) ;
			//Nurbsss[i][j].MakePeriodicUniformKnotVector() ;
			//Nurbsss[i][j] = ReconstructorUtility::fitNURBSforProfile2d( profiles2d[i][tempalteIdss[i][j]] ) ;

			int knotNum = Nurbsss[i][j].KnotCount() ;

		}
	}
}

void skeleton_whole::initSharpLable(){
	NurbsCVPixelsSharpLable.clear();
	NurbsCVPixelsSharpLable.resize( tempalteIdss[branchIdofTemplatetoDraw].size() ) ;
	
	int cvnum = Nurbsss[0][0].CVCount() ;
	for( int i=0; i<NurbsCVPixelsSharpLable.size(); ++i ){

		NurbsCVPixelsSharpLable[i] = std::vector<bool>( cvnum, false) ;

	}
}

void skeleton_whole::discretizeNurbs() {

	disNurbsss.clear() ;
	disNurbsss.resize( Nurbsss.size() ) ;

	disNurbsss3d.clear() ;
	disNurbsss3d.resize( Nurbsss.size() ) ;


	for( int i=0; i<Nurbsss.size(); ++i ){

		disNurbsss[i].resize( Nurbsss[i].size() ) ;
		disNurbsss3d[i].resize( Nurbsss[i].size() ) ;


		for( int j=0; j<Nurbsss[i].size(); ++j ){

			// discretize nurbs
			disNurbsss[i][j] = ReconstructorUtility::discretizeNurbs( Nurbsss[i][j],200 ) ;

			// convert it to 3d curve
			disNurbsss3d[i][j] = ReconstructorUtility::convert2dProfileTo3d(disNurbsss[i][j],cutplanes[i][tempalteIdss[i][j]], smoothedBrchPts[i][tempalteIdss[i][j]] ) ;

		}


	}



}



void skeleton_whole::convertTemplate2Pixles(){

	int bid = branchIdofTemplatetoDraw ;
	int n=tempalteIdss[bid].size() ;

	std::vector<Profile2D> profiles ;
	for( int i=0; i<n; ++i ){
		Profile2D profile =  profiles2d[bid][tempalteIdss[bid][i] ] ;
		profiles.push_back(profile) ;
	}

	Box2f box = GlobalFun::computeBox( profiles ) ;

	temPixels = ReconstructorUtility::cvtProf2Pixels( profiles ,box) ;

}
void skeleton_whole::convertDisNurbs2Pixles(){

	int bid = branchIdofTemplatetoDraw ;
	int n=tempalteIdss[bid].size() ;

	std::vector<Profile2D> profiles ;
	for( int i=0; i<n; ++i ){
		Profile2D profile =  profiles2d[bid][tempalteIdss[bid][i] ] ;
		profiles.push_back(profile) ;
	}


	Box2f box = GlobalFun::computeBox( profiles ) ;


	// convert discretized nurbs to pixels
	std::vector<Profile2D> curves ;
	for( int i=0; i<disNurbsss[bid].size(); ++i ){
		Profile2D curve =  disNurbsss[bid][i] ;
		curves.push_back(curve) ;
	}
	disNurbsssPixels = ReconstructorUtility::cvtProf2Pixels(curves, box) ;


	double halfwdth = (std::max)(  (std::max)( fabs( box.max.X() ), fabs( box.min.X() ) ) , (std::max)( fabs( box.max.Y() ), fabs( box.min.Y() ) )    ) ;  
	Box2f box2  ;
	box2.min = Point2f(-halfwdth,-halfwdth) * 1.3;
	box2.max = Point2f(halfwdth,halfwdth) * 1.3;


	// convert control points to pixels
	std::vector<Profile2D> ctrlpoinss ;
	for( int i=0; i<Nurbsss[bid].size(); ++i ){

		Profile2D  ctrlpoints ; 
		int cvnum = Nurbsss[bid][i].CVCount() - Nurbsss[bid][i].Degree() ;
		for( int cpid = 0; cpid<cvnum + Nurbsss[bid][i].Degree(); ++cpid ){
			ON_4dPoint p ;
			Nurbsss[bid][i].GetCV(cpid, p) ;
			Point2f cvp = Point2f(p[0]/p[3], p[1]/p[3]) ;
			
			// correct out out boundary control points
			if( cvp.X() < box2.min.X() ){ cvp.X() = box2.min.X() ; p[0] =  box2.min.X()* p[3] ; }
			if( cvp.X() > box2.max.X() ){ cvp.X() = box2.max.X() ; p[0] =  box2.max.X()* p[3] ; }

			if( cvp.Y() < box2.min.Y() ){ cvp.Y() = box2.min.Y() ; p[1] =  box2.min.Y() * p[3] ; }
			if( cvp.Y() > box2.max.Y() ){ cvp.Y() = box2.max.Y() ; p[1] =  box2.max.Y() * p[3] ; }

			Nurbsss[bid][i].SetCV(cpid, p) ;
			if( cpid < Nurbsss[bid][i].Degree() )
				Nurbsss[bid][i].GetCV(cpid+cvnum, p) ;

			ctrlpoints.push_back( cvp  ) ;
		}

		ctrlpoinss.push_back(ctrlpoints) ;
	}
	NurbsCVPixels = ReconstructorUtility::cvtProf2Pixels(ctrlpoinss, box) ;


}

void skeleton_whole::setSharpNurbsCV( int NurbsId, int CVId) {

	std::cout<<"void skeleton_whole::setSharpNurbsCV( int NurbsId, int CVId)"<<std::endl;
	std::cout<<"NurbsId="<<NurbsId<<std::endl;
	std::cout<<"CVId="<<CVId<<std::endl;
	
	// set sharp feature by weighting
	NurbsCVPixelsSharpLable[NurbsId][CVId] = true ;

	double w = ReconstructorPara::weight_sharpCV ;
	ON_4dPoint p ;
	Nurbsss[branchIdofTemplatetoDraw][NurbsId].GetCV(CVId, p);  p[0] *= w ; p[1] *= w ; p[2] *= w ; p[3] *= w ;
	Nurbsss[branchIdofTemplatetoDraw][NurbsId].SetCV( CVId, p) ;
	if( CVId < ReconstructorPara::NurbsDegree )
		Nurbsss[branchIdofTemplatetoDraw][NurbsId].SetCV( CVId + ReconstructorPara::cvNum, p) ;

	std::cout << "weight: " ;
	for( int i=0; i<Nurbsss[branchIdofTemplatetoDraw].size(); ++i )
		for( int j=0; j<Nurbsss[branchIdofTemplatetoDraw][i].CVCount(); ++j ){
			ON_4dPoint p ;
			Nurbsss[branchIdofTemplatetoDraw][i].GetCV(j, p); 
			std::cout << p.w <<" " ;

		}
		std::cout<<std::endl;

	discretizeNurbs() ;
	convertDisNurbs2Pixles();

	return ;



	// set sharp feature by knot multipilcity

	// set lable 
	NurbsCVPixelsSharpLable[NurbsId][CVId] = true ;

	std::cout<<"void skeleton_whole::setSharpNurbsCV( int NurbsId, int CVId)"<<std::endl;
	std::cout<<"NurbsId="<<NurbsId<<std::endl;
	std::cout<<"CVId="<<CVId<<std::endl;


	ON_NurbsCurve & nurbs0 = Nurbsss[branchIdofTemplatetoDraw][NurbsId] ;

	std::cout<<"degree = "<<nurbs0.Degree() <<std::endl;

	std::cout<<"before Set sharp feature:"<<std::endl;
	for( int k = 0; k<nurbs0.KnotCount(); ++k )
		std::cout <<  nurbs0.Knot(k) <<std::endl ;

	
	int cvnum = nurbs0.CVCount() ;

	if( CVId >= 3 && CVId < cvnum - 3){
		double midKnotValue = nurbs0.Knot(CVId+1) ;
		nurbs0.SetKnot(CVId+0, midKnotValue ) ;
		nurbs0.SetKnot(CVId+1, midKnotValue ) ;
		nurbs0.SetKnot(CVId+2, midKnotValue ) ;
	}
	else
	{
		if( CVId == 0 || CVId== cvnum-3 ){
			nurbs0.SetKnot(0, 0 ) ;
			nurbs0.SetKnot(1, 0 ) ;
			nurbs0.SetKnot(2, 0 ) ;
			nurbs0.SetKnot(cvnum-3, 1  ) ;
			nurbs0.SetKnot(cvnum-2, 1  ) ;
			nurbs0.SetKnot(cvnum-1, 1  ) ;
		}

		if( CVId == 1 || CVId== cvnum-2 ){
			nurbs0.SetKnot(1, 0 ) ;
			nurbs0.SetKnot(2, 0 ) ;
			nurbs0.SetKnot(3, 0 ) ;
			nurbs0.SetKnot(cvnum-2, 1  ) ;
			nurbs0.SetKnot(cvnum-1, 1  ) ;
			nurbs0.SetKnot(cvnum+0, 1  ) ;
		}
		if( CVId == 2 || CVId== cvnum-1 ){
			nurbs0.SetKnot(2, 0 ) ;
			nurbs0.SetKnot(3, 0 ) ;
			nurbs0.SetKnot(4, 0 ) ;
			nurbs0.SetKnot(cvnum-1, 1  ) ;
			nurbs0.SetKnot(cvnum+0, 1  ) ;
			nurbs0.SetKnot(cvnum+1, 1  ) ;
		}

	}




	discretizeNurbs() ;
	convertDisNurbs2Pixles();

	std::cout<<"\nafter Set sharp feature:"<<std::endl;
	for( int k = 0; k<nurbs0.KnotCount(); ++k )
		std::cout <<  nurbs0.Knot(k) <<std::endl ;

}



void skeleton_whole::cvtPixels2CV(){
	

	int bid = branchIdofTemplatetoDraw ;
	int n=tempalteIdss[bid].size() ;

	std::vector<Profile2D> profiles ;
	for( int i=0; i<n; ++i ){
		Profile2D profile =  profiles2d[bid][tempalteIdss[bid][i] ] ;
		profiles.push_back(profile) ;
	}
	Box2f box = GlobalFun::computeBox( profiles ) ;

	
	// convert pixels to control points
	std::vector<Profile2D> ctrlpoinss= ReconstructorUtility::cvtPixels2Prof( NurbsCVPixels, box ) ;

	for( int i=0; i<Nurbsss[bid].size(); ++i ){

		ON::point_style cvtype = Nurbsss[bid][i].CVStyle() ;

		int cvnum = Nurbsss[bid][i].CVCount() ;
		for( int cpid = 0; cpid<cvnum; ++cpid ){

			ON_4dPoint p ;
			Nurbsss[bid][i].GetCV(cpid, p);

			p[0] = ctrlpoinss[i][cpid].X()*p[3] ;
			p[1] = ctrlpoinss[i][cpid].Y()*p[3] ;
			p[2] = 0 ;

			Nurbsss[bid][i].SetCV( cpid,cvtype, p ) ;
		}
	}

}

void skeleton_whole::drawTemplates() {

	for( int i=0; i<temPixels.size(); ++i)
	GlobalFun::draw2dPointsOnScreen( temPixels[i], double3(0.5, 0.5, 1.0), 4 );

}

void skeleton_whole::draw3dNurbs(){

	if( disNurbsss3d.size() > 0){

		for( int i=0; i<disNurbsss3d[0].size(); ++i )
			if( i != activeTmpId )
				GlobalFun::draw3dCurves(std::vector<Curve3D>(1,disNurbsss3d[0][i]) , GLColor( 0.8,0, 0, 1), true, false, ReconstructorPara::nurbsDisplaySize ) ;

	}
	
}

void skeleton_whole::drawNurbs() {

	int bid = branchIdofTemplatetoDraw ;


	// draw sharp control point
	std::vector<Point2f> sharpcvpoint ;
	for( int i=0; i<NurbsCVPixels.size(); ++i){
		for( int j=0; j<NurbsCVPixels[i].size(); ++j )
			if( NurbsCVPixelsSharpLable[i][j] )
				sharpcvpoint.push_back( NurbsCVPixels[i][j] ) ;
	}

	GlobalFun::draw2dPointsOnScreen( sharpcvpoint, double3(0.0, 0.0, 0.0 ), 10 );

	// draw control points
	for( int i=0; i<NurbsCVPixels.size(); ++i){
		GlobalFun::draw2dPointsOnScreen( NurbsCVPixels[i], double3(1.0, 0.0, 0.0 ), 10 );
		GlobalFun::draw2dCurveOnScreen( NurbsCVPixels[i], double3(0.0, 0.0, 0.0 ), 1, false );
	}

	// draw nurbs
	for( int i=0; i<disNurbsssPixels.size(); ++i)
		GlobalFun::draw2dCurveOnScreen( disNurbsssPixels[i], double3(0.0, 1.0, 1.0 ), 4 );

}


void skeleton_whole::setActivBranch(int bid ){ 	
	branchIdofTemplatetoDraw = bid  ;


	convertTemplate2Pixles() ;

	convertDisNurbs2Pixles() ;
	initSharpLable() ;

}

/************************************************************************/
/*        skel_profEdit.cpp                                             */
/************************************************************************/
double cornerWidth = 400 ;
double planeWidth = 0.6 ;
double edgeWidth = 100 ;
Point2f profCenter ;

void skeleton_whole::selectActiveTmp( Point2f px ) {

	std::vector<Point3f> tmpCters ;
	for( int i=0; i<tempalteIdss[0].size(); ++i )
		tmpCters.push_back( smoothedBrchPts[0][ tempalteIdss[0][i] ] ) ;

	std::vector<Point2f> tmpCters_s2d = ReconstructorUtility::convertWorld2dToScreen2d( ReconstructorUtility::project3dPointsOntoScreen( tmpCters )) ;

	int nstId = ReconstructorUtility::NearstPoint( tmpCters_s2d, px ) ;
	if( (tmpCters_s2d[nstId]-px).Norm() < 40 ){
		if( activeTmpId == nstId ){
			activeTmpId = -1;
			activeCVID = -1 ;
		}
		else {
			activeTmpId = nstId ;
			activeCVID = -1 ;
		}
	}



	// update planeWidth & profileCter
	if( activeTmpId >= 0 ){

		if( profiles2d[0][tempalteIdss[0][activeTmpId]].size() >  5){
			Box2f box = GlobalFun::computeBox( profiles2d[0][tempalteIdss[0][activeTmpId]] ) ;
			planeWidth = 2 * std::max(  std::max( fabs(box.min.X()),  fabs(box.min.Y() ) ), std::max( fabs(box.max.X()),  fabs(box.max.Y() ) ) );
			planeWidth *= 1.3 ;

			profCenter = ( box.max + box.min ) * 0.5 ;
		}
		else {
			planeWidth = 0.3 ;
			profCenter *= 0 ;
		}
	}



	std::cout << "skeleton_whole::selectActiveTmp:\n" << "activeTmpId = " <<activeTmpId << "\nactiveCVID = "<<activeCVID<<std::endl;

}

void skeleton_whole::selectActiveTmpCV( Point2f px ){

	if( activeTmpId < 0 )
		return ;

	cutPlane ctpl = cutplanes[0][tempalteIdss[0][activeTmpId]] ;
	Point3f cter = smoothedBrchPts[0][ tempalteIdss[0][activeTmpId] ] ;

	std::vector<Point2f> cvs_l2d = ReconstructorUtility::getCvPts( Nurbsss[0][activeTmpId] ) ;
	std::vector<Point3f> cvs = ReconstructorUtility::convert2dProfileTo3d( cvs_l2d, ctpl, cter );

	std::vector<Point2f> cvs_s2d = ReconstructorUtility::convertWorld2dToScreen2d( ReconstructorUtility::project3dPointsOntoScreen( cvs )) ;

	int nstId = ReconstructorUtility::NearstPoint( cvs_s2d, px ) ;
	if( (cvs_s2d[nstId]-px).Norm() < 40 ){
		activeCVID = nstId ;
		cvSelectionType = 1 ;

	}
	//else
	//	activeCVID = -1 ;


	if( activeCVID == -1 ){ // select from corner

		cvs_s2d.clear();
		for( int i=0; i<cvs_l2d.size(); ++i )
			cvs_s2d.push_back(  (cvs_l2d[i] - profCenter) * cornerWidth/planeWidth + Point2f( 1,1) * (cornerWidth/2 + edgeWidth)  ) ;

		int nstId = ReconstructorUtility::NearstPoint( cvs_s2d, px ) ;
		if( (cvs_s2d[nstId]-px).Norm() < 40 ){
			activeCVID = nstId ;
			cvSelectionType = 2 ;
		}
	}


	std::cout << "skeleton_whole::selectActiveTmpCV:\n" << "activeTmpId = " <<activeTmpId << "\nactiveCVID = "<<activeCVID<<std::endl;
}

void skeleton_whole::moveActiveTmpCV( Point2f pixelPos ) {

	std::cout << " skeleton_whole::moveActiveTmpCV( "<<pixelPos.X()<<", "<<pixelPos.Y()<<")"<<std::endl;

	if( activeCVID<0)
		return ;

	cutPlane ctpl = cutplanes[0][tempalteIdss[0][activeTmpId]] ;
	Point3f cter = smoothedBrchPts[0][ tempalteIdss[0][activeTmpId] ] ;

	ON_NurbsCurve &nbs = Nurbsss[0][activeTmpId]  ;


	if( cvSelectionType == 1 ){

		Point3f pixelPos_l3d = ReconstructorUtility::convertWorld3dToLocal3d( ReconstructorUtility::convertScreen2dToWorld3d( std::vector<Point2f>(1,pixelPos) ) ) [0];
		Point3f original_l3d = ReconstructorUtility::convertWorld3dToLocal3d( std::vector<Point3f>(1, Point3f(0,0,0)) )[0] ;


		Point3f intersectPoint ;
		if( ReconstructorUtility::CalPlaneLineIntersectPoint(intersectPoint, ctpl.planeNormal, cter, original_l3d-pixelPos_l3d, pixelPos_l3d) ){

			Point2f dst = ReconstructorUtility::convert3dProfileTo2d( std::vector<Point3f>(1,intersectPoint),ctpl, cter )[0] ;

			std::vector<Point2f> cvs = ReconstructorUtility::getCvPts( Nurbsss[0][activeTmpId] ) ;
			cvs[activeCVID] = dst ;

			ReconstructorUtility::setCvOfNurbs( nbs, cvs) ;
			discretizeNurbs() ;
		}
	}else if( cvSelectionType == 2 ){


		std::vector<Point2f> cvs = ReconstructorUtility::getCvPts( Nurbsss[0][activeTmpId] ) ;

		Point2f dst =( pixelPos - Point2f( 1,1) * (cornerWidth/2 + edgeWidth) )  * ( planeWidth/cornerWidth)  +  profCenter ;
		cvs[activeCVID] = dst ;

		ReconstructorUtility::setCvOfNurbs( nbs, cvs) ;
		discretizeNurbs() ;
	}


}

void skeleton_whole::deselectActiveTmpCV() {
	activeCVID = -1 ;
}
void skeleton_whole::deselectActiveTmp() {
	activeTmpId = -1 ;
}

void skeleton_whole::drawActiveTmp() {

	if( activeTmpId<0)
		return ;

	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);




	cutPlane ctpl = cutplanes[0][tempalteIdss[0][activeTmpId]] ;
	Point3f cter = smoothedBrchPts[0][ tempalteIdss[0][activeTmpId] ] ;

	std::vector<Point3f> cvs = ReconstructorUtility::convert2dProfileTo3d( ReconstructorUtility::getCvPts( Nurbsss[0][activeTmpId] ), ctpl, cter );

	// draw control points
	for( int i=0; i<cvs.size(); ++i){
		if( i==activeCVID ){
			glColor4f(1,0,0, 1) ;
			glPointSize(30) ;
		}
		else{
			glColor4f(0,0,1, 1) ;
			glPointSize(20) ;
		}
		glBegin( GL_POINTS) ;
		glVertex3f( cvs[i][0],  cvs[i][1],  cvs[i][2] ) ;
		glEnd() ;
	}

	// draw line segment between control points
	glColor3f(0,0,0) ;
	glLineWidth(5) ;
	glBegin(GL_LINE_LOOP) ;
	for( int i=0; i<cvs.size(); ++i)
		glVertex3f( cvs[i][0],  cvs[i][1],  cvs[i][2] ) ;
	glEnd() ;


	if( activeTmpId >= 0 ){

		glColor4f(0.8,0,0, 0.7) ;
		glLineWidth(10) ;
		glBegin(GL_LINE_LOOP) ;
		cvs = disNurbsss3d[0][activeTmpId] ;
		for( int i=0; i<cvs.size(); ++i)
			glVertex3f( cvs[i][0],  cvs[i][1],  cvs[i][2] ) ;

		glEnd() ;
	}

}

void skeleton_whole::drawActiveTmp_cutplane(){

	if( activeTmpId<0)
		return ;

	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);




	cutPlane ctpl = cutplanes[0][tempalteIdss[0][activeTmpId]] ;
	Point3f cter = smoothedBrchPts[0][ tempalteIdss[0][activeTmpId] ] ;


	// draw plane
	std::vector<Point2f> vertices;
	vertices.push_back( Point2f(1,1) ) ;
	vertices.push_back( Point2f(1,-1) ) ;
	vertices.push_back( Point2f(-1,-1) ) ;
	vertices.push_back( Point2f(-1,1) ) ;

	for( int i=0; i<vertices.size(); ++i )
		vertices[i] *= planeWidth/2 ;

	std::vector<Point3f> v3d = ReconstructorUtility::convert2dProfileTo3d(vertices, ctpl, cter );


	glDisable( GL_CULL_FACE ) ;
	glDisable(GL_LIGHTING) ;
	glColor4f(1,0,1,0.3) ;


	if( appstate.path.find("cup") != string::npos)
		glColor4f(0,1,1,0.3) ;


	glBegin( GL_POLYGON ) ;
	for( int i=0; i<v3d.size(); ++i)
		glVertex3f( v3d[i][0],  v3d[i][1],  v3d[i][2] ) ;
	glEnd() ;
}

void skeleton_whole::drawActiveTmpOnCorner(){



	if( activeTmpId<0)
		return ;


	glEnable(GL_POINT_SMOOTH);   
	glEnable(GL_LINE_SMOOTH);   
	glEnable(GL_POLYGON_SMOOTH);

	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST); // Make round points, not square points   
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);  // Antialias the lines
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);


	double scaleRatio = cornerWidth / planeWidth ;


	glMatrixMode(GL_PROJECTION);

	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, glarea_global_ptr->width() ,0,glarea_global_ptr->height(), -1,1);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	glTranslatef( cornerWidth/2+edgeWidth, cornerWidth/2+edgeWidth, 0) ;

	glScalef(scaleRatio,scaleRatio, 1 ) ;

	glPushMatrix() ;
	glTranslatef( -profCenter.X(), -profCenter.Y(), 0) ;

	
	std::vector<Point2f> cvs2d =  ReconstructorUtility::getCvPts( Nurbsss[0][activeTmpId] ) ;
	std::vector<Point3f> cvs3d ; 
	for( int i=0; i<cvs2d.size(); ++i) 
		cvs3d.push_back( Point3f( cvs2d[i].X(), cvs2d[i].Y(), 0) ) ;

	// draw data points
	glColor3f(0,0,0) ;
	glPointSize(7) ;
	glBegin(GL_POINTS) ;
	for( int i=0; i<profiles2d[0][tempalteIdss[0][ activeTmpId] ].size(); ++i )
		glVertex3f( profiles2d[0][tempalteIdss[0][ activeTmpId]][i].X(), profiles2d[0][tempalteIdss[0][ activeTmpId]][i].Y(), 0 ) ;
	glEnd() ;

	// draw control points
	glDisable(GL_LIGHTING) ;
	for( int i=0; i<cvs2d.size(); ++i ){
		if( i==activeCVID )	glColor4f(1.0,0,0, 1) ;  
		else glColor4f(0,0,1,1) ;

		glPointSize( 20*(1 + 0.5*(int)( i==activeCVID )) ) ;
		glBegin(GL_POINTS);
		glVertex3f(cvs2d[i] .X(), cvs2d[i].Y(), 0.02  ) ;
		glEnd() ;
	}

	glDisable(GL_LIGHTING) ;
	glColor4f(0,0,1, 0.7);  glLineWidth(3);  
	glBegin(GL_LINE_LOOP) ;
	for( int i=0; i<cvs2d.size(); ++i )  glVertex3f(cvs2d[i].X(),cvs2d[i].Y(), 0.02 ) ;
	glEnd() ;


	// draw curves
	Curve2D  curve2d = ReconstructorUtility::discretizeNurbs(Nurbsss[0][activeTmpId],100) ;
	std::vector<Point3f> curve3d ; 
	for( int i=0; i<curve2d.size(); ++i) 
		curve3d.push_back( Point3f( curve2d[i].X(), curve2d[i].Y(), 0) ) ;

	//GlobalFun::draw3dCurves( std::vector<Curve3D>(1, curve3d ) ,GLColor( 0.8,0, 0, 0.7), true, false, ReconstructorPara::nurbsDisplaySize ) ;
	glDisable(GL_LIGHTING) ;
	glColor4f(0.8,0, 0, 0.7);  glLineWidth(20);  
	glBegin(GL_LINE_LOOP) ;
	for( int i=0; i<curve2d.size(); ++i )  glVertex3f(curve2d[i].X(),curve2d[i].Y(),0.01 ) ;
	glEnd() ;




	glPopMatrix() ;

	// draw plane
	std::vector<Point2f> vertices;
	vertices.push_back( Point2f(1,1) ) ;
	vertices.push_back( Point2f(1,-1) ) ;
	vertices.push_back( Point2f(-1,-1) ) ;
	vertices.push_back( Point2f(-1,1) ) ;

	for( int i=0; i<vertices.size(); ++i )
		vertices[i] *= planeWidth/2;

	std::vector<Point3f> v3d ;
	for( int i=0; i<vertices.size(); ++i) 
		v3d.push_back( Point3f( vertices[i].X(), vertices[i].Y(), 0) ) ;


	glDisable( GL_CULL_FACE ) ;
	glDisable(GL_LIGHTING) ;
	glColor4f(0,1,1,0.2) ;
	glBegin( GL_POLYGON ) ;
	for( int i=0; i<v3d.size(); ++i)
		glVertex3f( v3d[i][0],  v3d[i][1],  v3d[i][2] ) ;
	glEnd() ;

	// Closing 2D
	glPopMatrix(); // restore modelview
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);

}



void skeleton_whole::drawSpecialView(){

}



/************************************************************************/
/*                   skel_reconstruct.cpp                               */
/************************************************************************/



std::vector<Point3f> skeleton_whole::reconstruct_old() {
	// return point cloud
	std::vector<Point3f> points ;
	return points ;
	

	return points ;
}


std::vector<Point3f> skeleton_whole::reconstruct() {
	// return point cloud
	std::vector<Point3f> points ;


	if( solvers.size() == 0 || appstate.plyfilename.find("plant")!=string::npos ){

		if( tempalteIdss.size() == 0)
			return points ;

		solvers = generateJointSetting() ;

	}

	jointSettingNum = solvers.size() ;



	cutplanes_bp.clear() ;
	cutplanes_bp.resize( solvers.size() ) ;
	branches_bp.clear() ;
	branches_bp.resize( solvers.size() ) ;


	subpoints.clear() ;
	resultProf3d.clear() ;
	resultProf2d.clear() ;
	for( int iter =0; iter<solvers.size(); ++iter ){

		std::cout<<"iter="<<iter<<std::endl;

		completionSolver &comSolver = solvers[iter];

		char tmp[10] ;
		std::string xfilename = std::string("x") + std::string( itoa(iter, tmp,10) )+".txt" ;

		std::ifstream xifs(xfilename) ;
		std::vector<double> inx,outx ;
		//double tt ;
		//while( xifs >> tt )
		//	inx.push_back(tt);

		outx = comSolver.solve( inx ) ;

		std::ofstream xofs(xfilename) ;
		for( int i=0;i<outx.size(); ++i )
			xofs << outx[i] << std::endl;
		xofs.close() ;
		xifs.close() ;
		std::cout<<"mark--------1" <<std::endl;
		std::vector<cutPlane> cutpls = comSolver.cutplanes ;
		std::vector<Point3f> brch = comSolver.centers ;



		std::vector<std::vector<Point3f>> profiles3d ;
		std::vector<std::vector<Point2f>> profiles2d ;
		for( int i=0; i<comSolver.finalProfiles.size(); ++i ){

			Profile3D p3d =  ReconstructorUtility::convert2dProfileTo3d( comSolver.finalProfiles[i], cutpls[i], brch[i] ) ;
			profiles3d.push_back( p3d ) ;
			profiles2d.push_back( comSolver.finalProfiles[i] ) ;

			// store cutplane and plnormal of branch pairs for debuging
			cutplanes_bp[iter].push_back(cutpls[i] );
			branches_bp[iter].push_back(brch[i] );

		}

		resultProf3d.push_back( profiles3d ) ;
		resultProf2d.push_back( profiles2d ) ;


		std::cout<<"mark--------2" <<std::endl;

		std::vector<Point3f> subpoint = ReconstructorUtility::upsampleResult( profiles3d, int2(0, brch.size()-1 ), comSolver.firstIsTip, comSolver.lastIsTip ) ;
		points.insert(points.begin(), subpoint.begin(), subpoint.end() ) ;

		subpoints.push_back( subpoint ) ; 	// for blend

		std::cout<<"mark--------3" <<std::endl;

		// store control points
		extern std::vector<Point3f > globalPointsToDraw1  ; 
		for( int i=0; i<comSolver.finalCtrPoints.size(); ++i )
			for( int j=0; j<comSolver.finalCtrPoints[i].size(); ++j )
				globalPointsToDraw1.push_back( comSolver.finalCtrPoints[i][j] ) ;


		// store sharp features
		extern std::vector<Point3f > globalPointsToDraw2  ;

		for( int i=0; i<comSolver.finalCtrPoints.size(); ++i )
			for( int k=0; k<comSolver.sharpIds[0].size(); ++k)
					globalPointsToDraw2.push_back( comSolver.finalCtrPoints[i][comSolver.sharpIds[0][k]] ) ;

	}



	generateMesh() ;

	return points ;
}


bool mycompfunc (std::pair<double,int> a,std::pair<double,int> b) ;

std::vector<completionSolver> skeleton_whole::generateJointSetting(){


	std::vector<completionSolver> solvers ;



	// if there is only one branch
	if( smoothedBrchPts.size() == 1 ){

		extern std::vector<Curve3D> glonalOriSkel ;
		bool firstIsTip= ( ReconstructorUtility::isTip(*this, smoothedBrchPts[0][0]  ) || ReconstructorUtility::isMidNode( glonalOriSkel, smoothedBrchPts[0][0]) );
		bool lastIsTip = ( ReconstructorUtility::isTip(*this, smoothedBrchPts[0].back()  ) || ReconstructorUtility::isMidNode(glonalOriSkel, smoothedBrchPts[0].back()));


		if( firstMustBeTip && !firstMustBeNotTip ) firstIsTip = true ;
		if( !firstMustBeTip && firstMustBeNotTip ) firstIsTip = false ;
		if( lastMustBeTip && !lastMustBeNotTip ) lastIsTip = true ;
		if( !lastMustBeTip && lastMustBeNotTip ) lastIsTip = false ;


		if( brch0IsLoop )
			firstIsTip = lastIsTip = false ;

		solvers.push_back(completionSolver(*this, 0, firstIsTip, lastIsTip, corres_cues)  ) ;
		return solvers ;
	}else{

		std::cout << "multi-branch skeleton!" <<std::endl;
		system("pause") ;
		
	}

	return solvers ;
}

/************************************************************************/ 
/*              skel_segment.cpp                                        */
/************************************************************************/
std::vector<int2> globalBrchNodeInJoint ;
void skeleton_whole::intialSegment( std::vector<Point3f> pts ) {

	// find tips
	getTips( tips) ;

	refineSeg_curves.clear() ;

	points = pts ;

	std::cout<<" skeleton_whole::intialSegment(){}" <<std::endl;
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





void skeleton_whole::refineSegment( std::vector<Point2f> curve , bool ctrlPressed) {
	std::cout<<" skeleton_whole::refineSegment(){}" <<std::endl;

	if( !ptMap2BrchNode.size() )
		return ;


	Curve2D curve2d = ReconstructorUtility::convertScreen2dToWorld2d(curve ) ;
	std::vector<Point2f> points2d = ReconstructorUtility::project3dPointsOntoScreen(points);
	CurveArray2D brchPts2d( smoothedBrchPts.size() ) ;
	for( int i=0; i<smoothedBrchPts.size(); ++i )
		brchPts2d[i] =  ReconstructorUtility::project3dPointsOntoScreen(smoothedBrchPts[i]);




	// if the curve is a loop, update brchScopes
	if( (curve[0] - curve.back()).Norm() < 20 ){

		// find brchPts both in joints and in cicle
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

		// ---------------- do refine
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

void skeleton_whole::addPointInJoints(std::vector<Point3f> ptscld,std::vector<Point3f> &spts,std::vector<bool> &ptsSelected ,std::vector<Point3f> allJoints,std::vector<Point3f> nml){
	std::cout<<" skeleton_whole::addPointInJoints(){}" <<std::endl;
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


void skeleton_whole::drawSegmentation( int pointsize ){

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


bool myint2compfunc1 (std::pair<int,int> a,std::pair<int,int> b) {    return ( a.first < b.first ) ; }
bool myint2compfunc2 (std::pair<int,int> a,std::pair<int,int> b) {    return ( a.second < b.second ) ; }
void skeleton_whole::drawBrchesInDiffColor( GLColor skelColor){

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


		if( appstate.displayProfileConfidence ){

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


	// --------------------- draw branch in joint 
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

void skeleton_whole::recenter(){

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




void skeleton_whole::computeBrchNodeInJoint() {


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






void skeleton_whole::compute3dJoints( ) {
	std::cout<<"skeleton_whole::compute3dJoints"<<std::endl;

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
			//int nid ;
			//if( ( smoothedBrchPts[bid][0] -  joints[i] ).Norm() <  ( smoothedBrchPts[bid].back() -  joints[i] ).Norm() )
			//	nid = 0 ;
			//else
			//	nid = smoothedBrchPts[bid].size()-1 ;

			//brchEndForJoints.back().push_back( int2(bid, nid )) ;


			if( ( smoothedBrchPts[bid][0] -  joints[i] ).Norm() < 0.05 )
				brchEndForJoints.back().push_back( int2(bid, 0 )) ;

			if( ( smoothedBrchPts[bid].back() -  joints[i] ).Norm() < 0.05 )
				brchEndForJoints.back().push_back( int2(bid, smoothedBrchPts[bid].size()-1 )) ;

		}
	}



	// ------------------------------- compute joints ---------------------
	calculateCutplanes() ;

	//fragments = possionFragments ;
	//return ;

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


/************************************************************************/
/*     skel_tips.cpp                                                    */
/************************************************************************/

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


void skeleton_whole::snapTips(std::vector<Point3f> ptcld , std::vector<Point3f> nmls) {

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

		totalNormal = ReconstructorUtility::ComputeEigenvector( normalAroundTip );
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

			//Point3f centerOfTip = Point3f(0,0,0);
			/*int n=0;*/
			//for( int pid=0; pid<fra.size(); ++pid ){
			//	double dis = ( totalNormal * ( fra[pid] - farthest ) ) / totalNormal.Norm();
			//	if(  fabs(dis) < 0.01 ){
			//		centerOfTip += fra[pid];
			//		edges.push_back(fra[pid]);
			//		n++;
			//	}
			//}
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

		//double bsize = smoothedBrchPts[bid].size()-1 ;
		//double deformNsize = bsize/2;

		//for( int j=0; j<smoothedBrchPts[bid].size(); ++j ){
		//	if( abs(nid-j) < deformNsize )
		//		smoothedBrchPts[bid][j] += trans * (deformNsize-abs(nid-j))/deformNsize ; 
		//}

		//for(int j=0; trans.Norm()>0.00001; j++){
		//	smoothedBrchPts[bid][abs(nid - j)] += trans;
		//	trans = trans / 2;
		//}
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

void skeleton_whole::snapFaceTips( std::vector<Curve3D>  impSkel) {

	// tips are initialized in intialSegment (){}
	if( smoothedBrchPts.size() > 1  ) 
		return ;

	std::vector<int2> faceTips ;
	if( !ReconstructorUtility::isTip(*this, smoothedBrchPts[0][0] ) &&!ReconstructorUtility::isMidNode(impSkel, smoothedBrchPts[0][0] )  )
		faceTips.push_back(int2(0,0)) ;
	if( !ReconstructorUtility::isTip(*this, smoothedBrchPts[0].back() )&& !ReconstructorUtility::isMidNode(impSkel, smoothedBrchPts[0].back() ) )
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
	if( !ReconstructorUtility::isTip(*this, smoothedBrchPts[0][0] ) &&!ReconstructorUtility::isMidNode(impSkel, smoothedBrchPts[0][0] )  )
		faceTips.push_back(int2(0,0)) ;
	if( !ReconstructorUtility::isTip(*this, smoothedBrchPts[0].back() )&& !ReconstructorUtility::isMidNode(impSkel, smoothedBrchPts[0].back() ) )
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

void skeleton_whole::addTipSeed( Point2f seed){

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


void skeleton_whole::removeTips( Point2f seed){

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


/************************************************************************/
/*        skel_updateResult.cpp                                         */
/************************************************************************/

std::vector<Point3f> skeleton_whole::updateResult(int newProfId, ON_NurbsCurve newProf ){

	std::cout << "skeleton_whole::updateResult" <<std::endl;

	std::vector<Point3f> points_;


	solvers[0].resolve(newProfId, newProf ) ;


	subpointsLast = subpoints ;
	points.clear();
	subpoints.clear() ;
	resultProf3d.clear() ;
	resultProf2d.clear() ;
	for( int iter =0; iter<solvers.size(); ++iter ){
		completionSolver &comSolver = solvers[iter];


		std::vector<cutPlane> cutpls = comSolver.cutplanes ;
		std::vector<Point3f> brch = comSolver.centers ;

		std::vector<std::vector<Point3f>> profiles3d ;
		std::vector<std::vector<Point2f>> profiles2d ;
		for( int i=0; i<comSolver.finalProfiles.size(); ++i ){

			Profile3D p3d =  ReconstructorUtility::convert2dProfileTo3d( comSolver.finalProfiles[i], cutpls[i], brch[i] ) ;
			profiles3d.push_back( p3d ) ;
			profiles2d.push_back( comSolver.finalProfiles[i] ) ;

		}

		resultProf3d.push_back( profiles3d ) ;
		resultProf2d.push_back( profiles2d ) ;

		std::vector<Point3f> subpoint = ReconstructorUtility::upsampleResult( profiles3d, int2(0, brch.size()-1 ), comSolver.firstIsTip, comSolver.lastIsTip ) ;
		points.insert(points.begin(), subpoint.begin(), subpoint.end() ) ;

		subpoints.push_back( subpoint ) ; 	// for blend

	}




	generateMesh() ;

	return points_ ;
}
