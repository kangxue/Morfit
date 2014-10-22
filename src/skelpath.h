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


#ifndef _skelpath_h_ 
#define _skelpath_h_

#include "Types.h"
#include<vector>
#include "CMesh.h"

#include <lemon/list_graph.h>

#include <pcl/surface/on_nurbs/fitting_curve_2d_pdm.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_tdm.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_sdm.h>
#include <pcl/surface/on_nurbs/triangulation.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>
#include <pcl/surface/3rdparty/opennurbs/opennurbs_nurbscurve.h>

#include "smtf.h"
#include "nrigidtransform.h"
#include "Types.h"
#include "triangleList.h"
class completionSolver ;

class skelpath {
	//this class is used for select one skeleton
public:
	skelpath(){ branchIdofTemplatetoDraw = 0;  iditifier = clock(); }
	skelpath( std::vector<std::vector<Point3f>> & brchPts, bool dsample ) ;
	skelpath( completionSolver solver ) ;

	unsigned iditifier ;

	std::vector<Point3f> reconstruct() ;


	// data structure utilities
	void convertVerticesListToGraph() ;
	void convertGraphToVerticeList() ;
	void convertVeticesListToBranchPoints() ;
	void convertBranchPointsToVeticesList(  double disToMerge ) ;

	void smoothEachBranch() ;

	void getJoints( std::vector<Point3f> &joints, std::vector<std::vector<int>> &ajBrchId ) ;
	
	void getTips( std::vector<int2> &endpoints) ;

	void calculateCutplanes() ;

	void calculateProfiles( ) ;

	void initAfterClaculateProfiles();

	void calculateSharpness() ;

	void calculateProfileConfidence() ;

	void selectTemplates() ;

	void convertTemplate2Pixles() ;

	void fitInitialNurbs() ;

	void discretizeNurbs() ;

	void convertDisNurbs2Pixles() ;

	void initSharpLable() ;


	void switchTemplate_online(Point2f p) ;

	void cvtPixels2CV() ;
	void setSharpNurbsCV( int NurbsId, int CVId) ;

	void normalize( Box3f box ) ;

	void drawBrchesInDiffColor(  GLColor skelColor );

	void drawCutRadii() ;

	void drawTemplates() ;

	void drawNurbs() ;

	void draw3dNurbs() ;
	
	void snapTips(std::vector<Point3f> ptcld, std::vector<Point3f> nmls) ;
	void snapFaceTips(std::vector<Curve3D>  impSkel  ) ;

	//   segment
	void intialSegment(  std::vector<Point3f> pointsset) ;
	void refineSegment( std::vector<Point2f> curve , bool ctrlPressed ) ;
	void addPointInJoints( std::vector<Point3f> pointset,std::vector<Point3f> &subPointSet,std::vector<bool> &ptsSelected ,std::vector<Point3f> allJoints,std::vector<Point3f> nml);

	void compute3dJoints( ) ;

	void drawSegmentation( int pointsize ) ;

	void recenter() ;

	//*** segement data members  ***
	std::vector<Point3f> points ;
	std::vector<int2> ptMap2BrchNode ;
	std::vector<bool> isJoint ;
	std::vector<int2> brchNodeInJoint ;
	std::vector<int2> brchScopes ; // need initialization
	std::vector<std::vector<Point3f>> joints3d ;
	std::vector<std::vector<Point3f>> fragments ;
	std::vector<std::vector<int>> fragMapId; 
	std::vector<std::vector<int>> fragments_nid ;
	std::vector<std::vector<Point2f>> refineSeg_curves ;
	std::vector<std::vector<int2>> brchEndForJoints ;
	void computeBrchNodeInJoint() ;
	std::vector<completionSolver> generateJointSetting() ;
	int jointSettingNum ;
	//***    end                 ***


	// ---- graph representation --- //
	lemon::ListGraph  *skeletonGraph ;
	lemon::ListGraph::NodeMap<Point3f> *graphNodeP3f;

	// ---- vertices list representation --- //
	std::vector<Point3f> nodes ;  
	std::vector<std::vector<int>> branches ; 

	// ---- branch points representation --- //
	std::vector<std::vector<Point3f>>  brchPts ;

	// degree of each points
	std::vector<std::vector<int>>  brchNodeDegree ;


	// ---- smoothed branch points, cannot used for calculus purpurse --- //
	std::vector<std::vector<Point3f>>  smoothedBrchPts ;

	// profiling planes  //
	std::vector<std::vector<cutPlane>> cutplanes ;

	// profiles  //
	std::vector<std::vector<Profile3D>> profiles3d_full ; 
	std::vector<std::vector<Profile3D>> profiles3d ; 
	std::vector<std::vector<Profile2D>> profiles2d;
	std::vector<std::vector<std::vector<double>>> sharpness ;
	std::vector<std::vector<double>> profile_conf ;
	std::vector<double2> profconf_maxmin ;



	// final transformation of each profile
	std::vector<std::vector<ST>> FinalSimTrans ;

	// selected templates
	std::vector<std::vector<int>> tempalteIdss ;
	std::vector<int2> flipTemplatesId ;
	std::vector<std::vector<ON_NurbsCurve>> Nurbsss ;
	std::vector<std::vector<std::vector<Point2f>>> disNurbsss;
	std::vector<std::vector<std::vector<Point3f>>> disNurbsss3d;
	//std::vector<std::vector<std::vector<int>>> KnotCorresId ;
	std::vector<std::vector<Point2f>> disNurbsssPixels;
	std::vector<std::vector<Point2f>> NurbsCVPixels;
	std::vector<std::vector<bool>> NurbsCVPixelsSharpLable;
	std::vector<std::vector<Point2f>> temPixels ;


	// reconstructed 3d profiles
	std::vector<std::vector<Profile3D>> resultProf3d;
	std::vector<std::vector<Profile2D>> resultProf2d;
	std::vector<std::vector<Point3f>> ptsOfEachBrchPair ;
	std::vector<std::vector<cutPlane>> cutplanes_bp;
	std::vector<std::vector<Point3f>> branches_bp;


	//std::vector<std::vector<Profile3D>> reconstructed3d ; 

	int branchIdofTemplatetoDraw ;
	void setActivBranch(int bid );

	// tips or end plane
	std::vector<int2> tips ;
	std::vector<Point3f> tipSeed ;   // user selected dst for tips
	void addTipSeed(Point2f p) ;     // press alt + left button to add a seed point for tips
	void removeTips(Point2f p) ;	 // press shift + left button to remove a tip
	void switchTemplate(Point2f p) ; // press ctrl + right button to add a template

	// blend
	std::vector<std::vector<Point3f>> subpoints ;
	std::vector<std::vector<Point3f>> subpointsLast ;
	std::vector<std::vector<Point3f>> pwBrches ;
	std::vector<std::vector<cutPlane>> pwBrchesCtpl ;
	std::vector<int2> pwBrchesCmpnt ;
	std::vector<int> pwBrcheJid ;
	void removeInnerPoints();
	std::vector<Point3f> blendBranchBetweenJoints();

	// for debug
	std::vector<std::vector<Point3f>> jointBranches ;


	// for edit result
	std::vector<completionSolver> solvers ;
	std::vector<int> changedPwbranId ;
	std::vector<Point3f> updateResult(int newProfId, ON_NurbsCurve newProf  ) ;


	// get mesh
	std::vector<triangleList> resultMesh ;
	void generateMesh() ;


	// edit profile on 3d plane
	int activeTmpId ;
	int activeCVID ;
	int cvSelectionType; 
	void selectActiveTmp( Point2f px ) ;
	void selectActiveTmpCV( Point2f px ) ;
	void moveActiveTmpCV(  Point2f pixelPos ) ;
	void deselectActiveTmpCV() ;
	void deselectActiveTmp() ;
	void drawActiveTmp() ;
	void drawActiveTmp_cutplane() ;
	void drawActiveTmpOnCorner() ;
	void drawSpecialView() ;

	// for correspondence detection
	std::vector<int2> corres_cues; // user intervention to help detect correspondeces.   profile ID, cv ID. 


	// for tip correcting
	bool firstMustBeTip ;
	bool firstMustBeNotTip ;
	bool lastMustBeTip ;
	bool lastMustBeNotTip; 


	// for loop
	bool brch0IsLoop ;

};

#include "completionSolver_v2.h" 

#else

class sweeper ;

#endif