
#ifndef _skeleton_kx_h_ 
#define _skeleton_kx_h_

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
//#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/surface/3rdparty/opennurbs/opennurbs_nurbscurve.h>


typedef std::vector<Point2f> Profile2D ;
typedef std::vector<Point3f> Profile3D ;
typedef std::vector<PointPL> ProfilePL ;

#include "smtf.h"
#include "nrigidtransform.h"
#include "Types.h"
#include "triangleList.h"

struct cutPlane{
	Point3f cutRadius ;
	Point3f planeNormal ;
};

enum imgOutputType{ balckdot, withsharpness, withcompleteness };

class completionSolver ;

class skeleton_kx {

public:
	skeleton_kx(){ branchIdofTemplatetoDraw = 0;}
	skeleton_kx( std::vector<std::vector<Point3f>> & brchPts, bool dsample ) ;

	std::vector<Point3f> reconstruct() ;
	std::vector<Point3f> reconstruct_old() ;


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
	void initSharpLable() ;
	void discretizeNurbs() ;
	void convertDisNurbs2Pixles() ;
	void outProfilesAsImage( double len, int res, int imageWidth, imgOutputType OT = balckdot, int psize = 2 ) ;
	void cvtPixels2CV() ;
	void setSharpNurbsCV( int NurbsId, int CVId) ;

	void normalize( Box3f box ) ;
	void drawBrchesInDiffColor();
	void drawCutRadii() ;
	void drawTemplates() ;
	void drawNurbs() ;
	void draw3dNurbs() ;
	void drawGraph() ;
	void update() ;

	void snapTips() ;

	// ---------------- segment ------------------------- //
	void intialSegment(  std::vector<Point3f> pointsset ) ;
	void refineSegment( std::vector<Point2f> curve , bool ctrlPressed ) ;
	void compute3dJoints( ) ;
	void drawSegmentation( int pointsize ) ;
	void recenter() ;
	std::vector<Point3f> points ;
	std::vector<int2> ptMap2BrchNode ;
	std::vector<bool> isJoint ;
	std::vector<int2> brchNodeInJoint ;
	std::vector<int2> brchScopes ; // need initialization
	std::vector<std::vector<Point3f>> joints3d ;
	std::vector<std::vector<Point3f>> fragments ;
	std::vector<std::vector<int>> fragments_nid ;
	std::vector<std::vector<Point2f>> refineSeg_curves ;
	std::vector<std::vector<int2>> brchEndForJoints ;
	void computeBrchNodeInJoint() ;
	std::vector<completionSolver> generateJointSetting() ;
	int jointSettingNum ;
	// -------------------------------------------------- //


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

	// profiling planes --- //
	std::vector<std::vector<cutPlane>> cutplanes ;

	// profiles --- //
	std::vector<std::vector<Profile3D>> profiles3d_full ; 
	std::vector<std::vector<Profile3D>> profiles3d ; 
	std::vector<std::vector<Profile2D>> profiles2d ;
	std::vector<std::vector<std::vector<double>>> sharpness ;
	std::vector<std::vector<double>> profile_conf ;
	std::vector<double2> profconf_maxmin ;



	// final transformation of each profile
	std::vector<std::vector<ST>> FinalSimTrans ;

	// selected templates
	std::vector<std::vector<int>> tempalteIdss ;
	std::vector<int2> manualTemplates ;
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

	// tips
	std::vector<int2> tips ;
	std::vector<Point3f> tipSeed ;  // user selected dst for tips
	void addTipSeed(Point2f p) ;    // press alt + left button to add a seed point for tips
	void removeTips(Point2f p) ;	// press shift + left button to remove a tip
	void addTemplate(Point2f p) ;   // press ctrl + right button to add a template

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
	ON_NurbsCurve ctrTmp ;
	bool ctrTmpValid ;
	Box2f initBoxCtrTmp ;
	int2 ctrTmpIdx ;
	std::vector<Point2f> disCtrTmp2d;
	std::vector<Point3f> disCtrTmp3d;
	std::vector<Point2f> disCtrTmp2dPixel;
	std::vector<Point2f> ctrTmp2dCVPixel;
	std::vector<int> changedPwbranId ;
	void selectCtrTmp( Point2f px ) ;
	void moveCtrTmpCV( int cvid, Point2f pixelPos ) ;
	void drawCtrlTmp() ;
	std::vector<Point3f> updateResult( ) ;


	// get mesh
	std::vector<triangleList> resultMesh ;
	void generateMesh() ;


};

#include "completionSolver_v2.h" 

#endif