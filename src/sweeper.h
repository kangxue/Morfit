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


#ifndef _sweeper_h_
#define  _sweeper_h_



#include "skeleton_mul.h"
#include "skelpath.h"

class traceInfo{
public:
	std::vector<int> lastsubPSet ;


};

class GLArea ;

class sweeper{

public:

	sweeper(){}
	sweeper( std::vector<Point3f> &ptcld, std::vector<Point3f> &normal, std::vector<Curve3D> &skel ) ;
	sweeper( std::vector<Point3f> &ptcld, std::vector<Point3f> &normal, std::vector<Curve3D> &skel,std::vector<Point3f> &allTheJoints ,std::vector< std::vector<Point3f> > &branches) ;

	skeleton_mul oriSkel;
	std::vector<Point3f> allJoints;

	Curve2D ongoingSkelStroke ;
	//Curve3D ongoingSkelStroke3d ;

	void startSkelStroke(Point2f p ) ;
	void moveSkelStroke(Point2f p) ;
	bool endSkelStroke(Point2f p, float *mvmatrix) ;

	bool generateSetting() ;
	void reconLastStroke( ) ;    //controler.cpp(737): dataMgr.swp.reconLastStroke() ;

	// feature stroke
	Curve2D ongoingFeatureStroke ;
	Curve3D featureStroke3d ;
	int featureStrokeSkelId ;
	void startFeatureStroke(Point2f p ) ;
	void moveFeatureStroke(Point2f p) ;
	bool endFeatureStroke(Point2f p, float *mvmatrix) ;
	void performFeatureStroke() ;
	int getFeatureStrokeSkelId() ;
	void segment(std::vector<Point2f> segStroke, bool shiftPressed, bool ctrPressed, bool leftbutton ) ;


	std::vector<Curve2D> strokes ;
	std::vector<MVMatrix> mvmatrices ;// model view matrix when user draw the skeleton
	std::vector<Curve3D> branches ;
	std::vector<std::vector<int>> subPtsIds ;
	std::vector<skelpath> settings ;
	std::vector<triangleList> meshes ;
	std::vector<unsigned> meshLable ;
	
	skelpath ongoingSkel ;

	std::vector<int> lastsubPSet ;

	std::vector<Point3f> points ;
	std::vector<int> segLable ;
	std::vector<Point2f> points_w2d ;
	std::vector<bool> ptsSelected ;
	std::vector<Point3f> normals ;

	std::vector<std::vector<int>> kneibs ;

	std::vector<Curve3D>  impSkel ;



	// for update result
	int ctrlSkelId ;
	int ctrlProfId ;
	int ctrlCVId ;
	int ctrlCVSelectionType ;
	ON_NurbsCurve ctrNurbs ;
	void selectCtrlProf( Curve2D stroke ) ;
	void selectCtrlProfCV( Point2f px ) ;
	void moveCtrlProfCV(  Point2f pixelPos ) ;
	void drawCtrlNurbs() ;
	void drawCtrlNurbs_cutplane() ;
	void drawCtrlNurbsOnCorner() ;

	void startEditResultStroke( Point2f p ) ;
	void moveEditResultStroke( Point2f p ) ;
	void EndEditResultStroke( ) ;
	Curve2D editResultStroke ;

	void updateResult() ;


	triangleList mergedMesh ;
	void mergeMesh() ;

	void saveSettings() ;
	void loadSettings() ;


	// correspondences
	void generateCorrespondenceInterface() ;
	void recoverAfterCorrespondenceInterface() ;
	std::vector<int> trajSegsToDisplay ;
	std::vector<Curve3D> corresTrajectory ;

	std::vector<skelpath> settings_bk ;
	std::vector<triangleList> meshes_bk ;
	std::vector<unsigned> meshLable_bk ;

	//recenter the skeleton
	void recenterizeSkel( skeleton_mul &skel ) ;

	vector<Point3f> nml;

};

#else

class sweeper ;

#endif