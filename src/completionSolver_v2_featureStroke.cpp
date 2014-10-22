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



#include <fstream>
extern std::ofstream logout ;

#include "completionSolver_v2.h"
#include "reconstructionUtility.h"
#include "reconstructorPara.h"
#include <cmath>
#include "appstate.h"
#include "para.h"
void completionSolver::setFeatureStrokeCVId( int CVId, Curve2D stroke_w2d ,Curve3D stroke_l3d , int2 sharpScope) {

	std::cout << "completionSolver::setFeatureStrokeCVId( "<<CVId <<", ,("<<sharpScope.x<<","<<sharpScope.y<<")"<<std::endl;

	featureStrokeCVIds.push_back( CVId ) ;

	featureStroke_w2d.push_back( stroke_w2d) ;
	featureStroke_l3d.push_back( stroke_l3d) ;

	extern GLdouble globalModelviewMatrix[16] ;
	MVMatrix mvm ;
	for( int i=0; i<16; ++i )
		mvm.m[i] = globalModelviewMatrix[i] ;
	mvmtrices.push_back( mvm ) ;

	sharpStrokeScope.push_back(sharpScope) ;
}


std::vector<double> completionSolver::init_afterFeatureStroke( ) {

	int cvnum = ReconstructorPara::cvNum ;
	int xn = profiles2d.size() ;



	// init weights
	int strokeTimes = featureStrokeCVIds.size() ;
	assert( strokeTimes == featureStroke_w2d.size() && strokeTimes == mvmtrices.size() && strokeTimes ==  sharpStrokeScope.size() ) ;
	for( int i=0; i<strokeTimes; ++i ){

		int2 scope = sharpStrokeScope[i] ;

		for( int pid=0; pid<profiles2d.size(); ++pid ){

			if( pid>=scope.x && pid<=scope.y)
				WeightOfCV[pid][featureStrokeCVIds[i]] = ReconstructorPara::weight_sharpCV ; 
			else{

				double winsize = 10 ;
				double offset = std::min( abs(pid-scope.x), abs(pid-scope.y) ) ;
				double w = 1.0 + (ReconstructorPara::weight_sharpCV-1.0) * pow(2.71828, -offset*offset/(winsize*winsize) ) ;
				WeightOfCV[pid][featureStrokeCVIds[i]] = w ; 
			}
		}

	}


	return ReconstructorUtility::convertNRTransform2Vector( FinalNRTrans ) ;

	// init x0
	for( int i=0; i<strokeTimes; ++i ){
		int strokeCVID = featureStrokeCVIds[i] ;

		int2 scope = sharpStrokeScope[i] ;
		for( int pid=0; pid<profiles2d.size(); ++pid ){

			Point2f trans ;
			int id = pid ; 
			if( id < scope.x ) id = scope.x ;
			if( id > scope.y ) id = scope.y ;


			Curve2D cvs = ReconstructorUtility::convert3dProfileTo2d( finalCtrPoints[id],cutplanes[id], centers[id]) ;

			std::vector<std::vector<int>> nstId ;
			ReconstructorUtility::computeFlannKNN( featureStroke_l3d[i], std::vector<Point3f>(1, finalCtrPoints[id][strokeCVID]), nstId , 1) ;
	
			Point2f dst = ReconstructorUtility::convert3dProfileTo2d( std::vector<Point3f>( 1,featureStroke_l3d[i][nstId[0][0]]  ), cutplanes[id], centers[id])[0] ;
			Point2f src = ReconstructorUtility::convert3dProfileTo2d( std::vector<Point3f>( 1,finalCtrPoints[id][strokeCVID]  ), cutplanes[id], centers[id])[0] ;
			trans = dst-src ;

			for( int cvid=0; cvid<cvnum; ++cvid ){
				double winsize = 2 ;
				double offset = std::min( abs( cvid-strokeCVID ), cvnum - abs( cvid-strokeCVID ) )  ;
				double w =  pow(2.71828, -offset*offset/(winsize*winsize) ) ;
				Point2f dt = trans * w; 

				double verticalLength = dt *  Point2f(cvs[strokeCVID]).Normalize() ;
				double tangentLength =  dt * Point2f( cvs[strokeCVID].Y(), -cvs[strokeCVID].X()).Normalize() ;
				FinalNRTrans[pid].t[cvid] += Point2f(cvs[cvid]).Normalize() * verticalLength + Point2f( cvs[cvid].Y(), -cvs[cvid].X()).Normalize() * tangentLength ;

				//double h = 0.02 ;
				//double dis = ( finalCtrPoints[pid][cvid] - finalCtrPoints[pid][strokeCVID]).Norm()  ;
				//double w =  pow(2.71828, -dis*dis/(h*h) ) ;
				//FinalNRTrans[pid].t[cvid] += trans * w; 


			}

		}
	}

	return ReconstructorUtility::convertNRTransform2Vector( FinalNRTrans ) ;
}

double completionSolver::E_featrueStroke( std::vector<NRTSF> &x, int profId  ){

	int strokeTimes = featureStrokeCVIds.size() ;
	assert( strokeTimes == featureStroke_w2d.size() && strokeTimes == mvmtrices.size() && strokeTimes ==  sharpStrokeScope.size() ) ;

	double penSum = 0 ;

	for( int sid=0; sid<strokeTimes; ++sid ){
		int2 scope = sharpStrokeScope[sid] ;

		if( profId < scope.x || profId > scope.y )
			continue ;

		Point3f cvp = ReconstructorUtility::convert2dProfileTo3d( 
						ReconstructorUtility::getCvPts( ReconstructorUtility::deformNurbs(Tmps[0], x[profId], WeightOfCV[profId] ) ), 
						cutplanes[profId], centers[profId]  

						)   [ featureStrokeCVIds[sid] ] ;

		Point2f cvp_w2d = ReconstructorUtility::project3dPointsOntoScreen( std::vector<Point3f>(1,cvp), mvmtrices[sid].m )[0] ;

		penSum += ( featureStroke_w2d[sid][ ReconstructorUtility::NearstPoint( featureStroke_w2d[sid], cvp_w2d ) ] - cvp_w2d ).SquaredNorm() ;

	}


	return penSum ;
}



std::vector<ON_NurbsCurve> completionSolver::consolidateFinalNurbs( std::vector<ON_NurbsCurve> finalNurbs ) {

	// make the trajectory as smooth as possible by shifting the control points

	for( int i=1; i<finalNurbs.size(); ++i ){

		std::vector<Point2f> cv0 = ReconstructorUtility::getCvPts(finalNurbs[i-1]) ;
		std::vector<Point2f> cv1 = ReconstructorUtility::getCvPts(finalNurbs[i]) ;

		int bestOffset = 0 ;
		double minDis = 1e10 ;
		for( int offset=0; offset<cv1.size(); ++offset ){
			double dis = 0;
			for( int j=0; j<cv1.size(); ++j)
				dis += ( cv0[j] - cv1[(j+offset)%cv1.size()]).Norm() ;
			if( dis < minDis ){
				minDis = dis ;
				bestOffset = offset ;
			}
		}

		std::vector<Point2f> newcv1(cv1.size()) ;
		for( int j=0; j<cv1.size(); ++j)
			newcv1[j] = cv1[(j+bestOffset)%cv1.size()] ;


		ReconstructorUtility::setCvOfNurbs(finalNurbs[i], newcv1 ) ;
	}



	return finalNurbs ;
}

std::vector<ON_NurbsCurve> completionSolver::smoothFinalTrajectory( std::vector<ON_NurbsCurve> finalNurbs ) {


	// get control points
	std::vector<std::vector<Point3f>> ctrlPointss ;
	for( int id=0; id<finalNurbs.size(); ++id )
		ctrlPointss.push_back( ReconstructorUtility::convert2dProfileTo3d( ReconstructorUtility::getCvPts(finalNurbs[id]) , cutplanes[ id ], centers[ id ] )) ;

	for( int cvi=0; cvi<ReconstructorPara::cvNum; ++cvi ){

		Curve3D traj ;
		for( int i=0; i<ctrlPointss.size(); ++i )
			traj.push_back( ctrlPointss[i][cvi] ) ;

		Curve3D bspline, bezier;

		//if( completionSolverPara::smoothTrajactory  )
		//	ReconstructorUtility::getBezier(traj,bspline,1000) ;
		//else 
		//	ReconstructorUtility::getUniformCubicBSpline( traj, bspline, 1000 ) ;

		if( completionSolverPara::smoothTrajactory  )
			ReconstructorUtility::getNurbs(traj,bspline, 1000, completionSolverPara::trajectoryNurbsDegree ) ;
		else{
			ReconstructorUtility::getUniformCubicBSpline( traj, bspline, 1000 ) ;
			//bspline = traj ;
		}


		std::vector<std::vector<int>> nstId ;
		ReconstructorUtility::computeFlannKNN( bspline, traj, nstId, 1) ;

		for( int i=0; i<traj.size(); ++i )
			ctrlPointss[i][cvi] = traj[i] = bspline[ nstId[i][0] ] ;

	}


	for( int id=0; id<finalNurbs.size(); ++id )
		ReconstructorUtility::setCvOfNurbs(finalNurbs[id],  ReconstructorUtility::convert3dProfileTo2d( ctrlPointss[id] , cutplanes[ id ], centers[ id ] )) ;


	return finalNurbs ;
}