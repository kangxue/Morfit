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
#include "Types.h"
#include <fstream>
extern std::ofstream logout ;

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <algorithm>    // std::sort

#include <GL/GL.h>

//#include <cv.h>
//#include <highgui.h>

#include <string> ;

#include "reconstructionUtility.h"
#include "reconstructorPara.h"

#include "declarations.h"



#include "reconstructorPara.h" 
#include "appstate.h"

void skelpath::calculateCutplanes(){

	//calculate and get the member variable "cutplanes"
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

	// mark,  smooth cut plane for loop
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
}



void skelpath::calculateProfiles(  ) {

	//calculate and get the member variable "profiles"
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


	// down sample
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

void skelpath::initAfterClaculateProfiles(){

	//calling the following functions to initialize


	//calculate the sharpness and save it in the member variable "sharpness",and normalize sharpness
	calculateSharpness() ;

	//calculate the profile confidence ,using ReconstructorUtility::evaluateConfident()
	//save it in the member variable "profile_conf",smooth and normalize sharpness
	calculateProfileConfidence() ;

	//select templates with the help of the result "profile_conf"
	selectTemplates() ;

	//get the profile of the template and convert it to pixels with " ReconstructorUtility::cvtProf2Pixels()"
	convertTemplate2Pixles() ;

	//calculate the nurbsss and save it in the member variable "Nurbsss"
	fitInitialNurbs() ;

	//discretize nurbs and convert it to 3d curve
	discretizeNurbs() ;

	//convert discretized nurbs to pixels
	convertDisNurbs2Pixles() ;

	//initialization the member variable "NurbsCVPixelsSharpLable"
	initSharpLable() ;
}


void skelpath::calculateProfileConfidence(){

	//calculate the profile confidence ,using ReconstructorUtility::evaluateConfident()
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


void skelpath::selectTemplates() {

	//select templates with the help of the result "profile_conf"
	int maxTnum = ReconstructorPara::maxTemplateNum ;
	tempalteIdss.clear() ;
	tempalteIdss.resize( profile_conf.size() ) ;
	for( int i=0; i<profile_conf.size(); ++i ){

		std::vector<double> profconf_i = profile_conf[i];

		int windsize = std::max( 10, (int)(15 * profconf_i.size()/100) );
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



}

void skelpath::drawCutRadii() {

	//draw the smoothed branched branch points , the cut plane and profile
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


void skelpath::calculateSharpness(){

	//calculate the sharpness and save it in the member variable "sharpness",and normalize sharpness

	std::cout<<"skelpath::calculateSharpness() begin"<<std::endl;
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
				std::sort (dis.begin(), dis.end(), GlobalFun::mycompfunc );   
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

				//if( k==0 ){
				//	for( int i=0; i<n;++i )
				//		logout <<"[ "<<DataPoints(0,i) <<" " << DataPoints(1,i) <<" ] " <<"," ;
				//	logout<<std::endl;


				//}
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

					//logout<<"after mean"<<std::endl;
					//for( int i=0; i<n;++i )
					//	logout <<"[ "<<DataPoints(0,i) <<" " << DataPoints(1,i) <<" ] " <<"," ;
					//logout<<std::endl;

				}

				// get the covariance matrix
				MatrixXd Covariance = MatrixXd::Zero(m, m);
				Covariance = (1 / (double) n) * DataPoints * DataPoints.transpose();
				//  std::cout << Covariance ;	

				if( k==0 ){

					//logout<<"Covariance"<<std::endl;
					//for( int i=0; i<2;++i )
					//	for( int j=0; j<2;++j )
					//		logout << Covariance(i,j) <<", ";
					//logout<<std::endl;

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
						//logout <<  i << "-eigenvalue " << eigenvalues(i) << std::endl;

						//logout << i << "-eigenvector " << eigenVectors.row(i) << std::endl;

						//logout << std::endl ;
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

					std::cout<<"skelpath::calculateSharpness() end"<<std::endl;
}

void skelpath::switchTemplate( Point2f seed ){


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


void skelpath::switchTemplate_online(Point2f p) {

	//switch template and initialization after calculate profiles
	switchTemplate(p) ;
	initAfterClaculateProfiles() ;
}