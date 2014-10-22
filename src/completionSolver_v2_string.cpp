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




#include "completionSolver_v2.h"
#include "reconstructionUtility.h"
#include "reconstructorPara.h"
#include <cmath>
#include "appstate.h"
#include "para.h"

#include <sstream>

std::string completionSolver::cvtToString() {

	std::ostringstream oss ;

	
	//  std::vector<double> finalX ;
	oss << "## std::vector<double> finalX" <<std::endl;
	oss << finalX.size() <<std::endl;
	for( int i=0; i< finalX.size(); ++i )
		oss << finalX[i] <<" " ;
	oss << std::endl;


	// std::vector<Profile2D> profiles2d ;  
	oss << "## std::vector<Profile2D> profiles2d" <<std::endl;
	oss << profiles2d.size() <<std::endl;
	for( int i=0; i<profiles2d.size(); ++i ){
		oss << profiles2d[i].size() <<std::endl;
		for( int j=0; j<profiles2d[i].size(); ++j )
			oss << profiles2d[i][j].X() << " " << profiles2d[i][j].Y() << std::endl;
	}

	
	// std::vector<cutPlane> cutplanes ;  
	oss << "## std::vector<cutPlane> cutplanes " <<std::endl;
	oss << cutplanes.size() <<std::endl;
	for( int i=0; i<cutplanes.size(); ++i ){
		oss << cutplanes[i].cutRadius.X() <<" "<< cutplanes[i].cutRadius.Y() <<" "<< cutplanes[i].cutRadius.Z()  <<std::endl;
		oss << cutplanes[i].planeNormal.X() <<" "<< cutplanes[i].planeNormal.Y() <<" "<< cutplanes[i].planeNormal.Z()  <<std::endl;
	}

	// std::vector<Point3f> centers ; 
	oss << "## std::vector<Point3f> centers " <<std::endl;
	oss << centers.size() <<std::endl;
	for( int i=0; i<centers.size(); ++i )
		oss << centers[i].X() <<" " << centers[i].Y() <<" " << centers[i].Z()  <<std::endl;
	


	// std::vector<int> templatesIds ;  
	oss << "## std::vector<int> templatesIds" <<std::endl;
	oss << templatesIds.size() <<std::endl;
	for( int i=0; i< templatesIds.size(); ++i )
		oss << templatesIds[i] <<" " ;
	oss << std::endl;



	// std::vector<ON_NurbsCurve> Tmps ;
	oss << "## std::vector<ON_NurbsCurve> Tmps" <<std::endl;
	oss << Tmps.size() <<std::endl;
	for( int i=0; i<Tmps.size(); ++i ){
		Curve2D cv2d = ReconstructorUtility::getCvPts( Tmps[i] ) ;
		for( int j=0; j<cv2d.size(); ++j )
			oss << cv2d[j].X() << " " << cv2d[j].Y() << std::endl;
	}


	// std::vector<int2> TScope ;
	oss << "## std::vector<int2> TScope" <<std::endl;
	oss << TScope.size() <<std::endl;
	for( int i=0; i<TScope.size(); ++i ){
		oss << TScope[i].X() << " " << TScope[i].Y() << std::endl;
	}


	// std::vector< std::vector<int> > sharpIds ;
	oss << "## std::vector< std::vector<int> > sharpIds" <<std::endl;
	oss << sharpIds.size() <<std::endl;
	for( int i=0; i<sharpIds.size(); ++i ){
		oss << sharpIds[i].size() <<std::endl;
		for( int j=0; j<sharpIds[i].size(); ++j )
			oss << sharpIds[i][j] << " " ;
		oss <<  std::endl;
	}
	
	
	// std::vector< std::vector<int> > cvIdx ;
	oss << "## std::vector< std::vector<int> > cvIdx" <<std::endl;
	oss << cvIdx.size() <<std::endl;
	for( int i=0; i<cvIdx.size(); ++i ){
		oss << cvIdx[i].size() <<std::endl;
		for( int j=0; j<cvIdx[i].size(); ++j )
			oss << cvIdx[i][j] << " " ;
		oss <<  std::endl;
	}


	// bool firstIsTip ; 
	oss << "##  bool firstIsTip" <<std::endl;
	oss << firstIsTip <<std::endl;


	//bool lastIsTip ;
	oss << "##  bool lastIsTip" <<std::endl;
	oss << lastIsTip <<std::endl;

	//bool isLoop ; 
	oss << "##  bool isLoop" <<std::endl;
	oss << isLoop <<std::endl;



	return oss.str() ;

}


completionSolver::completionSolver( std::string str) {

	// remove comment
	int commentPos ;
	while( (commentPos = str.find( "##")) !=  string::npos ){

		int commentEnd = str.find( "\n", commentPos) ;

		str.erase( str.begin() + commentPos, str.begin()+commentEnd ) ;

	}


	std::istringstream iss( str ) ;


	int t, n, m, u, v;
	double x,y ,z ;

	int cvnum = ReconstructorPara::cvNum ;

	////  std::vector<double> finalX ;
	//oss << "## std::vector<double> finalX" <<std::endl;
	//oss << finalX.size() <<std::endl;
	//for( int i=0; i< finalX.size(); ++i )
	//	oss << finalX[i] <<" " ;
	//oss << std::endl;
	finalX.clear() ;
	iss >> n ;
	for( int i=0; i<n; ++i ){
		iss >> x ;   finalX.push_back( x ) ;
	}



	//// std::vector<Profile2D> profiles2d ;  
	//oss << "## std::vector<Profile2D> profiles2d" <<std::endl;
	//oss << profiles2d.size() <<std::endl;
	//for( int i=0; i<profiles2d.size(); ++i ){
	//	oss << profiles2d[i].size() <<std::endl;
	//	for( int j=0; j<profiles2d[i].size(); ++j )
	//		oss << profiles2d[i][j].X() << " " << profiles2d[i][j].Y() << std::endl;
	//}
	
	iss >> n ;
	profiles2d.clear() ;profiles2d.resize(n) ;
	for( int i=0; i<n; ++i ){
		iss >> m ;
		for( int j=0; j<m; ++j ){
			iss >> x >> y ;
			profiles2d[i].push_back( Point2f(x,y) ) ;
		}
	}


	//// std::vector<cutPlane> cutplanes ;  
	//oss << "## std::vector<cutPlane> cutplanes " <<std::endl;
	//oss << cutplanes.size() <<std::endl;
	//for( int i=0; i<cutplanes.size(); ++i ){
	//	oss << cutplanes[i].cutRadius.X() <<" "<< cutplanes[i].cutRadius.Y() <<" "<< cutplanes[i].cutRadius.Z()  <<std::endl;
	//	oss << cutplanes[i].planeNormal.X() <<" "<< cutplanes[i].planeNormal.Y() <<" "<< cutplanes[i].planeNormal.Z()  <<std::endl;
	//}
	cutplanes.clear() ;
	iss >> n ;
	cutplanes.resize(n) ;
	for( int i=0; i<n; ++i ){
		iss >> x >> y >> z;   cutplanes[i].cutRadius = Point3f(x,y,z) ;
		iss >> x >> y >> z;   cutplanes[i].planeNormal = Point3f(x,y,z) ;
	}





	//// std::vector<Point3f> centers ; 
	//oss << "## std::vector<Point3f> centers " <<std::endl;
	//oss << centers.size() <<std::endl;
	//for( int i=0; i<centers.size(); ++i )
	//	oss << centers[i].X() <<" " << centers[i].Y() <<" " << centers[i].Z()  <<std::endl;
	centers.clear() ;
	iss >> n ;
	centers.resize(n) ;
	for( int i=0; i<n; ++i ){
		iss >> x >> y >> z;   centers[i] = Point3f(x,y,z) ;
	}



	//// std::vector<int> templatesIds ;  
	//oss << "## std::vector<int> templatesIds" <<std::endl;
	//oss << templatesIds.size() <<std::endl;
	//for( int i=0; i< templatesIds.size(); ++i )
	//	oss << templatesIds[i] <<" " ;
	//oss << std::endl;
	templatesIds.clear() ;
	iss >> n ;
	templatesIds.resize(n) ;
	for( int i=0; i<n; ++i ){
		iss >> t;   templatesIds[i] = t  ;
	}



	// std::vector<ON_NurbsCurve> Tmps ;
	//oss << "## std::vector<ON_NurbsCurve> Tmps" <<std::endl;
	//oss << Tmps.size() <<std::endl;
	//for( int i=0; i<Tmps.size(); ++i ){
	//	Curve2D cv2d = ReconstructorUtility::getCvPts( Tmps[i] ) ;
	//	for( int j=0; j<cv2d.size(); ++j )
	//		oss << cv2d[j].X() << " " << cv2d[j].Y() << std::endl;
	//}
	Tmps.clear() ;
	iss >> n ;
	Tmps.resize(n) ;
	for( int i=0; i<n; ++i ){

		Curve2D cv2d ;
		for( int j=0; j<cvnum; ++j ){
			iss >> x >> y ; cv2d.push_back( Point2f(x,y) ) ;
		}

		ON_NurbsCurve nbs= ReconstructorUtility::fitNURBSforProfile2d( profiles2d[0]) ;
		ReconstructorUtility::setCvOfNurbs( nbs, cv2d, std::vector<double>(cvnum,1.0)) ;
		Tmps[i] = nbs ;
	}



	//// std::vector<int2> TScope ;
	//oss << "## std::vector<int2> TScope" <<std::endl;
	//oss << TScope.size() <<std::endl;
	//for( int i=0; i<TScope.size(); ++i ){
	//	oss << TScope[i].X() << " " << TScope[i].Y() << std::endl;
	//}
	iss >> n ;
	TScope.clear() ; TScope.resize(n) ;
	for( int i=0; i<n; ++i){
		iss >> u >> v ;
		TScope[i] = int2(u,v) ;
	}



	//// std::vector< std::vector<int> > sharpIds ;
	//oss << "## std::vector< std::vector<int> > sharpIds" <<std::endl;
	//oss << sharpIds.size() <<std::endl;
	//for( int i=0; i<sharpIds.size(); ++i ){
	//	oss << sharpIds[i].size() <<std::endl;
	//	for( int j=0; j<sharpIds[i].size(); ++j )
	//		oss << sharpIds[i][j] << " " ;
	//	oss <<  std::endl;
	//}
	iss >> n ;
	sharpIds.clear() ; sharpIds.resize(n) ;
	for( int i=0; i<n; ++i ){
		iss >> m ;
		for( int j=0; j<m; ++j ){
			iss >> t ;
			sharpIds[i].push_back( t) ;
		}
	}


	//// std::vector< std::vector<int> > cvIdx ;
	//oss << "## std::vector< std::vector<int> > cvIdx" <<std::endl;
	//oss << cvIdx.size() <<std::endl;
	//for( int i=0; i<cvIdx.size(); ++i ){
	//	oss << cvIdx[i].size() <<std::endl;
	//	for( int j=0; j<cvIdx[i].size(); ++j )
	//		oss << cvIdx[i][j] << " " ;
	//	oss <<  std::endl;
	//}
	iss >> n ;
	cvIdx.clear() ; cvIdx.resize(n) ;
	for( int i=0; i<n; ++i ){
		iss >> m ;
		for( int j=0; j<m; ++j ){
			iss >> t ;
			cvIdx[i].push_back( t) ;
		}
	}



	//// bool firstIsTip ; 
	//oss << "##  bool firstIsTip" <<std::endl;
	//oss << firstIsTip <<std::endl;

	////bool lastIsTip ;
	//oss << "##  bool lastIsTip" <<std::endl;
	//oss << lastIsTip <<std::endl;

	////bool isLoop ; 
	//oss << "##  bool isLoop" <<std::endl;
	//oss << isLoop <<std::endl;

	iss >> firstIsTip >> lastIsTip >> isLoop ;




	// compute the rest componets


	FinalNRTrans = ReconstructorUtility::convertVector2NRTransform(finalX, cvnum);  
	WeightOfCV = std::vector<std::vector<double>>(profiles2d.size(),std::vector<double>(cvnum,1) ) ;


	finalProfiles.clear() ;
	finalProfiles.resize( profiles2d.size() ) ;

	finalCtrPoints.clear() ;
	finalCtrPoints.resize( profiles2d.size() ) ;

	finalNurbs.clear() ;
	finalNurbs.resize( profiles2d.size() ) ;


	for( int pid = 0; pid< profiles2d.size();  ++pid )
		finalNurbs[pid] = ReconstructorUtility::deformNurbs(Tmps[0], FinalNRTrans[pid], WeightOfCV[pid] ) ;
	finalNurbs = consolidateFinalNurbs( finalNurbs ) ;
	finalNurbs = smoothFinalTrajectory(finalNurbs) ;

	for( int pid = 0; pid< profiles2d.size();  ++pid ){
		finalProfiles[pid] = ReconstructorUtility::discretizeNurbs( finalNurbs[pid], ReconstructorPara::finalProfilesSampleNum ) ;
		finalCtrPoints[pid] = ReconstructorUtility::convert2dProfileTo3d( ReconstructorUtility::getCvPts( finalNurbs[pid] ), cutplanes[pid], centers[pid] );
	}




}