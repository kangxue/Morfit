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
#include "skeleton_mul.h"
#include "bspline.h"

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

void skeleton_mul::calculateCutplanes(){


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
}
