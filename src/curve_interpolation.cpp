
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include <cv.h>
#include <highgui.h>

#include "E:\featureSnap\image\imgSalCurDetector\getCurves\LapVectorField.h"

#ifdef _DEBUG 
#pragma comment(lib, "opencv_core242d.lib") ;
#pragma comment(lib, "opencv_highgui242d.lib") ;
#pragma comment(lib, "opencv_imgproc242d.lib") ;
#else
#pragma comment(lib, "opencv_core242.lib") ;
#pragma comment(lib, "opencv_highgui242.lib") ;
#pragma comment(lib, "opencv_imgproc242.lib") ;
#endif

#include "E:\\SloppyPasting\\strucDrivenCompletion\\interactiveInpaint\\cubicInterpolation.h"
#include "E:\featureSnap\image\imgSalCurDetector\getCurves\declarations.h"

void curve_interplation( Edges &input, Edges &output ){

	std::vector<std::vector<CvPoint>> &ctrpoints = input.elems ;
	std::vector<std::vector<CvPoint>> &curves = output.elems ;
	curves.resize( ctrpoints.size() ) ;


	for( int i=0; i<ctrpoints.size(); ++i ){


		// cubic spline interpolation
		std::vector<double> x, y,t ;
		for( int j=0; j<ctrpoints[i].size(); ++j){
			x.push_back(ctrpoints[i][j].x);
			y.push_back(ctrpoints[i][j].y);
			t.push_back( (double)j/(ctrpoints[i].size()-1) ) ;
		}

		Spline3Interp poly1( t, x);
		poly1.calcCoefs();
		std::vector<double> px,py ;
		for( double t=0.0; t<1.0; t+=0.001)
			px.push_back(poly1.evaluate(t)) ;

		Spline3Interp poly2( t, y);
		poly2.calcCoefs();
		for( double t=0.0; t<1.0; t+=0.001)
			py.push_back(poly2.evaluate(t)) ;



		for( int id=0; id<px.size(); ++id)
			curves[i].push_back( cvPoint(px[id], py[id] ) );
		CvPoint prepoint = curves[i][0] ;

		// reduce
		for( int j=1; j<curves[i].size(); ++j ){

			if( sqrt( (double)((prepoint.x - curves[i][j].x) *(prepoint.x - curves[i][j].x) 
				+ (prepoint.y - curves[i][j].y) *(prepoint.y - curves[i][j].y)) ) >5  ){

					prepoint = curves[i][j] ;
			}
			else{
				curves[i].erase( curves[i].begin() + j ) ;
				j-- ;
			}

		}

		for( int j=1; j<curves[i].size(); ++j ){

			if( curves[i][j].x == curves[i][j-1].x &&  curves[i][j].y == curves[i][j-1].y ){

				curves[i].erase( curves[i].begin() + j ) ;
				j-- ;
			}

		}

		for( int j=0; j<curves[i].size(); ++j ){
			for( int k=j+1; k<curves[i].size(); ++k ){

				if( curves[i][j].x == curves[i][k].x &&  curves[i][j].y == curves[i][k].y ){

					curves[i].erase( curves[i].begin() + k ) ;
					k-- ;
				}
			}
		}



	}


}

void getBernsteinCoefficient( std::vector<double> &bc, int n   ){

	bc.resize(n);

	n=n-1;

	int j,k;
	for (k=0;k<=n;k++)
	{ //compute n! / (k!*(n-k)!)
		bc[k] = 1;
		for (j = n;j>=k+1;j--){
			bc[k] *= (double)j / (double)(j-k);
		}
	}



}

void getBezier( Edges &input, Edges &output, int step = 3 ){

	std::vector<std::vector<CvPoint>> &ctrpoints = input.elems ;
	std::vector<std::vector<CvPoint>> &curves = output.elems ;
	curves.resize( ctrpoints.size() ) ;


	for( int i=0; i<ctrpoints.size(); ++i ){

		// Coefficients of Bernstein Polynomials
		std::vector<double> bc ;


		const int n = ctrpoints[i].size() ;


		getBernsteinCoefficient(bc, n ) ;

		// store the Bezier curve
		curves[i].clear();
		for( double u=0; u<=1.0; u+=0.01 ){
			double2 p ;
			for( int k = 0; k<n; ++k ){
				p.x += bc[k] * pow(u, k) * pow( 1-u, n-1-k) * ctrpoints[i][k].x ;
				p.y += bc[k] * pow(u, k) * pow( 1-u, n-1-k) * ctrpoints[i][k].y ;
			}	

			curves[i].push_back( cvPoint(p.x,p.y) ) ;
		}

		CvPoint prepoint = curves[i][0] ;

		// reduce
		for( int j=1; j<curves[i].size()-1; ++j ){

			if( sqrt( (double)((prepoint.x - curves[i][j].x) *(prepoint.x - curves[i][j].x) 
				+ (prepoint.y - curves[i][j].y) *(prepoint.y - curves[i][j].y)) ) >step  ){

					prepoint = curves[i][j] ;
			}
			else{
				curves[i].erase( curves[i].begin() + j ) ;
				j-- ;
			}

		}


	}


}



double M[4][4] = {
	-1, 3, -3, 1,
	3,-6,3,0,
	-3,0,3,0,
	1,4,1,0
} ;

void getUniformCubicBSpline( std::vector<double2> &input, std::vector<double2> &output ){
//  refer http://en.wikipedia.org/wiki/B-spline

	output.clear() ;


	for( int i=1; i+2<input.size(); ++i  ){
		for( double t = 0; t<1.0; t+= 0.1 ){

			// index of control points


			double T[4] ; 
			for( int id = 0 ; id< 3; ++id )
				T[id] = pow( t, 3.0-id ) /6.0 ;
			T[3] = 1.0/6.0 ;


			double TM[4] ;
			for( int id = 0 ; id< 4; ++id ){
				TM[id] = 0;
				for( int j = 0; j<4; ++j)
					TM[id] += T[j] * M[j][id] ;
			}

			double2 St(0,0);
			for( int id = 0 ; id< 4; ++id ){
				St = St +  input[i-1 + id ] * TM[id] ;
			}

			output.push_back(St) ;

		}

	}


}


void getUniformCubicBSpline( Edges &input,Edges &output ){

	 std::vector<std::vector<CvPoint>> &ctrpoints = input.elems ;
	 std::vector<std::vector<CvPoint>> &curves = output.elems ;
	 curves.clear() ;
	 curves.resize( ctrpoints.size() ) ;


	 for( int i=0; i<ctrpoints.size(); ++i ){

		 std::vector<double2> ctrpoint_flt ; 
		 std::vector<double2> res ;
		 for( int j=0; j<ctrpoints[i].size(); ++j )
			 ctrpoint_flt.push_back( double2(ctrpoints[i][j].x, ctrpoints[i][j].y) ) ;
		 
		 getUniformCubicBSpline( ctrpoint_flt, res) ;


		 for( int j=0; j<res.size(); ++j )
			 curves[i].push_back( cvPoint( res[j].x,res[j].y)  ) ;
		 		
		 
		 CvPoint prepoint = curves[i][0] ;

		 // reduce
		 for( int j=1; j<curves[i].size(); ++j ){

			 if( sqrt( (double)((prepoint.x - curves[i][j].x) *(prepoint.x - curves[i][j].x) 
				 + (prepoint.y - curves[i][j].y) *(prepoint.y - curves[i][j].y)) ) >5  ){

					 prepoint = curves[i][j] ;
			 }
			 else{
				 curves[i].erase( curves[i].begin() + j ) ;
				 j-- ;
			 }

		 }

		 for( int j=1; j<curves[i].size(); ++j ){

			 if( curves[i][j].x == curves[i][j-1].x &&  curves[i][j].y == curves[i][j-1].y ){

				 curves[i].erase( curves[i].begin() + j ) ;
				 j-- ;
			 }

		 }

		 for( int j=0; j<curves[i].size(); ++j ){
			 for( int k=j+1; k<curves[i].size(); ++k ){

				 if( curves[i][j].x == curves[i][k].x &&  curves[i][j].y == curves[i][k].y ){

					 curves[i].erase( curves[i].begin() + k ) ;
					 k-- ;
				 }
			 }
		 }

	 }
}