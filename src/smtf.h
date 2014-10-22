

#ifndef _similarityTransform_
#define _similarityTransform_

#include "CMesh.h"


namespace transformSolverPara{

	extern double angleGraStep ;
	extern double scalingGraStep ;
	extern double transGraStep  ;
}

typedef class similarityTransform{
public:

	similarityTransform():rotateAngle(0),scaling(1.0), translate(Point2f(0,0)) {}
	similarityTransform(double r, double s, Point2f t):rotateAngle(r),scaling(s), translate(t) {}

	double rotateAngle ;
	double scaling ;
	Point2f translate ;

	inline double operator%( similarityTransform & st1 ){
		double angleDiff = (rotateAngle - st1.rotateAngle) * transformSolverPara::transGraStep / transformSolverPara::angleGraStep;
		double  scalDiff = (scaling - st1.scaling) * transformSolverPara::transGraStep / transformSolverPara::scalingGraStep;
		double transDiff = (translate - st1.translate).Norm() ;

		return angleDiff*angleDiff +  scalDiff*scalDiff + transDiff*transDiff;
	}

	//double operator=( similarityTransform & st1  ){
	//	rotateAngle = st1.rotateAngle ;
	//	scaling = st1.scaling ;
	//	translate = st1.translate ;
	//}


}ST;



#endif