#include <fstream>
extern std::ofstream logout ;

#include "transformSolver.h"
#include "GlobalFunction.h"

#include "..\nlopt-2.3\nlopt.hpp"
#include "reconstructionUtility.h"

namespace transformSolverPara{
	int winside = 7; 

	double leastQuareWeight = 1.0e10;

	double lsw_incspeed = 1.5 ;

	double angleGraStep = 0.01 ;
	double scalingGraStep = 0.01 ;
	double transGraStep = 0.01 ;

	double angleUpperBoundary = 3.14 ;
	double angleLowerBoundary = -3.14 ;

	double scalingUpperBoundary = 4.0 ;
	double scalingLowerBoundary = 0.2 ;


	double transUpperBoundary = 0.1 ;
	double transLowerBoundary = -0.1 ;


	bool diagnosticMessage = true ;

}


//std::vector<Point2f> transformSolver::transformProfile( int proId,  ST stf) {
//
//	std::vector<Point2f> res( profiles2d[proId ] ) ;
//
//	double theta = stf.rotateAngle ;
//	for( int i=0; i< profiles2d[proId].size(); ++i ){
//
//		double x = profiles2d[proId][i][0] ;
//		double y = profiles2d[proId][i][1] ;
//
//		res[i] = Point2f( x*cos(theta) - y*sin(theta) ,  x*sin(theta) + y*cos(theta) )  * stf.scaling + stf.translate ;
//	}
//	return res ;
//}

double transformSolver::profDis( std::vector<Point2f> &prof1, std::vector<Point2f> &prof2 , int id1, int id2  ) {

	int size1 = prof1.size() ;
	int size2 = prof2.size() ;

	std::vector<std::vector<int>> neiId ;
	
	std::vector<double> dis1to2( size1 ) ;

	GlobalFun::computeKNN(prof2, prof1, neiId, 1 ) ;
	for( int i=0; i<size1; ++i )
		dis1to2[i] = ( prof1[i] - prof2[ neiId[i][0] ] ).Norm() ;

	//std::vector<double> dis2to1(  size2  ) ;

	//GlobalFun::computeKNN(prof1, prof2, neiId, 1 ) ;
	//for( int i=0; i<size2; ++i )
	//	dis2to1[i] = ( prof2[i] - prof1[ neiId[i][0] ] ).Norm() ;

	double dis = 0;
	for( int i=0; i<dis1to2.size(); ++i )
		dis += dis1to2[i];
	//for( int i=0; i<dis2to1.size(); ++i )
	//	dis += dis2to1[i] * sharpness[id2][i] ;

	return dis ;
}

double transformSolver::objFunc(  std::vector<ST> &stf , double LSWeight) {


	int winside = transformSolverPara::winside ;

	double sigma = 0.3 * ( winside / 2.0 - 1.0 )  + 0.8  ;


	int npf = profiles2d.size() ;

	std::vector<std::vector<Point2f>> tsfPros( npf ) ;  // transformed profiles

	for( int i=0; i<npf; ++i )
		tsfPros[i] = ReconstructorUtility::transformProfile( profiles2d[currentTemplateId], stf[i] ) ;

	double totalDis = 0.0 ;
	for( int i=0; i<npf; ++i ){
		totalDis += profDis(profiles2d[i],tsfPros[i],  i, i )  ;
	}

	
	double leastSquare = 0.0 ;
	for( int i=0; i<npf; ++i )
		for( int offset = -winside/2 ;  offset<winside/2; ++offset ) {
			if( i+offset>=0 &&  i+offset<npf && offset!=0){
				double gaussWeight = pow( 2.71828, (-offset*offset/(2*sigma*sigma) ) ) ;
				leastSquare += (stf[i]%stf[i+offset])  *  gaussWeight ;
			}

		}


	if( transformSolverPara::diagnosticMessage) {
		
		std::vector<double> x = ReconstructorUtility::convertTransform2Vector( stf ) ;
		
		std::cout << "------------------------   ObjFunc -------------------------------\n x = [ " ;
		for( int i=0; i<x.size(); ++i ){
			std::cout <<x[i] <<" " ;
			if( i%4 == 3 ) std::cout <<"\n" ;
		}
		std::cout<<"]"<<std::endl;

		std::cout <<" totalDis = " <<totalDis <<",  leastSquare = " << leastSquare <<std::endl;
		std::cout << "f(x) = " << totalDis + leastSquare * LSWeight <<std::endl;
		std::cout << "------------------------------------------------------- "<<std::endl ;


		logout << "------------------------   ObjFunc -------------------------------\n x = [ " ;
		for( int i=0; i<x.size(); ++i ){
			logout <<x[i] <<" " ;
			if( i%4 == 3 ) logout <<"\n" ;
		}
		logout<<"]"<<std::endl;
		logout <<" totalDis = " <<totalDis <<",  leastSquare = " << leastSquare <<std::endl;
		logout << "f(x) = " << totalDis + leastSquare * LSWeight <<std::endl;
		logout << "-------------------------------------------------------"<<std::endl ;

	}

	return totalDis + leastSquare * LSWeight ;
}


double vfunc_ts(const std::vector<double> &x, std::vector<double> &grad, void* f_data){


	//if( transformSolverPara::diagnosticMessage) {
	//	std::cout << "-------------------------------------------------------\n x = [ " ;
	//	for( int i=0; i<x.size(); ++i )
	//		std::cout <<x[i] <<" " ;
	//	std::cout<<"]"<<std::endl;


	//	logout << "-------------------------------------------------------\n x = [ " ;
	//	for( int i=0; i<x.size(); ++i )
	//		logout <<x[i] <<" " ;
	//	logout<<"]"<<std::endl;

	//}


	transformSolver *tfs = (transformSolver *)f_data ;
	double fx = tfs->myvfunc( x, grad, tfs ) ;


	//if( transformSolverPara::diagnosticMessage) {
	//	std::cout << "------------------------------------------------------- "<<std::endl ;
	//	logout << "------------------------------------------------------- "<<std::endl ;
	//}

	return  fx;
}


double transformSolver::gradient( std::vector<ST> &stf, int xdiff_id,  double LSWeight ) {

	// init some parameters
	int winside = transformSolverPara::winside ;
	double sigma = 0.3 * ( winside / 2.0 - 1.0 )  + 0.8  ;
	int npf = profiles2d.size() ;
	int dstProfileId = xdiff_id/4 ;

	// transformation of x0 , and x1 on i
	ST st0 = stf[dstProfileId] ;
	ST st1 = st0;
	double step ;
	if     ( xdiff_id%4 == 0 ) { step = transformSolverPara::angleGraStep ; st0.rotateAngle -= step ; st1.rotateAngle += step ; }  
	else if( xdiff_id%4 == 1 ) { step = transformSolverPara::scalingGraStep ; st0.scaling -= step ;  st1.scaling += step ; }
	else { step = transformSolverPara::transGraStep ;   st0.translate[xdiff_id%2] -= step;  st1.translate[xdiff_id%2] += step; } 

	// profiles after transformation
	std::vector<Point2f> dstprofile0 = ReconstructorUtility::transformProfile( profiles2d[currentTemplateId],  st0 )  ;
	std::vector<Point2f> dstprofile1 = ReconstructorUtility::transformProfile( profiles2d[currentTemplateId],  st1 )  ;



	double totalDis = profDis(profiles2d[dstProfileId] , dstprofile0, dstProfileId, currentTemplateId )   ;

	double leastSquare = 0.0 ;
	for( int offset = -winside/2 ;  offset<winside/2; ++offset ) {
		if( dstProfileId+offset>=0 &&  dstProfileId+offset<npf &&  offset!= 0){
			double gaussWeight = pow( 2.71828, (-offset*offset/(2*sigma*sigma) ) ) ;
			leastSquare += (st0%stf[dstProfileId+offset])  *  gaussWeight ;
		}
	}
	double sum0 = totalDis + LSWeight * leastSquare ;


	totalDis = profDis(profiles2d[dstProfileId] , dstprofile1, dstProfileId, currentTemplateId )   ;
	leastSquare = 0.0 ;
	for( int offset = -winside/2 ;  offset<winside/2; ++offset ) {
		if( dstProfileId+offset>=0 &&  dstProfileId+offset<npf &&  offset!= 0){
			double gaussWeight = pow( 2.71828, (-offset*offset/(2*sigma*sigma) ) ) ;
			leastSquare += (st1%stf[dstProfileId+offset])  *  gaussWeight ;
		}
	}

	double sum1 = totalDis + LSWeight * leastSquare ;


	//if( transformSolverPara::diagnosticMessage) {

	//	std::cout << "g["<<xdiff_id<<"] = " << (sum1-sum0) / step <<"\n";
	//	logout << "g[ "<<xdiff_id<<" ] = " << (sum1-sum0) / step <<"\n";

	//}

	return (sum1-sum0) / (step*2)  ;

}

double transformSolver::myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){

	std::vector<ST> stf = ReconstructorUtility::convertVector2Transform(x) ;

	double fx = objFunc(stf,currentLSWeight) ;

	if (!grad.empty()) {
		
		//if( transformSolverPara::diagnosticMessage){
		//	std::cout << "------------------------   Gradient -------------------------------\n  " ;
		//	logout << "------------------------   Gradient -------------------------------\n  " ;
		//}

		for( int i=0; i<x.size(); ++i )
			grad[i] =  gradient( stf, i, currentLSWeight) ;

		//if( transformSolverPara::diagnosticMessage){
		//	std::cout << "-------------------------------------------------------"<<std::endl ;
		//	logout <<  "------------------------------------------------------- " <<std::endl ;
		//}

	}

	return fx;
}

#include "reconstructorPara.h" 
std::vector<ST> transformSolver::oneIterationSolver(  const std::vector<ST> &stf,  double LSWeight ) {
		
	// write least square weight
	currentLSWeight = LSWeight ;

	int pn = profiles2d.size() ;

	// convert stf to vector
	std::vector<double> x0 = ReconstructorUtility::convertTransform2Vector(stf ) ;

	int xn = x0.size() ;


	nlopt::opt opt(nlopt::LD_MMA, xn);

	std::vector<double> lb(xn);
	std::vector<double> ub(xn);
	for( int i=0; i<xn; ++i ){
		if( i%4 == 0 ){
			lb[i] = transformSolverPara::angleLowerBoundary ;
			ub[i] = transformSolverPara::angleUpperBoundary ;
		} else if(  i%4 == 1 ) {
			lb[i] = transformSolverPara::scalingLowerBoundary ;
			ub[i] = transformSolverPara::scalingUpperBoundary ;
		}else{
			lb[i] = transformSolverPara::transLowerBoundary ;
			ub[i] = transformSolverPara::transUpperBoundary ;
		}
	}

	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);

	opt.set_min_objective(vfunc_ts, this );


	opt.set_ftol_abs(ReconstructorPara::MMA_FTOL);

	double minf ;
	nlopt::result result ;
	
	try{
		result = opt.optimize(x0, minf);
	}catch (std::exception& e){  std::cerr << "exception caught: " << e.what() << '\n'; }

	std::cout << "result code: " << result <<std::endl;

	// convert result vector to transform
	return ReconstructorUtility::convertVector2Transform(x0)  ;


}


void transformSolver::solve( double LSW, int templateId ) {



	// if caller does not supply a positive templateId, choose the most confident template
	if( templateId < 0 ){
		int maxId = 0 ; 
		double maxConf = 0 ;
		for( int i=0; i<profile_conf.size(); ++i )
			if( profile_conf[i] > maxConf ) { maxConf=profile_conf[i] ;  maxId =i ; }
		templateId = maxId ;
	}

	// write templateId, to make objfunc can read it
	currentTemplateId = templateId ;


	// enlarge the confidence of template to make sure it match it self
	double conf_bk = profile_conf[templateId] ;
	profile_conf[templateId] = 2.0 ;


	std::vector<double> x0( profiles2d.size() * 4 ) ;
	for( int i=0; i<x0.size() ; ++i )
		if( i%4 == 1 )  x0[i] = 1.0 ; else x0[i] = 0.0 ;

	std::vector<ST> res = oneIterationSolver(  ReconstructorUtility::convertVector2Transform(x0), LSW ) ;

	x0 =  ReconstructorUtility::convertTransform2Vector(res ) ;
	
	//std::cout <<"result:\n" ;
	//for( int i=0; i<x0.size(); ++i ){
	//	std::cout << x0[i] <<" ";
	//}

	//system("pause") ;

	//res = oneIterationSolver( res, 1.0 ) ;

	FinalSimTrans = res ;

	x0 =  ReconstructorUtility::convertTransform2Vector(FinalSimTrans ) ;


	std::cout <<"result:\n" ;
	for( int i=0; i<x0.size(); ++i ){
		std::cout << x0[i] <<" ";
	}
	std::cout<<std::endl;


	profile_conf[templateId] = conf_bk ;



}