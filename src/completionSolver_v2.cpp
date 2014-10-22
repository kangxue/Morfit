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
#include "GlobalFunction.h"

#include "..\nlopt-2.3\nlopt.hpp"
#include "reconstructionUtility.h"
#include "reconstructorPara.h"
#include <cmath>
#include "declarations.h"
#include "omp.h"
#include "para.h"
#include "appstate.h"


completionSolver::completionSolver(const skelpath &skel, int branId, bool isTip_first, bool IsTip_last, std::vector<int2> cc ){ 
	profiles2d = skel.profiles2d[branId]; //sharpness = skel.sharpness[branId] ; 
	cutplanes = skel.cutplanes[branId];   centers = skel.smoothedBrchPts[branId] ;
	//profile_conf = skel.profile_conf[branId] ; 
	templatesIds = skel.tempalteIdss[branId] ;  Tmps = skel.Nurbsss[branId] ;

	firstIsTip = isTip_first ;
	lastIsTip = IsTip_last ;

	corres_cues = cc ;

	// sort template Ids
	std::vector<int> newTemplatesIds = templatesIds;     
	std::vector<ON_NurbsCurve> newTmps( Tmps.size() ); 
	std::sort( newTemplatesIds.begin(), newTemplatesIds.end() ) ;
	for( int i=0; i<newTemplatesIds.size(); ++i )
		for( int j=0; j<templatesIds.size(); ++j )
			if( templatesIds[j] == newTemplatesIds[i] )
				newTmps[i] = Tmps[j] ;
	templatesIds = newTemplatesIds ;
	Tmps = newTmps ;


	// set sharpIds
	sharpIds.clear() ;
	for( int i=0; i<Tmps.size(); ++i )
		sharpIds.push_back( ReconstructorUtility::getSharpCVofNurbs(Tmps[i]) ) ;



	// set Tscope
	TScope.resize( Tmps.size() ) ;
	for( int Tid=0; Tid<TScope.size(); ++Tid ){
		// compute sweeping scope of TEMPLATE Tid
		int Tnum =  templatesIds.size()   ;
		int Pid_begin, Pid_end ;

		if( Tid == 0 ) 	Pid_begin = 0 ;
		//else Pid_begin = ( templatesIds[Tid-1] + templatesIds[Tid] ) / 2;
		//else Pid_begin =  templatesIds[Tid-1];
		else Pid_begin = ( templatesIds[Tid-1] + templatesIds[Tid] ) / 2;

		if( Tid == Tnum-1) Pid_end = profiles2d.size()-1 ;
		//else Pid_end = ( templatesIds[Tid+1] + templatesIds[Tid] ) / 2;
		//else Pid_end = templatesIds[Tid+1];
		else Pid_end = ( templatesIds[Tid+1] + templatesIds[Tid] ) / 2 - 1;

		TScope[Tid].X() = Pid_begin ;
		TScope[Tid].Y() = Pid_end ;

	}

	// initialize weight of CV
	int cvNum = ReconstructorPara::cvNum ;
	std::vector<std::vector<double>>  w_cv( profiles2d.size(), std::vector<double>( cvNum,0 ) ) ;

	for( int pid=0; pid<profiles2d.size(); ++pid ){
		for( int cvid =0; cvid<cvNum; ++cvid ){
			ON_4dPoint p;  Tmps[0].GetCV(cvid, p ) ;
			w_cv[pid][cvid] = p.w ;
		}
	}

	WeightOfCV = w_cv ;

	tmpConsWeights.resize(Tmps.size() ) ;
	for( int i=0;i<Tmps.size(); ++i)
		tmpConsWeights[i] = completionSolverPara::w_tmpconstraint ;

	wdataScaled = false ;

	// detect loop

	isLoop = skel.brch0IsLoop ;
} 

int WhichStage = 1 ;


bool count_time = 0 ;
double timecount1 = 0.0 ;
double timecount2 = 0.0 ;
double timecount3 = 0.0 ;
double timecount4 = 0.0 ;
double timecount5 = 0.0 ;

double lastclock ;

	
Point2f centerOfProfile;
Point2f centerOfx;
Profile2D pointForEval;
float radiusOfProfile;


double vfunc_nrts(const std::vector<double> &x, std::vector<double> &grad, void* f_data){

	completionSolver *tfs = (completionSolver *)f_data ;
	double fx = tfs->myvfunc( x, grad, tfs ) ;

	return  fx;
}

double completionSolver::myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {

	int cvNum = ReconstructorPara::cvNum  ;
	int nrtsf_size = ( ReconstructorPara::cvNum * 2 + 3 ) ;
	
	double e =  objFunc( x);
	
	if( !grad.empty()){

		timecount1 = timecount2 = timecount3 = timecount4 = timecount5 = 0;
		
		std::vector<double> origEnergies( profiles2d.size() ) ;

	
#ifdef _use_openmp_
#define TNUM _Solver_TNUM_
			#pragma omp parallel num_threads(TNUM)
			{
				int i = omp_get_thread_num();

				for( int j=0; j<x.size()/TNUM+1; ++j ){
					int xid = j*TNUM+i  ;
					if( xid< x.size() )
						grad[xid] = gradient( x, xid, origEnergies[xid/nrtsf_size] ) ;
				}

			}

#else
			for( int i=0; i<x.size(); ++i ){
				grad[i] = gradient( x, i) ;
			}
		
#endif


		if( count_time )
			std::cout<<"grad time: "<<timecount1/CLOCKS_PER_SEC<<" "<<timecount2/CLOCKS_PER_SEC<<" "<<timecount3/CLOCKS_PER_SEC<<" "<<timecount4/CLOCKS_PER_SEC<<" "<<timecount5/CLOCKS_PER_SEC<<std::endl;
	
		if( completionSolverPara::diagnosticMessage) {
			double norm = 0.0 ;
			for( int i=0; i<grad.size(); ++i )
				norm += grad[i]*grad[i] ;
			std::cout<<"norm(g) = " << sqrt( norm/grad.size() ) <<std::endl;
			std::cout<<std::endl;
		}	
	}
	return e ;
}

std::vector<double>  completionSolver::oneIterationSolver(  std::vector<double> x0  ) {



	oneIterationSolverCalled = true ;

	lastX = x0 ;

	int pNum = profiles2d.size() ;
	int cvNum = ReconstructorPara::cvNum ;

	int nrtsf_size = ( ReconstructorPara::cvNum * 2 + 3 ) ;

	int xn = x0.size() ;
	nlopt::opt opt(nlopt::LD_MMA, xn);

	std::vector<double> lb(xn);
	std::vector<double> ub(xn);
	for( int i=0; i<xn; ++i ){
		
		int innerId = i % nrtsf_size ;
		
		if( innerId <cvNum*2 ){
			lb[i] = completionSolverPara::t_lb ;
			ub[i] = completionSolverPara::t_ub ;
		}else if( innerId == cvNum*2 ){
			lb[i] = completionSolverPara::A_lb ;
			ub[i] = completionSolverPara::A_ub ;
		}else{
			lb[i] = completionSolverPara::b_lb ;
			ub[i] = completionSolverPara::b_ub ;
		}

	}

	// check x0 
	for( int i=0; i<x0.size(); ++i )
		if( _isnan(x0[i]) || !_finite( x0[i]) || x0[i]<lb[i] || x0[i]>ub[i] ){
			std::cout << "bad x0"<<std::endl;
			std::cout << "x0["<<i<<"] = " << x0[i] <<std::endl;
			system("pause") ;
		}


	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);

	opt.set_min_objective(vfunc_nrts, this );


	opt.set_ftol_rel(ReconstructorPara::MMA_FTOL);

	double minf ;
	nlopt::result result ;

	try{
		result = opt.optimize(x0, minf);
	}catch (std::exception& e){  std::cerr << "exception caught: " << e.what() << '\n'; }

	std::cout << "result code: " << result <<std::endl;

	// convert result vector to transform
	return x0 ;
}


std::vector<double>  completionSolver::solveOneSegment(  std::vector<double> x, int lfPfid  ) {

	std::cout << "Left Profile ID = " << lfPfid <<std::endl;


	oneIterationSolverCalled = false ;

	int pNum = profiles2d.size() ;
	int cvNum = ReconstructorPara::cvNum ;
	int nrtsf_size = ( ReconstructorPara::cvNum * 2 + 3 ) ;


	// information used in objective function
	LeftProfId = lfPfid ;
	lastX = x ;


	int nT = Tmps.size() ;

	int pid0 = LeftProfId<0?0: templatesIds[LeftProfId];
	int pid1 = LeftProfId+1<nT? templatesIds[LeftProfId+1] :  profiles2d.size();

	int xn = (pid1 - pid0) *   nrtsf_size;


	// extract x0
	int startX = pid0 * nrtsf_size ;
	std::vector<double> x0(xn) ;
	for( int i=0; i<xn; ++i )
		x0[i] = x[startX+i] ;

	// solver one segment 
	nlopt::opt opt(nlopt::LD_MMA, xn);
	std::vector<double> lb(xn);
	std::vector<double> ub(xn);
	for( int i=0; i<xn; ++i ){
		int innerId = i % nrtsf_size ;
		if( innerId <cvNum*2 ){
			lb[i] = completionSolverPara::t_lb ;
			ub[i] = completionSolverPara::t_ub ;
		}else if( innerId == cvNum*2 ){
			lb[i] = completionSolverPara::A_lb ;
			ub[i] = completionSolverPara::A_ub ;
		}else{
			lb[i] = completionSolverPara::b_lb ;
			ub[i] = completionSolverPara::b_ub ;
		}
	}

	opt.set_lower_bounds(lb);
	opt.set_upper_bounds(ub);
	opt.set_min_objective(vfunc_nrts, this );
	opt.set_ftol_rel(ReconstructorPara::MMA_FTOL);

	double minf ;
	nlopt::result result ;

	try{
		result = opt.optimize(x0, minf);
	}catch (std::exception& e){  std::cerr << "exception caught: " << e.what() << '\n'; }

	//std::cout << "result code: " << result <<std::endl;

	// convert result to x of whole branch
	for( int i=0; i<xn; ++i )
		x[startX+i] = x0[i] ;

	return x ;
}




std::vector<double> completionSolver::morphing( ) {


	// modified for loop

	// we have only the first seed curve to transform , so there is  no need to morph weight

	int cvnum = ReconstructorPara::cvNum ;
	int xn = profiles2d.size() ;

	if( Tmps.size()==1 ){

		// initialize position
		std::vector<NRTSF> x(xn, NRTSF(cvnum) ) ;
		return ReconstructorUtility::convertNRTransform2Vector(x) ;
	}

	for( int i=0; i<corres.size(); ++i )
		if( corres[i].size() <cvnum ){ std::cout<<"correspondence not enough for morphing!"<<std::endl;   system("pause");  exit(1) ; } 

	//  position and weight
	std::vector<NRTSF> x(xn, NRTSF(cvnum) ) ;
	std::vector<std::vector<double>>  w_cv( xn, std::vector<double>(cvnum,0 ) ) ;


	// initialize x for each template
	std::vector<NRTSF> tmpX( Tmps.size() ) ;
	for( int tid=0; tid<Tmps.size(); ++tid ){
		tmpX[tid] = NRTSF(cvnum) ;
		for(int i=0; i<cvnum; ++i){
			ON_4dPoint p1,p2;
			Tmps[tid].GetCV( cvIdx[tid][i] , p1 ) ; p1.x/=p1.w; p1.y/=p1.w ;
			Tmps[0].GetCV( i, p2 ) ;  p2.x/=p2.w; p2.y/=p2.w ;
			tmpX[tid].t[i].X() = p1.x - p2.x;
			tmpX[tid].t[i].Y() = p1.y - p2.y;

		}
	}
	// initialize w_cv for each template
	std::vector<std::vector<double>> tmpW( Tmps.size() ) ;
	for( int tid=0; tid<Tmps.size(); ++tid ){
		tmpW[tid] = std::vector<double>(cvnum,0 ) ;
		for(int i=0; i<cvnum; ++i){
			ON_4dPoint p1,p2;
			Tmps[tid].GetCV( cvIdx[tid][i] , p1 ) ; 
			tmpW[tid][i] = p1.w ;
		}
	}

	// Morphing of position and weight
	for( int pid=0; pid<profiles2d.size(); ++pid ){

		int il = -1;
		for( int i=0; i<Tmps.size(); ++i )
			if( templatesIds[i] <= pid )  il = i ;

		if( il!=-1 && il!=Tmps.size()-1){
			
			int ir = il+1 ;
			int nstep = templatesIds[ir] - templatesIds[il] ;

			// position
			std::vector<Point2f> p_inc(cvnum) ;
			for( int i=0; i<cvnum; ++i )
				p_inc[i] = Point2f( tmpX[ir].t[i].X() - tmpX[il].t[i].X(), tmpX[ir].t[i].Y() - tmpX[il].t[i].Y()  ) / nstep ;

			int n = pid - templatesIds[il] ;
			for( int cvid=0; cvid<cvnum; ++cvid )
				x[pid].t[cvid] = tmpX[il].t[cvid] + p_inc[cvid] * n ;

			// weight
			std::vector<double> w_inc(cvnum) ;
			for(int i=0; i<cvnum; ++i){
				ON_4dPoint p1,p2;
				Tmps[il].GetCV( cvIdx[il][i], p1 ) ; 
				Tmps[ir].GetCV( cvIdx[ir][i], p2 ) ;  
				w_inc[i] = (p2.w-p1.w) / nstep ;
			}

			//int n = pid - templatesIds[il] ;
			for( int cvid=0; cvid<cvnum; ++cvid ){
				w_cv[pid][cvid] =  tmpW[il][cvid] + w_inc[cvid] * n ;
			}

		}else{
			if( il==-1){
				for( int cvid=0; cvid<cvnum; ++cvid ){
					x[pid].t[cvid] = Point2f(0,0) ;
					w_cv[pid][cvid] = tmpW[0][cvid] ;
				}
			}else{
				for( int cvid=0; cvid<cvnum; ++cvid ){
					x[pid].t[cvid] = tmpX.back().t[cvid] ;
					w_cv[pid][cvid] = tmpW.back()[cvid] ;
				}
			}
		}

	}

	// morphing of tips
	if( firstIsTip ){
		int nstep = templatesIds[0] - 0 ;
		std::vector<Point2f> p_inc(cvnum) ;
		std::vector<double> w_inc(cvnum) ;
		for(int i=0; i<cvnum; ++i){
			ON_4dPoint p1,p2;

			p1.x = 0; p1.y = 0 ; p1.z=0 ; p1.w = 1.0 ;
			Tmps[0].GetCV( i, p2 ) ;  p2.x/=p2.w; p2.y/=p2.w ;
			
			p_inc[i] = Point2f(p2.x-p1.x, p2.y-p1.y ) / nstep ;
			w_inc[i] = (p2.w-p1.w) / nstep ;
		}

		for( int i=0; i<nstep; ++i ){
			// morphing of position
			for( int j=0;j<cvnum; ++j )
				x[i].t[j] -= p_inc[j] * (nstep-i) ;

			// morphing of weight of CV
			for( int j=0; j<cvnum; ++j )
				w_cv[i][j] -= w_inc[j] * (nstep-i) ;
		}

	}

	if( lastIsTip ){
		int nstep = xn-1 - templatesIds.back() ;
		int offset = templatesIds.back()  ;

		std::vector<Point2f> p_inc(cvnum) ;
		std::vector<double> w_inc(cvnum) ;
		for(int i=0; i<cvnum; ++i){
			ON_4dPoint p1,p2;

			Tmps.back().GetCV( cvIdx.back()[i], p1 ) ;  p1.x/=p1.w; p1.y/=p1.w ;
			p2.x = 0; p2.y = 0 ; p2.z=0 ; p2.w = 1.0 ;

			p_inc[i] = Point2f(p2.x-p1.x, p2.y-p1.y ) / nstep ;
			w_inc[i] = (p2.w-p1.w) / nstep ;
		}

		for( int i=1; i<=nstep; ++i ){
			// morphing of position
			for( int j=0;j<cvnum; ++j )
				x[offset+i].t[j] += p_inc[j] * i ;

			// morphing of weight of CV
			for( int j=0; j<cvnum; ++j )
				w_cv[offset+i][j] += w_inc[j] * i ;
		}

	}


	// morphing of loop

	if( isLoop ){

		Curve2D cvs0 = ReconstructorUtility::getCvPts( Tmps[0]) ;
		Curve2D cvse = ReconstructorUtility::getCvPts( Tmps.back()) ;

		int steps = profiles2d.size() - ( templatesIds.back() - templatesIds[0]) ;
		int n = profiles2d.size() ;
		Curve2D delta ;
		for( int i=0; i<cvs0.size(); ++i )
			delta.push_back( (cvse[i] - cvs0[i])/ steps)  ;


		int pid0 = templatesIds[0] ;
		for( int i=1; i<steps; ++i ){

			int pid = (pid0 - i + n) % n ;

			for( int j=0; j<cvs0.size(); ++j )
				x[pid].t[j] = x[pid0].t[j] + delta[j] * i;

		}
	}
	//added by huajie
	//adjustCtrlPoints(x);
	initMofit(x);
	
	// we have only the first seed curve to transform , so there is  no need to morph weight
	//WeightOfCV = w_cv ;

	return ReconstructorUtility::convertNRTransform2Vector(x) ;
}

void loadParaConfig(){


	//std::string dir = appstate.plyfilename ;
	//dir.erase( dir.begin() + dir.find_last_of("\/"), dir.end()  );
	//std::string parafname = dir + "\\config.txt" ;

	//std::ifstream paraifs(parafname) ;
	//if( paraifs.good() ) paraifs >> completionSolverPara::w_data ;
	//if( paraifs.good() ) paraifs >> completionSolverPara::w_rigid ;
	//if( paraifs.good() ) paraifs >> completionSolverPara::w_smooth;
	//if( paraifs.good() ) paraifs >> completionSolverPara::w_tmpconstraint ;


	std::cout<<"-------------------------------------------------"<<std::endl;
	std::cout<< "loadParaConfig() v2 called!" 
		<<"\nw_data = " << completionSolverPara::w_data 
		<<"\nw_rigid = " << completionSolverPara::w_rigid 
		<<"\nw_smooth = " << completionSolverPara::w_smooth 
		<<"\nw_tmpconstraint = " << completionSolverPara::w_tmpconstraint  
		<<std::endl;
	std::cout<<"-------------------------------------------------"<<std::endl;

	//paraifs.close() ;
}


void completionSolver::detectCorrespondence () {


	int cvnum = ReconstructorPara::cvNum ;
	corres.clear() ;  cvIdx.clear() ;


	bool corres_cues_valid = true ;
	for( int i=0; i<Tmps.size(); ++i ){
		bool found = false ;
		for( int id=0; id<corres_cues.size(); ++id )
			if( corres_cues[id].x == templatesIds[i]  )
				found = true ;

		if( !found )
			corres_cues_valid = false ;
	}

	if( Tmps.size()!= corres_cues.size()) corres_cues_valid = false ;
	
	if( !corres_cues_valid){
		std::cout <<"corres_cues not valid. " <<std::endl;
		int offset = 0;
		for( int i=0; i<Tmps.size()-1; ++i ){
			offset += TScope[i].Y() -  TScope[i].X() + 1 ;
			corres.push_back( ReconstructorUtility::detectCorrespondenceForTmps( 
				ReconstructorUtility::deformNurbs( Tmps[i], NRTSF(cvnum), WeightOfCV[offset-1] ) ,
				ReconstructorUtility::deformNurbs( Tmps[i+1], NRTSF(cvnum), WeightOfCV[offset] ) ,
				sharpIds[i], sharpIds[i+1] ) );
		}

	}else{

		std::cout <<"corres_cues  valid. " <<std::endl;

		std::vector<int>  ProfIds, CVIds ;
		for( int i=0; i<corres_cues.size();++i){ ProfIds.push_back(corres_cues[i].x);  CVIds.push_back( corres_cues[i].y); }

		for( int i=0; i<Tmps.size()-1; ++i ){
			std::vector<int2> corr;

			int CVid_i = CVIds[ ReconstructorUtility::FindInVerctor(ProfIds, templatesIds[i] ) ] ;
			int CVid_i1= CVIds[ ReconstructorUtility::FindInVerctor(ProfIds, templatesIds[i+1] ) ] ;

			for( int j=0; j<cvnum; ++j)
				corr.push_back( int2( (CVid_i+j)%cvnum, (CVid_i1+j)%cvnum ) ) ;

			corres.push_back( corr );
		}

	}

	cvIdx.resize( Tmps.size() ) ;
	for( int i=0; i<cvnum; ++i )
		cvIdx[0].push_back( i ) ;
	for( int i=1; i<cvIdx.size(); ++i )
		cvIdx[i].resize(cvnum) ;

	for( int Tid=0; Tid<Tmps.size()-1; ++Tid ){
		for( int i=0; i<cvnum; ++i )
			for( int j=0; j<corres[Tid].size(); ++j )
				if( corres[Tid][j].x == cvIdx[Tid][i] )
					cvIdx[Tid+1][i] = corres[Tid][j].y ;
	}


	std::cout<<"correspondence:"<<std::endl;
	for( int i=0; i<corres.size(); ++i ){
		std::cout << i<<":"<<std::endl;
		for( int j=0; j<corres[i].size(); ++j)
			std::cout<<corres[i][j].x <<","<<corres[i][j].y <<std::endl;
	}



}
std::vector<double> completionSolver::solve( std::vector<double> x0_from_caller ){

	std::cout << "isLoop = " << isLoop <<std::endl;

	loadParaConfig() ;
	tmpConsWeights.resize(Tmps.size() ) ;
	for( int i=0;i<Tmps.size(); ++i)
		tmpConsWeights[i] = completionSolverPara::w_tmpconstraint ;


	double time_begin = clock() ;
	int cvnum = ReconstructorPara::cvNum ;
	int xn = profiles2d.size() ;

	xn *= ReconstructorPara::cvNum * 2 + 3  ;

	int nrtsf_size =  ReconstructorPara::cvNum * 2 + 3  ;

	std::vector<double> x0( xn, 0.0 ) ;

	wdataScaled = false ;
	

	// ------------------------- Warm-Up --------------------------------------------------------
	corres.clear() ;
	WhichStage = 1 ;



	// ------------------------- solve correspondences -----------------------------------------
	detectCorrespondence() ;


	// ------------------------- Morphing ------------------------------------------------------
	if( Warm0_Morph1_Fit2 >0 ){

		if( featureStroke_w2d.size() == 0 ){
			x0 = morphing() ;
		}
		else
			x0 = init_afterFeatureStroke() ;

		FinalNRTrans = ReconstructorUtility::convertVector2NRTransform( x0, ReconstructorPara::cvNum ) ;
	}


#ifdef After_solve_DEBUG
	FinalNRTrans = ReconstructorUtility::convertVector2NRTransform( x0, ReconstructorPara::cvNum ) ;
	goto solved;
#endif

	//if( sharpStrokeScope.size() )  // check initializing
	//	goto solved ;


	//// add new func here shihuajie
	//initMofit(FinalNRTrans);

	// ------------------------- Fit -------------------------------------------------------


	if( _Fitting_Only )
		x0 = std::vector<double>( xn, 0.0 ) ;

	if( Warm0_Morph1_Fit2 >1  ){

		WhichStage = 2 ;

		if( x0.size() == x0_from_caller.size() )
			x0 = x0_from_caller ;

		x0 = oneIterationSolver(x0) ;

		//for( int lftId=-1; lftId<(int)(Tmps.size()); ++lftId )
		//	x0 = solveOneSegment(x0, lftId) ;
		
		FinalNRTrans = ReconstructorUtility::convertVector2NRTransform( x0, ReconstructorPara::cvNum ) ;

	}


solved:

	// ------------------------- calculate finalProfiles and finalCtrPoints ------------------------------------------

	// ---------make tip be tip
	NRTSF trsf(cvnum) ;
	std::vector<Point2f> ctrpts = ReconstructorUtility::getCvPts( Tmps[0]) ;
	for( int i=0; i<cvnum; ++i)
		trsf.t[i] = -ctrpts[i] ;
	if( firstIsTip )
		FinalNRTrans[0] = trsf ;
	if( lastIsTip )
		FinalNRTrans.back() = trsf ;
	x0 = ReconstructorUtility::convertNRTransform2Vector(FinalNRTrans) ;

	// ---------make tip be tip  *end





	finalX = x0 ;


	finalProfiles.clear() ;
	finalProfiles.resize( profiles2d.size() ) ;

	finalCtrPoints.clear() ;
	finalCtrPoints.resize( profiles2d.size() ) ;

	finalNurbs.clear() ;
	finalNurbs.resize( profiles2d.size() ) ;


	for( int pid = 0; pid< profiles2d.size();  ++pid )
		finalNurbs[pid] = ReconstructorUtility::deformNurbs(Tmps[0], FinalNRTrans[pid], WeightOfCV[pid] ) ;

	// if no sharp stroke, then consolidate the nurbs
	if( featureStroke_w2d.size()==0 )
		finalNurbs = consolidateFinalNurbs( finalNurbs ) ;

	finalNurbs = smoothFinalTrajectory(finalNurbs) ;


	for( int pid = 0; pid< profiles2d.size();  ++pid ){
		finalProfiles[pid] = ReconstructorUtility::discretizeNurbs( finalNurbs[pid], ReconstructorPara::finalProfilesSampleNum ) ;
		finalCtrPoints[pid] = ReconstructorUtility::convert2dProfileTo3d( ReconstructorUtility::getCvPts( finalNurbs[pid] ), cutplanes[pid], centers[pid] );
	}


	if( _show_correspondence &&0){
		//for( int i=0; i<xn; ++i )
		//	gradient(x0, i) ;

		extern std::vector<std::pair<Point3f, Point3f> > globalVectorsToDraw  ; 
		for( int i=0; i<corres.size(); ++i ){


			std::vector<Point2f> CVs1, CVs2 ;
			std::vector<Point3f> CVs3d1, CVs3d2 ;

			ON_NurbsCurve Tl, Tr ;
			Tl = ReconstructorUtility::deformNurbs(Tmps[i], FinalNRTrans[TScope[i].Y()], WeightOfCV[TScope[i].Y()]  ) ;
			Tr = ReconstructorUtility::deformNurbs(Tmps[i+1], FinalNRTrans[TScope[i+1].X()], WeightOfCV[TScope[i+1].X()]  ) ;
			for( int j=0; j<corres[i].size(); ++j ){

				ON_4dPoint p1,p2 ;
				Tl.GetCV( corres[i][j].first, p1 ) ;
				Tr.GetCV( corres[i][j].second, p2 ) ;

				CVs1.push_back( Point2f(p1.x/p1.w, p1.y/p1.w));
				CVs2.push_back( Point2f(p2.x/p2.w, p2.y/p2.w));
			}

			CVs3d1 = ReconstructorUtility::convert2dProfileTo3d( CVs1, cutplanes[ TScope[i].Y()], centers[TScope[i].Y()] ) ;
			CVs3d2 = ReconstructorUtility::convert2dProfileTo3d( CVs2, cutplanes[ TScope[i+1].X()], centers[TScope[i+1].X()] ) ;

			for(int id=0; id<CVs3d1.size(); ++id )
				globalVectorsToDraw.push_back( std::pair<Point3f,Point3f>( CVs3d1[id], CVs3d2[id] )  ) ;
		}

		std::cout<<"globalVectorsToDraw:"<<std::endl;
		for( int i=0; i<globalVectorsToDraw.size(); ++i ){
			std::cout<<"("<<globalVectorsToDraw[i].first[0]<<","<<globalVectorsToDraw[i].first[1]<<","<<globalVectorsToDraw[i].first[2]<<"), " ;
			std::cout<<"("<<globalVectorsToDraw[i].second[0]<<","<<globalVectorsToDraw[i].second[1]<<","<<globalVectorsToDraw[i].second[2]<<")" <<std::endl;
		}

	}

	if( !_show_correspondence  &&1){
		extern std::vector<std::pair<Point3f, Point3f> > globalVectorsToDraw  ; 
		globalVectorsToDraw.clear() ;
	}

	//std::cout << "weight: " ;
	//for( int i=0; i<finalNurbs.size(); ++i )
	//	for( int j=0; j<finalNurbs[i].CVCount(); ++j ){
	//		ON_4dPoint p ;
	//		finalNurbs[i].GetCV(j, p); 
	//		std::cout << p.w <<" " ;

	//	}
	std::cout<<std::endl;


	double time_end = clock() ;

	std::cout<<"-------------------------------------------------------   TIME: "<< (time_end -time_begin)/(double)CLOCKS_PER_SEC  <<" sec."<<std::endl; 

	return x0 ;
}




std::vector<double> completionSolver::resolve( int newTmpId, ON_NurbsCurve newTmp ) {


	loadParaConfig() ;

	int cvnum = ReconstructorPara::cvNum ;

	//// add new template
	//bool  newTmpInsertFirst = false ;
	//bool  newTmpWriteFirst = false ;
	std::vector<int> templatesIds_bk = templatesIds ;                  
	std::vector<ON_NurbsCurve> Tmps_bk = Tmps ;
	//for( int i=0; i<templatesIds.size(); ++i ){

	//	if( newTmpId == templatesIds[i] ){
	//		Tmps[i] = newTmp ;
	//		if( i== 0 ) newTmpWriteFirst = true ;
	//		break;
	//	}

	//	if( newTmpId < templatesIds[i] ){
	//		templatesIds.insert( templatesIds.begin()+i, newTmpId ) ;
	//		Tmps.insert( Tmps.begin()+i, newTmp ) ;
	//		if( i== 0 ) newTmpInsertFirst = true ;
	//		break;
	//	}
	//}

	templatesIds.resize(1) ;  templatesIds[0] = newTmpId ;
	Tmps.resize(1) ; Tmps[0] = newTmp ;


	///////////////////////////  part of the construction function ////////////////////
	// set sharpIds
	sharpIds.clear() ;
	for( int i=0; i<Tmps.size(); ++i )
		sharpIds.push_back( ReconstructorUtility::getSharpCVofNurbs(Tmps[i]) ) ;


	// set Tscope
	TScope.resize( Tmps.size() ) ;
	for( int Tid=0; Tid<TScope.size(); ++Tid ){
		// compute sweeping scope of TEMPLATE Tid
		int Tnum =  templatesIds.size()   ;
		int Pid_begin, Pid_end ;

		if( Tid == 0 ) 	Pid_begin = 0 ;
		else Pid_begin = ( templatesIds[Tid-1] + templatesIds[Tid] ) / 2;

		if( Tid == Tnum-1) Pid_end = profiles2d.size()-1 ;
		//else Pid_end = ( templatesIds[Tid+1] + templatesIds[Tid] ) / 2;
		//else Pid_end = templatesIds[Tid+1];
		else Pid_end = ( templatesIds[Tid+1] + templatesIds[Tid] ) / 2 - 1;

		TScope[Tid].X() = Pid_begin ;
		TScope[Tid].Y() = Pid_end ;

	}




	//////////////////////////////////////////////////////////////////////////////

	// update tmpConsWeights
	tmpConsWeights.clear() ;
	tmpConsWeights.resize( templatesIds.size() ) ;
	for( int i=0; i<tmpConsWeights.size(); ++i )
		if( templatesIds[i]== newTmpId)
			tmpConsWeights[i] = 1e5;
		else
			tmpConsWeights[i] = 0;




	// calculate initial x
	std::vector<double> x0 ;

	std::vector<NRTSF> tsfm0 = FinalNRTrans ;
	std::vector<Point2f> cvsSrc = ReconstructorUtility::getCvPts(newTmp) ;
	std::vector<Point2f> cvsCtrlTmpOld =  ReconstructorUtility::getCvPts( finalNurbs[newTmpId] );
	for( int pid=0; pid<profiles2d.size(); ++pid ){
		std::vector<Point2f> cvsDst =  ReconstructorUtility::getCvPts( finalNurbs[pid] );
		for( int i=0; i<cvsSrc.size(); ++i ){
			tsfm0[pid].t[i] =  ( cvsDst[i] - cvsSrc[i]);
			tsfm0[pid].A = 0 ;
			tsfm0[pid].b = Point2f(0,0) ;
		}


		// -----  deform
		int winside = completionSolverPara::edit_winside/2 ;
		double sigma = 0.3 * ( winside / 2.0 - 1.0 )  + 0.8  ;
		int offset = abs(templatesIds[0]-pid) ; 
		double  weight =  pow(2.718, (-offset*offset/(2*sigma*sigma) ) ) ;
		for( int i=0; i<cvsSrc.size(); ++i ){
			tsfm0[pid].t[i] +=  ( cvsSrc[i] - cvsCtrlTmpOld[i])*weight;
		}


	}





	x0 = ReconstructorUtility::convertNRTransform2Vector(tsfm0) ;
	FinalNRTrans = tsfm0 ;


	if( featureStroke_l3d.size() ){
		x0 = init_afterFeatureStroke() ;
		FinalNRTrans = ReconstructorUtility::convertVector2NRTransform(x0,cvnum) ;
	}
	
	// resolve
	WhichStage = 2 ;
	//completionSolverPara::w_data *= 0.01 ;
	wdataScaled = true ;
	x0 = oneIterationSolver(x0) ;

	FinalNRTrans = ReconstructorUtility::convertVector2NRTransform( x0, ReconstructorPara::cvNum ) ;

	finalX = x0 ;


	finalProfiles.clear() ;
	finalProfiles.resize( profiles2d.size() ) ;

	finalCtrPoints.clear() ;
	finalCtrPoints.resize( profiles2d.size() ) ;

	finalNurbs.clear() ;
	finalNurbs.resize( profiles2d.size() ) ;


	for( int pid = 0; pid< profiles2d.size();  ++pid )
		finalNurbs[pid] = ReconstructorUtility::deformNurbs(Tmps[0], FinalNRTrans[pid], WeightOfCV[pid] ) ;

	// if no sharp stroke, then consolidate the nurbs
	if( featureStroke_w2d.size()==0 )
		finalNurbs = consolidateFinalNurbs( finalNurbs ) ;

	finalNurbs = smoothFinalTrajectory(finalNurbs) ;


	for( int pid = 0; pid< profiles2d.size();  ++pid ){
		finalProfiles[pid] = ReconstructorUtility::discretizeNurbs( finalNurbs[pid], ReconstructorPara::finalProfilesSampleNum ) ;
		finalCtrPoints[pid] = ReconstructorUtility::convert2dProfileTo3d( ReconstructorUtility::getCvPts( finalNurbs[pid] ), cutplanes[pid], centers[pid] );
	}





	return x0 ;
}

double completionSolver::xmult(Point2f p1,Point2f p2,Point2f p0)
{
	return (p1.X() -p0.X())*(p2.Y()-p0.Y()) - (p2.X()-p0.X())*(p1.Y()-p0.Y());
}

void completionSolver::MiniDisWith2pointss(Point2f p,Point2f q,int n)
{
	centerOfProfile = (p+q)/2.0;
	radiusOfProfile = (p-q).Norm()/2.0;
	int i;
	float c1,c2,t1,t2,t3;
	for(i = 0;i <= n;i ++)
	{
		if((pointForEval[i]-centerOfProfile).Norm() <= radiusOfProfile)
			continue;
		if(xmult(p,q,pointForEval[i]) != 0)
		{
			c1 = (p.X()*p.X()+p.Y()*p.Y()-q.X()*q.X()-q.Y()*q.Y())/2.0;
			c2 = (p.X()*p.X()+p.Y()*p.Y()-pointForEval[i].X()*pointForEval[i].X()-pointForEval[i].Y()*pointForEval[i].Y())/2.0;
			centerOfProfile.X() =(c1*(p.Y()-pointForEval[i].Y())-c2*(p.Y()-q.Y()))/((p.X()-q.X())*(p.Y()-pointForEval[i].Y())-(p.X()-pointForEval[i].X())*(p.Y()-q.Y()));
			centerOfProfile.Y() =(c1*(p.X()-pointForEval[i].X())-c2*(p.X()-q.X()))/((p.Y()-q.Y())*(p.X()-pointForEval[i].X())-(p.Y()-pointForEval[i].Y())*(p.X()-q.X()));
			radiusOfProfile=(centerOfProfile-pointForEval[i]).Norm();
		}
		else
		{
			t1 = (p-q).Norm();
			t2 = (p-pointForEval[i]).Norm();
			t3 = (q-pointForEval[i]).Norm();
			if(t1>=t2&&t1>=t3)
			{
				centerOfProfile=(p+q)/2.0;
				radiusOfProfile=(p-q).Norm()/2.0;
			}
			else if(t2>=t1&&t2>=t3)
			{
				centerOfProfile=(pointForEval[i]+q)/2.0;
				radiusOfProfile=(pointForEval[i]-q).Norm()/2.0;
			}
			else
			{
				centerOfProfile=(pointForEval[i]+p)/2.0;
				radiusOfProfile=(pointForEval[i]-p).Norm()/2.0;
			}
		}
	}
}

void completionSolver::MiniDisWithpointss(Point2f pi,int n)
{
	centerOfProfile = (pi+pointForEval[0])/2.0;
	radiusOfProfile = (pointForEval[0]-pi).Norm()/2.0;
	for(int j = 1;j <= n;j ++)
	{
		if((pointForEval[j]-centerOfProfile).Norm() <= radiusOfProfile)
			continue;
		MiniDisWith2pointss(pi,pointForEval[j],j-1);
	}
}


//************************************
// 函数名称: initMofit
// 函数说明：morph to fit the raw point
// 作    者：huajie.Shi
// 日    期：2014.8.2
// 返 回 值: void
// 参    数: std::vector<NRTSF> &x
//************************************
void completionSolver::initMofit(std::vector<NRTSF> &x){

	int cvnum = ReconstructorPara::cvNum ;
	int xn = profiles2d.size();
	std::vector<NRTSF> x0(xn, NRTSF(cvnum));
	for( int i=0; i<x.size(); i++ ){
		for(int j=0; j < cvnum; ++j){
			x0[i].t[j] = x[i].t[j];
		}
		x0[i].A = x[i].A;
		x0[i].b = x[i].b;
		//std::cout<<"A=="<<x[i].A<<"   b=="<<x[i].b.X()<<"  "<<x[i].b.Y()<<endl;
	}

	for( int i=0; i<x.size(); i++ ){
		for(int j=0; j < cvnum; ++j){
			ON_4dPoint p;
			Tmps[0].GetCV( j, p ) ;  p.x/=p.w; p.y/=p.w ;
			x[i].t[j].X() += p.x;
			x[i].t[j].Y() += p.y;
		}
	}

	//smooth the center of each ctrlpoint in curve
	Profile3D center;
	Profile3D output;
	//Profile2D centerOfx2D;
	center.clear();
	output.clear();
	//According focus Computing x Center
	for(int i=0;i<x.size();++i){	
		Point3f centerOfProfile3D(0.0,0.0,0.0);
		if(profiles2d[i].size()<5){
			center.push_back(centerOfProfile3D);
			continue;
		}
		Profile3D curve3D = ReconstructorUtility::convert2dProfileTo3d(profiles2d[i], cutplanes[i],  centers[i] ) ;
		for(int j=0;j<curve3D.size();j++){
			centerOfProfile3D = centerOfProfile3D + curve3D[j];
		}
		centerOfProfile3D = centerOfProfile3D / curve3D.size();
		center.push_back(centerOfProfile3D);
	}
	//ReconstructorUtility::getNurbs( center, output, center.size(), 3 );
	ReconstructorUtility::getBezier( center,  output, center.size()) ;

	//centerOfx2D.push_back(ReconstructorUtility::convert3dProfileTo2d( output,cutplanes[i],centers[i] ).at(i)) ;


	std::vector<Point2f> centerOfX; 
	std::vector<double> ratio;
	ratio.clear();
	centerOfX.clear();
	//calculate the translate
	for(int i=0;i<x.size();++i){
	//Filtered Correspondence
/*		bool flag=false;
		for( int j=0; j<Tmps.size(); ++j )
			if( templatesIds[j] == i ){
				flag=true;
				break;
			}
		if(flag) continue;*/

	//According focus Computing x Center
		/////////////////////
	/*	if(profiles2d[i].size()<=10) return ;
		pointForEval.clear() ;
		pointForEval.resize( profiles2d[i].size()) ;
		for(int k=0;k<profiles2d[i].size();k++){
			pointForEval[k]=profiles2d[i][k];
		}
		//random points 
		random_shuffle(pointForEval.begin(), pointForEval.end());
		centerOfProfile = (pointForEval[0]+pointForEval[1])/2.0;
		radiusOfProfile = (pointForEval[0]-pointForEval[1]).Norm()/2.0;
		for(int k = 2;k < pointForEval.size();k++)
		{
			if((centerOfProfile-pointForEval[k]).Norm() <= radiusOfProfile)
				continue;
			MiniDisWithpointss(pointForEval[k],k-1);
		}
		std::cout<<"radiusOfProfile"<<radiusOfProfile<<endl;
		*/////////////////////
		std::vector<Point2f> curve = ReconstructorUtility::discretizeNurbs( ReconstructorUtility::deformNurbs( Tmps[0], x0[i],WeightOfCV[i] ),  100 );
		centerOfx=Point2f(0,0);
		for(int j=0;j<curve.size();j++){
			centerOfx += curve[j];
		}
		centerOfx /= curve.size();
		
		//center.push_back(centerOfx);
	/*	centerOfx=Point2f(0,0);
		for(int j=0;j<x[i].t.size();j++){
			centerOfx += x[i].t[j];
		}
		centerOfx /=x[i].t.size();
		*/

		//centerOfProfile=ReconstructorUtility::convert3dProfileTo2d( output,cutplanes[i],centers[i] ).at(i);
	//According focus Computing profile Center
	/*	for(int j=0;j<profiles2d[i].size();j++){
			centerOfProfile += profiles2d[i][j];
		}
		if(profiles2d[i].size()!=0)
			centerOfProfile /=profiles2d[i].size();
		else centerOfProfile=Point2f(0,0);
	*/
		

	//Calculated according to the minimum point of the center surrounded by round
	/*	pointForEval.clear() ;
		pointForEval.resize( x[i].t.size()) ;
		for(int j=0;j<x[i].t.size();j++){
			pointForEval[j] = x[i].t[j];
		}
		//random points 
		random_shuffle(pointForEval.begin(), pointForEval.end());
		centerOfProfile = (pointForEval[0]+pointForEval[1])/2.0;
		radiusOfProfile = (pointForEval[0]-pointForEval[1]).Norm()/2.0;
		for(int j = 2;j < pointForEval.size();j++)
		{
			if((centerOfProfile-pointForEval[j]).Norm() <= radiusOfProfile)
				continue;
			MiniDisWithpointss(pointForEval[j],j-1);
		}
	*/
		centerOfProfile = Point2f(0.0,0.0);
		//centerOfProfile =ReconstructorUtility::convert3dProfileTo2d( output,cutplanes[i],centers[i] ).at(i);
		//std::cout<<"size=="<<ReconstructorUtility::convert3dProfileTo2d( output,cutplanes[i],centers[i] ).size()<<endl;;
		//Point2f trans = centerOfProfile - centerOfx;
		//for(int j=0;j<x[i].t.size();j++){
		//	x[i].t[j] += trans;
		//}

		//centerOfx = centerOfProfile;
		centerOfX.push_back(centerOfx);

		double cvLength=0.0;
		double profileLength=0.0;
		//Filtered Correspondence
/*		bool flag=false;
		for( int j=0; j<Tmps.size(); ++j )
			if( templatesIds[j] == i ){
				flag=true;
				break;
			}
		if(flag){
			ratio.push_back( 1 );
			continue;
		}
*/
		for( int j=0; j < x[i].t.size(); ++j ){
			cvLength += ( x[i].t[j]-centerOfx).Norm();
		}
		cvLength /=  x[i].t.size();
		if(cvLength <= 0.00001){
			ratio.push_back( 1 );
			continue;
		}
		for(int k=0;k<profiles2d[i].size();++k){
			profileLength += (profiles2d[i][k]-centerOfProfile).Norm();
		}
		if(profiles2d[i].size()==0){
			ratio.push_back( 1 );
			continue;
		}
		profileLength /= profiles2d[i].size();
		//profileLength = radiusOfProfile;
		if( (profileLength/cvLength) < 3  ) 
			ratio.push_back( profileLength / cvLength );
		else
			ratio.push_back( 3 );
	}

	//std::cout<<"ratio.size=="<<ratio.size()<<endl;
	//std::cout<<"x.size=="<<x.size()<<endl;

	for( int i=0; i<x.size(); i++ ){	
		for(int j=0; j< x[i].t.size(); ++j){
			x[i].t[j] -= centerOfX[i];
		}
	}
	
	ReconstructorUtility::smoothVectorData( ratio,30) ;
	for( int i=0; i<x.size(); i++ ){
		for( int j=0; j < x[i].t.size(); ++j ){
			x[i].t[j] *= ratio[i];
		}
	}

	for( int i=0; i<x.size(); i++ ){	
		for(int j=0; j< x[i].t.size(); ++j){
			x[i].t[j] += centerOfX[i];
		}
	}
	
	for( int i=0; i<x.size(); i++ ){
		for(int j=0; j<cvnum; ++j){
			ON_4dPoint p;
			Tmps[0].GetCV( j, p ) ;  p.x/=p.w; p.y/=p.w ;
			x[i].t[j].X() -= p.x;
			x[i].t[j].Y() -= p.y;
		}
	}
}
//************************************
// 函数名称: adjustCtrlPoints
// 函数说明：adjust the control point
// 作    者：huajie.Shi
// 日    期：2014.8.6
// 返 回 值: void
// 参    数: std::vector<NRTSF> &x
//************************************
void completionSolver::adjustCtrlPoints(std::vector<NRTSF> &x){
	int cvnum = ReconstructorPara::cvNum ;
	std::cout<<"cvnum=="<<cvnum<<endl;
	int xn = profiles2d.size();
	std::vector<NRTSF> x0(xn, NRTSF(cvnum));
	for( int i=0; i<x.size(); i++ ){
		for(int j=0; j < cvnum; ++j){
			ON_4dPoint p;
			Tmps[0].GetCV( j, p ) ;  p.x/=p.w; p.y/=p.w ;
			x[i].t[j].X() += p.x;
			x[i].t[j].Y() += p.y;
		}
	}
	//////////////////////////////////
	Point2f d1,d2;
	for( int i=0; i<x.size(); ++i ){
		//int i = templatesIds[k];
		int s=x[i].t.size();
		for(int j=0;j<s;j++){
			d1 = x[i].t[   j   ] - x[i].t[(j-1+s)%s];
			d2 = x[i].t[(j+1)%s] - x[i].t[(j-1+s)%s];
			if( (d1.X() * d2.Y() - d1.Y() * d2.X()) > 0 ){
				x[i].t[j] = (x[i].t[(j-1+s)%s] + x[i].t[(j+1)%s]) / 2;
			}
		}
			
		for(int j=0; j < cvnum; ++j){
			x0[i].t[j] = x[i].t[j];
			ON_4dPoint p;
			Tmps[0].GetCV( j, p ) ;  p.x/=p.w; p.y/=p.w ;
			x0[i].t[j].X() -= p.x;
			x0[i].t[j].Y() -= p.y;
		}
		x0[i].A = x[i].A;
		x0[i].b = x[i].b;

		std::vector<Point2f> curve = ReconstructorUtility::discretizeNurbs( ReconstructorUtility::deformNurbs( Tmps[0], x0[i],WeightOfCV[i] ),  100 );
		Profile2D profile;
		int size = curve.size() / cvnum;
		for(int j=0; j< cvnum; ++j){
			profile.push_back(curve[j*size]);		
		}
		x[i].t = ReconstructorUtility::discretizeNurbs( ReconstructorUtility::fitNURBSforProfile2d( profile ), cvnum);
/*		std::vector<double> alpha;
		alpha.clear();
		for(int j=0; j<curve.size(); ++j){
			Point2f vector1 = curve[(j+curve.size()) % curve.size()] - curve[(j-1+curve.size()) % curve.size()];
			Point2f vector2 = curve[(j+1+curve.size()) % curve.size()] - curve[(j+curve.size()) % curve.size()];
			alpha.push_back(ReconstructorUtility::vectorIntersectionAngle( vector1, vector2 ));
		}
		for(int j=0; j<curve.size()/cvnum; ++j){
			x[i].t[j] = curve[max_element(curve.begin() + j*curve.size(), curve.begin() + std::min((j+1)*curve.size()/cvnum,curve.size())) - curve.begin()];
		}
*/
	}
	
	//////////////////////////////////

	for( int i=0; i<x.size(); i++ ){
		for(int j=0; j < cvnum; ++j){
			ON_4dPoint p;
			Tmps[0].GetCV( j, p ) ;  p.x/=p.w; p.y/=p.w ;
			x[i].t[j].X() -= p.x;
			x[i].t[j].Y() -= p.y;
		}
	}
}