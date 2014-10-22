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

double completionSolver::E_data(std::vector<NRTSF> &x, int profId ) {



	double disSum = 0.0 ;

	ON_NurbsCurve Tmp = Tmps[0] ;
	NRTSF ngTrans = x[profId] ;
	ON_NurbsCurve Ttid = ReconstructorUtility::deformNurbs(Tmp, ngTrans,  WeightOfCV[profId] ) ;
	std::vector<Point2f> disnbs = ReconstructorUtility::discretizeNurbs(Ttid,completionSolverPara::pwCurveSampleNum, true )  ;  // 3/7

	if( _use_360dis )
		disSum += ReconstructorUtility::profileDis_360deg( profiles2d[profId ],  disnbs ) ; //
	else
		disSum += ReconstructorUtility::profileDis( profiles2d[profId ],  disnbs ) ;											//  3/7


	// scaling for edit the sweeping result interactively
	//if( wdataScaled ){
	//	int winside = completionSolverPara::edit_winside ;
	//	double sigma = 0.3 * ( winside / 2.0 - 1.0 )  + 0.8  ;
	//	int offset = abs(templatesIds[0]-profId) ; 
	//	double  scale =  1 - pow(2.718, (-offset*offset/(2*sigma*sigma) ) ) ;
	//	return disSum * scale ;
	//}

	return disSum ;
}


double completionSolver::E_rigid( std::vector<NRTSF> &x, int profId  )  {


	if( profId == x.size()-1 )
		return 0 ;

	for( int i=0; i<templatesIds.size(); ++i )
		if( profId == templatesIds[i] - 1 )
			return 0 ;

	double penSum = 0.0 ;

	NRTSF ngTrans0 = x[profId] ;
	NRTSF ngTrans1 = x[profId+1] ;

	for( int i=0; i<ngTrans0.t.size(); ++i )
		penSum += ( ngTrans0.t[i] - ngTrans1.t[i] ).SquaredNorm() ; 



	return penSum ;

}

double completionSolver::E_smooth( std::vector<NRTSF> &x, int profId, int changedCvID   ) {



	// modified for loop


	if( !oneIterationSolverCalled ){
		for( int i=0; i<templatesIds.size(); ++i )
			if( profId == templatesIds[i] - 1 )
				return 0 ;
	}


	int cvnum = x[0].t.size() ;
	int profnum = x.size() ;



	int winside = completionSolverPara::smooth_winside ;

	winside = 3 ;
	std::vector<std::vector<Point3f>> ctrlPointss(winside) ;

	for( int id=0; id<winside; ++id ){

		int offset = id-winside/2 ;

		int pid = profId+offset ;

		if( pid <0 || pid >= profiles2d.size() ){
			pid = (pid + profnum)%profnum ;

		}

		ON_NurbsCurve nbs = ReconstructorUtility::deformNurbs( Tmps[0], x[pid], WeightOfCV[pid] ) ;

		if( changedCvID < 0 ){
			std::vector<Point2f> CVs2d (cvnum);
			for( int cvi=0; cvi<cvnum; ++cvi ){
				ON_4dPoint p;  nbs.GetCV( cvi, p) ; p.x/=p.w; p.y/=p.w ;
				CVs2d[ cvi ] =  Point2f(p.x, p.y)  ;  
			}
			ctrlPointss[id] = ReconstructorUtility::convert2dProfileTo3d( CVs2d, cutplanes[ pid ], centers[ pid] ) ;
		}else{

			std::vector<Point2f> CVs2d (1);
			int cvi = changedCvID ;
			ON_4dPoint p;  nbs.GetCV( cvi, p) ; p.x/=p.w; p.y/=p.w ;
			CVs2d[ 0 ] =  Point2f(p.x, p.y)  ;  

			ctrlPointss[id] = ReconstructorUtility::convert2dProfileTo3d( CVs2d, cutplanes[ pid ], centers[ pid ] ) ;
		}

	}

	double penSum = 0.0 ;

	double totalWeight = 0.0 ;

	double sigma = 0.3 * ( winside / 2.0 - 1.0 )  + 0.8  ;

	for( int offset=1;  offset<=winside/2; ++offset ){


		if( !isLoop  ){
			if( profId - offset <0 || profId + offset >= profiles2d.size() )
				continue ;
		}


		int mid = winside/2 ;

		//double weight = pow(2.718, (-(offset-1)*(offset-1)/(2*sigma*sigma) ) ) ;
		double weight = 1.0 ;

		for( int i=0; i<ctrlPointss[0].size(); ++i )
			penSum +=  (1 - cos(ReconstructorUtility::vectorIntersectionAngle( ctrlPointss[mid+offset][i]-ctrlPointss[mid][i], ctrlPointss[mid][i]-ctrlPointss[mid-offset][i] ) ) ) * weight  ;

		totalWeight += weight ;

	}

	//if( profId == 1 || profId == profiles2d.size()-2 )
	//	penSum*=2 ;

	if( totalWeight == 0 )
		return 0 ;


	return penSum /totalWeight;
}

double completionSolver::E_template_constrain( std::vector<NRTSF> &x, int profId ) {


	int cvnum = x[0].t.size() ;

	double penSum = 0; 

	int tIdx = 0 ;
	for( int i=0; i<templatesIds.size(); ++i  ){

		int tid = templatesIds[i] ;

		if( tid!=profId)
			continue ;

		tIdx = i ;

		ON_NurbsCurve nbs1 = ReconstructorUtility::deformNurbs(Tmps[0], x[tid], WeightOfCV[tid] ) ;
		ON_NurbsCurve nbs2 = Tmps[i] ;
		double avedis = 0;
		for( int cvi=0; cvi<cvnum; ++cvi ){
			ON_4dPoint p;
			nbs1.GetCV( cvi, p) ; p.x/=p.w; p.y/=p.w ;
			Point2f cv1 =  Point2f(p.x, p.y)  ;   
			nbs2.GetCV( cvIdx[i][cvi], p) ; p.x/=p.w; p.y/=p.w ;
			Point2f cv2 =  Point2f(p.x, p.y)  ;   
			avedis += (cv1-cv2).SquaredNorm() ;
		}
		penSum += avedis/cvnum ;

	}

	return penSum  * tmpConsWeights[tIdx];
}

double completionSolver::E_tip_constrain( std::vector<NRTSF> &x, int profId ) {

	int cvnum = x[0].t.size() ;
	double penSum = 0; 
	// constraint of tips
	if( firstIsTip &&  profId == 0){
		std::vector<Point2f> cvs ;
		ON_NurbsCurve nbs1 = ReconstructorUtility::deformNurbs(Tmps[0], x[0], WeightOfCV[0] ) ;		
		double avedis = 0;
		for( int cvi=0; cvi<cvnum; ++cvi ){
			ON_4dPoint p;
			nbs1.GetCV( cvi, p) ; p.x/=p.w; p.y/=p.w ;
			Point2f cv1 =  Point2f(p.x, p.y)  ;   
			cvs.push_back( cv1 ) ;
		}

		//Point2f center = GlobalFun::centerOfPoints( cvs ) ; 
		for( int i=0; i<cvs.size(); ++i )
			//avedis += (cvs[i]-center).SquaredNorm() ;
			avedis += (cvs[i]).SquaredNorm() ;
		penSum += avedis/cvnum ;
	}
	if( lastIsTip &&  profId == profiles2d.size()-1){
		std::vector<Point2f> cvs ;
		ON_NurbsCurve nbs1 = ReconstructorUtility::deformNurbs(Tmps[0], x.back(), WeightOfCV.back() ) ;		
		double avedis = 0;
		for( int cvi=0; cvi<cvnum; ++cvi ){
			ON_4dPoint p;
			nbs1.GetCV( cvi, p) ; p.x/=p.w; p.y/=p.w ;
			Point2f cv1 =  Point2f(p.x, p.y)  ;   
			cvs.push_back( cv1 ) ;
		}

		//Point2f center = GlobalFun::centerOfPoints( cvs ) ;
		for( int i=0; i<cvs.size(); ++i )
			//avedis += (cvs[i]-center).SquaredNorm() ;
			avedis += (cvs[i]).SquaredNorm() ;
		penSum += avedis/cvnum ;
	}


	return penSum ;
}
double completionSolver::objFunc(  std::vector<double> x0_  ) {


	int pNum = profiles2d.size() ;
	int cvNum = ReconstructorPara::cvNum ;
	int nrtsf_size = ( ReconstructorPara::cvNum * 2 + 3 ) ;
	int nT = Tmps.size() ;

	int pid0,pid1,startX;
	if( oneIterationSolverCalled ){
		pid0 = 0 ; 
		pid1 = profiles2d.size() ;
		startX = pid0 * nrtsf_size ;
	}else{
		pid0 = LeftProfId<0?0: templatesIds[LeftProfId];
		pid1 = LeftProfId+1<nT? templatesIds[LeftProfId+1] :  profiles2d.size();
		startX = pid0 * nrtsf_size ;
	}


	// get x
	std::vector<double> x0 = lastX ;
	for( int i=0; i<x0_.size(); ++i )
		x0[startX+i] = x0_[i];

	std::vector<NRTSF> x = ReconstructorUtility::convertVector2NRTransform(x0, cvNum) ;



	double edata = 0  ;
	double erigid = 0  ;
	double esmooth = 0  ;
	double e_tmpconstrain = 0  ;
	double e_tipconstrain = 0  ;
	double e_feaconstrain = 0  ;



#ifdef _use_openmp_
#define TNUM _Solver_TNUM_

	double edataArray[TNUM] ;
	double erigidArray[TNUM] ;
	double esmoothArray[TNUM] ;
	double etmpconstArray[TNUM] ;
	double etipconstArray[TNUM] ;
	double efeaconstArray[TNUM] ;

	for( int i=0; i<TNUM;++i)
		edataArray[i]=erigidArray[i]=esmoothArray[i]=etmpconstArray[i]=etipconstArray[i]=efeaconstArray[i]=0;

#pragma omp parallel num_threads(TNUM)
	{
		int tid = omp_get_thread_num();
		for( int j=0; j<(pid1-pid0)/TNUM+1 ; ++j ){
			int pid = pid0 + j*TNUM+tid ;
			if( pid < pid1 ){
				edataArray[tid]   += E_data(x,pid ) * completionSolverPara::w_data ;
				erigidArray[tid]  += E_rigid(x,pid )* completionSolverPara::w_rigid ;
				esmoothArray[tid] += E_smooth(x,pid )* completionSolverPara::w_smooth  ;
				etmpconstArray[tid] += E_template_constrain(x, pid)* completionSolverPara::one  ;
				etipconstArray[tid] += E_tip_constrain(x, pid)* completionSolverPara::w_tipconstraint  ;
				efeaconstArray[tid] += E_featrueStroke(x, pid)* completionSolverPara::w_featureconstraint  ;
			}
		}


	}

	for( int i=0; i<TNUM;++i){
		edata += edataArray[i];
		erigid += erigidArray[i];
		esmooth += esmoothArray[i];
		e_tmpconstrain += etmpconstArray[i];
		e_tipconstrain += etipconstArray[i];
		e_feaconstrain += efeaconstArray[i];
	}
#else
	for( int i=pid0; i<pid1; ++i ){
		edata += E_data(x,i ) * completionSolverPara::w_data ;
		erigid += E_rigid(x,i )* completionSolverPara::w_rigid ;
		esmooth += E_smooth(x,i )* completionSolverPara::w_smooth  ;
		e_tmpconstrain += E_template_constrain(x, i)* completionSolverPara::one  ;
		e_tipconstrain += E_tip_constrain(x, i)* completionSolverPara::w_tipconstraint  ;
		e_feaconstrain += E_featrueStroke(x, i)* completionSolverPara::w_featureconstraint  ;
	}
#endif


	double energy = edata 
		+ erigid 
		+ esmooth
		+ e_tmpconstrain
		+ e_tipconstrain 
		+ e_feaconstrain ;

	if( completionSolverPara::diagnosticMessage) {
		std::cout <<edata<<" + "<<  erigid <<" + " <<  esmooth  <<" + " << e_tmpconstrain<<" + " << e_tipconstrain<<" + " << e_feaconstrain << " = " << energy <<" | " <<std::flush;
	}

	return energy;
}


double completionSolver::gradient( const std::vector<double> &x_ , int xdiff_id, double oriProfEnergy ) {
	
	int pNum = profiles2d.size() ;
	int cvNum = ReconstructorPara::cvNum ;
	int nrtsf_size = ( ReconstructorPara::cvNum * 2 + 3 ) ;
	int nT = Tmps.size() ;

	int pid0,pid1,startX;
	if( oneIterationSolverCalled ){
		pid0 = 0 ; 
		pid1 = profiles2d.size() ;
		startX = pid0 * nrtsf_size ;
	}else{
		pid0 = LeftProfId<0?0: templatesIds[LeftProfId];
		pid1 = LeftProfId+1<nT? templatesIds[LeftProfId+1] :  profiles2d.size();
		startX = pid0 * nrtsf_size ;
	}


	// get x
	std::vector<double> x = lastX ;
	for( int i=0; i<x_.size(); ++i )
		x[startX+i] = x_[i];

	xdiff_id += startX ;
	
	



	std::vector<NRTSF> nrtsf = ReconstructorUtility::convertVector2NRTransform(x, cvNum) ;

	int nrtsfId = xdiff_id / nrtsf_size ;
	int innerId = xdiff_id % nrtsf_size ;
	
	// get step
	double step ;
	if( innerId < cvNum * 2 )
		step = completionSolverPara::tStep ;
	else if( innerId == cvNum * 2 )
		step = completionSolverPara::AStep ;
	else
		step = completionSolverPara::bStep ;



	std::vector<double> xr = x ;  xr[xdiff_id] += step ;
	std::vector<double> xl = x ;  xl[xdiff_id] -= step ;
/*
	std::vector<NRTSF> nrtsf_r = ReconstructorUtility::convertVector2NRTransform( xr, cvNum ) ;
	std::vector<NRTSF> nrtsf_l = ReconstructorUtility::convertVector2NRTransform( xl, cvNum ) ;*/


	std::vector<NRTSF> nrtsf_r = nrtsf ;
	std::vector<NRTSF> nrtsf_l = nrtsf ;

	std::vector<double> item_r,item_l;
	item_r.insert(item_r.begin(), xr.begin()+nrtsfId*nrtsf_size, xr.begin()+nrtsfId*nrtsf_size+nrtsf_size ) ;
	item_l.insert(item_l.begin(), xl.begin()+nrtsfId*nrtsf_size, xl.begin()+nrtsfId*nrtsf_size+nrtsf_size ) ;
	nrtsf_r[nrtsfId] = ReconstructorUtility::convertVector2NRTransform( item_r, cvNum )[0] ;
	nrtsf_l[nrtsfId] = ReconstructorUtility::convertVector2NRTransform( item_l, cvNum )[0];


	int changedId = -1 ;
	if( innerId < cvNum * 2 )
		changedId = innerId/2 ;

	double edata = E_data( nrtsf_r, nrtsfId ) ;
	double erigid = E_rigid( nrtsf_r, nrtsfId ) ;
	double esmooth = E_smooth( nrtsf_r, nrtsfId,changedId ) ;
	double etmpconstrain = E_template_constrain( nrtsf_r, nrtsfId ) ;
	double etipconstrain = E_tip_constrain( nrtsf_r, nrtsfId ) ;
	double efeaconstrain = E_featrueStroke( nrtsf_r, nrtsfId ) ;

	double edata1 = E_data( nrtsf_l, nrtsfId ) ;
	double erigid1 = E_rigid( nrtsf_l, nrtsfId ) ;
	double esmooth1 = E_smooth( nrtsf_l, nrtsfId,changedId ) ;
	double etmpconstrain1 = E_template_constrain( nrtsf_l, nrtsfId ) ;
	double etipconstrain1 = E_tip_constrain( nrtsf_l, nrtsfId ) ;
	double efeaconstrain1 = E_featrueStroke( nrtsf_l, nrtsfId ) ;

	double energy = edata * completionSolverPara::w_data  
				+ erigid * completionSolverPara::w_rigid  
				+ esmooth * completionSolverPara::w_smooth  
				+ etmpconstrain * completionSolverPara::one
				+ etipconstrain * completionSolverPara::w_tipconstraint
				+ efeaconstrain * completionSolverPara::w_featureconstraint ;



	double energy1 = edata1 * completionSolverPara::w_data  
		+ erigid1 * completionSolverPara::w_rigid  
		+ esmooth1 * completionSolverPara::w_smooth  
		+ etmpconstrain1 * completionSolverPara::one
		+ etipconstrain1 * completionSolverPara::w_tipconstraint
		+ efeaconstrain1 * completionSolverPara::w_featureconstraint ;

	
	// compute neighbour smooth terms
	double delta_smooth_nei = 0 ;
	if( nrtsfId != 0)
		delta_smooth_nei += E_smooth(nrtsf_r, nrtsfId-1, changedId)  - E_smooth(nrtsf_l, nrtsfId-1, changedId)  ;
	if( nrtsfId != profiles2d.size()-1 )
		delta_smooth_nei += E_smooth(nrtsf_r, nrtsfId+1, changedId)  - E_smooth(nrtsf_l, nrtsfId+1,changedId )  ;
	delta_smooth_nei *= completionSolverPara::w_smooth  ;


	// compute neighbour rigid terms
	double delta_rigid_nei = 0 ;
	if( nrtsfId != 0)
		delta_rigid_nei += E_rigid(nrtsf_r, nrtsfId-1)  - E_rigid(nrtsf_l, nrtsfId-1)  ;
	delta_rigid_nei *= completionSolverPara::w_rigid  ;



	return (energy-energy1+delta_smooth_nei+delta_rigid_nei ) / (step*2) ;

}

