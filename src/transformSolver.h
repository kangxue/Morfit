

#ifndef _transform_solver_
#define  _transform_solver_


#include "smtf.h"
#include "skeleton_kx.h"


class transformSolver{
public:

	transformSolver( skeleton_one skel, int pro_id ){ 
		profiles2d = skel.profiles2d[pro_id]; sharpness = skel.sharpness[pro_id] ; 
		cutplanes = skel.cutplanes[pro_id];   profile_conf = skel.profile_conf[pro_id] ; 
	} 

	void solve(  double LSW = 1.0, int templateId = -1) ;


	std::vector<ST> FinalSimTrans ;                                    // similarity transform
	//std::vector<std::vector<Point2f>> finalProfiles ;				   // final profiles

	int currentTemplateId ;

private:


	std::vector<Profile2D> profiles2d ;   // 2d profiles
	std::vector<std::vector<double>> sharpness ;     // sharpness of each points

	std::vector<cutPlane> cutplanes ;                // profiling planes of each profile
	std::vector<double> profile_conf ;               // confidence of each profile


	//std::vector<Point2f> transformProfile( int proId,  ST stf) ;

	double profDis( std::vector<Point2f> &prof1, std::vector<Point2f> &prof2 , int id1, int id2  ) ;

	double objFunc(  std::vector<ST> &x , double LSWeight) ;

	std::vector<ST> oneIterationSolver(  const std::vector<ST> &x0,  double LSWeight ) ;


	double gradient( std::vector<ST> &stf, int xdiff_id,  double LSWeight ) ;

	double currentLSWeight ;

public:

	double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) ;

};


#endif