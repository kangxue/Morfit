
#ifndef _bspline_h_
#define  _bspline_h_ 

#include <vector>

template<class  _Point>
void getUniformCubicBSpline( std::vector<_Point> input, std::vector<_Point> &output, int pnumperbranch = 20 ){
	//  refer http://en.wikipedia.org/wiki/B-spline

	output.clear() ;


	double M[4][4] = {
		-1, 3, -3, 1,
		3,-6,3,0,
		-3,0,3,0,
		1,4,1,0
	} ;


	for( int i=0; i<input.size(); ++i  ){
		for( double t = 0; t<1.0; t+= 0.01 ){

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

			
			float zero[3] = {0.0, 0.0, 0.0 } ;
			_Point St( zero );
			//St = St - St ;

			double totoWeight = 0.0 ;
			for( int id = 0 ; id< 4; ++id ){
				if( i-1 + id >= 0 &&i-1 + id <input.size() ){

					//if( i-1 + id ==0 ||  i-1 + id == input.size()-1 ){
					//	St = St +  input[i-1 + id ] * TM[id]*2 ;
					//	totoWeight += TM[id] * 2;


					//}else{
						St = St +  input[i-1 + id ] * TM[id] ;
						totoWeight += TM[id] ;
					//}
				}
			}

			//std::cout<<"totoWeight = " <<totoWeight<<std::endl;
			St = St / totoWeight ;

			output.push_back(St) ;

		}

	}

	// calculate  branch length
	double stepLength = 0 ;
	for( int i=0; i< input.size()-1; ++i )
		stepLength += (input[i] - input[i+1]).Norm();
	stepLength /= pnumperbranch ;

	// reduce
	_Point prepoint = output[0] ;
	for( int j=1; j<output.size(); ++j ){

		if( (output[j] - prepoint).Norm() > stepLength  || j==output.size()-1 ){

			prepoint = output[j] ;
		}
		else{
			output.erase( output.begin() + j ) ;
			j-- ;
		}

	}


	// correct head and tail

	output[0] = input[0] ;
	output.back() = input.back() ;
}

#endif