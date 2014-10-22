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


#include <flann/flann.hpp>

#include "Types.h"
#include <vector>

#include "reconstructionUtility.h"

namespace ReconstructorUtility{

	void computeFlannKNN( const std::vector<Point2f> &dataPoints, const std::vector<Point2f> &queryPoints, std::vector<std::vector<int>> &nearestIds, int nnk ) {


		int nn = nnk;
		flann::Matrix<float> dataset(new float[dataPoints.size()*2],  dataPoints.size(), 2);
		flann::Matrix<float> query(new float[queryPoints.size()*2],  queryPoints.size(), 2);

		for( int i=0; i<dataPoints.size(); ++i ){
			dataset[i][0] = dataPoints[i][0] ;
			dataset[i][1] = dataPoints[i][1] ;
		}

		for( int i=0; i<queryPoints.size(); ++i ){
			query[i][0] = queryPoints[i][0] ;
			query[i][1] = queryPoints[i][1] ;
		}


		flann::Matrix<int> indices(new int[query.rows*nn], query.rows, nn);
		flann::Matrix<float> dists(new float[query.rows*nn], query.rows, nn);
		// construct an randomized kd-tree index using 4 kd-trees
		flann::Index<flann::L2<float> > index(dataset, flann::KDTreeIndexParams(1));
		index.buildIndex();
		// do a knn search, using 128 checks
		index.knnSearch(query, indices, dists, nn, flann::SearchParams(128));

		//flann::save_to_file(indices,"result.hdf5","result");


		nearestIds.clear() ;
		nearestIds.resize( queryPoints.size() ) ;
		for( int i=0; i<queryPoints.size(); ++i ){
			for( int j=0; j<nn; ++j )
				nearestIds[i].push_back( indices[i][j] ) ;
		}

		delete[] dataset.ptr();
		delete[] query.ptr();
		delete[] indices.ptr();
		delete[] dists.ptr();

	}



	void computeFlannKNN( const std::vector<Point3f> &dataPoints, const std::vector<Point3f> &queryPoints, std::vector<std::vector<int>> &nearestIds, int nnk ) {


		int nn = nnk;
		flann::Matrix<float> dataset(new float[dataPoints.size()*3],  dataPoints.size(), 3);
		flann::Matrix<float> query(new float[queryPoints.size()*3],  queryPoints.size(), 3);

		for( int i=0; i<dataPoints.size(); ++i ){
			dataset[i][0] = dataPoints[i][0] ;
			dataset[i][1] = dataPoints[i][1] ;
			dataset[i][2] = dataPoints[i][2] ;
		}

		for( int i=0; i<queryPoints.size(); ++i ){
			query[i][0] = queryPoints[i][0] ;
			query[i][1] = queryPoints[i][1] ;
			query[i][2] = queryPoints[i][2] ;
		}


		flann::Matrix<int> indices(new int[query.rows*nn], query.rows, nn);
		flann::Matrix<float> dists(new float[query.rows*nn], query.rows, nn);
		// construct an randomized kd-tree index using 4 kd-trees
		flann::Index<flann::L2<float> > index(dataset, flann::KDTreeIndexParams(1));
		index.buildIndex();
		// do a knn search, using 128 checks
		index.knnSearch(query, indices, dists, nn, flann::SearchParams(128));

		//flann::save_to_file(indices,"result.hdf5","result");


		nearestIds.clear() ;
		nearestIds.resize( queryPoints.size() ) ;
		for( int i=0; i<queryPoints.size(); ++i ){
			for( int j=0; j<nn; ++j )
				nearestIds[i].push_back( indices[i][j] ) ;
		}

		delete[] dataset.ptr();
		delete[] query.ptr();
		delete[] indices.ptr();
		delete[] dists.ptr();

	}

	void MycomputeFlannKNN(std::vector<Point3f> &dataPoints,std::vector<Point3f> &queryPoints,std::vector<Point3f> &nmls, std::vector<int> &nearestId){
		//add the influence of the vector
	float  k=0.5;
	nearestId.resize(queryPoints.size());
	
	for(int pid =0 ; pid<queryPoints.size() ; pid ++){
		float tmpDis=1;
		float minDis=1000;
		
		for(int spid = 0; spid < dataPoints.size() ; spid++){
			Point3f vecPtoSp = queryPoints[pid]-dataPoints[spid];
			float vcos=(vecPtoSp*nmls[pid])/(vecPtoSp.Norm()*nmls[pid].Norm());
			/*float vcos_2;
			if(vecPtoSp*nmls[pid]>=0){
				vcos_2=-sqrtf((vcos+1)/2);
			}else{
				vcos_2=sqrtf((vcos+1)/2);
			}*/
			tmpDis = (queryPoints[pid]-dataPoints[spid]).Norm()*(1-k*vcos);;//*(1+(1-vcos*vcos))
			if(tmpDis < minDis){
				minDis = tmpDis;
				nearestId[pid]=spid;
			}
		}
	}
	}

}