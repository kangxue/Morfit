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



#include "skelpath.h"

void skelpath::generateMesh() {

	resultMesh.clear() ;

	int psize = ReconstructorPara::finalProfilesSampleNum ;
	for( int pwbid = 0; pwbid<resultProf3d.size(); ++pwbid ) {

		//if( pwbid!=0 )
		//	continue ;

		std::vector<Point3f> vertices ;
		std::vector<Point3f> normals ;
		std::vector<int3> faces ;

		std::vector<Profile3D> &profs = resultProf3d[pwbid] ;

		for( int i=0; i<psize; ++i )
			vertices.push_back( profs[0][i] ) ;


		bool firstIsTip = solvers[pwbid].firstIsTip ;
		bool lastIsTip = solvers[pwbid].lastIsTip ;

		for( int i=0; i<profs.size()-1; ++i  ){

			// estimate offset
			std::vector<int> offsetCand ;
			for( int id=0; id<psize; id+=psize/10 ){
				int nstId = ReconstructorUtility::NearstPoint(profs[i+1], profs[i][id] ) ;
				offsetCand.push_back( (nstId - id + psize ) % psize ) ;
			}
			std::vector<int> diffSums ;
			for( int i=0; i<psize; ++i){
				double sum = 0;
				for( int i=0; i<offsetCand.size(); ++i )
					sum += (std::min)( (offsetCand[i] - i + psize) % psize , psize - (offsetCand[i] - i + psize) % psize ) ;
				diffSums.push_back(sum) ;
			}
			int maxSum = 0 ;  int offset=0;
			for( int i=0; i<diffSums.size(); ++i ){
				if( diffSums[i]>maxSum){ maxSum= diffSums[i]; offset = i;}
			}


			offset = 0 ;

			// generates triangles
			//  ------------ add vertices
			for( int id=0; id<psize; ++id )
				vertices.push_back( profs[i+1][id] ) ;

			while( normals.size()<vertices.size() )
				normals.push_back(Point3f(0,0,0)) ;



			for( int id=0; id<psize; ++id ){
				int id1 = (id+1)%psize ;

				int startId = i * psize ; 
				int startId1 = i * psize + psize ; 

				int3 t1 = int3(startId+id, startId + id1,                startId1 + (id1+offset)%psize ) ;
				int3 t2 = int3(startId+id, startId1 + (id1+offset)%psize, startId1 + (id+offset)%psize ) ;
				faces.push_back( t1 ) ;
				faces.push_back( t2 ) ;

				Point3f n1 = -(vertices[t1[1]] - vertices[t1[0]])^(vertices[t1[2]] - vertices[t1[1]]) ;
				Point3f n2 = -(vertices[t2[1]] - vertices[t2[0]])^(vertices[t2[2]] - vertices[t2[1]]) ;
				n1.Normalize() ;
				n2.Normalize() ;

				if( !(i==0&&firstIsTip) ){
					normals[t1[0]] += n1 ;     
					normals[t1[1]] += n1 ;
					normals[t1[2]] += n1 ;
				}

				if( !(i==profs.size()-1&&lastIsTip) ){
					normals[t2[0]] += n2 ;
					normals[t2[1]] += n2 ;
					normals[t2[2]] += n2 ;
				}


			}


		}

		if( brch0IsLoop ){
			// connect the first and last profile

			// estimate offset
			std::vector<int> offsetCand ;
			for( int id=0; id<psize; id+=psize/10 ){
				int nstId = ReconstructorUtility::NearstPoint(profs[0], profs.back()[id] ) ;
				offsetCand.push_back( (nstId - id + psize ) % psize ) ;
			}
			std::vector<int> diffSums ;
			for( int i=0; i<psize; ++i){
				double sum = 0;
				for( int i=0; i<offsetCand.size(); ++i )
					sum += (std::min)( (offsetCand[i] - i + psize) % psize , psize - (offsetCand[i] - i + psize) % psize ) ;
				diffSums.push_back(sum) ;
			}
			int maxSum = 0 ;  int offset=0;
			for( int i=0; i<diffSums.size(); ++i ){
				if( diffSums[i]>maxSum){ maxSum= diffSums[i]; offset = i;}
			}


			// generates triangles
			for( int id=0; id<psize; ++id ){
				int id1 = (id+1)%psize ;

				int startId  = (profs.size()-1) * psize ; 
				int startId1 = 0; 

				int3 t1 = int3(startId+id, startId + id1,                startId1 + (id1+offset)%psize ) ;
				int3 t2 = int3(startId+id, startId1 + (id1+offset)%psize, startId1 + (id+offset)%psize ) ;
				faces.push_back( t1 ) ;
				faces.push_back( t2 ) ;

			}


		}

		// calculate normal
		for( int i=0; i<profs.size(); ++i  ){
			for( int j=0; j<profs[i].size(); ++j ){

				if( i>0 && i<profs.size()-1 ){
					Point3f v1 = ( profs[i][(j+1)%psize] - profs[i][(j-1+psize)%psize] ).Normalize();
					Point3f v2 = ( profs[i+1][j] - profs[i-1][j] ).Normalize();

					normals[i * psize + j] = -v1^v2 ;
				}

				if( i==0&&firstIsTip ){
					Point3f n = solvers[pwbid].centers[0] -  solvers[pwbid].centers[1] ;
					normals[i * psize + j] = n ;
				}

				if( i==profs.size()-1&&lastIsTip ){
					Point3f n = solvers[pwbid].centers.back() -  solvers[pwbid].centers[ solvers[pwbid].centers.size()-2 ] ;
					normals[i * psize + j] = n ;
				}

			}

		}


		// update tip normal
		if( firstIsTip ){
			Point3f n = solvers[pwbid].centers[0] -  solvers[pwbid].centers[1] ;
			for( int i=0; i<psize; ++i)
				normals[i] = n ;
		}

		if( lastIsTip ){
			Point3f n = solvers[pwbid].centers.back() -  solvers[pwbid].centers[ solvers[pwbid].centers.size()-2 ] ;
			for( int i=0; i<psize; ++i)
				normals[ normals.size()-i-1 ] = n ;
		}


		// delete redundent tip vertices
		if( firstIsTip ){
			for( int i=0; i<psize-1; ++i) { normals.erase(normals.begin()); vertices.erase( vertices.begin()); }

			for( int i=0; i< faces.size(); ++i){
				if( faces[i].X() >= psize )  faces[i].X() -= psize-1 ;
				else  faces[i].X() = 0 ;

				if( faces[i].Y() >= psize )  faces[i].Y() -= psize-1 ;
				else  faces[i].Y() = 0 ;

				if( faces[i].Z() >= psize )  faces[i].Z() -= psize-1 ;
				else  faces[i].Z() = 0 ;
			}
		}
		if( lastIsTip ){
			for( int i=0; i<psize-1; ++i) { normals.erase(normals.begin() + normals.size()-1); vertices.erase(vertices.begin() + vertices.size()-1); }
			for( int i=0; i< faces.size(); ++i){
				if( faces[i].X() >= normals.size()-1 )  faces[i].X() =  normals.size()-1 ;
				if( faces[i].Y() >= normals.size()-1 )  faces[i].Y() =  normals.size()-1 ;
				if( faces[i].Z() >= normals.size()-1 )  faces[i].Z() =  normals.size()-1 ;
			}
		}


		// add endface ,,  do not share vertices with the extisting surface 
		if( !firstIsTip && !brch0IsLoop){
			Point3f cter = GlobalFun::centerOfPoints( profs[0] ) ;
			Point3f n = solvers[pwbid].centers[0] -  solvers[pwbid].centers[1] ;
			vertices.push_back( cter ) ;  int posmark = vertices.size();
			normals.push_back( n ) ;
			for( int i=0; i<psize; ++i ){
				vertices.push_back( vertices[i] ) ;
				normals.push_back( n ) ;
			}
			for( int i=0; i<psize; ++i )
				faces.push_back( int3(posmark-1, posmark+(i+1)%psize , posmark+i ) ) ;

		}
		if( !lastIsTip &&!brch0IsLoop ){
			Point3f cter = GlobalFun::centerOfPoints( profs.back() ) ;
			Point3f n = solvers[pwbid].centers.back() -  solvers[pwbid].centers[ solvers[pwbid].centers.size()-2 ] ;
			vertices.push_back( cter ) ;  int posmark = vertices.size();
			normals.push_back( n ) ;
			int pid = profs.size()-1 ;
			for( int i=0; i<psize; ++i ){
				vertices.push_back( vertices[pid*psize + i] ) ;
				normals.push_back( n ) ;
			}
			for( int i=0; i<psize; ++i )
				faces.push_back( int3(posmark-1 , posmark+i , posmark+(i+1)%psize) ) ;

		}

		// normalize normals
		for( int i=0; i<normals.size(); ++i ){
			normals[i].Normalize() ;
			//normals[i] *= -1 ;
		}

		// convert the vertices and faces to  c-style array

		triangleList TL( vertices, normals,faces ) ;

		resultMesh.push_back(TL) ;
	}

}
