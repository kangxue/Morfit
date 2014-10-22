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



#ifndef triangleList_h
#define  triangleList_h


#include "Types.h"
#include <GL/GL.h>
#include <GL/GLU.h>

class triangleList{
public:

	triangleList(){
		nvertices = nfaces = 0;
		vertices = normals  = NULL ;
		indices = NULL ;
	}
	triangleList( int nvertices, int nfaces) ;
	triangleList( std::vector<Point3f> ver,  std::vector<Point3f> nor,std::vector<int3> faces ) ;
	triangleList( CMesh &mesh);

	void operator=(triangleList &tl) ;

	// draw the triangle list
	void drawItself();
	void drawItself( std::vector<double3> colors ) ;

	//free the room
	void destroy() ;

	//add another triangleList in this mesh
	void addMesh( triangleList &a) ;

	//sort depth and calculate and buffer the result in data member "indices"
	void sortDepth() ;

	//convert the data members mesh 
	CMesh cvtToCmesh() ;

	//return vertices in the form of vector<Point3f>
	std::vector<Point3f> getVertices() ; 

	CMesh cmeshFormat ;

	float *vertices  ;
	float *normals ;
	unsigned *indices  ;
	int nvertices ;
	int nfaces ;
} ;
#endif