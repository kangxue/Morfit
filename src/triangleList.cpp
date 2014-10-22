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



#include <gl/glew.h>
#include "triangleList.h"
//#include <QOpenGLShader>
#include "GlobalFunction.h"


#include <fstream>
extern std::ofstream logout ;

triangleList::triangleList( int nvertices, int nfaces){
	vertices = new float[ nvertices * 3] ;
	indices = new unsigned[  nfaces*3 ] ;

	this->nvertices = nvertices ;
	this->nfaces = nfaces ;
}

triangleList::triangleList( std::vector<Point3f> ver,std::vector<Point3f> nor,  std::vector<int3> faces ){
	nvertices = ver.size() ;
	nfaces = faces.size() ;

	vertices = new float[ nvertices * 3] ;
	normals = new float[ nvertices * 3] ;
	indices = new unsigned[  nfaces * 3] ;

	// set vertices
	for( int i=0; i<ver.size(); ++i ){
		vertices[i*3+0] = ver[i][0] ;
		vertices[i*3+1] = ver[i][1] ;
		vertices[i*3+2] = ver[i][2] ;
	}
	// set normals
	for( int i=0; i<nor.size(); ++i ){
		normals[i*3+0] = nor[i][0] ;
		normals[i*3+1] = nor[i][1] ;
		normals[i*3+2] = nor[i][2] ;
	}
	// set indices
	for( unsigned i=0; i<faces.size(); ++i ){
		indices[i*3+0] = faces[i][0] ;
		indices[i*3+1] = faces[i][1] ;
		indices[i*3+2] = faces[i][2] ;
	}


}

triangleList::triangleList( CMesh &mesh) {

	nvertices = mesh.vert.size() ;
	nfaces = mesh.face.size() ;

	for( int i=0; i<mesh.vert.size(); ++i )
		mesh.vert[i].m_index = i ;

	vertices = new float[ nvertices * 3] ;
	normals = new float[ nvertices * 3] ;
	indices = new unsigned[  nfaces * 3] ;

	// set vertices
	for( int i=0; i<nvertices; ++i ){
		vertices[i*3+0] = mesh.vert[i].P().X() ;
		vertices[i*3+1] = mesh.vert[i].P().Y() ;
		vertices[i*3+2] = mesh.vert[i].P().Z() ;
	}
	// set normals
	for( int i=0; i<nvertices; ++i ){
		normals[i*3+0] = mesh.vert[i].N().X() ;
		normals[i*3+1] = mesh.vert[i].N().Y() ;
		normals[i*3+2] = mesh.vert[i].N().Z() ;
	}
	// set indices
	for( unsigned i=0; i<nfaces; ++i ){

		indices[i*3+2] = mesh.face[i].V(0)->m_index ; 
		indices[i*3+1] = mesh.face[i].V(1)->m_index ; 
		indices[i*3+0] = mesh.face[i].V(2)->m_index ; 
	}


	//for( int i=0; i<100; ++i )
	//	std::cout << vertices[i] <<" " ;
	//for( int i=0; i<100; ++i )
	//	std::cout << normals[i] <<" " ;
	//for( int i=0; i<100; ++i )
	//	std::cout << indices[i] <<" " ;


	

}

CMesh triangleList::cvtToCmesh() {

		
	CMesh m;
	for( int i=0; i<nvertices; ++i ){
		CMesh::VertexIterator vi = vcg::tri::Allocator<CMesh>::AddVertices(m,1);
		vi->P()=CMesh::CoordType ( vertices[i* 3 + 0] , vertices[ i * 3 + 1] , vertices[ i * 3 + 2] ); 
		vi->N()=CMesh::CoordType ( normals [i* 3 + 0] , normals[ i * 3 + 1] , normals[ i * 3 + 2] ); 
	}

	for( int i=0; i<nfaces; ++i ){

		CMesh::FaceIterator fi = vcg::tri::Allocator<CMesh>::AddFaces(m,1);
		CMesh::VertexPointer ivp[3];
		ivp[2]= &( m.vert[indices[i*3+0]] ) ;
		ivp[1]= &( m.vert[indices[i*3+1]] ) ;
		ivp[0]= &( m.vert[indices[i*3+2]] ) ;
		fi->V(0)=ivp[0];
		fi->V(1)=ivp[1];
		fi->V(2)=ivp[2];
	}

	m.vn = m.vert.size() ;
	ReconstructorUtility::CopyCmesh(cmeshFormat, m) ;

	return m ;

}


void triangleList::drawItself(){
	if( nfaces == 0 ) 
		return ;
//	glDisable(GL_POLYGON_SMOOTH);

	glEnable(GL_LIGHTING) ;
	glEnable(GL_LIGHT0) ;
	glEnable(GL_LIGHT1) ;

	glEnableClientState(GL_VERTEX_ARRAY);
//	glEnableClientState(GL_COLOR_ARRAY );
	glEnableClientState(GL_NORMAL_ARRAY);

	glVertexPointer(3,GL_FLOAT, 0, vertices ) ;
	//glColorPointer(3, GL_FLOAT, 0, colors ) ;
	glNormalPointer( GL_FLOAT, 0, normals ) ;
	glDrawElements( GL_TRIANGLES, nfaces, GL_UNSIGNED_INT, indices) ;

	glDisableClientState(GL_VERTEX_ARRAY);
//	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glColor4f(0,0,0,1) ;
	glEnable(GL_LIGHTING) ;
	glEnable(GL_LIGHT0) ;
	glDisable(GL_LIGHT1) ;
	glBegin(GL_TRIANGLES) ;
	for( int i=0; i<nfaces*3; ++i ){		
		glNormal3f(normals[ indices[i]*3 ],normals[ indices[i]*3+1 ],normals[ indices[i]*3+2 ] ) ;
		glVertex3f(vertices[ indices[i]*3 ],vertices[ indices[i]*3+1 ],vertices[ indices[i]*3+2 ] ) ;
	}
	glEnd() ;

	glDisable(GL_LIGHTING) ;
}



void triangleList::drawItself( std::vector<double3> colors ){
	if( nfaces == 0 ) 
		return ;

	if( colors.size() != nvertices )
		return ;




	glEnable(GL_LIGHTING) ;
	glEnable(GL_LIGHT0) ;
	glEnable(GL_LIGHT1) ;
	glEnable(GL_LIGHT2) ;
	//glEnable(GL_LIGHT3) ;

	glBegin(GL_TRIANGLES) ;
	for( int i=0; i<nfaces*3; ++i ){		
		glNormal3f(normals[ indices[i]*3 ],normals[ indices[i]*3+1 ],normals[ indices[i]*3+2 ] ) ;
		glColor3f ( colors[ indices[i] ].x, colors[ indices[i] ].y,  colors[ indices[i] ].z ) ;
		glVertex3f(vertices[ indices[i]*3 ],vertices[ indices[i]*3+1 ],vertices[ indices[i]*3+2 ] ) ;
	}
	glEnd() ;

	
	glDisable(GL_LIGHT0) ;
	glDisable(GL_LIGHT1) ;
	glDisable(GL_LIGHT2) ;
	glDisable(GL_LIGHT3) ;

	glDisable(GL_LIGHTING) ;

	//CMesh mesh1 ;

	//for( int i=0; i<nvertices; ++i )
	//	mesh1.vert.


}


void triangleList::destroy(){

	if( nvertices!=0){
		delete[] vertices ;
		delete[] normals ;
	}
	if( nfaces!= 0)
		delete[] indices ;


	nvertices = 0;
	nfaces = 0;
	
	cvtToCmesh() ;
}



void triangleList::addMesh( triangleList &a) {


	float *newvertices = new float[( nvertices + a.nvertices) * 3] ;
	float *newnormals = new float[( nvertices + a.nvertices) * 3] ;
	unsigned *newindices = new unsigned[  ( nfaces + a.nfaces)*3 ] ;

	for( int i=0; i<nvertices*3; ++i ){
		newvertices[i] = vertices[i] ;
		newnormals[i] = normals[i] ;
	}
	for( int i=0; i<nfaces*3; ++i )
		newindices[i] = indices[i] ;

	for( int i=0; i<a.nvertices*3; ++i ){
		newvertices[nvertices*3+i] = a.vertices[i] ;
		newnormals[nvertices*3+i] = a.normals[i] ;
	}
	for( int i=0; i<a.nfaces*3; ++i )
		newindices[nfaces*3+i] = a.indices[i] + nvertices;

	if( nvertices > 0 ){
		delete[] vertices ;
		delete[] indices ;
	}
	if( nfaces  > 0 )	
		delete[] normals ;

	vertices = newvertices ;
	normals = newnormals ;
	indices = newindices ;

	nvertices += a.nvertices ;
	nfaces += a.nfaces ;

	cvtToCmesh() ;
}



void triangleList::sortDepth() {

	//sort depth and calculate and buffer the result in data member "indices"
	if( nfaces == 0 || nvertices==0)
		return ;

	//std::cout << "triangleList::sortDepth(){ "<< std::endl;

	extern GLdouble globalModelviewMatrix[16] ;
	extern GLdouble globalProjectMatrix[16] ;
	float mvmatrix[16];   for( int i=0; i<16;++i ) mvmatrix[i]=globalModelviewMatrix[i] ;
	float pjmatrix[16];   for( int i=0; i<16;++i ) pjmatrix[i]=globalProjectMatrix[i] ;

	
	std::vector<double> vdepths(nvertices) ;
	for( int i=0; i<nvertices; ++i )
		vdepths[i] =  GlobalFun::vecMulMatrix( Point3f(vertices[i*3],vertices[i*3+1],vertices[i*3+2] ), mvmatrix ).Z() ;


	std::vector<std::pair<double,int>> triDepthsPair(nfaces) ;
	for( int i=0; i<nfaces; ++i ){
		triDepthsPair[i].first = -( vdepths[indices[i*3]] + vdepths[indices[i*3+1]] + vdepths[indices[i*3+2]] )/3;
		triDepthsPair[i].second = i ;
	}


	std::sort( triDepthsPair.begin(), triDepthsPair.end() ) ;

	//for( int i=0; i<triDepthsPair.size() ; ++i)
	//	logout << triDepthsPair[i].first <<std::endl;
	//system("pause" );

	unsigned *newindices = new unsigned[ nfaces*3 ] ;
	for( int i=0; i<nfaces; ++i ){
		newindices[i*3]   = indices[triDepthsPair[i].second *3 ] ;
		newindices[i*3+1] = indices[triDepthsPair[i].second *3 +1] ;
		newindices[i*3+2] = indices[triDepthsPair[i].second *3 +2] ;
	}

	delete indices ;
	indices = newindices ;
}



void triangleList::operator=(triangleList &tl) {



	ReconstructorUtility::CopyCmesh( cmeshFormat,tl.cmeshFormat) ;
	vertices = tl.vertices ;
	normals = tl.normals ;
	indices = tl.indices ;
	nvertices = tl.nvertices ;
	nfaces = tl.nfaces ;

}



std::vector<Point3f> triangleList::getVertices() {

	std::vector<Point3f> vlst ;

	for( int i=0; i<nvertices; ++i  )
		vlst.push_back( Point3f(  vertices[i*3],  vertices[i*3+1],  vertices[i*3+2] ) ) ;
	

	return vlst ;
}