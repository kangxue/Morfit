
#include "GLDrawer.h"

#include <fstream>
extern std::ofstream logout ;




bool samewithprecolor( GLColor randcolor[32], int i ){

	for( int j=0; j<i; ++j ){
		if( (Point3f(randcolor[i].r, randcolor[i].g,randcolor[i].b ) - Point3f(randcolor[j].r, randcolor[j].g,randcolor[j].b )  ).Norm() < 0.2 )
			return true ;
	}

	return false ; 
}
GLDrawer::GLDrawer(RichParameterSet* _para)
{
	para = _para;

	randColor[0] = GLColor( 0.1,1.0, 0.5) ;
	randColor[1] = GLColor( 0.1,0.5,0.5) ;
	randColor[2] = GLColor( 0.1,0.5,1.0 ) ;
	randColor[3] = GLColor( 0.4,0.0,1.0) ;
	randColor[4] = GLColor( 0.4,1.0,0.5) ;
	randColor[5] = GLColor( 0.4,0,1.0) ;
	randColor[6] = GLColor( 0.7,0.5,0.5) ;
	randColor[7] = GLColor( 0.7,1.0,0.0) ;
	randColor[8] = GLColor( 0.7,0.0,1.0) ;
	randColor[9] = GLColor( 1.0,1.0,0.1) ;

	srand(time(0)) ;
	for( int i=10; i<32; ++i ){
		randColor[i] = GLColor( (rand()%1000)/1000.0, (rand()%1000)/1000.0, (rand()%1000)/1000.0  ) ;



		GLColor color0 = para->getColor("Pick Point Color"); 
		GLColor color1 =  para->getColor("Sample Point Color") ;

		Point3f p0 ( color0.r, color0.g, color0.b ) ;
		Point3f p1 ( color1.r, color1.g, color1.b ) ;
		Point3f p2 ( 0.8, 0.8,  0.8 ) ;

		while( (Point3f(randColor[i].r, randColor[i].g,randColor[i].b ) - p0  ).Norm() < 0.2 
			|| (Point3f(randColor[i].r, randColor[i].g,randColor[i].b ) - p1 ).Norm() < 0.2 
			|| (Point3f(randColor[i].r, randColor[i].g,randColor[i].b ) - p2 ).Norm() < 0.2  
			|| samewithprecolor( randColor, i )  
			|| randColor[i].r > 0.7 ){


				randColor[i] = GLColor( (rand()%1000)/1000.0, (rand()%1000)/1000.0, (rand()%1000)/1000.0  ) ;


		}



	}

}


GLDrawer::~GLDrawer(void)
{
}

void GLDrawer::updateDrawer(vector<int>& pickList)
{
	bCullFace = para->getBool("Need Cull Points");

	original_draw_width = para->getDouble("Original Draw Width");
	sample_draw_width = para->getDouble("Sample Draw Width");
	Quade_draw_width = para->getDouble("Quade Draw Width") ;

	original_color = para->getColor("Original Point Color");
	sample_color = para->getColor("Sample Point Color");
	normal_width = para->getDouble("Normal Line Width");
	normal_length  = para->getDouble("Normal Line Length");
	sample_dot_size = para->getDouble("Sample Dot Size");
	original_dot_size = para->getDouble("Original Dot Size");
	normal_color = para->getColor("Normal Line Color");
	feature_color = para->getColor("Feature Color");
	pick_color = para->getColor("Pick Point Color");

	skel_bone_color = para->getColor("Skeleton Bone Color");
	skel_node_color = para->getColor("Skeleton Node Color");
	skel_branch_color = para->getColor("Skeleton Branch Color");
	skel_bone_width = para->getDouble("Skeleton Bone Width") / 10000.;
	skel_node_size = para->getDouble("Skeleton Node Size") / 10000.;
	skel_branch_size = para->getDouble("Skeleton Branch Size") / 10000.;


	if (!pickList.empty())
	{
		curr_pick_indx = pickList[0];
	}



}



void GLDrawer::draw(DrawType type, CMesh* _mesh, bool isSamples, std::vector<bool> &sampleIsSelected )
{
	if (!_mesh)
	{
		return;
	}

	if( sampleIsSelected.size() && isSamples )
		drawSelectedPoint(_mesh,  sampleIsSelected, true ) ;

	bool doPick = para->getBool("Doing Pick");

	int qcnt = 0;
	CMesh::VertexIterator vi;
	int i=-1;
	for(vi = _mesh->vert.begin(); vi != _mesh->vert.end(); ++vi) 
	{
		i++ ;

		if (doPick)
		{
			glLoadName(qcnt);
		}

		Point3f& p = vi->P();      
		Point3f& normal = vi->N();

		if(!bCullFace ||  isCanSee(p, normal) || !isSamples )
		{
			switch(type)
			{
			case DOT:
				if( isSamples && (sampleIsSelected.size()==_mesh->vert.size()) && sampleIsSelected[i] ){
				}
				else{

					drawDot(*vi);
				}

				break;
			case CIRCLE:
				drawCircle(*vi);
				break;
			case QUADE:
				drawQuade(*vi);
				break;
			case NORMAL:
				drawNormal(*vi);
				break;
			case SPHERE:
				drawSphere(*vi);
				break;
			default:
				break;
			}
		}

		if (doPick) 
		{
			qcnt++;
		}
	}

	para->setValue("Doing Pick", BoolValue(false));
}

void GLDrawer::DrawHeightMap(  CMesh* _mesh, std::vector<bool> &sampleIsSelected ){

	int size = sample_dot_size;


	Box3f bbox = _mesh->bbox ;

	glPointSize(size);

	for( int i=0; i< _mesh->vert.size(); ++i ){


		Point3f p = _mesh->vert[i].P();
		double height = (p[1]-bbox.min.Y() ) / bbox.DimY() ;
		double3 color = GlobalFun::scalar2color( height ) ;
		glColor3f(color.x,  color.y,  color.z ) ;
		glBegin(GL_POINTS);
		glVertex3f(p[0], p[1], p[2]);
		glEnd(); 

	}

}

void  GLDrawer::drawSampleDiskOnScreen(  CMesh* _mesh, std::vector<bool> &sampleIsSelected, std::vector<bool> &sampleIsBoudary,  std::vector<int> &sampleClusterId ){


	// inverse matrix of modle view matrix
	float modelviewMatrix[16] ;
	float invmvMatrix[16] ;
	glGetFloatv(GL_MODELVIEW_MATRIX, modelviewMatrix);
	GlobalFun::gluInvertMatrix( modelviewMatrix, invmvMatrix );

	glEnable(GL_LIGHTING) ;
	glEnable(GL_LIGHT0) ;
	glDisable(GL_LIGHT1) ;
	glColor4f(0.8, 0.8,0.8, 1.0 ) ;

	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY );

	int npoints = _mesh->vert.size() ;
	int nTrianglePerPoint = 16 ;
	int nvertices =  npoints * nTrianglePerPoint * 3  ;
	float *vertices = new float[ nvertices * 3] ;
	float *normals = new float[ nvertices * 3] ;
	float *colors = new float[ nvertices * 3] ;
	unsigned *indices = new unsigned[  nvertices ] ;

	for( int i=0; i<npoints; ++i){
		int firstVer_offset = i * nTrianglePerPoint * 3 * 3 ;

		Point3f center = _mesh->vert[i].P() ;
		Point3f normal = _mesh->vert[i].N() ;

		// writeNormal for each point
		for( int offset = firstVer_offset; offset<firstVer_offset+nTrianglePerPoint * 3*3; ++offset )
			normals[offset] = normal[(offset-firstVer_offset)%3] ;



		//vertices[offset]
		Point3f c_world = GlobalFun::vecMulMatrix(center, modelviewMatrix ) ;
		Point3f c_scrn = c_world * fabs( 0.1/c_world.Z() ) ;


		Point3f arcp1_scrn = c_scrn ;  arcp1_scrn.Y() += 0.1155 / 200 ;
		Point3f arcp1_world = arcp1_scrn * fabs( c_world.Z() / 0.1 ) ;

		//Point3f arcp2_scrn = c_scrn ;  arcp2_scrn.X() += 0.1155 / 50 ;
		//Point3f arcp2_world = arcp2_scrn * fabs( c_world.Z() / 0.1 ) ;


		Point3f arcp1 =  GlobalFun::vecMulMatrix( arcp1_world, invmvMatrix ) ;

		// write vertices for each point
		for( int tid =0; tid<nTrianglePerPoint; ++tid ){
			int offset = firstVer_offset + tid * 3 * 3 ;
			
			
			Point3f V1 , V2;
			GlobalFun::Rotate_Point3D(tid/(double)nTrianglePerPoint * 3.14159 * 2, Point3f(0,0,1), arcp1_scrn-c_scrn, V1 );
			GlobalFun::Rotate_Point3D(  ((tid+1)%nTrianglePerPoint)/(double)nTrianglePerPoint * 3.14159 * 2, Point3f(0,0,1), arcp1_scrn-c_scrn, V2 );


			V1 = (c_scrn + V1) *  fabs( c_world.Z() / 0.1 ) ;
			V2 = (c_scrn + V2) *  fabs( c_world.Z() / 0.1 ) ;



			V1 =  GlobalFun::vecMulMatrix( V1, invmvMatrix ) ;
			V2 =  GlobalFun::vecMulMatrix( V2, invmvMatrix ) ;



			vertices[offset+0] = center[0] ;
			vertices[offset+1] = center[1] ;
			vertices[offset+2] = center[2] ;
			vertices[offset+3] = V2[0] ;
			vertices[offset+4] = V2[1] ;
			vertices[offset+5] = V2[2] ;
			vertices[offset+6] = V1[0] ;
			vertices[offset+7] = V1[1] ;
			vertices[offset+8] = V1[2] ;

		}

	}


	// set indices
	for( unsigned i=0; i<nvertices; ++i )
		indices[i] = i ;


	// set colors
	if( sampleIsSelected.size() == npoints ){
		for( int i=0; i<nvertices; ++i ){
			if( sampleIsSelected[i/(nTrianglePerPoint * 3) ]   && !sampleIsBoudary[i/(nTrianglePerPoint * 3)]  ){
				colors[i*3] = 1.0 ; colors[i*3+1] = 0.0 ; colors[i*3+2] = 1.0 ;
			}else if( !sampleIsSelected[i/(nTrianglePerPoint * 3) ]   && !sampleIsBoudary[i/(nTrianglePerPoint * 3)]  ){
				colors[i*3] = 0.8 ; colors[i*3+1] = 0.8 ; colors[i*3+2] = 0.8 ;
			}else{

				if( sampleClusterId.size() == npoints ){
					int cid = sampleClusterId[ i/(nTrianglePerPoint * 3) ] %32;
					colors[i*3  ] = randColor[cid].r ;
					colors[i*3+1] = randColor[cid].g ;
					colors[i*3+2] = randColor[cid].b ;

				}
				else{
					colors[i*3] = 0 ; colors[i*3+1] = 1 ; colors[i*3+2] = 1 ;
				}
			}
		}
	}
	else{

		for( int i=0; i<nvertices*3; ++i )
			colors[i] = 0.8 ;
	}



	glVertexPointer(3,GL_FLOAT, 0, vertices ) ;
	glNormalPointer( GL_FLOAT, 0, normals ) ;
	glColorPointer(3, GL_FLOAT, 0, colors ) ;

	glDrawElements( GL_TRIANGLES, nvertices, GL_UNSIGNED_INT, indices) ;


	// release memory
	delete []vertices ;
	delete []normals;
	delete []indices ;
	delete []colors ;

	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);


}


bool drawPolygonOnScreenNeedUpdate = true ;
void  GLDrawer::drawPolygonOnScreen(  CMesh* _mesh, double length , int edgeNum ){

	if( edgeNum < 3 ) edgeNum = 3 ;
	if( edgeNum > 16 ) edgeNum = 16 ;


	bool drawPolygonOnScreenNeedUpdate = true ;

	// inverse matrix of modle view matrix
	float modelviewMatrix[16] ;
	float invmvMatrix[16] ;
	glGetFloatv(GL_MODELVIEW_MATRIX, modelviewMatrix);
	GlobalFun::gluInvertMatrix( modelviewMatrix, invmvMatrix );

	glEnable(GL_LIGHTING) ;
	glDisable(GL_LIGHT0) ;
	glEnable(GL_LIGHT1) ;
	glColor4f(1, 1,1, 1.0 ) ;

	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY );

	int npoints = _mesh->vert.size() ;
	int nTrianglePerPoint = edgeNum ;
	int nvertices =  npoints * nTrianglePerPoint * 3  ;
	
	static float *vertices = NULL ;
	static float *normals  = NULL;
	static float *colors = NULL ;
	static unsigned *indices = NULL ;

	if( drawPolygonOnScreenNeedUpdate ){

		//std::cout<<"drawPolygonOnScreenNeedUpdate = true";

		// release memory
		if( vertices!=NULL) 	delete []vertices ;
		if( normals!=NULL) 	delete []normals ;
		if( colors!=NULL) 	delete []colors ;
		if( indices!=NULL) 	delete []indices ;

		vertices = new float[ nvertices * 3] ;
		normals = new float[ nvertices * 3] ;
		colors = new float[ nvertices * 3] ;
		indices = new unsigned[  nvertices ] ;

		for( int i=0; i<npoints; ++i){
			int firstVer_offset = i * nTrianglePerPoint * 3 * 3 ;

			Point3f center = _mesh->vert[i].P() ;
			Point3f normal = _mesh->vert[i].N() ;

			// writeNormal for each point
			for( int offset = firstVer_offset; offset<firstVer_offset+nTrianglePerPoint * 3*3; ++offset )
				normals[offset] = normal[(offset-firstVer_offset)%3] ;



			//vertices[offset]
			Point3f c_world = GlobalFun::vecMulMatrix(center, modelviewMatrix ) ;
			Point3f c_scrn = c_world * fabs( 0.1/c_world.Z() ) ;


			Point3f arcp1_scrn = c_scrn ;  arcp1_scrn.Y() += length / 200 ;
			Point3f arcp1_world = arcp1_scrn * fabs( c_world.Z() / 0.1 ) ;

			//Point3f arcp2_scrn = c_scrn ;  arcp2_scrn.X() += 0.1155 / 50 ;
			//Point3f arcp2_world = arcp2_scrn * fabs( c_world.Z() / 0.1 ) ;


			Point3f arcp1 =  GlobalFun::vecMulMatrix( arcp1_world, invmvMatrix ) ;

			// write vertices for each point
			for( int tid =0; tid<nTrianglePerPoint; ++tid ){
				int offset = firstVer_offset + tid * 3 * 3 ;


				Point3f V1 , V2;
				GlobalFun::Rotate_Point3D(tid/(double)nTrianglePerPoint * 3.14159 * 2, Point3f(0,0,1), arcp1_scrn-c_scrn, V1 );
				GlobalFun::Rotate_Point3D(  ((tid+1)%nTrianglePerPoint)/(double)nTrianglePerPoint * 3.14159 * 2, Point3f(0,0,1), arcp1_scrn-c_scrn, V2 );


				V1 = (c_scrn + V1) *  fabs( c_world.Z() / 0.1 ) ;
				V2 = (c_scrn + V2) *  fabs( c_world.Z() / 0.1 ) ;



				V1 =  GlobalFun::vecMulMatrix( V1, invmvMatrix ) ;
				V2 =  GlobalFun::vecMulMatrix( V2, invmvMatrix ) ;



				vertices[offset+0] = center[0] ;
				vertices[offset+1] = center[1] ;
				vertices[offset+2] = center[2] ;
				vertices[offset+3] = V2[0] ;
				vertices[offset+4] = V2[1] ;
				vertices[offset+5] = V2[2] ;
				vertices[offset+6] = V1[0] ;
				vertices[offset+7] = V1[1] ;
				vertices[offset+8] = V1[2] ;
			}
		}

		// set indices
		for( unsigned i=0; i<nvertices; ++i )
			indices[i] = i ;


		// set colors
		for( int i=0; i<nvertices*3; ++i )
			colors[i] = 0.8 ;

		//drawPolygonOnScreenNeedUpdate = false ;
	}



	//for( int i=0; i<nvertices*3; ++i )
	//	normals[i] *=-1 ;

	glVertexPointer(3,GL_FLOAT, 0, vertices ) ;
	glNormalPointer( GL_FLOAT, 0, normals ) ;
	glColorPointer(3, GL_FLOAT, 0, colors ) ;

	glDrawElements( GL_TRIANGLES, nvertices, GL_UNSIGNED_INT, indices) ;



	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);


}
void  GLDrawer::drawDisk(  CMesh* _mesh ){




	glEnable(GL_LIGHTING) ;
	glEnable(GL_LIGHT0) ;
	glDisable(GL_LIGHT1) ;
	glColor4f(0.8, 0.8,0.8, 1.0 ) ;

	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);

	int npoints = _mesh->vert.size() ;
	int nTrianglePerPoint = 16 ;
	int nvertices =  npoints * nTrianglePerPoint * 3  ;
	float *vertices = new float[ nvertices * 3] ;
	float *normals = new float[ nvertices * 3] ;
	unsigned *indices = new unsigned[  nvertices ] ;

	for( int i=0; i<npoints; ++i){
		int firstVer_offset = i * nTrianglePerPoint * 3 * 3 ;

		Point3f center = _mesh->vert[i].P() ;
		Point3f normal = _mesh->vert[i].N() ;

		// writeNormal for each point
		for( int offset = firstVer_offset; offset<firstVer_offset+nTrianglePerPoint * 3*3; ++offset )
			normals[offset] = normal[(offset-firstVer_offset)%3] ;


		// write vertices for each point
		for( int tid =0; tid<nTrianglePerPoint; ++tid ){
			int offset = firstVer_offset + tid * 3 * 3 ;
			//vertices[offset]



			Point3f V0 = _mesh->vert[i].eigen_vector0;
			V0.Normalize();
			V0 = V0*0.02 ;


			Point3f V1 , V2;
			GlobalFun::Rotate_Point3D(tid/(double)nTrianglePerPoint * 3.14159 * 2, normal, V0, V1 );
			GlobalFun::Rotate_Point3D(  ((tid+1)%nTrianglePerPoint)/(double)nTrianglePerPoint * 3.14159 * 2, normal, V0, V2 );

			V1 = center + V1 ;
			V2 = center + V2 ;

			vertices[offset+0] = center[0] ;
			vertices[offset+1] = center[1] ;
			vertices[offset+2] = center[2] ;
			vertices[offset+3] = V2[0] ;
			vertices[offset+4] = V2[1] ;
			vertices[offset+5] = V2[2] ;
			vertices[offset+6] = V1[0] ;
			vertices[offset+7] = V1[1] ;
			vertices[offset+8] = V1[2] ;

		}

	}

	for( unsigned i=0; i<nvertices; ++i )
		indices[i] = i ;


	glVertexPointer(3,GL_FLOAT, 0, vertices ) ;
	glNormalPointer( GL_FLOAT, 0, normals ) ;
	
	glDrawElements( GL_TRIANGLES, nvertices, GL_UNSIGNED_INT, indices) ;


	// release memory
	delete []vertices ;
	delete []normals;
	delete []indices ;

	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);


}


bool GLDrawer::isCanSee(const Point3f& pos, const Point3f& normal)
{
	return  ( (view_point - pos) * normal >= 0 );
}

GLColor GLDrawer::getColorByType(const CVertex& v)
{
	//return original_color;

	if (v.bIsOriginal)
	{
		return original_color;
	}
	else
	{
		if (v.is_fixed_sample)
		{
			return feature_color;
		}
		return sample_color;
	}
}

void GLDrawer::drawDot(const CVertex& v)
{
	int size;
	if (v.bIsOriginal)
	{
		size = original_dot_size;
	}
	else
	{
		size = sample_dot_size;
	}

	GLColor color = getColorByType(v);
	glColor4f(color.r, color.g, color.b, 1);

	glPointSize(size);
	glBegin(GL_POINTS);

	Point3f p = v.P();
	glVertex3f(p[0], p[1], p[2]);
	glEnd(); 
}

void GLDrawer::drawSphere(const CVertex& v)
{

	Point3f p = v.P(); 

	glDrawSphere( p, GLColor(1.0, 1.0, 0.0, 1.0), 0.01, 5 ) ;

}

void GLDrawer::drawCircle(const CVertex& v)
{

}

void GLDrawer::drawQuade(const CVertex& v)
{
	if (v.is_skel_ignore)
	{
		return;
	}


	double h = Quade_draw_width;
	GLColor color = getColorByType(v);
	glColor4f(color.r, color.g, color.b, 1);

	Point3f p = v.P(); 
	Point3f normal = v.cN();
	Point3f V0 = v.eigen_vector0;
	Point3f V1 = v.eigen_vector1;

	Point3f p0 = p + V0 * h;
	Point3f p1 = p + V1 * h;
	Point3f p2 = p + (-V0) * h;
	Point3f p3 = p + (-V1) * h;


	//const double *m = v.m;
	glBegin(GL_QUADS);

	glNormal3d(normal.X(), normal.Y(), normal.Z());
	glVertex3d(p0.X(), p0.Y(), p0.Z());
	glVertex3d(p1.X(), p1.Y(), p1.Z());
	glVertex3d(p2.X(), p2.Y(), p2.Z());
	glVertex3d(p3.X(), p3.Y(), p3.Z());

	glEnd();  

}

void GLDrawer::drawNormal(const CVertex& v)
{
	double width = normal_width;
	double length = normal_length;
	QColor qcolor = normal_color;

	glDisable(GL_LIGHTING);

	glLineWidth(width); 
	GLColor color(qcolor);
	glColor4f(color.r, color.g, color.b, 0.8);  

	Point3f p = v.P(); 
	Point3f m = v.cN();

	glBegin(GL_LINES);	
	glVertex3d(p[0], p[1], p[2]);
	glVertex3f(p[0] + m[0]*length, p[1]+m[1]*length, p[2]+m[2]*length);
	glEnd(); 


	//glBegin(GL_LINES);	
	//glVertex3d(p[0], p[1], p[2]);
	//glVertex3f(p[0] - m[0]*length, p[1] - m[1]*length, p[2] - m[2]*length);
	//glEnd(); 

	glEnable(GL_LIGHTING);
}


void GLDrawer::drawSelectedPoint(CMesh* samples, vector<bool>& sampleIsSelected, bool bShow_as_dot)
{

	double width = para->getDouble("Sample Draw Width");
	GLColor pick_color = para->getColor("Pick Point Color");
	glColor3f(pick_color.r, pick_color.g, pick_color.b);

	for(int i = 0; i < sampleIsSelected.size(); i++) 
	{
		if( !sampleIsSelected[i] )
			continue ;

		CVertex &v = samples->vert[i];
		Point3f &p = v.P();     

		if(bShow_as_dot)
		{
			glPointSize(sample_dot_size * 2);
			glBegin(GL_POINTS);

			GLColor color = pick_color;
			glColor4f(color.r, color.g, color.b, 1);

			glVertex3d(p[0], p[1], p[2]);

			glEnd(); 
		}
		else
		{
			glDrawSphere(p, pick_color, sample_draw_width * 2, 40);
		}
	}    
}


void GLDrawer::drawPickPoint(CMesh* samples, vector<int>& pickList, bool bShow_as_dot)
{
	double width = para->getDouble("Sample Draw Width");
	GLColor pick_color = para->getColor("Pick Point Color");
	glColor3f(pick_color.r, pick_color.g, pick_color.b);

	for(int ii = 0; ii < pickList.size(); ii++) 
	{
		int i = pickList[ii];

		if(i < 0 || i >= samples->vert.size())
			continue;

		CVertex &v = samples->vert[i];
		Point3f &p = v.P();     

		if(bShow_as_dot)
		{
			glPointSize(sample_dot_size * 2);
			glBegin(GL_POINTS);

			GLColor color = pick_color;
			glColor4f(color.r, color.g, color.b, 1);

			glVertex3d(p[0], p[1], p[2]);

			glEnd(); 
		}
		else
		{
			glDrawSphere(p, pick_color, sample_draw_width * 2, 40);
		}
	}    
}

void GLDrawer::glDrawLine(Point3f& p0, Point3f& p1, GLColor color, double width)
{
	glColor3f(color.r, color.g, color.b);
	glLineWidth(width);
	glBegin(GL_LINES);
	glVertex3f(p0[0], p0[1], p0[2]);
	glVertex3f(p1[0], p1[1], p1[2]);
	glEnd();
}


void RenderBone(float x0, float y0, float z0, float x1, float y1, float z1, double width = 20)  
{  
	GLdouble  dir_x = x1 - x0;  
	GLdouble  dir_y = y1 - y0;  
	GLdouble  dir_z = z1 - z0;  
	GLdouble  bone_length = sqrt( dir_x*dir_x + dir_y*dir_y + dir_z*dir_z );  
	static GLUquadricObj *  quad_obj = NULL;  
	if ( quad_obj == NULL )  
		quad_obj = gluNewQuadric();  
	gluQuadricDrawStyle( quad_obj, GLU_FILL );  
	gluQuadricNormals( quad_obj, GLU_SMOOTH );  
	glPushMatrix();  
	// 平移到起始点   
	glTranslated( x0, y0, z0 );  
	// 计算长度   
	double  length;  
	length = sqrt( dir_x*dir_x + dir_y*dir_y + dir_z*dir_z );  
	if ( length < 0.0001 ) {   
		dir_x = 0.0; dir_y = 0.0; dir_z = 1.0;  length = 1.0;  
	}  
	dir_x /= length;  dir_y /= length;  dir_z /= length;  
	GLdouble  up_x, up_y, up_z;  
	up_x = 0.0;  
	up_y = 1.0;  
	up_z = 0.0;  
	double  side_x, side_y, side_z;  
	side_x = up_y * dir_z - up_z * dir_y;  
	side_y = up_z * dir_x - up_x * dir_z;  
	side_z = up_x * dir_y - up_y * dir_x;  
	length = sqrt( side_x*side_x + side_y*side_y + side_z*side_z );  
	if ( length < 0.0001 ) {  
		side_x = 1.0; side_y = 0.0; side_z = 0.0;  length = 1.0;  
	}  
	side_x /= length;  side_y /= length;  side_z /= length;  
	up_x = dir_y * side_z - dir_z * side_y;  
	up_y = dir_z * side_x - dir_x * side_z;  
	up_z = dir_x * side_y - dir_y * side_x;  
	// 计算变换矩阵   
	GLdouble  m[16] = { side_x, side_y, side_z, 0.0,  
		up_x,   up_y,   up_z,   0.0,  
		dir_x,  dir_y,  dir_z,  0.0,  
		0.0,    0.0,    0.0,    1.0 };  
	glMultMatrixd( m );  
	// 圆柱体参数   
	GLdouble radius= width;        // 半径   
	GLdouble slices = 40.0;      //  段数   
	GLdouble stack = 3.0;       // 递归次数   
	gluCylinder( quad_obj, radius, radius, bone_length, slices, stack );   
	glPopMatrix();  
}  

void GLDrawer::glDrawCylinder(Point3f& p0, Point3f& p1, GLColor color, double width)
{
	glColor3f(color.r, color.g, color.b);
	RenderBone(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], width);
}


void GLDrawer::glDrawSphere(Point3f& p, GLColor color, double radius, int slide)
{
	if (radius < 0.0001)
		radius = 0.01;

	glColor3f(color.r, color.g, color.b);
	glPushMatrix();      
	glTranslatef(p[0], p[1], p[2]);
	glutSolidSphere(radius, slide, slide);
	glPopMatrix();
}


void GLDrawer::cleanPickPoint()
{
	curr_pick_indx = 0;
	prevPickIndex = 0;
}

//
//void GLDrawer::glDrawBranches(vector<Branch>& branches, GLColor gl_color)
//{
//	for(int i = 0; i < branches.size(); i++)
//	{
//		Curve& curve = branches[i].curve;
//
//		if (curve.empty())
//		{
//			cout << "curve empty!" << endl;
//			continue;
//		}
//
//		if (curve.size() == 1)
//		{
//			Point3f p0 = curve[0];
//			glDrawSphere(p0, GLColor(cBlue), skel_node_size, 40);
//
//		}
//		Point3f p0, p1;
//		double size_scale = 0.5;
//		//draw nodes
//		for (int j = 0; j < curve.size(); j++)
//		{
//			Point3f p0 = p0 = curve[j].P();
//
//			if (curve[j].is_skel_virtual)
//			{
//				glDrawSphere(p0, GLColor(feature_color), skel_node_size * 1.1, 40);
//			}
//			else
//			{
//				glDrawSphere(p0, GLColor(skel_node_color), skel_node_size, 40);
//			}
//		}
//
//		for (int j = 0; j < curve.size()-1; j++)
//		{
//			p0 = curve[j].P();
//			p1 = curve[j+1].P();
//
//			if (skel_bone_width > 0.0003)
//			{
//				glDrawCylinder(p0, p1, gl_color, skel_bone_width);
//			}
//		}
//
//	}
//}


void GLDrawer::drawSkeleton(skeleton_mul &skel, GLColor gl_color ){

	glEnable(GL_COLOR_MATERIAL) ;


	//GLfloat mat_ambient[4] = {0.6, 0.6, 0.6,1.0}; 
	//GLfloat mat_diffuse[4] = {0.0, 0.6, 0.6, 1.0 };
	//GLfloat mat_specular[] = {0, 1.0, 1.0, 1.0 };
	//GLfloat shininess = 0.5*128;

	//glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_ambient);
	//glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse); 
	//glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	//glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);



	glEnable(GL_LIGHTING) ;
	glEnable(GL_LIGHT1) ;

	if(  skel.branches.size() == 0 || skel.nodes.size() == 0 )
		return ;

	for(int i = 0; i < skel.branches.size(); i++)
	{

		std::vector<Point3f> curve ;
		for( int pbid=0; pbid<skel.branches[i].size(); ++pbid  )
			curve.push_back( skel.nodes[ skel.branches[i][pbid] ] ) ;



		if (curve.size() == 1)
		{
			Point3f p0 = curve[0];
			glDrawSphere(p0, GLColor(cBlue), skel_node_size, 40);


		}
		Point3f p0, p1;
		double size_scale = 0.5;
		//draw nodes
		for (int j = 0; j < curve.size(); j++)
		{
			Point3f p0 = p0 = curve[j];

			if( skel.brchNodeDegree[i][j] > 2 )
				// glDrawSphere(p0, GLColor(feature_color), skel_node_size * 1.1, 40);
				glDrawSphere(p0, GLColor(skel_node_color), skel_node_size, 40);

		}

		for (int j = 0; j < curve.size()-1; j++)
		{
			p0 = curve[j];
			p1 = curve[j+1];

			if (skel_bone_width > 0.0003)
			{
				glDrawCylinder(p0, p1, gl_color, skel_bone_width);
			}
		}

	}


	GLfloat mat_ambient_[4] = {0.2, 0.2, 0.2,1.0}; 
	GLfloat mat_diffuse_[4] = {0.8, 0.8, 0.8, 1.0 };
	GLfloat mat_specular_[] = {0.0, 0.0, 0.0, 1.0 };
	GLfloat shininess_ = 0;

	//glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_ambient_);
	//glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse_); 
	//glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular_);
	//glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess_);

	glDisable(GL_LIGHT1) ;

}

//void GLDrawer::glDrawCurves(vector<Curve>& curves, GLColor gl_color)
//{
//	for(int i = 0; i < curves.size(); i++)
//	{
//		Curve& curve = curves[i];
//
//		if (curve.empty())
//		{
//			cout << "curve empty!" << endl;
//			continue;
//		}
//
//		if (curve.size() == 1)
//		{
//			Point3f p0 = curve[0];
//			glDrawSphere(p0, GLColor(cBlue), skel_node_size, 40);
//
//		}
//		Point3f p0, p1;
//		double size_scale = 0.5;
//		//draw nodes
//		for (int j = 0; j < curve.size(); j++)
//		{
//			Point3f p0 = p0 = curve[j].P();
//
//			if (curve[j].is_skel_virtual)
//			{
//				glDrawSphere(p0, GLColor(feature_color), skel_node_size * 1.1, 40);
//			}
//			else
//			{
//				glDrawSphere(p0, GLColor(skel_node_color), skel_node_size, 40);
//			}
//		}
//
//		for (int j = 0; j < curve.size()-1; j++)
//		{
//			p0 = curve[j].P();
//			p1 = curve[j+1].P();
//
//			if (skel_bone_width > 0.0003)
//			{
//				glDrawCylinder(p0, p1, gl_color, skel_bone_width);
//			}
//		}
//
//	}
//}


//void GLDrawer::drawCurveSkeleton(Skeleton& skeleton)
//{
//	if (para->getBool("Skeleton Light"))
//	{
//		glEnable(GL_LIGHTING);
//	}
//
//	glDrawBranches(skeleton.branches, skel_bone_color);
//
//	if (para->getBool("Skeleton Light"))
//	{
//		//glDisable(GL_LIGHTING);
//	}
//}


void GLDrawer::drawSamplesAsDot(CMesh* _mesh, std::vector<bool> &sampleIsSelected,std::vector<bool> &sampleIsBoundary, std::vector<int> &sampleClusterId ){

	int n = _mesh->vert.size() ;
	if( sampleIsSelected.size() != n || sampleIsBoundary.size() != n ){

		int size = sample_dot_size ;
		GLColor color = sample_color ;
		glColor4f(color.r, color.g, color.b, 1);
		glPointSize(size);
		glBegin(GL_POINTS);
		for( int i=0; i< n; ++i ){
				Point3f p = _mesh->vert[i].P();
				glVertex3f(p[0], p[1], p[2]);
		}
		glEnd(); 
				
		return ;

	}

	// draw unselected points except for boundary
	int size = sample_dot_size ;
	GLColor color = sample_color ;
	glColor4f(color.r, color.g, color.b, 1);
	glPointSize(size);
	glBegin(GL_POINTS);
	for( int i=0; i< n; ++i ){
		if(  !sampleIsSelected[i] && !sampleIsBoundary[i] ){
			Point3f p = _mesh->vert[i].P();
			glVertex3f(p[0], p[1], p[2]);
		}
	}
	glEnd(); 


	// draw selected points except for boundary
	color = para->getColor("Pick Point Color"); ;
	glColor4f(color.r, color.g, color.b, 1);
	glColor4f( 0.8, 0.8, 0.8, 1 ) ;
	glPointSize(size);
	glBegin(GL_POINTS);
	for( int i=0; i< n; ++i ){
		if(  sampleIsSelected[i] &&  !sampleIsBoundary[i]  ){
			Point3f p = _mesh->vert[i].P();
			glVertex3f(p[0], p[1], p[2]);
		}
	}
	glEnd(); 


	// draw boundary
	glPointSize(size);
	glColor4f( 0, 1, 1, 1 ) ;
	for( int i=0; i< n; ++i ){
		if(  sampleIsBoundary[i]  ){
			Point3f p = _mesh->vert[i].P();

			if( sampleClusterId.size() == n ){
				int cid =  sampleClusterId[i] %32;
				glColor4f( randColor[cid].r,randColor[cid].g, randColor[cid].b, 1  ) ;
			}

			glBegin(GL_POINTS);
			glVertex3f(p[0], p[1], p[2]);
			glEnd(); 
		}
	}


} 