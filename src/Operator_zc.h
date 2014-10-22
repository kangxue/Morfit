
#ifndef _Operator_zc_h_ 
#define _Operator_zc_h_ 

#include "GLArea.h"
#include <GL/glu.h>
#include "skeleton_mul.h"
#include <QCursor>
#include "bspline.h"
#include <math.h>
#include "reconstructorPara.h"
#include "reconstructionUtility.h"

#define PI 3.14159265

extern GLdouble globalProjectMatrix[16] ;
extern GLdouble globalModelviewMatrix[16] ;

static const int maxPixelsForPoint=10;
static const int maxPixelsForEdge=10;
static const double minCurvatureForCircle=0.02;
static const int maxPixelsFor3DegPoint=10;


class Operator_zc
{
private:

	static enum MoveTendEnum{tendToDefault,tendToCutSegment,tendToCutCrossPoint,tendToJointPoints,tendToDelBranch,tendToAttach,tendToExtend};

	typedef lemon::ListGraph::Node Node ;
	typedef lemon::ListGraph::NodeIt NodeIt ;
	typedef lemon::ListGraph::Edge Edge ;
	typedef lemon::ListGraph::EdgeIt EdgeIt ;
	typedef lemon::ListGraph::IncEdgeIt IncEdgeIt ;

	typedef lemon::ListGraph::NodeMap<Point3f> P3fNodeMap ;

	vector<lemon::ListGraph*> undoGraph;
	vector<lemon::ListGraph*> redoGraph;
	vector<P3fNodeMap*> redoNodeMap;
	vector<P3fNodeMap*> undoNodeMap;


	class M22{
	public:
		M22( double a, double b, double c, double d){
			m[0][0]=a; m[0][1] = b;  m[1][0] = c; m[1][1] = d ;
		}
		M22(){m[0][0]=0.0; m[0][1] = 0.0;  m[1][0] = 0.0; m[1][1] = 0.0 ;}

		double m[2][2] ;

	};

	Point2f calFootPoint(Point2f a,Point2f b,Point2f c)
	{
		float abx = b.X() - a.X();
		float aby = b.Y() - a.Y();
		float acx = c.X() - a.X();
		float acy = c.Y() - a.Y();
		float f = (abx*acx+aby*acy )/(abx*abx+aby*aby);  // 注意ab必须是直线上的两个不同点
		return Point2f(a.X() + f*abx,a.Y() + f*aby);
	}


	// 功能：判断点是否在多边形内 
	// 方法：求解通过该点的水平线与多边形各边的交点 
	// 结论：单边交点为奇数，成立!
	//参数： 
	// POINT p 指定的某个点 
	// LPPOINT ptPolygon 多边形的各个顶点坐标（首末点可以不一致） 
	// int nCount 多边形定点的个数

	bool PtInPolygon (Point2f p, std::vector<Point2f> &ptPolygon ) 
	{ 
		int nCross = 0;
		int nCount = ptPolygon.size() ;

		for (int i = 0; i < nCount; i++) 
		{ 
			Point2f p1 = ptPolygon[i]; 
			Point2f p2 = ptPolygon[(i + 1) % nCount];
			// 求解 y=p.y 与 p1p2 的交点
			if ( p1.Y() == p2.Y()) // p1p2 与 y=p0.y平行 
				continue;
			if ( p.Y()< min(p1.Y(), p2.Y()) ) // 交点在p1p2延长线上 
				continue; 
			if ( p.Y()>= max(p1.Y(), p2.Y()) ) // 交点在p1p2延长线上 
				continue;
			// 求交点的 X 坐标 -------------------------------------------------------------- 
			double x = (double)(p.Y()- p1.Y()) * (double)(p2.X() - p1.X()) / (double)(p2.Y() - p1.Y()) + p1.X();
			if ( x > p.X() ) 
				nCross++; // 只统计单边交点 
		}
		// 单边交点为偶数，点在多边形之外 --- 
		return (nCross % 2 == 1); 
	}
	void getBernsteinCoefficient( std::vector<double> &bc, int n   ){

		bc.resize(n);

		n=n-1;

		int j,k;
		for (k=0;k<=n;k++)
		{ //compute n! / (k!*(n-k)!)
			bc[k] = 1;
			for (j = n;j>=k+1;j--){
				bc[k] *= (double)j / (double)(j-k);
			}
		}
	}

	// calculate curvatures 
	std::vector<double> calAverCurvature( std::vector<std::vector<Point2f>> curves ) {

		std::vector<double> res ;
		for( int cid=0; cid<curves.size(); ++cid ){

			double total_curvature = 0;

			std::vector<double2> curve ;
			for( int i=0; i<curves[cid].size(); ++i)
				curve.push_back( double2( curves[cid][i].X(),curves[cid][i].Y() )  ) ;

			int n = curve.size();
			std::vector<double> bc ;
			getBernsteinCoefficient(bc, n ) ;

			// store the Bezier curve
			std::vector<double2> bezierCurves ;
			for( double u=0; u<=1.0; u+=0.01 ){
				double2 p ;
				for( int k = 0; k<n; ++k ){
					p.x += bc[k] * pow(u, k) * pow( 1-u, n-1-k) * curve[k].x ;
					p.y += bc[k] * pow(u, k) * pow( 1-u, n-1-k) * curve[k].y ;
				}	

				bezierCurves.push_back( p ) ;
			}

			Curve2D curveP2f,bezierCurvesP2f ;
			for( int i=0; i<curve.size(); ++i )
				curveP2f.push_back( Point2f( curve[i].x , curve[i].y ) ) ;

			ReconstructorUtility::getUniformCubicBSpline(curveP2f, bezierCurvesP2f, 100) ;

			bezierCurves.clear() ;
			for( int i=0; i<bezierCurvesP2f.size(); ++i )
				bezierCurves.push_back( double2(bezierCurvesP2f[i].X(), bezierCurvesP2f[i].Y() ) ) ;


			for( int id =1; id<bezierCurves.size()-1; ++id ){
				double2 p1 = bezierCurves[id-1] ;
				double2 p2 = bezierCurves[id  ] ;
				double2 p3 = bezierCurves[id+1] ;

				total_curvature +=  fabs( acos( (std::min)( (std::max)( (p3-p2)*(p2-p1)/( (p3-p2).norm()*(p2-p1).norm()), -1.0), 1.0 )  ) / (p3-p1).norm()  ) ;


				if( _isnan(total_curvature) || !_finite( total_curvature ) ){


					std::cout<<"calAverCurvature NaN:\n"<<p1.x<<" "<<p1.y<<"\n"
						<<p2.x<<" "<<p2.y<<"\n"
						<<p3.x<<" "<<p3.y<<std::endl;
					// system("pause");
				}

			}

			res.push_back(total_curvature/curve.size()) ;

		}

		return res ;

	}

	double cross(const Point2f& a, const Point2f& b, const Point2f& c)
	{
		return (a.X()- c.X()) * (b.Y() - c.Y()) - (a.Y() - c.Y()) * (b.X() - c.X());
	}

	double ex(double a) 
	{ 
		const double eps = 1e-8;
		return (a < eps && a > -eps) ? 0 : (a > 0 ? 1 : -1);
	}

	bool between(const Point2f& a, const Point2f& c, const Point2f& b)
	{
		return (c.X()>= min(a.X(), b.X()) && c.X() <= max(a.X(), b.X())) && (c.Y() >= min(a.Y(), c.Y()) && c.Y() <= max(a.Y(), c.Y()));
	}

	int intersect(const Point2f& a, const Point2f& b, const Point2f& c, const Point2f& d, Point2f &res)
	{
		double s1,s2,s3,s4;
		int d1 = ex(s1=cross(a, b, c));
		int d2 = ex(s2=cross(a, b, d));
		int d3 = ex(s3=cross(c, d, a));
		int d4 = ex(s4=cross(c, d, b)); 
		if(((d1 ^ d2) == -2 && (d3 ^ d4) == -2)) // 规范相交 -1 ^ 1 = -2
		{
			res=Point2f((c.X()*s2-d.X()*s1)/(s2-s1),(c.Y()*s2-d.Y()*s1)/(s2-s1));
			return 1;
		}
		// 非规范相交 
		else if(d1 == 0 && between(a, c, b))
		{
			res=c;
			if(d2==0 && c!=a && c!=b)return 3;
			else return 2;
		}
		else if(d2 == 0 && between(a, d, b))
		{
			res=d;
			if(d1 == 0 && d!=a && d!=b)return 3;
			else return 2;
		}
		else if(d3 == 0 && between(c, a, d))
		{
			res=a;
			if(d4 == 0 && a!=c && a!=d)return 3;
			else return 2;
		}
		else if(d4 == 0 && between(c, b, d))
		{
			res=b;
			if(d3 == 0 && b!=c && b!=d)return 3;
			else return 2;
		}
		else return 0;
	}

	float GetPointDistance(Point2f p1, Point2f p2)    
	{   
		return sqrt((p1.X()-p2.X())*(p1.X()-p2.X())+(p1.Y()-p2.Y())*(p1.Y()-p2.Y()));   
	}   

	float GetNearestDistance(Point2f PA, Point2f PB, Point2f P3)   
	{   

		//----------图2--------------------   
		float a,b,c;   
		a=GetPointDistance(PB,P3);   
		if(a<=0.00001)   
			return 0.0f;   
		b=GetPointDistance(PA,P3);   
		if(b<=0.00001)   
			return 0.0f;   
		c=GetPointDistance(PA,PB);   
		if(c<=0.00001)   
			return a;//如果PA和PB坐标相同，则退出函数，并返回距离   
		//------------------------------   


		if(a*a>=b*b+c*c)//--------图3--------   
			return b;   
		if(b*b>=a*a+c*c)//--------图4-------   
			return a;    

		//图1   
		float l=(a+b+c)/2;     //周长的一半   
		float s=sqrt(l*(l-a)*(l-b)*(l-c));  //海伦公式求面积   
		return 2*s/c;   
	}

	Point2f myProject(Point3f p3f)
	{
		GLint viewport[4];
		GLdouble winX, winY,winZ;
		glGetIntegerv(GL_VIEWPORT, viewport);
		gluProject(p3f.X(),p3f.Y(),p3f.Z(),globalModelviewMatrix,globalProjectMatrix,viewport,&winX,&winY,&winZ);
		return Point2f(winX,winY);
	}

	Point2f myProject(Point3f p3f,GLdouble &winZ)
	{
		GLint viewport[4];
		GLdouble winX, winY;
		glGetIntegerv(GL_VIEWPORT, viewport);
		gluProject(p3f.X(),p3f.Y(),p3f.Z(),globalModelviewMatrix,globalProjectMatrix,viewport,&winX,&winY,&winZ);
		return Point2f(winX,winY);
	}
	Point2f myProject(skeleton_mul skel,Node node,GLdouble &winZ)
	{
		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		Point3f p3f=nodeMap[node];

		GLint viewport[4];
		GLdouble winX, winY;
		glGetIntegerv(GL_VIEWPORT, viewport);
		gluProject(p3f.X(),p3f.Y(),p3f.Z(),globalModelviewMatrix,globalProjectMatrix,viewport,&winX,&winY,&winZ);
		return Point2f(winX,winY);
	}
	Point3f myUnProject(GLdouble winX,GLdouble winY,GLdouble winZ)
	{
		GLdouble posX,posY,posZ;
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);
		gluUnProject(winX, winY, winZ,globalModelviewMatrix,globalProjectMatrix,viewport, &posX, &posY, &posZ); 
		return Point3f(posX,posY,posZ);
	}
	bool belongToSameBranch(skeleton_mul skel,Edge e1,Edge e2)
	{
		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		if(e1==e2)return true;

		Edge lastEdge=e1;
		Node lastNode=g.u(lastEdge);
		while(lemon::countIncEdges(g,lastNode)==2)
		{
			for(IncEdgeIt iei=IncEdgeIt(g,lastNode);iei!=lemon::INVALID;iei++)
			{
				if(lastEdge!=iei)
				{
					if(e2==iei)return true;
					lastEdge=iei;
					lastNode=g.oppositeNode(lastNode,iei);
					break;
				}
			}
		}

		lastEdge=e1;
		lastNode=g.v(lastEdge);
		while(lemon::countIncEdges(g,lastNode)==2)
		{
			for(IncEdgeIt iei=IncEdgeIt(g,lastNode);iei!=lemon::INVALID;iei++)
			{
				if(lastEdge!=iei)
				{
					if(e2==iei)return true;
					lastEdge=iei;
					lastNode=g.oppositeNode(lastNode,iei);
					break;
				}
			}
		}

		return false;
	}
	Edge isOnBranch(Point2f p2f,skeleton_mul skel)
	{
		float minDis;
		Edge nearestEdge=Edge(lemon::INVALID);

		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		for(EdgeIt eit=EdgeIt(g);eit!=lemon::INVALID;eit++)
		{
			float distmp=GetNearestDistance(myProject(nodeMap[g.u(eit)]),myProject(nodeMap[g.v(eit)]),p2f);
			if(nearestEdge==Edge(lemon::INVALID)||minDis>distmp)
			{
				minDis=distmp;
				nearestEdge=Edge(eit);
			}
		}
		if(nearestEdge!=Edge(lemon::INVALID)&&minDis<maxPixelsForPoint)
			return nearestEdge;
		else return lemon::INVALID;
	}
	Node isOnEndPoint(Point2f p2f,skeleton_mul skel)
	{
		float minDis;
		Node nearestNode=Node(lemon::INVALID);

		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		for(NodeIt nit=NodeIt(g);nit!=lemon::INVALID;nit++)
		{
			if(lemon::countIncEdges(g,nit)==1||lemon::countIncEdges(g,nit)>2)
			{
				float distmp=GetPointDistance(myProject(nodeMap[nit]),p2f);
				if(nearestNode==Node(lemon::INVALID)||minDis>distmp)
				{
					minDis=distmp;
					nearestNode=Node(nit);
				}
			}
		}
		if(nearestNode!=Node(lemon::INVALID)&&minDis<maxPixelsForEdge)
			return nearestNode;
		else return lemon::INVALID;
	}
	Node isCut3DegNode(skeleton_mul skel,vector<Point2f> track)
	{
		float minDis;
		Node nearestNode=Node(lemon::INVALID);

		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		for(NodeIt nit=NodeIt(g);nit!=lemon::INVALID;nit++)
		{
			if(lemon::countIncEdges(g,nit)>2)
			{
				for(int i=0;i<track.size()-1;i++)
				{
					float distmp=GetNearestDistance(track[i],track[i+1],myProject(nodeMap[nit]));
					if(nearestNode==Node(lemon::INVALID)||minDis>distmp)
					{
						minDis=distmp;
						nearestNode=Node(nit);
					}
				}
			}
		}
		if(nearestNode!=Node(lemon::INVALID)&&minDis<maxPixelsFor3DegPoint)
			return nearestNode;
		else return lemon::INVALID;
	}
	int justOneBranch(skeleton_mul skel,vector<Point2f> track)
	{
		int hasOne=-1;
		for(int i=0;i<skel.branches.size();i++)
		{
			for(int j=0;j<skel.branches[i].size();j++)
			{
				if(PtInPolygon(myProject(skel.nodes[skel.branches[i][j]]),track))
				{
					if(hasOne!=-1)return -1;
					else hasOne=i;
					break;
				}
			}
		}
		return hasOne;
	}
	bool moreThanOneEndPoints(skeleton_mul skel,vector<Point2f> track)
	{
		Node firstEndPoint=lemon::INVALID;

		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		for(NodeIt nit=NodeIt(g);nit!=lemon::INVALID;nit++)
		{
			if(lemon::countIncEdges(g,nit)==1||lemon::countIncEdges(g,nit)>2)
			{
				if(PtInPolygon(myProject(nodeMap[nit]),track))
				{
					if(firstEndPoint==lemon::INVALID)
						firstEndPoint=nit;
					else if(nit!=firstEndPoint)return true;
					else continue;
				}
			}
		}
		return false;
	}
	Edge justCutOne(skeleton_mul skel,vector<Point2f> track)
	{
		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		Edge toBeDel=Edge(lemon::INVALID);

		for(EdgeIt e=EdgeIt(g);e!=lemon::INVALID;e++)
		{
			Point3f p1=nodeMap[g.u(e)],p2=nodeMap[g.v(e)];
			GLint viewport[4];   
			GLdouble winX[2], winY[2], winZ[2];   
			glGetIntegerv(GL_VIEWPORT, viewport);   
			gluProject(p1.X(),p1.Y(),p1.Z(),globalModelviewMatrix,globalProjectMatrix,viewport,winX,winY,winZ);
			gluProject(p2.X(),p2.Y(),p2.Z(),globalModelviewMatrix,globalProjectMatrix,viewport,winX+1,winY+1,winZ+1);

			//winY[0]=viewport[3]-winY[0];
			//winY[1]=viewport[3]-winY[1];

			for(int k=0;k<track.size()-1;k++)
			{
				Point2f res;
				int res_status;
				res_status=intersect(Point2f(winX[0],winY[0]),Point2f(winX[1],winY[1]),track[k],track[k+1],res);
				if(res_status)
				{
					if(toBeDel!=lemon::INVALID)return lemon::INVALID;
					else
					{
						toBeDel=e;
						break;
					}
				}
			}
		}
		return toBeDel;
	}
	vector<Node> pointsInCircle(skeleton_mul skel,vector<Point2f> track)
	{
		vector<Node> res;
		res.clear();

		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		for(NodeIt nit=NodeIt(g);nit!=lemon::INVALID;nit++)
		{
			if(lemon::countIncEdges(g,nit)==1||lemon::countIncEdges(g,nit)>2)
			{
				if(PtInPolygon(myProject(nodeMap[nit]),track))
				{
					for(int i=0;i<res.size();i++)
						if(nit==res[i])continue;
					res.push_back(Node(nit));
				}
			}
		}
		return res;
	}
	Node getPointOnEdge(skeleton_mul &skel,Edge e,Point2f p)
	{
		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		Node res;

		GLdouble d1,d2;
		Point2f u2d,v2d;
		Point2f footPoint;

		u2d=myProject(nodeMap[g.u(e)],d1);
		v2d=myProject(nodeMap[g.v(e)],d2);
		footPoint=calFootPoint(u2d,v2d,p);
		if(GetPointDistance(footPoint,u2d)<maxPixelsForEdge)
			res=g.u(e);
		else if(GetPointDistance(footPoint,v2d)<maxPixelsForEdge)
			res=g.v(e);
		else
		{
			nodeMap[res=g.addNode()]=myUnProject(footPoint.X(),footPoint.Y(),d1 + (d2-d1) * (footPoint.X()-u2d.X()) / (v2d.X()-u2d.X()) );
			g.addEdge(g.u(e),res);
			g.addEdge(res,g.v(e));
			g.erase(e);
		}
		return res;
	}
	MoveTendEnum moveTend_Left(skeleton_mul skel,vector<Point2f> track)
	{
		if(track.size()>2)
		{
			vector<vector<Point2f>> vvp2f;
			vector<double> vdres;
			vvp2f.clear();
			vvp2f.push_back(track);
			vdres.clear();
			vdres=calAverCurvature(vvp2f);

			Edge e1=isOnBranch(track[0],skel),e2=isOnBranch(track.back(),skel);
			if(e1!=lemon::INVALID&&e2!=lemon::INVALID)
				return tendToAttach;

			if(e1!=lemon::INVALID)
				return tendToExtend;

			if(vdres.back()>=minCurvatureForCircle&&justOneBranch(skel,track)!=-1)
				return tendToDelBranch;

			if(vdres.back()>=minCurvatureForCircle&&moreThanOneEndPoints(skel,track))
				return tendToJointPoints;

			Node node=isCut3DegNode(skel,track);
			if(node!=lemon::INVALID)
				return tendToCutCrossPoint;
			if(justCutOne(skel,track)!=lemon::INVALID)
				return tendToCutSegment;
		}

		return tendToDefault;

	}
	void attachBranches(skeleton_mul &skel,vector<Point2f> track)
	{
		printf("       attachBranches!\n");

		//vector<Point2f> track;
		//getUniformCubicBSpline(oriTrack,track);

		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		GLdouble startPointDepth,endPointDepth;
		Node startPoint,endPoint;

		if((startPoint=isOnEndPoint(track[0],skel))==lemon::INVALID)
			startPoint=getPointOnEdge(skel,isOnBranch(track[0],skel),track[0]);

		if((endPoint=isOnEndPoint(track.back(),skel))==lemon::INVALID)
			endPoint=getPointOnEdge(skel,isOnBranch(track.back(),skel),track.back());

		myProject(skel, startPoint,startPointDepth);
		myProject(skel, endPoint,endPointDepth);

		Node last=startPoint;
		Point3f sp=nodeMap[startPoint],ep=nodeMap[endPoint],np;

		for(int i=1;i<track.size()-1;i++)
		{
			Node ni;
			np=nodeMap[ni=g.addNode()]=myUnProject(track[i].X(),track[i].Y(),startPointDepth+(endPointDepth-startPointDepth)*i/(track.size()-1.0));
			g.addEdge(last,ni);
			last=ni;
		}
		g.addEdge(last,endPoint);
		skel.convertGraphToVerticeList();
		skel.convertVeticesListToBranchPoints();
		skel.smoothEachBranch();
	}
	void jointPoints(skeleton_mul &skel,vector<Point2f> track)
	{
		printf("       jointPoints!\n");

		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		vector<Node> points=pointsInCircle(skel,track);

		Point3f newpoint=nodeMap[points[0]];
		for(int i=1;i<points.size();i++)
			newpoint+=nodeMap[points[i]];
		newpoint/=points.size();
		Node newnode;
		nodeMap[newnode=g.addNode()]=newpoint;
		for(int i=0;i<points.size();i++)
			g.addEdge(newnode,points[i]);

		skel.convertGraphToVerticeList();
		skel.convertVeticesListToBranchPoints();
		skel.smoothEachBranch();
	}
	void deleteBranch(skeleton_mul &skel,vector<Point2f> track)
	{
		printf("       deleteBranch!\n");

		skel.convertGraphToVerticeList();
		skel.convertVeticesListToBranchPoints();
		skel.smoothEachBranch();
		skel.branches.erase(skel.branches.begin()+justOneBranch(skel,track));
		skel.convertVerticesListToGraph();
		skel.convertGraphToVerticeList();
		skel.convertVeticesListToBranchPoints();
		skel.smoothEachBranch();
	}
	void extendbranch(skeleton_mul &skel,vector<Point2f> track)
	{
		printf("       extendbranch!\n");

		//vector<Point2f> track;
		//getUniformCubicBSpline(oriTrack,track);

		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		int times=0;
		GLdouble startPointDepth,averRate=0;
		Node startPoint;
		Point2f startPoint2d;

		if((startPoint=isOnEndPoint(track[0],skel))==lemon::INVALID)
			startPoint=getPointOnEdge(skel,isOnBranch(track[0],skel),track[0]);

		startPoint2d=myProject(skel, startPoint,startPointDepth);

		for(IncEdgeIt e(g,startPoint);e!=lemon::INVALID;e++)
		{
			Point2f tmp2d;
			double tmpDepth;
			tmp2d=myProject(skel,g.oppositeNode(startPoint,e),tmpDepth);
			averRate+=(startPointDepth-tmpDepth)/(tmp2d-startPoint2d).Norm();
			times++;
		}
		averRate/=times;

		Node last=startPoint;
		Point2f last2d=startPoint2d;
		GLdouble lastDepth=startPointDepth;

		for(int i=1;i<track.size();i++)
		{
			Node ni;
			nodeMap[ni=g.addNode()]=myUnProject(track[i].X(),track[i].Y(),lastDepth+(track[i]-last2d).Norm()*averRate);
			g.addEdge(last,ni);
			last=ni;
		}
		skel.convertGraphToVerticeList();
		skel.convertVeticesListToBranchPoints();
		skel.smoothEachBranch() ;
	}
	void cutCrossPoint(skeleton_mul &skel,vector<Point2f> track)
	{
		printf("       cutCrossPoint!\n");

		Node node=isCut3DegNode(skel,track);
		skel.skeletonGraph->erase(node);
		skel.convertGraphToVerticeList();
		skel.convertVeticesListToBranchPoints();
		skel.smoothEachBranch() ;
	}
	void cutBranch(skeleton_mul &skel,vector<Point2f> track)
	{
		printf("       cutBranch!\n");

		Edge e=justCutOne(skel,track);
		skel.skeletonGraph->erase(e);
		skel.convertGraphToVerticeList();
		skel.convertVeticesListToBranchPoints();
		skel.smoothEachBranch();
	}
	void backup(skeleton_mul &skel)
	{
		while(!redoGraph.empty())
		{
			lemon::ListGraph *gi=redoGraph.back();
			redoGraph.pop_back();
			P3fNodeMap *mapi=redoNodeMap.back();
			redoNodeMap.pop_back();
			gi->clear();
			delete gi;
			delete mapi;
		}

		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		lemon::ListGraph *backup=new lemon::ListGraph;
		lemon::ListGraph::NodeMap<Node> gtb(g);

		for(NodeIt ni(g);ni!=lemon::INVALID;ni++)
			gtb[ni]=backup->addNode();

		P3fNodeMap *bnodeMap=new P3fNodeMap(*backup);
		for(NodeIt ni(g);ni!=lemon::INVALID;ni++)
			(*bnodeMap)[gtb[ni]]=nodeMap[ni];

		for(EdgeIt ei(g);ei!=lemon::INVALID;ei++)
			backup->addEdge(gtb[g.u(ei)],gtb[g.v(ei)]);

		undoGraph.push_back(backup);
		undoNodeMap.push_back(bnodeMap);
	}


	vector<Point2f> GussianDeform(vector<Point2f> pts, vector<int> steps, Point2f translation , int max_gussian_winside)	{

		int winside = max_gussian_winside ;
		if( winside > ReconstructorPara::winside_skel_deform )
			winside = ReconstructorPara::winside_skel_deform  ;

		double sigma = 0.3 * ( winside / 2.0 - 1.0 )  + 0.8  ;

		for( int i=0; i<pts.size(); ++i ){
			pts[i] += translation * pow(2.718, -steps[i]*steps[i]/(2.0*sigma*sigma) ) ;
		}


		return pts ;
	}


	vector<Point2f> MLSDeform(vector<Point2f> v,vector<Point2f> handless,vector<Point2f> dstpoints)
	{
		assert(handless.size() == dstpoints.size()) ;

		const double alpha = 1.0 ;

		vector<Point2f> p,q;
		for( int i=0; i<handless.size(); ++i){
			p.push_back(Point2f(handless[i].X(), handless[i].Y())) ;
			q.push_back(Point2f(dstpoints[i].X(), dstpoints[i].Y())) ;
		}


		double hsize = handless.size() ;

		// traverse all pixel in v, and store the result into vdst ;
		std::vector<Point2f> vdst(v.size()) ;
		for( int vid=0; vid<v.size(); ++vid ){

			// weight
			std::vector<double> w(hsize) ;
			for( int i=0; i<hsize; ++i ){
				w[i] = 1.0/( pow((double) ((p[i] -v[vid]).Norm()),2*alpha) );
			}


			// p*, q*
			Point2f p_as  = Point2f(0,0) ;
			Point2f q_as  = Point2f(0,0) ;
			for( int i=0; i<hsize; ++i ){

				p_as = p_as +  p[i] * w[i] ;
				q_as = q_as +  q[i] * w[i] ;
			}
			double wSum = 0;
			for( int i=0; i<hsize; ++i){
				wSum += w[i] ;
			}
			p_as = p_as * (1/wSum) ;
			q_as = q_as * (1/wSum) ;


			// p^ , q^
			std::vector<Point2f> p_tip( hsize );
			std::vector<Point2f> q_tip(  hsize );
			for( int i=0; i<hsize; ++i){
				p_tip[i] = p[i] - p_as ;
				q_tip[i] = q[i] - q_as ;
			}

			// A

			std::vector<M22> A(hsize) ;
			Point2f v_pas = v[vid] - p_as ; 
			for( int i=0; i<hsize; ++i){
				A[i].m[0][0] = w[i] * ( p_tip[i].X() * v_pas.X() + p_tip[i].Y() * v_pas.Y() ) ;
				A[i].m[0][1] = w[i] * ( p_tip[i].X() * v_pas.Y() - p_tip[i].Y() * v_pas.X() ) ;
				A[i].m[1][0] = w[i] * ( p_tip[i].Y() * v_pas.X() - p_tip[i].X() * v_pas.Y() ) ;
				A[i].m[1][1] = w[i] * ( p_tip[i].Y() * v_pas.Y() + p_tip[i].X() * v_pas.X() ) ;
			}


			// fr->
			Point2f fr_ar = Point2f(0,0);
			for( int i=0; i<hsize; ++i ){
				fr_ar = fr_ar + Point2f( q_tip[i] * Point2f(A[i].m[0][0],A[i].m[1][0] ), q_tip[i] * Point2f(A[i].m[0][1],A[i].m[1][1] ) ) ;
			}

			// fr(v)

			Point2f fr_v = fr_ar *( v_pas.Norm() / fr_ar.Norm() ) + q_as ; 


			vdst[vid] = fr_v ;


		}
		return vdst;
	}
	void deform(skeleton_mul &skel,Node endpoint)        //deform using MLS
	{
		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		vector<Point2f> oldHandles;
		vector<Point2f> newHandles;

		vector<Point2f> toTransform;  // points to transformed
		vector<int> steps;            // steps to the draged point
		vector<Point2f> transformed;
		vector<Node> nodelist;
		vector<GLdouble> depths;

		GLdouble depthtmp;

		toTransform.push_back(myProject(skel,endpoint,depthtmp));  
		steps.push_back( 0 ) ;
		depths.push_back(depthtmp);
		nodelist.push_back(endpoint);

		oldHandles.push_back(mouseTrack[0]);
		newHandles.push_back(mouseTrack[1]);


		int max_gussian_winside = 10000;

		for(IncEdgeIt iei(g,endpoint);iei!=lemon::INVALID;iei++)
		{
			Node lastNode=g.oppositeNode(endpoint,iei);
			Edge lastEdge=iei;

			int step_count = 1 ;
			toTransform.push_back(myProject(skel,lastNode,depthtmp));
			steps.push_back( step_count ) ;

			depths.push_back(depthtmp);
			nodelist.push_back(lastNode);

			while(lemon::countIncEdges(g,lastNode)==2)
			{
				step_count ++ ;

				for(IncEdgeIt iei2=IncEdgeIt(g,lastNode);iei2!=lemon::INVALID;iei2++)
				{
					if(lastEdge!=iei2)
					{
						lastEdge=iei2;
						lastNode=g.oppositeNode(lastNode,iei2);

						toTransform.push_back(myProject(skel,lastNode,depthtmp));
						steps.push_back( step_count ) ;
						nodelist.push_back(lastNode);
						depths.push_back(depthtmp);

						break;
					}
				}
			}

			if( step_count-1 < max_gussian_winside )
				max_gussian_winside = step_count-1 ;

			toTransform.pop_back();
			steps.pop_back() ;
			nodelist.pop_back();
			depths.pop_back();
			if(lastNode!=endpoint)
			{
				oldHandles.push_back(myProject(skel,lastNode,depthtmp));
				newHandles.push_back(oldHandles.back());
			}
		}

		transformed = MLSDeform(toTransform,oldHandles,newHandles);

		//transformed = GussianDeform(toTransform,steps, mouseTrack[1]-mouseTrack[0], max_gussian_winside );



		for(int i=0;i<toTransform.size();i++)
			nodeMap[nodelist[i]]=myUnProject(transformed[i].X(),transformed[i].Y(),depths[i]);

		skel.convertGraphToVerticeList();
		skel.convertVeticesListToBranchPoints();
		skel.smoothEachBranch();
	}

	double GlgetGaussianVal(double s)
	{
		const double SIGMA=1,MIU=0;
		return 1/sqrt(2*PI)/SIGMA*exp(-(s-MIU)*(s-MIU)/2/SIGMA/SIGMA);
	}

	void deformLocally(skeleton_mul &skel,Node endpoint)        //deform using Gaussian distribution (unfinished)
	{
		lemon::ListGraph & g = *(skel.skeletonGraph);
		P3fNodeMap &nodeMap=*(skel.graphNodeP3f);

		vector<Point2f> toTransform;
		vector<Point2f> transformed;
		vector<Node> nodelist;
		vector<GLdouble> depths;
		vector<double> svec;

		toTransform.clear();
		svec.clear();
		transformed.clear();
		nodelist.clear();
		depths.clear();

		Point2f varVec=mouseTrack[1]-mouseTrack[0];
		double maxGaussianVal=GlgetGaussianVal(0);
		GLdouble depthtmp;
		double minS=-1;

		toTransform.push_back(myProject(skel,endpoint,depthtmp));
		depths.push_back(depthtmp);
		nodelist.push_back(endpoint);
		svec.push_back(0.0);

		transformed.push_back(toTransform.back()+varVec);

		for(IncEdgeIt iei(g,endpoint);iei!=lemon::INVALID;iei++)
		{
			Node lastNode=g.oppositeNode(endpoint,iei);
			Edge lastEdge=iei;

			toTransform.push_back(myProject(skel,lastNode,depthtmp));
			depths.push_back(depthtmp);
			nodelist.push_back(lastNode);
			svec.push_back(svec[0]+(toTransform.back()-toTransform[0]).Norm());

			while(lemon::countIncEdges(g,lastNode)==2)
			{
				for(IncEdgeIt iei2=IncEdgeIt(g,lastNode);iei2!=lemon::INVALID;iei2++)
				{
					if(lastEdge!=iei2)
					{
						lastEdge=iei2;
						lastNode=g.oppositeNode(lastNode,iei2);

						toTransform.push_back(myProject(skel,lastNode,depthtmp));
						nodelist.push_back(lastNode);
						depths.push_back(depthtmp);
						svec.push_back(svec.back()+(toTransform.back()-toTransform[toTransform.size()-2]).Norm());

						break;
					}
				}
			}

			if(minS<0||svec.back()<minS)minS=svec.back();
		}

		for(int i=transformed.size();i<toTransform.size();i++)
		{

		}

		for(int i=0;i<toTransform.size();i++)
			nodeMap[nodelist[i]]=myUnProject(transformed[i].X(),transformed[i].Y(),depths[i]);

		skel.convertGraphToVerticeList();
		skel.convertVeticesListToBranchPoints();
		skel.smoothEachBranch();
	}
	
public:

	vector<Point2f> mouseTrack;

	/*void cutBranch(skeleton_mul &skel,vector<Point2f> track)
	{
	int branch_no=-1,node_no;
	double depth;
	for(int i=0;i<skel.brchPts.size();i++)
	{
	for(int j=0;j<skel.brchPts[i].size()-1;j++)
	{
	Point3f p1=skel.brchPts[i][j],p2=skel.brchPts[i][j+1];

	GLint viewport[4];   
	GLdouble modelview[16];    
	GLdouble projection[16];   
	GLdouble winX[2], winY[2], winZ[2];   
	glGetIntegerv(GL_VIEWPORT, viewport);   
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);    
	glGetDoublev(GL_PROJECTION_MATRIX, projection);    
	gluProject(p1.X(),p1.Y(),p1.Z(),modelview,projection,viewport,winX,winY,winZ);
	gluProject(p2.X(),p2.Y(),p2.Z(),modelview,projection,viewport,winX+1,winY+1,winZ+1);

	winY[0]=viewport[3]-winY[0];
	winY[1]=viewport[3]-winY[1];

	for(int k=0;k<track.size()-1;k++)
	{
	Point2f res;
	int res_status;
	res_status=intersect(Point2f(winX[0],winY[0]),Point2f(winX[1],winY[1]),track[k],track[k+1],res);
	if(res_status&&(branch_no==-1||depth>(res.X()-winX[0])/(winX[1]-winX[0])*(winZ[1]-winZ[0])+winZ[0]))
	{
	depth=(res.X()-winX[0])/(winX[1]-winX[0])*(winZ[1]-winZ[0])+winZ[0];
	branch_no=i;
	node_no=j;
	break;
	}
	}
	}
	}
	if(branch_no!=-1)
	{
	if(node_no==0)
	{
	skel.brchPts[branch_no].erase(skel.brchPts[branch_no].begin());
	skel.convertBranchPointsToVeticesList();
	skel.convertVerticesListToGraph();
	return;
	}
	else if(node_no==skel.brchPts[branch_no].size()-1)
	{
	skel.brchPts[branch_no].pop_back();
	skel.convertBranchPointsToVeticesList();
	skel.convertVerticesListToGraph();
	return;
	}
	else
	{
	vector<Point3f> new_branch;
	for(int i=node_no+1;i<skel.brchPts[branch_no].size();i++)
	new_branch.push_back(skel.brchPts[branch_no][i]);

	for(int i=0;i<new_branch.size();i++)skel.brchPts[branch_no].pop_back();
	skel.brchPts.push_back(new_branch);
	skel.convertBranchPointsToVeticesList();
	skel.convertVerticesListToGraph();
	return ;
	}
	}
	return ;
	}*/

	void undo(skeleton_mul &skel)
	{
		if(!undoGraph.empty())
		{
			redoGraph.push_back(skel.skeletonGraph);
			redoNodeMap.push_back(skel.graphNodeP3f);

			skel.skeletonGraph=undoGraph.back();
			undoGraph.pop_back();
			skel.graphNodeP3f=undoNodeMap.back();
			undoNodeMap.pop_back();

			skel.convertGraphToVerticeList();
			skel.convertVeticesListToBranchPoints();
			skel.smoothEachBranch();
		}
	}

	void redo(skeleton_mul &skel)
	{
		if(!redoGraph.empty())
		{
			undoGraph.push_back(skel.skeletonGraph);
			undoNodeMap.push_back(skel.graphNodeP3f);

			skel.skeletonGraph=redoGraph.back();
			redoGraph.pop_back();
			skel.graphNodeP3f=redoNodeMap.back();
			redoNodeMap.pop_back();

			skel.convertGraphToVerticeList();
			skel.convertVeticesListToBranchPoints();
			skel.smoothEachBranch();
		}
	}
	void doOp_Left(skeleton_mul &skel)
	{
		printf("   Do my job!\n");

		switch(moveTend_Left(skel,mouseTrack))
		{
		case(tendToCutCrossPoint):
			backup(skel);cutCrossPoint(skel,mouseTrack);break;
		case(tendToCutSegment):
			backup(skel);cutBranch(skel,mouseTrack);break;
		case(tendToAttach):
			backup(skel);attachBranches(skel,mouseTrack);break;
		case(tendToJointPoints):
			backup(skel);jointPoints(skel,mouseTrack);break;
		case(tendToDelBranch):
			backup(skel);deleteBranch(skel,mouseTrack);break;
		case(tendToExtend):
			backup(skel);extendbranch(skel,mouseTrack);break;
		case(tendToDefault):
			printf("       Default!\n");
			break;
		}
	}
	void doOp_Right(skeleton_mul &skel)
	{
		Node ntmp;
		Edge etmp;
		if(mouseTrack.size()==2)
		{
			if((ntmp=isOnEndPoint(mouseTrack[0],skel))!=lemon::INVALID)
			{
				printf("       Deform!\n");
				backup(skel);
				deform(skel,ntmp);
			}
			else if((etmp=isOnBranch(mouseTrack[0],skel))!=lemon::INVALID)
			{
				printf("       Deform!\n");
				backup(skel);
				ntmp=getPointOnEdge(skel,etmp,mouseTrack[0]);
				deform(skel,ntmp);
			}
		}
		return;
	}
	void saveSkeletonAsSkel(skeleton_mul skel,QString filename)
	{
		ofstream outfile;
		outfile.open(filename.toStdString().c_str());

		ostringstream strStream; 

		strStream << "CN " << skel.smoothedBrchPts.size() << endl;
		for (int i = 0; i < skel.smoothedBrchPts.size(); i++)
		{
			vector<Point3f>& branch = skel.smoothedBrchPts[i];
			strStream << "CNN " << branch.size() << endl;
			for (int j = 0; j < branch.size(); j++) {
				strStream <<skel.smoothedBrchPts[i][j].X() << "	" << skel.smoothedBrchPts[i][j].Y() << "	" << skel.smoothedBrchPts[i][j].Z() << "	" << endl;
			}
		}
		strStream << endl;

		outfile.write( strStream.str().c_str(), strStream.str().size() ); 
		outfile.close();
	}

	//void saveSkeletonAsSkel(skeleton_mul skel,QString filename)
	//{
	//	ofstream outfile;
	//	outfile.open(filename.toStdString().c_str());

	//	ostringstream strStream; 

	//	skel.convertGraphToVerticeList() ;
	//	skel.convertVeticesListToBranchPoints() ;
	//	strStream << "CN " << skel.branches.size() << endl;
	//	for (int i = 0; i < skel.branches.size(); i++)
	//	{
	//		vector<int>& branch = skel.branches[i];
	//		strStream << "CNN " << branch.size() << endl;
	//		for (int j = 0; j < branch.size(); j++)
	//		{
	//			strStream << skel.nodes[branch[j]].X() << "	" << skel.nodes[branch[j]].Y() << "	" << skel.nodes[branch[j]].Z() << "	" << endl;
	//		}
	//	}
	//	strStream << endl;

	//	outfile.write( strStream.str().c_str(), strStream.str().size() ); 
	//	outfile.close();
	//}

	void init()
	{
		while(!redoGraph.empty())
		{
			lemon::ListGraph *gi=redoGraph.back();
			redoGraph.pop_back();
			P3fNodeMap *mapi=redoNodeMap.back();
			redoNodeMap.pop_back();
			gi->clear();
			delete gi;
			delete mapi;
		}
		while(!undoGraph.empty())
		{
			lemon::ListGraph *gi=undoGraph.back();
			undoGraph.pop_back();
			P3fNodeMap *mapi=undoNodeMap.back();
			undoNodeMap.pop_back();
			gi->clear();
			delete gi;
			delete mapi;
		}
	}
};

#endif