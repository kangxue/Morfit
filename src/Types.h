/*
	This file is part of the software described in the paper,

	@ARTICLE{Morfit,
	  title = {Morfit: Interactive Surface Reconstruction from Incomplete Point Clouds with Curve-Driven Topology and Geometry Control},
	  author = {Kangxue Yin and Hui Huang and Hao Zhang and Minglun Gong and Daniel Cohen-Or and Baoquan Chen},
	  journal = {ACM Transactions on Graphics(Proc. of SIGGRAPH Asia 2014)},
	  volume = {33},
	  number = {6},
	  pages = {41:1--41:12},
	  year = {2014},
	}

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
#ifndef _types_h_
#define _types_h_


#include <vector>
 #include <QColor>
#include "cmesh.h"
//#include <vcg\space\deprecated_point3.h>

typedef std::vector<Point3f> Curve3D ; 
typedef std::vector<Curve3D> CurveArray3D ; 
typedef std::vector<Point2f> Curve2D ; 
typedef std::vector<Curve2D> CurveArray2D ; 

typedef std::vector<Point2f> Profile2D ;
typedef std::vector<Point3f> Profile3D ;


class GLColor
{
public:
	GLColor(const float& _r = 0, const float& _g = 0, const float& _b = 0, const float& _a = 1.0):r(_r),g(_g),b(_b),a(_a){}
	GLColor(const QColor & qcolor)
	{
		int _r;
		int _g;
		int _b;
		qcolor.getRgb(&_r, &_g, &_b);
		r = _r / 255.0;
		g = _g / 255.0;
		b = _b / 255.0;
	}
	float r;
	float g;
	float b;
	float a;
};


#define MyMax(a,b) (((a) > (b)) ? (a) : (b))  
#define MyMin(a,b) (((a) < (b)) ? (a) : (b))  


class double2{
public:
	double x ;
	double y ;

	double norm(){
		return sqrt( x*x+y*y) ;
	}
	double2 normalize(){
		double n = norm() ;
		x /= n;
		y /= n ;

		return *this ;
	}
	double2(){ x = y = 0.0;}
	double2(double _x, double _y ): x(_x), y(_y) {}

	double2 operator+(double2 &a){
		return double2(x+a.x, y+a.y) ;
	}
	double operator*(double2 &a){
		return x*a.x + y*a.y ;
	}
	double2 operator*(double a){
		return double2(a*x, a*y) ;
	}
	double2 operator/(double a){
		return double2(x/a, y/a) ;
	}
	double2 operator-(double2 &a){
		return double2(x-a.x, y-a.y) ;
	}

	double2 operator-(){
		return double2(-x, -y) ;
	}

	bool operator==(double2 a){
		return x==a.x&&y==a.y ;
	}
};


class int2{
public:
	int2():first(0),second(0){}
	int2(int a, int b):first(a), second(b) {}

	union{
		int first ;
		int x ;
	} ;

	union{

		int second ;
		int y ;
	};

	bool operator==( int2 t){
		return x==t.x && y==t.y ;
	}

	int &operator[](int a){
		if( a == 0)
			return x ;
		else if( a==1)
			return y;
		else{
			std::cerr<<"invalid index of int2" ;
			exit(1) ;
		}
	}

	inline int &X(){ return x; } 
	inline int &Y(){ return y; } 

};

class int3{
public:
	int3():first(0),second(0),third(0) {}
	int3(int a, int b, int c):first(a), second(b), third(c) {}

	union{  int first ;  int x ; };
	union{  int second;  int y ; };
	union{  int third;   int z ; };
	bool operator==( int3 t){
		return x==t.x && y==t.y &&z==t.z;
	}

	int &operator[](int a){
		if( a == 0)
			return x ;
		else if( a==1)
			return y;
		else if( a==2)
			return z;
		else{
			std::cerr<<"invalid index of int2" ;
			system("pause") ;
		}
	}

	inline int &X(){ return x; } 
	inline int &Y(){ return y; } 
	inline int &Z(){ return z; } 

};

class double3{
public:
	union{
		struct{
			double x, y, z;
		};

		struct{
			double r, g, b;
		};
	};
	double3(){x=y=z=0; }
	double3(double x0, double y0, double z0):x(x0), y(y0), z(z0){}

	double3 operator/(double a){
		return double3(x/a, y/a, z/a) ;
	}
	double3 operator*(double a){
		return double3(x*a, y*a, z*a) ;
	}
	double3 operator-(double3 &a){
		return double3(x-a.x, y-a.y, z-a.z) ;
	}
	double3 operator+(double3 &a){
		return double3(x+a.x, y+a.y, z+a.z) ;
	}
	double norm(){
		return sqrt( x*x+y*y+z*z) ;
	}

	double3 normalize(){
		double n = norm() ;
		x /= n ;
		y /= n ;
		z /= n ;
		return *this ;
	}

};



class PointPL: public Point2f{
public:

	PointPL(double theta, double radius){
		X() = theta ;
		Y() = radius ;
	}


	float &th(){ return X() ; }  // theta
	float &rd(){ return Y() ; }  // radius
};

typedef std::vector<PointPL> ProfilePL ;
//class triangleList{
//public:
//	std::vector<Point3f> vertices ;
//	std::vector<int3> face ;
//}

class uniformValue{
public:
	uniformValue(){}
	uniformValue(QString n, float v):vname(n),val(v){}

	QString vname ;
	float val ;
};

class MVMatrix{
public:
	float m[16] ;
};

enum endType{
	ET_TipEnd,
	ET_FaceEnd,
	ET_MidNode
};

enum imgOutputType
{
	balckdot, 
	withsharpness, 
	withcompleteness 
};

struct cutPlane{
	Point3f cutRadius ;
	Point3f planeNormal ;
};
#endif