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

#include <QDir>
#include <QMessageBox>

#include "skeleton_mul.h"
#include "sweeper.h"

#include "appstate.h"
extern  AppState appstate ;

#include "reconstructionUtility.h"
void sweeper::saveSettings() {

	std::cout << "sweeper::saveSettings() {"<<std::endl;

	char buffer[1024] ;


	std::string path =  appstate.path + "\/settings" ;

	QDir dir ;
	dir.mkpath ( QString(path.c_str()) );


	std::string settingNum_fname = path + "\/settingNum.txt" ;
	// save setting number
	std::ofstream ofs( settingNum_fname ) ;

	if( ofs.fail() ){
		std::cout << "fail to open "<<settingNum_fname<<std::endl;
		return ;
		//system("pause") ;
	}

	ofs << settings.size() ;
	ofs.close() ;

	// save meshes
	assert( settings.size() == meshes.size() ) ;
	for( int i=0; i<settings.size(); ++i ){

		std::string submesh_fname = path+"\/submesh_"+std::string(itoa(i,buffer,10))+".ply" ;
		meshes[i].cvtToCmesh() ;
		int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;
		mask += tri::io::Mask::IOM_VERTCOLOR;
		mask += tri::io::Mask::IOM_BITPOLYGONAL;
		tri::io::ExporterPLY<CMesh>::Save(meshes[i].cmeshFormat,  submesh_fname.c_str(), mask, false);

	}

	// save  solvers
	for( int i=0; i<settings.size(); ++i ){
		std::string brch_fname =  path+"\/solver_"+std::string(itoa(i,buffer,10)) + ".txt" ;

		ofs.open( brch_fname ) ;

		ofs<<settings[i].solvers[0].cvtToString() ;
		
		ofs.close() ;
	}

	// save profiles3d
	for( int i=0; i<settings.size(); ++i ){
		std::string p3d_fname =  path+"\/profile3d_"+std::string(itoa(i,buffer,10)) + ".txt" ;

		ofs.open( p3d_fname ) ;

		std::vector<Curve3D> prof3d = settings[i].profiles3d[0] ;

		ofs << prof3d.size() <<std::endl;
		for( int j=0; j<prof3d.size(); ++j ){
			ofs << prof3d[j].size() <<std::endl ;
			for( int k=0; k<prof3d[j].size(); ++k  )
				ofs << prof3d[j][k].X() << " " <<prof3d[j][k].Y() << " " <<prof3d[j][k].Z() << std::endl;
		}

		ofs.close() ;
	}



	std::cout << "}"<<std::endl;

}


void sweeper::loadSettings() {

	std::cout << "sweeper::loadSettings() {"<<std::endl;

	char buffer[1024] ;
	std::string path =  appstate.path + "\/settings" ;
	std::string settingNum_fname = path + "\/settingNum.txt" ;


	int m,n , u,v,t;
	double x,y,z ;
	// load setting number
	int settingNum = 0 ;
	std::ifstream ifs(settingNum_fname) ;


	if( ifs.fail() ){
		std::cout << "fail to open "<<settingNum_fname<<std::endl;
		
		QMessageBox msgBox;
		msgBox.setText("setting file does not exist!");
		msgBox.exec();

		return ;
		//system("pause") ;
	}


	if( !ifs.fail() ){
		ifs >> settingNum ;
	}
	ifs.close() ;


	// load meshes
	meshes.clear() ;
	for( int i=0; i<settingNum; ++i ){
		
		std::string submesh_fname = path+"\/submesh_"+std::string(itoa(i,buffer,10))+".ply" ;

		CMesh m ;
		int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;
		int err = tri::io::Importer<CMesh>::Open( m, submesh_fname.c_str(), mask);  
		if(err) {
			cout << "Failed reading mesh: " << err << "\n";
			return;
		}  

		meshes.push_back( triangleList(m));

	}

	// load solvers
	settings.clear() ;
	for( int i=0; i<settingNum; ++i ){
		std::string brch_fname =  path+"\/solver_"+std::string(itoa(i,buffer,10)) + ".txt" ;

		ifs.open( brch_fname ) ;
		string str((istreambuf_iterator<char>(ifs)), istreambuf_iterator<char>());
	
		completionSolver slver = completionSolver(str) ;



		settings.push_back( skelpath(  slver )) ;

		ifs.close() ;
	}


	// load profiles3d
	for( int i=0; i<settings.size(); ++i ){
		std::string p3d_fname =  path+"\/profile3d_"+std::string(itoa(i,buffer,10)) + ".txt" ;

		ifs.open( p3d_fname ) ;

		ifs >> n ;
		std::vector<Curve3D> prof3d(n)  ;
		for( int j=0; j<prof3d.size(); ++j ){
			ifs >> m ;
			for( int k=0; k<m; ++k  ){
				ifs >> x>>y>>z ;
				prof3d[j].push_back( Point3f(x,y,z) ) ;
			}
		}

		settings[i].profiles3d = std::vector<std::vector<Curve3D>>(1, prof3d ) ;

		ifs.close() ;
	}


	// write iditifier
	meshLable.clear() ;
	for( int i=0; i<settings.size(); ++i )
		meshLable.push_back( settings[i].iditifier ) ;

	// write meshes
	for( int i=0; i<settingNum; ++i )
		settings[i].resultMesh = std::vector<triangleList>(1, meshes[i]) ;



	mergeMesh() ;







	std::cout << "}"<<std::endl;

}
