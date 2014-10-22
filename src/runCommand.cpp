
#include <QMessageBox>
#include "GLArea.h"
#include "appstate.h"
extern  AppState appstate ;





//#include <QTextStream>
//
//#include <stdio.h>
//

//QTextStream cin(stdin, QIODevice::ReadOnly);
//QTextStream cout(stdout, QIODevice::WriteOnly);
//QTextStream cerr(stderr, QIODevice::WriteOnly);


//if( cin.status() != QTextStream::Ok ){
//	std::cout <<"cin fail"<<std::endl;
//	system("pause") ;
//}

triangleList geckoBoundary ;  //  for labeling interaction area


void GLArea::runCommand( std::string cmd ){

	std::cout << "run command: "<<cmd<<std::endl;
	
	char buffer[128] ;

	if( cmd == "save og nbs" ) {

		if( get_operating_state() != onesweepReady )
			return ;

		std::ofstream ofs( appstate.path + "\/ognbs.txt" ) ;
			
		std::vector<int> tmpIdx =  dataMgr.swp.ongoingSkel.tempalteIdss[0] ;
		std::vector<ON_NurbsCurve> nbss =  dataMgr.swp.ongoingSkel.Nurbsss[0] ;
		Curve3D  brch = dataMgr.swp.ongoingSkel.smoothedBrchPts[0] ;
		std::vector<double>  conf = dataMgr.swp.ongoingSkel.profile_conf[0] ;

		if( tmpIdx.size() != nbss.size() ){
			std::cout <<__FILE__<<__LINE__<<std::endl;
			system("error") ;
		}

		// save skeleton
		ofs << brch.size() <<std::endl ;
		for( int i=0; i<brch.size(); ++i )
			ofs << brch[i].X() << " " <<  brch[i].Y()<< " " <<  brch[i].Z() <<std::endl;

		// save confidence
		ofs << conf.size() <<std::endl ;
		for( int i=0; i<conf.size(); ++i )
			ofs << conf[i] <<std::endl;


		// save tmpId 
		ofs<< tmpIdx.size() <<std::endl;
		for( int i=0; i<tmpIdx.size(); ++i )
			ofs << tmpIdx[i] <<std::endl;

		// save nurbs
		ofs<< nbss.size() <<std::endl;
		for( int i=0; i<nbss.size(); ++i ){
			Curve2D cvs = ReconstructorUtility::getCvPts( nbss[i] ) ;
			for( int j=0; j<cvs.size(); ++j )
				ofs << cvs[j].X() << " " <<  cvs[j].Y() <<std::endl;
		}

		ofs.close() ;

		int ogSkelId = dataMgr.swp.settings.size() ;
		if( ogSkelId!=0 && dataMgr.swp.settings.back().iditifier == dataMgr.swp.ongoingSkel.iditifier )
			ogSkelId-- ;

		std::string cpyCmd = std::string( "copy " ) + "\"" + appstate.path + "\/ognbs.txt\"" + " \"" + appstate.path + "\/ognbs_" + itoa( ogSkelId, buffer,10 ) + ".txt\"" ;
		
		for( int i=0; i<cpyCmd.size(); ++i )
			if( cpyCmd[i] == '\/' ) 
				 cpyCmd[i] = '\\' ;

		system(cpyCmd.c_str()) ;

		std::cout << "cpyCmd = "<< cpyCmd <<std::endl;

	} else if( cmd == "load og nbs" ){

		if( get_operating_state() != onesweepReady )
			return ;
		
		int ogSkelId = dataMgr.swp.settings.size() ;
		if( ogSkelId!=0 && dataMgr.swp.settings.back().iditifier == dataMgr.swp.ongoingSkel.iditifier )
			ogSkelId-- ;


		std::string fname1 = appstate.path + "\/ognbs_" + itoa( ogSkelId, buffer,10 ) + ".txt" ;
		std::string fname2 = appstate.path + "\/ognbs.txt" ;
		
		std::ifstream ifs(fname1) ;

		if( ifs.fail() ){
			std::cout << "didn't find " << fname1 <<std::endl;
			ifs.clear() ;
			ifs.open( fname2 ) ;

			if( ifs.fail() ){
				std::cout << "didn't find " << fname2 <<std::endl;
				return ;
			}

		}


		std::vector<int> tmpIdx ;
		std::vector<ON_NurbsCurve> nbss;
		Curve3D  brch ;
		std::vector<double>  conf ;

		// load skeleton
		int brchSize ;
		ifs >> brchSize ;
		for( int i=0; i<brchSize; ++i ){
			double x,y,z ;  
			ifs >> x >> y >> z;  
			brch.push_back( Point3f(x,y,z) ) ;
		}


		// load confidence
		int confSize ;
		ifs >> confSize ;
		for( int i=0; i<confSize; ++i ){
			double t;  ifs >>t ;   conf.push_back( t) ;
		}


		// load tmpId 
		int tmpNum ;
		ifs >> tmpNum ;
		for( int i=0; i<tmpNum; ++i ){
			int id;
			ifs >> id ;
			tmpIdx.push_back( id ) ;
		}


		// load nurbs

		int nbsNum ;
		ifs >> nbsNum ;

		for( int i=0; i<nbsNum; ++i ){
			Curve2D cvs  ;
				
			for( int j=0; j< ReconstructorPara::cvNum; ++j ){
				Point2f p ;
				ifs >> p.X() >> p.Y() ;
				cvs.push_back(p) ;
			}

			ON_NurbsCurve nbs = ReconstructorUtility::fitNURBSforProfile2d( dataMgr.swp.ongoingSkel.profiles2d[0][0] ) ;
			ReconstructorUtility::setCvOfNurbs(nbs, cvs) ;
			nbss.push_back( nbs ) ;
		}

		ifs.close() ;


		
		if( tmpIdx.size() != nbss.size() ){
			std::cout <<__FILE__<<__LINE__<<std::endl;
			system("error") ;
		}

		// rank tmpId and nbs
		std::vector<int> newTmpId = tmpIdx ;
		std::sort( newTmpId.begin(), newTmpId.end() ) ;

		std::vector<ON_NurbsCurve> newNbss;
		for( int i=0; i<newTmpId.size(); ++i )
			newNbss.push_back(  nbss[ReconstructorUtility::FindInVerctor(tmpIdx, newTmpId[i] ) ] ) ;
		
		tmpIdx = newTmpId ;
		nbss = newNbss ;




		dataMgr.swp.ongoingSkel.smoothedBrchPts = std::vector<Curve3D>(1, brch) ;
		std::vector<int> subptsMapId;
		subptsMapId.clear();
		dataMgr.swp.ongoingSkel.intialSegment( dataMgr.swp.ongoingSkel.points ) ;
		dataMgr.swp.ongoingSkel.calculateCutplanes() ;
		dataMgr.swp.ongoingSkel.calculateProfiles() ;
		dataMgr.swp.ongoingSkel.initAfterClaculateProfiles() ;

		dataMgr.swp.ongoingSkel.profile_conf[0] = conf ;
		dataMgr.swp.ongoingSkel.tempalteIdss[0] = tmpIdx;
		dataMgr.swp.ongoingSkel.Nurbsss[0] = nbss;
		dataMgr.swp.ongoingSkel.discretizeNurbs() ;

		//dataMgr.swp.ongoingSkel.firstMustBeTip = false ;
		//dataMgr.swp.ongoingSkel.firstMustBeNotTip = true ;
		//dataMgr.swp.ongoingSkel.lastMustBeTip = false ;
		//dataMgr.swp.ongoingSkel.lastMustBeNotTip = true ;


	}else if( cmd == "load mesh"){

		CMesh m ;
		int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;
		int err = tri::io::Importer<CMesh>::Open( m, (appstate.path+"\/mesh.ply").c_str(), mask);  
		if(err) {
			cout << "Failed reading mesh: " << err << "\n";


		QMessageBox msgBox;
		msgBox.setText( QString("fail to load ") + QString(  (appstate.path+"\/mesh.ply").c_str() ) );
		msgBox.exec();


			return;
		}  

		m.bbox = dataMgr.boxBeforeNormalize ;
		
		dataMgr.normalizeROSA_Mesh(m) ;

		dataMgr.swp.mergedMesh = triangleList(m) ;
		dataMgr.swp.mergedMesh.cvtToCmesh() ;

	
	}else if( cmd == "save mesh"){

		CMesh m = dataMgr.swp.mergedMesh.cvtToCmesh() ;
		int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;
		mask += tri::io::Mask::IOM_VERTCOLOR;
		mask += tri::io::Mask::IOM_BITPOLYGONAL;
		tri::io::ExporterPLY<CMesh>::Save(m,  (appstate.path+"\/mesh.ply").c_str(), mask, false);

	}else if(cmd == "save brch mesh"){
	
		for( int i=0; i<dataMgr.swp.meshes.size();++i){
			CMesh m = dataMgr.swp.meshes[i].cvtToCmesh() ;
			int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;
			mask += tri::io::Mask::IOM_VERTCOLOR;
			mask += tri::io::Mask::IOM_BITPOLYGONAL;
			char buffer[1024] ;
			tri::io::ExporterPLY<CMesh>::Save(m,  (appstate.path+"\/mesh_brch"+std::string(itoa(i,buffer,10))+".ply").c_str(), mask, false);
		}

	}else if(cmd == "save settings"){
		dataMgr.swp.saveSettings() ;
	}else if(cmd == "load settings"){
		dataMgr.swp.loadSettings() ;
	}else if( cmd == "recenter"){

		dataMgr.swp.recenterizeSkel( dataMgr.skel ) ;

	}else if(cmd=="h"){
		//std::cout << "command list:"<<std::endl; 
		//std::cout << "save  og nbs"<<std::endl; 
		//std::cout << "load  og nbs"<<std::endl; 

	} else {
		std::cout << "command: " << cmd << " not found " <<std::endl;
	}
	
	
	
	/*else {

		cin.ignore(numeric_limits <streamsize> ::max(), '\n');
		std::cout << " no operation matched, input again:"<<std::endl;

		getline(cin, cmd) ;

		if( cmd.size()==1 && cmd[0]=='\n' )
			return;
	}*/


}