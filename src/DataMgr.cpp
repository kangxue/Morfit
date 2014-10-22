
#include "DataMgr.h"
#include <QMessageBox>

#include <fstream>
extern std::ofstream logout ;
#include "appstate.h"
extern  AppState appstate ;

DataMgr::DataMgr(RichParameterSet* _para)
{
	para = _para;
	//skeleton = new Skeleton ;
}


DataMgr::~DataMgr(void)
{
	//delete skeleton ;
}

void DataMgr::clearCMesh(CMesh& mesh)
{

	mesh.face.clear();
	mesh.fn = 0;
	mesh.vert.clear();
	mesh.vn = 0;
	mesh.bbox = Box3f();

	mesh.Clear() ;

}

bool DataMgr::isSamplesEmpty()
{
	return samples.vert.empty();
}


bool DataMgr::isOriginalEmpty()
{
	return original.vert.empty();
}

void DataMgr::loadPlyToOriginal(QString fileName)
{
	//load the file of ply and buffer it in CMesh original
	std:: cout << "begin DataMgr::loadPlyToOriginal\n";

	//clear the variable original
	clearCMesh(original);
	curr_file_name = fileName;

	int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;

	int err = tri::io::Importer<CMesh>::Open( original, curr_file_name.toStdString().c_str(), mask);  
	if(err) 
	{
		QMessageBox msgBox;
		msgBox.setText(QString("fail to load ") + fileName+"\nNote that this program doesn't support unicode file path.");
		msgBox.exec();

		return;
	}  


	CMesh::VertexIterator vi;
	int idx = 0;
	for(vi = original.vert.begin(); vi != original.vert.end(); ++vi)
	{
		vi->bIsOriginal = true;
		vi->m_index = idx++;
		original.bbox.Add(vi->P());
	}
	original.vn = original.vert.size();

	downSamplesByNum();
	
	std:: cout << "end" <<std::endl;
}

void DataMgr::loadPlyToPossion(QString fileName)
{
	std:: cout << "begin DataMgr::loadPlyToPossion\n";

	clearCMesh(possion);
	curr_file_name = fileName;

	int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;

	int err = tri::io::Importer<CMesh>::Open(possion, curr_file_name.toStdString().c_str(), mask);  
	if(err) 
	{
		cout << "Failed reading mesh: " << err << "\n";
		return;
	}  


	CMesh::VertexIterator vi;
	int idx = 0;
	for(vi = possion.vert.begin(); vi != possion.vert.end(); ++vi)
	{
		vi->bIsOriginal = true;
		vi->m_index = idx++;
		possion.bbox.Add(vi->P());
	}
	possion.vn = possion.vert.size();


	std:: cout << "end" <<std::endl;
}


void DataMgr::loadPlyToSample(QString fileName)
{
	std:: cout << "begin DataMgr::loadPlyToSam\n";

	clearCMesh(samples);
	curr_file_name = fileName;

	int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;
	mask += tri::io::Mask::IOM_VERTCOLOR;
	mask += tri::io::Mask::IOM_BITPOLYGONAL;

	int err = tri::io::Importer<CMesh>::Open(samples, curr_file_name.toStdString().c_str(), mask);  
	if(err) 
	{
		cout << "Failed reading mesh: " << err << "\n";
		return;
	}  

	CMesh::VertexIterator vi;
	int idx = 0;
	for(vi = samples.vert.begin(); vi != samples.vert.end(); ++vi)
	{
		vi->bIsOriginal = false;
		vi->m_index = idx++;
		samples.bbox.Add(vi->P());
	}
	samples.vn = samples.vert.size();


	std:: cout << "end" <<std::endl;
}

#include "reconstructorPara.h"
void DataMgr::loadSkeletonWithPly(QString fileName ){

	//load the skeleton file of skel and buffer it in skeleton_mul skel
	std::cout << "begin DataMgr::loadSkeletonWithPly\n";

	std::string fname = fileName.toStdString() +".skel" ;

	ifstream infile;
	infile.open( fname.c_str());

	if( infile.fail() ){
		QMessageBox msgBox;
		msgBox.setText(QString("fail to load ") + fileName + ".skel" );
		msgBox.exec();


		std::cout<<"end" <<std::endl;
		return ;
	}
		
	stringstream sem; 
	sem << infile.rdbuf(); 

	string str;
	int num;
	int num2;


	while( sem >> str ){

		if (str == "CN")
		{

			// get number of branches
			sem >> num;

			std::vector< std::vector<Point3f> > branches(num) ; 

			for (int i = 0; i < num; i++)
			{
				sem >> str;
				sem >> num2;
				for(int j = 0; j < num2; j++)
				{
					Point3f p;
					sem >> p[0] >> p[1] >> p[2];
					branches[i].push_back(p) ;
				}

			}

			// delete too short branches
			for( int i=0; i< branches.size(); ++i )
				if( branches[i].size() <= ReconstructorPara::minimumBranchPointNum  )
				{	branches.erase( branches.begin()+i) ; i--; }

			// convert branches to skeleton_mul
			skel = skeleton_mul(branches, true) ;

			//save the branches also
			theBranches=branches;
			break ;
		}

	}



	infile.close() ;

	

}

void DataMgr::intializeSwp(){

	// initialize sweeper
	std::vector<Point3f> pointset ;
	for( int i=0; i<original.vert.size(); ++i )
		pointset.push_back(original.vert[i].P() ) ;
	std::vector<Point3f> normals ;
	for( int i=0; i<original.vert.size(); ++i )
		normals.push_back(original.vert[i].N() ) ;

	//swp = sweeper(pointset,normals, skel.smoothedBrchPts ) ;
    

	//swp = sweeper(pointset,normals, skel.smoothedBrchPts) ;
	//swp = sweeper(pointset,normals,skel,allJoints);
	swp = sweeper(pointset,normals,skel.smoothedBrchPts);
}

CMesh* DataMgr::getCurrentSamples()
{
	if(samples.vert.empty())
	{
		//cout << "DataMgr::getCurrentSamples samples.vert.empty()!!" <<endl;
		return NULL;
	}
	if(&samples == NULL)
	{
		//cout << "DataMgr::getCurrentSamples samples = NULL!!" <<endl;
		return NULL;
	}
	return & samples;
}

CMesh* DataMgr::getCurrentOriginal()
{
	if(original.vert.empty())
	{
		//cout << "DataMgr::getCurrentOriginal() original.vert.empty()!!" <<endl;
		return NULL;
	}
	if(&original == NULL)
	{
		//cout << "DataMgr::getCurrentOriginal() samples = NULL!!" <<endl;
		return NULL;
	}
	return & original;
}



void DataMgr::recomputeBox()
{
	//compute the box of the model
	std:: cout << "begin DataMgr::recomputeBox\n";

	//if( appstate.plyfilename.find( "recSkirt") != string::npos ){
	//	return ;
	//}

	samples.bbox.SetNull();
	original.bbox.SetNull();

	CMesh::VertexIterator vi;
	for(vi = samples.vert.begin(); vi != samples.vert.end(); ++vi) 
	{
		if (vi->is_skel_ignore)
		{
			continue;
		}
		samples.bbox.Add(vi->P());
	}

	for(vi = original.vert.begin(); vi != original.vert.end(); ++vi) 
	{
		original.bbox.Add(vi->P());
	}

	std::cout << original.bbox.min.X() << " "<< original.bbox.min.Y() << " "<< original.bbox.min.Z() << std::endl;
	std::cout << original.bbox.max.X() << " "<< original.bbox.max.Y() << " "<< original.bbox.max.Z() << std::endl;

	std::cout << samples.bbox.min.X() << " "<< samples.bbox.min.Y() << " "<< samples.bbox.min.Z() << std::endl;
	std::cout << samples.bbox.max.X() << " "<< samples.bbox.max.Y() << " "<< samples.bbox.max.Z() << std::endl;

	std:: cout << "end" <<std::endl;
}

double DataMgr::getInitRadiuse()
{
	//compute the radius of the model and return the result of double
	std:: cout << "begin DataMgr::getInitRadiuse\n";

	double init_para = para->getDouble("Init Radius Para");
	if (!isOriginalEmpty())
	{
		//use the data bbox of original,should call after recomputeBox() 
		Box3f box = original.bbox;
		if ( abs(box.min.X() - box.max.X()) < 1e-5 ||   
			abs(box.min.Y() - box.max.Y()) < 1e-5 ||   
			abs(box.min.Z() - box.max.Z()) < 1e-5 )
		{
			double diagonal_length = sqrt((box.min - box.max).SquaredNorm());
			double original_size = sqrt(double(original.vn));
			init_radius = 2 * init_para * diagonal_length / original_size;
		}
		else
		{
			double diagonal_length = sqrt((box.min - box.max).SquaredNorm());
			double original_size = pow(double(original.vn), 0.333);
			init_radius = init_para * diagonal_length / original_size;
		}
	}

	std:: cout << "end" <<std::endl;
	return init_radius;
}


void DataMgr::downSamplesByNum()
{
	//vcg::tri::Append::Mesh( samples,original,false) ;

	//just CopyCmesh(samples, original), get the value of original to samples
	ReconstructorUtility::CopyCmesh(samples, original) ;

	return ;

	//****the following not execute actually*****
	std:: cout << "begin DataMgr::downSamplesByNum\n";

	if (isOriginalEmpty() && !isSamplesEmpty())
	{
		subSamples();
		return;
	}

	if (isOriginalEmpty())
	{
		return;
	}

	int want_sample_num = para->getDouble("Down Sample Num");

	want_sample_num = original.vert.size() * 0.1;
	//if( want_sample_num > 10000 )
		want_sample_num = 10000 ;

	if (want_sample_num > original.vn)
	{
		want_sample_num = original.vn;
	}

	clearCMesh(samples);
	samples.vn = want_sample_num;

	vector<int> nCard = GlobalFun::GetRandomCards(original.vert.size());
	for(int i = 0; i < samples.vn; i++) 
	{
		int index = nCard[i]; //1-18 not random!

		CVertex& v = original.vert[index];
		samples.vert.push_back(v);
		samples.bbox.Add(v.P());
	}

	CMesh::VertexIterator vi;
	for(vi = samples.vert.begin(); vi != samples.vert.end(); ++vi)
	{
		vi->bIsOriginal = false;
	}

	std:: cout << "end" <<std::endl;
}

void DataMgr::subSamples()
{
	std:: cout << "begin DataMgr::subSamples\n";

	clearCMesh(original);

	CMesh::VertexIterator vi;
	original.vn = samples.vert.size();
	original.bbox.SetNull();
	for(vi = samples.vert.begin(); vi != samples.vert.end(); ++vi)
	{
		CVertex v = (*vi);
		v.bIsOriginal = true;
		original.vert.push_back(v);
		original.bbox.Add(v.P());
	}

	downSamplesByNum();
	std:: cout << "end" <<std::endl;
}


void DataMgr::savePly(QString fileName, CMesh& mesh)
{
	int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;
	mask += tri::io::Mask::IOM_VERTCOLOR;
	mask += tri::io::Mask::IOM_BITPOLYGONAL;


	std::string fname = fileName.toStdString() ;
	if( fname.find(".ply") == string::npos )
		fname = fname + ".ply" ;

	std::cout << "saveFile: " << fname <<std::endl;

	if( tri::io::ExporterPLY<CMesh>::Save(mesh, fname.c_str(), mask, false) != 0);
		std::cout << "save fail" <<std::endl;
}

void DataMgr::normalizeROSA_Mesh(CMesh& mesh)
{
	//this function can normalize the mesh 
	Box3f box_bk ;
	//if( appstate.plyfilename.find( "recSkirt") != string::npos ){
	//	box_bk = mesh.bbox ;
	//}


	std:: cout << "begin DataMgr::normalizeROSA_Mesh\n";


	if (mesh.vert.empty())
	{
		return;
	}
	Box3f box = mesh.bbox;

	std::cout << "bbox = "<<box.min.X() << " " <<box.min.Y() << " "<<box.min.Z() << "\n"
						  <<box.max.X() << " "<<box.max.Y() << " "<<box.max.Z() << std::endl;


	mesh.bbox.SetNull();
	float max_x = abs((box.min - box.max).X());
	float max_y = abs((box.min - box.max).Y());
	float max_z = abs((box.min - box.max).Z());
	float max_length = max_x > max_y ? max_x : max_y;
	max_length = max_length > max_z ? max_length : max_z;

	for(int i = 0; i < mesh.vert.size(); i++)
	{
		Point3f& p = mesh.vert[i].P();

		p -= box.min;
		p /= max_length;

		p = (p - Point3f( 0.5*max_x/max_length,  0.5*max_y/max_length,  0.5*max_z/max_length));

		mesh.vert[i].N().Normalize(); 
		mesh.bbox.Add(p);
	}


	//if( appstate.plyfilename.find( "recSkirt") != string::npos ){
	//	mesh.bbox = box_bk;
	//}

	std:: cout << "end" <<std::endl;
}

void DataMgr::normalizeSkeletonkx( skeleton_mul &skel, Box3f box ) {


	std:: cout << "begin DataMgr::normalizeSkeletonkx\n";

	skel.normalize( box ) ;


	std:: cout << "end" <<std::endl;

}

Box3f DataMgr::normalizeAllMesh()
{


	std:: cout << "begin DataMgr::normalizeAllMesh\n";
	Box3f box;
	if (!isSamplesEmpty())
	{
		for (int i = 0; i < samples.vert.size(); i++)
		{
			box.Add(samples.vert[i].P());
		}
		samples.bbox = box;
	}
	if (!isOriginalEmpty())
	{
		for (int i = 0; i < original.vert.size(); i++)
		{
			box.Add(original.vert[i].P());
		}
		original.bbox =box;
	}

	/*if( appstate.plyfilename.find( "recSkirt" ) != string::npos  ){
		

		Box3f box_x;
		
		std::string fname = appstate.path + "\/recSkirt2_old.ply" ;
		
		int mask= tri::io::Mask::IOM_VERTCOORD + tri::io::Mask::IOM_VERTNORMAL ;

		int err = tri::io::Importer<CMesh>::Open( samples, fname.c_str(), mask);  
		if(err) 
		{
			cout << "Failed reading mesh: err = " << err << "\n" << "fname = " << fname <<std::endl;;

			system("pause") ;
		}  

		CMesh::VertexIterator vi;
		for(vi = samples.vert.begin(); vi != samples.vert.end(); ++vi){
			box_x.Add(vi->P());
		}
		samples.vn = samples.vert.size();



		box = original.bbox = samples.bbox = box_x;

	}*/

	//normalize the mesh samples and original use normalizeROSA_Mesh()
	normalizeROSA_Mesh(samples);
	normalizeROSA_Mesh(original);

	boxBeforeNormalize = box ;
	normalizeSkeletonkx( skel, box ) ;

	//compute the box and radius
	recomputeBox();
	getInitRadiuse();

	std:: cout << "end" <<std::endl;

	return samples.bbox;
}


void DataMgr::eraseRemovedSamples()
{
	std:: cout << "begin DataMgr::eraseRemovedSamples\n";
	int cnt = 0;
	vector<CVertex> temp_mesh;
	for (int i = 0; i < samples.vert.size(); i++)
	{
		CVertex& v = samples.vert[i];
		if (!v.is_skel_ignore)
		{
			temp_mesh.push_back(v);
		}
	}

	samples.vert.clear();
	samples.vn = temp_mesh.size();
	for (int i = 0; i < temp_mesh.size(); i++)
	{
		temp_mesh[i].m_index = i;
		samples.vert.push_back(temp_mesh[i]);
	}

	std:: cout << "end" <<std::endl;
}

void DataMgr::clearData()
{
	std:: cout << "begin DataMgr::clearData\n";
	clearCMesh(original);
	clearCMesh(samples);

	curves2d.clear() ;
	curves3d.clear() ;
	snappedCurves3d.clear() ;
	sampleNormDiff.clear() ;

	std:: cout << "end" <<std::endl;
}


void DataMgr::updateSampEigenVectors() {
	//for  vector of each sample "recompute_m_render"
	for( int i=0; i<samples.vert.size();++i )
		samples.vert[i].recompute_m_render() ;
}