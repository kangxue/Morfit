
#pragma once

#include "Types.h"

#include "Parameter.h"
//#include "Algorithm/Skeleton.h"
#include "GlobalFunction.h"
#include "sweeper.h"
#include "skeleton_mul.h"

#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <sstream>
#include <fstream>
#include <set>


using namespace vcg;
using namespace std;
using namespace tri;


class sweeper ;



class  DataMgr
{
public:
	DataMgr(RichParameterSet* _para);
	~DataMgr(void);

	void loadPlyToOriginal(QString fileName);

	void loadSkeletonWithPly(QString fileName ) ;

	void loadPlyToPossion(QString fileName);

	void loadPlyToSample(QString fileName);

	void savePly(QString fileName, CMesh& mesh);

	bool isSamplesEmpty();
	bool isOriginalEmpty();

	CMesh* getCurrentSamples();
	CMesh* getCurrentOriginal();


	void recomputeBox();

	double getInitRadiuse();

	void downSamplesByNum();
	void subSamples();

	void normalizeROSA_Mesh(CMesh& mesh);

	Box3f normalizeAllMesh();
	
	void normalizeSkeletonkx( skeleton_mul &skel, Box3f box ) ;
	
	void eraseRemovedSamples();
	void clearData();

	void updateSampEigenVectors() ;

	void clearCMesh(CMesh& mesh);

	void intializeSwp();

public:
	CMesh original;

	std::vector< std::vector<Point3f> > theBranches;
    skeleton_mul skel;

	sweeper swp ;

	CMesh samples;
	CMesh possion;

	//Skeleton *skeleton;

	CurveArray3D curves3d ;
	CurveArray3D snappedCurves3d ;
	CurveArray2D curves2d ;
	std::vector<double> sampleNormDiff ;

	RichParameterSet* para;
	double init_radius;
	QString curr_file_name;

	Box3f boxBeforeNormalize ;
};

