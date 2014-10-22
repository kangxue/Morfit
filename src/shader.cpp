

#include "GLArea.h"
#include <QString>
#include <QMessageBox>

void PopupError( const QString& info )
{
	QString header = QObject::tr( "Sorry but there are some errors in using shaders.\n"
		"Here are some useful information:\n" );
	QMessageBox::warning(
		Q_NULLPTR,
		QObject::tr( "Shader program error" ),
		header + info );
}

#define CHECK( stmt, info ) \
	if ( !stmt )\
{\
	PopupError( info );\
	return;\
}

QString GetShaderFromFile( const QString& fileName )
{
	QString retStr;
	QFile file( fileName );
	if ( file.open( QIODevice::ReadOnly ) )
	{
		retStr = file.readAll( );
		file.close( );
	}
	return retStr;
}
void GLArea::SetupShaders( QString shaderName, std::vector<uniformValue> univals ){

	QString g_VertexShaderText = GetShaderFromFile( QString("") + shaderName + ".vert" );
	QString g_FragmentShaderText = GetShaderFromFile( QString("")+shaderName+".frag" );

	QOpenGLShaderProgram *m_ShaderProgram = new QOpenGLShaderProgram() ;

	QOpenGLShader *m_pVertexShader, *m_pFragmentShader ;

	m_pVertexShader = new QOpenGLShader( QOpenGLShader::Vertex,
		m_ShaderProgram );
	CHECK( m_pVertexShader->compileSourceCode(g_VertexShaderText ),
		m_pVertexShader->log( ) );

	m_pFragmentShader = new QOpenGLShader( QOpenGLShader::Fragment,
		m_ShaderProgram );

	CHECK( m_pFragmentShader->compileSourceCode( g_FragmentShaderText ),
		m_pFragmentShader->log( ) );

	CHECK( m_ShaderProgram->addShader( m_pVertexShader ),
		m_ShaderProgram->log( ) );
	CHECK( m_ShaderProgram->addShader( m_pFragmentShader ),
		m_ShaderProgram->log( ) );

	CHECK( m_ShaderProgram->link( ),
		m_ShaderProgram->log( ) );




	for( int i=0; i<univals.size(); ++i )
		m_ShaderProgram->setUniformValue( m_ShaderProgram->uniformLocation(univals[i].vname),univals[i].val ) ;

	bool found = false ;
	for( int i=0; i<shaderNames.size(); ++i ){
		if( shaderNames[i] == shaderName ){
			shaderPrograms[i] = m_ShaderProgram ;
			found = true ;
		}
	}
	if( !found ){
		shaderNames.push_back( shaderName) ;
		shaderPrograms.push_back( m_ShaderProgram ) ;
	}



}


void GLArea::bindShader( QString shaderName ) {

	for( int i=0; i<shaderNames.size(); ++i )
		if( shaderNames[i] == shaderName )
			shaderPrograms[i]->bind() ;

}

void GLArea::releaseShader() {
	for( int i=0; i<shaderPrograms.size(); ++i )
		shaderPrograms[i]->release() ;
}



void GLArea::initShader(){


	std::vector<uniformValue> univalues ;
	//univalues.push_back( uniformValue("ambient", 0.1) ) ;
	//univalues.push_back( uniformValue("edgefalloff", 0.1) ) ;
	//univalues.push_back( uniformValue("intensity", 0.1) ) ;
	//SetupShaders( "xray",univalues ) ;
	univalues.clear() ;


	//SetupShaders( "glass",univalues ) ;

	//SetupShaders( "gooch",univalues ) ;
	//SetupShaders( "gooch_light",univalues ) ;
	//SetupShaders( "phong",univalues ) ;
	//SetupShaders( "Cook-Torrance",univalues ) ;
	//SetupShaders( "slicingplane",univalues ) ;
	SetupShaders( "Oren-Nayar",univalues ) ;
	//SetupShaders( "minnaert",univalues ) ;
	//SetupShaders( "polkadot3d",univalues ) ;
	//SetupShaders( "dimple",univalues ) ;
	//SetupShaders( "pink",univalues ) ;

	


}
