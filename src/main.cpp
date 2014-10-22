
#include "mainwindow.h"
#include <QtWidgets/QApplication>
#include <QStyleFactory>
#include "Console.h"

//函数入口，一般不做修改
//#include <GL/glew.h>
//#include <GL/gl.h>
//#include <GL/glu.h>
//#include <GL/glut.h>
#include <fstream>
std::ofstream logout ;


int main(int argc, char *argv[])
{	

	//logout.open("log.txt") ;

	CConsoleOutput::Instance();
	//QApplication app(argc, argv);
	QApplication::setStyle(QStyleFactory::create("cleanlooks"));
	/* 
	"windows", "motif", "cde", "plastique" and "cleanlooks". Depending on the platform, "windowsxp", "windowsvista" and "macintosh" may be available.
	Note that keys are case insensitive.
	*/
	QApplication::addLibraryPath("qtplugins");


	QApplication a(argc, argv);
	MainWindow *mainWindow = new MainWindow(NULL);
	mainWindow->showMaximized();
	mainWindow->show();

	return a.exec();
}
