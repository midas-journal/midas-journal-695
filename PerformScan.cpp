#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "vtkLidarScanner.h"
#include "vtkLidarPoint.h"

#include "vtkMath.h"
#include "vtkPolyData.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkTransform.h"


int main(int argc, char* argv[])
{
	if(argc != 15)
	{
		std::cout << "Incorrect arguments! Required args:" << std::endl
				<< "InputFilename OutputFilename CreateMesh? Tx Ty Tz Rx Ry Rz (R in degrees - dictated by VTK) NumThetaPoints NumPhiPoints ThetaSpan(radians) PhiSpan(radians) StoreRays? " << std::endl;
		exit(-1);
	}
	//get the input/output filenames from the command line
	std::string InputFilename = argv[1];
	std::string OutputFilename = argv[2];
	std::string strCreateMesh = argv[3];
	std::string strTx = argv[4];
	std::string strTy = argv[5];
	std::string strTz = argv[6];
	std::string strRx = argv[7];
	std::string strRy = argv[8];
	std::string strRz = argv[9];
	std::string strThetaPoints = argv[10];
	std::string strPhiPoints = argv[11];
	std::string strThetaSpan = argv[12];
	std::string strPhiSpan = argv[13];
	std::string strStoreRays = argv[14];
	
	std::cout << "InputFilename: " << InputFilename << std::endl;
	std::cout << "OutputFilename: " << OutputFilename << std::endl;
	
	//convert string to bool
	std::stringstream ssCreateMesh(strCreateMesh);
	bool CreateMesh;
	ssCreateMesh >> CreateMesh;
	
	std::cout << "Create Mesh? " << CreateMesh << std::endl;
	
	//convert strings to doubles
	double Tx, Ty, Tz, Rx, Ry, Rz;
	
	std::stringstream ssTx(strTx);
	ssTx >> Tx;
	
	std::stringstream ssTy(strTy);
	ssTy >> Ty;
	
	std::stringstream ssTz(strTz);
	ssTz >> Tz;
	
	std::stringstream ssRx(strRx);
	ssRx >> Rx;
	
	std::stringstream ssRy(strRy);
	ssRy >> Ry;
	
	std::stringstream ssRz(strRz);
	ssRz >> Rz;
	
	std::cout << "Tx: " << Tx << std::endl;
	std::cout << "Ty: " << Ty << std::endl;
	std::cout << "Tz: " << Tz << std::endl;
	std::cout << "Rx: " << Rx << std::endl;
	std::cout << "Ry: " << Ry << std::endl;
	std::cout << "Rz: " << Rz << std::endl;
	
	
	unsigned int ThetaPoints;
	unsigned int PhiPoints;
	
	std::stringstream ssThetaPoints(strThetaPoints);
	ssThetaPoints >> ThetaPoints;
	
	std::stringstream ssPhiPoints(strPhiPoints);
	ssPhiPoints >> PhiPoints;
		
	std::cout << "Theta points: " << ThetaPoints << std::endl;
	std::cout << "Phi points: " << PhiPoints << std::endl;
	
	
	double ThetaSpan, PhiSpan;
	
	std::stringstream ssThetaSpan(strThetaSpan);
	ssThetaSpan >> ThetaSpan;
	
	std::stringstream ssPhiSpan(strPhiSpan);
	ssPhiSpan >> PhiSpan;
	
	//convert string to bool
	std::stringstream ssStoreRays(strStoreRays);
	bool StoreRays;
	ssStoreRays >> StoreRays;
	
	std::cout << "Store Rays? " << StoreRays << std::endl;
	
	//read the input vtp file
	vtkXMLPolyDataReader* reader = vtkXMLPolyDataReader::New();
	reader->SetFileName(InputFilename.c_str());
	reader->Update();
	
	//construct a vtkLidarScanner and set all of its parameters
	vtkLidarScanner* Scanner = vtkLidarScanner::New();
	
	//Scanner->WriteScanner("scanner_original.vtp");
			
	//double testAngle = vtkMath::Pi()/4.0;
	Scanner->SetPhiSpan(PhiSpan);
	Scanner->SetThetaSpan(ThetaSpan);
			
	Scanner->SetNumberOfThetaPoints(ThetaPoints);
	Scanner->SetNumberOfPhiPoints(PhiPoints);
	
	Scanner->SetStoreRays(StoreRays);
	
	//"aim" the scanner.  This is a very simple translation, but any transformation will work
	vtkTransform* transform = vtkTransform::New();
	transform->PostMultiply();
	
	transform->RotateX(Rx);
	transform->RotateY(Ry);
	transform->RotateZ(Rz);
	transform->Translate(Tx, Ty, Tz);
	
	Scanner->SetTransform(transform);
	Scanner->WriteScanner("scanner_transformed.vtp");
	
	Scanner->MakeSphericalGrid(); //indicate to use uniform spherical spacing
	
	Scanner->SetCreateMesh(CreateMesh);
	
	Scanner->SetInput(reader->GetOutput());
	Scanner->Update();
	
	vtkPolyData* poly = Scanner->GetOutput();
	
	//create a writer and write the output vtp file
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(OutputFilename.c_str());
	writer->SetInput(poly);
	writer->Write();

	
	return EXIT_SUCCESS;
}
