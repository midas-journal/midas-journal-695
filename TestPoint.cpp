#include <iostream>
#include <string>

#include "vtkLidarPoint.h"
#include "vtkRay.h"

int main(int argc, char* argv[])
{
	//delcare
	vtkLidarPoint* P = vtkLidarPoint::New();
	
	//set attributes
	P->SetHit(true);
	
	double coord[3] = {1.1, 2.1, 3.1};
	P->SetCoordinate(coord);
	
	double normal[3] = {5.1, 6.1, 7.1};
	P->SetNormal(normal);
	
	double origin[3] = {0.0, 0.0, 0.0};
	double direction[3] = {10.1, 11.1, 12.1};
	vtkRay* R = vtkRay::New();
	R->SetDirection(direction);
	R->SetOrigin(origin);
	P->SetRay(R);
	
	//get attributes
	std::cout << "P: " << *P << std::endl;
	
	bool gothit = P->GetHit();
	std::cout << "Hit: " << gothit << std::endl;
	
	double* gotcoord = P->GetCoordinate();
	std::cout << gotcoord[0] << " " << gotcoord[1] << " " << gotcoord[2] << std::endl;
	
	if((coord[0] != gotcoord[0]) || (coord[1] != gotcoord[1]) || (coord[2] != gotcoord[2]) )
	{
		return EXIT_FAILURE;
	}
	
	double* gotdirection = P->GetRay()->GetDirection();
	std::cout << gotdirection[0] << " " << gotdirection[1] << " " << gotdirection[2] << std::endl;
	if((direction[0] != gotdirection[0]) || (direction[1] != gotdirection[1]) || (direction[2] != gotdirection[2]) )
	{
		return EXIT_FAILURE;
	}
	
	double* gotnormal = P->GetNormal();
	std::cout << gotnormal[0] << " " << gotnormal[1] << " " << gotnormal[2] << std::endl;
	if((normal[0] != gotnormal[0]) || (normal[1] != gotnormal[1]) || (normal[2] != gotnormal[2]) )
	{
		return EXIT_FAILURE;
	}
	
	
	return EXIT_SUCCESS;
}
