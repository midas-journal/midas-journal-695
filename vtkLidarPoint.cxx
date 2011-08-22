#include "vtkLidarPoint.h"
#include "vtkRay.h"

#include "vtkObjectFactory.h" //for new() macro

vtkStandardNewMacro(vtkLidarPoint);

vtkLidarPoint::vtkLidarPoint()
{
	//until specified, this LidarPoint is invalid
	Hit = false;
}

vtkLidarPoint::~vtkLidarPoint() 
{
}


//vtkCxxGetObjectMacro(vtkLidarPoint, Ray, vtkRay);
vtkRay* vtkLidarPoint::GetRay()
{
	return this->Ray;
}

//vtkCxxSetObjectMacro(vtkLidarPoint, Ray, vtkRay);
void vtkLidarPoint::SetRay(vtkRay* R)
{
	Ray = R;
}

double* vtkLidarPoint::GetNormal()
{
	//there is only a normal if this is a valid LidarPoint
	if(this->Hit)
		return Normal;
	else
		return NULL;
}

void vtkLidarPoint::PrintSelf(ostream &os, vtkIndent indent)
{
	//if the LidarPoint is valid, print its information
	if(this->Hit)
	{
		os << "Coordinate: " << Coordinate[0] << " " << Coordinate[1] << " " << Coordinate[2] << std::endl
			<< "Ray: " << *Ray << std::endl
			<< "Normal: " << Normal[0] << " " << Normal[1] << " " << Normal[2] << std::endl;
	}
	else
	{
		os << "Invalid!" << std::endl;
	}

}

