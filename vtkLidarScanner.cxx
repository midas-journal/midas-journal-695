#include <limits>

#include "vtkRay.h"
#include "vtkLidarPoint.h"
#include "vtkLidarScanner.h"

#include "vtkObjectFactory.h" //for new() macro
#include "vtkIdList.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkPointData.h"
#include "vtkTriangle.h"
#include "vtkLine.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkDelaunay2D.h"
#include "vtkModifiedBSPTree.h"
#include "vtkMath.h"

vtkStandardNewMacro(vtkLidarScanner);

//The scanner transform is relative to the coordinate system described at the top of the header file. We need only the Forward vector.
const double vtkLidarScanner::Forward[3] = {0.0, 1.0, 0.0};
double vtkLidarScanner::Origin[3] = {0.0, 0.0, 0.0};

vtkLidarScanner::vtkLidarScanner()
{
	//initialze everything to NULL/zero values
	Transform = vtkTransform::New(); //vtkTransform is identity by default
	
	NumberOfThetaPoints = 0;
	NumberOfPhiPoints = 0;
	MinThetaAngle = 0.0;
	MaxThetaAngle = 0.0;
	MinPhiAngle = 0.0;
	MaxPhiAngle = 0.0;
	
	//Noise mode: noiseless unless specified
	LOSVariance = 0.0;
	OrthogonalVariance = 0.0;
			
	CreateMesh = false;
	StoreRays = false;
	MaxMeshEdgeLength = std::numeric_limits<double>::infinity();; //don't throw away any edges unless this is specified
}

vtkLidarScanner::~vtkLidarScanner()
{
}

double* vtkLidarScanner::GetPosition() const 
{
	//Return the 3D coordinate of the scanners location. If the scanner transform has not been set, this value is not valid.
	if(Transform)
	{
		return this->Transform->GetPosition();
	}
	else 
	{
		return NULL;
	}
}


void vtkLidarScanner::SetScene(vtkPolyData* scene)
{
	//set the polydata that we are going to scan
	this->Scene = scene;
	this->Tree = vtkSmartPointer<TreeType>::New();
	this->Tree->SetDataSet(Scene);
	this->Tree->BuildLocator();
}


vtkCxxSetObjectMacro(vtkLidarScanner, Transform, vtkTransform);

int vtkLidarScanner::RequestData(vtkInformation *vtkNotUsed(request),
		vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
	//This function calls the scanners input and output to allow it to 
	//work in the vtk algorithm pipeline
	
  // get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and ouptut
	vtkPolyData *input = vtkPolyData::SafeDownCast(
			inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData *output = vtkPolyData::SafeDownCast(
			outInfo->Get(vtkDataObject::DATA_OBJECT()));

	this->SetScene(input);
	
	this->PerformScan();
	
	if(this->CreateMesh)
	{
		this->GetOutputMesh(output);	
	}
	else
	{
		this->GetOutputPoints(output);	
	}
	
	
	return 1;
}
bool vtkLidarScanner::AcquirePoint(const unsigned int ThetaIndex, const unsigned int PhiIndex)
{
	//This function performs a ray/mesh intersection of the ray at the specified grid location
	
	//The grid location must be in bounds.
	if(ThetaIndex >= NumberOfThetaPoints || PhiIndex >= NumberOfPhiPoints)
	{
		std::cout << "Out of range!" << std::endl;
		exit(-1);
	}
	
	//We have computed the grid of rays (in MakeSphericalGrid()) relative to a "default" scanner (i.e. Forward = (0,1,0))
	//so we have to apply the scanner's transform to the ray before casting it
	vtkRay* Ray = OutputGrid[ThetaIndex][PhiIndex]->GetRay();
	Ray->ApplyTransform(this->Transform);
	
	double t;
	double x[3];
	double pcoords[3];
	int subId;
	vtkIdType cellId;
	int hit = this->Tree->IntersectWithLine(Ray->GetOrigin(), Ray->GetPointAlong(1000.0), .01, t, x, pcoords, subId, cellId);
	/*
	std::cout << "t: " << t << std::endl;
	std::cout << "x: " << x[0] << " " << x[1] << " " << x[2] << std::endl;
	std::cout << "pcoords: " << pcoords[0] << " " << pcoords[1] << " " << pcoords[2] << std::endl;
	std::cout << "subId: " << subId << std::endl;
	std::cout << "cellId: " << cellId << std::endl;
	*/
	
	// the ray does not intersect the mesh at all, we can stop here
	if(!hit)
	{
		return false;
	}
	
	vtkTriangle* Tri = vtkTriangle::SafeDownCast(this->Scene->GetCell(cellId));
	vtkPoints* TriPoints = Tri->GetPoints();

	double n[3];
	double t0[3];
	double t1[3];
	double t2[3];

	TriPoints->GetPoint(0, t0);
	TriPoints->GetPoint(1, t1);
	TriPoints->GetPoint(2, t2);
	vtkTriangle::ComputeNormal(t0, t1, t2, n);
	
	//save the coordinate of the intersection and the normal of the triangle that was intersected in the grid
	OutputGrid[ThetaIndex][PhiIndex]->SetCoordinate(x);
	OutputGrid[ThetaIndex][PhiIndex]->SetNormal(n);
	
	//set the flag for this point indicating that there was a valid intersection
	OutputGrid[ThetaIndex][PhiIndex]->SetHit(true);

	this->AddNoise(OutputGrid[ThetaIndex][PhiIndex]);
	
	//return if this ray had a valid intersection with the scene
	return OutputGrid[ThetaIndex][PhiIndex]->GetHit();
}


void vtkLidarScanner::PerformScan()
{
	std::cout << "Performing scan..." << std::endl;
	//loop through the grid and intersect each ray with the scene.
	for(unsigned int theta = 0; theta < this->NumberOfThetaPoints; theta++)
	{
		for(unsigned int phi = 0; phi < this->NumberOfPhiPoints; phi++)
		{
			AcquirePoint(theta, phi);
		}
	}

}


void vtkLidarScanner::MakeSphericalGrid()
{

	//Make a uniformly spaced spherical grid assuming a scanner position of (0,0,0) and facing Forward (0, 1, 0).
	//This grid will later be traversed and the ray at each position will be cast into the scene to look for valid intersections.
	//(the rays will be transformed using the scanner Transform so that they come from the scanner in its current orientation)
	
	//size the grid 
	//set the number of columns
	this->OutputGrid.clear();
	this->OutputGrid.resize(this->NumberOfThetaPoints);

	//declare a column
	std::vector<vtkLidarPoint*> Column;

	for(unsigned int thetaCounter = 0; thetaCounter < NumberOfThetaPoints; thetaCounter++)
	{
		//set the number of rows
		Column.clear();
		Column.resize(NumberOfPhiPoints);

		for(unsigned int phiCounter = 0; phiCounter < NumberOfPhiPoints; phiCounter++)
		{
			//compute the (phi,theta) coordinate for the current grid location
			double phi = this->MinPhiAngle + phiCounter * this->GetPhiStep();
			double theta = this->MinThetaAngle + thetaCounter * this->GetThetaStep();

			//convert the spherical coordinates into a cartesian direction
			vtkTransform* Transform = vtkTransform::New();
			//caution - these VTK functions take parameters in degrees!
			Transform->RotateZ(-theta*180./vtkMath::Pi()); //the negative is to obtain the coordinate system we defined
			Transform->RotateX(phi*180./vtkMath::Pi());
			double* RayDir = Transform->TransformPoint(this->Forward);

			//construct a ray
			vtkRay* Ray = vtkRay::New();
			
			Ray->SetOrigin(this->Origin);
			Ray->SetDirection(RayDir);
			
			//construct a vtkLidarPoint to store in the grid and set its ray to the ray that was just computed
			vtkLidarPoint* LidarPoint = vtkLidarPoint::New();
			LidarPoint->SetRay(Ray);
			
			//store the point in the grid
			Column[phiCounter] = LidarPoint;

		}

		//save the column in the grid
		OutputGrid[thetaCounter] = Column;
	}

}

void vtkLidarScanner::GetOutputMesh(vtkPolyData* output)
{
	std::cout << "GetOutputMesh" << std::endl;
	//This function connects the raw points output into a mesh using the connectivity 
	//information from the known point acquisition order. The result is stored in 'output'.
	
	//create a polydata of the points only
	GetOutputPoints(output);
		
	//create a grid of theta/phi coordinates (keeps only connectivity, not geometry)
	vtkSmartPointer<vtkPoints> points2D = vtkSmartPointer<vtkPoints>::New();

	unsigned int PointCounter = 0;
	//for(unsigned int theta = 0; theta < OutputGrid.size(); theta++ )
	for(unsigned int theta = 0; theta < this->NumberOfThetaPoints; theta++ )
	{
		//for(unsigned int phi = 0; phi < OutputGrid[0].size(); phi++ )
		for(unsigned int phi = 0; phi < this->NumberOfPhiPoints; phi++ )
		{
			if(this->OutputGrid[theta][phi]->GetHit())
			{
				points2D->InsertNextPoint(theta, phi, 0);
				PointCounter++;
			}
		}
	}
	
	std::cout << PointCounter << " points were added." << std::endl;
		
	//add the 2d grid points to a polydata object
	vtkSmartPointer<vtkPolyData> polydata2d = vtkSmartPointer<vtkPolyData>::New();
	polydata2d->SetPoints(points2D);
	
	//triangulate the grid points
	vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();
	delaunay->SetInput(polydata2d);
	delaunay->Update();

	//get the resulting triangles from the triangulation
	vtkSmartPointer<vtkCellArray> cells = delaunay->GetOutput()->GetPolys();
	
	//create the 3d triangle array
	vtkSmartPointer<vtkCellArray> Triangles3D = vtkSmartPointer<vtkCellArray>::New();
	
	//initialize some variables
	vtkIdType npts; // the number of points in a cell
	vtkIdType* pts; //indexes to the points

	//go through all the triangles of the Delaunay triangulation and add them to the 3d polydata if they are shorter than MaxMeshEdgeLength
	cells->InitTraversal();
	while (cells->GetNextCell(npts,pts))
	{
		//get the 3 points of the current triangle
		double p0[3];
		double p1[3];
		double p2[3];
		
		//If the rays are stored, point 0 is the scanner position. If the rays are not stored, point 0 is the first scan point
		unsigned int offset;
		if(this->StoreRays)
		{
			offset = 1;
		}
		else
		{
			offset = 0;
		}
		
		unsigned int TriP0 = pts[0] + offset;
		unsigned int TriP1 = pts[1] + offset;
		unsigned int TriP2 = pts[2] + offset;
		//std::cout << "pts: " << pts[0] << " " << pts[1] << " " << pts[2] << std::endl;
		std::cout << "pts: " << TriP0 << " " << TriP1 << " " << TriP2 << std::endl;
		output->GetPoint(TriP0, p0);
		output->GetPoint(TriP1, p1);
		output->GetPoint(TriP2, p2);
		
		//throw away triangles that are bigger than a threshold
		double dist1 = vtkMath::Distance2BetweenPoints(p0, p1);
		if(dist1 > this->MaxMeshEdgeLength)
		{
			continue;
		}
		
		double dist2 = vtkMath::Distance2BetweenPoints(p1, p2);
		if(dist2 > this->MaxMeshEdgeLength)
		{
			continue;
		}
		
		double dist3 = vtkMath::Distance2BetweenPoints(p0, p2);
		if(dist3 > this->MaxMeshEdgeLength)
		{
			continue;
		}
		
     	//add the triangle to the 3d polydata
		vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		/*
		triangle->GetPointIds()->SetId(0,pts[0]);
		triangle->GetPointIds()->SetId(1,pts[1]);
		triangle->GetPointIds()->SetId(2,pts[2]);
		*/
		triangle->GetPointIds()->SetId(0,TriP0);
		triangle->GetPointIds()->SetId(1,TriP1);
		triangle->GetPointIds()->SetId(2,TriP2);
		Triangles3D->InsertNextCell(triangle);

	}
	
	//save the 3d triangles in the output polydata
	output->SetPolys(Triangles3D);
	
}

void vtkLidarScanner::GetOutputPoints(vtkPolyData* output)
{
	std::cout << "GetOutputPoints" << std::endl;
	
	//declare a geometry and topology array
	vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> Vertices = vtkSmartPointer<vtkCellArray>::New();
	
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	
	//declare an array to store the normals of the triangles that were intersected
	vtkSmartPointer<vtkDoubleArray> norms = vtkSmartPointer<vtkDoubleArray>::New();
	norms->SetNumberOfComponents(3);
	norms->SetName("Normals");

	//store the scanner position as a point so rays can be drawn between scanner and scan points
	vtkIdType ScannerLocationPid[1];
	if(StoreRays)
	{
		ScannerLocationPid[0] = Points->InsertNextPoint(this->GetLocation());
		double OriginNormal[3] = {0.0, 0.0, 1.0}; //must be a valid normal so tube filter will work (!!! should change this to the scanner's 'up' vector)
		norms->InsertNextTupleValue(OriginNormal); // have to insert a normal for this new point or it will complain "the array is too short"
	}
	
	//traverse the grid, storing valid scene intersections in the geometry/topology/normal arrays
	for(unsigned int thetaCounter = 0; thetaCounter < NumberOfThetaPoints; thetaCounter++)
	{
		for(unsigned int phiCounter = 0; phiCounter < NumberOfPhiPoints; phiCounter++)
		{
			//if the ray in this grid location had a valid scene intersection
			if(this->OutputGrid[thetaCounter][phiCounter]->GetHit())
			{
				//set the next element in the geometry/topology/normal vector
				double* p = OutputGrid[thetaCounter][phiCounter]->GetCoordinate();
				vtkIdType pid[1];
				pid[0] = Points->InsertNextPoint(p);
				Vertices->InsertNextCell(1,pid);
				
				norms->InsertNextTupleValue(OutputGrid[thetaCounter][phiCounter]->GetNormal());
				
				if(StoreRays)
				{
					vtkSmartPointer<vtkLine> Line = vtkSmartPointer<vtkLine>::New();
					Line->GetPointIds()->SetId(0,ScannerLocationPid[0]);
					Line->GetPointIds()->SetId(1,pid[0]);
					lines->InsertNextCell(Line);
				}
			}
			else
			{
				//if StoreRays is on, we will save a unit distance point in every direction so the rays can be visualized even if the target is missed
				if(StoreRays)
				{
					//set the next element in the geometry/topology/normal vector
					vtkRay* R = OutputGrid[thetaCounter][phiCounter]->GetRay();
					double* p = R->GetPointAlong(1.0);
					vtkIdType pid[1];
					pid[0] = Points->InsertNextPoint(p);
					Vertices->InsertNextCell(1,pid);
				
					double n[3] = {1.0, 0.0, 0.0}; //the normal of a miss point is not defined, so we set it to an arbitrary (1,0,0)
					norms->InsertNextTupleValue(n);
									
					vtkSmartPointer<vtkLine> Line = vtkSmartPointer<vtkLine>::New();
					Line->GetPointIds()->SetId(0,ScannerLocationPid[0]);
					Line->GetPointIds()->SetId(1,pid[0]);
					lines->InsertNextCell(Line);
				}
			}
		}
	}
	
	//save these arrays in the polydata
	output->SetPoints(Points);
	output->SetVerts(Vertices);
	if(StoreRays)
	{
		output->SetLines(lines);
	}
	output->GetPointData()->SetNormals(norms);

}


bool vtkLidarScanner::WritePTX(const std::string &Filename) const
{
	std::ofstream fout(Filename.c_str(), ios::out);

	fout << this->NumberOfPhiPoints << std::endl
			<< this->NumberOfThetaPoints << std::endl
			<< "0 0 0" << std::endl
			<< "1 0 0" << std::endl
			<< "0 1 0" << std::endl
			<< "0 0 1" << std::endl
			<< "1 0 0 0" << std::endl
			<< "0 1 0 0" << std::endl
			<< "0 0 1 0" << std::endl
			<< "0 0 0 1" << std::endl;
	
	for(unsigned int thetaCounter = 0; thetaCounter < NumberOfThetaPoints; thetaCounter++)
	{
		for(unsigned int phiCounter = 0; phiCounter < NumberOfPhiPoints; phiCounter++)
		{
			//if the ray in this grid location had a valid scene intersection
			if(this->OutputGrid[thetaCounter][phiCounter]->GetHit())
			{
				vtkLidarPoint* LP = OutputGrid[thetaCounter][phiCounter];
				double* coord = LP->GetCoordinate();
				fout << coord[0] << " " << coord[1] << " " << coord[2] << " .5 0 0 0" << std::endl;
			}
			else
			{
				fout << "0 0 0 0 0 0 0" << std::endl;
			}
			
		}
	}
	
	fout.close();
		
	return true;//write successful
}

void vtkLidarScanner::WriteScanner(const std::string &Filename) const
{
	vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
	
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->InsertNextPoint(0.0, 0.0, 0.0);
	points->InsertNextPoint(1.0, 0.0, 0.0);
	points->InsertNextPoint(0.0, 1.0, 0.0);
	points->InsertNextPoint(0.0, 0.0, 1.0);
	
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	
	vtkSmartPointer<vtkLine> Line = vtkSmartPointer<vtkLine>::New();
	//x axis
	Line->GetPointIds()->SetId(0,0);
	Line->GetPointIds()->SetId(1,1);
	lines->InsertNextCell(Line);
	
	//y axis
	Line->GetPointIds()->SetId(0,0);
	Line->GetPointIds()->SetId(1,2);
	lines->InsertNextCell(Line);
	
	//z axis
	Line->GetPointIds()->SetId(0,0);
	Line->GetPointIds()->SetId(1,3);
	lines->InsertNextCell(Line);
	
	//setup colors
	unsigned char Red[3] = {255, 0, 0};
	unsigned char Yellow[3] = {255, 255, 0};
	unsigned char Green[3] = {0, 255, 0};

	vtkSmartPointer<vtkUnsignedCharArray> Colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	Colors->SetNumberOfComponents(3);
	Colors->SetName("Colors");
	Colors->InsertNextTupleValue(Red);
	Colors->InsertNextTupleValue(Yellow);
	Colors->InsertNextTupleValue(Green);

	//add points and lines to polydata
	poly->SetPoints(points);
	poly->SetLines(lines);
	
	//add colors to lines
	poly->GetCellData()->SetVectors(Colors);
	
	vtkSmartPointer<vtkTransformPolyDataFilter> TranslateFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	TranslateFilter->SetInput(poly);
	TranslateFilter->SetTransform(this->Transform);
	TranslateFilter->Update();
 
	vtkSmartPointer<vtkPolyData> Transformed = TranslateFilter->GetOutput();

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(Filename.c_str());
	writer->SetInput(Transformed);
	writer->Write();
}

void vtkLidarScanner::AddNoise(vtkLidarPoint* Point)
{
	if(!Point->GetHit())
	{
		return; //don't change anything if the point is not valid
	}
			
	//get the original point information
	double* OriginalPoint = Point->GetCoordinate();
	vtkRay* OriginalRay = Point->GetRay();
			
	//create line of sight noise (a vector to add to the point to affect the distance that was seen)
	double* LOS = OriginalRay->GetDirection();
	double LOSX = LOS[0];
	double LOSY = LOS[1];
	double LOSZ = LOS[2];
	double LOSNoise[3];
	if(this->LOSVariance > 0.0)
	{
		double LOSLength = vtkMath::Gaussian(0.0, this->LOSVariance); //LOS noise should be zero mean
		std::cout << "LOSLength: " << LOSLength << std::endl;
		
		LOSNoise[0] = LOSLength * LOSX;
		LOSNoise[1] = LOSLength * LOSY;
		LOSNoise[2] = LOSLength * LOSZ;
	}
	else
	{
		for(unsigned int i = 0; i < 3; i++)
			LOSNoise[i] = 0.0;
	}
			
	//create orthogonal noise
	double OrthogonalNoise[3];
	if(this->OrthogonalVariance > 0.0)
	{
		double OriginalDirection[3];
		OriginalDirection[0] = LOSX;
		OriginalDirection[1] = LOSY;
		OriginalDirection[2] = LOSZ;
		GetOrthogonalVector(OriginalDirection, OrthogonalNoise);
		double OrthogonalLength = vtkMath::Gaussian(0.0, this->OrthogonalVariance); //LOS noise should be zero mean
		std::cout << "OrthogonalLength: " << OrthogonalLength << std::endl;
		for(unsigned int i = 0; i < 3; i++)
			OrthogonalNoise[i] = OrthogonalLength * OrthogonalNoise[i];
	}	
	else
	{
		for(unsigned int i = 0; i < 3; i++)
			OrthogonalNoise[i] = 0.0;
	}			
	
	//combine the noise and add it to the point
	double NoiseVector[3];
	for(unsigned int i = 0; i < 3; i++)
		NoiseVector[i] = LOSNoise[i] + OrthogonalNoise[i];
			
	double NewPoint[3];
	for(unsigned int i = 0; i < 3; i++)
		NewPoint[i] = OriginalPoint[i] + NoiseVector[i];
			
	//set the noisy point
	std::cout << "Adding noise..." << std::endl;
	Point->SetCoordinate(NewPoint);
}

void GetOrthogonalVector(const double* V, double* OrthogonalVector)
{
	//Gram Schmidt Orthogonalization
	
	//create a random vector
	double RandomVector[3] = {vtkMath::Random(0.0, 1.0), vtkMath::Random(0.0, 1.0), vtkMath::Random(0.0, 1.0)};
	vtkMath::Normalize(RandomVector);
	
	double ProjectionOfRandOnV[3];
	Project(RandomVector, V, ProjectionOfRandOnV);
	
	for(unsigned int i = 0; i < 3; i++)
		OrthogonalVector[i] = RandomVector[i] - ProjectionOfRandOnV[i];
	
	vtkMath::Normalize(OrthogonalVector);
	
}

void Project(const double* A, const double* B, double* Projection)
{
	//the projection of A on B
	double scale = vtkMath::Dot(A,B)/pow(vtkMath::Norm(B), 2);
	for(unsigned int i = 0; i < 3; i++)
		Projection[i] = scale * B[i];
	
}


////////// External Operators /////////////

void vtkLidarScanner::PrintSelf(ostream &os, vtkIndent indent)
{
	//print the scanners information when << is called
  os << "Scanner" << std::endl
		  << "-----------" << std::endl;
  if(Transform)
  {
	  std::cout << "Transform: " << *(this->Transform) << std::endl
			  << "Location: " << this->GetPosition()[0] << " " << this->GetPosition()[1] << " " << this->GetPosition()[2] << std::endl;
  }
  else
  {
	std::cout << "Transform (and hence Location) are NULL" << std::endl;
  }

   std::cout << "NumberOfThetaPoints: " << this->NumberOfThetaPoints << std::endl
		  << "NumberOfPhiPoints: " << this->NumberOfPhiPoints << std::endl
		  << "MinThetaAngle: " << this->MinThetaAngle << std::endl
		  << "MaxThetaAngle: " << this->MaxThetaAngle << std::endl
		  << "MinPhiAngle: " << this->MinPhiAngle << std::endl
		  << "MaxPhiAngle: " << this->MaxPhiAngle << std::endl << std::endl;
}
