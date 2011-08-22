#ifndef __vtkLidarScanner_h
#define __vtkLidarScanner_h

/* 
This Lidar scanner class models the output of the Leica HDS 3000. It acquires data in a series of vertical strips, from low to high, left to right. The output is either a pointcloud (with CreateMesh == 0) or a connected mesh (with CreateMesh == 1).

Coordinate System:

The scanner coordinate system is as follows:
z = up
y = forward
therefore (to be right handed), x = right

Theta:
The angle in the "XY" (forward-right) plane (a rotation around Z), measured from +y (Forward). It's range is -pi to pi. -pi/2 is left, pi/2 is right. This is obtained by rotating around the "up" axis.

Phi:
The elevation angle, in the YZ (forward-up) plane (a rotation around X), measured from +y. It's range is -pi/2 (down) to pi/2 (up). This is obtained by rotating around the "right" axis (AFTER the new right axis is obtained by setting Theta).
*/


#include <vector> //The main output grid is a vector of vectors.

#include "vtkSmartPointer.h" // compiler errors if this is forward declared

#include "vtkPolyDataAlgorithm.h" //superclass

class vtkPolyData;
class vtkTransform;
class vtkInformation;
class vtkInformationVector;
class vtkModifiedBSPTree;

class vtkRay;
class vtkLidarPoint;

class vtkLidarScanner : public vtkPolyDataAlgorithm
{

	public:
		static vtkLidarScanner *New();
		void PrintSelf(ostream &os, vtkIndent indent);
		
		vtkGetMacro(NumberOfThetaPoints, unsigned int);
		vtkGetMacro(NumberOfPhiPoints, unsigned int);
		vtkGetMacro(MinPhiAngle, double);
		vtkGetMacro(MaxPhiAngle, double);
		vtkGetMacro(MinThetaAngle, double);
		vtkGetMacro(MaxThetaAngle, double);
		vtkGetMacro(LOSVariance, double);
		vtkGetMacro(OrthogonalVariance, double);
		
		///////// Calculated Accessors /////////////////
		double GetPhiStep() const {return fabs(MaxPhiAngle - MinPhiAngle)/static_cast<double> (NumberOfPhiPoints - 1);}
		double GetThetaStep() const {return fabs(MaxThetaAngle - MinThetaAngle)/static_cast<double>(NumberOfThetaPoints - 1);}
		double GetNumberOfTotalPoints() const {return NumberOfPhiPoints * NumberOfThetaPoints;}
		
		vtkSetMacro(NumberOfThetaPoints, unsigned int);
		vtkSetMacro(NumberOfPhiPoints, unsigned int);
		vtkSetMacro(MinPhiAngle, double);
		vtkSetMacro(MaxPhiAngle, double);
		vtkSetMacro(MinThetaAngle, double);
		vtkSetMacro(MaxThetaAngle, double);
		vtkSetMacro(CreateMesh, bool);
		vtkSetMacro(MaxMeshEdgeLength, double);
		vtkSetMacro(StoreRays, bool);
		vtkSetMacro(LOSVariance, double);
		vtkSetMacro(OrthogonalVariance, double);
		
		void SetScene(vtkPolyData* scene);
		void SetTransform(vtkTransform* transform);
		
		void SetThetaSpan(const double theta) // (radians)
		{
			//this is a convenience function that simply divides the span by two and evenly splits the span across zero
			MinThetaAngle = -fabs(theta)/2.0;
			MaxThetaAngle = fabs(theta)/2.0;
		}

		void SetPhiSpan(const double phi)  // (radians)
		{
			//this is a convenience function that simply divides the span by two and evenly splits the span across zero
			MinPhiAngle = -fabs(phi)/2.0;
			MaxPhiAngle = fabs(phi)/2.0;
		}
		
		///////////// Accessors ////////////
		double* GetPosition() const;
		double* GetLocation() const {return this->GetPosition();} //convenience function
		vtkRay* GetRay(const double Theta, const double Phi) const;
		
		//////////////Functions///////////
		bool AcquirePoint(const unsigned int ThetaIndex, const unsigned int PhiIndex); //do a single ray/scene intersection
		
		void PerformScan(); //actually do all of the ray/scene intersections
				
		//std::vector<std::vector<vtkRay*> > MakeLinearGrid(const double DistanceToGrid); //some scanners try to set the spacing as if it were sampling a plane at a certain distance uniformly 
		
		void MakeSphericalGrid(); //use a uniform spherical spacing
		
		void GetOutputPoints(vtkPolyData* output); //put all of the valid scene intersections into a PolyData
		void GetOutputMesh(vtkPolyData* output); //put all of the valid scene intersections into a PolyData and connect them using Delaunay triangulation
		
		void AddNoise(vtkLidarPoint* Point);
		
		bool WritePTX(const std::string &Filename) const; //write a ptx file of the acquired scan
		void WriteScanner(const std::string &Filename) const; //write a vtp file of a coordinate system indicating the scanner's location and orientation
		
	protected:
		vtkLidarScanner();
		~vtkLidarScanner();
		int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *); //the function that makes this class work with the vtk pipeline
		
	private:
		static const double Forward[3]; //the direction of the "default" scanner
		static double Origin[3];
		
		unsigned int NumberOfThetaPoints; //the number of strips
		unsigned int NumberOfPhiPoints; //the number of points per strip
		double MinPhiAngle; //phi angle of the first point in each strip (radians)
		double MaxPhiAngle; //phi angle of the last point in each strip (radians)
		double MinThetaAngle; //theta angle of the first strip (radians)
		double MaxThetaAngle; //theta angle of the last strip (radians)
		
		std::vector<double> PhiAngles; // a list of the phi angles
		std::vector<double> ThetaAngles;// a list of the theta angles

		vtkTransform* Transform; //the transformation to take the scanner from its default orientation and position to the correction orientation and position
		
		bool StoreRays;
		bool CreateMesh;
		double MaxMeshEdgeLength;
	//public:
		std::vector<std::vector<vtkLidarPoint*> > OutputGrid; //outer vector is size NumberOfThetaPoints, inner vector is size NumberOfPhiPoints
		
		double LOSVariance;
		double OrthogonalVariance;
	private:
		vtkSmartPointer<vtkPolyData> Scene; //the mesh that is to be intersected
		
		typedef vtkModifiedBSPTree TreeType;
		vtkSmartPointer<TreeType> Tree; //an efficient storage of the scene
		
};


void Project(const double* A, const double* B, double* Projection);
void GetOrthogonalVector(const double* V, double* OrthogonalVector);

#endif

