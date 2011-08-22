cmake_minimum_required(VERSION 2.6)

PROJECT(vtkLidarScanner)

#################################
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

#################################

set(Sources
vtkLidarPoint.cxx
vtkLidarScanner.cxx
vtkRay.cxx
)

set(Libraries
vtkHybrid
)

ADD_LIBRARY(MyLibs ${Sources} )
TARGET_LINK_LIBRARIES(MyLibs ${Libraries} )


ENABLE_TESTING()

##########################################
ADD_EXECUTABLE(PerformScan PerformScan.cpp)
TARGET_LINK_LIBRARIES(PerformScan MyLibs)


##########################################
ADD_EXECUTABLE(TestPoint TestPoint.cpp)
TARGET_LINK_LIBRARIES(TestPoint MyLibs)
ADD_TEST(TestPoint TestPoint)

############################################
ADD_EXECUTABLE(TestScanner TestScanner.cpp)
TARGET_LINK_LIBRARIES(TestScanner MyLibs)
ADD_TEST(TestScanner TestScanner ${CMAKE_SOURCE_DIR}/XZPlane.vtp XZPlane_scan.vtp 1 0 -2 0 0 0 0 2 3 .6 .6 0)
