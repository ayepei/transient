cmake_minimum_required(VERSION 3.0)
project(OpenMOC-transient LANGUAGES CXX )
set(CMAKE_CXX_STANDARD 14)
 
set(CMAKE_BUILD_TYPE "Debug")
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/run )
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

 
 
find_package(OpenMP)
add_definitions( -DOPENMP  -fopenmp)
set(MOC_LIBS ${MOC_LIBS} OpenMP::OpenMP_CXX)

add_definitions( -DFP_PRECISION=double -DCMFD_PRECISION=double -DLINALG_TOL=1E-7 -DVEC_ALIGNMENT=16 -DVEC_LENGTH=8 )

# ---[ Includes
set( MOC_INCLUDE_DIR ${MOC_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/include/moc  )
                    
# ---[ lib
set(moc_lib "libmoc")
set(MOC_LIBS ${MOC_LIBS} )

add_subdirectory(src)
 

