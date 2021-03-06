cmake_minimum_required(VERSION 2.8)
project(fccCaloAnalysis CXX)

# Set up C++14
set(CPP14FLAGS "-std=c++14")
if (${APPLE})
    set(CPP14FLAGS "-std=c++14\ -stdlib=libc++")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DDROP_CGAL ${CPP14FLAGS} -Wall -Wextra -Wpedantic -Wno-unused-variable -Wno-unused-parameter")

# Make sure we find the Find*.cmake functions distributed with this package
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_PREFIX_PATH $ENV{PODIO} $ENV{FCCEDM})

find_package(podio REQUIRED HINTS $ENV{PODIO})
message(STATUS "Found podio: ${podio_DIR}")

find_library(PODIOLIB podio)
if (NOT PODIOLIB)
    message(FATAL_ERROR "libpodio.so(dylib) cannot be found dynamically. Make sure you have sourced PODIO init*.sh file to set up your environment to use PODIO")
endif()
link_directories(${podio_LIBRARY_DIR})

find_package(fccedm REQUIRED)
message(STATUS "Found fcc-edm: ${fccedm_DIR}")

find_package(ROOT REQUIRED)
message(STATUS "Found ROOT: ${ROOT_DIR}")

# set up include-directories
include_directories(${PROJECT_SOURCE_DIR}/include/
  ${ROOT_INCLUDE_DIR}
  ${podio_INCLUDE_DIRS}
  ${FCCEDM_INCLUDE_DIRS}
  )

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
ROOT_GENERATE_DICTIONARY(CaloAnalysisDict
   TreeAnalysis.h
   BaseAnalysis.h
   MultiFileAnalysis.h
   SimpleAnalysis.h
   SimpleAnalysis_standalone.h
   ShowerProfiles.h
   CombinedShowerProfiles.h
   CombinedShowerProfilesPbSpacer.h
   ChiSquareMinimisationClusterBarrel.h
   ChiSquareMinimisationBarrel.h
   ChiSquareMinimisationEdepBarrel.h
   ChiSquareMinimisationEndcap.h
   CellAnalysis.h
   CellSignificances.h
   BasicClusterAnalysis.h
   ClusterAnalysis.h
   ClusterAnalysisTree.h
   CombinedCellAnalysis.h
   CombinedCellAnalysisPbSpacer.h
   CombinedCellAnalysisBenchmarkBarrel.h
   CombinedCellAnalysisEndcap.h
   CombinedCellAnalysisFwd.h
   ReconstructionExample.h
   SingleParticleRecoMonitors.h
   Decoder.h
   LINKDEF LinkDef.h)

add_library(CaloAnalysis SHARED ${sources} ${headers} CaloAnalysisDict.cxx)

target_link_libraries(CaloAnalysis
   datamodel datamodelDict podio
   ${ROOT_LIBRARIES} ${ROOT_COMPONENT_LIBRARIES})

install(TARGETS CaloAnalysis DESTINATION lib)

# Install the dictionary
if (${ROOT_VERSION} GREATER 6)
  install(FILES
      "${PROJECT_BINARY_DIR}/CaloAnalysisDict_rdict.pcm"
      DESTINATION lib COMPONENT dev)
endif()

