#Setting Version Number, Project Name
cmake_minimum_required (VERSION 2.8)
project(PHASTACHEF Fortran C CXX)
enable_language(Fortran)
enable_language(CXX)
enable_language(C)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

execute_process(
  COMMAND git rev-parse HEAD
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_FILE "${CMAKE_BINARY_DIR}/VERSION.phastaChef"
)

set(PHASTA_CHEF_ENABLED ON)

enable_testing()
include(CTest)
set(MPIRUN "mpirun"
  CACHE STRING
  "the mpirun or srun executable")
set(MPIRUN_PROCFLAG "-np"
  CACHE STRING
  "the command line flag to give process count to MPIRUN")
set(PHASTA_SRC_DIR phasta CACHE FILEPATH "path to phasta source code")
add_subdirectory(${PHASTA_SRC_DIR} ${CMAKE_BINARY_DIR}/phasta)

#test to see if simmetrix models are supported
find_package(SCOREC 2.1.0 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
if(TARGET SCOREC::gmi_sim)
  set(GMI_SIM_FOUND 1)
  message(STATUS "Found gmi_sim package - enabling Simmetrix modeler support")
  add_definitions(-DGMI_SIM_FOUND)
endif()

find_package(phastaIC PATHS ${CMAKE_BINARY_DIR})
message(STATUS "PHASTAIC_LIBS ${PHASTAIC_LIBS}")
find_package(phastaC PATHS ${CMAKE_BINARY_DIR})
message(STATUS "PHASTAC_LIBS ${PHASTAC_LIBS}")

if(phastaIC_FOUND)
set(PHASTA_INCLUDE_DIRS ${PHASTAIC_INCLUDE_DIRS})
set(PHASTA_LIBS ${PHASTAIC_LIBS})
elseif(phastaC_FOUND)
set(PHASTA_INCLUDE_DIRS ${PHASTAC_INCLUDE_DIRS})
set(PHASTA_LIBS ${PHASTAC_LIBS})
else()
 message(FATAL_ERROR "Neither the incompressible or compressible solver was found")
endif()

find_library(ACUSOLVE_LIB libles)

find_package(Armadillo REQUIRED)

add_library(phastaChef STATIC
  pcWriteFiles.cc
  pcUpdateMesh.cc
  pcAdapter.cc
  pcTimeDepMesh.cc
  pcSmooth.cc
  pcError.cc
  pcShock.cc
)
set_target_properties(phastaChef PROPERTIES LINKER_LANGUAGE HAS_CXX)
set_target_properties(phastaChef PROPERTIES LINKER_LANGUAGE Fortran)

#chef
target_link_libraries(phastaChef PUBLIC SCOREC::core)

#phasta
include_directories(${PHASTA_INCLUDE_DIRS})
target_link_libraries(phastaChef PUBLIC ${PHASTA_LIBS})

#armadillo
include_directories(${ARMADILLO_INCLUDE_DIRS})
target_link_libraries(phastaChef PUBLIC ${ARMADILLO_LIBRARIES})

macro(setup_exe exename srcname)
add_executable(${exename} ${srcname})
# Link to PHASTA_LIBS again to fix circular dependency between common/compressible.
target_link_libraries(${exename} PUBLIC phastaChef ${PHASTA_LIBS})
endmacro(setup_exe)

setup_exe(chefPhasta_posix chef_phasta_posix.cc)
setup_exe(chefPhasta_stream chef_phasta_stream.cc)
setup_exe(chefPhastaLoop_stream_ur chef_phasta_loop_stream_ur.cc)
setup_exe(chefPhastaLoop_stream_adapt chef_phasta_adaptLoop.cc)
setup_exe(chefPhastaLoop_files_adapt chef_phasta_adaptLoop_files.cc)
setup_exe(chefPhastaLoop_sam_stream_adapt chef_phasta_sam_adaptLoop.cc)
setup_exe(loopChefPhasta loopChefPhasta.cc)
setup_exe(transferAndAdapter transferAndAdapter.cc)
setup_exe(solutionProjection solutionProjection.cc)
setup_exe(calcEfficiency calcEfficiency.cc)
setup_exe(meshGrading meshGrading.cc)

add_subdirectory(test)
