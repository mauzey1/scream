# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-v100.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
set(NetCDF_C_PATHS /usr/gdata/climdat/lassen_builds/libs/netcdf/netcdf-c/netcdf-c-install/ CACHE STRING "")
set(NetCDF_Fortran_PATHS /usr/gdata/climdat/lassen_builds/libs/netcdf/netcdf-f/netcdf-f-install/ CACHE STRING "")
set(SCREAM_MPIRUN_EXE "jsrun" CACHE STRING "")
set(SCREAM_MPI_NP_FLAG "--np" CACHE STRING "")
