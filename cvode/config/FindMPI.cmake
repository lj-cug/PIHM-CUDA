# - Find MPI
# This module looks for MPI (Message Passing Interface) support
# it will define the following values
#  MPI_INCLUDE_PATH = cached location of mpi.h
#  MPI_LIBRARIES    = cached list of libraries to link in (mpi mpich etc)

FIND_PATH(MPI_INCLUDE_PATH mpi.h
  PATHS /usr/local/include 
  /usr/include 
  /usr/include/mpi
  /usr/local/mpi/include
  "$ENV{ProgramFiles}/MPICH/SDK/Include"
  "$ENV{ProgramFiles}/MPICH2/include"
  "C:/Program Files/MPICH/SDK/Include"
  )

FIND_LIBRARY(MPI_LIBRARIES
  NAMES mpich2 mpi mpich 
  PATHS /usr/lib /usr/local/lib /usr/local/mpi/lib
  "$ENV{ProgramFiles}/MPICH/SDK/Lib"
  "$ENV{ProgramFiles}/MPICH2/Lib"
  "C:/Program Files/MPICH/SDK/Lib" 
  )

FIND_LIBRARY(MPI_EXTRA_LIBRARIES 
  NAMES mpi++
  PATHS /usr/lib /usr/local/lib /usr/local/mpi/lib 
  "$ENV{ProgramFiles}/MPICH/SDK/Lib"
  "C:/Program Files/MPICH/SDK/Lib" 
  DOC "If a second mpi library is necessary, specify it here.")
MARK_AS_ADVANCED(MPI_EXTRA_LIBRARIES)

IF(MPI_EXTRA_LIBRARIES)
  SET(MPI_LIBRARIES ${MPI_LIBRARIES} ${MPI_EXTRA_LIBRARIES})
ENDIF(MPI_EXTRA_LIBRARIES)

