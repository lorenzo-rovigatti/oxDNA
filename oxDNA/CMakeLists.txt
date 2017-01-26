CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

OPTION(CUDA "Set to ON to compile with CUDA support" OFF)
OPTION(MPI "Set to ON to compile with MPI support" OFF)
OPTION(Debug "Set to ON to compile with debug symbols" OFF)
OPTION(G "Set to ON to compile with optimisations and debug symbols" OFF)
OPTION(INTEL "Use the Intel compiler" OFF)
OPTION(MOSIX "Make oxDNA compatible with MOSIX" OFF)
OPTION(SIGNAL "Enable SignalManager - set to OFF for OSX compatibility" OFF) 

# these operations have to be performed before PROJECT(oxDNA) or we will have
# problems at linking time
if(INTEL)
	FIND_PROGRAM(CMAKE_C_COMPILER NAMES icc)
	FIND_PROGRAM(CMAKE_CXX_COMPILER NAMES icpc)
	FIND_PROGRAM(CMAKE_AR NAMES xiar)
	FIND_PROGRAM(CMAKE_LINKER NAMES xild)

	if(CMAKE_C_COMPILER MATCHES CMAKE_C_COMPILER-NOTFOUND OR
	CMAKE_CXX_COMPILER MATCHES CMAKE_CXX_COMPILER-NOTFOUND OR
	CMAKE_AR MATCHES CMAKE_AR-NOTFOUND OR
	CMAKE_LINKER MATCHES CMAKE_LINKER-NOTFOUND)
	    MESSAGE(FATAL_ERROR "Cannot find Intel compiler")
	ENDIF()

	if(CUDA)
		MESSAGE(FATAL_ERROR "oxDNA with CUDA support cannot be compiled with the Intel compiler")
	ENDIF(CUDA)
ENDIF(INTEL) 

PROJECT(oxDNA)

# The default behaviour is to compile with optimisations
SET(CMAKE_BUILD_TYPE Release)

IF(Debug)
	SET(CMAKE_BUILD_TYPE Debug)
	MESSAGE(STATUS "Compiling in Debug mode")
ENDIF(Debug)

IF(G)
	SET(CMAKE_BUILD_TYPE RelWithDebInfo)
	MESSAGE(STATUS "Compiling in Release+Debug symbols mode")
ENDIF(G)

IF(MOSIX)
	ADD_DEFINITIONS(-DMOSIX)
	MESSAGE(STATUS "Disabling timings")
ENDIF(MOSIX)

IF(SIGNAL)
	ADD_DEFINITIONS(-DSIGNAL)
ENDIF(SIGNAL)

ADD_DEFINITIONS(-Wall -Wshadow)

IF(NOT CMAKE_BUILD_TYPE MATCHES Debug)
	IF(INTEL)
		ADD_DEFINITIONS(-fast)
	ELSE()
		ADD_DEFINITIONS(-ffast-math)
	ENDIF(INTEL)
ENDIF()

IF(CUDA)
	CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
	FIND_PACKAGE(CUDA)
	
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES")

	INCLUDE_DIRECTORIES(${CUDA_INSTALL_PREFIX}/include/)
	LINK_DIRECTORIES(${CUDA_INSTALL_PREFIX}/lib)
	
	IF(CUDA_VERSION_MAJOR GREATER 7)
		SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -use_fast_math --compiler-options "-ffast-math -O3" -Xcudafe "--diag_suppress=unrecognized_gcc_pragma")
	ELSE()
		SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -arch=sm_20 -use_fast_math --compiler-options "-ffast-math -O3" -Xcudafe "--diag_suppress=unrecognized_gcc_pragma")
	ENDIF(CUDA_VERSION_MAJOR GREATER 7)
ELSE()
	ADD_DEFINITIONS(-DNOCUDA)
ENDIF(CUDA)

IF(MPI)
	ADD_DEFINITIONS(-DHAVE_MPI)
ENDIF(MPI)

# get the current svn version, if svn is installed. Avoid warnings if it isn't
FIND_PACKAGE(Subversion)
IF(Subversion_FOUND)
	EXECUTE_PROCESS(COMMAND "svnversion" WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE revision OUTPUT_STRIP_TRAILING_WHITESPACE)
ELSE()
	SET(revision SVN_NOT_FOUND)
ENDIF()

# return the current date as day/month/year
MACRO(CURRENT_DATE result)
    IF(WIN32)
        EXECUTE_PROCESS(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE ${result})
        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\1/\\2/\\3" ${result} ${${result}})
    ELSEIF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%d/%m/%Y" OUTPUT_VARIABLE ${result})
        string(REGEX REPLACE "(..)/(..)/..(..).*" "\\1/\\2/\\3" ${result} ${${result}})
    ELSE(WIN32)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${result} 000000)
    ENDIF (WIN32)
ENDMACRO(CURRENT_DATE)

CURRENT_DATE(today)
ADD_DEFINITIONS(-D'SVN_VERSION="${revision}"' -D'BUILD_TIME="${today}"')

# take care of the 'make docs' target in the makefile
find_package(Doxygen)
if(DOXYGEN_FOUND)
	CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/DOCS/Doxyfile_oxDNA ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/DOCS/oxDNA.css ${CMAKE_CURRENT_BINARY_DIR}/oxDNA.css @ONLY)
	ADD_CUSTOM_TARGET(docs 
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
	)
ENDIF(DOXYGEN_FOUND)

ADD_CUSTOM_TARGET(test_run
    ${CMAKE_CURRENT_SOURCE_DIR}/TEST_LR/TestSuite.py test_folder_list.txt ${PROJECT_BINARY_DIR}/bin/oxDNA run
    COMMAND ${TEST_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/TEST_LR/TestSuite.py test_folder_list.txt ${PROJECT_BINARY_DIR}/bin/DNAnalysis analysis
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/TEST_LR
    COMMENT "Running build integration tests" VERBATIM
)

ADD_CUSTOM_TARGET(test_quick
    ${CMAKE_CURRENT_SOURCE_DIR}/TEST_LR/TestSuite.py test_folder_list.txt ${PROJECT_BINARY_DIR}/bin/oxDNA quick
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/TEST_LR
    COMMENT "Running quick tests" VERBATIM
)

ADD_CUSTOM_TARGET(test_scientific
    ${CMAKE_CURRENT_SOURCE_DIR}/TEST_LR/TestSuite.py test_folder_list.txt ${PROJECT_BINARY_DIR}/bin/oxDNA scientific
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/TEST_LR
    COMMENT "Running scientific tests" VERBATIM
)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

ADD_SUBDIRECTORY(src)
