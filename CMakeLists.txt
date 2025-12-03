CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

if(POLICY CMP0054)
	cmake_policy(SET CMP0054 NEW)
endif()

set(CMAKE_ENABLE_EXPORTS ON)
project(oxDNA LANGUAGES CXX)

# this is required when working with newer gcc versions (> 5.4)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES")

OPTION(CUDA "Set to ON to compile with CUDA support" OFF)
OPTION(MPI "Set to ON to compile with MPI support" OFF)
OPTION(Debug "Set to ON to compile with debug symbols" OFF)
OPTION(G "Set to ON to compile with optimisations and debug symbols" OFF)
OPTION(INTEL "Use the Intel compiler" OFF)
OPTION(MOSIX "Make oxDNA compatible with MOSIX" OFF)
OPTION(SIGNAL "Enable SignalManager - set to OFF for OSX compatibility" OFF)
OPTION(DOUBLE "Set the numerical precision to double" ON)
OPTION(CUDA_DOUBLE "Set the numerical precision for the CUDA backend to double" OFF)
OPTION(CUDA_COMMON_ARCH "Set to OFF to autodetect the GPU and compile for its architecture, set to ON (default) to compile for the most common architectures" ON)
OPTION(NATIVE_COMPILATION "Set to OFF to compile without the -march=native flag. This may be required when compiling binaries to be used elsewhere" ON)
OPTION(JSON_ENABLED "Set to OFF to compile without JSON support, disabling the possibility of initialising external forces and/or observables from a JSON file. Disabling JSON can sometimes help with compiling with older compilers" ON)

OPTION(Python "Enable Python bindings" OFF)
OPTION(OxpySystemInstall "Install oxpy as a system-wide package" OFF)

if(DEFINED PYTHON)
	set(Python ${PYTHON})
elseif(DEFINED python)
	set(Python ${python})
endif()

# these operations have to be performed before PROJECT(oxDNA) or we will have problems at linking time
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
	SET(CMAKE_BUILD_TYPE "Debug")
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

IF(NOT DOUBLE)
	ADD_DEFINITIONS(-DFLOAT_PRECISION)
ENDIF()

include_directories(src/extern)

ADD_DEFINITIONS(-Wall -Wshadow -Wextra -Wno-unused-parameter)

IF(NOT Debug)
	IF(INTEL)
		ADD_DEFINITIONS(-fast)
	ELSE()
		ADD_DEFINITIONS(-ffast-math)
		if(NATIVE_COMPILATION)
			include(CheckCXXCompilerFlag)
            unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
            CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
            if(COMPILER_SUPPORTS_MARCH_NATIVE)
                    ADD_DEFINITIONS(-march=native)
			else()
				message(STATUS "The '-march=native' flag is not supported, compiling without it")
            endif()
		else()
			message(STATUS "Compiling without the '-march=native' flag")
		endif()
	ENDIF(INTEL)
ENDIF()

IF(MPI)
	ADD_DEFINITIONS(-DHAVE_MPI)
ENDIF(MPI)

if(CUDA)
	find_package(CUDAToolkit REQUIRED)
	if(CUDAToolkit_VERSION VERSION_GREATER_EQUAL "13.0")
		message(STATUS "CUDA ${CUDAToolkit_VERSION} detected, setting C++17 as the standard, as required for CCCL")
    	set(USE_CXX17_FOR_CUDA ON)
	endif()

	if(CMAKE_VERSION VERSION_LESS 3.18)
		FIND_PACKAGE("CUDA" REQUIRED)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_FORCE_INLINES")

		INCLUDE_DIRECTORIES(${CUDA_INSTALL_PREFIX}/include/)
		LINK_DIRECTORIES(${CUDA_INSTALL_PREFIX}/lib)

		INCLUDE(../select_compute_arch.cmake)
		
		IF(CUDA_COMMON_ARCH)
			CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS Common)
		ELSE(CUDA_COMMON_ARCH)
			CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS Auto)
		ENDIF(CUDA_COMMON_ARCH)
		
		MESSAGE(STATUS "Generating code for the following CUDA architectures: ${ARCH_FLAGS_readable}")

		LIST(APPEND CUDA_NVCC_FLAGS ${ARCH_FLAGS})
		set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};${CUDA_COMMON_FLAGS}")
	else()
		include(CheckLanguage)
		check_language(CUDA)
		if(NOT CMAKE_CUDA_COMPILER)
			message(FATAL_ERROR 
				"CUDA requested but no CUDA compiler found.\n"
				"Please install CUDA or pass -DCMAKE_CUDA_COMPILER=/path/to/nvcc"
			)
		endif()

		enable_language(CUDA)

		if(CUDA_COMMON_ARCH)
			# Query NVCC for supported architectures
			execute_process(
				COMMAND ${CMAKE_CUDA_COMPILER} --list-gpu-arch
				OUTPUT_VARIABLE NVCC_ARCH_LIST
				OUTPUT_STRIP_TRAILING_WHITESPACE
			)

			# Convert newlines into semicolons, which apparently is the way to produce a proper CMake list
			string(REPLACE "\n" ";" NVCC_ARCH_LIST "${NVCC_ARCH_LIST}")

			set(FINAL_ARCH_LIST "")
			foreach(arch ${NVCC_ARCH_LIST})
				# Match compute_XX or compute_XXX
				if(arch MATCHES "compute_([0-9]+)")
					list(APPEND FINAL_ARCH_LIST "${CMAKE_MATCH_1}")
				endif()
			endforeach()

			# Remove duplicates and sort
			list(REMOVE_DUPLICATES FINAL_ARCH_LIST)
			list(SORT FINAL_ARCH_LIST)

			# Set globally for all CUDA targets
			set(CMAKE_CUDA_ARCHITECTURES "${FINAL_ARCH_LIST}")

			message(STATUS "CUDA architectures detected via NVCC: ${CMAKE_CUDA_ARCHITECTURES}")
		else()
			if(CMAKE_VERSION VERSION_LESS 3.24)
				# Auto GPU detection through nvidia-smi, since we cannot rely on 'native' support in CMake < 3.24
				execute_process(
					COMMAND nvidia-smi --query-gpu=compute_cap --format=csv,noheader
					OUTPUT_VARIABLE GPU_ARCH
					OUTPUT_STRIP_TRAILING_WHITESPACE
					ERROR_QUIET
				)
				if(GPU_ARCH STREQUAL "")
					message(WARNING "Could not detect GPU arch, defaulting to 7.0")
					set(GPU_ARCH 70)
				else()
					# If nvidia-smi works, remove the dot: e.g. "6.1" → "61"
					MESSAGE(STATUS "Compiling for the following detected CUDA architecture: ${GPU_ARCH}")
					if(GPU_ARCH MATCHES "([0-9]+)\\.([0-9]+)")
						set(GPU_ARCH "${CMAKE_MATCH_1}${CMAKE_MATCH_2}")
					endif()
				endif()
				set(CMAKE_CUDA_ARCHITECTURES "${GPU_ARCH}")
			else()
				message(STATUS "Selecting AUTO (native) architectures")
				set(CMAKE_CUDA_ARCHITECTURES "native")
			endif()
		endif()
	endif()
endif(CUDA)

if(USE_CXX17_FOR_CUDA)
    # CUDA 13+ requires C++17 for CCCL
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
else()
    # Fallback to C++14
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CUDA_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

if(JSON_ENABLED)
	add_definitions(-DJSON_ENABLED)
else()
	message(STATUS "Disabling JSON support")
endif(JSON_ENABLED)

# get the current git revision, if git is installed. Avoid warnings if it isn't
FIND_PACKAGE(Git)
IF(Git_FOUND)
	EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} rev-parse --git-dir WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE git_dir OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
	IF(git_dir MATCHES ".git$")
		SET(IS_GIT_REPO TRUE) 
	ELSE()
		SET(IS_GIT_REPO FALSE)
	ENDIF()
ELSE()
	SET(IS_GIT_REPO FALSE)
ENDIF(Git_FOUND)

IF(IS_GIT_REPO)
	EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE revision OUTPUT_STRIP_TRAILING_WHITESPACE)
	EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0 HEAD WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE tag OUTPUT_STRIP_TRAILING_WHITESPACE)
ELSE()
	get_filename_component(dir_name ${PROJECT_SOURCE_DIR} NAME)
	string(REGEX MATCH "-([0-9]+\\.[0-9]+\\.[0-9]+.*)$" version_match ${dir_name})
	SET(revision "none")
	if(version_match)
		SET(tag "v${CMAKE_MATCH_1}")
	else()
		SET(tag "unknown")
	endif()
ENDIF(IS_GIT_REPO)

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
    ENDIF(WIN32)
ENDMACRO(CURRENT_DATE)

CURRENT_DATE(today)
ADD_DEFINITIONS(-D'RELEASE="${tag}"' -D'GIT_COMMIT="${revision}"' -D'BUILD_TIME="${today}"')

# take care of the 'make docs' target in the makefile
find_package(Doxygen)
if(DOXYGEN_FOUND)
	CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile_oxDNA ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
	CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/docs/oxDNA.css ${CMAKE_CURRENT_BINARY_DIR}/oxDNA.css @ONLY)
	ADD_CUSTOM_TARGET(docs 
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM
	)
ENDIF(DOXYGEN_FOUND)

add_custom_target(test_run
    ${CMAKE_CURRENT_SOURCE_DIR}/test/TestSuite.py test_folder_list.txt ${PROJECT_BINARY_DIR}/bin/oxDNA run
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
    COMMENT "Running build integration tests" VERBATIM
)

add_custom_target(test_quick
    ${CMAKE_CURRENT_SOURCE_DIR}/test/TestSuite.py test_folder_list.txt ${PROJECT_BINARY_DIR}/bin/oxDNA quick
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test
    COMMENT "Running quick tests" VERBATIM
)

add_custom_target(test
	COMMENT "Running all tests" VERBATIM
)

add_dependencies(test test_quick test_run)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

ADD_SUBDIRECTORY(src)
#ADD_SUBDIRECTORY(contrib)
ADD_SUBDIRECTORY(src/oxpy)
