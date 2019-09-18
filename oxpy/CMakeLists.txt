SET(oxpy_SOURCES
	OxpyContext.cpp
	OxpyManager.cpp
)

IF(Python)
	ADD_DEFINITIONS(-DPYTHON_BINDINGS)
	FIND_PACKAGE( PythonLibs 3 REQUIRED )
	INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} )
	
	ADD_SUBDIRECTORY(pybind11)
	
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/oxpy/pybind11/include )
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/src )
	
	IF(CMAKE_VERSION VERSION_LESS "3.1")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
	ELSE()
		SET(CMAKE_CXX_STANDARD 14)
		SET(CMAKE_CXX_STANDARD_REQUIRED ON)
	ENDIF()
	
	MESSAGE(STATUS "Enabling Python bindings")
	
	SET(CMAKE_SHARED_LIBRARY_PREFIX "")
	
	# we need the oxDNA "common" library to be compiled so that the resulting code is position-independent and the library can be linked dynamically
	SET_TARGET_PROPERTIES(common PROPERTIES POSITION_INDEPENDENT_CODE ON)
	
	ADD_LIBRARY(_oxpy STATIC ${oxpy_SOURCES})
	TARGET_LINK_LIBRARIES(_oxpy ${PYTHON_LIBRARIES} common)
	SET_TARGET_PROPERTIES(_oxpy PROPERTIES POSITION_INDEPENDENT_CODE ON)
	pybind11_add_module(oxpy SHARED bindings.cpp)
	TARGET_LINK_LIBRARIES(oxpy PRIVATE _oxpy)
	
	# the library should be in a folder with the same name to make it easier to use Python's distutils
	SET_TARGET_PROPERTIES(oxpy 
		PROPERTIES 
			SUFFIX ".so"
			LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/oxpy"
	)
	
	# copy the files required to install the library with Python's distutils
	CONFIGURE_FILE(__init__.py
		${CMAKE_CURRENT_BINARY_DIR}/oxpy/__init__.py
	)

	SET(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/oxpy/setup.py")
	CONFIGURE_FILE(setup.py
		${SETUP_PY}
	)

	# setup the make install target to run setup.py (in user mode, no root permissions required)
	FIND_PROGRAM(PYTHON "python3")
	INSTALL(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_PY} install --user)")
ENDIF(Python)
