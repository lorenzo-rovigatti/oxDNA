SET(oxpy_SOURCES
	OxpyContext.cpp
	OxpyManager.cpp
)

IF(Python)
	ADD_SUBDIRECTORY(pybind11)
	ADD_DEFINITIONS(-DPYTHON_BINDINGS)
	FIND_PACKAGE( Python 3 REQUIRED )
	INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} )
	MESSAGE(STATUS ${PYTHON_INCLUDE_DIRS} )
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/oxpy/pybind11/include )
	INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/src )
	
	MESSAGE(STATUS "Enabling Python bindings")
	
	SET(CMAKE_SHARED_LIBRARY_PREFIX "")
	
	# we need the oxDNA "common" library to be compiled so that the resulting code is position-independent and the library can be linked dynamically
	SET_TARGET_PROPERTIES(common PROPERTIES POSITION_INDEPENDENT_CODE ON)
	
	ADD_LIBRARY(_oxpy_lib STATIC ${oxpy_SOURCES})
	TARGET_LINK_LIBRARIES(_oxpy_lib ${PYTHON_LIBRARIES} common)
	SET_TARGET_PROPERTIES(_oxpy_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
	pybind11_add_module(core SHARED bindings.cpp)
	TARGET_LINK_LIBRARIES(core PRIVATE _oxpy_lib)
	
	# the library should be in a folder with the same name to make it easier to use Python's distutils
	SET_TARGET_PROPERTIES(core 
		PROPERTIES 
			SUFFIX ".so"
			LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/oxpy"
	)
	
	# copy the files required to install the library with Python's distutils
	CONFIGURE_FILE(__init__.py
		${CMAKE_CURRENT_BINARY_DIR}/oxpy/__init__.py
	)
	
	CONFIGURE_FILE(utils.py
		${CMAKE_CURRENT_BINARY_DIR}/oxpy/utils.py
	)
	
	# read the whole defs.h file
	FILE(READ ${PROJECT_SOURCE_DIR}/src/defs.h all_defs)
	# extract the release version
	STRING(REGEX MATCH "#define RELEASE \"([^\"]+)\"" _ ${all_defs})
	SET(release ${CMAKE_MATCH_1})
	CONFIGURE_FILE(__version__.py
		${CMAKE_CURRENT_BINARY_DIR}/oxpy/__version__.py
	)
	MESSAGE(STATUS "Parsed release: ${release}")

	SET(OXPY_DIR "${CMAKE_CURRENT_BINARY_DIR}/oxpy")
	CONFIGURE_FILE(setup.py
		"${OXPY_DIR}/setup.py"
	)

	# setup the make install target to run setup.py (in user mode, no root permissions required)
	INSTALL(CODE "execute_process(COMMAND ${Python_EXECUTABLE} -m pip install ${OXPY_DIR})")
ENDIF(Python)
