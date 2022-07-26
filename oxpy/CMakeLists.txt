SET(oxpy_SOURCES
	OxpyContext.cpp
	OxpyManager.cpp
)

IF(Python)
	# tell pybind11 to use python3
	set(PYBIND11_PYTHON_VERSION 3 CACHE STRING "Python version to use for compiling modules")
	ADD_SUBDIRECTORY(pybind11)
	ADD_DEFINITIONS(-DPYTHON_BINDINGS)
	
	INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} )
	MESSAGE(STATUS "The detected Python's include dirs are: ${PYTHON_INCLUDE_DIRS}" )
	
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
	
	set(OXPY_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/python/oxpy")
	
	# the library should be in a folder with the same name to make it easier to use Python's distutils
	SET_TARGET_PROPERTIES(core 
		PROPERTIES 
			SUFFIX ".so"
			LIBRARY_OUTPUT_DIRECTORY "${OXPY_OUTPUT_DIR}"
	)
	
	# copy the files required to install the library with Python's distutils
	CONFIGURE_FILE(__init__.py
		${OXPY_OUTPUT_DIR}/__init__.py
	)
	
	CONFIGURE_FILE(utils.py
		${OXPY_OUTPUT_DIR}/utils.py
	)
	
	# add a target to test oxpy
	add_custom_target(test_oxpy
	    ${PROJECT_SOURCE_DIR}/test/TestSuite.py test_folder_list.txt ${PYTHON_EXECUTABLE} oxpy
	    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
	    COMMENT "Running oxpy tests" VERBATIM
	)
	add_dependencies(test test_oxpy)
	
	CONFIGURE_FILE(make_install_setup.py
		${OXPY_OUTPUT_DIR}/setup.py
	)

	# install the package
	if(OxpySystemInstall)
		INSTALL(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install ${OXPY_OUTPUT_DIR})")
		INSTALL(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install ${PROJECT_SOURCE_DIR}/analysis)")
	else()
		INSTALL(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --user ${OXPY_OUTPUT_DIR})")
		INSTALL(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --user ${PROJECT_SOURCE_DIR}/analysis)")
	endif()
	
ENDIF(Python)
