from distutils.core import setup

import sys
if sys.version_info < (3,0):
    sys.exit('Error, Python < 3.0 is not supported')

setup(
    name = 'cmake_cpp_pybind11',
    version = '${revision}',
    packages = [ 'oxpy' ],
    package_dir = {
        '': '${CMAKE_CURRENT_BINARY_DIR}'
    },
    package_data = {
        '': ['core.so']
    }
)
