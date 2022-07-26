from setuptools import setup

import sys
if sys.version_info < (3,0):
    sys.exit('Error, Python < 3.0 is not supported')

setup(
    name = 'oxpy',
    use_scm_version = {
        "root": "${CMAKE_SOURCE_DIR}",
        "fallback_version": "3.3",
        },
    setup_requires = ['setuptools_scm'],
    packages = [ 'oxpy' ],
    package_dir = {
        'oxpy': '${OXPY_OUTPUT_DIR}'
    },
    package_data = {
        'oxpy': ['core.so']
    }
)
