from setuptools import setup
import os

import sys
if sys.version_info < (3,8):
    sys.exit('Error, Python < 3.8 is not supported')

setup(
    name = 'oxpy',
    use_scm_version = {
        "root": "${CMAKE_SOURCE_DIR}",
        "fallback_version": "3.3",
        },
    setup_requires = ['setuptools-scm'],
    install_requires = [
        f"oxDNA_analysis_tools @ file://localhost/${CMAKE_SOURCE_DIR}/analysis/"
    ],
    packages = [ 'oxpy' ],
    package_dir = {
        'oxpy': '${OXPY_OUTPUT_DIR}'
    },
    package_data = {
        'oxpy': ['core.so']
    }
)
