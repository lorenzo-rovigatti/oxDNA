#!/bin/bash

# Autoincrement the version, remove old files, rebuild, and reinstall
awk -F'.' '{printf("%s.%s.%s'\''",$1,$2,$3+1)}' ./src/oxDNA_analysis_tools/__init__.py > ./src/oxDNA_analysis_tools/__init__.py.tmp
mv ./src/oxDNA_analysis_tools/__init__.py.tmp ./src/oxDNA_analysis_tools/__init__.py
rm dist/*
python -m build

#Can't upload system-specific wheels to PyPi.  Use auditwheel to make a generic Linux wheel and remove my own wheel.
auditwheel repair dist/*cp39-cp39-linux_x86_64.whl
rm dist/*.whl 
mv wheelhouse/* dist/
rm -r wheelhouse

#reinstall
python -m pip install .