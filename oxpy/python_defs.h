/*
 * python_defs.h
 *
 *  Created on: Mar 9, 2019
 *      Author: lorenzo
 */

#ifndef PYTHON_DEFS_H_
#define PYTHON_DEFS_H_

// gcc spits out *a lot* of warnings when compiling pybind11 code with -Wshadow
#pragma GCC diagnostic ignored "-Wshadow"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;

#endif
