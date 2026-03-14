/*
 * OrderParameters.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: sulc
 */

#include "Weights.h"

#include "OrderParameters.h"

#include <fstream>
#include <sstream>

// functions for Weights class
// Constructor
Weights::Weights () {
	_dim = -1;
	_ndim = -1;
	_w = nullptr;
	_sizes = nullptr;
}

Weights::~Weights () {
	delete [] _w;
	//if (_ndim > 1) delete [] _sizes;
	free (_sizes);
}

/**
 * loads weights from a file
 * @param filename filename from which to load weights
 * @param op order parameters object for which weights will be loaded
 * @param safe if True, check that the number of weights in the file match the number of possible states defined by the order parameter.
 * @param default_weight default weight to use for states not explicitly assigned weights, if param `safe` is set to false
 */
void Weights::init (const char * filename, OrderParameters * op, bool safe, double default_weight) {
	ifstream inp;
	inp.open(filename, ios::in);

	OX_LOG(Logger::LOG_INFO, "(Weights.cpp) Weight parser: parsing `%s\'", filename);

	if (!inp.good()) throw oxDNAException ("(Weights.cpp) parser: Can't read file `%s\'. Aborting", filename);

	// initialise the array
	_ndim = op->get_all_parameters_count ();
	_sizes = (int *) calloc ((size_t)_ndim, sizeof (int));

    memcpy (_sizes, op->get_state_sizes(), ((size_t)_ndim) * sizeof (int));

    _dim = 1;
    for (int i = 0; i < _ndim; i ++) _dim *= _sizes[i];

    _w = new double[_dim];
    for (int i = 0; i < _dim; i ++)  _w[i] = default_weight;

    OX_LOG (Logger::LOG_INFO, "(Weights.cpp) weights found; O.P. dim: %d, tot size: %d", _ndim, _dim);

    std::string line;
    int lineno = 1;
    while (inp.good()) {
        getline (inp, line);

        if (line.empty()) continue;

        // Strip leading whitespace
        size_t pos = line.find_first_not_of(" \t\r\n");
        if (pos == std::string::npos) continue;
        if (line[pos] == '#') continue;

        // Parse using istringstream instead of manual sscanf/pointer walking
        std::istringstream iss(line);
        std::vector<int> tmp(_ndim);
        double tmpf;
        int check = 0;

        for (int i = 0; i < _ndim; i++) {
            if (iss >> tmp[i]) check++;
        }
        if (iss >> tmpf) check++;

        if (check != _ndim + 1) throw oxDNAException ("(Weights.cpp) weight parser: error parsing line %d in %s. Not enough numbers in line. Aborting", lineno, filename);

        // Check that we're within the order parameter boundries
        for (int i = 0; i < _ndim; i ++) {
            if (tmp[i] < 0 || tmp[i] > (_sizes[i] + 2)) {
                throw oxDNAException ("(Weights.cpp) parser: error parsing line %d of `%s`': index %d out of OrderParameters bounds. Aborting\n", lineno, filename, tmp[i]);
            }
        }
    	if (tmpf < -DBL_EPSILON) {
    		throw oxDNAException ("(Weights.cpp) parser: error parsing line %d of `%s`': weight %lf < 0. Cowardly refusing to proceed. Aborting\n", lineno, filename, tmpf);
    	}

        int index = 0;
        for (int i = 0; i < _ndim; i ++) {
            int pindex = 1;
            for (int k = 0; k < i; k ++) {
                pindex *= _sizes[k];
            }
            index += tmp[i] * pindex;
        }

        if (index < _dim) _w[index] = tmpf;
        else OX_LOG (Logger::LOG_WARNING, "(Weights.cpp) Trying to assign weight to non-existent index the order parameter. Weight file too long/inconsistent?");

        lineno ++;
    }

    lineno --;
    if (lineno != _dim) {
        if (safe) throw oxDNAException ("(Weights.cpp) number of lines in the weights file do not match the dimensions of the order parameter. Expecting %d lines, got %d.\n\tUse safe_weights = False and default_weight = <float> to assume a default weight. Aborting", _dim, lineno);
        else OX_LOG (Logger::LOG_INFO, "(Weights.cpp) number of lines in the weights file do not match the dimensions of the order parameter. Expecting %d lines, got %d. Using default weight %g", _dim, lineno, default_weight);
    }

    OX_LOG (Logger::LOG_INFO, "(Weights.cpp) parser: parsing done");
}

void Weights::print() {
	printf ("######## weights ##################\n");
	int tmp[_ndim];
	for (int i = 0; i < _dim; i ++) {
		for (int j = 0; j < _ndim; j ++) {
			int pindex = 1;
			for (int k = 0; k < j; k ++) {
				pindex *= _sizes[k];
			}
			tmp[j] = (i / pindex) % _sizes[j];
		}
		printf ("%d %lf %lf\n", i, _w[i], get_weight(tmp));
	}
	printf ("###################################\n");
}

double Weights::get_weight_by_index (int index) {
	return _w[index];
}

double Weights::get_weight(int * arg, int * ptr) {
	int index = 0;
	for (int i = 0; i < _ndim; i ++) {
		assert (arg[i] < _sizes[i]);
		int pindex = 1;
		for (int k = 0; k < i; k ++) {
			pindex *= _sizes[k];
		}
		index += arg[i] * pindex;
	}
	if (index >= _dim) {
		printf ("index > dim: %i > %i\n", index, _dim);
		for (int k = 0; k<_ndim; k ++) {
			printf ("%d ",arg[k]);
		}
		for (int k = 0; k<_ndim; k ++) {
			printf ("%d ",arg[k]);
		}
		printf ("\n");
	}
	* ptr = index;
	return _w[index];
}

double Weights::get_weight(int * arg) {
	int index = 0;
	for (int i = 0; i < _ndim; i ++) {
		assert (arg[i] < _sizes[i]);
		int pindex = 1;
		for (int k = 0; k < i; k ++) {
			pindex *= _sizes[k];
		}
		index += arg[i] * pindex;
	}
	if (index >= _dim) {
		printf ("index > dim: %i > %i\n", index, _dim);
		for (int k = 0; k<_ndim; k ++) {
			printf ("%d ",arg[k]);
		}
		for (int k = 0; k<_ndim; k ++) {
			printf ("%d ",arg[k]);
		}
		printf ("\n");
	}

	return _w[index];
}

double Weights::get_weight(OrderParameters * arg) {
	//return get_weight (arg->get_hb_states());
	return get_weight (arg->get_all_states());
}

double Weights::get_weight(OrderParameters * arg, int * ptr) {
	//return get_weight (arg->get_hb_states(), ptr);
	return get_weight (arg->get_all_states(), ptr);
}
