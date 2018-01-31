/*
 * OrderParameters.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: sulc
 */

#include <cfloat>

#include "Weights.h"

// functions for Weights class
// Constructor
Weights::Weights () {
	_dim = -1;
	_ndim = -1;
	_w = NULL;
	_sizes = NULL;
}

Weights::~Weights () {
	delete [] _w;
	//if (_ndim > 1) delete [] _sizes;
	free (_sizes);
}

void Weights::init (const char * filename, OrderParameters * op, bool safe, double default_weight) {
	// aprire file
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
	//for (int i = 0; i < _dim; i ++)  _w[i] = 1.;
	for (int i = 0; i < _dim; i ++)  _w[i] = default_weight;
	
	OX_LOG (Logger::LOG_INFO, "(Weights.cpp) weights found; O.P. dim: %d, tot size: %d", _ndim, _dim);

	// read the file linewise...
	char line[1024];
	std::string cpp_line;
	int lineno = 1;
	while (inp.good()) {
		getline (inp, cpp_line);

		if (cpp_line.length() > 1000) throw oxDNAException ("(Weights.cpp) weight parser: error parsing line %d in %s. Lines cannot exceed 1000 characters.", lineno, filename);
		strcpy (line, cpp_line.c_str());

		if (strlen (line) == 0)
			continue;

		char * aux = (char *) line;
		while (isspace (* aux))
			aux ++;

		if (* aux == '#') continue;

		int tmp[_ndim], check = 0;
		double tmpf;
		for (int i = 0; i < _ndim; i ++) {
			check += sscanf (aux, "%d", &(tmp[i]));
			while (isalnum(* aux))
				aux ++;
			while (isspace(* aux))
				aux ++;
		}
		check += sscanf (aux, "%lf", &tmpf);

		if (check != _ndim + 1) throw oxDNAException  ("(Weights.cpp) weight parser: error parsing line %d in %s. Not enough numbers in line. Aborting", lineno, filename);

		// now we check that we are within the boundaries;
		for (int i = 0; i < _ndim; i ++) {
			if (tmp[i] < 0 || tmp[i] > (_sizes[i] + 2)) {
				throw oxDNAException ("(Weights.cpp) parser: error parsing line %d of `%s`': index %d out of OrderParameters bounds. Aborting\n", lineno, filename, tmp[i]);
			}
			if (tmpf < -DBL_EPSILON) {
				throw oxDNAException ("(Weights.cpp) parser: error parsing line %d of `%s`': weight %lf < 0. Cowardly refusing to proceed. Aborting\n", lineno, filename, tmpf);
			}
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

	// to get the real line number, we have to remove 1 since lineno was
	// initialised at 1 (not at 0 as usual)
	lineno --;
	if (lineno != _dim) {
		if (safe) throw oxDNAException ("(Weights.cpp) number of lines in the weights file do not match the dimensions of the order parameter. Expecting %d lines, got %d.\n\tUse safe_weights = False and default_weight = <float> to assume a default weight. Aborting", _dim, lineno);
		else OX_LOG (Logger::LOG_INFO, "(Weights.cpp) number of lines in the weights file do not match the dimensions of the order parameter. Expecting %d lines, got %d. Using default weight %g", _dim, lineno, default_weight);
	}
	
	OX_LOG (Logger::LOG_INFO, "(Weights.cpp) parser: parsing done"); 

	return;
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
	return;
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
