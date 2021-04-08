/*
 * OrderParameters.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: sulc
 */

#include "Weights.h"

//using namespace std;

// functions for Weights class
// Constructor
Weights::Weights () {
	_dim = 0;
	_ndim = 0;
	_w = NULL;
	_sizes = NULL;
}

Weights::~Weights () {
	delete [] _w;
	delete [] _sizes;
}

void Weights::init (const char * filename, OrderParameters * op) {
	// aprire file
	ifstream inp;
	inp.open(filename, ios::in);

	fprintf (stderr, "Weight parser: parsing `%s\'\n", filename);
	
	if (!inp.good()) {
		fprintf (stderr, "Can't read file %s. Aborting\n", filename);
		abort ();
	}
	
	// initialise the array
	_ndim = op->get_hb_parameters_count ();
	_sizes = new int[_ndim];
	memcpy (_sizes, op->get_max_hb_states(), _ndim * sizeof (int));
	for (int i = 0; i < _ndim; i ++)
		_sizes[i] += 1;  // +1 to include max_value
	
	_dim = 1;
	for (int i = 0; i < _ndim; i ++) {
		_dim *= _sizes[i];
	}
	_w = new double[_dim];
	for (int i = 0; i < _dim; i ++) {
		_w[i] = 1.;
	}
	fprintf (stderr, "Weight found: %d, %d\n", _ndim, _dim);
	
	// read the file linewise...
	char line[512];
	int lineno = 1;
	while (inp.good()) {
		inp.getline(line, 512);
		if (strlen (line) == 0)
			continue;
		
		char * aux = (char *) line;
		while (isspace (* aux))
			aux ++;
		
		if (* aux == '#')
			continue;
		
		//fprintf (stderr, "Weights parser: parsing line: \"%s\"\n", aux);

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

		if (check != _ndim + 1) {
			fprintf (stderr, "Weights parser: error parsing line %d in %s. Not enough numbers in line. Aborting\n", lineno, filename);
			abort ();
		}
		
		// now we check that we are within the boundaries;
		for (int i = 0; i < _ndim; i ++) {
			if (tmp[i] < 0 || tmp[i] > (_sizes[i] + 2)) {
				fprintf (stderr, "Weights parser: error parsing line %d of `%s`': index %d out of OrderParameters bounds. Aborting\n", lineno, filename, tmp[i]);
				abort ();
			}
			if (tmpf < -DBL_EPSILON) {
				fprintf (stderr, "Weights parser: error parsing line %d of `%s`': weight %lf < 0. Cowardly refusing to proceed. Aborting\n", lineno, filename, tmpf);
				abort ();
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

		assert (index < _dim);
		_w[index] = tmpf;
		//fprintf (stderr, "w[%d]: %lf\n", index, tmpf);

		lineno ++;
	}
	
	fprintf (stderr, "Weight parser: parsing done\n");
	
	//this->print();
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
	return get_weight (arg->get_hb_states());
}

double Weights::get_weight(OrderParameters * arg, int * ptr) {
	return get_weight (arg->get_hb_states(), ptr);
}
