/*
 * Histogram.cpp
 * author: flavio
 */

#include <cfloat>
#include <sstream>

#include "Histogram.h"
#include "OrderParameters.h"

Histogram::Histogram() {
	_data = NULL;
	_rdata = NULL;
	_sizes = NULL;
	_etemps = NULL;
	_erdata = NULL;
	_simtemp = 1.;
	_ntemps = 0;
	_ndim = 0;
	_dim = 0;
}

Histogram::~Histogram() {
	delete[] _data;
	delete[] _rdata;
	delete[] _sizes;
	delete[] _etemps;
	for(int i = 0; i < _ntemps; i++) {
		delete[] _erdata[i];
	}
	delete[] _erdata;
}

void Histogram::init(OrderParameters *op, double *temps, int ntemps) {
	_ndim = op->get_all_parameters_count();
	_sizes = new int[_ndim];
	memcpy(_sizes, op->get_state_sizes(), _ndim * sizeof(int));
	_dim = 1;
	for(int i = 0; i < _ndim; i++) {
		_dim *= _sizes[i];
	}
	_data = new double[_dim]();
	_rdata = new double[_dim]();

	_ntemps = ntemps;
	_etemps = new double[_ntemps];
	if(ntemps > 0) {
		if(temps != NULL) memcpy(_etemps, temps, ntemps * sizeof(double));
		_erdata = new double *[_ntemps];
		for(int i = 0; i < _ntemps; i++) {
			_erdata[i] = new double[_dim]();
		}
	}
}

void Histogram::init(OrderParameters * op, std::vector<number> temps) {
	double *dtemps = (double *) malloc(temps.size() * sizeof(double));
	for(unsigned int i = 0; i < temps.size(); i++)
		dtemps[i] = (double) temps[i];
	init(op, dtemps, temps.size());
	free(dtemps);
}

void Histogram::init(OrderParameters *op) {
	_ndim = op->get_all_parameters_count();
	_sizes = new int[_ndim];
	memcpy(_sizes, op->get_state_sizes(), _ndim * sizeof(int));
	_dim = 1;
	for(int i = 0; i < _ndim; i++) {
		_dim *= _sizes[i];
	}
	_data = new double[_dim]();
	_rdata = new double[_dim]();
}

void Histogram::init(const char *filename, OrderParameters *op, double *temps, int ntemps) {
	// aprire file
	ifstream inp;
	inp.open(filename, ios::in);

	OX_LOG(Logger::LOG_INFO, "Histogram parser: parsing `%s\'\n", filename);

	if(!inp.good()) {
		OX_LOG(Logger::LOG_INFO, "Can't read file %s. Aborting\n", filename);
		abort ();
	}

	// initialise the array
	_ndim = op->get_all_parameters_count();
	_sizes = new int[_ndim];
	memcpy(_sizes, op->get_state_sizes(), _ndim * sizeof(int));

	_dim = 1;
	for(int i = 0; i < _ndim; i++) {
		_dim *= _sizes[i];
	}

	_data = new double[_dim];
	_rdata = new double[_dim];
	for(int i = 0; i < _dim; i++) {
		_data[i] = 0.;
		_rdata[i] = 0.;
	}

	_ntemps = ntemps;
	_etemps = new double[_ntemps];
	if(ntemps > 0) {
		if(temps != NULL) memcpy(_etemps, temps, ntemps * sizeof(double));
		_erdata = new double *[_ntemps];
		for(int i = 0; i < _ntemps; i++) {
			_erdata[i] = new double[_dim]();
		}
	}

	OX_LOG(Logger::LOG_INFO, "Histogram found: %d, %d %d\n", _ndim, _dim, _ntemps);

	// read the file linewise...
	char line[2048];
	std::string cpp_line;
	int lineno = 1;
	inp.getline(line, 2048); // to remove the first line
	while(inp.good()) {
		getline(inp, cpp_line);

		if(cpp_line.length() == 0) continue;
		if(cpp_line.length() > 2000) throw oxDNAException("(Histogram.cpp) Histogram parser: error parsing line %d in %s. Too many characters in line (2000 max). Aborting\n", lineno, filename);

		strcpy(line, cpp_line.c_str());

		char * aux = (char *) line;
		while(isspace(*aux))
			aux++;

		if(*aux == '#') continue;

		//OX_LOG(Logger::LOG_INFO, "Histogram parser: parsing line: \"%s\"\n", aux);

		int tmp[_ndim], check = 0;
		double tmpf, tmpf2, tmpdata[_ntemps + 2];
		for(int i = 0; i < _ndim; i++) {
			check += sscanf(aux, "%d", &(tmp[i]));
			while(isalnum(*aux))
				aux++;
			while(isspace(*aux))
				aux++;
		}
		for(int i = 0; i < _ntemps + 2; i++) {
			check += sscanf(aux, "%lf", &(tmpdata[i]));
			//while (isalnum(* aux))
			while(!isspace(*aux))
				aux++;
			while(isspace(*aux))
				aux++;
		}
		/*
		 OX_LOG(Logger::LOG_INFO, "got : ");
		 for (int i = 0; i < _ntemps + 2; i ++) {
		 OX_LOG(Logger::LOG_INFO, "%lf ", tmpdata[i]);
		 }		OX_LOG(Logger::LOG_INFO, "\n");*/

		tmpf = tmpdata[0];
		tmpf2 = tmpdata[1];

		if(check != _ndim + 2 + _ntemps) {
			throw oxDNAException("Histogram parser: error parsing line %d in %s. Not enough numbers in line. Aborting\n", lineno, filename);
		}

		// now we check that we are within the boundaries;
		for(int i = 0; i < _ndim; i++) {
			if(tmp[i] < 0 || tmp[i] > (_sizes[i] + 2)) {
				throw oxDNAException("Histogram parser: error parsing line %d of `%s\': index %d out of OrderParameters bounds. Aborting\n", lineno, filename, tmp[i]);
			}
			if(tmpf < -DBL_EPSILON) {
				throw oxDNAException("Histogram parser: error parsing line %d of `%s\': weight %lf < 0. Cowardly refusing to proceed. Aborting\n", lineno, filename, tmpf);
			}
		}

		int index = 0;
		for(int i = 0; i < _ndim; i++) {
			int pindex = 1;
			for(int k = 0; k < i; k++) {
				pindex *= _sizes[k];
			}
			index += tmp[i] * pindex;
		}

		assert(index < _dim);
		_data[index] = tmpf;
		_rdata[index] = tmpf2;

		for(int i = 2; i < _ntemps + 2; i++) {
			_erdata[i - 2][index] = tmpdata[i];
		}
		//OX_LOG(Logger::LOG_INFO, "w[%d]: %lf\n", index, tmpf);

		lineno++;
	}

	OX_LOG(Logger::LOG_INFO, "Histogram parser: parsing done\n");

	//this->print();
	return;
}

void Histogram::init(const char *filename, OrderParameters *op) {
	// aprire file
	ifstream inp;
	inp.open(filename, ios::in);

	OX_LOG(Logger::LOG_INFO, "Histogram parser: parsing `%s\'\n", filename);

	if(!inp.good()) {
		OX_LOG(Logger::LOG_INFO, "Can't read file %s. Aborting\n", filename);
		abort ();
	}

	// initialise the array
	_ndim = op->get_all_parameters_count();
	_sizes = new int[_ndim];
	memcpy(_sizes, op->get_state_sizes(), _ndim * sizeof(int));

	_dim = 1;
	for(int i = 0; i < _ndim; i++) {
		_dim *= _sizes[i];
	}

	_data = new double[_dim];
	_rdata = new double[_dim];
	for(int i = 0; i < _dim; i++) {
		_data[i] = 0.;
		_rdata[i] = 0.;
	}

	//OX_LOG(Logger::LOG_INFO, "Histogram found: %d, %d\n", _ndim, _dim);
	char line[2048];
	std::string cpp_line;
	int lineno = 1;
	inp.getline(line, 2048); // to remove the first line
	while(inp.good()) {
		getline(inp, cpp_line);

		if(cpp_line.length() == 0) continue;
		if(cpp_line.length() > 2000) throw oxDNAException("(Histogram.cpp) Histogram parser: error parsing line %d in %s. Too many characters in line (2000 max). Aborting\n", lineno, filename);

		strcpy(line, cpp_line.c_str());

		char * aux = (char *) line;
		while(isspace(*aux))
			aux++;

		if(*aux == '#') continue;

		OX_LOG(Logger::LOG_INFO, "Histogram parser: parsing line: \"%s\"\n", aux);

		int tmp[_ndim], check = 0;
		double tmpf, tmpf2;
		for(int i = 0; i < _ndim; i++) {
			check += sscanf(aux, "%d", &(tmp[i]));
			while(isalnum(*aux))
				aux++;
			while(isspace(*aux))
				aux++;
		}
		check += sscanf(aux, "%lf %lf", &tmpf, &tmpf2);

		if(check != _ndim + 2) {
			throw oxDNAException("Histogram parser: error parsing line %d in %s. Not enough numbers in line. Aborting\n", lineno, filename);
		}

		// now we check that we are within the boundaries;
		for(int i = 0; i < _ndim; i++) {
			if(tmp[i] < 0 || tmp[i] > (_sizes[i] + 2)) {
				throw oxDNAException("Histogram parser: error parsing line %d of `%s\': index %d out of OrderParameters bounds. Aborting\n", lineno, filename, tmp[i]);
			}
			if(tmpf < -DBL_EPSILON) {
				throw oxDNAException("Histogram parser: error parsing line %d of `%s\': weight %lf < 0. Cowardly refusing to proceed. Aborting\n", lineno, filename, tmpf);
			}
		}

		int index = 0;
		for(int i = 0; i < _ndim; i++) {
			int pindex = 1;
			for(int k = 0; k < i; k++) {
				pindex *= _sizes[k];
			}
			index += tmp[i] * pindex;
		}

		assert(index < _dim);
		_data[index] = tmpf;
		_rdata[index] = tmpf2;
		//OX_LOG(Logger::LOG_INFO, "w[%d]: %lf\n", index, tmpf);

		lineno++;
	}

	OX_LOG(Logger::LOG_INFO, "Histogram parser: parsing done\n");

	//this->print();
	return;
}

void Histogram::reset() {
	for(int k = 0; k < _dim; k++) {
		_data[k] = _rdata[k] = 0.;
	}
	return;
}

/*
 void Histogram::add (int i, double w, double e_state, double e_stack) {

 double e0, et;
 // etot = e0 + T * et;
 et = e_stack * STCK_FACT_EPS / (STCK_BASE_EPS + _simtemp * STCK_FACT_EPS);
 e0 = e_state - _simtemp * et;

 _data[i] += 1.;
 if (w > 0) {
 _rdata[i] += 1./w;
 for (int k = 0; k < _ntemps; k ++)
 _erdata[k][i] += exp (- (e0 + _etemps[k] * et) / _etemps[k] + (e0 + _simtemp * et) / _simtemp) / w;
 }
 else {
 _rdata[i] += 1.;
 for (int k = 0; k < _ntemps; k ++)
 _erdata[k][i] += exp (- (e0 + _etemps[k] * et) / _etemps[k] + (e0 + _simtemp * et) / _simtemp);
 }
 }
 */

void Histogram::add(int i, double w, double e_state, double e_stack, double e_ext) {

	double e0, et;
	if(_oxDNA2_stacking) et = e_stack * STCK_FACT_EPS_OXDNA2 / (STCK_BASE_EPS_OXDNA2 + _simtemp * STCK_FACT_EPS_OXDNA2);
	else et = e_stack * STCK_FACT_EPS_OXDNA / (STCK_BASE_EPS_OXDNA + _simtemp * STCK_FACT_EPS_OXDNA);
	e0 = e_state - _simtemp * et;

	_data[i] += 1.;
	if(w > 0) {
		_rdata[i] += 1. / w;
		for(int k = 0; k < _ntemps; k++)
			_erdata[k][i] += exp(-(e0 + e_ext + _etemps[k] * et) / _etemps[k] + (e0 + e_ext + _simtemp * et) / _simtemp) / w;
	}
	else {
		_rdata[i] += 1.;
		for(int k = 0; k < _ntemps; k++)
			_erdata[k][i] += exp(-(e0 + e_ext + _etemps[k] * et) / _etemps[k] + (e0 + e_ext + _simtemp * et) / _simtemp);
	}
}

void Histogram::add(int i, double am, double w, double e_state, double e_stack, double e_ext) {
	double e0, et;
	// etot = e0 + T * et;
	if(_oxDNA2_stacking) et = e_stack * STCK_FACT_EPS_OXDNA2 / (STCK_BASE_EPS_OXDNA2 + _simtemp * STCK_FACT_EPS_OXDNA2);
	else et = e_stack * STCK_FACT_EPS_OXDNA / (STCK_BASE_EPS_OXDNA + _simtemp * STCK_FACT_EPS_OXDNA);
	e0 = e_state - _simtemp * et;

	_data[i] += am;
	if(w > 0) {
		_rdata[i] += 1. / w;
		for(int k = 0; k < _ntemps; k++)
			_erdata[k][i] += am * exp(-(e0 + e_ext + _etemps[k] * et) / _etemps[k] + (e0 + e_ext + _simtemp * et) / _simtemp) / w;
	}
	else {
		_rdata[i] += 1.;
		for(int k = 0; k < _ntemps; k++)
			_erdata[k][i] += am * exp(-(e0 + e_ext + _etemps[k] * et) / _etemps[k] + (e0 + e_ext + _simtemp * et) / _simtemp);
	}
}

void Histogram::add(int i, int am, double w, double e_state, double e_stack, double e_ext) {
	this->add(i, (double) am, w, e_state, e_stack, e_ext);
}

// fix here, this function does not work
std::string Histogram::print_to_string(bool skip_zeros) {
	std::stringstream my_stream;

	int tmp[_ndim];
	for(int i = 0; i < _dim; i++) {
		if(_data[i] > 0 || skip_zeros == false) {
			for(int j = 0; j < _ndim; j++) {
				int pindex = 1;
				for(int k = 0; k < j; k++) {
					pindex *= _sizes[k];
				}
				tmp[j] = (i / pindex) % _sizes[j];
				my_stream << tmp[j] << " ";
			}
			my_stream << _data[i] << " " << _rdata[i] << " ";
			for(int k = 0; k < _ntemps; k++)
				my_stream << _erdata[k][i] << " ";
			my_stream << endl;
		}
	}

	return my_stream.str();
}

void Histogram::print_to_file(const char * filename, long long int time, bool only_last, bool skip_zeros) {
	FILE * outfile;

	if(only_last) {
		outfile = fopen(filename, "w");
	}
	else {
		outfile = fopen(filename, "a");
	}

	if(outfile == NULL) throw oxDNAException("could not open histogram file `%s\' for writing. Aborting.\n", filename);
	fprintf(outfile, "#t = %llu; extr. Ts: ", time);
	for(int k = 0; k < _ntemps; k++)
		fprintf(outfile, "%g ", _etemps[k]);
	fprintf(outfile, "\n");
	// the format is: index index ... index weight reweight
	int tmp[_ndim];
	for(int i = 0; i < _dim; i++) {
		if(_data[i] > 0 || skip_zeros == false) {
			for(int j = 0; j < _ndim; j++) {
				int pindex = 1;
				for(int k = 0; k < j; k++) {
					pindex *= _sizes[k];
				}
				tmp[j] = (i / pindex) % _sizes[j];
				fprintf(outfile, "%d ", tmp[j]);
			}
			fprintf(outfile, "%g %g ", _data[i], _rdata[i]);
			for(int k = 0; k < _ntemps; k++)
				fprintf(outfile, "%g ", _erdata[k][i]);
			fprintf(outfile, "\n");
		}
	}

	fclose(outfile);

	return;
}

void Histogram::read_interaction(input_file &inp) {
	// need this to set the stacking strength
	_oxDNA2_stacking = false;
	std::string inter_type("");
	getInputString(&inp, "interaction_type", inter_type, 0);
	if(inter_type.compare("DNA2") == 0) _oxDNA2_stacking = true;
}
