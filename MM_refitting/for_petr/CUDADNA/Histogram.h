/*
 * Histogram class; designed by the poor object-oriented thinking of a
 * (happy) C programmer
 */

#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <cstdlib>
#include <cstdio>

class OrderParameters;

class Histogram {
	private:
		double * _data;
		double * _rdata;
		int _dim;
		int _ndim;
		int * _sizes;

		int _ntemps;
		double * _etemps;
		double ** _erdata;
		double _simtemp;
	
	
	public:
		Histogram();
		~Histogram();

		void init (OrderParameters *);
		void init (OrderParameters *, double *, int);
		void init (const char *, OrderParameters *);
		void init (const char *, OrderParameters *, double *, int);

		void reset();
		void set_simtemp (double arg) { _simtemp = arg; }
		
		void add(int index, double w, double e_state, double e_stack);
		void add(int index, int amount, double w, double e_state, double e_stack);
		void add(int index, double amount, double w, double e_state, double e_stack);
		void print();
		void print_to_file (const char * filename, long long int time, bool only_last);
		void load_from_file (const char * filename);
};

#endif

