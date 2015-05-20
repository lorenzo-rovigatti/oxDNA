/*
 * AnalysisBackend.h
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#ifndef ANALYSISBACKEND_H_
#define ANALYSISBACKEND_H_

#include "SimBackend.h"

/**
 * @brief Backend used by DNAnalysis to perfor analysis on trajectories.
 *
 * @verbatim
[analysis_confs_to_skip = <int> (number of configurations that should be excluded from the analysis.)]
analysis_data_output_<n> = {\nObservableOutput\n} (specify an analysis output stream. <n> is an integer number and should start from 1. The setup and usage of output streams are documented in the ObservableOutput class.)
@endverbatim
 */
class AnalysisBackend : public SimBackend<double> {
protected:
	bool _done;
	llint _n_conf;

public:
	AnalysisBackend();
	virtual ~AnalysisBackend();

	/**
	 * @brief This method does nothing, exactly as it is supposed to do.
	 *
	 * It is only implemented because it is a purely virtual method in ISimBackend.
	 * Yes, it is a design flaw, I am aware of that :-)
	 * @param cur_step
	 */
	void sim_step(llint cur_step) {}
	void print_conf(llint cur_step, bool reduced, bool only_last) {}

	void analyse();
	bool done() { return _done; };

	void get_settings(input_file &inp);
	void init();
};

#endif /* ANALYSISBACKEND_H_ */
