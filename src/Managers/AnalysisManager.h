/*
 * AnalysisManager.h
 *
 *  Created on: Feb 13, 2013
 *      Author: rovigatti
 */

#ifndef ANALYSISMANAGER_H_
#define ANALYSISMANAGER_H_

#include "../defs.h"
#include "../Backends/AnalysisBackend.h"

/**
 * @brief Managing class for DNAnalysis
 *
 * This very basic class runs the analysis. It internally uses an instance of AnalysisBackend
 * to perform all the necessary analysis on the input trajectory.
 */
class AnalysisManager {
protected:
	input_file _input;
	std::shared_ptr<AnalysisBackend> _backend;

public:
	AnalysisManager(input_file input);
	virtual ~AnalysisManager();

	void load_options();
	void init();
	void analysis();
};

#endif /* ANALYSISMANAGER_H_ */
