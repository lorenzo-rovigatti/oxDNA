/*
 * DeltaManager.h
 *
 *  Created on: 20/set/2011
 *      Author: lorenzo
 */

#ifndef DELTAMANAGER_H_
#define DELTAMANAGER_H_

#include "../../SimManager.h"
#include "DeltaBackend.h"

class DeltaManager : public SimManager {
protected:
	int _append;
	char _E_output[256];
	int _change_every;
	int _every_change_skip;
	int _confs_to_skip;
	double _threshold;
	DeltaBackend *_dbackend;

	virtual void _get_options();

public:
	DeltaManager(IOManager &IO, const char *input_file);
	virtual ~DeltaManager();

	virtual void run();
};

#endif /* DELTAMANAGER_H_ */
