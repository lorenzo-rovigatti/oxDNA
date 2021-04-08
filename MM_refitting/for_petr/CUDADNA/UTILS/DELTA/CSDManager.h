/*
 * CSDManager.h
 *
 *  Created on: 26/set/2011
 *      Author: lorenzo
 */

#ifndef CSDMANAGER_H_
#define CSDMANAGER_H_

#include "../../SimManager.h"
#include "CSDBackend.h"

class CSDManager : public SimManager {
protected:
	int _N_confs;
	int _N_cyls;
	int *_csd;
	char _csd_output[256];
	int _confs_to_skip;
	CSDBackend *_csdbackend;

	virtual void _get_options();
	virtual void _load_backend(int confs_to_skip);

public:
	CSDManager(IOManager &IO, const char *input_file);
	virtual ~CSDManager();

	virtual void init();
	virtual void run();
};

#endif /* CSDMANAGER_H_ */
