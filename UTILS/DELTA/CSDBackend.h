/*
 * CSDBackend.h
 *
 *  Created on: 26/set/2011
 *      Author: lorenzo
 */

#ifndef CSDBACKEND_H_
#define CSDBACKEND_H_

#include "../../MC_CPUBackend.h"
#include "Strand.h"

class CSDBackend : public MC_CPUBackend<double> {
	friend class CSDManager;
protected:
	Strand **_cyls;
	int *_clusters;
	int *_csd;
	int _N_cyls;

	void _flip_neighs(int cyl);

public:
	CSDBackend(IOManager *IO);
	virtual ~CSDBackend();

	virtual void set_confs_to_skip(int n) { _confs_to_skip = n; };
	virtual void init(ifstream &conf_input);

	virtual void set_csd();
};

#endif /* CSDBACKEND_H_ */
