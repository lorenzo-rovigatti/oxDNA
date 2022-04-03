#ifndef GAUSSTRAPMETARATIO1D_H_
#define GAUSSTRAPMETARATIO1D_H_

#include "BaseForce.h"
#include <vector>

template<typename number>
class GaussTrapMetaRatio1D : public BaseForce<number> {
private:
	std::vector<int> _p1a;
	std::vector<int> _p2a;
	std::vector<int> _p1b;
	std::vector<int> _p2b;

public:
    int _ref_id;
	std::vector<BaseParticle<number> *>  _p1a_ptr;
	std::vector<BaseParticle<number> *>  _p2a_ptr;
	std::vector<BaseParticle<number> *>  _p1b_ptr;
	std::vector<BaseParticle<number> *>  _p2b_ptr;
	number xmin;
	number xmax;
    int N_grid; 
    number dX;
    std::vector<double>  potential_grid;
    int _mode;
	bool PBC;
	BaseBox<number> * _box_ptr;
	GaussTrapMetaRatio1D ();
	virtual ~GaussTrapMetaRatio1D() {}

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, BaseBox<number> *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential(llint step, LR_vector<number> &pos);

    int cutoff_int;

protected:
	LR_vector<number> _distance(LR_vector<number> u, LR_vector<number> v);
};

#endif // GAUSSTRAPMETARATIO1D_H
