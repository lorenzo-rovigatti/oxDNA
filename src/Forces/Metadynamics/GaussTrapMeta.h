#ifndef GAUSSTRAPMETA_H_
#define GAUSSTRAPMETA_H_

#include "../BaseForce.h"
#include <vector>

class GaussTrapMeta: public BaseForce {
private:
	std::vector<int> _p1a;
	std::vector<int> _p2a;
	std::vector<int> _p1b;
	std::vector<int> _p2b;

public:
	std::vector<BaseParticle*> _p1a_ptr;
	std::vector<BaseParticle*> _p2a_ptr;
	std::vector<BaseParticle*> _p1b_ptr;
	std::vector<BaseParticle*> _p2b_ptr;
	number xmin;
	number xmax;
	int N_grid;
	number dX;
	std::vector<std::vector<number> > potential_grid;
	int _mode = 0;
	bool PBC = false;
	BaseBox *_box_ptr = nullptr;
	GaussTrapMeta();
	virtual ~GaussTrapMeta() {

	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp, BaseBox *box_ptr) override;

	LR_vector value(llint step, LR_vector &pos) override;
	number potential(llint step, LR_vector &pos) override;

protected:
	LR_vector _distance(LR_vector u, LR_vector v);
};

#endif // GAUSSTRAPMETA_H
