//
// Created by jonah on 2/25/21.
//

#include "ANMTParticle.h"


LR_vector const ANMTParticle::principal_axis(1, 0, 0);
LR_vector const ANMTParticle::second_axis(0, 0, 1);
LR_vector const ANMTParticle::third_axis(0, 1, 0);

ANMTParticle::ANMTParticle() : BaseParticle()  {
    this->bonded_neighs={};
}


ANMTParticle::~ANMTParticle() {

}


void ANMTParticle::add_bonded_neighbor(BaseParticle *nn) {
    if(!is_bonded(nn)) {
        auto *Cq = dynamic_cast<ANMTParticle *>(nn);
        bonded_neighs.push_back(nn->index);
        Cq->bonded_neighs.push_back(this->index);

        ParticlePair new_pair(this, nn);
        this->affected.push_back(new_pair);
        nn->affected.push_back(new_pair);
    }
}



bool ANMTParticle::is_bonded(BaseParticle *q) {
    std::vector<int>::iterator it;
    it = find (bonded_neighs.begin(), bonded_neighs.end(), q->index);
    return it != bonded_neighs.end();
}
