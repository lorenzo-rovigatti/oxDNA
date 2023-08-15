//
// Created by jonah on 2/25/21.
//

#include "ANMParticle.h"


ANMParticle::ANMParticle() : BaseParticle()  {
    this->bonded_neighs = {};
}


ANMParticle::~ANMParticle() {

}


void ANMParticle::add_bonded_neighbor(BaseParticle *nn) {
    if(!is_bonded(nn)) {
        auto *Cq = dynamic_cast<ANMParticle *>(nn);
        bonded_neighs.push_back(Cq->index);

        Cq->bonded_neighs.push_back(this->index);

        ParticlePair new_pair(this, nn);
        this->affected.push_back(new_pair);
        nn->affected.push_back(new_pair);
    }
}



bool ANMParticle::is_bonded(BaseParticle *q) {
    std::vector<int>::iterator it;
    it = find (bonded_neighs.begin(), bonded_neighs.end(), q->index);
    return it != bonded_neighs.end();
}
