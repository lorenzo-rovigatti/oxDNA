//
// Created by jonah on 2/25/21.
//

#include "GSParticle.h"


GSParticle::GSParticle() : BaseParticle()  {
    this->bonded_neighs = {};
}


GSParticle::~GSParticle() {

}


void GSParticle::add_bonded_neighbor(BaseParticle *nn) {
    if(!is_bonded(nn)) {
        auto *Cq = dynamic_cast<GSParticle *>(nn);
        bonded_neighs.push_back(Cq->index);

        Cq->bonded_neighs.push_back(this->index);

        ParticlePair new_pair(this, nn);
        this->affected.push_back(new_pair);
        nn->affected.push_back(new_pair);
    }
}



bool GSParticle::is_bonded(BaseParticle *q) {
    std::vector<int>::iterator it;
    it = find (bonded_neighs.begin(), bonded_neighs.end(), q->index);
    return it != bonded_neighs.end();
}
