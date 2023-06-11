//CCGParticle.cpp
// Subhajit

#include "CCGParticle.h"

CCGParticle::CCGParticle():BaseParticle(){
    spring_neighbours={};
    Bfactor={};
    ro={};
    strength=1.f;
    // lockedTo=0; // patchy particle with index 0 should not exists
    // multipatch=true; // if multiple patches could attach or not
}

CCGParticle::~CCGParticle(){
    
}

void CCGParticle::add_neighbour(BaseParticle *p){
    if(!has_bond(p)){
        auto *CCGp = static_cast<CCGParticle*>(p);
        spring_neighbours.push_back(p->index);
        CCGp->spring_neighbours.push_back(this->index);
        ParticlePair pair(this,p);
        this->affected.push_back(pair);
        p->affected.push_back(pair);
    }
}

bool CCGParticle::has_bond(BaseParticle *p){
    if(std::find(spring_neighbours.begin(),spring_neighbours.end(),p->index)!=spring_neighbours.end()){
        return true;
    }else{
        return false;
    }
}

double CCGParticle::return_bfactor(int particleIndex){
    std::vector<int>::iterator itr = std::find(spring_neighbours.begin(), spring_neighbours.end(), particleIndex);
    if(itr!=spring_neighbours.cend()){
        return Bfactor[std::distance(spring_neighbours.begin(),itr)];
    }else{
        return 0.f;
    }
}

void CCGParticle::return_kro(int particleIndex,double *k,double *r0){
    std::vector<int>::iterator itr = std::find(spring_neighbours.begin(), spring_neighbours.end(), particleIndex);
    if(itr!=spring_neighbours.cend()){
        int jindex = std::distance(spring_neighbours.begin(),itr);
        *k = 1.0/Bfactor[jindex];
        *r0 = ro[jindex];
    }else{
        *k=0;*r0=0;
    }
}