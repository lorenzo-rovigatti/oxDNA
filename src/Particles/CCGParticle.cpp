//CCGParticle.cpp
// Subhajit

#include "CCGParticle.h"

CCGParticle::CCGParticle():BaseParticle(){
    spring_neighbours={};
    Bfactor={};
    ro={};
    strength=1.f;
    lockedTo=0; // patchy particle with index 0 should not exists
}

CCGParticle::~CCGParticle(){
    
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