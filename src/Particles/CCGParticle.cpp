//CCGParticle.cpp
// Subhajit
// Major breaking changes k = Bfactor;

#include "CCGParticle.h"

CCGParticle::CCGParticle():BaseParticle(){
    spring_neighbours={};
    Bfactor={};
    ro={};
    strength=1;
    // lockedTo=0; // patchy particle with index 0 should not exists
    // multipatch=true; // if multiple patches could attach or not
}

CCGParticle::~CCGParticle(){
    
}

void CCGParticle::add_neighbour(BaseParticle *nn,double bfact, double tempro){
    if(!has_bond(nn)){
        auto *Cq = dynamic_cast<CCGParticle *>(nn);
        spring_neighbours.push_back(Cq->index);
        Bfactor.push_back(bfact);
        ro.push_back(tempro);
        Cq->spring_neighbours.push_back(index);
        Cq->Bfactor.push_back(bfact);
        Cq->ro.push_back(tempro);

        ParticlePair newPair(this,nn);
        affected.push_back(newPair);
        nn->affected.push_back(newPair);
    }
}

// void CCGParticle::add_Bfactor()

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
        // std::cout << index << "\t"<< particleIndex<<"\t"<<(double) Bfactor[std::distance(spring_neighbours.begin(),itr)]<<"\n";
        return Bfactor[std::distance(spring_neighbours.begin(),itr)];
    }else{
        std::cout<<"error \n";
        return 0.f;
    }
}

void CCGParticle::return_kro(int particleIndex,double *k,double *r0){
    std::vector<int>::iterator itr = std::find(spring_neighbours.begin(), spring_neighbours.end(), particleIndex);
    if(itr!=spring_neighbours.cend()){
        int jindex = (int)std::distance(spring_neighbours.begin(),itr);
        *k = (double)Bfactor[jindex];
        *r0 = (double)ro[jindex];
        // std::cout<<index<<"\t"<<particleIndex<<"\t"<<*k<<"\t"<<*r0<<"\n";
    }else{
        *k=0;*r0=0;
    }
}