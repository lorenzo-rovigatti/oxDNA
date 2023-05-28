//CCGParticle.cpp
// Subhajit

#include "CCGParticle.h"

CCGParticle::CCGParticle():BaseParticle(){
    this->spring_neighbours={};
    // this->color=0;
    Bfactor={};
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