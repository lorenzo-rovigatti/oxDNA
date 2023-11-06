// Subhajit 

#include "PHBParticle.h"

PHBParticle::PHBParticle():BaseParticle(){
    
}

PHBParticle::~PHBParticle(){

}

bool PHBParticle::is_bonded(const PHBParticle *q){
    if(type==1){//helix bundle
        return (n3==q->n5||n5==q->n3);// 3' connected 5' and vice versa
    }else{
        return false; // no connection
    }
}