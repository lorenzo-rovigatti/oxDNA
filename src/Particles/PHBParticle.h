// Subhajit

//btype = color and 100 is colorless

#ifndef PHBParticle_H
#define PHBParticle_H

#include "BaseParticle.h"

class PHBParticle: public BaseParticle {
public:
    PHBParticle();
    virtual ~PHBParticle();

    //variables
    double radius=1;

    //functions
    virtual bool is_bonded(const PHBParticle *q);
};

#endif