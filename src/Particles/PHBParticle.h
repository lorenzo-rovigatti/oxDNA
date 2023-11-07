// Subhajit

//btype = color and 100 is colorless

#ifndef PHBParticle_H
#define PHBParticle_H

#include "CCGParticle.h"

class PHBParticle: public CCGParticle {
public:

    //strange variables
    number th_b_0=0,beta_b_0=0;
    PHBParticle();
    virtual ~PHBParticle();
};

#endif