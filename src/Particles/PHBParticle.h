// Subhajit

//btype = color and 100 is colorless

#ifndef PHBParticle_H
#define PHBParticle_H

#include "CCGParticle.h"

class PHBParticle: public CCGParticle {
public:

    //strange variables
    number th_b_0=0,beta_b_0=0,xu_bending=0.952319757,xk_bending=1.14813301,kb1=0,kb2=0.80,kt_pref=1;
    PHBParticle();
    virtual ~PHBParticle();
};

#endif