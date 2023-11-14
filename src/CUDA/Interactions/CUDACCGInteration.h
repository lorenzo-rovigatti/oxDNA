#ifndef CUDACCGINTERACTION_H_
#define CUDACCGINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../CCGInteraction.h"

class CUDACCGInteraction: public CUDABaseInteraction, public CCGInteraction{
public:
    CUDACCGInteraction();
    ~CUDACCGInteraction();

protected:
}

#endif