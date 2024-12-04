//
// Created by josh on 12/3/24.
//

#ifndef ESTIMATE_TM_PY_CUDARASPBERRYINTERACTION_H
#define ESTIMATE_TM_PY_CUDARASPBERRYINTERACTION_H

#include "RaspberryInteraction.h"
#include "CUDA/Interactions/CUDABaseInteraction.h"

#endif //ESTIMATE_TM_PY_CUDARASPBERRYINTERACTION_H

//TODO: i will probably want to make an additional base-class for cpp and cuda versions
class CUDARaspberryInteraction : public CUDABaseInteraction, public RaspberryInteraction {

};
