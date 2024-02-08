/**
 * @file AttractionPlane.h
 * @brief Header file for the AttractionPlane class copied form Michael's file.
 * @author CopySubho
*/

#ifndef ATTRACTIONPLANE_H
#define ATTRACTIONPLANE_H

#include "BaseForce.h"

class AttractionPlane : public BaseForce {
public:
    int _particle;
    number _position;
    AttractionPlane();
    virtual ~AttractionPlane() {}

    std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos);
	virtual number potential(llint step, LR_vector &pos);
};

#endif