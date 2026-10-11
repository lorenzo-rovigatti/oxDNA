/*
 * MetalBox.h
 *
 *  Created for Metal backend
 */

#ifndef METALBOX_H_
#define METALBOX_H_

#include "../../Boxes/BaseBox.h"
#include "../metal_defs.h"

class MetalBox {
protected:
	m_number _Lx, _Ly, _Lz;
	bool _cubic;

public:
	MetalBox() : _Lx(0), _Ly(0), _Lz(0), _cubic(false) {

	}

	MetalBox(const MetalBox &b) {
		_cubic = b._cubic;
		_Lx = b._Lx;
		_Ly = b._Ly;
		_Lz = b._Lz;
	}

	~MetalBox() {

	}
    
    // Host-side methods
	void set_Metal_from_CPU(BaseBox *box) {
		LR_vector sides = box->box_sides();
		change_sides(sides.x, sides.y, sides.z);
	}

	void set_CPU_from_Metal(BaseBox *box) {
		box->init(_Lx, _Ly, _Lz);
	}

	void change_sides(m_number nLx, m_number nLy, m_number nLz) {
		_Lx = nLx;
		_Ly = nLy;
		_Lz = nLz;
		if(nLx == nLy && nLy == nLz) _cubic = true;
	}

	m_number V() {
		return _Lx * _Ly * _Lz;
	}
    
    // Data access for copying to GPU
    struct BoxData {
        m_number box_sides[3];
        m_number inv_sides[3];
    };
    
    BoxData get_box_data() const {
        BoxData data;
        data.box_sides[0] = _Lx;
        data.box_sides[1] = _Ly;
        data.box_sides[2] = _Lz;
        data.inv_sides[0] = 1.0 / _Lx;
        data.inv_sides[1] = 1.0 / _Ly;
        data.inv_sides[2] = 1.0 / _Lz;
        return data;
    }
};

#endif /* METALBOX_H_ */
