/**
 * @file    MetalNoList.h
 * @date    Created for Metal backend
 * @author  antigravity
 *
 */

#ifndef METALNOLIST_H_
#define METALNOLIST_H_

#include "MetalBaseList.h"

/**
 * @brief Implements a O(N^2) type of simulation (each particle interact with each other) with Metal.
 */

class MetalNoList: public MetalBaseList {
public:
	MetalNoList();
	virtual ~MetalNoList();

	void update(id<MTLBuffer> poss, id<MTLBuffer> list_poss, id<MTLBuffer> bonds) override {}
    
    // Just for override compliance
    void metal_init(int N, m_number rcut, MetalBox *h_metal_box, id<MTLBuffer> d_metal_box, id<MTLDevice> device, id<MTLLibrary> library) override {
        MetalBaseList::metal_init(N, rcut, h_metal_box, d_metal_box, device, library);
    }
	void clean() override {}

	void get_settings(input_file &inp) override;
};

#endif /* METALNOLIST_H_ */
