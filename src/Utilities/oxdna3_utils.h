/**
 * @file    oxdna3_utils.h
 * @date    05/jul/2025
 * @author  lorenzo
 *
 *
 */

#ifndef OXDNA3_UTILS_H_
#define OXDNA3_UTILS_H_

#include <array>
#include <cassert>
#include <cstddef>

#include "../defs.h"

template<size_t D0, size_t D1, size_t D2, size_t D3>
struct MultiDimArray {
    static constexpr size_t total_size = D0 * D1 * D2 * D3;
    static constexpr size_t S0 = D1 * D2 * D3;
    static constexpr size_t S1 = D2 * D3;
    static constexpr size_t S2 = D3;

    number data[total_size];

    inline number& operator()(size_t i0, size_t i1, size_t i2, size_t i3) noexcept {
        return data[i0*S0 + i1*S1 + i2*S2 + i3];
    }

    inline const number& operator()(size_t i0, size_t i1, size_t i2, size_t i3) const noexcept {
        return const_cast<MultiDimArray*>(this)->operator()(i0, i1, i2, i3);
    }

    void fill(const number& v) {
        std::fill(data, data + total_size, v);
    }

    number get_average_par(int i, int j, int k, int l) const noexcept {
        number average = 0.;
        number den = 4.;
        if(i == 5 && l == 5) den = 16.;
        if(i == 5 && l != 5) {
            for(int m = 0; m < 4; m++) {
                average += this->operator()(m, j, k, l);
            }
        }
        else if(l == 5 && i != 5){
            for(int m = 0; m < 4; m++) {
                average += this->operator()(i, j, k, m);
            }
        }
        else if(l == 5 && i == 5){
            for(int m = 0; m < 4; m++) {
                for(int n = 0; n < 4; n++) average += this->operator()(m, j, k, n);
            }
        }
        return average / den;
    }
};

#endif /* OXDNA3_UTILS_H_ */
