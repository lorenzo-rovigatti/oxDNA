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
class MultiDimArray {
    static constexpr size_t total = D0 * D1 * D2 * D3;
    static constexpr size_t S0 = D1 * D2 * D3;
    static constexpr size_t S1 = D2 * D3;
    static constexpr size_t S2 = D3;

    number data[total];

  public:
    inline number& operator()(size_t i0, size_t i1, size_t i2, size_t i3) noexcept {
        return data[i0*S0 + i1*S1 + i2*S2 + i3];
    }

    inline const number& operator()(size_t i0, size_t i1, size_t i2, size_t i3) const noexcept {
        return const_cast<MultiDimArray*>(this)->operator()(i0, i1, i2, i3);
    }

    void fill(const number& v) {
        std::fill(data, data + total, v);
    }

    number get_average_par() const noexcept {
        number average = 0.;
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                for(int k = 0; k < 4; k++) {
                    for(int l = 0; l < 4; l++) {
                        average += this->operator()(i, j, k, l);
                    }
                }
            }
        }
        return average / 256.;
    }
};

#endif /* OXDNA3_UTILS_H_ */