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

// Helper to compute total size at compile time
template<size_t... Dims>
struct SizeHelper;

template<size_t First, size_t... Rest>
struct SizeHelper<First, Rest...> {
    static constexpr size_t value = First * SizeHelper<Rest...>::value;
};

template<>
struct SizeHelper<> {
    static constexpr size_t value = 1;
};

template<size_t... Dims>
class MultiDimArray {
public:
    static constexpr size_t N = sizeof...(Dims);
    static constexpr size_t total_size = SizeHelper<Dims...>::value;

    std::array<number, total_size> data;
    
private:
    std::array<size_t, N> _strides;
    std::array<size_t, N> _sizes = {{Dims...}};

    void compute_strides() {
        _strides[N - 1] = 1;
        for (int i = static_cast<int>(N) - 2; i >= 0; --i) {
            _strides[i] = _strides[i + 1] * _sizes[i + 1];
        }
    }

    size_t compute_index(const std::array<size_t, N>& idx) const {
        size_t offset = 0;
        for (size_t i = 0; i < N; ++i) {
            assert(idx[i] < _sizes[i]);
            offset += idx[i] * _strides[i];
        }
        return offset;
    }

public:
    MultiDimArray() {
        compute_strides();
    }

    template<typename... Indices>
    number &operator()(Indices... indices) {
        static_assert(sizeof...(Indices) == N, "Wrong number of indices");
        std::array<size_t, N> idx = {{static_cast<size_t>(indices)...}};
        return data[compute_index(idx)];
    }

    template<typename... Indices>
    const number &operator()(Indices... indices) const {
        static_assert(sizeof...(Indices) == N, "Wrong number of indices");
        std::array<size_t, N> idx = {{static_cast<size_t>(indices)...}};
        return data[compute_index(idx)];
    }

    number &operator()(const std::array<size_t, N>& idx) {
        return data[compute_index(idx)];
    }

    const number &operator()(const std::array<size_t, N>& idx) const {
        return data[compute_index(idx)];
    }

    number *raw_data() { return data.data(); }
    const number *raw_data() const { return data.data(); }

    const std::array<size_t, N>& shape() const { return _sizes; }

    void fill(number value) {
        data.fill(value);
    }

    number get_average_par() const {
        static_assert(N == 4, "Cannot take the average of an array with a number of indices != 4");
        number average = 0.;
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                for(int k = 0; k < 4; k++) {
                    for(int l = 0; l < 4; l++) {
                        average += this->operator()(i, j, k, l) / 256.;
                    }
                }
            }
        }
        return average;
    }
};

#endif /* OXDNA3_UTILS_H_ */