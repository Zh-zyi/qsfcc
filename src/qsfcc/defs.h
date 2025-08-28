#ifndef DEFS_H
#define DEFS_H

#include <cstdint>
#include <vector>
#include <array>
#include <cmath>

template <size_t D>
using mdPoint = std::array<double, D>;

template <size_t D>
using QuantizedPoint = std::array<uint64_t, D>;

struct QuantizationParams {
    double delta;   // 量化步长, 通常是 2 * epsilon
};

enum class SFCType : uint8_t {
    HILBERT,
    Z_ORDER
};

struct FileHeader {
    uint8_t dimensions;
    uint16_t block_size;
    uint8_t sfc_precision; // SFC的精度/层级。
    uint64_t total_points;
    SFCType sfc_type;
};

#endif //DEFS_H
