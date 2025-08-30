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
    uint64_t total_points;
    SFCType sfc_type;
    uint8_t bits;
};

struct BlockMetaData {
    uint64_t offset;    // 对应的块数据相对于数据区的起始偏移量

};

#endif //DEFS_H
