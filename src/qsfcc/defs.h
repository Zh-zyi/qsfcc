#ifndef DEFS_H
#define DEFS_H

#include <cstdint>
#include <vector>
#include <array>
#include <cmath>

#define SUCCESS 1
#define DIMENSION_ERROR -100
#define QUANTIZED_ERROR -200
#define SFC_ERROR -300
#define LOSSLESS_COMPRESS_ERROR -400

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
    uint8_t bits;
    uint64_t compressed_size;
};

#endif //DEFS_H
