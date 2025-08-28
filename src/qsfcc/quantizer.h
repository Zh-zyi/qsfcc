#ifndef QUANTIZER_H
#define QUANTIZER_H

#include "defs.h"


template <size_t D>
QuantizedPoint<D> quantize_point(const mdPoint<D>& point, const std::vector<QuantizationParams>& q_params);

template <size_t D>
mdPoint<D> dequantize_point(const QuantizedPoint<D>& q_point, const std::vector<QuantizationParams>& q_params);

inline uint64_t zigzag_encode(int64_t n) {
    return (static_cast<uint64_t>(n) << 1) ^ (n >> 63);
}

inline int64_t zigzag_decode(uint64_t n) {
    return (n >> 1) ^ (-(n & 1));
}

// template <size_t D>
// QuantizedPoint<D> quantize_point(const Point<D>& point, const std::vector<QuantizationParams>& q_params) {
//     QuantizedPoint<D> q_point;
//     for (size_t i = 0; i < D; ++i) {
//         auto signed_quantized_index = static_cast<int64_t>(round(point[i] / q_params[i].delta));
//         q_point[i] = zigzag_encode(signed_quantized_index);
//     }
//     return q_point;
// }
//
// template <size_t D>
// Point<D> dequantize_point(const QuantizedPoint<D>& q_point, const std::vector<QuantizationParams>& q_params) {
//     Point<D> point;
//     for (size_t i = 0; i < D; ++i) {
//         int64_t signed_quantized_index = zigzag_decode(q_point[i]);
//         point[i] = static_cast<double>(signed_quantized_index) * q_params[i].delta;
//     }
//     return point;
// }

inline std::vector<uint64_t> quantize_point(const std::vector<double>& point,
                                     const std::vector<QuantizationParams>& q_params) {
    size_t D = point.size();
    std::vector<uint64_t> q_point(D);

    for (size_t i = 0; i < D; ++i) {
        int64_t signed_quantized_index = static_cast<int64_t>(round(point[i] / q_params[i].delta));
        q_point[i] = zigzag_encode(signed_quantized_index);
    }

    return q_point;
}

inline std::vector<double> dequantize_point(const std::vector<uint64_t>& q_point,
                                     const std::vector<QuantizationParams>& q_params) {
    size_t D = q_point.size();
    std::vector<double> point(D);

    for (size_t i = 0; i < D; ++i) {
        int64_t signed_quantized_index = zigzag_decode(q_point[i]);
        point[i] = static_cast<double>(signed_quantized_index) * q_params[i].delta;
    }

    return point;
}

#endif //QUANTIZER_H
