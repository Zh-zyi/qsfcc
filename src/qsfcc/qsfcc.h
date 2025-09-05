#ifndef QSFCC_H
#define QSFCC_H

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <cstring>
#include <sstream>
#include <filesystem>
#include <zstd.h>

#include "defs.h"
#include "sfc.h"
#include "quantizer.h"

// 量化降维后使用无损压缩，Zstd
inline std::vector<char> lossless_compress_zstd(const std::vector<int64_t> &data) {
    size_t src_size = data.size() * sizeof(int64_t);
    const void *src = reinterpret_cast<const void *>(data.data());

    size_t max_dst_size = ZSTD_compressBound(src_size);
    std::vector<char> compressed_data(max_dst_size);

    size_t compressed_size = ZSTD_compress(
        compressed_data.data(), max_dst_size,
        src, src_size,
        6
    );

    if (ZSTD_isError(compressed_size)) {
        throw std::runtime_error(std::string("ZSTD compression error: ") + ZSTD_getErrorName(compressed_size));
    }

    compressed_data.resize(compressed_size);
    return compressed_data;
}

inline std::vector<int64_t> lossless_decompress_zstd(const std::vector<char> &compressed_data, size_t original_size) {
    // original_size 是压缩前数据的字节数 (== data.size() * sizeof(int64_t))
    std::vector<int64_t> data(original_size / sizeof(int64_t));

    size_t decompressed_size = ZSTD_decompress(
        data.data(), original_size,
        compressed_data.data(), compressed_data.size()
    );

    if (ZSTD_isError(decompressed_size)) {
        throw std::runtime_error(std::string("ZSTD decompression error: ") + ZSTD_getErrorName(decompressed_size));
    }

    return data;
}

class TimeSeriesCompressor {
public:
    TimeSeriesCompressor(size_t dim_count, uint8_t block_size,
                         std::vector<QuantizationParams> q_params,
                         SFCType sfc_type)
        : D(dim_count), block_size_(block_size),
          q_params_(std::move(q_params)),
          sfc_type_(sfc_type) {
    }

    bool quantize_data(const std::vector<std::vector<double> > &input,
                          std::vector<std::pair<uint32_t, uint32_t> > &points,
                          uint32_t &max_coord);

    bool sfc_encode(std::vector<std::pair<uint32_t, uint32_t> > &quantized_points,
                       uint32_t max_coord,
                       std::vector<uint64_t> &out_sfc_sequence);

    ssize_t lossless_compress(const std::vector<uint64_t> &sfc_sequence,
                              const std::string &output_filepath);

    ssize_t qsfcc_compress(const std::vector<std::vector<double> > &input, const std::string &output_filepath_data);

private:
    uint8_t D;
    uint16_t block_size_;
    std::vector<QuantizationParams> q_params_;
    SFCType sfc_type_;
};

#endif // QSFCC_H
