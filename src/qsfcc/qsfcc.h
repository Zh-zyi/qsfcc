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
        10
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
    TimeSeriesCompressor(size_t dim_count, uint8_t block_size, uint8_t sfc_precision,
                         std::vector<QuantizationParams> q_params,
                         SFCType sfc_type)
        : D(dim_count), block_size_(block_size), sfc_precision_(sfc_precision),
          q_params_(std::move(q_params)),
          sfc_type_(sfc_type) {
    }

    ssize_t qsfcc_compress(const std::vector<std::vector<double> > &input, const std::string &output_filepath);

    // bool qsfcc_compress(const std::vector<mdPoint<D> > &data, const std::string &output_filename) {
    //     std::vector<uint64_t> sfc_sequence;
    //     sfc_sequence.reserve(data.size());
    //
    //     for (const auto &point: data) {
    //         QuantizedPoint<D> q_point = quantize_point(point, q_params_);
    //
    //         if (q_point[0] < 0 || q_point[1] < 0 ||
    //             q_point[0] > std::numeric_limits<uint32_t>::max() ||
    //             q_point[1] > std::numeric_limits<uint32_t>::max()) {
    //             std::cerr << "错误: 量化后的索引超出uint32_t范围。" << std::endl;
    //             return false;
    //         }
    //
    //         auto qx = static_cast<uint32_t>(q_point[0]);
    //         auto qy = static_cast<uint32_t>(q_point[1]);
    //
    //         switch (sfc_type_) {
    //             case SFCType::HILBERT:
    //                 sfc_sequence.push_back(map_to_hilbert_2d(qx, qy, sfc_precision_));
    //                 break;
    //             case SFCType::Z_ORDER:
    //                 sfc_sequence.push_back(map_to_z_curve_2d(qx, qy));
    //                 break;
    //         }
    //     }
    //
    //     std::ofstream out_file(output_filename, std::ios::binary);
    //     if (!out_file) {
    //         std::cerr << "错误: 无法打开输出文件 " << output_filename << std::endl;
    //         return false;
    //     }
    //
    //     FileHeader header = {D, block_size_, sfc_precision_, static_cast<uint64_t>(data.size()), sfc_type_};
    //     out_file.write(reinterpret_cast<const char *>(&header), sizeof(FileHeader));
    //     out_file.write(reinterpret_cast<const char *>(q_params_.data()), q_params_.size() * sizeof(QuantizationParams));
    //
    //     uint64_t block_index_offset = out_file.tellp();
    //     size_t num_blocks = (sfc_sequence.size() + block_size_ - 1) / block_size_;
    //     std::vector<uint64_t> block_offsets(num_blocks);
    //     out_file.seekp(sizeof(uint64_t) * num_blocks, std::ios::cur);
    //
    //     uint64_t data_body_offset = out_file.tellp();
    //
    //     for (size_t i = 0; i < num_blocks; ++i) {
    //         block_offsets[i] = (uint64_t) out_file.tellp() - data_body_offset;
    //         size_t start_idx = i * block_size_;
    //         size_t end_idx = std::min(start_idx + block_size_, sfc_sequence.size());
    //         std::vector<int64_t> block_data_to_compress;
    //         block_data_to_compress.reserve(end_idx - start_idx);
    //         block_data_to_compress.push_back(sfc_sequence[start_idx]);
    //         for (size_t j = start_idx + 1; j < end_idx; ++j) {
    //             block_data_to_compress.push_back(sfc_sequence[j] - sfc_sequence[j - 1]);
    //         }
    //         auto compressed_block = lossless_compress_zstd(block_data_to_compress);
    //         uint64_t compressed_size = compressed_block.size();
    //         out_file.write(reinterpret_cast<const char *>(&compressed_size), sizeof(compressed_size));
    //         out_file.write(compressed_block.data(), compressed_block.size());
    //     }
    //
    //     out_file.seekp(block_index_offset);
    //     out_file.write(reinterpret_cast<const char *>(block_offsets.data()), block_offsets.size() * sizeof(uint64_t));
    //     out_file.close();
    //     return true;
    // }

    // bool qsfcc_query(const std::string &input_filename, size_t timestamp, Point<D> &result) {
    //     std::ifstream in_file(input_filename, std::ios::binary);
    //     if (!in_file) {
    //         std::cerr << "错误: 无法打开输入文件 " << input_filename << std::endl;
    //         return false;
    //     }
    //
    //     FileHeader header;
    //     in_file.read(reinterpret_cast<char *>(&header), sizeof(FileHeader));
    //     if (timestamp >= header.total_points) {
    //         std::cerr << "错误: 时间戳超出范围" << std::endl;
    //         return false;
    //     }
    //
    //     std::vector<QuantizationParams> q_params(header.dimensions);
    //     in_file.read(reinterpret_cast<char *>(q_params.data()), q_params.size() * sizeof(QuantizationParams));
    //
    //     size_t num_blocks = (header.total_points + header.block_size - 1) / header.block_size;
    //     std::vector<uint64_t> block_offsets(num_blocks);
    //     in_file.read(reinterpret_cast<char *>(block_offsets.data()), block_offsets.size() * sizeof(uint64_t));
    //
    //     uint64_t data_body_offset = in_file.tellg();
    //
    //     size_t block_idx = timestamp / header.block_size;
    //     size_t idx_in_block = timestamp % header.block_size;
    //     if (block_idx >= block_offsets.size()) {
    //         std::cerr << "错误: 块索引计算错误" << std::endl;
    //         return false;
    //     }
    //
    //     in_file.seekg(data_body_offset + block_offsets[block_idx]);
    //
    //     uint64_t compressed_size;
    //     in_file.read(reinterpret_cast<char *>(&compressed_size), sizeof(compressed_size));
    //     std::vector<char> compressed_block(compressed_size);
    //     in_file.read(compressed_block.data(), compressed_size);
    //
    //     auto decompressed_data = lossless_decompress_zstd(compressed_block);
    //
    //     uint64_t sfc_value = decompressed_data[0];
    //     for (size_t i = 1; i <= idx_in_block; ++i) {
    //         sfc_value += decompressed_data[i];
    //     }
    //
    //     QuantizedPoint<D> q_point;
    //     uint32_t qx, qy;
    //     switch (header.sfc_type) {
    //         case SFCType::HILBERT:
    //             map_from_hilbert_2d(sfc_value, header.sfc_precision, qx, qy);
    //             break;
    //         case SFCType::Z_ORDER:
    //             map_from_z_curve_2d(sfc_value, qx, qy);
    //             break;
    //     }
    //     q_point[0] = qx;
    //     q_point[1] = qy;
    //
    //     result = dequantize_point(q_point, q_params);
    //     return true;
    // }

private:
    uint8_t D;
    uint16_t block_size_;
    uint8_t sfc_precision_;
    std::vector<QuantizationParams> q_params_;
    SFCType sfc_type_;
};

#endif // QSFCC_H
