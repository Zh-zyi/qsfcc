#include "qsfcc.h"

#include <map>

bool TimeSeriesCompressor::quantize_data(const std::vector<std::vector<double> > &input,
                                            std::vector<std::pair<uint32_t, uint32_t> > &points,
                                            uint32_t &max_coord) {
    for (const auto &coords: input) {
        if (coords.size() != D) {
            std::cerr << "Error: input dim wrong, should be: " << D << " dimensions" << std::endl;
            return false;
        }
        std::vector<int16_t> integer_part(D);
        std::vector<double> fractional_part(D);

        std::vector<double> point = coords;
        std::vector<uint64_t> q_point = quantize_point(point, q_params_);

        // // 划分成整数和小数部分
        // for (size_t i = 0; i < D; i++) {
        //     double int_val;
        //     fractional_part[i] = std::modf(coords[i], &int_val);
        //     integer_part[i] = static_cast<int16_t>(int_val);
        //
        //     if (static_cast<long long>(int_val) % 2 != 0) {
        //         // 如果整数部分是奇数，则反转小数部分
        //         fractional_part[i] = 1.0 - fractional_part[i];
        //     }
        // }

        // std::vector<uint64_t> q_point = quantize_point(fractional_part, q_params_);

        if (q_point[0] < 0 || q_point[1] < 0 ||
            q_point[0] > std::numeric_limits<uint32_t>::max() ||
            q_point[1] > std::numeric_limits<uint32_t>::max()) {
            std::cerr << "Error: Quantified index exceeds the range of uint32_t." << std::endl;
            return false;
        }

        auto qx = static_cast<uint32_t>(q_point[0]);
        auto qy = static_cast<uint32_t>(q_point[1]);

        points.push_back({qx, qy});
        max_coord = std::max({max_coord, qx, qy});
    }

    return true;
}

bool TimeSeriesCompressor::sfc_encode(std::vector<std::pair<uint32_t, uint32_t> > &quantized_points,
                                         uint32_t max_coord,
                                         std::vector<uint64_t> &out_sfc_sequence) {
    out_sfc_sequence.clear();
    out_sfc_sequence.reserve(quantized_points.size());

    uint8_t bits = (max_coord == 0) ? 1 : static_cast<uint8_t>(std::floor(std::log2(max_coord))) + 1;
    if (bits > 32) {
        std::cerr << "Error: Quantified index exceeds the range of uint32_t." << std::endl;
        return false;
    }

    for (const auto &q_pair: quantized_points) {
        switch (sfc_type_) {
            case SFCType::Z_ORDER:
                out_sfc_sequence.push_back(map_to_z_curve_2d(q_pair.first, q_pair.second, bits));
                break;
            case SFCType::HILBERT:
                out_sfc_sequence.push_back(map_to_hilbert_curve_2d(q_pair.first, q_pair.second, bits));
                break;
        }
    }

    return true;
}

ssize_t TimeSeriesCompressor::lossless_compress(const std::vector<uint64_t> &sfc_sequence,
                                                const std::string &output_filepath) {
    std::string output_filepath_data = output_filepath + std::to_string(q_params_[0].delta / 2) + ".bin";
    std::ofstream out_file(output_filepath_data, std::ios::binary);
    if (!out_file) {
        std::cerr << "error: cannot open the file " << output_filepath_data << std::endl;
        return -1;
    }

    // 写入量化参数等元数据
    out_file.write(reinterpret_cast<const char*>(q_params_.data()), q_params_.size() * sizeof(QuantizationParams));

    // 准备数据进行压缩 (注意：Zstd可以处理任何连续内存，这里直接使用uint64_t即可)
    // 如果您的lossless_compress_zstd需要特定类型，请在此处转换。
    // 为了与原始代码保持一致，这里转为int64_t。
    std::vector<int64_t> block_data_to_compress(sfc_sequence.begin(), sfc_sequence.end());

    auto compressed_block = lossless_compress_zstd(block_data_to_compress);

    // 写入压缩数据块
    uint64_t compressed_size = compressed_block.size();
    out_file.write(reinterpret_cast<const char*>(&compressed_size), sizeof(compressed_size));
    out_file.write(compressed_block.data(), compressed_block.size());

    out_file.close();

    // 检查文件是否成功写入并获取大小
    if (!std::filesystem::exists(output_filepath_data)) {
        return -1;
    }
    return static_cast<ssize_t>(std::filesystem::file_size(output_filepath_data));
}

ssize_t TimeSeriesCompressor::qsfcc_compress(const std::vector<std::vector<double> > &input,
                                             const std::string &output_filepath) {
    // 记录中间数据的文件
    // std::string record_filepath = "D:/University/qsfcc/result/quantized_" + std::to_string(q_params_[0].delta / 2) + ".txt";
    // std::ofstream record_file(record_filepath, std::ios::out);
    // if (!record_file) {
    //     std::cerr << "error: cannot open the record_file " << std::endl;
    //     return false;
    // }

    // step 1: quantize data
    std::vector<std::pair<uint32_t, uint32_t> > quantized_points;
    uint32_t max_coord = 0;
    if (!quantize_data(input, quantized_points, max_coord)) {
        return QUANTIZED_ERROR;
    }

    // step 2: sfc encode
    std::vector<uint64_t> sfc_sequence;
    if (!sfc_encode(quantized_points, max_coord, sfc_sequence)) {
        return SFC_ERROR;
    }

    // record_file.close();

    // 步骤 3: 无损编码并写入文件
    ssize_t file_size = lossless_compress(sfc_sequence, output_filepath);
    if (file_size == -1) {
        return LOSSLESS_COMPRESS_ERROR;
    }

    return file_size;
}
