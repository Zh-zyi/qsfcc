#include "qsfcc.h"

#include <map>

ssize_t TimeSeriesCompressor::qsfcc_compress(const std::vector<std::vector<double> > &input,
                                             const std::string &output_filepath) {
    std::map<std::vector<int16_t>, uint16_t> integer_part_dict;
    uint16_t next_dict_index = 0;

    // std::vector<uint16_t> dict_indices;
    // dict_indices.reserve(input.size());

    // 记录中间数据的文件
    // std::string record_filepath = "D:/University/qsfcc/result/quantized_" + std::to_string(q_params_[0].delta / 2) + ".txt";
    // std::ofstream record_file(record_filepath, std::ios::out);
    // if (!record_file) {
    //     std::cerr << "error: cannot open the record_file " << std::endl;
    //     return false;
    // }

    uint32_t max_coord = 0;
    std::vector<std::pair<uint32_t, uint32_t> > quantized_points;
    quantized_points.reserve(input.size());

    for (const auto &coords: input) {
        if (coords.size() != D) {
            std::cerr << "Error: input dim wrong, should be: " << D << " dimensions" << std::endl;
            return -1;
        }
        std::vector<int16_t> integer_part(D);
        std::vector<double> fractional_part(D);

        // std::vector<double> point = coords;
        // std::vector<uint64_t> q_point = quantize_point(point, q_params_);

        // 划分成整数和小数部分
        for (size_t i = 0; i < D; i++) {
            double int_val;
            fractional_part[i] = std::modf(coords[i], &int_val);
            integer_part[i] = static_cast<int16_t>(int_val);

            if (static_cast<long long>(int_val) % 2 != 0) {
                // 如果整数部分是奇数，则反转小数部分
                fractional_part[i] = 1.0 - fractional_part[i];
            }
        }

        std::vector<uint64_t> q_point = quantize_point(fractional_part, q_params_);

        // record_file << point[0] << ", " << point[1] << std::endl;

        if (q_point[0] < 0 || q_point[1] < 0 ||
            q_point[0] > std::numeric_limits<uint32_t>::max() ||
            q_point[1] > std::numeric_limits<uint32_t>::max()) {
            std::cerr << "Error: Quantified index exceeds the range of uint32_t." << std::endl;
            return -1;
        }
        // record_file << q_point[0] << ", " << q_point[1] << std::endl;

        auto qx = static_cast<uint32_t>(q_point[0]);
        auto qy = static_cast<uint32_t>(q_point[1]);

        quantized_points.push_back({qx, qy});
        max_coord = std::max({max_coord, qx, qy});
    }
    // 确定每个值的位数
    uint8_t bits = (max_coord == 0) ? 1 : static_cast<uint8_t>(std::floor(std::log2(max_coord))) + 1;
    if (bits > 32) {
        // 安全检查
        std::cerr << "Error: Required bits exceed 32." << std::endl;
        return -1;
    }

    std::vector<uint64_t> sfc_fractional_sequence;
    sfc_fractional_sequence.reserve(input.size());

    for (const auto &q_pair: quantized_points) {
        switch (sfc_type_) {
            case SFCType::Z_ORDER:
                sfc_fractional_sequence.push_back(map_to_z_curve_2d(q_pair.first, q_pair.second, bits));
            // record_file << map_to_z_curve_2d(qx, qy) << std::endl << std::endl;
                break;
            case SFCType::HILBERT: // 假设您已在 SFCType 枚举中添加 HILBERT
                sfc_fractional_sequence.push_back(map_to_hilbert_curve_2d(q_pair.first, q_pair.second, bits));
                break;
        }
    }

    // record_file.close();

    // std::string output_filepath_dict = output_filepath + "dict.bin";
    // std::ofstream dict_file(output_filepath_dict, std::ios::binary);
    // if (!dict_file) {
    //     std::cerr << "error: cannot open the dict_file " << std::endl;
    //     return false;
    // }
    //
    // uint16_t dict_size = integer_part_dict.size();
    // dict_file.write(reinterpret_cast<const char*>(&dict_size), sizeof(uint16_t));
    //
    // for(const auto& pair : integer_part_dict) {
    //     // 写入键：整数向量 (D * 8字节)
    //     dict_file.write(reinterpret_cast<const char*>(pair.first.data()), D * sizeof(int16_t));
    // }
    // dict_file.close();

    std::string output_filepath_data = output_filepath + std::to_string(q_params_[0].delta / 2) + ".bin";
    std::ofstream out_file(output_filepath_data, std::ios::binary);
    if (!out_file) {
        std::cerr << "error: cannot open the data_file " << std::endl;
        return false;
    }

    FileHeader header = {D, block_size_, static_cast<uint64_t>(input.size()), sfc_type_, bits};
    out_file.write(reinterpret_cast<const char *>(&header), sizeof(FileHeader));
    out_file.write(reinterpret_cast<const char *>(q_params_.data()), q_params_.size() * sizeof(QuantizationParams));

    std::vector<int64_t> block_data_to_compress;
    block_data_to_compress.reserve(sfc_fractional_sequence.size());

    block_data_to_compress.push_back(sfc_fractional_sequence[0]);

    for (size_t j = 1; j < sfc_fractional_sequence.size(); ++j) {
        // block_data_to_compress.push_back(sfc_fractional_sequence[j] - sfc_fractional_sequence[j - 1]);
        block_data_to_compress.push_back(sfc_fractional_sequence[j]);
    }

    auto compressed_block = lossless_compress_zstd(block_data_to_compress);

    uint64_t compressed_size = compressed_block.size();
    out_file.write(reinterpret_cast<const char *>(&compressed_size), sizeof(compressed_size));
    out_file.write(compressed_block.data(), compressed_block.size());

    out_file.close();

    ssize_t file_size = static_cast<ssize_t>(
        std::filesystem::file_size(output_filepath_data)
    );
    return file_size;
}
