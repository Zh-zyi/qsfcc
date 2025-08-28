#include "qsfcc.h"

ssize_t TimeSeriesCompressor::qsfcc_compress(const std::vector<std::vector<double>> &input,
    const std::string &output_filepath)  {
        std::vector<uint64_t> sfc_sequence;
        sfc_sequence.reserve(input.size());

        // std::ofstream record_file("D:/University/myCompression/src/QSFCC/data/record1.txt", std::ios::app);
        // if (!record_file) {
        //     std::cerr << "error: cannot open the record_file " << std::endl;
        //     return false;
        // }
        for (const auto &coords: input) {
            if (coords.size() != D) {
                std::cerr << "Error: input dim wrong, should be: " << D << " dimensions" << std::endl;
                return -1;
            }

            std::vector<double> point = coords;
            std::vector<uint64_t> q_point = quantize_point(point, q_params_);

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

            switch (sfc_type_) {
                case SFCType::HILBERT:
                    sfc_sequence.push_back(map_to_hilbert_2d(qx, qy, sfc_precision_));
                // record_file << map_to_hilbert_2d(qx, qy, sfc_precision_) << std::endl << std::endl;
                    break;
                case SFCType::Z_ORDER:
                    sfc_sequence.push_back(map_to_z_curve_2d(qx, qy));
                // record_file << map_to_z_curve_2d(qx, qy) << std::endl << std::endl;
                    break;
            }
        }
        // record_file.close();

        std::ofstream out_file(output_filepath, std::ios::binary);
        if (!out_file) {
            std::cerr << "error: cannot open the file " << std::endl;
            return false;
        }

        FileHeader header = {D, block_size_, sfc_precision_, static_cast<uint64_t>(input.size()), sfc_type_};
        out_file.write(reinterpret_cast<const char *>(&header), sizeof(FileHeader));
        out_file.write(reinterpret_cast<const char *>(q_params_.data()), q_params_.size() * sizeof(QuantizationParams));

        uint64_t block_index_offset = out_file.tellp();
        size_t num_blocks = (sfc_sequence.size() + block_size_ - 1) / block_size_;
        std::vector<uint64_t> block_offsets(num_blocks);
        out_file.seekp(sizeof(uint64_t) * num_blocks, std::ios::cur);

        uint64_t data_body_offset = out_file.tellp();

        for (size_t i = 0; i < num_blocks; ++i) {
            block_offsets[i] = (uint64_t) out_file.tellp() - data_body_offset;
            size_t start_idx = i * block_size_;
            size_t end_idx = std::min(start_idx + block_size_, sfc_sequence.size());
            std::vector<int64_t> block_data_to_compress;
            block_data_to_compress.reserve(end_idx - start_idx);
            block_data_to_compress.push_back(sfc_sequence[start_idx]);
            for (size_t j = start_idx + 1; j < end_idx; ++j) {
                block_data_to_compress.push_back(sfc_sequence[j] - sfc_sequence[j - 1]);
            }
            auto compressed_block = lossless_compress_zstd(block_data_to_compress);
            uint64_t compressed_size = compressed_block.size();
            out_file.write(reinterpret_cast<const char *>(&compressed_size), sizeof(compressed_size));
            out_file.write(compressed_block.data(), compressed_block.size());
        }

        out_file.seekp(block_index_offset);
        out_file.write(reinterpret_cast<const char *>(block_offsets.data()), block_offsets.size() * sizeof(uint64_t));
        out_file.close();

        ssize_t file_size = static_cast<ssize_t>(
            std::filesystem::file_size(output_filepath)
        );
        return file_size;
    }
