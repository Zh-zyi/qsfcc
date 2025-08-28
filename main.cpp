#include <fstream>
#include <string>
#include <utility>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <iomanip>
#include <iostream>

#include "src/snappy/snappy.h"
#include "src/machete/machete.h"
#include "src/sz2/sz/include/sz.h"
#include "src/sim_piece/sim_piece.h"
#include "src/qsfcc/qsfcc.h"

// Remember to change this if you run single precision experiment
const static size_t kDoubleSize = 64;
const static std::string kExportExprTablePrefix = "D:/University/qsfcc/src/";
const static std::string kExportExprTableFileName = "perf_table.csv";
const static std::string kDataSetDirPrefix = "D:/University/qsfcc/dataset/";

struct DatasetInfo {
    std::string fileName;
    int dim_count;
};

// const static std::string kDataSetList[] = {
//     "T-Drive-data.csv",
//     "barometric_pressure.csv"
// };
// const static int kDatasetDimensionList[] = {2, 5};
const static DatasetInfo kDataSetList[] = {
    {"T-Drive-data.csv", 2},
    // {"barometric_pressure.csv", 5},
};
const static std::string kMethodList[] = {
    "LZ77", "Zstd", "Snappy", "SZ2", "Machete", "SimPiece", "Deflate", "LZ4", "FPC", "Gorilla", "Chimp128",
    "Elf"
};
const static std::string kMethodList32[] = {
    "LZ77", "Zstd", "Snappy", "SZ2", "Deflate", "LZ4", "Chimp128", "Elf"
};
const static std::string kAbbrList[] = {
    "AP", "AS", "BM", "BT", "BW", "CT", "DT", "IR", "PM10", "SDE", "SUK", "SUSA", "WS"
};
const static std::string kAbbrList32[] = {
    "CT", "DT", "SDE", "SUK", "SUSA"
};
//const static std::string kAbbrList32[] = {
//    "AP", "AS", "BM", "BT", "BW", "CT", "DT", "IR", "PM10", "SDE", "SUK", "SUSA", "WS"
//};

//constexpr static int kBlockSizeList[] = {50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
constexpr static int kBlockSizeList[] = {1000};
//constexpr static double kMaxDiffList[] = {1.0E-1, 1.0E-2, 1.0E-3, 1.0E-4, 1.0E-5, 1.0E-6, 1.0E-7, 1.0E-8};
constexpr static double kMaxDiffList[] = {1.0E-2, 1.0E-3, 1.0E-4};

static uint16_t global_block_size = 0;

std::vector<double> ReadBlock(std::ifstream &file_input_stream_ref, int block_size) {
    std::vector<double> ret;
    ret.reserve(block_size);
    int entry_count = 0;
    double buffer;
    while (!file_input_stream_ref.eof() && entry_count < block_size) {
        file_input_stream_ref >> buffer;
        ret.emplace_back(buffer);
        ++entry_count;
    }
    return ret;
}

std::vector<float> ReadBlock32(std::ifstream &file_input_stream_ref, int block_size) {
    std::vector<float> ret;
    ret.reserve(block_size);
    int entry_count = 0;
    float buffer;
    while (!file_input_stream_ref.eof() && entry_count < block_size) {
        file_input_stream_ref >> buffer;
        ret.emplace_back(buffer);
        ++entry_count;
    }
    return ret;
}

void ResetFileStream(std::ifstream &data_set_input_stream_ref) {
    data_set_input_stream_ref.clear();
    data_set_input_stream_ref.seekg(0, std::ios::beg);
}

static std::string double_to_string_with_precision(double val, size_t precision) {
    std::ostringstream stringBuffer;
    stringBuffer << std::fixed << std::setprecision(precision) << val;
    return stringBuffer.str();
}

class PerfRecord {
private:
    std::chrono::microseconds compression_time_ = std::chrono::microseconds::zero();
    std::chrono::microseconds decompression_time_ = std::chrono::microseconds::zero();
    long compressed_size_in_bits_ = 0;
    int block_count_ = 0;

    std::vector<long> per_dimension_compressed_bits_;
    int dimension_count_ = 0;

public:
    PerfRecord() = default;

    void InitPerDimension(int D) {
        dimension_count_ = D;
        per_dimension_compressed_bits_.assign(D, 0);
    }

    int dimension_count() {
        return dimension_count_;
    }

    void IncreaseCompressionTime(std::chrono::microseconds &duration) {
        compression_time_ += duration;
    }

    auto &compression_time() {
        return compression_time_;
    }

    auto AvgCompressionTimePerBlock() {
        return static_cast<double>(compression_time_.count()) / block_count_;
    }

    void IncreaseDecompressionTime(std::chrono::microseconds &duration) {
        decompression_time_ += duration;
    }

    auto &decompression_time() {
        return decompression_time_;
    }

    auto AvgDecompressionTimePerBlock() {
        return static_cast<double>(decompression_time_.count()) / block_count_;
    }

    long compressed_size_in_bits() {
        return compressed_size_in_bits_;
    }

    long CompressedSizeForDim(int dim) {
        return (dim >= 0 && dim < dimension_count_) ? per_dimension_compressed_bits_[dim] : 0;
    }

    void AddCompressedSize(long size) {
        compressed_size_in_bits_ += size;
    }

    void AddCompressedSizeForDim(int dim, long size) {
        if (dim >= 0 && dim < dimension_count_) {
            per_dimension_compressed_bits_[dim] += size;
        }
    }

    void set_block_count(int blockCount_) {
        block_count_ = blockCount_;
    }

    int block_count() {
        return block_count_;
    }

    double CalCompressionRatio() const {
        return (double) (block_count_ * global_block_size * kDoubleSize * dimension_count_) / (double)
               compressed_size_in_bits_;
    }

    double CalCompressionRatioForDim(int dim) const {
        if (dim < 0 || dim >= dimension_count_) return 0.0;
        return (double) (block_count_ * global_block_size * kDoubleSize) / (double)
               per_dimension_compressed_bits_[dim];
    }

    float CalCompressionRatio_32() {
        return (float) (block_count_ * global_block_size * kDoubleSize) / (float) compressed_size_in_bits_;
    }

    friend std::ostream &operator<<(std::ostream &os, const PerfRecord &record);
};

std::ostream &operator<<(std::ostream &os, const PerfRecord &record) {
    os << std::fixed << std::setprecision(6);

    os << "======== PerfRecord Summary ========\n";
    os << std::left << std::setw(30) << "Block count: " << record.block_count_ << "\n";

    os << std::left << std::setw(30) << "Total comp time: " << record.compression_time_.count() << " us\n";
    os << std::left << std::setw(30) << "Avg comp time per block: "
            << (record.block_count_ ? (double) record.compression_time_.count() / record.block_count_ : 0.0) << " us\n";

    os << std::left << std::setw(30) << "Total decomp time: " << record.decompression_time_.count() << " us\n";
    os << std::left << std::setw(30) << "Avg decomp time per block: "
            << (record.block_count_ ? (double) record.decompression_time_.count() / record.block_count_ : 0.0) <<
            " us\n";

    os << std::left << std::setw(30) << "Original size: " << record.block_count_ * global_block_size * kDoubleSize <<
            " bits\n";
    os << std::left << std::setw(30) << "Total compressed size: " << record.compressed_size_in_bits_ << " bits\n";
    os << std::left << std::setw(30) << "Overall CR: " << record.CalCompressionRatio() << "\n";

    if (!record.per_dimension_compressed_bits_.empty()) {
        os << "Per-dim CR:\n";
        for (int i = 0; i < record.dimension_count_; ++i) {
            os << "  Dim " << i << ": " << record.CalCompressionRatioForDim(i) << "\n";
        }
    }

    os << "=====================================\n\n";
    return os;
}

class ExprConf {
public:
    struct hash {
        std::size_t operator()(const ExprConf &conf) const {
            return std::hash<std::string>()(conf.method_ + conf.data_set_ + conf.max_diff_);
        }
    };

    ExprConf() = delete;

    ExprConf(std::string method, std::string data_set, double max_diff) : method_(method),
                                                                          data_set_(data_set),
                                                                          max_diff_(double_to_string_with_precision(
                                                                              max_diff, 8)) {
    }

    bool operator==(const ExprConf &otherConf) const {
        return method_ == otherConf.method_ && data_set_ == otherConf.data_set_ && max_diff_ == otherConf.max_diff_;
    }

    std::string method() const {
        return method_;
    }

    std::string data_set() const {
        return data_set_;
    }

    std::string max_diff() const {
        return max_diff_;
    }

private:
    const std::string method_;
    const std::string data_set_;
    const std::string max_diff_;
};

std::unordered_map<ExprConf, PerfRecord, ExprConf::hash> expr_table;

void ExportTotalExprTable() {
    std::ofstream expr_table_output_stream(kExportExprTablePrefix + kExportExprTableFileName);
    if (!expr_table_output_stream.is_open()) {
        std::cerr << "Failed to export performance data." << std::endl;
        exit(-1);
    }
    // Write header
    expr_table_output_stream
            << "Method,DataSet,D,MaxDiff,CompressionTime(AvgPerBlock),DecompressionTime(AvgPerBlock),CompressionRatio"
            << std::endl;

    // Write record
    for (const auto &conf_record: expr_table) {
        auto conf = conf_record.first;
        auto record = conf_record.second;
        expr_table_output_stream << conf.method() << "," << conf.data_set() << ","
                << record.dimension_count() << "," << conf.max_diff() << ","
                << record.AvgCompressionTimePerBlock() << ","
                << record.AvgDecompressionTimePerBlock() << ","
                << record.CalCompressionRatio();
        for (int i = 0; i < record.dimension_count(); i++) {
            expr_table_output_stream << "," << record.CalCompressionRatioForDim(i);
        }
        expr_table_output_stream << std::endl;
    }
    // Go!!
    expr_table_output_stream.flush();
    expr_table_output_stream.close();
}

void ExportExprTableWithCompressionRatio() {
    std::ofstream expr_table_output_stream(kExportExprTablePrefix + kExportExprTableFileName);
    if (!expr_table_output_stream.is_open()) {
        std::cerr << "Failed to export performance data." << std::endl;
        exit(-1);
    }
    // Write header
    expr_table_output_stream
            << "Method,DataSet,MaxDiff,CompressionRatio"
            << std::endl;
    // Write record
    for (const auto &conf_record: expr_table) {
        auto conf = conf_record.first;
        auto record = conf_record.second;
        expr_table_output_stream << conf.method() << "," << conf.data_set() << "," << conf.max_diff() << ","
                << record.CalCompressionRatio() << std::endl;
    }
    // Go!!
    expr_table_output_stream.flush();
    expr_table_output_stream.close();
}

void ExportExprTableWithCompressionTime() {
    std::ofstream expr_table_output_stream(kExportExprTablePrefix + kExportExprTableFileName);
    if (!expr_table_output_stream.is_open()) {
        std::cerr << "Failed to export performance data." << std::endl;
        exit(-1);
    }
    // Write header
    expr_table_output_stream
            << "Method,DataSet,MaxDiff,CompressionTime(Total)"
            << std::endl;
    // Write record
    for (const auto &conf_record: expr_table) {
        auto conf = conf_record.first;
        auto record = conf_record.second;
        expr_table_output_stream << conf.method() << "," << conf.data_set() << "," << conf.max_diff() << ","
                << record.compression_time().count() << std::endl;
    }
    // Go!!
    expr_table_output_stream.flush();
    expr_table_output_stream.close();
}

void ExportExprTableWithDecompressionTime() {
    std::ofstream expr_table_output_stream(kExportExprTablePrefix + kExportExprTableFileName);
    if (!expr_table_output_stream.is_open()) {
        std::cerr << "Failed to export performance data." << std::endl;
        exit(-1);
    }
    // Write header
    expr_table_output_stream
            << "Method,DataSet,MaxDiff,DecompressionTime(Total)"
            << std::endl;
    // Write record
    for (const auto &conf_record: expr_table) {
        auto conf = conf_record.first;
        auto record = conf_record.second;
        expr_table_output_stream << conf.method() << "," << conf.data_set() << "," << conf.max_diff() << ","
                << record.decompression_time().count() << std::endl;
    }
    // Go!!
    expr_table_output_stream.flush();
    expr_table_output_stream.close();
}

std::string trim(const std::string& str) {
    const std::string whitespace = " \t\n\r\f\v";
    size_t first = str.find_first_not_of(whitespace);
    if (std::string::npos == first) {
        return ""; // 字符串全是空白
    }
    size_t last = str.find_last_not_of(whitespace);
    return str.substr(first, (last - first + 1));
}

std::vector<std::vector<double> > ReadBlockMultiDim(std::ifstream &ifs, int block_size, int D) {
    std::vector<std::vector<double> > result;
    std::string line;
    int count = 0;

    while (count < block_size && std::getline(ifs, line)) {
        std::vector<double> values;
        std::stringstream ss(line);
        std::string token;
        while (std::getline(ss, token, ',')) {
            std::string trimmed_token = trim(token); // 清理 token
            if (trimmed_token.empty()) {
                continue; // 如果清理后为空，则跳过
            }
            values.push_back(std::stod(token));
        }
        if (values.size() == D) {
            result.push_back(values);
            ++count;
        }
    }
    return result;
}

PerfRecord PerfSnappy(std::ifstream &data_set_input_stream_ref, int block_size, size_t D) {
    PerfRecord perf_record;
    perf_record.InitPerDimension(D);

    int block_count = 0;

    while (true) {
        std::vector<std::vector<double> > multi_dim_data = ReadBlockMultiDim(data_set_input_stream_ref, block_size, D);
        if (multi_dim_data.size() < block_size) break;
        ++block_count;

        // 初始化每个维度的数据
        std::vector<std::vector<double> > dim_data(D, std::vector<double>(block_size));
        for (int i = 0; i < block_size; ++i)
            for (int d = 0; d < D; ++d)
                dim_data[d][i] = multi_dim_data[i][d];

        for (int d = 0; d < D; ++d) {
            std::string compression_output;
            std::string decompression_output;
            auto compression_start_time = std::chrono::steady_clock::now();
            size_t compression_output_len = snappy::Compress(reinterpret_cast<const char *>(dim_data[d].data()),
                                                             dim_data[d].size() * sizeof(double), &compression_output);
            auto compression_end_time = std::chrono::steady_clock::now();

            perf_record.AddCompressedSize(compression_output_len * 8);
            perf_record.AddCompressedSizeForDim(d, compression_output_len * 8);

            auto decompression_start_time = std::chrono::steady_clock::now();
            snappy::Uncompress(compression_output.data(), compression_output.size(), &decompression_output);
            auto decompression_end_time = std::chrono::steady_clock::now();

            auto compression_time_in_a_block = std::chrono::duration_cast<std::chrono::microseconds>(
                compression_end_time - compression_start_time);
            auto decompression_time_in_a_block = std::chrono::duration_cast<std::chrono::microseconds>(
                decompression_end_time - decompression_start_time);

            perf_record.IncreaseCompressionTime(compression_time_in_a_block);
            perf_record.IncreaseDecompressionTime(decompression_time_in_a_block);
        }
    }

    perf_record.set_block_count(block_count);
    std::cout << perf_record;
    return perf_record;
}

PerfRecord PerfMachete(std::ifstream &data_set_input_stream_ref, double max_diff, uint16_t block_size, size_t D) {
    PerfRecord perf_record;
    perf_record.InitPerDimension(D);

    int block_count = 0;

    while (true) {
        std::vector<std::vector<double> > multi_dim_data = ReadBlockMultiDim(data_set_input_stream_ref, block_size, D);
        if (multi_dim_data.size() < block_size) break;
        ++block_count;

        // 初始化每个维度的数据
        std::vector<std::vector<double> > dim_data(D, std::vector<double>(block_size));
        for (int i = 0; i < block_size; ++i)
            for (int d = 0; d < D; ++d)
                dim_data[d][i] = multi_dim_data[i][d];

        for (int d = 0; d < D; ++d) {
            auto *compression_buffer = new uint8_t[100000];
            auto *decompression_buffer = new double[block_size];

            auto compression_start_time = std::chrono::steady_clock::now();
            ssize_t compression_output_len = machete_compress<lorenzo1, hybrid>(
                dim_data[d].data(), block_size, &compression_buffer, max_diff);
            auto compression_end_time = std::chrono::steady_clock::now();

            perf_record.AddCompressedSize(compression_output_len * 8);
            perf_record.AddCompressedSizeForDim(d, compression_output_len * 8);

            auto decompression_start_time = std::chrono::steady_clock::now();
            ssize_t decompression_output_len = machete_decompress<lorenzo1, hybrid>(
                compression_buffer, compression_output_len, decompression_buffer);
            auto decompression_end_time = std::chrono::steady_clock::now();

            auto compression_time = std::chrono::duration_cast<std::chrono::microseconds>(
                compression_end_time - compression_start_time);
            auto decompression_time = std::chrono::duration_cast<std::chrono::microseconds>(
                decompression_end_time - decompression_start_time);

            perf_record.IncreaseCompressionTime(compression_time);
            perf_record.IncreaseDecompressionTime(decompression_time);

            delete[] compression_buffer;
            delete[] decompression_buffer;
        }
    }

    perf_record.set_block_count(block_count);
    std::cout << perf_record;
    return perf_record;
}

PerfRecord PerfSZ2(std::ifstream &data_set_input_stream_ref, double max_diff, int block_size, size_t D) {
    PerfRecord perf_record;
    perf_record.InitPerDimension(D);

    int block_count = 0;
    while (true) {
        std::vector<std::vector<double> > multi_dim_data = ReadBlockMultiDim(data_set_input_stream_ref, block_size, D);
        if (multi_dim_data.size() < block_size) break;
        ++block_count;

        // 初始化每个维度的数据
        std::vector<std::vector<double> > dim_data(D, std::vector<double>(block_size));
        for (int i = 0; i < block_size; ++i)
            for (int d = 0; d < D; ++d)
                dim_data[d][i] = multi_dim_data[i][d];

        for (int d = 0; d < D; ++d) {
            size_t compression_output_len;
            auto decompression_output = new double[block_size];

            auto compression_start_time = std::chrono::steady_clock::now();
            auto compression_output = SZ_compress_args(SZ_DOUBLE, dim_data[d].data(), &compression_output_len,
                                                       ABS, max_diff * 0.99, 0, 0, 0, 0, 0, 0, dim_data[d].size());
            auto compression_end_time = std::chrono::steady_clock::now();

            perf_record.AddCompressedSize(compression_output_len * 8);
            perf_record.AddCompressedSizeForDim(d, compression_output_len * 8);

            auto decompression_start_time = std::chrono::steady_clock::now();
            size_t decompression_output_len = SZ_decompress_args(SZ_DOUBLE, compression_output,
                                                                 compression_output_len, decompression_output, 0, 0,
                                                                 0, 0, block_size);
            auto decompression_end_time = std::chrono::steady_clock::now();

            auto compression_time_in_a_block = std::chrono::duration_cast<std::chrono::microseconds>(
                compression_end_time - compression_start_time);
            auto decompression_time_in_a_block = std::chrono::duration_cast<std::chrono::microseconds>(
                decompression_end_time - decompression_start_time);

            perf_record.IncreaseCompressionTime(compression_time_in_a_block);
            perf_record.IncreaseDecompressionTime(decompression_time_in_a_block);

            delete[] decompression_output;
        }
    }

    perf_record.set_block_count(block_count);
    std::cout << perf_record;
    return perf_record;
}

PerfRecord PerfSimPiece(std::ifstream &data_set_input_stream_ref, double max_diff, int block_size, size_t D) {
    PerfRecord perf_record;
    perf_record.InitPerDimension(D);

    int block_count = 0;
    while (true) {
        std::vector<std::vector<double> > multi_dim_data = ReadBlockMultiDim(data_set_input_stream_ref, block_size, D);
        if (multi_dim_data.size() < block_size) break;
        ++block_count;

        // 初始化每个维度的数据
        std::vector<std::vector<double> > dim_data(D, std::vector<double>(block_size));
        for (int i = 0; i < block_size; ++i)
            for (int d = 0; d < D; ++d)
                dim_data[d][i] = multi_dim_data[i][d];

        for (int d = 0; d < D; ++d) {
            std::vector<Point> input_points;
            for (int i = 0; i < dim_data[d].size(); ++i) {
                input_points.emplace_back(i, dim_data[d][i]);
            }

            char *compression_output = new char [dim_data[d].size() * 8];
            int compression_output_len = 0;
            int timestamp_store_size;
            auto compression_start_time = std::chrono::steady_clock::now();
            SimPiece sim_piece_compress(input_points, max_diff);
            compression_output_len = sim_piece_compress.toByteArray(compression_output, true, &timestamp_store_size);
            auto compression_end_time = std::chrono::steady_clock::now();

            perf_record.AddCompressedSize((compression_output_len - timestamp_store_size) * 8);
            perf_record.AddCompressedSizeForDim(d, (compression_output_len - timestamp_store_size) * 8);

            auto decompression_start_time = std::chrono::steady_clock::now();
            SimPiece sim_piece_decompress(compression_output, compression_output_len, true);
            sim_piece_decompress.decompress();
            auto decompression_end_time = std::chrono::steady_clock::now();

            auto compression_time_in_a_block = std::chrono::duration_cast<std::chrono::microseconds>(
                compression_end_time - compression_start_time);
            auto decompression_time_in_a_block = std::chrono::duration_cast<std::chrono::microseconds>(
                decompression_end_time - decompression_start_time);

            perf_record.IncreaseCompressionTime(compression_time_in_a_block);
            perf_record.IncreaseDecompressionTime(decompression_time_in_a_block);
        }
    }

    perf_record.set_block_count(block_count);
    std::cout << perf_record;
    return perf_record;
}

PerfRecord PerfQSFCC(std::ifstream &data_set_input_stream_ref, double max_diff, uint16_t block_size, size_t D) {
    PerfRecord perf_record;
    // const size_t D = 2;
    perf_record.InitPerDimension(D);

    const uint8_t sfc_precision = 16;
    std::vector<QuantizationParams> q_params = {
        {2.0 * max_diff}, {2.0 * max_diff}
    };
    TimeSeriesCompressor qsfcc_compressor(D, block_size, sfc_precision, q_params, SFCType::Z_ORDER);

    uint32_t block_count = 0;

    while (true) {
        std::vector<std::vector<double> > multi_dim_data = ReadBlockMultiDim(data_set_input_stream_ref, block_size, D);
        if (multi_dim_data.size() < block_size) break;
        ++block_count;

        auto *compression_buffer = new uint8_t[100000];
        auto *decompression_buffer = new double[block_size];
        std::string compress_result_filepath = "D:/University/qsfcc/src/QSFCC/data/compress_result" +
                                               std::to_string(max_diff) + ".bin";

        auto compression_start_time = std::chrono::steady_clock::now();
        ssize_t compression_output_len = qsfcc_compressor.qsfcc_compress(multi_dim_data, compress_result_filepath);
        auto compression_end_time = std::chrono::steady_clock::now();

        perf_record.AddCompressedSize(compression_output_len * 8);

        auto compression_time = std::chrono::duration_cast<std::chrono::microseconds>(
            compression_end_time - compression_start_time);
        perf_record.IncreaseCompressionTime(compression_time);

        delete[] compression_buffer;
        delete[] decompression_buffer;
    }

    perf_record.set_block_count(block_count);
    std::cout << perf_record;
    return perf_record;
}

int main() {
    global_block_size = kBlockSizeList[0];
    for (const auto &data_set: kDataSetList) {
        std::string data_set_name = data_set.fileName;
        size_t D = data_set.dim_count;
        std::ifstream data_set_input_stream(kDataSetDirPrefix + data_set_name);
        if (!data_set_input_stream.is_open()) {
            std::cerr << "Failed to open the file [" << kDataSetDirPrefix + data_set_name << "]" << std::endl;
        }
        std::cout << "-------- dataset: " << data_set_name << "   D: " << D << " --------" << std::endl;
        for (const auto &max_diff: kMaxDiffList) {
            std::cout << "-------- max_diff: " << max_diff << "--------" << std::endl;

            std::cout << "--------- Machete ---------\n";
            expr_table.insert(std::make_pair(ExprConf("Machete", data_set_name, max_diff),
                                             PerfMachete(data_set_input_stream, max_diff, global_block_size, D)));
            ResetFileStream(data_set_input_stream);

            std::cout << "--------- SimPiece ---------\n";
            expr_table.insert(std::make_pair(ExprConf("SimPiece", data_set_name, max_diff),
                                             PerfSimPiece(data_set_input_stream, max_diff, global_block_size, D)));
            ResetFileStream(data_set_input_stream);

            std::cout << "--------- SZ2 ---------\n";
            expr_table.insert(std::make_pair(ExprConf("SZ2", data_set_name, max_diff),
                                             PerfSZ2(data_set_input_stream, max_diff, global_block_size, D)));
            ResetFileStream(data_set_input_stream);

            std::cout << "---------- qsfcc ----------\n";
            expr_table.insert(std::make_pair(ExprConf("qsfcc", data_set_name, max_diff),
                                             PerfQSFCC(data_set_input_stream, max_diff, global_block_size, D)));
            ResetFileStream(data_set_input_stream);
        }

        // lossless
        std::cout << "--------- Snappy ---------\n";
        expr_table.insert(std::make_pair(ExprConf("Snappy", data_set_name, 0),
                                         PerfSnappy(data_set_input_stream, global_block_size, D)));
        ResetFileStream(data_set_input_stream);

        data_set_input_stream.close();
    }

    ExportTotalExprTable();
    return 0;
}
