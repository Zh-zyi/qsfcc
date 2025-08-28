#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <array>

// 代表一个D维数据点
template <size_t D>
using Point = std::array<double, D>;

// 代表一个量化后的D维整数格点
template <size_t D>
using QuantizedPoint = std::array<int64_t, D>;

// 每个维度的量化参数
struct QuantizationParams {
    double min_val;
    double delta; // 步长, 通常是 2 * epsilon
};

// 空间填充曲线类型枚举
enum class SFCType : uint8_t {
    HILBERT,
    Z_ORDER
};

// 压缩文件的头部信息
struct FileHeader {
    uint32_t version = 1;
    uint32_t dimensions;
    uint32_t block_size;
    uint32_t sfc_precision; // SFC的精度/层级
    uint64_t total_points; // 总点数，用于解压时确定块数量
    SFCType sfc_type;
};

// --- 希尔伯特曲线函数 ---
uint64_t map_to_hilbert_2d(uint32_t x, uint32_t y, uint32_t precision) {
    uint64_t d = 0;
    for (int s = precision - 1; s >= 0; s--) {
        uint32_t rx = (x >> s) & 1;
        uint32_t ry = (y >> s) & 1;
        d += (uint64_t)(3 * rx) ^ ry;
        d <<= 2;
        if (ry == 0) {
            if (rx == 1) {
                x = (1u << precision) - 1 - x;
                y = (1u << precision) - 1 - y;
            }
            std::swap(x, y);
        }
    }
    return d >> 2;
}

void map_from_hilbert_2d(uint64_t d, uint32_t precision, uint32_t& x, uint32_t& y) {
    x = y = 0;
    for (int s = 0; s < precision; ++s) {
        uint32_t rx = (d >> 1) & 1;
        uint32_t ry = d & 1;
        if (ry == 0) {
            if (rx == 1) {
                x = (1u << s) - 1 - x;
                y = (1u << s) - 1 - y;
            }
            std::swap(x, y);
        }
        x += rx << s;
        y += ry << s;
        d >>= 2;
    }
}

// --- Z曲线 ---
uint64_t map_to_z_curve_2d(uint32_t x, uint32_t y) {
    uint64_t z = 0;
    for (int i = 0; i < 32; ++i) { // 假设坐标是32位无符号整数
        z |= ((uint64_t)(x & (1u << i))) << i;
        z |= ((uint64_t)(y & (1u << i))) << (i + 1);
    }
    return z;
}

// 从1D Z阶索引恢复2D坐标
void map_from_z_curve_2d(uint64_t z, uint32_t& x, uint32_t& y) {
    x = 0;
    y = 0;
    for (int i = 0; i < 32; ++i) {
        x |= ((z >> (2 * i)) & 1) << i;
        y |= ((z >> (2 * i + 1)) & 1) << i;
    }
}

// --- 占位符：模拟真实的无损压缩库 ---
// 在实际应用中，这里应该调用像 Zstandard, zlib, lzma 等库
std::vector<char> lossless_compress_placeholder(const std::vector<int64_t>& data) {
    // 简单的序列化：直接拷贝内存
    std::vector<char> compressed_data(data.size() * sizeof(int64_t));
    memcpy(compressed_data.data(), data.data(), compressed_data.size());
    return compressed_data;
}

std::vector<int64_t> lossless_decompress_placeholder(const std::vector<char>& compressed_data) {
    std::vector<int64_t> data(compressed_data.size() / sizeof(int64_t));
    memcpy(data.data(), compressed_data.data(), compressed_data.size());
    return data;
}


template <size_t D>
class TimeSeriesCompressor {
public:
    TimeSeriesCompressor(uint32_t block_size, uint32_t sfc_precision, std::vector<QuantizationParams> q_params)
        : block_size_(block_size), sfc_precision_(sfc_precision), q_params_(std::move(q_params)) {
        if (q_params_.size() != D) {
            throw std::invalid_argument("Quantization parameters must match data dimensions.");
        }
        if (D != 2) {
            throw std::invalid_argument("This example implementation only supports 2D data for SFC.");
        }
    }

    bool compress(const std::vector<Point<D>>& data, const std::string& output_filename) {
        // --- 步骤 1 & 2: 量化并映射到SFC ---
        std::vector<uint64_t> sfc_sequence;
        sfc_sequence.reserve(data.size());

        for (const auto& point : data) {
            QuantizedPoint<D> q_point;
            for (size_t i = 0; i < D; ++i) {
                q_point[i] = static_cast<int64_t>(floor((point[i] - q_params_[i].min_val) / q_params_[i].delta));
            }
            sfc_sequence.push_back(map_to_z_curve_2d(q_point[0], q_point[1]));
        }

        // --- 步骤 3: 分块与多阶段压缩 ---
        std::ofstream out_file(output_filename, std::ios::binary);
        if (!out_file) {
            std::cerr << "错误: 无法打开输出文件 " << output_filename << std::endl;
            return false;
        }

        // 写入文件头
        FileHeader header = {1, D, block_size_, sfc_precision_, static_cast<uint64_t>(data.size())};
        out_file.write(reinterpret_cast<const char*>(&header), sizeof(FileHeader));

        // 写入量化参数
        out_file.write(reinterpret_cast<const char*>(q_params_.data()), q_params_.size() * sizeof(QuantizationParams));
        
        // 块索引
        uint64_t block_index_offset = out_file.tellp();
        size_t num_blocks = (sfc_sequence.size() + block_size_ - 1) / block_size_;
        std::vector<uint64_t> block_offsets(num_blocks);
        out_file.seekp(sizeof(uint64_t) * num_blocks, std::ios::cur); // 为索引表预留空间

        uint64_t data_body_offset = out_file.tellp();

        // 逐块压缩
        for (size_t i = 0; i < num_blocks; ++i) {
            block_offsets[i] = (uint64_t)out_file.tellp() - data_body_offset;
            
            size_t start_idx = i * block_size_;
            size_t end_idx = std::min(start_idx + block_size_, sfc_sequence.size());

            std::vector<int64_t> block_data_to_compress;
            block_data_to_compress.reserve(end_idx - start_idx);

            // 1. 存储检查点
            block_data_to_compress.push_back(sfc_sequence[start_idx]);

            // 2. Delta编码
            for (size_t j = start_idx + 1; j < end_idx; ++j) {
                block_data_to_compress.push_back(sfc_sequence[j] - sfc_sequence[j - 1]);
            }
            
            // 3. (可选的RLE) - 此处省略
            
            // 4. 最终无损压缩
            auto compressed_block = lossless_compress_placeholder(block_data_to_compress);
            
            // 写入块大小和数据
            uint64_t compressed_size = compressed_block.size();
            out_file.write(reinterpret_cast<const char*>(&compressed_size), sizeof(compressed_size));
            out_file.write(compressed_block.data(), compressed_block.size());
        }

        // 回填块索引表
        out_file.seekp(block_index_offset);
        out_file.write(reinterpret_cast<const char*>(block_offsets.data()), block_offsets.size() * sizeof(uint64_t));

        out_file.close();
        return true;
    }

    bool query(const std::string& input_filename, size_t timestamp, Point<D>& result) {
        std::ifstream in_file(input_filename, std::ios::binary);
        if (!in_file) {
            std::cerr << "错误: 无法打开输入文件 " << input_filename << std::endl;
            return false;
        }

        // 读取文件头和量化参数
        FileHeader header;
        in_file.read(reinterpret_cast<char*>(&header), sizeof(FileHeader));
        
        if (timestamp >= header.total_points) {
             std::cerr << "错误: 时间戳超出范围" << std::endl;
            return false;
        }

        std::vector<QuantizationParams> q_params(header.dimensions);
        in_file.read(reinterpret_cast<char*>(q_params.data()), q_params.size() * sizeof(QuantizationParams));
        
        // 读取块索引
        size_t num_blocks = (header.total_points + header.block_size - 1) / header.block_size;
        std::vector<uint64_t> block_offsets(num_blocks);
        in_file.read(reinterpret_cast<char*>(block_offsets.data()), block_offsets.size() * sizeof(uint64_t));
        
        uint64_t data_body_offset = in_file.tellg();

        // 1. 定位块
        size_t block_idx = timestamp / header.block_size;
        size_t idx_in_block = timestamp % header.block_size;

        if (block_idx >= block_offsets.size()) {
            std::cerr << "错误: 块索引计算错误" << std::endl;
            return false;
        }

        in_file.seekg(data_body_offset + block_offsets[block_idx]);
        
        // 2. 解压块
        uint64_t compressed_size;
        in_file.read(reinterpret_cast<char*>(&compressed_size), sizeof(compressed_size));
        std::vector<char> compressed_block(compressed_size);
        in_file.read(compressed_block.data(), compressed_size);
        
        auto decompressed_data = lossless_decompress_placeholder(compressed_block);

        // 3. 重建SFC值
        uint64_t sfc_value = decompressed_data[0]; // 检查点
        for (size_t i = 1; i <= idx_in_block; ++i) {
            sfc_value += decompressed_data[i];
        }

        // 4. 逆向映射
        QuantizedPoint<D> q_point;
        uint32_t qx, qy;
        map_from_z_curve_2d(sfc_value, qx, qy);
        q_point[0] = qx;
        q_point[1] = qy;

        std::cout << "逆映射后: (" << qx << ", " << qy << ")" << std::endl;
        
        // 5. (可选) 逆量化
        for (size_t i = 0; i < D; ++i) {
            // 使用区间中点作为重建值
            result[i] = q_params[i].min_val + (q_point[i] + 0.5) * q_params[i].delta;
        }

        return true;
    }

private:
    uint32_t block_size_;
    uint32_t sfc_precision_;
    std::vector<QuantizationParams> q_params_;
};

int main() {
    const size_t D = 2; // 2维数据
    const size_t NUM_POINTS = 5000;
    const uint32_t BLOCK_SIZE = 1024;
    const uint32_t SFC_PRECISION = 16; // 16位精度，每个维度可表示 2^16 个级别

    // 1. 生成模拟数据 (一个在2D平面上缓慢移动的螺旋线)
    std::vector<Point<D>> original_data;
    original_data.reserve(NUM_POINTS);
    for (size_t i = 0; i < NUM_POINTS; ++i) {
        double angle = 0.1 * i;
        double radius = 100 + 0.5 * i;
        original_data.push_back({
            500 + radius * cos(angle),
            500 + radius * sin(angle)
        });
    }

    // 2. 定义量化参数
    double epsilon_x = 0.5;
    double epsilon_y = 0.5;
    std::vector<QuantizationParams> q_params = {
        {-1000.0, 2.0 * epsilon_x}, // min_val_x, delta_x
        {-1000.0, 2.0 * epsilon_y}  // min_val_y, delta_y
    };

    // 3. 压缩
    TimeSeriesCompressor<D> compressor(BLOCK_SIZE, SFC_PRECISION, q_params);
    std::string filename = "compressed_timeseries.bin";
    std::cout << "正在压缩 " << NUM_POINTS << " 个数据点到 " << filename << " ..." << std::endl;
    if (compressor.compress(original_data, filename)) {
        std::cout << "压缩成功!" << std::endl;

        // --- 新增：计算并显示压缩率 ---
        size_t original_size = original_data.size() * sizeof(Point<D>);
        
        std::ifstream compressed_file(filename, std::ios::binary | std::ios::ate);
        if(compressed_file.is_open()){
            size_t compressed_size = compressed_file.tellg();
            compressed_file.close();

            if (compressed_size > 0) {
                double compression_ratio = static_cast<double>(original_size) / compressed_size;
                std::cout << "------------------------------------" << std::endl;
                std::cout << "原始数据大小: " << original_size << " 字节" << std::endl;
                std::cout << "压缩后文件大小: " << compressed_size << " 字节" << std::endl;
                std::cout << "压缩率: " << compression_ratio << " : 1" << std::endl;
                std::cout << "------------------------------------" << std::endl;
            } else {
                 std::cerr << "错误: 压缩文件大小为0。" << std::endl;
            }
        } else {
            std::cerr << "错误: 无法打开压缩文件以计算大小。" << std::endl;
        }
        // --- 结束新增代码 ---

    } else {
        std::cerr << "压缩失败!" << std::endl;
        return 1;
    }

    // 4. 查询
    size_t query_timestamp = 2500;
    Point<D> queried_point;
    std::cout << "\n正在查询时间戳 " << query_timestamp << " 的数据..." << std::endl;
    if (compressor.query(filename, query_timestamp, queried_point)) {
        const auto& original_point = original_data[query_timestamp];
        std::cout << "查询成功!" << std::endl;
        std::cout << "原始数据: (" << original_point[0] << ", " << original_point[1] << ")" << std::endl;
        std::cout << "量化数据：("<< static_cast<int64_t>(floor((original_point[0] - q_params[0].min_val) / q_params[0].delta)) <<
         ", " << static_cast<int64_t>(floor((original_point[1] - q_params[1].min_val) / q_params[1].delta)) <<")" <<std::endl;
        std::cout << "重建数据: (" << queried_point[0] << ", " << queried_point[1] << ")" << std::endl;
        
        double error_x = std::abs(original_point[0] - queried_point[0]);
        double error_y = std::abs(original_point[1] - queried_point[1]);
        std::cout << "重建误差: (" << error_x << ", " << error_y << ")" << std::endl;
        std::cout << "允许误差: (" << epsilon_x << ", " << epsilon_y << ")" << std::endl;
    } else {
        std::cerr << "查询失败!" << std::endl;
        return 1;
    }

    return 0;
}
