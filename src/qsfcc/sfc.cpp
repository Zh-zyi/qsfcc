#include "sfc.h"
#include <utility>

// 将2D坐标 (x, y) 转换为 1D希尔伯特曲线值
// x, y: 坐标值，必须小于 (1 << bits)
// bits: 每个坐标的位数，决定了希尔伯特曲线的阶数 (order)
uint64_t map_to_hilbert_curve_2d(uint32_t x, uint32_t y, uint8_t bits) {
    uint64_t h_idx = 0;
    // 从最高位（最大子网格）开始迭代
    for (int s = bits - 1; s >= 0; --s) {
        // 提取当前处理层级的 x 和 y 的位
        uint64_t rx = (x >> s) & 1;
        uint64_t ry = (y >> s) & 1;

        // 当前 2x2 子网格的象限索引 (0, 1, 2, or 3)
        uint64_t quadrant = (ry << 1) | (rx ^ ry);

        // 将象限索引添加到希尔伯特值中
        h_idx = (h_idx << 2) | quadrant;

        // --- 这是希尔伯特曲线的核心：根据象限对坐标进行旋转/翻转 ---
        if (ry == 0) {
            if (rx == 1) {
                // 右上象限的翻转和旋转
                x = (1u << s) - 1 - x;
                y = (1u << s) - 1 - y;
            }
            // 左上象限的旋转 (交换 x 和 y)
            std::swap(x, y);
        }
    }
    return h_idx;
}

// 将 1D希尔伯特曲线值 h 转换回 2D坐标 (x, y)
// h: 希尔伯特曲线值
// bits: 每个坐标的位数
// x, y: 用于接收结果的引用
void map_from_hilbert_curve_2d(uint64_t h, uint32_t& x, uint32_t& y, uint8_t bits) {
    x = 0;
    y = 0;
    // 从最低位（最小子网格）开始迭代
    for (int s = 0; s < bits; ++s) {
        // 提取当前希尔伯特值的最低两位作为象限索引
        uint64_t quadrant = h & 3;
        h >>= 2;

        uint64_t rx = (quadrant == 1 || quadrant == 2);
        uint64_t ry = (quadrant == 2 || quadrant == 3);

        // --- 这是反向变换的核心：根据象限对坐标进行逆向旋转/翻转 ---
        if (ry == 0) {
            if (rx == 1) {
                x = (1u << s) - 1 - x;
                y = (1u << s) - 1 - y;
            }
            // 左上象限的逆向旋转 (交换 x 和 y)
            std::swap(x, y);
        }

        // 将当前计算出的位添加到 x 和 y
        x |= rx << s;
        y |= ry << s;
    }
}

// --- Z阶曲线 ---
uint64_t map_to_z_curve_2d(uint32_t x, uint32_t y, uint8_t bits) {
    uint64_t z = 0;
    for (int i = 0; i < bits; ++i) { // 假设坐标是32位无符号整数
        z |= ((uint64_t)(x & (1u << i))) << i;
        z |= ((uint64_t)(y & (1u << i))) << (i + 1);
    }
    return z;
}

void map_from_z_curve_2d(uint64_t z, uint32_t& x, uint32_t& y, uint8_t bits) {
    x = 0;
    y = 0;
    for (int i = 0; i < bits; ++i) {
        x |= ((z >> (2 * i)) & 1) << i;
        y |= ((z >> (2 * i + 1)) & 1) << i;
    }
}
