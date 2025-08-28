#include "sfc.h"
#include <utility>

// --- 希尔伯特曲线 ---
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

// --- Z阶曲线 ---
uint64_t map_to_z_curve_2d(uint32_t x, uint32_t y) {
    uint64_t z = 0;
    for (int i = 0; i < 32; ++i) { // 假设坐标是32位无符号整数
        z |= ((uint64_t)(x & (1u << i))) << i;
        z |= ((uint64_t)(y & (1u << i))) << (i + 1);
    }
    return z;
}

void map_from_z_curve_2d(uint64_t z, uint32_t& x, uint32_t& y) {
    x = 0;
    y = 0;
    for (int i = 0; i < 32; ++i) {
        x |= ((z >> (2 * i)) & 1) << i;
        y |= ((z >> (2 * i + 1)) & 1) << i;
    }
}
