#ifndef SFC_H
#define SFC_H

#include <cstdint> // For uint32_t, uint64_t

uint64_t map_to_hilbert_curve_2d(uint32_t x, uint32_t y, uint8_t bits);

void map_from_hilbert_curve_2d(uint64_t h, uint32_t& x, uint32_t& y, uint8_t bits);

uint64_t map_to_z_curve_2d(uint32_t x, uint32_t y, uint8_t bits);

void map_from_z_curve_2d(uint64_t z, uint32_t &x, uint32_t &y, uint8_t bits);

#endif // SFC_H
