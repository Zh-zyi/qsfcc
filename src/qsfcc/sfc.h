#ifndef SFC_H
#define SFC_H

#include <cstdint> // For uint32_t, uint64_t

uint64_t map_to_hilbert_2d(uint32_t x, uint32_t y, uint32_t precision);
void map_from_hilbert_2d(uint64_t d, uint32_t precision, uint32_t& x, uint32_t& y);

uint64_t map_to_z_curve_2d(uint32_t x, uint32_t y);
    void map_from_z_curve_2d(uint64_t z, uint32_t& x, uint32_t& y);

#endif // SFC_H
