#ifndef Z_UTILS_H_INCLUDED
#define Z_UTILS_H_INCLUDED

// Utility functions to create and manipulate z-values.

#include <stdint.h>
#include <vector>

namespace zorder {

// Returns the Z-order value of the cell containing (lat, lon).
uint_fast32_t zvalue_from_coordinates(const double &lon, const double &lat);

// (out_lat, out_lon) makes the center point of the cell with zvalue.
void coordinates_from_zvalue(const uint_fast32_t &zvalue, double &out_lon,
                             double &out_lat);

uint_fast32_t neighbor_north(const uint_fast32_t &zvalue);
uint_fast32_t neighbor_east(const uint_fast32_t &zvalue);
uint_fast32_t neighbor_south(const uint_fast32_t &zvalue);
uint_fast32_t neighbor_west(const uint_fast32_t &zvalue);

// Returns the set of cells which intersects the query window.
std::vector<uint_fast32_t> qw_decomposition(const double &ll_lon,
                                            const double &ll_lat,
                                            const double &ur_lon,
                                            const double &ur_lat);

}  // namespace zorder

#endif
