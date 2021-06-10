#include <bitset>

#include "./libmorton/morton.h"
#include "sql/gis/z_utils.h"

namespace zorder {
// NB, modifies lower or upper.
char get_bit(const double c, double &lower, double &upper) {
  char res;
  double middle = (lower + upper) / 2;

  if (c < middle) {
    res = '0';
    upper = middle;
  } else {
    res = '1';
    lower = middle;
  }

  return res;
}

uint_fast32_t zvalue_from_coordinates(const double lon, const double lat) {
  constexpr int fidelity = 32;
  std::string zstring = "";
  double lat_lower = -90;
  double lat_upper = 90;
  double lon_lower = -180;
  double lon_upper = 180;

  for (int i = 0; i < fidelity; i += 2) {
    zstring.push_back(get_bit(lon, lon_lower, lon_upper));
    zstring.push_back(get_bit(lat, lat_lower, lat_upper));
  }

  uint_fast32_t z = std::bitset<32>(zstring).to_ulong();
  return z;
}

void coordinates_from_zvalue(const uint_fast32_t zvalue, double &out_lon,
                             double &out_lat) {
  std::string zstring = std::bitset<32>(zvalue).to_string();
  double lat_lower = -90;
  double lat_upper = 90;
  double lon_lower = -180;
  double lon_upper = 180;
  constexpr char one = '1';

  for (uint i = 0; i < zstring.size(); i += 2) {
    if (zstring[i] == one) {
      lon_lower = (lon_lower + lon_upper) / 2;
    } else {
      lon_upper = (lon_lower + lon_upper) / 2;
    }

    if (zstring[i + 1] == one) {
      lat_lower = (lat_lower + lat_upper) / 2;
    } else {
      lat_upper = (lat_lower + lat_upper) / 2;
    }
  }

  out_lon = (lon_lower + lon_upper) / 2;
  out_lat = (lat_lower + lat_upper) / 2;
}

uint_fast32_t neighbor_north(const uint_fast32_t zvalue) {
  // NB, adds to the opposite coordinate to accomodate the different
  // rotation of the grid (z vs n shaped curve).
  uint_fast16_t lon;
  uint_fast16_t lat;
  libmorton::morton2D_32_decode(zvalue, lon, lat);
  return libmorton::morton2D_32_encode(lon + 1, lat);
}
uint_fast32_t neighbor_east(const uint_fast32_t zvalue) {
  // NB, adds to the opposite coordinate to accomodate the different
  // rotation of the grid (z vs n shaped curve).
  uint_fast16_t lon;
  uint_fast16_t lat;
  libmorton::morton2D_32_decode(zvalue, lon, lat);
  return libmorton::morton2D_32_encode(lon, lat + 1);
}
uint_fast32_t neighbor_south(const uint_fast32_t zvalue) {
  // NB, subtracts from the opposite coordinate to accomodate the different
  // rotation of the grid (z vs n shaped curve).
  uint_fast16_t lon;
  uint_fast16_t lat;
  libmorton::morton2D_32_decode(zvalue, lon, lat);
  return libmorton::morton2D_32_encode(lon - 1, lat);
}
uint_fast32_t neighbor_west(const uint_fast32_t zvalue) {
  // NB, subtracts from the opposite coordinate to accomodate the different
  // rotation of the grid (z vs n shaped curve).
  uint_fast16_t lon;
  uint_fast16_t lat;
  libmorton::morton2D_32_decode(zvalue, lon, lat);
  return libmorton::morton2D_32_encode(lon, lat - 1);
}

std::vector<uint_fast32_t> qw_decomposition(const double ll_lon,
                                            const double ll_lat,
                                            const double ur_lon,
                                            const double ur_lat) {
  std::vector<uint_fast32_t> result;

  uint_fast32_t lower_left = zorder::zvalue_from_coordinates(ll_lon, ll_lat);
  uint_fast32_t upper_right = zorder::zvalue_from_coordinates(ur_lon, ur_lat);

  // Only one cell --> early termination
  if (lower_left == upper_right) {
    result.push_back(lower_left);
    return result;
  }

  uint_fast32_t northern_bound =
      zorder::neighbor_north(zorder::zvalue_from_coordinates(ll_lon, ur_lat));
  uint_fast32_t eastern_bound =
      zorder::neighbor_east(zorder::zvalue_from_coordinates(ur_lon, ll_lat));
  uint_fast32_t eastwards = lower_left;
  uint_fast32_t northwards;

  // Adding cells to result column by column, starting in lower left corner.
  while (eastwards != eastern_bound) {
    northwards = eastwards;
    while (northwards != northern_bound) {
      result.push_back(northwards);
      northwards = zorder::neighbor_north(northwards);
    }
    eastwards = zorder::neighbor_east(eastwards);
    northern_bound = zorder::neighbor_east(northern_bound);
  }

  std::sort(result.begin(), result.end());
  return result;
}

}  // namespace zorder
