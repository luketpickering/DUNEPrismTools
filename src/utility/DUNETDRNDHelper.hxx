#ifndef DUNETDRNDHELPER_SEEN
#define DUNETDRNDHELPER_SEEN

#include "CAFReader.hxx"

inline bool ND_FHC_Select(CAFReader const &ev) {
  return (ev.reco_q == -1) && (ev.muon_exit == 0) && (ev.Ehad_veto < 30);
}
inline bool ND_RHC_Select(CAFReader const &ev) {
  return (ev.reco_q == 1) && (ev.muon_exit == 0) && (ev.Ehad_veto < 30);
}
inline bool IsCathode_Select(double const &vtx_x) {
  constexpr double half_gap = 1.5;
  for (int i = -3; i < 4; ++i) {
    double cathode_center = i * 102.1; // cm
    if ((vtx_x > (cathode_center - half_gap)) &&
        (vtx_x < (cathode_center + half_gap))) {
      return true;
    }
  }
  return false;
}
inline bool IsWall_Select(double const &vtx_x) {
  constexpr double half_gap = 1.5;
  for (int i = -3; i < 4; ++i) {
    double wall_center = -51.05 + i * 102.1; // cm
    if ((vtx_x > (wall_center - half_gap)) &&
        (vtx_x < (wall_center + half_gap))) {
      return true;
    }
  }
  return false;
}
inline bool XFV_Select(double const &vtx_x) {
  return (std::abs(vtx_x) < 300) && !IsWall_Select(vtx_x) &&
         !IsCathode_Select(vtx_x);
}
inline bool FV_Select(CAFReader const &ev) {
  return XFV_Select(ev.vtx_x) && (std::abs(ev.vtx_y) < 100) &&
         (ev.vtx_z > 50) && (ev.vtx_z < 350);
}

#endif
