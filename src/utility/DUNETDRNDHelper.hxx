#ifndef DUNETDRNDHELPER_SEEN
#define DUNETDRNDHELPER_SEEN

#include "CAFReader.hxx"

inline bool FHC_Numu_True_Select(CAFReader const &ev) {
  return ev.isCC && (ev.nuPDG == 14);
}
inline bool RHC_Numu_True_Select(CAFReader const &ev) {
  return ev.isCC && (ev.nuPDG == -14);
}

inline bool ND_FHC_Numu_Select(CAFReader const &ev) {
  return (ev.reco_q == -1) && (ev.muon_exit == 0) && (ev.Ehad_veto < 30);
}
inline bool ND_RHC_Numu_Select(CAFReader const &ev) {
  return (ev.reco_q == 1) && (ev.muon_exit == 0) && (ev.Ehad_veto < 30);
}

inline bool ND_IsCathode_Select(double const &vtx_x) {
  constexpr double half_gap = 2;
  for (int i = -3; i < 4; ++i) {
    double cathode_center = i * 102; // cm
    if ((vtx_x > (cathode_center - half_gap)) &&
        (vtx_x < (cathode_center + half_gap))) {
      return true;
    }
  }
  return false;
}
inline bool ND_IsWall_Select(double const &vtx_x) {
  constexpr double half_gap = 2;
  for (int i = -3; i < 4; ++i) {
    double wall_center = -51 + i * 102; // cm
    if ((vtx_x > (wall_center - half_gap)) &&
        (vtx_x < (wall_center + half_gap))) {
      return true;
    }
  }
  return false;
}

inline bool FD_IsCathode_Select(double const &vtx_x) { return false; }
inline bool FD_IsWall_Select(double const &vtx_x) { return false; }

inline bool PRISMVertex_Select(double const &vtx_x) {
  return std::abs(vtx_x) < 200;
}

inline bool ND_XFV_Select(double const &vtx_x) {
  return PRISMVertex_Select(vtx_x) && !ND_IsWall_Select(vtx_x) &&
         !ND_IsCathode_Select(vtx_x);
}

inline bool FD_XFV_Select(double const &vtx_x) { return std::abs(vtx_x) < 310; }

inline bool FV_Select(CAFReader const &ev) {
  if (ev.isFD) {
    return FD_XFV_Select(ev.vtx_x) && (std::abs(ev.vtx_y) < 550) &&
           (ev.vtx_z > 50) && (ev.vtx_z < 1244);
  } else {
    return ND_XFV_Select(ev.vtx_x) && (std::abs(ev.vtx_y) < 100) &&
           (ev.vtx_z > 50) && (ev.vtx_z < 350);
  }
}

#endif
