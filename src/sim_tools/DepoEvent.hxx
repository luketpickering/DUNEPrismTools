#ifndef DEPOEVENT_HXX_SEEN
#define DEPOEVENT_HXX_SEEN

#include "DepoParticle.hxx"

#include "TVector3.h"
#include "TObjString.h"

#include <vector>
#include <map>
#include <string>

struct PrimaryParticle {
  PrimaryParticle();

  bool IsFinalState;
  int PDG;
  double EKin;
  double EMass;
  TVector3 ThreeMom;

  std::string ToString();
};

struct DepoEvent {
  TVector3 VertexPosition;
  TObjString *RooTrackerInteractionCode;
  Int_t ev_id;

  std::vector<PrimaryParticle> PrimaryParticles;

  void PrintGENIEPassthrough();

  template <size_t N>
  PrimaryParticle GetFirstPrimaryWithPDG(int (&pdg)[N], bool FinalState) {
    for (auto p : PrimaryParticles) {
      for (size_t p_it = 0; p_it < N; ++p_it) {
        if (p.IsFinalState != FinalState) {
          continue;
        }
        if (p.PDG == pdg[p_it]) {
          return p;
        }
      }
    }
    return PrimaryParticle();
  }

  template <size_t N>
  size_t CountPrimaryWithPDG(int (&pdg)[N]) {
    size_t Count = 0;
    for (auto p : PrimaryParticles) {
      for (size_t p_it = 0; p_it < N; ++p_it) {
        if (p.PDG == pdg[p_it]) {
          Count++;
        }
      }
    }
    return Count;
  }

  std::vector<DepoTracked> TrackedDeposits;
  std::vector<DepoParticle> TotalDeposits;

  std::map<size_t, DepoParticle *> RollupPrimaryParticle;
  std::map<size_t, DepoTracked *> TrackedParticleMap;
  std::map<size_t, bool> IsPrimary;

  DepoParticle *GetPrimaryParticle(size_t parent_id, size_t TrackID);

  bool GetIsPrimary(size_t TrackID);

  DepoTracked *GetTrackedDeposit(size_t TrackID);
};

#endif
