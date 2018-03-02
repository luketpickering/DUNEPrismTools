#include "TH3D.h"
#include "TObjString.h"
#include "TVector3.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>

struct DepoParticle {
  int PDG;
  size_t TrackID;

  TH3D *Deposits;
  TH3D *DaughterDeposits;

  TH3D *Deposits_timesep;
  TH3D *DaughterDeposits_timesep;

  bool TrackTime;
  TH3D *Deposits_ChrgWSumTime;
  TH3D *DaughterDeposits_ChrgWSumTime;

  double timesep_us;

  DepoParticle() {
    PDG = 0;
    TrackID = 0;
    Deposits = nullptr;
    DaughterDeposits = nullptr;
    Deposits_timesep = nullptr;
    DaughterDeposits_timesep = nullptr;
    TrackTime = false;
    Deposits_ChrgWSumTime = nullptr;
    DaughterDeposits_ChrgWSumTime = nullptr;

    timesep_us = 0xdeadbeef;
  }

  void SetTrackTime() {
    TrackTime = true;

    if (!Deposits_ChrgWSumTime && Deposits) {
      Deposits_ChrgWSumTime = static_cast<TH3D *>(Deposits->Clone());
      Deposits_ChrgWSumTime->SetDirectory(nullptr);
    }
    if (!DaughterDeposits_ChrgWSumTime && DaughterDeposits) {
      DaughterDeposits_ChrgWSumTime =
          static_cast<TH3D *>(DaughterDeposits->Clone());
      DaughterDeposits_ChrgWSumTime->SetDirectory(nullptr);
    }
  }

  DepoParticle(int PDG, size_t trackID, double timesep_us = 0xdeadbeef) {
    this->PDG = PDG;
    TrackID = trackID;
    this->timesep_us = timesep_us;

    Deposits = nullptr;
    DaughterDeposits = nullptr;

    Deposits_timesep = nullptr;
    DaughterDeposits_timesep = nullptr;

    TrackTime = false;
    Deposits_ChrgWSumTime = nullptr;
    DaughterDeposits_ChrgWSumTime = nullptr;
  }

  DepoParticle(DepoParticle const &) = delete;
  DepoParticle &operator=(DepoParticle const &) = delete;

  void Swap(DepoParticle &&other) {
    PDG = other.PDG;
    TrackID = other.TrackID;
    timesep_us = other.timesep_us;

    Deposits = other.Deposits;
    DaughterDeposits = other.DaughterDeposits;

    other.Deposits = nullptr;
    other.DaughterDeposits = nullptr;

    Deposits_timesep = other.Deposits_timesep;
    DaughterDeposits_timesep = other.DaughterDeposits_timesep;

    other.Deposits_timesep = nullptr;
    other.DaughterDeposits_timesep = nullptr;

    TrackTime = other.TrackTime;

    Deposits_ChrgWSumTime = other.Deposits_ChrgWSumTime;
    DaughterDeposits_ChrgWSumTime = other.DaughterDeposits_ChrgWSumTime;

    other.Deposits_ChrgWSumTime = nullptr;
    other.DaughterDeposits_ChrgWSumTime = nullptr;
  }

  DepoParticle(DepoParticle &&other) { Swap(std::move(other)); }

  DepoParticle &operator=(DepoParticle &&other) {
    Swap(std::move(other));
    return (*this);
  }

  void LendDepositMaps(TH3D *map1, TH3D *map2) {
    Deposits = map1;
    DaughterDeposits = map2;

    // Clear any deposits from before
    Deposits->Reset();
    DaughterDeposits->Reset();

    if (timesep_us != 0xdeadbeef) {
      Deposits_timesep = static_cast<TH3D *>(map1->Clone());
      Deposits_timesep->SetDirectory(nullptr);
      DaughterDeposits_timesep = static_cast<TH3D *>(map2->Clone());
      DaughterDeposits_timesep->SetDirectory(nullptr);
    }

    if (TrackTime) {
      Deposits_ChrgWSumTime = static_cast<TH3D *>(map1->Clone());
      Deposits_ChrgWSumTime->SetDirectory(nullptr);
      DaughterDeposits_ChrgWSumTime = static_cast<TH3D *>(map2->Clone());
      DaughterDeposits_ChrgWSumTime->SetDirectory(nullptr);
    }
  }

  void DisownDepositMaps() {
    Deposits = nullptr;
    DaughterDeposits = nullptr;
  }
  void DeleteDepositMaps() {
    delete Deposits;
    delete DaughterDeposits;
  }

  void AddDeposit(double *x, double edep, bool IsPrimary) {
    TH3D *dephist = IsPrimary ? Deposits : DaughterDeposits;

    if (timesep_us != 0xdeadbeef) {
      if ((x[3] * 1E-3) > timesep_us) {
        dephist = IsPrimary ? Deposits_timesep : DaughterDeposits_timesep;

        // std::cout << "[INFO]: Filling other " << (x[3] * 1E-3) << std::endl;
      }
    }

    dephist->Fill(x[0], x[1], x[2], edep);

    if (TrackTime) {
      TH3D *timehist =
          IsPrimary ? Deposits_ChrgWSumTime : DaughterDeposits_ChrgWSumTime;

      timehist->Fill(x[0], x[1], x[2], x[3] * edep);
    }
  }

  void AddDeposit(DepoParticle &other) {
    Deposits->Add(other.Deposits);
    DaughterDeposits->Add(other.DaughterDeposits);

    if ((timesep_us != 0xdeadbeef) && (other.timesep_us != 0xdeadbeef)) {
      Deposits_timesep->Add(other.Deposits_timesep);
      DaughterDeposits_timesep->Add(other.DaughterDeposits_timesep);
    }

    if (TrackTime && other.TrackTime) {
      Deposits_ChrgWSumTime->Add(other.Deposits_ChrgWSumTime);
      DaughterDeposits_ChrgWSumTime->Add(other.DaughterDeposits_ChrgWSumTime);
    }
  }

  void Reset() {
    Deposits->Reset();
    DaughterDeposits->Reset();

    if (timesep_us != 0xdeadbeef) {
      Deposits_timesep->Reset();
      DaughterDeposits_timesep->Reset();
    }

    if (TrackTime) {
      Deposits_ChrgWSumTime->Reset();
      DaughterDeposits_ChrgWSumTime->Reset();
    }
  }

  void FinalizeTime() {
    if (TrackTime) {
      Deposits_ChrgWSumTime->Divide(Deposits);
      DaughterDeposits_ChrgWSumTime->Divide(DaughterDeposits);
    }
  }

  virtual ~DepoParticle() {
    DisownDepositMaps();

    if (timesep_us != 0xdeadbeef) {
      delete Deposits_timesep;
      delete DaughterDeposits_timesep;
    }

    if (TrackTime) {
      delete Deposits_ChrgWSumTime;
      delete DaughterDeposits_ChrgWSumTime;
    }
  }
};

struct DepoTracked : public DepoParticle {
  size_t kMaxTrackedSteps;
  double *_Position;
  double *_Momentum;

  double **Position;
  double **Momentum;
  size_t NSteps;

  DepoTracked(int PDG, size_t trackID, int NMaxTrackSteps = 1000,
              double timesep_us = 0xdeadbeef)
      : DepoParticle(PDG, trackID, timesep_us),
        kMaxTrackedSteps(NMaxTrackSteps) {
    _Position = new double[kMaxTrackedSteps * 3];
    _Momentum = new double[kMaxTrackedSteps * 3];

    if (NMaxTrackSteps == 0) {
      throw;
    }

    std::fill_n(_Position, kMaxTrackedSteps * 3, 0xdeadbeef);
    std::fill_n(_Momentum, kMaxTrackedSteps * 3, 0xdeadbeef);

    Position = new double *[kMaxTrackedSteps];
    Momentum = new double *[kMaxTrackedSteps];

    for (size_t st_it = 0; st_it < kMaxTrackedSteps; ++st_it) {
      Position[st_it] = &_Position[st_it * 3];
      Momentum[st_it] = &_Momentum[st_it * 3];
    }
    NSteps = 0;
  }

  DepoTracked(DepoTracked const &) = delete;
  DepoTracked &operator=(DepoTracked const &) = delete;

  void Swap(DepoTracked &&other) {
    kMaxTrackedSteps = other.kMaxTrackedSteps;
    _Position = other._Position;
    _Momentum = other._Momentum;
    Position = other.Position;
    Momentum = other.Momentum;
    NSteps = other.NSteps;

    other._Position = nullptr;
    other._Momentum = nullptr;
    other.Position = nullptr;
    other.Momentum = nullptr;
  }

  DepoTracked(DepoTracked &&other) : DepoParticle(std::move(other)) {
    Swap(std::move(other));
  }

  DepoTracked &operator=(DepoTracked &&other) {
    DepoParticle::Swap(std::move(other));
    Swap(std::move(other));
    return (*this);
  }

  void AddStep(double *x, double *p) {
    if (NSteps == kMaxTrackedSteps) {
      std::cout << "[WARN]: Tried to add step number " << NSteps
                << ", but this exceeds the maximum number of steps ("
                << kMaxTrackedSteps
                << "). Please "
                   "increase DepoTracked::kMaxTrackedSteps and re-run."
                << std::endl;
      return;
    }

    Position[NSteps][0] = x[0];
    Position[NSteps][1] = x[1];
    Position[NSteps][2] = x[2];

    Momentum[NSteps][0] = p[0];
    Momentum[NSteps][1] = p[1];
    Momentum[NSteps][2] = p[2];

    NSteps++;
  }

  virtual ~DepoTracked() {
    delete[] _Position;
    delete[] _Momentum;
    delete[] Position;
    delete[] Momentum;
  }
};

struct PrimaryParticle {
  PrimaryParticle()
      : IsFinalState(false), PDG(0), EKin(0), EMass(0), ThreeMom(0, 0, 0) {}

  bool IsFinalState;
  int PDG;
  double EKin;
  double EMass;
  TVector3 ThreeMom;

  std::string ToString() {
    std::stringstream ss("");

    ss << "PDG = " << PDG << ", E = " << (EKin + EMass) << ", 3Mom = {"
       << ThreeMom.X() << ", " << ThreeMom.Y() << ", " << ThreeMom.Z() << "}"
       << ", IsFinalState: " << IsFinalState << std::flush;
    return ss.str();
  }
};

struct Event {
  TVector3 VertexPosition;
  TObjString *RooTrackerInteractionCode;
  Int_t ev_id;

  std::vector<PrimaryParticle> PrimaryParticles;

  void PrintGENIEPassthrough() {
    std::cout << "[GENIE p/t] Ev id = " << ev_id
              << ", Interaction : " << RooTrackerInteractionCode->GetString()
              << std::endl;
    std::cout << "\tInteraction position: {" << VertexPosition.X() << ", "
              << VertexPosition.Y() << ", " << VertexPosition.Z() << "}"
              << std::endl;
    std::cout << "\tParticle stack: " << std::endl;
    for (auto &p : PrimaryParticles) {
      std::cout << "\t\t" << p.ToString() << std::endl;
    }
  }

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

  DepoParticle *GetPrimaryParticle(size_t parent_id, size_t TrackID) {
    if (RollupPrimaryParticle.count(TrackID)) {
      return RollupPrimaryParticle[TrackID];
    }

    if (RollupPrimaryParticle.count(parent_id)) {
      RollupPrimaryParticle[TrackID] = RollupPrimaryParticle[parent_id];
      IsPrimary[TrackID] = false;
      return GetPrimaryParticle(parent_id, TrackID);
    }

    for (DepoTracked &trk : TrackedDeposits) {
      if (parent_id == trk.TrackID) {
        RollupPrimaryParticle[TrackID] = static_cast<DepoParticle *>(&trk);
        IsPrimary[TrackID] = false;
        return GetPrimaryParticle(parent_id, TrackID);
      }
      if (TrackID == trk.TrackID) {
        RollupPrimaryParticle[TrackID] = static_cast<DepoParticle *>(&trk);
        IsPrimary[TrackID] = true;
        TrackedParticleMap[TrackID] = &trk;
        return GetPrimaryParticle(parent_id, TrackID);
      }
    }
    for (DepoParticle &dep : TotalDeposits) {
      if (parent_id == dep.TrackID) {
        RollupPrimaryParticle[TrackID] = &dep;
        IsPrimary[TrackID] = false;
        return GetPrimaryParticle(parent_id, TrackID);
      }
      if (TrackID == dep.TrackID) {
        RollupPrimaryParticle[TrackID] = static_cast<DepoParticle *>(&dep);
        IsPrimary[TrackID] = true;
        return GetPrimaryParticle(parent_id, TrackID);
      }
    }

    std::cout << "[WARN]: Ev # " << ev_id << " Particle parent: " << parent_id
              << " not found for particle: " << TrackID << std::endl;
    return nullptr;
  }

  bool GetIsPrimary(size_t TrackID) { return IsPrimary[TrackID]; }

  DepoTracked *GetTrackedDeposit(size_t TrackID) {
    if (TrackedParticleMap.count(TrackID)) {
      return TrackedParticleMap[TrackID];
    }
    return nullptr;
  }
};
