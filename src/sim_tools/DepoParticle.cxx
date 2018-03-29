#include "DepoParticle.hxx"

#include <iostream>

DepoParticle::DepoParticle() {
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

void DepoParticle::SetTrackTime() {
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

DepoParticle::DepoParticle(int PDG, size_t trackID, double timesep_us) {
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
void DepoParticle::Swap(DepoParticle &&other) {
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
DepoParticle::DepoParticle(DepoParticle &&other) { Swap(std::move(other)); }
DepoParticle &DepoParticle::operator=(DepoParticle &&other) {
  Swap(std::move(other));
  return (*this);
}
void DepoParticle::LendDepositMaps(TH3D *map1, TH3D *map2) {
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

void DepoParticle::DisownDepositMaps() {
  Deposits = nullptr;
  DaughterDeposits = nullptr;
}
void DepoParticle::DeleteDepositMaps() {
  delete Deposits;
  delete DaughterDeposits;
}

void DepoParticle::AddDeposit(double *x, double edep, bool IsPrimary) {
  TH3D *dephist = IsPrimary ? Deposits : DaughterDeposits;

  if (timesep_us != 0xdeadbeef) {
    if ((x[3] * 1E-3) > timesep_us) {
      dephist = IsPrimary ? Deposits_timesep : DaughterDeposits_timesep;
    }
  }
#ifdef DEBUG
  std::cout << "[INFO]: Add dep = " << edep << " at time " << (x[3] * 1E-3)
            << std::endl;
#endif
  dephist->Fill(x[0], x[1], x[2], edep);

  if (TrackTime) {
    TH3D *timehist =
        IsPrimary ? Deposits_ChrgWSumTime : DaughterDeposits_ChrgWSumTime;

    timehist->Fill(x[0], x[1], x[2], x[3] * edep);
  }
}

void DepoParticle::AddDeposit(DepoParticle &other) {
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

void DepoParticle::Reset() {
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

void DepoParticle::FinalizeTime() {
  if (TrackTime) {
    Deposits_ChrgWSumTime->Divide(Deposits);
    DaughterDeposits_ChrgWSumTime->Divide(DaughterDeposits);
  }
}

DepoParticle::~DepoParticle() {
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


DepoTracked::DepoTracked(int PDG, size_t trackID, int NMaxTrackSteps,
            double timesep_us)
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


void DepoTracked::Swap(DepoTracked &&other) {
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

DepoTracked::DepoTracked(DepoTracked &&other) : DepoParticle(std::move(other)) {
  Swap(std::move(other));
}

DepoTracked &DepoTracked::operator=(DepoTracked &&other) {
  DepoParticle::Swap(std::move(other));
  Swap(std::move(other));
  return (*this);
}

void DepoTracked::AddStep(double *x, double *p) {
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

DepoTracked::~DepoTracked() {
  delete[] _Position;
  delete[] _Momentum;
  delete[] Position;
  delete[] Momentum;
}
