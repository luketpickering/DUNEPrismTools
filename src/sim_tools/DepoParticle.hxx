#ifndef DEPOPARTICLE_HXX_SEEN
#define DEPOPARTICLE_HXX_SEEN

#include "TH3D.h"

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

  DepoParticle();

  void SetTrackTime();

  DepoParticle(int PDG, size_t trackID, double timesep_us = 0xdeadbeef);

  DepoParticle(DepoParticle const &) = delete;
  DepoParticle &operator=(DepoParticle const &) = delete;

  void Swap(DepoParticle &&other);

  DepoParticle(DepoParticle &&other);

  DepoParticle &operator=(DepoParticle &&other);

  void LendDepositMaps(TH3D *map1, TH3D *map2);

  void DisownDepositMaps();
  void DeleteDepositMaps();

  void AddDeposit(double *x, double edep, bool IsPrimary);

  void AddDeposit(DepoParticle &other);

  void Reset();

  void FinalizeTime();

  virtual ~DepoParticle();
};

struct DepoTracked : public DepoParticle {
  size_t kMaxTrackedSteps;
  double *_Position;
  double *_Momentum;

  double **Position;
  double **Momentum;
  size_t NSteps;

  DepoTracked(int PDG, size_t trackID, int NMaxTrackSteps = 1000,
              double timesep_us = 0xdeadbeef);

  DepoTracked(DepoTracked const &) = delete;
  DepoTracked &operator=(DepoTracked const &) = delete;

  void Swap(DepoTracked &&other);

  DepoTracked(DepoTracked &&other);

  DepoTracked &operator=(DepoTracked &&other);

  void AddStep(double *x, double *p);

  virtual ~DepoTracked();
};

#endif
