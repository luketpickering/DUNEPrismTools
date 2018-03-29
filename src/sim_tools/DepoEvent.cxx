#include "DepoEvent.hxx"

#include <iostream>
#include <sstream>

PrimaryParticle::PrimaryParticle()
    : IsFinalState(false), PDG(0), EKin(0), EMass(0), ThreeMom(0, 0, 0) {}

std::string PrimaryParticle::ToString() {
  std::stringstream ss("");

  ss << "PDG = " << PDG << ", E = " << (EKin + EMass) << ", 3Mom = {"
     << ThreeMom.X() << ", " << ThreeMom.Y() << ", " << ThreeMom.Z() << "}"
     << ", IsFinalState: " << IsFinalState << std::flush;
  return ss.str();
}

void DepoEvent::PrintGENIEPassthrough() {
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
DepoParticle *DepoEvent::GetPrimaryParticle(size_t parent_id, size_t TrackID) {
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
bool DepoEvent::GetIsPrimary(size_t TrackID) { return IsPrimary[TrackID]; }

DepoTracked *DepoEvent::GetTrackedDeposit(size_t TrackID) {
  if (TrackedParticleMap.count(TrackID)) {
    return TrackedParticleMap[TrackID];
  }
  return nullptr;
}
