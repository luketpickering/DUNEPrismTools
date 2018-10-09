#ifndef FLUXCOMBINER_SEEN
#define FLUXCOMBINER_SEEN

#ifdef USE_FHICL

#include "TH1D.h"

#include "fhiclcpp/ParameterSet.h"

#include <memory>

std::unique_ptr<TH1D> GetCombinedFlux(fhicl::ParameterSet const &,
                                      std::vector<std::unique_ptr<TH1D>> &);

inline std::unique_ptr<TH1D> GetCombinedFlux(fhicl::ParameterSet const &ps) {
  std::vector<std::unique_ptr<TH1D>> dummy;
  return GetCombinedFlux(ps, dummy);
}
#endif

#endif
