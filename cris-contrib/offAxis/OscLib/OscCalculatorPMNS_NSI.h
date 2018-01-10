#include "OscCalculatorPMNS_NSI.h"

#include <cassert>
#include <cstdlib>

namespace osc
{
  OscCalculatorPMNS_NSI::OscCalculatorPMNS_NSI()
    : fMixDirty(true), fDmDirty(true), fPropDirty(true), fEpsDirty(true), fPrevAnti(0)
  {
  }

  OscCalculatorPMNS_NSI::~OscCalculatorPMNS_NSI()
  {
  }

  double OscCalculatorPMNS_NSI::P(int flavBefore, int flavAfter, double E)
  {
    const int anti = (flavBefore > 0) ? +1 : -1;
    assert(flavAfter/anti > 0);
    if(anti != fPrevAnti) fPropDirty = true;

    int i = -1, j = -1;
    if(abs(flavBefore) == 12) i = 0;
    if(abs(flavBefore) == 14) i = 1;
    if(abs(flavBefore) == 16) i = 2;
    if(abs(flavAfter) == 12) j = 0;
    if(abs(flavAfter) == 14) j = 1;
    if(abs(flavAfter) == 16) j = 2;
    assert(i >= 0 && j >= 0);

    if(fMixDirty){
      fPMNS_NSI.SetMix(fTh12, fTh23, fTh13, fdCP);
      fMixDirty = false;
    }
    if(fDmDirty){
      fPMNS_NSI.SetDeltaMsqrs(fDmsq21, fDmsq32);
      fDmDirty = false;
    }
    if(fEpsDirty){
      fPMNS_NSI.SetNSI(fEps_ee,    fEps_emu,    fEps_etau,
                       fEps_mumu,  fEps_mutau,  fEps_tautau,
                       fDelta_emu, fDelta_etau, fDelta_mutau);
      fEpsDirty = false;
    }


    fPMNS_NSI.ResetToFlavour(i);
    // Assume Z/A=0.5
    const double Ne = fRho/2;
    fPMNS_NSI.PropMatter(fL, E, Ne, anti);
    return fPMNS_NSI.P(j);
  }
} // namespace
