#include "OscCalculatorPMNSOpt.h"

#include <cassert>
#include <cstdlib>

namespace osc
{
  OscCalculatorPMNSOpt::OscCalculatorPMNSOpt()
    : fMixDirty(true), fDmDirty(true), fPropDirty(true), fPrevAnti(0)
  {
  }

  OscCalculatorPMNSOpt::~OscCalculatorPMNSOpt()
  {
  }

  double OscCalculatorPMNSOpt::P(int flavBefore, int flavAfter, double E)
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
      fPMNSOpt.SetMix(fTh12, fTh23, fTh13, fdCP);
      fMixDirty = false;
    }
    if(fDmDirty){
      fPMNSOpt.SetDeltaMsqrs(fDmsq21, fDmsq32);
      fDmDirty = false;
    }


    fPMNSOpt.ResetToFlavour(i);
    // Assume Z/A=0.5
    const double Ne = fRho/2;
    fPMNSOpt.PropMatter(fL, E, Ne, anti);
    return fPMNSOpt.P(j);
  }
} // namespace


